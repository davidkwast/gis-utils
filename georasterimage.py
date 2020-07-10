import shutil
import json
import math

from pathlib import Path

import rasterio
import rasterio.windows
import rasterio.mask
from osgeo import gdal;gdal.UseExceptions()
from osgeo import ogr;ogr.UseExceptions()
from osgeo import osr;osr.UseExceptions()


wgs84_spatial_ref = osr.SpatialReference()
wgs84_spatial_ref.ImportFromEPSG(4326)

def get_wgs84_transform(crs_epsg):
    source = osr.SpatialReference()
    source.ImportFromEPSG(crs_epsg)
    
    return osr.CreateCoordinateTransformation(source, wgs84_spatial_ref)

def get_from_wgs84_transform(crs_epsg):
    target = osr.SpatialReference()
    target.ImportFromEPSG(crs_epsg)
    
    return osr.CreateCoordinateTransformation(wgs84_spatial_ref, target)


def convert_wgs_to_utm(lat: float, lon: float) -> int:
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
    else:
        epsg_code = '327' + utm_band
    return int(epsg_code)


class GeoRasterImageException(Exception): pass

class GeoRasterImageOpenError(GeoRasterImageException): pass

class GeoRasterImageReadError(GeoRasterImageOpenError): pass

class GeoRasterImageCropError(GeoRasterImageOpenError): pass


class GeoRasterImage:
    
    def __init__(self, file_path=''):
        
        if isinstance(file_path, Path):
            self.file_path = file_path
        elif isinstance(file_path, str):
            self.file_path = Path(file_path)
        else:
            raise ValueError('"file_path" must be "str" or "pathlib.Path" instance')
        
        self._dataset = None
        self.profile = None
        self.transform = None
        self.crs = None
        self.bands_count = 0
        self.dtypes = (None,)
        self.bounds = None
    
    def open(self, mode='r', **kw):
        if self._dataset:
            raise GeoRasterImageOpenError('dataset already opened')
        if mode not in ('r','w'):
            raise GeoRasterImageOpenError('invalid mode')
        # GTiff performance and memory care
        if mode == 'w' and 'driver' in kw and kw['driver'] == 'GTiff':
            kw['compress'] = 'deflate'
            kw['zlevel'] = '1'
            kw['num_threads'] = '2'
            kw['interleave'] = 'band'
        
        # loads attrs - TODO: reload method
        self._dataset = rasterio.open(self.file_path, mode, **kw)
        self.profile = self._dataset.profile
        self.transform = self._dataset.transform
        self.crs = self._dataset.crs
        self.bands_count = self._dataset.count
        self.dtypes =  self._dataset.dtypes
        self.bounds =  self._dataset.bounds
    
    def write(self, bands_iterable):
        for band, array in enumerate(bands_iterable):
            self._dataset.write(array, band+1)
    
    def read(self, **kwargs):
        if not self._dataset:
            raise GeoRasterImageReadError('dataset not opened')
        return self._dataset.read(**kwargs)
    
    def close(self):
        self._dataset.close()
        self._dataset = None
    
    
    def get_xy_size(self):
        return self._dataset.width, self._dataset.height
    
    def get_array_shape(self):
        return self._dataset.height, self._dataset.width
    
    def get_pixel_size(self):
        return abs(self.profile['transform'][0])
    
    def get_spatial_coords(self, array_x, array_y):
        return self._dataset.xy(array_x, array_y)
    
    def get_array_coords(self, spatial_x, spatial_y, roundind_method='floor'):
        'ref: https://rasterio.readthedocs.io/en/stable/api/rasterio.io.html#rasterio.io.DatasetReader.index'
        if roundind_method == 'floor':
            return self._dataset.index(spatial_x, spatial_y)
        elif roundind_method == 'round':
            return self._dataset.index(spatial_x, spatial_y, op=round)
    
    
    def generate_windows(self, block_size=(5000,5000), step_ratio=(0.9,0.9)):
        block_size_x, block_size_y = block_size
        size_x, size_y = self.get_xy_size()
        step_x, step_y = (int(s*v) for s,v in zip(block_size,step_ratio))
        for x in range(0, size_x, step_x):
            for y in range(0, size_y, step_y):
                yield rasterio.windows.Window(x,y,block_size_x, block_size_y)
    
    
    def process_profile__from_window(self, window):
        new_profile = self.profile
        new_profile.update({
            'height': window.height,
            'width': window.width,
            'transform': rasterio.windows.transform(window, self.transform),
        })
        return new_profile
    
    
    def crop__by_shape(self, shape_geojson: str, new_img_path: Path):
        
        dataset = self._dataset
        profile = self.profile
        
        shape_geom = ogr.CreateGeometryFromJson(shape_geojson)
        if shape_geom is None:
            raise GeoRasterImageCropError('Invalid shape:', shape_geojson[:40], '...')
        shape_geom.Transform( get_from_wgs84_transform( self.crs.to_epsg() ) )
        shape_geom_str = shape_geom.ExportToJson()
        shape_geojson_dict = json.loads(shape_geom_str)
        
        indexes = [x+1 for x in range(self.bands_count)]
        
        try:
            arrays, affine = rasterio.mask.mask(
                dataset,
                [shape_geojson_dict],
                indexes = indexes,
                crop = True,
                all_touched = False,
                invert = False,
                nodata = profile['nodata'],#0,#np.nan,
                filled = True,
                pad = False,
            )
        
        except ValueError as error:
            raise GeoRasterImageCropError(str(error))
        
        new_img = GeoRasterImage(new_img_path)
        new_profile = profile.copy()
        
        new_profile['height'],new_profile['width'] = arrays[0].shape
        new_profile['transform'] = affine
        
        new_img.open('w', **new_profile)
        new_img.write(arrays)
        new_img.close()
        
        return new_img
    
    
    def polygonize(self, shape_dir_path: Path, select_value=0, buffer_treatment=None):
        if shape_dir_path.exists():
            shutil.rmtree(shape_dir_path)
        shape_dir_path.mkdir(parents=True, exist_ok=False)
        layer_name = shape_dir_path.name
        
        osr_srs = osr.SpatialReference()
        osr_srs.ImportFromEPSG(self.crs.to_epsg())
        
        src_ds = gdal.Open(str(self.file_path))
        srcband = src_ds.GetRasterBand(1)
        
        mem_drv = ogr.GetDriverByName('MEMORY')
        dst_ds = mem_drv.CreateDataSource('')
        dst_layer = dst_ds.CreateLayer("L1", srs=osr_srs)
        dst_layer.CreateFields([ogr.FieldDefn('DN')])
        
        gdal.Polygonize(srcband, None, dst_layer, 0)
        dst_layer.ResetReading()
        
        filtered_polygons = []
        for feature in dst_layer:
            value = feature.GetFieldAsInteger(0)
            if value == select_value:
                geom = feature.GetGeometryRef()#.Clone()
                
                if buffer_treatment is not None:
                    geom = geom.Buffer(distance= -1*buffer_treatment, quadsecs=1)
                    geom = geom.Buffer(distance= buffer_treatment, quadsecs=1)
                else:
                    geom = geom.Clone()
                
                filtered_polygons.append(geom)
        
        shp_drv = ogr.GetDriverByName('ESRI Shapefile')
        shp_ds = shp_drv.CreateDataSource(str(shape_dir_path))
        shp_lyr = shp_ds.CreateLayer(layer_name, srs=osr_srs)
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        shp_lyr.CreateField(idField)
        
        for index, polygon in enumerate(filtered_polygons):
            feat_def = shp_lyr.GetLayerDefn()
            feat = ogr.Feature(feat_def)
            feat.SetGeometry(polygon)
            feat.SetField("id", index)
            shp_lyr.CreateFeature(feat)
        
        shp_ds.SyncToDisk()
        
        geojson_dict = {'type': 'GeometryCollection', 'geometries': []} # type: ignore
        for polygon in filtered_polygons:
            geom_wgs84 = polygon.Clone()
            geom_wgs84.Transform( get_wgs84_transform( self.crs.to_epsg() ) )
            geom_wgs84_str = geom_wgs84.ExportToJson()
            geojson_dict['geometries'].append(json.loads(geom_wgs84_str)) # type: ignore
        
        with open( shape_dir_path/(layer_name+'.json'), 'w' ) as fd:
            fd.write(json.dumps(geojson_dict, indent=1))
        
        return filtered_polygons, geojson_dict
