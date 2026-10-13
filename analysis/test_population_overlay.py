"""Check cached fractions against direct geometry on offset native grids."""
from pathlib import Path
import tempfile,os
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
import numpy as np,rasterio
from rasterio.transform import from_origin
from shapely.geometry import box
from population_overlay import geometry,union_population_support,matrix_for_area,apply_counts
from spatial_weights import exposure_geometry,allocate_population
root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(dir=root/'qa') as temp:
    cache=Path(temp)
    fine={'driver':'GTiff','height':3,'width':4,'count':1,'dtype':'float32','crs':'EPSG:4326','transform':from_origin(30,-19,.01,.01),'nodata':-9999}
    coarse={'height':3,'width':3,'crs':'EPSG:4326','transform':from_origin(29.995,-18.995,.017,.017)}
    populations=[np.array([[0,8,12,0],[6,0,20,5],[0,9,0,11]],dtype='float32'),np.array([[4,9,10,0],[5,7,15,8],[1,6,3,9]],dtype='float32')]
    paths=[]
    for i,p in enumerate(populations):
        path=cache/f'population_{i}.tif'
        with rasterio.open(path,'w',**fine) as dst:dst.write(p,1)
        paths.append(path)
    boundary=box(30.006,-19.028,30.037,-19.003)
    geos=geometry(coarse,{'test':boundary});support=union_population_support(paths,cache)
    matrix=matrix_for_area(support,geos['test'],coarse,boundary,cache,'test',chunk=4)
    direct=exposure_geometry(coarse,boundary)
    for p in populations:
        w,c,qa=apply_counts(matrix,p)
        expected,expected_c,_=allocate_population(p,fine,direct)
        np.testing.assert_allclose(w,expected,rtol=0,atol=1e-10)
        np.testing.assert_allclose(c,expected_c,rtol=0,atol=1e-10)
        assert qa['allocation_absolute_error']<1e-10
    cached=matrix_for_area(support,geos['test'],coarse,boundary,cache,'test',chunk=4)
    assert cached==matrix
print('Cached population overlay matches direct geometry for offset grids, partial borders and changing support',flush=True)
