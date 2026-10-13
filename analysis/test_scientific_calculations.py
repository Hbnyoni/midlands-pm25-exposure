"""Analytic edge-case checks for consequential numerical operations."""
import os
os.environ.pop('PROJ_LIB',None)
os.environ.pop('PROJ_DATA',None)
import numpy as np
from rasterio.transform import from_origin
from shapely.geometry import box
from risk_core import who_band,lag01_reduction,transferred_rr,longest_exceedance_run
from spatial_weights import exposure_geometry,allocate_population

assert np.array_equal(who_band([5,5.001,10,15,25,35,35.001,np.nan],[5,10,15,25,35]),[0,1,1,2,3,4,5,-1])
x=np.array([10,15,20,30,np.nan,40,15])[:,None]
expected=np.array([np.nan,0,2.5,10,np.nan,np.nan,12.5])[:,None]
np.testing.assert_allclose(lag01_reduction(x),expected,equal_nan=True)
assert longest_exceedance_run(x)[0]==2
assert transferred_rr(10,1.136)==1.136
assert transferred_rr(0,1.136)==1
# One population cell spanning two equal-width exposure cells: 100 people
# must split 50:50; clipping the study to one half must retain exactly 50.
pr={'height':1,'width':2,'transform':from_origin(30,-19,.01,.01),'crs':'EPSG:4326'}
pp={'height':1,'width':1,'transform':from_origin(30,-19,.02,.01),'crs':'EPSG:4326'}
g=exposure_geometry(pr,box(30,-19.01,30.02,-19))
w,_,q=allocate_population(np.array([[100.]]),pp,g)
np.testing.assert_allclose(w,[50,50],atol=1e-7)
g2=exposure_geometry(pr,box(30,-19.01,30.01,-19))
w2,_,q2=allocate_population(np.array([[100.]]),pp,g2)
np.testing.assert_allclose(w2,[50,0],atol=1e-7)
print('PASS: exact WHO cutpoints, lag missingness, RR reference, exceedance runs, population splitting and boundary clipping')
