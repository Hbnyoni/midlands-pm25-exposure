"""Reusable, count-conserving population overlays for multiple years and areas.

Sparse matrices store geometric fractions, not reprojected population densities.
All six years share a source grid. Each matrix is built once on their union of
positive cells and then applied to actual year-specific counts.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import rasterio,shapely
from shapely.geometry import shape
from shapely.ops import transform
from pyproj import Transformer
from scipy.sparse import csr_matrix,save_npz,load_npz
from spatial_weights import native_pixel_polygons,clip_native_cells

def geometry(profile,boundaries):
    rr,cc=np.indices((profile['height'],profile['width']))
    polygons=native_pixel_polygons(profile['transform'],rr.ravel(),cc.ravel(),profile['crs'])
    tree=shapely.STRtree(polygons)
    tx=Transformer.from_crs(4326,6933,always_xy=True).transform
    out={}
    for name,boundary in boundaries.items():
        b=transform(tx,boundary);clip=clip_native_cells(polygons,b);area=shapely.area(clip)
        out[name]={'polygons':polygons,'tree':tree,'boundary':b,'inside':area>0,'area_m2':area}
    return out

def union_population_support(paths,cache):
    """Read the public rasters one at a time to keep memory bounded."""
    mask=None;reference=None
    for path in paths:
        with rasterio.open(path) as src:
            a=src.read(1)
            if reference is None:
                reference=src.profile.copy();mask=np.zeros(a.shape,dtype=bool)
            assert src.transform==reference['transform'] and a.shape==mask.shape
            assert src.crs==reference['crs']
            mask|=np.isfinite(a)&(a>0)
    rr,cc=np.where(mask)
    support={'rows':rr.astype('int32'),'cols':cc.astype('int32'),'profile':reference,
             'hash':hashlib.sha256(mask.tobytes()).hexdigest()}
    np.savez_compressed(cache/'population_positive_support.npz',rows=support['rows'],cols=support['cols'])
    return support

def matrix_for_area(support,exposure,exposure_profile,boundary_wgs84,cache,label,chunk=100_000):
    safe=''.join(c if c.isalnum() else '_' for c in label)
    primary=cache/(safe+'_fractional.npz');centre=cache/(safe+'_centre.npz')
    boundary_file=cache/(safe+'_boundary_fractions.npz');meta=cache/(safe+'_metadata.json')
    profile=support['profile'];rr=support['rows'];cc=support['cols']
    signature={'population_transform':list(profile['transform'])[:6],'population_shape':[profile['height'],profile['width']],
               'exposure_transform':list(exposure_profile['transform'])[:6],'exposure_shape':[exposure_profile['height'],exposure_profile['width']],
               'positive_support_sha256':support['hash'],'boundary_wkb_sha256':hashlib.sha256(boundary_wgs84.wkb).hexdigest(),
               'method':'EPSG6933 native footprint intersections v1'}
    if all(f.exists() for f in [primary,centre,boundary_file,meta]) and json.loads(meta.read_text())==signature:
        return {'primary':primary,'centre':centre,'fractions':boundary_file}
    t=profile['transform'];x0,y0,x1,y1=boundary_wgs84.bounds
    cmin=int(np.floor((x0-t.c)/t.a))-1;cmax=int(np.ceil((x1-t.c)/t.a))+1
    rmin=int(np.floor((y1-t.f)/t.e))-1;rmax=int(np.ceil((y0-t.f)/t.e))+1
    select=(rr>=rmin)&(rr<=rmax)&(cc>=cmin)&(cc<=cmax)
    rr=rr[select];cc=cc[select]
    rid=[];cid=[];vals=[];cent_rows=[];cent_cols=[];fractions=[]
    expected_cols=[];expected_fraction=[]
    ncell=len(exposure['polygons']);nsource=profile['height']*profile['width']
    for start in range(0,len(rr),chunk):
        r=rr[start:start+chunk];c=cc[start:start+chunk]
        indices=r.astype('int64')*profile['width']+c
        full=native_pixel_polygons(t,r,c,profile['crs']);area=shapely.area(full)
        clipped=clip_native_cells(full,exposure['boundary']);cliparea=shapely.area(clipped)
        use=cliparea>0
        full,clipped,area,cliparea,indices=full[use],clipped[use],area[use],cliparea[use],indices[use]
        if not len(indices):continue
        left,right=exposure['tree'].query(clipped,predicate='intersects')
        overlap=shapely.area(shapely.intersection(clipped[left],exposure['polygons'][right]))
        value=overlap/area[left];positive=value>0
        rid.append(right[positive]);cid.append(indices[left[positive]]);vals.append(value[positive])
        frac=cliparea/area
        near=exposure['tree'].nearest(shapely.centroid(full))
        cent_rows.append(near);cent_cols.append(indices);fractions.append(frac)
        expected_cols.append(indices);expected_fraction.append(frac)
    assert vals,label
    matrix=csr_matrix((np.concatenate(vals),(np.concatenate(rid),np.concatenate(cid))),shape=(ncell,nsource))
    cmatrix=csr_matrix((np.concatenate(fractions),(np.concatenate(cent_rows),np.concatenate(cent_cols))),shape=(ncell,nsource))
    indices=np.concatenate(expected_cols);fraction=np.concatenate(expected_fraction)
    # Geometry conservation is tested column by column without assuming counts.
    by_source=np.asarray(matrix.sum(axis=0)).ravel()[indices]
    assert np.max(np.abs(by_source-fraction))<1e-7,(label,np.max(np.abs(by_source-fraction)))
    assert np.all(np.asarray(matrix.sum(axis=1)).ravel()[~exposure['inside']]<1e-10)
    save_npz(primary,matrix);save_npz(centre,cmatrix)
    np.savez_compressed(boundary_file,source_index=indices,area_fraction=fraction)
    meta.write_text(json.dumps(signature,indent=2))
    print('Cached population overlay',label,'source cells',len(indices),'intersections',matrix.nnz,flush=True)
    return {'primary':primary,'centre':centre,'fractions':boundary_file}

def apply_counts(paths,population):
    flat=np.nan_to_num(population.astype('float64',copy=False),nan=0,copy=False).ravel();flat[flat<0]=0
    w=load_npz(paths['primary'])@flat;c=load_npz(paths['centre'])@flat
    with np.load(paths['fractions']) as f:
        expected=float(np.sum(flat[f['source_index']]*f['area_fraction']))
        positive=int(np.sum(flat[f['source_index']]>0))
    error=abs(w.sum()-expected)
    assert error<max(.01,expected*1e-8),error
    assert abs(c.sum()-expected)<max(.01,expected*1e-8)
    return w,c,{'population_in_boundary':expected,'allocated_population':float(w.sum()),'allocation_absolute_error':float(error),
                'positive_source_cells_intersecting_boundary':positive,'method':'Cached EPSG6933 native-footprint fractions applied to year-specific counts',
                'within_cell_assumption':'Uniform distribution within each native population cell'}
