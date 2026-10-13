"""Conserve population counts when overlaying different native raster grids."""
import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.ops import transform as transform_geometry
from pyproj import Transformer


def clip_native_cells(polygons, boundary):
    """Intersect only edge cells; prepared predicates spare interior cells.

    This is mathematically the same overlay as intersecting every footprint,
    but avoids repeatedly processing a detailed boundary for interior cells.
    """
    shapely.prepare(boundary)
    touches=shapely.intersects(boundary,polygons)
    contained=shapely.contains(boundary,polygons)
    clipped=polygons.copy()
    clipped[~touches]=Polygon()
    edge=touches&~contained
    clipped[edge]=shapely.intersection(polygons[edge],boundary)
    return clipped


def native_pixel_polygons(transform, rows, columns, crs):
    """Pixel footprints in a common equal-area CRS, EPSG:6933."""
    rows, columns = np.asarray(rows), np.asarray(columns)
    corners = np.stack([
        np.column_stack((columns, rows)), np.column_stack((columns+1, rows)),
        np.column_stack((columns+1, rows+1)), np.column_stack((columns, rows+1)),
        np.column_stack((columns, rows)),
    ], axis=1)
    x = transform.a * corners[:, :, 0] + transform.b * corners[:, :, 1] + transform.c
    y = transform.d * corners[:, :, 0] + transform.e * corners[:, :, 1] + transform.f
    project = Transformer.from_crs(crs, 6933, always_xy=True)
    xx, yy = project.transform(x, y)
    return shapely.polygons(np.stack((xx, yy), axis=-1))


def exposure_geometry(profile, boundary_wgs84):
    rr, cc = np.indices((profile['height'], profile['width']))
    polygons = native_pixel_polygons(profile['transform'], rr.ravel(), cc.ravel(), profile['crs'])
    tx = Transformer.from_crs(4326, 6933, always_xy=True).transform
    boundary = transform_geometry(tx, boundary_wgs84)
    clipped = clip_native_cells(polygons, boundary)
    area = shapely.area(clipped)
    inside = area > 0
    return {'polygons': polygons, 'clipped': clipped, 'area_m2': area, 'inside': inside,
            'tree': shapely.STRtree(polygons), 'boundary': boundary}


def allocate_population(population, population_profile, exposure):
    """Area-apportion each source count to the study boundary and exposure cells.

    Positive population cells are treated as uniformly occupied within their native
    footprint. This assumption is explicit and retained at administrative edges.
    Return conservative weights and a centre-assignment sensitivity result.
    """
    rr, cc = np.where(np.isfinite(population) & (population > 0))
    counts = population[rr, cc].astype('float64')
    footprints = native_pixel_polygons(population_profile['transform'], rr, cc, population_profile['crs'])
    areas = shapely.area(footprints)
    clipped = clip_native_cells(footprints, exposure['boundary'])
    valid = shapely.area(clipped) > 0
    counts, areas, footprints, clipped = counts[valid], areas[valid], footprints[valid], clipped[valid]
    fraction = shapely.area(clipped) / areas
    expected = float(np.sum(counts * fraction))
    weights = np.zeros(len(exposure['polygons']), dtype='float64')
    tree = exposure['tree']
    left, right = tree.query(clipped, predicate='intersects')
    overlap = shapely.area(shapely.intersection(clipped[left], exposure['polygons'][right]))
    contribution = counts[left] * overlap / areas[left]
    np.add.at(weights, right, contribution)
    # Each whole or boundary-clipped count is assigned to the nearest native-cell
    # centroid only for sensitivity; primary results use polygon intersections.
    centre_weights = np.zeros_like(weights)
    centres = shapely.centroid(footprints)
    nearest = tree.nearest(centres)
    np.add.at(centre_weights, nearest, counts*fraction)
    error = abs(weights.sum() - expected)
    if error > max(0.01, expected*1e-8):
        raise AssertionError(f'Population mass conservation failed: error={error}')
    return weights, centre_weights, {
        'population_in_boundary': expected,
        'allocated_population': float(weights.sum()),
        'allocation_absolute_error': float(error),
        'positive_source_cells_intersecting_boundary': int(valid.sum()),
        'method': 'Equal-area intersections of native pixel footprints; EPSG:6933',
        'within_cell_assumption': 'Uniform distribution within each native population cell',
    }
