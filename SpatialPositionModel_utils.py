# -*- coding: utf-8 -*-
"""
SpatialPositionModel Utils
"""
import numpy as np
from qgis.core import QgsGeometry, QgsPoint
from matplotlib.pyplot import contourf
from matplotlib.mlab import griddata


def parse_expression(expr):
    nexpr = []
    dico = {
        '*': np.multiply,
        '-': np.subtract,
        '/': np.divide,
        '+': np.sum
        }
    for i, char in enumerate(expr):
        nexpr.append(char)
        if char in ('/', '+', '*', '-', '(', ')'):
            try:
                if expr[i+1] != ' ':
                    nexpr.append(' ')
            except IndexError:
                pass
    nexpr = ''.join(nexpr).split(' ')
    if len(nexpr) == 3:
        fields = nexpr[:3:2]
        if "\"" in fields[0][0] or "'" in fields[0][0]:
            fields[0] = fields[0][1:-1]
        if "\"" in fields[1][0] or "'" in fields[1][0]:
            fields[1] = fields[1][1:-1]
        return fields, dico[nexpr[1]]
    else:
        return -1
        

def hav_dist(locs1, locs2, k=np.pi/180):
    # (lat, lon)
    locs1 = locs1 * k
    locs2 = locs2 * k
    cos_lat1 = np.cos(locs1[..., 0])
    cos_lat2 = np.cos(locs2[..., 0])
    cos_lat_d = np.cos(locs1[..., 0] - locs2[..., 0])
    cos_lon_d = np.cos(locs1[..., 1] - locs2[..., 1])
    return 6367 * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))


def make_dist_mat(xy1, xy2, longlat=False):
    if not longlat:
        d0 = np.subtract.outer(xy1[:,0], xy2[:,0])
        d1 = np.subtract.outer(xy1[:,1], xy2[:,1])
        return np.hypot(d0, d1)
    else:
        return hav_dist(xy1, xy2)


def compute_interact_density(matdist, typefun, beta, span):
    if 'pareto' in typefun:
        alpha  = (2 ** (1 / beta) - 1) / span
        matDens = (1 + alpha * matdist) ** (-beta)
    elif 'exponential' in typefun:
        alpha  = np.log(2) / span ** beta
        matDens = np.exp(- alpha * matdist ** beta)
    else:
        raise ValueError('Bad interaction function argument: {}'
                         .format(typefun))
    return matDens.round(8)


def compute_opportunity(pts_values, matdens):
    matOpport = pts_values[:, np.newaxis] * matdens
    return matOpport.round(8)

def compute_potentials(matopport):
    return matopport.sum(axis=0)

def compute_reilly(mattoportes):
    return np.argmax(mattoportes, axis=0)

def compute_huff(unknownpts, matopport):
    sum_lines = matopport.sum(1)
    sum_lines = sum_lines[np.where(sum_lines > 0)[0]]
    matopportPct = np.array([100]) * (matopport[np.where(sum_lines > 0)[0]] / sum_lines[:, np.newaxis]).T
    return matopportPct.max(axis=1)

def make_regular_points(bounds, resolution):
    """
    Return a grid of regular points.
    """
    xmin, ymin, xmax, ymax = bounds
    nb_x = int(round((xmax - xmin) / resolution))
    nb_y = int(round((ymax - ymin) / resolution))
    w = (xmax - xmin) + (xmax - xmin) / 5
    h = (ymax - ymin) + (ymax - ymin) / 5
    prog_x = [(xmin - (xmax - xmin) / 5) + i * (w / nb_x) for i in xrange(nb_x+1)]
    prog_y = [(ymin - (ymax - ymin) / 5) + i * (h / nb_y) for i in xrange(nb_y+1)]
    return np.array([(x, y) for x in prog_x for y in prog_y])


def render_stewart(pot, unknownpts, nb_class):
    x = np.array([c[0] for c in unknownpts])
    y = np.array([c[1] for c in unknownpts])
    xi = np.linspace(np.nanmin(x), np.nanmax(x), len(x))
    yi = np.linspace(np.nanmin(y), np.nanmax(y), len(y))
    zi = griddata(x, y, pot, xi, yi, interp='linear')
    collec_poly = contourf(
        xi, yi, zi, nb_class, vmax=abs(zi).max(), vmin=-abs(zi).max())
    _ = [0] + [100/i for i in range(1, nb_class)][::-1]
    break_values = np.percentile(pot[pot.nonzero()[0]], q=_)
    break_values = np.append(np.array([0]), break_values)
    res_poly = qgsgeom_from_mpl_collec(collec_poly.collections, levels=break_values.tolist())
#    if mask is not None:
#        res.geometry = res.geometry.apply(lambda x: x.intersection(unary_union(mask.geometry)))
#        res = GeoDataFrame(geometry = [i for i in mask.geometry.values] + [i for i in res.geometry.values])
    return res_poly, break_values

def render_reilly(reilly_values, unknownpts):
    x = np.array([c[0] for c in unknownpts])
    y = np.array([c[1] for c in unknownpts])
    xi = np.linspace(np.nanmin(x), np.nanmax(x), round(np.sqrt(len(x))))
    yi = np.linspace(np.nanmin(y), np.nanmax(y), round(np.sqrt(len(y))))
    zi = griddata(x, y, reilly_values, xi, yi, interp='linear')
    nb_class = len({ft for ft in reilly_values})
    collec_poly = contourf(
        xi, yi, zi, nb_class, vmax=abs(zi).max(), vmin=-abs(zi).max())
    _ = [0] + [100/i for i in range(1, nb_class)][::-1]
    break_values = np.percentile(reilly_values[reilly_values.nonzero()[0]], q=_)
    break_values = np.append(np.array([0]), break_values)
    res_poly = qgsgeom_from_mpl_collec(collec_poly.collections,
                                       levels=break_values.tolist())
#    if mask is not None:
#        res.geometry = res.geometry.apply(lambda x: x.intersection(unary_union(mask.geometry)))
#        res = GeoDataFrame(geometry = [i for i in mask.geometry.values] + [i for i in res.geometry.values])
    return res_poly, break_values


def qgsgeom_from_mpl_collec(collections, levels):
    polygons = []
    for i, polygon in enumerate(collections):
        mpoly = []
        for path in polygon.get_paths():
            path.should_simplify = False
            poly = path.to_polygons()
            if len(poly) > 0 and len(poly[0]) > 4:
                exterior = [QgsPoint(*p.tolist()) for p in poly[0]]
                holes = [
                    [QgsPoint(*p.tolist()) for p in h]
                    for h in poly[1:] if len(h) > 4
                    ]
                if len(holes) == 1:
                    mpoly.append([exterior, holes[0]])
                elif len(holes) > 1:
                    mpoly.append([exterior] + [h for h in holes])
                else:
                    mpoly.append([exterior, holes])
        if len(mpoly) > 1:
            polygons.append(QgsGeometry.fromMultiPolygon(mpoly))
        elif len(mpoly) == 1:
            polygons.append(QgsGeometry.fromPolygon(mpoly[0]))
    return polygons
