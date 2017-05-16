# -*- coding: utf-8 -*-
"""
SpatialPositionModel Utils
"""
import os
import numpy as np
from osgeo import gdal, osr
from PyQt4.QtGui import QColor
from qgis.core import (
    QgsGeometry, QgsPoint, QgsVectorGradientColorRampV2, QgsFeature,
    QgsGraduatedSymbolRendererV2, QgsFillSymbolV2, QgsRendererRangeV2,
    QgsRasterShader, QgsColorRampShader, QgsSingleBandPseudoColorRenderer,
    QgsRasterBandStats)
from shapely.wkb import loads
from shapely.ops import cascaded_union
from tempfile import mkdtemp
from uuid import uuid4


def parse_expression(expr):
    nexpr = []
    dico = {
        '*': np.multiply, '-': np.subtract,
        '/': np.divide, '+': np.sum
        }
    for i, char in enumerate(expr):
        nexpr.append(char)
        if char in ('/', '+', '*', '-', '(', ')'):
            try:
                if nexpr[-2] != ' ':
                    _ix = len(nexpr) - 1
                    nexpr.insert(_ix, ' ')
                if expr[i + 1] != ' ':
                    nexpr.append(' ')
            except IndexError:
                pass
    nexpr = ''.join(nexpr).split(' ')
    for i in nexpr:
        if i.isdigit() and i not in ('*', '+', '-', '/'):
            return -1
    if len(nexpr) == 3:
        fields = nexpr[:3:2]
        if "\"" in fields[0][0] or "'" in fields[0][0]:
            fields[0] = fields[0][1:-1]
        if "\"" in fields[1][0] or "'" in fields[1][0]:
            fields[1] = fields[1][1:-1]
        print(fields, dico[nexpr[1]])
        return fields, dico[nexpr[1]]
    else:
        return -1


def hav_dist(locs1, locs2):
    """
    Return a distance matrix between two set of coordinates.
    Use geometric distance (default) or haversine distance (if longlat=True).

    Parameters
    ----------
    locs1 : numpy.array
        The first set of coordinates as [(long, lat), (long, lat)].
    locs2 : numpy.array
        The second set of coordinates as [(long, lat), (long, lat)].

    Returns
    -------
    mat_dist : numpy.array
        The distance matrix between locs1 and locs2
    """
    locs1 = np.radians(locs1)
    locs2 = np.radians(locs2)
    cos_lat1 = np.cos(locs1[..., 0])
    cos_lat2 = np.cos(locs2[..., 0])
    cos_lat_d = np.cos(locs1[..., 0] - locs2[..., 0])
    cos_lon_d = np.cos(locs1[..., 1] - locs2[..., 1])
    return 6367 * np.arccos(
        cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))


def make_dist_mat(xy1, xy2, longlat=False):
    if not longlat:
        d0 = np.subtract.outer(xy1[:, 0], xy2[:, 0])
        d1 = np.subtract.outer(xy1[:, 1], xy2[:, 1])
        return np.hypot(d0, d1) / 1000
    elif longlat:
        return hav_dist(xy1[:, None], xy2)


def compute_interact_density(matdist, typefun, beta, span):
    if 'pareto' in typefun:
        alpha = (2 ** (1 / beta) - 1) / span
        matDens = (1 + alpha * matdist) ** (-beta)
    elif 'exponential' in typefun:
        alpha = np.log(2) / span ** beta
        matDens = np.exp(- alpha * matdist ** beta)
    else:
        raise ValueError('Bad interaction function argument: {}'
                         .format(typefun))
    return matDens.round(8)


def gen_unknownpts(pts_layer, mask_layer, resolution, longlat):
    if mask_layer:
        ext = mask_layer.extent()
        bounds = (ext.xMinimum(), ext.yMinimum(),
                  ext.xMaximum(), ext.yMaximum())
    else:
        ext = pts_layer.extent()
        bounds = (ext.xMinimum(), ext.yMinimum(),
                  ext.xMaximum(), ext.yMaximum())
        tmp = ((bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
        bounds = (bounds[0] - tmp, bounds[1] - tmp,
                  bounds[2] + tmp, bounds[3] + tmp)

    return make_regular_points(bounds, resolution, longlat)


def parse_class_breaks(class_breaks):
    try:
        values = [float(i.strip().replace(',', '.'))
                  for i in class_breaks.split('-')]
        last = -float('inf')
        for i in values:
            assert i > last
            last = i
        return values
    except:
        return None


def compute_opportunity(pts_values, matdens):
    matOpport = pts_values[:, np.newaxis] * matdens
    return matOpport.round(8)


def compute_potentials(matopport):
    return matopport.sum(axis=0)


def get_height_width(bounds, longlat):
    minlon, minlat, maxlon, maxlat = bounds

    if longlat:
        height = hav_dist(
                np.array([(maxlon + minlon) / 2, minlat]),
                np.array([(maxlon + minlon) / 2, maxlat])
                )
        width = hav_dist(
                np.array([minlon, (maxlat + minlat) / 2]),
                np.array([maxlon, (maxlat + minlat) / 2])
                )
    else:
        height = np.linalg.norm(
            np.array([(maxlon + minlon) / 2, minlat]) -
            np.array([(maxlon + minlon) / 2, maxlat])) / 1000
        width = np.linalg.norm(
            np.array([minlon, (maxlat + minlat) / 2]) -
            np.array([maxlon, (maxlat + minlat) / 2])) / 1000

    return height, width


def make_regular_points(bounds, resolution, longlat=True):
    """
    Return a regular grid of points within `bounds` with the specified
    resolution.

    Parameters
    ----------
    bounds : 4-floats tuple
        The bbox of the grid, as xmin, ymin, xmax, ymax.
    resolution : int
        The resolution to use, in the same unit as `bounds`

    Returns
    -------
    points : numpy.array
        An array of coordinates
    shape : 2-floats tuple
        The number of points on each dimension (width, height)
    """
#    xmin, ymin, xmax, ymax = bounds
    minlon, minlat, maxlon, maxlat = bounds
    offset_lon = (maxlon - minlon) / 8
    offset_lat = (maxlat - minlat) / 8
    minlon -= offset_lon
    maxlon += offset_lon
    minlat -= offset_lat
    maxlat += offset_lat

    height, width = get_height_width((minlon, minlat, maxlon, maxlat), longlat)

    nb_x = int(round(width / resolution))
    nb_y = int(round(height / resolution))

    if nb_y * 0.6 > nb_x:
        nb_x = int(nb_x + nb_x / 3)
    elif nb_x * 0.6 > nb_y:
        nb_y = int(nb_y + nb_y / 3)

    if nb_y * nb_x > 200000:
        raise ProbableMemoryError(
            ("Please choose a lower resolution "
             "(by raising the value of the resolution parameter)"))
    return ((nb_x, nb_y),
            np.array([(x, y) for x in np.linspace(minlon, maxlon, nb_x)
                      for y in np.linspace(minlat, maxlat, nb_y)]))


class ProbableMemoryError(Exception):
    pass


def render_stewart(polygons, pot_layer, levels, nb_class, mask_layer):
    if mask_layer:
        try:
            renderer = _render_stewart_mask(
                polygons, pot_layer, levels, nb_class, mask_layer)
            return (False, renderer)
        except:
            renderer = _render_stewart(
                polygons, pot_layer, levels, nb_class)
            return (True, renderer)
    else:
        renderer = _render_stewart(
            polygons, pot_layer, levels, nb_class)
        return (False, renderer)


def _render_stewart(polygons, pot_layer, levels, nb_class):
    data_provider = pot_layer.dataProvider()
    colorRamp = QgsVectorGradientColorRampV2.create({
        'color1': '#ffffff',
        'color2': '#0037ff',
        'stops': '0.5;#72b2d7'})
    ranges = []
    features = []
    for i, poly in enumerate(polygons):
        if i == 0:
            last_level = 0
        else:
            last_level = float(levels[i - 1])
        current_level = float(levels[i])
        ft = QgsFeature()
        ft.setGeometry(poly)
        ft.setAttributes([i, last_level, current_level])
        features.append(ft)
        symbol = QgsFillSymbolV2()
        symbol.setColor(colorRamp.color(float(i) / len(polygons)))
        label = "{} - {}".format(last_level, current_level)
        rng = QgsRendererRangeV2(last_level, current_level, symbol, label)
        ranges.append(rng)
    data_provider.addFeatures(features[::-1])
    renderer = QgsGraduatedSymbolRendererV2('level_max', ranges)
    return renderer


def _render_stewart_mask(polygons, pot_layer, levels, nb_class, mask_layer):
    data_provider = pot_layer.dataProvider()
    colorRamp = QgsVectorGradientColorRampV2.create({
        'color1': '#ffffff',
        'color2': '#0037ff',
        'stops': '0.5;#72b2d7'})
    ranges = []
    features = []
    geoms = [loads(f.geometry().asWkb()) for f in mask_layer.getFeatures()]
    clip_geom = QgsGeometry.fromWkt(cascaded_union(geoms).wkt)

    for i, poly in enumerate(polygons):
        geom = poly.intersection(clip_geom.buffer(0, 16))
        if i == 0:
            last_level = 0
        else:
            last_level = float(levels[i - 1])
        current_level = float(levels[i])
        if geom.area() > 0:
            ft = QgsFeature()
            ft.setGeometry(geom)
            ft.setAttributes([i, last_level, current_level])
            features.append(ft)
            symbol = QgsFillSymbolV2()
            symbol.setColor(colorRamp.color(float(i) / len(polygons)))
            label = "{} - {}".format(last_level, current_level)
            rng = QgsRendererRangeV2(last_level, current_level, symbol, label)
            ranges.append(rng)
    data_provider.addFeatures(features[::-1])
    renderer = QgsGraduatedSymbolRendererV2('level_max', ranges)
    return renderer


def save_to_raster(pot, shape, bounds, proj4_value):
    minlon, minlat, maxlon, maxlat = bounds
    pixel_size_x = (maxlon - minlon) / shape[0]
    pixel_size_y = (maxlat - minlat) / shape[1]
    driver = gdal.GetDriverByName("GTiff")
    folder = mkdtemp()
    name = uuid4().hex
    path = os.path.sep.join([folder, name + '.geotiff'])
    dataset = driver.Create(path, shape[0], shape[1], 1, gdal.GDT_Float64)
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj4_value.encode("utf-8"))
    dataset.SetProjection(srs.ExportToWkt())
    dataset.SetGeoTransform(
        (minlon, pixel_size_x, 0, minlat, 0, pixel_size_y))
    dataset.GetRasterBand(1).WriteArray(pot.reshape((shape)).T)
    dataset = None
    return path


def color_raster(layer):
    provider = layer.dataProvider()
    extent = layer.extent()
    stats = provider.bandStatistics(1, QgsRasterBandStats.All, extent, 0)
    value_range = stats.maximumValue - stats.minimumValue

    value_list = [0] + [(value_range/i) for i in xrange(1, 5)][::-1]
    color_ramp_items = [
        QgsColorRampShader.ColorRampItem(value_list[0], QColor('#2c7bb6')),
        QgsColorRampShader.ColorRampItem(value_list[1], QColor('#abd9e9')),
        QgsColorRampShader.ColorRampItem(value_list[2], QColor('#ffffbf')),
        QgsColorRampShader.ColorRampItem(value_list[3], QColor('#fdae61')),
        QgsColorRampShader.ColorRampItem(value_list[4], QColor('#d7191c'))
        ]

    myRasterShader = QgsRasterShader()
    myColorRamp = QgsColorRampShader()
    myColorRamp.setColorRampItemList(color_ramp_items)
    myColorRamp.setColorRampType(QgsColorRampShader.INTERPOLATED)
    myRasterShader.setRasterShaderFunction(myColorRamp)
    myPseudoRenderer = QgsSingleBandPseudoColorRenderer(
        provider, layer.type(), myRasterShader)
    layer.setRenderer(myPseudoRenderer)
    layer.triggerRepaint()


def qgsgeom_from_mpl_collec(collections):
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

#def get_matdist_user(matdist, dim1, dim2):
#    try:
#        mat_dist = np.array([
#             map(int, feat.attributes()[1:])
#             for feat in matdist.getFeatures()
#             ])
#    except ValueError:
#        mat_dist = np.array([
#            map(float, feat.attributes()[1:])
#            for feat in matdist.getFeatures()
#            ])
#    assert dim1 in mat_dist.shape \
#        and dim2 in mat_dist.shape
#    return mat_dist
#
