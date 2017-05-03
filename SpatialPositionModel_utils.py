# -*- coding: utf-8 -*-
"""
SpatialPositionModel Utils
"""
import numpy as np
from qgis.core import (
    QgsGeometry, QgsPoint, QGis, QgsCoordinateTransform,
    QgsCoordinateReferenceSystem, QgsVectorGradientColorRampV2, QgsFeature,
    QgsGraduatedSymbolRendererV2, QgsFillSymbolV2, QgsVectorLayer,
    QgsSpatialIndex, QgsMapLayerRegistry,
    )
#from matplotlib.pyplot import contourf
from shapely.wkb import loads
from shapely.ops import cascaded_union

#import json
#from sys import version_info
#from httplib import HTTPConnection


#def check_host(url):
#    """ Helper function to get the hostname in desired format """
#    if len(url) < 4:
#        raise ValueError('Probably empty/non-valable url')
#    if not ('http' in url and '//' in url) and url[-1] == '/':
#        host = url[:-1]
#    elif not ('http:' in url and '//' in url):
#        host = url
#    elif 'http://' in url[:7] and url[-1] == '/':
#        host = url[7:-1]
#    elif 'http://' in url[:7]:
#        host = url[7:]
#    else:
#        host = url
#    return host


#def rectangular_light_table(src_coords, dest_coords, conn):
#    """
#    Function wrapping new OSRM 'table' function in order to get a rectangular
#    matrix of time distance as a numpy array
#
#    Params :
#
#        src_coords: list
#            A list of coord as (x, y) , like :
#                 list_coords = [(21.3224, 45.2358),
#                                (21.3856, 42.0094),
#                                (20.9574, 41.5286)] (coords have to be float)
#        dest_coords: list
#            A list of coord as (x, y) , like :
#                 list_coords = [(21.3224, 45.2358),
#                                (21.3856, 42.0094),
#                                (20.9574, 41.5286)] (coords have to be float)
#        conn: httplib.HTTPConnection object
#        headers: dict
#            headers (as dict) to be transmited to the
#            httplib.HTTPConnection.request function.
#
#    Output:
#        - a numpy array containing the time in tenth of seconds
#            (where 2147483647 means not-found route)
#        - a numpy array of snapped source coordinates (lat, lng)
#        - a numpy array of snapped destination coordinates (lat, lng)
#
#        ValueError is raised in case of any error
#            (wrong list of coords/ids, unknow host,
#            wrong response from the host, etc.)
#    """
#    headers = {
#        'connection': 'keep-alive',
#        'User-Agent': ' '.join(
#            ['QGIS-desktop', QGis.QGIS_VERSION, '/',
#             'Python-httplib', str(version_info[:3])[1:-1].replace(', ', '.')])
#        }
#    query = ['/table?src=']
#    # If only one source code (not nested) :
#    if len(src_coords) == 2 and not isinstance(src_coords[0],
#                                               (list, tuple, QgsPoint)):
#        query.append(''.join(
#            [str(src_coords[1]), ',', str(src_coords[0]), '&dst=']))
#    else: # Otherwise :
#        for coord in src_coords:  # Preparing the query
#            if coord is not None:
#                tmp = ''.join([str(coord[1]), ',', str(coord[0]), '&src='])
#                query.append(tmp)
#        query[-1] = query[-1][:-5] + '&dst='
#
#    if len(dest_coords) == 2 and not isinstance(dest_coords[0],
#                                                (list, tuple, QgsPoint)):
#        tmp = ''.join([str(dest_coords[1]), ',', str(dest_coords[0]), '&dst='])
#        query.append(tmp)
#    else:
#        for coord in dest_coords:  # Preparing the query
#            if coord is not None:
#                tmp = ''.join([str(coord[1]), ',', str(coord[0]), '&dst='])
#                query.append(tmp)
#
#    query = (''.join(query))[:-5]
#    try:  # Querying the OSRM instance
#        conn.request('GET', query, headers=headers)
#        parsed_json = json.loads(conn.getresponse().read().decode('utf-8'))
#    except Exception as err:
#        raise ValueError('Error while contacting OSRM instance : \n{}'
#                         .format(err))
#
#    if 'distance_table' in parsed_json.keys():  # Preparing the result matrix
#        mat = np.array(parsed_json['distance_table'], dtype='int32')
#        src_snapped = \
#            np.array(parsed_json['source_coordinates'], dtype=float)
#        dest_snapped = \
#            np.array(parsed_json['destination_coordinates'], dtype=float)
#        return mat, src_snapped, dest_snapped
#    else:
#        raise ValueError('No distance table return by OSRM instance')
#

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
                if expr[i+1] != ' ':
                    nexpr.append(' ')
            except IndexError:
                pass
    nexpr = ''.join(nexpr).split(' ')
    for i in nexpr:
        if i.isdigit() and not i in ('*', '+', '-', '/'):
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

def prepare_raster(mode, crs, values, unknownpts, mask_layer, resolution):
    pts_layer = QgsVectorLayer(
        "Point?crs={}&field=id:integer&field=level:double".format(crs),
        "{}_pts".format(mode), "memory")
    data_provider = pts_layer.dataProvider()
    features = []
    if not mask_layer:
        for i in xrange(len(unknownpts)):
            ft = QgsFeature()
            ft.setGeometry(QgsGeometry.fromPoint(
                QgsPoint(unknownpts[i][0], unknownpts[i][1])))
            ft.setAttributes([i, float(values[i])])
            features.append(ft)
        data_provider.addFeatures(features)

    elif mask_layer: # TODO : Améliorer le découpage qd il y a un mask
        index = QgsSpatialIndex()
        mask_features = {ft.id(): ft for ft in mask_layer.getFeatures()}
        map(index.insertFeature, mask_features.values())
        for i in xrange(len(unknownpts)):
            tmp_pt = QgsGeometry.fromPoint(
                QgsPoint(unknownpts[i][0], unknownpts[i][1]))
            ids = index.intersects(tmp_pt.buffer(resolution, 4).boundingBox())
            for _id in ids:
                mask_ft = mask_features[_id]
                if tmp_pt.buffer(resolution, 4).intersects(mask_ft.geometry()):
                    ft = QgsFeature()
                    ft.setGeometry(tmp_pt)
                    ft.setAttributes([i, float(values[i])])
                    features.append(ft)
        data_provider.addFeatures(features)
    QgsMapLayerRegistry.instance().addMapLayer(pts_layer)
    ext = pts_layer.extent()
    offset = resolution / 2
    grass_region_size = str(ext.xMinimum()-offset) + ',' + str(ext.xMaximum()+offset) \
        + ',' + str(ext.yMinimum()-offset) + ',' + str(ext.yMaximum()+offset)
    return pts_layer, grass_region_size


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
        d0 = np.subtract.outer(xy1[:,0], xy2[:,0])
        d1 = np.subtract.outer(xy1[:,1], xy2[:,1])
        return np.hypot(d0, d1)
    elif longlat:
        return hav_dist(xy1[:, None], xy2)


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


def gen_unknownpts(pts_layer, mask_layer, resolution, longlat):
    if mask_layer:
        crs_mask = int(
            mask_layer.dataProvider().crs().authid().split(':')[1])
        assert int(pts_layer.dataProvider().crs().authid().split(':')[1]) \
            == crs_mask

        ext = mask_layer.extent()
        bounds = (ext.xMinimum(), ext.yMinimum(),
                  ext.xMaximum(), ext.yMaximum())
    else:
        ext = pts_layer.extent()
        bounds = (ext.xMinimum(), ext.yMinimum(),
                  ext.xMaximum(), ext.yMaximum())

        tmp = (
            (bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
#                tmp = tmp if tmp > span else span

        bounds = (bounds[0] - tmp, bounds[1] - tmp,
                  bounds[2] + tmp, bounds[3] + tmp)

    return make_regular_points(bounds, resolution, longlat)


#def get_osrm_matdist(host, pts_coords, unknownpts, crs):
#    conn = HTTPConnection(host)
#    xform = QgsCoordinateTransform(QgsCoordinateReferenceSystem(crs),
#                                   QgsCoordinateReferenceSystem(4326))
#    pts_coords_temp = \
#        [xform.transform(QgsPoint(*point)) for point in pts_coords]
#    unknownpts_temp = \
#        [xform.transform(QgsPoint(*point)) for point in unknownpts]
#    mat_dist, _, _ = \
#        rectangular_light_table(pts_coords_temp, unknownpts_temp, conn)
#
#    return mat_dist


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

def compute_opportunity(pts_values, matdens):
    matOpport = pts_values[:, np.newaxis] * matdens
    return matOpport.round(8)


def compute_potentials(matopport):
    return matopport.sum(axis=0)


def compute_reilly(mattoportes):
    return np.argmax(mattoportes, axis=0)


def compute_huff(matopport):
    sum_lines = matopport.sum(1)
    sum_lines = sum_lines[np.where(sum_lines > 0)[0]]
    matopportPct = np.array([100]) \
        * (matopport[np.where(sum_lines > 0)[0]] / sum_lines[:, np.newaxis]).T
    return matopportPct.max(axis=1)


#def make_regular_points(bounds, reso, skip_limit=False, longlat=True):
#    """
#    Return a grid of regular points.
#    """
#    xmin, ymin, xmax, ymax = bounds
#    nb_x = int(round((xmax - xmin) / reso + ((xmax - xmin) / reso) /10))
#    nb_y = int(round((ymax - ymin) / reso + ((ymax - ymin) / reso) /10))
#    try:
#        prog_x = \
#            [(xmin - (xmax - xmin) / 20) + reso * i for i in range(nb_x + 1)]
#        prog_y = \
#            [(ymin - (ymax - ymin) / 20) + reso * i for i in range(nb_y + 1)]
#    except ZeroDivisionError:
#        raise ZeroDivisionError(
#            'Please choose a finest resolution (by lowering the value of the '
#            'resolution argument and/or providing an appropriate mask layer')
#    print(len(prog_x) * len(prog_y))
#    if not skip_limit:            
#        if len(prog_x) * len(prog_y) > 500000:
#            raise ProbableMemoryError("Please choose a lower resolution"
#                                      " (by raising the value of the resolutio"
#                                      "n parameter)")
#    return (np.array([(x, y) for x in prog_x for y in prog_y]),
#            (len(prog_x), len(prog_y)))
#

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
            np.array([(maxlon + minlon) / 2, minlat])
            - np.array([(maxlon + minlon) / 2, maxlat]))
        width = np.linalg.norm(
            np.array([minlon, (maxlat + minlat) / 2])
            - np.array([maxlon, (maxlat + minlat) / 2]))

    nb_x = int(round(width / resolution))
    nb_y = int(round(height / resolution))
    if nb_y * 0.6 > nb_x:
        nb_x = int(nb_x + nb_x / 3)
    elif nb_x * 0.6 > nb_y:
        nb_y = int(nb_y + nb_y / 3)
        if nb_y * nb_x > 60000:
            raise ProbableMemoryError(
                "Please choose a lower resolution (by raising the value of the resolution parameter)")
    return (
        np.array([(x, y) for x in np.linspace(minlon, maxlon, nb_x) for y in np.linspace(minlat, maxlat, nb_y)]),
        (nb_x, nb_y)
        )


class ProbableMemoryError(Exception):
    pass


def render_stewart(polygons, pot_layer, levels, nb_class, mask_layer):
    data_provider = pot_layer.dataProvider()
    if mask_layer:
        features = []
        geoms = [loads(f.geometry().asWkb()) for f in mask_layer.getFeatures()]
        clip_geom = QgsGeometry.fromWkt(cascaded_union(geoms).wkt)

        for i, poly in enumerate(polygons):
            geom = poly.intersection(clip_geom.buffer(0, 16))
            if i == 0:
                last_level = 0
            else:
                last_level = float(levels[i-1])
            if geom.area() > 0:
                ft = QgsFeature()
                ft.setGeometry(geom)
                ft.setAttributes([i, last_level, float(levels[i])])
                features.append(ft)
        data_provider.addFeatures(features[::-1])

    else:
        features = []
        for i, poly in enumerate(polygons):
            if i == 0:
                last_level = 0
            else:
                last_level = float(levels[i-1])
            ft = QgsFeature()
            ft.setGeometry(poly)
            ft.setAttributes([i, last_level, float(levels[i])])
            features.append(ft)
        data_provider.addFeatures(features[::-1])

    symbol = QgsFillSymbolV2()
    colorRamp = QgsVectorGradientColorRampV2.create(
        {'color1': '#ffffff',
         'color2': '#0037ff',
         'stops': '0.5;#72b2d7'})

    renderer = QgsGraduatedSymbolRendererV2.createRenderer(
        pot_layer, 'id', nb_class,
        QgsGraduatedSymbolRendererV2.EqualInterval,
        symbol, colorRamp)

    return renderer


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
