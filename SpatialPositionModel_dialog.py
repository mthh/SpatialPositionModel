# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SpatialPositionModelDialog
                                 A QGIS plugin
 Compute spatial interaction models
                             -------------------
        begin                : 2015-12-10
        git sha              : $Format:%H$
        copyright            : (C) 2015 by #H
        email                : mth@#!.org
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os
import numpy as np
from PyQt4 import QtGui, uic
#from qgis.core import (
#    QgsVectorLayer, QgsFeature, QgsFillSymbolV2,
#    QgsVectorGradientColorRampV2, QgsGraduatedSymbolRendererV2,
#    QgsMapLayerRegistry)
from qgis.core import *
from .SpatialPositionModel_utils import *
import os.path
import processing

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'SpatialPositionModel_dialog_base.ui'))


class SpatialPositionModelDialog(QtGui.QTabWidget, FORM_CLASS):
    def __init__(self, iface, parent=None):
        """Constructor."""
        super(SpatialPositionModelDialog, self).__init__(parent)
        self.setupUi(self)
        self.iface = iface
        self.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.StewartComboBox_field.setLayer(x))
        self.HuffComboBox_pts.layerChanged.connect(
            lambda x: self.HuffComboBox_field.setLayer(x))
        self.ReillyComboBox_pts.layerChanged.connect(
            lambda x: self.ReillyComboBox_field.setLayer(x))
        self.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.StewartpushButton.setEnabled(True)
            if len(str(x)) > 0 else None)
        self.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.mFieldExpressionWidget.setLayer(x))
        self.ReillyComboBox_pts.layerChanged.connect(
            lambda x: self.ReillypushButton.setEnabled(True)
            if len(str(x)) > 0 else None)
        self.HuffComboBox_pts.layerChanged.connect(
            lambda x: self.HuffpushButton.setEnabled(True)
            if len(str(x)) > 0 else None)

        self.StewartpushButton.clicked.connect(self.run_stewart)
        self.ReillypushButton.clicked.connect(self.run_reilly)
        self.HuffpushButton.clicked.connect(self.run_huff)
        self.buttonBox_close1.clicked.connect(self.close)
        self.buttonBox_close2.clicked.connect(self.close)
        self.buttonBox_close3.clicked.connect(self.close)
        self.StewartpushButton_clear.clicked.connect(self.clear_stewart_fields)
        self.ReillypushButton_clear.clicked.connect(self.clear_reilly_fields)
        self.HuffpushButton_clear.clicked.connect(self.clear_huff_fields)
        self.pushButton_data.clicked.connect(self.load_dataset)

    def clean_name_factory(self, name):
        assert name in {'Stewart', 'Reilly', 'Huff'}
        base = '        self.'
        actions = [
            'ComboBox_pts.setCurrentIndex(-1)',
            'ComboBox_mask.setCurrentIndex(-1)',
            'ComboBox_field.setCurrentIndex(-1)',
            'ComboBox_unknwPts.setCurrentIndex(-1)',
            'ComboBox_matdist.setCurrentIndex(-1)',
            'comboBox_function.setCurrentIndex(0)',
            'doubleSpinBox_beta.setValue(1.0)',
            'doubleSpinBox_span.setValue(1.0)',
            'doubleSpinBox_resolution.setValue(1.0)',
            'pushButton.setEnabled(False)'
            ]
        for action in actions:
            eval('{}{}{}'.format(base, name, action))

    def clear_stewart_fields(self):
        self.clean_name_factory('Stewart')
        self.StewartspinBox_class.setValue(10)
        self.mFieldExpressionWidget.setLayer(None)

    def clear_reilly_fields(self):
        self.clean_name_factory('Reilly')

    def clear_huff_fields(self):
        self.clean_name_factory('Huff')

    def run_reilly(self):
        pts_layer = self.ReillyComboBox_pts.currentLayer()
        mask_layer = self.ReillyComboBox_mask.currentLayer()
        unknownpts_layer = self.ReillyComboBox_unknwPts.currentLayer()
        shape = None
        matdist = self.ReillyComboBox_matdist.currentLayer()
        function = (self.ReillycomboBox_function.currentText()).lower()

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])
        self.crs = int(pts_layer.dataProvider().crs().authid().split(':')[1])

        pts_values_field = self.ReillyComboBox_field.currentField()
        if pts_values_field:
            pts_values = np.array(
                [f.attribute(pts_values_field) for f in pts_layer.getFeatures()])
        else:
            pts_values = np.array([i for i in xrange(len(pts_coords))])

        beta = self.ReillydoubleSpinBox_beta.value()
        span = self.ReillydoubleSpinBox_span.value()
        resolution = self.ReillydoubleSpinBox_resolution.value()

        assert resolution != 0 and beta != 0 and span != 0

        if unknownpts_layer:
            unknownpts_coords = \
                [f.geometry().asPoint() for f in unknownpts_layer.getFeatures()]
            crs_unknpts = int(
                pts_layer.dataProvider().crs().authid().split(':')[1])
            assert self.crs == crs_unknpts
            unknownpts = np.array(unknownpts_coords)
        else:
            if mask_layer:
                crs_mask = int(
                    mask_layer.dataProvider().crs().authid().split(':')[1])
                assert self.crs == crs_mask
                ext = mask_layer.extent()
                bounds = (ext.xMinimum(), ext.yMinimum(),
                          ext.xMaximum(), ext.yMaximum())
            else:
                ext = pts_layer.extent()
                bounds = (ext.xMinimum(), ext.yMinimum(),
                          ext.xMaximum(), ext.yMaximum())
                tmp = ((bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
                tmp = tmp if tmp > span else span
                bounds = (bounds[0] - tmp, bounds[1] - tmp,
                          bounds[2] + tmp, bounds[3] + tmp)
            try:
                unknownpts, shape = make_regular_points(bounds, resolution)
            except ProbableMemoryError:
                QtGui.QMessageBox.information(
                    self.iface.mainWindow(), 'Error',
                    "The computation have been aborted, "
                    "please choose a larger resolution")
                return -1

        if not matdist:
            mat_dist = make_dist_mat(pts_coords, unknownpts, longlat=False)
        else:
            mat_dist = None

        mat_dens = compute_interact_density(mat_dist, function, beta, span)
        mat_opport = compute_opportunity(pts_values, mat_dens)
        reilly_val = compute_reilly(mat_opport)

        reilly_layer = QgsVectorLayer(
            "Point?crs=epsg:{}&field=id:integer"
            "&field=level:double".format(self.crs),
            "reilly_pts", "memory")
        data_provider = reilly_layer.dataProvider()
        features = []
        if not mask_layer:
            for i in xrange(len(unknownpts)):
                ft = QgsFeature()
                ft.setGeometry(QgsGeometry.fromPoint(
                    QgsPoint(unknownpts[i][0], unknownpts[i][1])))
                ft.setAttributes([i, float(reilly_val[i])])
                features.append(ft)
            data_provider.addFeatures(features)

        elif mask_layer:  # TODO : Améliorer le découpage qd il y a un mask
            index = QgsSpatialIndex()
            mask_features = {ft.id(): ft for ft in mask_layer.getFeatures()}
            map(index.insertFeature,mask_features.values())
            for i in xrange(len(unknownpts)):
                tmp_pt = QgsGeometry.fromPoint(
                    QgsPoint(unknownpts[i][0], unknownpts[i][1]))
                ids = index.intersects(tmp_pt.buffer(resolution, 4).boundingBox())
                for _id in ids:
                    mask_ft = mask_features[_id]
                    if tmp_pt.buffer(resolution, 4).intersects(mask_ft.geometry()):
                        ft = QgsFeature()
                        ft.setGeometry(tmp_pt)
                        ft.setAttributes([i, float(reilly_val[i])])
                        features.append(ft)
            data_provider.addFeatures(features)
        QgsMapLayerRegistry.instance().addMapLayer(reilly_layer)
        ext = reilly_layer.extent()
        offset = resolution / 2
        size = str(ext.xMinimum()-offset) + ',' + str(ext.xMaximum()+offset) \
            + ',' + str(ext.yMinimum()-offset) + ',' + str(ext.yMaximum()+offset)
        processing.runandload(
            'grass:v.to.rast.attribute', reilly_layer, 0,
            "level", size, resolution, -1, 0.0001, "Rasterized")
        QgsMapLayerRegistry.instance().removeMapLayer(reilly_layer.id())
        rast = QgsMapLayerRegistry.instance().mapLayersByName("Rasterized")[0]
        rast.setLayerName('Reilly_areas_span_{}_beta_{}'.format(span, beta))
        self.render_raster(rast, 'reilly')
        self.iface.setActiveLayer(rast)
        self.iface.zoomToActiveLayer()

    def render_raster(self, layer, mode):
        provider = layer.dataProvider()
        extent = layer.extent()
        stats = provider.bandStatistics(1, QgsRasterBandStats.All, extent, 0)
        value_range = stats.maximumValue - stats.minimumValue

        if 'huff' in mode:
            value_list = [0] + [(value_range/i) for i in xrange(1, 5)][::-1]
            color_ramp_items = [
                QgsColorRampShader.ColorRampItem(value_list[0], 
                                                 QtGui.QColor('#2c7bb6')), 
                QgsColorRampShader.ColorRampItem(value_list[1], 
                                                 QtGui.QColor('#abd9e9')), 
                QgsColorRampShader.ColorRampItem(value_list[2], 
                                                 QtGui.QColor('#ffffbf')),
                QgsColorRampShader.ColorRampItem(value_list[3], 
                                                 QtGui.QColor('#fdae61')),
                QgsColorRampShader.ColorRampItem(value_list[4], 
                                                 QtGui.QColor('#d7191c'))
                ]

        elif 'reilly' in mode:
            add = value_range // 2
            interval = stats.minimumValue + add
            valueList =[stats.minimumValue, interval, stats.maximumValue]
            color_ramp_items = [
                QgsColorRampShader.ColorRampItem(valueList[0], 
                                                 QtGui.QColor('#ff0000')), 
                QgsColorRampShader.ColorRampItem(valueList[1], 
                                                 QtGui.QColor('#ffff00')), 
                QgsColorRampShader.ColorRampItem(valueList[2], 
                                                 QtGui.QColor('#0000ff'))
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

    def run_huff(self):
        pts_layer = self.HuffComboBox_pts.currentLayer()
        mask_layer = self.HuffComboBox_mask.currentLayer()
        unknownpts_layer = self.HuffComboBox_unknwPts.currentLayer()
        shape = None
        matdist = self.HuffComboBox_matdist.currentLayer()
        function = (self.HuffcomboBox_function.currentText()).lower()

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])
        self.crs = int(pts_layer.dataProvider().crs().authid().split(':')[1])

        pts_values_field = self.HuffComboBox_field.currentField()
        if pts_values_field:
            pts_values = np.array(
                [f.attribute(pts_values_field) for f in pts_layer.getFeatures()])
        else:
            pts_values = np.array([1 for i in xrange(len(pts_coords))])

        beta = self.HuffdoubleSpinBox_beta.value()
        span = self.HuffdoubleSpinBox_span.value()
        resolution = self.HuffdoubleSpinBox_resolution.value()

        assert resolution != 0 and beta != 0 and span != 0

        if unknownpts_layer:
            unknownpts_coords = \
                [f.geometry().asPoint() for f in unknownpts_layer.getFeatures()]
            crs_unknpts = int(
                pts_layer.dataProvider().crs().authid().split(':')[1])
            assert self.crs == crs_unknpts
            unknownpts = np.array(unknownpts_coords)
        else:
            if mask_layer:
                crs_mask = int(
                    mask_layer.dataProvider().crs().authid().split(':')[1])
                assert self.crs == crs_mask
                ext = mask_layer.extent()
                bounds = (ext.xMinimum()-resolution, ext.yMinimum()-resolution,
                          ext.xMaximum()+resolution, ext.yMaximum()+resolution)

            else:
                ext = pts_layer.extent()
                bounds = (ext.xMinimum(), ext.yMinimum(),
                          ext.xMaximum(), ext.yMaximum())
                tmp = ((bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
                tmp = tmp if tmp > span else span
                bounds = (bounds[0] - tmp, bounds[1] - tmp,
                          bounds[2] + tmp, bounds[3] + tmp)
            try:
                unknownpts, shape = make_regular_points(bounds, resolution)
            except ProbableMemoryError:
                QtGui.QMessageBox.information(
                    self.iface.mainWindow(), 'Error',
                    "The computation have been aborted, "
                    "please choose a larger resolution")
                return -1

        if not matdist:
            mat_dist = make_dist_mat(pts_coords, unknownpts, longlat=False)
        else:
            mat_dist = None

        mat_dens = compute_interact_density(mat_dist, function, beta, span)
#        mat_opport = compute_opportunity(pts_values, mat_dens)
        huff_vals = compute_huff(compute_opportunity(pts_values, mat_dens))

        huff_layer = QgsVectorLayer(
            "Point?crs=epsg:{}&field=id:integer"
            "&field=level:double".format(self.crs),
            "huff_pts", "memory")
        data_provider = huff_layer.dataProvider()
        features = []
        if not mask_layer:
            for i in xrange(len(unknownpts)):
                ft = QgsFeature()
                ft.setGeometry(QgsGeometry.fromPoint(
                    QgsPoint(unknownpts[i][0], unknownpts[i][1])))
                ft.setAttributes([i, float(huff_vals[i])])
                features.append(ft)
            data_provider.addFeatures(features)

        elif mask_layer: # TODO : Améliorer le découpage qd il y a un mask
            index = QgsSpatialIndex()
            mask_features = {ft.id(): ft for ft in mask_layer.getFeatures()}
            map(index.insertFeature,mask_features.values())
            for i in xrange(len(unknownpts)):
                tmp_pt = QgsGeometry.fromPoint(
                    QgsPoint(unknownpts[i][0], unknownpts[i][1]))
                ids = index.intersects(tmp_pt.buffer(resolution, 4).boundingBox())
                for _id in ids:
                    mask_ft = mask_features[_id]
                    if tmp_pt.buffer(resolution, 4).intersects(mask_ft.geometry()):
                        ft = QgsFeature()
                        ft.setGeometry(tmp_pt)
                        ft.setAttributes([i, float(huff_vals[i])])
                        features.append(ft)
            data_provider.addFeatures(features)
        QgsMapLayerRegistry.instance().addMapLayer(huff_layer)
        ext = huff_layer.extent()
        offset = resolution / 2
        size = str(ext.xMinimum()-offset) + ',' + str(ext.xMaximum()+offset) \
            + ',' + str(ext.yMinimum()-offset) + ',' + str(ext.yMaximum()+offset)
        processing.runandload(
            'grass:v.to.rast.attribute', huff_layer, 0,
            "level", size, resolution, -1, 0.0001, "Rasterized")
        QgsMapLayerRegistry.instance().removeMapLayer(huff_layer.id())
        rast = QgsMapLayerRegistry.instance().mapLayersByName("Rasterized")[0]
        rast.setLayerName('Huff_potential_catchment_area_span_{}_beta_{}'
                          .format(span, beta))
        self.render_raster(rast, 'huff')
        self.iface.setActiveLayer(rast)
        self.iface.zoomToActiveLayer()

    def run_stewart(self):
        pts_layer = self.StewartComboBox_pts.currentLayer()
        mask_layer = self.StewartComboBox_mask.currentLayer()
        unknownpts_layer = self.StewartComboBox_unknwPts.currentLayer()
        shape = None
        matdist = self.StewartComboBox_matdist.currentLayer()
        function = (self.StewartcomboBox_function.currentText()).lower()

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])
        self.crs = int(pts_layer.dataProvider().crs().authid().split(':')[1])

        pts_values_field = self.StewartComboBox_field.currentField()
        other_values_fields = self.mFieldExpressionWidget.currentField()

        beta = self.StewartdoubleSpinBox_beta.value()
        span = self.StewartdoubleSpinBox_span.value()
        resolution = self.StewartdoubleSpinBox_resolution.value()
        nb_class = self.StewartspinBox_class.value()

        try:
            assert resolution != 0 and beta != 0 and span != 0
        except AssertionError:
            QtGui.QMessageBox.information(
                self.iface.mainWindow(), 'Error',
                "Span, resolution and beta should not be set to 0")
            return -1

        if unknownpts_layer:
            unknownpts_coords = [f.geometry().asPoint() for f
                                 in unknownpts_layer.getFeatures()]
            crs_unknpts = int(
                pts_layer.dataProvider().crs().authid().split(':')[1])
            assert self.crs == crs_unknpts
            unknownpts = np.array(unknownpts_coords)

        else:
            if mask_layer:
                crs_mask = int(
                    mask_layer.dataProvider().crs().authid().split(':')[1])
                assert self.crs == crs_mask
                ext = mask_layer.extent()
                bounds = (ext.xMinimum(), ext.yMinimum(),
                          ext.xMaximum(), ext.yMaximum())
            else:
                ext = pts_layer.extent()
                bounds = (ext.xMinimum(), ext.yMinimum(),
                          ext.xMaximum(), ext.yMaximum())

                tmp = ((bounds[2] - bounds[0]) / 10
                       + (bounds[3] - bounds[1]) / 10) / 2
                tmp = tmp if tmp > span else span

                bounds = (bounds[0] - tmp, bounds[1] - tmp,
                          bounds[2] + tmp, bounds[3] + tmp)

            try:
                unknownpts, shape = make_regular_points(bounds, resolution)
            except ProbableMemoryError:
                QtGui.QMessageBox.information(
                    self.iface.mainWindow(), 'Error',
                    "The computation have been aborted, "
                    "please choose a larger resolution value")
                return -1

        if not matdist:
            mat_dist = make_dist_mat(pts_coords, unknownpts, longlat=False)
        else:
            try:
                mat_dist = np.array([
                     map(int, feat.attributes()[1:])
                     for feat in matdist.getFeatures()
                     ])
            except ValueError:
                try:
                    mat_dist = np.array([
                        map(float, feat.attributes()[1:])
                        for feat in matdist.getFeatures()
                        ])
                except Exception as err:
                    print('Error when loading the matrice : {}'.format(err))

        assert len(unknownpts) in mat_dist.shape \
            and len(pts_coords) in mat_dist.shape

        if not other_values_fields[0]:
            if pts_values_field:
                pts_values = np.array(
                    [f.attribute(pts_values_field)
                     for f in pts_layer.getFeatures()])
            else:
                pts_values = np.array([1 for i in xrange(len(pts_coords))])
            mat_dens = compute_interact_density(mat_dist, function, beta, span)
#            mat_opport = compute_opportunity(pts_values, mat_dens)
            pot = compute_potentials(compute_opportunity(pts_values, mat_dens))

        else:
            fields, operator = parse_expression(other_values_fields[0])
            mat_dens = compute_interact_density(mat_dist, function, beta, span)
            pts_values1 = np.array(
                [f.attribute(fields[0]) for f in pts_layer.getFeatures()])
            pts_values2 = np.array(
                [f.attribute(fields[1]) for f in pts_layer.getFeatures()])
            mat_opport1 = compute_opportunity(pts_values1, mat_dens)
            mat_opport2 = compute_opportunity(pts_values2, mat_dens)
            pot1 = compute_potentials(mat_opport1)
            pot2 = compute_potentials(mat_opport2)
            pot = operator(pot1, pot2)
            # pts_values = [[np.array(
            #     [f.attribute(fields[i]) for f in pts_layer.getFeatures()])]
            #     for i in xrange(len(fields))]
            # pots = [
            #     compute_potentials(compute_opportunity(pts_values[i], mat_dens))
            #     for i in xrange(len(pts_values))
            #     ]
            # pot = 
        polygons, levels = render_stewart(pot, unknownpts, nb_class, shape)
#        print("Levels :\n{}".format(levels))
        pot_layer = QgsVectorLayer(
            "MultiPolygon?crs=epsg:{}&field=id:integer"
            "&field=level_min:double"
            "&field=level_max:double".format(self.crs),
            "stewart_potentials_span_{}_beta_{}".format(span, beta), "memory")
        data_provider = pot_layer.dataProvider()
        if mask_layer:
            features = []
            mask_geom = [f.geometry() for f in mask_layer.getFeatures()][0]
            for i, poly in enumerate(polygons):
                geom = poly.intersection(mask_geom.buffer(0, 16))
                if geom.area() > 0:
                    ft = QgsFeature()
                    ft.setGeometry(geom)
                    ft.setAttributes([i, float(levels[i]), float(levels[i+1])])
                    features.append(ft)
            data_provider.addFeatures(features[::-1])

        else:
            features = []
            for i, poly in enumerate(polygons):
                ft = QgsFeature()
                ft.setGeometry(poly)
                ft.setAttributes([i, float(levels[i]), float(levels[i+1])])
                features.append(ft)
            data_provider.addFeatures(features[::-1])

        symbol = QgsFillSymbolV2()
        colorRamp = QgsVectorGradientColorRampV2.create(
            {'color1' : '#ffffff',
             'color2' : '#0037ff',
             'stops' : '0.5;#72b2d7'}
             )
        renderer = QgsGraduatedSymbolRendererV2.createRenderer(
            pot_layer, 'id', nb_class, QgsGraduatedSymbolRendererV2.EqualInterval,
            symbol, colorRamp)
        pot_layer.setRendererV2(renderer)

        QgsMapLayerRegistry.instance().addMapLayer(pot_layer)
        self.iface.setActiveLayer(pot_layer)
        self.iface.zoomToActiveLayer()

    def load_dataset(self):
        home_path = \
            os.getenv('HOMEPATH') or os.getenv('HOME') or os.getenv('USERPROFILE')
        _  = self.iface.addVectorLayer(os.sep.join(
            [home_path, '.qgis2', 'python', 'plugins', 'SpatialPositionModel',
            'test_data', 'paris_mask.geojson']), 'paris_mask', 'ogr')
        _ = self.iface.addVectorLayer(os.sep.join(
            [home_path, '.qgis2', 'python', 'plugins', 'SpatialPositionModel',
            'test_data', 'paris_hospitals.geojson']), 'paris_hospitals', 'ogr')
        QtGui.QMessageBox.information(
            self.iface.mainWindow(), 'Dataset loaded',
            "Set these values and give a first try to the plug-in:\n\n"
            "Point layer = \"paris_hospitals\"\n"
            "Field = \"Capacite\"\n"
            "Mask layer = \"paris_mask\"\n"
            "resolution = 100\n"
            "beta = 3\n"
            "span = 1250\n"
            "\n")