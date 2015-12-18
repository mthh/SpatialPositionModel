# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SpatialPositionModel
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
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt4.QtGui import QAction, QIcon
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from qgis.core import (
    QgsVectorLayer, QgsFeature, QgsFillSymbolV2,
    QgsVectorGradientColorRampV2, QgsGraduatedSymbolRendererV2,
    QgsMapLayerRegistry)
from SpatialPositionModel_dialog import SpatialPositionModelDialog
from .SpatialPositionModel_utils import *
import os.path
import numpy as np

class SpatialPositionModel:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'SpatialPositionModel_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = SpatialPositionModelDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&SpatialPositionModel')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'SpatialPositionModel')
        self.toolbar.setObjectName(u'SpatialPositionModel')

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('SpatialPositionModel', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/SpatialPositionModel/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'SpatialPositionModel'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&SpatialPositionModel'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def name_factory(self, name):
        assert name in {'Stewart', 'Reilly', 'Huff'}
        base = '        self.dlg.'
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
        self.name_factory('Stewart')
        self.dlg.StewartspinBox_class.setValue(10)
        self.dlg.mFieldExpressionWidget.setLayer(None)

    def clear_reilly_fields(self):
        self.name_factory('Reilly')

    def clear_huff_fields(self):
        self.name_factory('Huff')

    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.StewartComboBox_field.setLayer(x))
        self.dlg.HuffComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.HuffComboBox_field.setLayer(x))
        self.dlg.ReillyComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.ReillyComboBox_field.setLayer(x))
        self.dlg.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.StewartpushButton.setEnabled(True)
            if len(str(x)) > 0 else None
            )
        self.dlg.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.mFieldExpressionWidget.setLayer(x))
        self.dlg.ReillyComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.ReillypushButton.setEnabled(True)
            if len(str(x)) > 0 else None
            )
        self.dlg.HuffComboBox_pts.layerChanged.connect(
            lambda x: self.dlg.HuffpushButton.setEnabled(True)
            if len(str(x)) > 0 else None
            )
        
        self.dlg.StewartpushButton.clicked.connect(self.run_stewart)
        self.dlg.ReillypushButton.clicked.connect(self.run_reilly)
        self.dlg.HuffpushButton.clicked.connect(self.run_huff)
        self.dlg.buttonBox_close1.clicked.connect(self.dlg.close)
        self.dlg.buttonBox_close2.clicked.connect(self.dlg.close)
        self.dlg.buttonBox_close3.clicked.connect(self.dlg.close)
        self.dlg.StewartpushButton_clear.clicked.connect(self.clear_stewart_fields)
        self.dlg.ReillypushButton_clear.clicked.connect(self.clear_reilly_fields)
        self.dlg.HuffpushButton_clear.clicked.connect(self.clear_huff_fields)
        self.dlg.show()

    def run_reilly(self):
        pts_layer = self.dlg.ReillyComboBox_pts.currentLayer()
        mask_layer = self.dlg.ReillyComboBox_mask.currentLayer()
        unknownpts_layer = self.dlg.ReillyComboBox_unknwPts.currentLayer()
        matdist = self.dlg.ReillyComboBox_matdist.currentLayer()
        function = (self.dlg.ReillycomboBox_function.currentText()).lower()

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])
        self.crs = int(pts_layer.dataProvider().crs().authid().split(':')[1])


        pts_values_field = self.dlg.ReillyComboBox_field.currentField()
        if pts_values_field:
            pts_values = np.array(
                [f.attribute(pts_values_field) for f in pts_layer.getFeatures()])
        else:
            pts_values = np.array([i for i in xrange(len(pts_coords))])

        beta = self.dlg.ReillydoubleSpinBox_beta.value()
        span = self.dlg.ReillydoubleSpinBox_span.value()
        resolution = self.dlg.ReillydoubleSpinBox_resolution.value()

        assert resolution != 0
        assert beta != 0
        assert span != 0

        if unknownpts_layer:
            pts_coords = [f.geometry().asPoint() for f in unknownpts_layer.getFeatures()]
            crs_unknpts = int(pts_layer.dataProvider().crs().authid().split(':')[1])
            assert self.crs == crs_unknpts
        else:
            if mask_layer:
                crs_mask = int(mask_layer.dataProvider().crs().authid().split(':')[1])
                assert self.crs == crs_mask
                bounds = (mask_layer.extent().xMinimum(),
                          mask_layer.extent().yMinimum(),
                          mask_layer.extent().xMaximum(),
                          mask_layer.extent().yMaximum())
            else:
                bounds = (pts_layer.extent().xMinimum(),
                          pts_layer.extent().yMinimum(),
                          pts_layer.extent().xMaximum(),
                          pts_layer.extent().yMaximum())
                tmp = ((bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
                tmp = tmp if tmp > span else span
                bounds = (bounds[0] - tmp, bounds[1] - tmp,
                          bounds[2] + tmp, bounds[3] + tmp)
                
            unknownpts = make_regular_points(bounds, resolution)

        if not matdist:
            mat_dist = make_dist_mat(pts_coords, unknownpts, longlat=False)
        else:
            mat_dist = None

        mat_dens = compute_interact_density(mat_dist, function, beta, span)
        mat_opport = compute_opportunity(pts_values, mat_dens)
        reilly_val = compute_reilly(mat_opport)

#        polygons, levels = render_reilly(reilly_val, unknownpts)

        reilly_layer = QgsVectorLayer(
            "Point?crs=epsg:{}&field=id:integer"
            "&field=level:integer(10)".format(self.crs),
            "reilly_pts", "memory")
        data_provider = reilly_layer.dataProvider()
#        if mask_layer:
#            features = []
#            mask_geom = [f.geometry() for f in mask_layer.getFeatures()][0]
#            for i, poly in enumerate(polygons):
#                ft = QgsFeature()
#                ft.setGeometry(poly.intersection(mask_geom.buffer(0, 16)))
#                ft.setAttributes([i, levels[i]])
#                features.append(ft)
#            data_provider.addFeatures(features[::-1])
#        else:
#            features = []
#            for i, poly in enumerate(polygons):
#                ft = QgsFeature()
#                ft.setGeometry(poly)
#                ft.setAttributes([i, levels[i]])
#                features.append(ft)
#            data_provider.addFeatures(features[::-1])
#
#        symbol =  QgsFillSymbolV2()
#        colorRamp = QgsVectorGradientColorRampV2.create(
#            {'color1' : '#ffffff',
#             'color2' : '#0037ff',
#             'stops' : '0.5;#72b2d7'})
#        renderer = QgsGraduatedSymbolRendererV2.createRenderer(
#            reilly_layer, 'id', 12, QgsGraduatedSymbolRendererV2.EqualInterval,
#            symbol, colorRamp)
#        reilly_layer.setRendererV2(renderer)
#
#        QgsMapLayerRegistry.instance().addMapLayer(reilly_layer)
#        self.iface.setActiveLayer(reilly_layer)
#        self.iface.zoomToActiveLayer()

    def run_huff(self):
        pts_layer = self.dlg.HuffComboBox_pts.currentLayer()
        mask_layer = self.dlg.HuffComboBox_mask.currentLayer()
        unknownpts_layer = self.dlg.HuffComboBox_unknwPts.currentLayer()
        matdist = self.dlg.HuffComboBox_matdist.currentLayer()
        function = (self.dlg.HuffcomboBox_function.currentText()).lower()

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])
        self.crs = int(pts_layer.dataProvider().crs().authid().split(':')[1])


        pts_values_field = self.dlg.ReillyComboBox_field.currentField()
        if pts_values_field:
            pts_values = np.array(
                [f.attribute(pts_values_field) for f in pts_layer.getFeatures()])
        else:
            pts_values = np.array([i for i in xrange(len(pts_coords))])

        beta = self.dlg.ReillydoubleSpinBox_beta.value()
        span = self.dlg.ReillydoubleSpinBox_span.value()
        resolution = self.dlg.ReillydoubleSpinBox_resolution.value()

        assert resolution != 0
        assert beta != 0
        assert span != 0

        if unknownpts_layer:
            unknwownpts_coords = \
                [f.geometry().asPoint() for f in unknownpts_layer.getFeatures()]
            crs_unknpts = int(pts_layer.dataProvider().crs().authid().split(':')[1])
            assert self.crs == crs_unknpts
        else:
            if mask_layer:
                crs_mask = int(mask_layer.dataProvider().crs().authid().split(':')[1])
                assert self.crs == crs_mask
                bounds = (mask_layer.extent().xMinimum(),
                          mask_layer.extent().yMinimum(),
                          mask_layer.extent().xMaximum(),
                          mask_layer.extent().yMaximum())
            else:
                bounds = (pts_layer.extent().xMinimum(),
                          pts_layer.extent().yMinimum(),
                          pts_layer.extent().xMaximum(),
                          pts_layer.extent().yMaximum())
                tmp = ((bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
                tmp = tmp if tmp > span else span
                bounds = (bounds[0] - tmp, bounds[1] - tmp,
                          bounds[2] + tmp, bounds[3] + tmp)
                
            unknwownpts_coords = make_regular_points(bounds, resolution)

        if not matdist:
            mat_dist = make_dist_mat(pts_coords, unknwownpts_coords, longlat=False)
        else:
            mat_dist = None

        mat_dens = compute_interact_density(mat_dist, function, beta, span)
        mat_opport = compute_opportunity(pts_values, mat_dens)
        huff_vals = compute_huff(mat_opport)

        polygons, levels = render_stewart(huff_vals, unknownpts)

        pot_layer = QgsVectorLayer(
            "MultiPolygon?crs=epsg:{}&field=id:integer"
            "&field=level:integer(10)".format(self.crs),
            "huff_catchments_areas", "memory")
        data_provider = pot_layer.dataProvider()
        if mask_layer:
            features = []
            mask_geom = [f.geometry() for f in mask_layer.getFeatures()][0]
            for i, poly in enumerate(polygons):
                ft = QgsFeature()
                ft.setGeometry(poly.intersection(mask_geom.buffer(0, 16)))
                ft.setAttributes([i, levels[i]])
                features.append(ft)
            data_provider.addFeatures(features[::-1])
        else:
            features = []
            for i, poly in enumerate(polygons):
                ft = QgsFeature()
                ft.setGeometry(poly)
                ft.setAttributes([i, levels[i]])
                features.append(ft)
            data_provider.addFeatures(features[::-1])

        symbol =  QgsFillSymbolV2()
        colorRamp = QgsVectorGradientColorRampV2.create(
            {'color1' : '#ffffff',
             'color2' : '#0037ff',
             'stops' : '0.5;#72b2d7'})
        renderer = QgsGraduatedSymbolRendererV2.createRenderer(
            pot_layer, 'id', 12, QgsGraduatedSymbolRendererV2.EqualInterval,
            symbol, colorRamp)
        pot_layer.setRendererV2(renderer)

        QgsMapLayerRegistry.instance().addMapLayer(pot_layer)
        self.iface.setActiveLayer(pot_layer)
        self.iface.zoomToActiveLayer()

    def run_stewart(self):
        pts_layer = self.dlg.StewartComboBox_pts.currentLayer()
        mask_layer = self.dlg.StewartComboBox_mask.currentLayer()
        unknownpts_layer = self.dlg.StewartComboBox_unknwPts.currentLayer()
        matdist = self.dlg.StewartComboBox_matdist.currentLayer()
        function = (self.dlg.StewartcomboBox_function.currentText()).lower()

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])
        self.crs = int(pts_layer.dataProvider().crs().authid().split(':')[1])

        pts_values_field = self.dlg.StewartComboBox_field.currentField()
        other_values_fields = self.dlg.mFieldExpressionWidget.currentField()
        
        beta = self.dlg.StewartdoubleSpinBox_beta.value()
        span = self.dlg.StewartdoubleSpinBox_span.value()
        resolution = self.dlg.StewartdoubleSpinBox_resolution.value()
        nb_class = self.dlg.StewartspinBox_class.value()

        assert resolution != 0
        assert beta != 0
        assert span != 0

        if unknownpts_layer:
            unknownpts_coords = [f.geometry().asPoint() for f in unknownpts_layer.getFeatures()]
            crs_unknpts = int(pts_layer.dataProvider().crs().authid().split(':')[1])
            assert self.crs == crs_unknpts
            unknownpts = np.array(unknownpts_coords)

        else:
            if mask_layer:
                crs_mask = int(mask_layer.dataProvider().crs().authid().split(':')[1])
                assert self.crs == crs_mask
                bounds = (mask_layer.extent().xMinimum(),
                          mask_layer.extent().yMinimum(),
                          mask_layer.extent().xMaximum(),
                          mask_layer.extent().yMaximum())
            else:
                bounds = (pts_layer.extent().xMinimum(),
                          pts_layer.extent().yMinimum(),
                          pts_layer.extent().xMaximum(),
                          pts_layer.extent().yMaximum())
                tmp = ((bounds[2] - bounds[0]) / 10 + (bounds[3] - bounds[1]) / 10) / 2
                tmp = tmp if tmp > span else span
                bounds = (bounds[0] - tmp, bounds[1] - tmp,
                          bounds[2] + tmp, bounds[3] + tmp)
                
            unknownpts = make_regular_points(bounds, resolution)

        if len(unknownpts) > 2000000:
            QMessageBox.information(
                self.iface.mainWindow(), 'Error',
                "The computation have been aborted, please choose a larger resolution")
            return -1

        if not matdist:
            mat_dist = make_dist_mat(pts_coords, unknownpts, longlat=False)
        else:
            try:
                 mat_dist = np.array([
                    map(int,feat.attributes()[1:]) 
                    for feat in matdist.getFeatures()
                    ])
            except ValueError:
                try:
                    mat_dist = np.array([
                        map(float,feat.attributes()[1:]) 
                        for feat in matdist.getFeatures()
                        ])
                except Exception as err:
                    print('Error when loading the matrice : {}'.format(err))

        assert len(unknownpts) in mat_dist.shape
        assert len(pts_coords) in mat_dist.shape

        if not other_values_fields[0]:
            if pts_values_field:
                pts_values = np.array(
                    [f.attribute(pts_values_field) for f in pts_layer.getFeatures()])
            else:
                pts_values = np.array([i for i in xrange(len(pts_coords))])
            mat_dens = compute_interact_density(mat_dist, function, beta, span)
            mat_opport = compute_opportunity(pts_values, mat_dens)
            pot = compute_potentials(mat_opport)
        else:
            fields, operator = parse_expression(other_values_fields[0])
            pts_values1 = np.array(
                    [f.attribute(fields[0]) for f in pts_layer.getFeatures()])
            pts_values2 = np.array(
                    [f.attribute(fields[1]) for f in pts_layer.getFeatures()])
            mat_dens = compute_interact_density(mat_dist, function, beta, span)
            mat_opport1 = compute_opportunity(pts_values1, mat_dens)
            mat_opport2 = compute_opportunity(pts_values2, mat_dens)
            pot1 = compute_potentials(mat_opport1)
            pot2 = compute_potentials(mat_opport2)
            pot = operator(pot1, pot2)

        polygons, levels = render_stewart(pot, unknownpts, nb_class)

        pot_layer = QgsVectorLayer(
            "MultiPolygon?crs=epsg:{}&field=id:integer"
            "&field=level:integer(10)".format(self.crs),
            "stewart_potentials_span_{}_beta_{}".format(span, beta), "memory")
        data_provider = pot_layer.dataProvider()
        if mask_layer:
            features = []
            mask_geom = [f.geometry() for f in mask_layer.getFeatures()][0]
            for i, poly in enumerate(polygons):
                ft = QgsFeature()
                ft.setGeometry(poly.intersection(mask_geom.buffer(0, 16)))
                ft.setAttributes([i, levels[i]])
                features.append(ft)
            data_provider.addFeatures(features[::-1])
        else:
            features = []
            for i, poly in enumerate(polygons):
                ft = QgsFeature()
                ft.setGeometry(poly)
                ft.setAttributes([i, levels[i]])
                features.append(ft)
            data_provider.addFeatures(features[::-1])

        symbol =  QgsFillSymbolV2()
        colorRamp = QgsVectorGradientColorRampV2.create(
            {'color1' : '#ffffff',
             'color2' : '#0037ff',
             'stops' : '0.5;#72b2d7'})
        renderer = QgsGraduatedSymbolRendererV2.createRenderer(
            pot_layer, 'id', 12, QgsGraduatedSymbolRendererV2.EqualInterval,
            symbol, colorRamp)
        pot_layer.setRendererV2(renderer)

        QgsMapLayerRegistry.instance().addMapLayer(pot_layer)
        self.iface.setActiveLayer(pot_layer)
        self.iface.zoomToActiveLayer()
