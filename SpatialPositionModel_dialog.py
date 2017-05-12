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
from qgis.core import *
from .SpatialPositionModel_utils import (
    parse_expression, make_dist_mat, gen_unknownpts,
    compute_interact_density, compute_opportunity, compute_potentials,
    render_stewart, ProbableMemoryError, qgsgeom_from_mpl_collec)
from matplotlib.pyplot import contourf

import os.path


FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'SpatialPositionModel_dialog_base.ui'))


class SpatialPositionModelDialog(QtGui.QTabWidget, FORM_CLASS):
    def __init__(self, iface, parent=None):
        """Constructor."""
        super(SpatialPositionModelDialog, self).__init__(parent)
        self.setupUi(self)
        self.iface = iface
        self.host = None
        self.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.StewartComboBox_field.setLayer(x))
        self.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.StewartpushButton.setEnabled(True)
            if len(str(x)) > 0 else None)
        self.StewartComboBox_pts.layerChanged.connect(
            lambda x: self.mFieldExpressionWidget.setLayer(x))

        self.StewartpushButton.clicked.connect(self.run_stewart)
        self.buttonBox_close1.clicked.connect(self.close)
        self.StewartpushButton_clear.clicked.connect(self.clear_stewart_fields)
        self.pushButton_data.clicked.connect(self.load_dataset)

    def clean_name_factory(self, name):
        self.StewartComboBox_pts.setCurrentIndex(-1)
        self.StewartComboBox_mask.setCurrentIndex(-1)
        self.StewartComboBox_field.setCurrentIndex(-1)
        self.StewartComboBox_unknwPts.setCurrentIndex(-1)
        self.StewartComboBox_matdist.setCurrentIndex(-1)
        self.StewartcomboBox_function.setCurrentIndex(0)
        self.StewartdoubleSpinBox_beta.setValue(1.0)
        self.StewartdoubleSpinBox_span.setValue(1.0)
        self.StewartdoubleSpinBox_resolution.setValue(1.0)
        self.StewartpushButton.setEnabled(False)

    def clear_stewart_fields(self):
        self.clean_name_factory('Stewart')
        self.StewartspinBox_class.setValue(10)
        self.mFieldExpressionWidget.setLayer(None)

    def run_stewart(self):
        pts_layer = self.StewartComboBox_pts.currentLayer()
        mask_layer = self.StewartComboBox_mask.currentLayer()
#        unknownpts_layer = self.StewartComboBox_unknwPts.currentLayer()
        shape = None
#        matdist = self.StewartComboBox_matdist.currentLayer()
        function = (self.StewartcomboBox_function.currentText()).lower()

        self.crs = pts_layer.dataProvider().crs()
        self.longlat = self.crs.geographicFlag()

        if mask_layer:
            if self.crs.authid() != mask_layer.dataProvider().crs().authid():
                self.display_log_error("Crs mismatch : {} != {}".format(
                   self.crs.authid(),
                   mask_layer.dataProvider().crs().authid()
                   ), 5)
                mask_layer = None

        pts_coords = np.array(
            [f.geometry().asPoint() for f in pts_layer.getFeatures()])

        pts_values_field = self.StewartComboBox_field.currentField()
        other_values_fields = self.mFieldExpressionWidget.currentField()

        beta = self.StewartdoubleSpinBox_beta.value()
        span = self.StewartdoubleSpinBox_span.value()
        resolution = self.StewartdoubleSpinBox_resolution.value()
        nb_class = self.StewartspinBox_class.value()

        if (resolution == 0) or (span == 0) or (beta == 0):
            self.display_log_error(err, 4)
            return -1

        try:
            shape, unknownpts = \
                gen_unknownpts(
                    pts_layer, mask_layer, resolution, self.longlat)
        except ProbableMemoryError as err:
            self.display_log_error(err, 2)
            return -1
        except Exception as er:
            self.display_log_error(er, 3)
            return -1

        mat_dist = make_dist_mat(pts_coords, unknownpts, longlat=self.longlat)

        if not other_values_fields[0]:
            if pts_values_field:
                pts_values = np.array(
                    [f.attribute(pts_values_field)
                     for f in pts_layer.getFeatures()])
            else:
                pts_values = np.array([1 for i in xrange(len(pts_coords))])
            mat_dens = compute_interact_density(mat_dist, function, beta, span)
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

        x = np.array(list(set([c[0] for c in unknownpts])))
        y = np.array(list(set([c[1] for c in unknownpts])))
        xi = np.linspace(np.nanmin(x), np.nanmax(x), shape[0])
        yi = np.linspace(np.nanmin(y), np.nanmax(y), shape[1])

        collec_poly = contourf(
            xi, yi,
            pot.reshape((shape)).T,
            nb_class,
            vmax=abs(pot).max(),
            vmin=-abs(pot).max())
        levels = collec_poly.levels[1:]
        levels[-1] = np.nanmax(pot)
        levels = levels.tolist()

        res_poly = qgsgeom_from_mpl_collec(collec_poly.collections)

        pot_layer = QgsVectorLayer(
            "MultiPolygon?crs={}&field=id:integer"
            "&field=level_min:double"
            "&field=level_max:double".format(self.crs.authid()),
            "stewart_potentials_span_{}_beta_{}".format(span, beta), "memory")
        renderer = render_stewart(
            res_poly, pot_layer,
            levels, nb_class, mask_layer)

        pot_layer.setRendererV2(renderer)
        QgsMapLayerRegistry.instance().addMapLayer(pot_layer)
        self.iface.setActiveLayer(pot_layer)
        self.iface.zoomToActiveLayer()

    ## Todo : display better information/error message + show progression
    def display_log_error(self, error, msg_nb):
        error_msg = {
#            1: "Error when loading the choosen matrix. See QGis log for error traceback.",
            2: "The computation have been aborted, please choose a larger resolution value.",
            3: "Error when generating the grid of unknownpts.",
            4: "Span, resolution and beta should not be set to 0.",
            5: "Crs mismatch between point and mask layers. Mask not used."
            }
        QtGui.QMessageBox.information(
            self.iface.mainWindow(), 'Error', error_msg[msg_nb])
        QgsMessageLog.logMessage(
            'SpatialPositionModel plugin error report :\n {}'.format(error),
            level=QgsMessageLog.WARNING)

    def load_dataset(self):
        home_path = \
            os.getenv('HOMEPATH') or os.getenv('HOME') or os.getenv('USERPROFILE')
        self.iface.addVectorLayer(os.sep.join(
            [home_path, '.qgis2', 'python', 'plugins', 'SpatialPositionModel',
             'test_data', 'paris_mask.geojson']),
            'paris_mask', 'ogr')
        self.iface.addVectorLayer(os.sep.join(
            [home_path, '.qgis2', 'python', 'plugins', 'SpatialPositionModel',
             'test_data', 'paris_hospitals.geojson']),
            'paris_hospitals', 'ogr')
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
