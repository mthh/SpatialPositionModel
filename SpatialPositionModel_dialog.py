# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SpatialPositionModelDialog
                                 A QGIS plugin
 Compute spatial interaction models
                              -------------------
        begin                : 2017-05-15
        git sha              : $Format:%H$
        copyright            : (C) 2015 by mthh
        email                : matthieu.viry@cnrs.fr
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
    render_stewart, ProbableMemoryError, qgsgeom_from_mpl_collec,
    parse_class_breaks, get_height_width, save_to_raster, color_raster
    )
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

        self.StewartComboBox_pts.layerChanged.connect(self.on_change_layer)
        self.StewartpushButton.clicked.connect(self.run_stewart)
        self.buttonBox_close1.clicked.connect(self.close)
        self.StewartpushButton_clear.clicked.connect(self.clear_stewart_fields)
        self.pushButton_data.clicked.connect(self.load_dataset)
        self.radioButton_vector.clicked.connect(
            lambda _: self.toggleVectorRaster('vector'))
        self.radioButton_raster.clicked.connect(
            lambda _: self.toggleVectorRaster('raster'))

    def toggleVectorRaster(self, value):
        if 'vector' in value:
            self.StewarttextEdit_breaks.setEnabled(True)
            self.StewartspinBox_class.setEnabled(True)
            self.StewartComboBox_mask.setEnabled(True)
            self.StewarttextEdit_breaks.setVisible(True)
            self.StewartspinBox_class.setVisible(True)
            self.StewartComboBox_mask.setVisible(True)
            self.label_25.setVisible(True)
            self.label_38.setVisible(True)
            self.label_3.setVisible(True)
        else:
            self.StewarttextEdit_breaks.setEnabled(False)
            self.StewartspinBox_class.setEnabled(False)
            self.StewartComboBox_mask.setEnabled(False)
            self.StewarttextEdit_breaks.setVisible(False)
            self.StewartspinBox_class.setVisible(False)
            self.StewartComboBox_mask.setVisible(False)
            self.label_25.setVisible(False)
            self.label_38.setVisible(False)
            self.label_3.setVisible(False)

    def on_change_layer(self, layer):
        self.StewartComboBox_field.setLayer(layer)
        self.mFieldExpressionWidget.setLayer(layer)
        if not layer or not layer.extent() \
                or not layer.dataProvider() \
                or not layer.dataProvider().crs():
                return
        try:
            self.StewartpushButton.setEnabled(True)
            ext = layer.extent()
            bounds = (ext.xMinimum(), ext.yMinimum(),
                      ext.xMaximum(), ext.yMaximum())
            height, width = get_height_width(
                bounds, layer.dataProvider().crs().geographicFlag())
            reso = max([(height / 90), (width / 90)])
            reso += reso * 0.2
            self.StewartdoubleSpinBox_resolution.setValue(round(reso))
            self.StewartdoubleSpinBox_span.setValue(round(reso * 2.5))
        except TypeError:
            pass

    def clear_stewart_fields(self):
        self.StewartComboBox_pts.setCurrentIndex(-1)
        self.StewartComboBox_mask.setCurrentIndex(-1)
        self.StewartComboBox_field.setCurrentIndex(-1)
        self.StewartComboBox_unknwPts.setCurrentIndex(-1)
        self.StewartComboBox_matdist.setCurrentIndex(-1)
        self.StewartcomboBox_function.setCurrentIndex(0)
        self.StewartdoubleSpinBox_beta.setValue(1.0)
        self.StewartdoubleSpinBox_span.setValue(1.0)
        self.StewartdoubleSpinBox_resolution.setValue(1.0)
        self.StewarttextEdit_breaks.setPlainText("")
        self.StewartpushButton.setEnabled(False)
        self.StewartspinBox_class.setValue(7)
        self.mFieldExpressionWidget.setLayer(None)
        self.radioButton_vector.setChecked(True)

    def run_stewart(self):
        pts_layer = self.StewartComboBox_pts.currentLayer()
        mask_layer = self.StewartComboBox_mask.currentLayer() \
            if self.radioButton_vector.isChecked() else None
        shape = None
        function = (self.StewartcomboBox_function.currentText()).lower()
#        unknownpts_layer = self.StewartComboBox_unknwPts.currentLayer()
#        matdist = self.StewartComboBox_matdist.currentLayer()

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

        if self.radioButton_vector.isChecked():
            nb_class = self.StewartspinBox_class.value()
            class_breaks_txt = self.StewarttextEdit_breaks.toPlainText()

            if len(class_breaks_txt) > 0:
                class_breaks = parse_class_breaks(class_breaks_txt)
                if class_breaks:
                    nb_class = len(class_breaks)
                else:
                    nb_class = 7
                    self.display_log_error(None, 7)
            else:
                class_breaks = None

            if (resolution == 0) or (span == 0) or (beta == 0):
                self.display_log_error(None, 4)
                return -1

            xi = np.linspace(np.nanmin(x), np.nanmax(x), shape[0])
            yi = np.linspace(np.nanmin(y), np.nanmax(y), shape[1])

            if not class_breaks:
                class_breaks = np.percentile(
                    pot, np.linspace(0.0, 100.0, nb_class + 1))
            try:
                collec_poly = contourf(
                    xi, yi,
                    pot.reshape((shape)).T,
                    class_breaks,
                    vmax=abs(pot).max(),
                    vmin=-abs(pot).max())
            except:
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
                "stewart_potentials_span_{}_beta_{}".format(span, beta),
                "memory")
            err, renderer = render_stewart(
                res_poly, pot_layer,
                levels, nb_class, mask_layer)
            if err:
                self.display_log_error(None, 6)
            pot_layer.setRendererV2(renderer)
            QgsMapLayerRegistry.instance().addMapLayer(pot_layer)
            self.iface.setActiveLayer(pot_layer)
            self.iface.zoomToActiveLayer()
        else:
            bounds = (np.nanmin(x), np.nanmin(y), np.nanmax(x), np.nanmax(y))
            path = save_to_raster(pot, shape, bounds, self.crs.toProj4())
            raster_layer = self.iface.addRasterLayer(
                path, "stewart_potentials_span_{}_beta_{}".format(span, beta))
            raster_layer.setCrs(self.crs)
        # color_raster(raster_layer)

    # Todo : display better information/error message + show progression
    def display_log_error(self, error, msg_nb):
        error_msg = {
            2: ("The computation have been aborted, "
                "please choose a larger resolution value."),
            3: "Error when generating the grid of unknownpts.",
            4: "Span, resolution and beta should not be set to 0.",
            5: "Crs mismatch between point and mask layers. Mask not used.",
            6: "Mask contains invalid geometries. Mask not used.",
            7: ("Unable to parse class breaks (levels must be increasing and "
                "separated by a dash). Using default value of 7 class.")
            }
        QtGui.QMessageBox.information(
            self.iface.mainWindow(), 'Error', error_msg[msg_nb])

    def load_dataset(self):
        home_path = os.getenv('HOMEPATH') \
            or os.getenv('HOME') \
            or os.getenv('USERPROFILE')
        pts_layer = self.iface.addVectorLayer(os.sep.join(
            [home_path, '.qgis2', 'python', 'plugins', 'SpatialPositionModel',
             'test_data', 'paris_mask.geojson']),
            'paris_mask', 'ogr')
        mask_layer = self.iface.addVectorLayer(os.sep.join(
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
        self.StewartComboBox_pts.setLayer(pts_layer)
        self.StewartComboBox_mask.setLayer(mask_layer)
        self.StewartComboBox_field.setCurrentIndex(0)
        self.StewartcomboBox_function.setCurrentIndex(0)
        self.StewartdoubleSpinBox_beta.setValue(3.0)
        self.StewartdoubleSpinBox_span.setValue(1.250)
        self.StewartdoubleSpinBox_resolution.setValue(0.100)
        self.StewarttextEdit_breaks.setPlainText("")
        self.StewartpushButton.setEnabled(True)
