# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SpatialPositionModel_dialog_base.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
from qgis import gui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_SpatialPositionModel_dialog_base(object):
    def setupUi(self, SpatialPositionModel_dialog_base):
        SpatialPositionModel_dialog_base.setObjectName(_fromUtf8("SpatialPositionModel_dialog_base"))
        SpatialPositionModel_dialog_base.resize(760, 628)
        self.tabWidget = QtGui.QTabWidget(SpatialPositionModel_dialog_base)
        self.tabWidget.setGeometry(QtCore.QRect(0, 0, 751, 621))
        font = QtGui.QFont()
        font.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.tabWidget.setFont(font)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.layoutWidget_2 = QtGui.QWidget(self.tab)
        self.layoutWidget_2.setGeometry(QtCore.QRect(20, 0, 711, 391))
        self.layoutWidget_2.setObjectName(_fromUtf8("layoutWidget_2"))
        self.gridLayout_2 = QtGui.QGridLayout(self.layoutWidget_2)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_27 = QtGui.QLabel(self.layoutWidget_2)
        self.label_27.setObjectName(_fromUtf8("label_27"))
        self.gridLayout.addWidget(self.label_27, 12, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.gridLayout.addItem(spacerItem, 10, 2, 1, 1)
        self.mFieldExpressionWidget = gui.QgsFieldExpressionWidget(self.layoutWidget_2)
        self.mFieldExpressionWidget.setFilters(gui.QgsFieldProxyModel.Numeric)
        self.mFieldExpressionWidget.setObjectName(_fromUtf8("mFieldExpressionWidget"))
        self.gridLayout.addWidget(self.mFieldExpressionWidget, 2, 2, 1, 2)
        self.StewartdoubleSpinBox_resolution = QtGui.QDoubleSpinBox(self.layoutWidget_2)
        self.StewartdoubleSpinBox_resolution.setDecimals(4)
        self.StewartdoubleSpinBox_resolution.setMinimum(0.0001)
        self.StewartdoubleSpinBox_resolution.setMaximum(150000.0)
        self.StewartdoubleSpinBox_resolution.setProperty("value", 10.0)
        self.StewartdoubleSpinBox_resolution.setObjectName(_fromUtf8("StewartdoubleSpinBox_resolution"))
        self.gridLayout.addWidget(self.StewartdoubleSpinBox_resolution, 5, 3, 1, 1)
        self.radioButton_vector = QtGui.QRadioButton(self.layoutWidget_2)
        self.radioButton_vector.setChecked(True)
        self.radioButton_vector.setObjectName(_fromUtf8("radioButton_vector"))
        self.gridLayout.addWidget(self.radioButton_vector, 12, 2, 1, 1)
        self.radioButton_raster = QtGui.QRadioButton(self.layoutWidget_2)
        self.radioButton_raster.setObjectName(_fromUtf8("radioButton_raster"))
        self.gridLayout.addWidget(self.radioButton_raster, 12, 3, 1, 1)
        self.StewartComboBox_unknwPts = gui.QgsMapLayerComboBox(self.layoutWidget_2)
        self.StewartComboBox_unknwPts.setEnabled(False)
        self.StewartComboBox_unknwPts.setFilters(gui.QgsMapLayerProxyModel.PointLayer)
        self.StewartComboBox_unknwPts.setObjectName(_fromUtf8("StewartComboBox_unknwPts"))
        self.gridLayout.addWidget(self.StewartComboBox_unknwPts, 8, 2, 1, 2)
        self.label_29 = QtGui.QLabel(self.layoutWidget_2)
        self.label_29.setEnabled(False)
        self.label_29.setObjectName(_fromUtf8("label_29"))
        self.gridLayout.addWidget(self.label_29, 9, 0, 1, 1)
        self.label_30 = QtGui.QLabel(self.layoutWidget_2)
        self.label_30.setEnabled(False)
        self.label_30.setObjectName(_fromUtf8("label_30"))
        self.gridLayout.addWidget(self.label_30, 8, 0, 1, 1)
        self.label_37 = QtGui.QLabel(self.layoutWidget_2)
        self.label_37.setObjectName(_fromUtf8("label_37"))
        self.gridLayout.addWidget(self.label_37, 2, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.gridLayout.addItem(spacerItem1, 3, 2, 1, 1)
        self.StewartdoubleSpinBox_span = QtGui.QDoubleSpinBox(self.layoutWidget_2)
        self.StewartdoubleSpinBox_span.setDecimals(4)
        self.StewartdoubleSpinBox_span.setMinimum(0.0001)
        self.StewartdoubleSpinBox_span.setMaximum(250000.0)
        self.StewartdoubleSpinBox_span.setProperty("value", 10.0)
        self.StewartdoubleSpinBox_span.setObjectName(_fromUtf8("StewartdoubleSpinBox_span"))
        self.gridLayout.addWidget(self.StewartdoubleSpinBox_span, 6, 3, 1, 1)
        self.StewartComboBox_pts = gui.QgsMapLayerComboBox(self.layoutWidget_2)
        self.StewartComboBox_pts.setMinimumContentsLength(0)
        self.StewartComboBox_pts.setFilters(gui.QgsMapLayerProxyModel.PointLayer)
        self.StewartComboBox_pts.setObjectName(_fromUtf8("StewartComboBox_pts"))
        self.gridLayout.addWidget(self.StewartComboBox_pts, 0, 2, 1, 2)
        self.StewartComboBox_field = gui.QgsFieldComboBox(self.layoutWidget_2)
        self.StewartComboBox_field.setMinimumContentsLength(0)
        self.StewartComboBox_field.setFilters(gui.QgsFieldProxyModel.Numeric)
        self.StewartComboBox_field.setObjectName(_fromUtf8("StewartComboBox_field"))
        self.gridLayout.addWidget(self.StewartComboBox_field, 1, 2, 1, 2)
        self.label_6 = QtGui.QLabel(self.layoutWidget_2)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout.addWidget(self.label_6, 7, 2, 1, 1)
        self.label = QtGui.QLabel(self.layoutWidget_2)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_4 = QtGui.QLabel(self.layoutWidget_2)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 5, 2, 1, 1)
        self.StewartdoubleSpinBox_beta = QtGui.QDoubleSpinBox(self.layoutWidget_2)
        self.StewartdoubleSpinBox_beta.setDecimals(2)
        self.StewartdoubleSpinBox_beta.setMinimum(0.0)
        self.StewartdoubleSpinBox_beta.setMaximum(100.0)
        self.StewartdoubleSpinBox_beta.setProperty("value", 1.0)
        self.StewartdoubleSpinBox_beta.setObjectName(_fromUtf8("StewartdoubleSpinBox_beta"))
        self.gridLayout.addWidget(self.StewartdoubleSpinBox_beta, 7, 3, 1, 1)
        self.label_5 = QtGui.QLabel(self.layoutWidget_2)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 6, 2, 1, 1)
        self.label_2 = QtGui.QLabel(self.layoutWidget_2)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.label_24 = QtGui.QLabel(self.layoutWidget_2)
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.gridLayout.addWidget(self.label_24, 4, 2, 1, 1)
        self.StewartcomboBox_function = QtGui.QComboBox(self.layoutWidget_2)
        self.StewartcomboBox_function.setObjectName(_fromUtf8("StewartcomboBox_function"))
        self.StewartcomboBox_function.addItem(_fromUtf8(""))
        self.StewartcomboBox_function.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.StewartcomboBox_function, 4, 3, 1, 1)
        self.StewartComboBox_matdist = gui.QgsMapLayerComboBox(self.layoutWidget_2)
        self.StewartComboBox_matdist.setEnabled(False)
        self.StewartComboBox_matdist.setFilters(gui.QgsMapLayerProxyModel.NoGeometry|gui.QgsMapLayerProxyModel.PluginLayer)
        self.StewartComboBox_matdist.setObjectName(_fromUtf8("StewartComboBox_matdist"))
        self.gridLayout.addWidget(self.StewartComboBox_matdist, 9, 2, 1, 2)
        self.clean_field = QtGui.QToolButton(self.layoutWidget_2)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Monospace"))
        font.setBold(True)
        font.setWeight(75)
        self.clean_field.setFont(font)
        self.clean_field.setText(_fromUtf8(""))
        self.clean_field.setIconSize(QtCore.QSize(16, 16))
        self.clean_field.setObjectName(_fromUtf8("clean_field"))
        self.gridLayout.addWidget(self.clean_field, 1, 1, 1, 1)
        self.clean_custom_expr = QtGui.QToolButton(self.layoutWidget_2)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Monospace"))
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.clean_custom_expr.setFont(font)
        self.clean_custom_expr.setText(_fromUtf8(""))
        self.clean_custom_expr.setIconSize(QtCore.QSize(16, 16))
        self.clean_custom_expr.setObjectName(_fromUtf8("clean_custom_expr"))
        self.gridLayout.addWidget(self.clean_custom_expr, 2, 1, 1, 1)
        self.clean_point_layer = QtGui.QToolButton(self.layoutWidget_2)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Monospace"))
        font.setBold(True)
        font.setWeight(75)
        self.clean_point_layer.setFont(font)
        self.clean_point_layer.setText(_fromUtf8(""))
        self.clean_point_layer.setIconSize(QtCore.QSize(16, 16))
        self.clean_point_layer.setObjectName(_fromUtf8("clean_point_layer"))
        self.gridLayout.addWidget(self.clean_point_layer, 0, 1, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout, 4, 0, 2, 1)
        self.label_7 = QtGui.QLabel(self.layoutWidget_2)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_2.addWidget(self.label_7, 0, 0, 1, 1)
        self.layoutWidget = QtGui.QWidget(self.tab)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 390, 711, 149))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout_3 = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_38 = QtGui.QLabel(self.layoutWidget)
        self.label_38.setObjectName(_fromUtf8("label_38"))
        self.gridLayout_3.addWidget(self.label_38, 2, 0, 1, 1)
        self.label_3 = QtGui.QLabel(self.layoutWidget)
        self.label_3.setEnabled(True)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_3.addWidget(self.label_3, 1, 0, 1, 1)
        self.label_25 = QtGui.QLabel(self.layoutWidget)
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.gridLayout_3.addWidget(self.label_25, 0, 0, 1, 1)
        self.StewartspinBox_class = QtGui.QSpinBox(self.layoutWidget)
        self.StewartspinBox_class.setMinimum(3)
        self.StewartspinBox_class.setMaximum(25)
        self.StewartspinBox_class.setProperty("value", 7)
        self.StewartspinBox_class.setObjectName(_fromUtf8("StewartspinBox_class"))
        self.gridLayout_3.addWidget(self.StewartspinBox_class, 0, 3, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem2, 0, 2, 1, 1)
        self.StewartComboBox_mask = gui.QgsMapLayerComboBox(self.layoutWidget)
        self.StewartComboBox_mask.setEnabled(True)
        self.StewartComboBox_mask.setFilters(gui.QgsMapLayerProxyModel.PolygonLayer)
        self.StewartComboBox_mask.setObjectName(_fromUtf8("StewartComboBox_mask"))
        self.gridLayout_3.addWidget(self.StewartComboBox_mask, 1, 2, 1, 2)
        self.clean_mask_layer = QtGui.QToolButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.clean_mask_layer.sizePolicy().hasHeightForWidth())
        self.clean_mask_layer.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Monospace"))
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.clean_mask_layer.setFont(font)
        self.clean_mask_layer.setText(_fromUtf8(""))
        self.clean_mask_layer.setCheckable(False)
        self.clean_mask_layer.setAutoRepeat(False)
        self.clean_mask_layer.setAutoExclusive(False)
        self.clean_mask_layer.setObjectName(_fromUtf8("clean_mask_layer"))
        self.gridLayout_3.addWidget(self.clean_mask_layer, 1, 1, 1, 1)
        self.StewarttextEdit_breaks = QtGui.QTextEdit(self.layoutWidget)
        self.StewarttextEdit_breaks.setMaximumSize(QtCore.QSize(16777215, 98))
        self.StewarttextEdit_breaks.setTabChangesFocus(False)
        self.StewarttextEdit_breaks.setObjectName(_fromUtf8("StewarttextEdit_breaks"))
        self.gridLayout_3.addWidget(self.StewarttextEdit_breaks, 2, 2, 1, 2)
        self.layoutWidget1 = QtGui.QWidget(self.tab)
        self.layoutWidget1.setGeometry(QtCore.QRect(18, 540, 711, 41))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.layoutWidget1)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.StewartpushButton_clear = QtGui.QPushButton(self.layoutWidget1)
        self.StewartpushButton_clear.setMinimumSize(QtCore.QSize(120, 23))
        self.StewartpushButton_clear.setMaximumSize(QtCore.QSize(200, 23))
        font = QtGui.QFont()
        font.setItalic(True)
        self.StewartpushButton_clear.setFont(font)
        self.StewartpushButton_clear.setObjectName(_fromUtf8("StewartpushButton_clear"))
        self.horizontalLayout.addWidget(self.StewartpushButton_clear)
        self.button_box = QtGui.QDialogButtonBox(self.layoutWidget1)
        self.button_box.setOrientation(QtCore.Qt.Horizontal)
        self.button_box.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.button_box.setObjectName(_fromUtf8("button_box"))
        self.horizontalLayout.addWidget(self.button_box)
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.textBrowser_2 = QtGui.QTextBrowser(self.tab_2)
        self.textBrowser_2.setGeometry(QtCore.QRect(10, 10, 731, 281))
        self.textBrowser_2.setObjectName(_fromUtf8("textBrowser_2"))
        self.layoutWidget_3 = QtGui.QWidget(self.tab_2)
        self.layoutWidget_3.setGeometry(QtCore.QRect(10, 340, 721, 213))
        self.layoutWidget_3.setObjectName(_fromUtf8("layoutWidget_3"))
        self.gridLayout_7 = QtGui.QGridLayout(self.layoutWidget_3)
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.label_26 = QtGui.QLabel(self.layoutWidget_3)
        self.label_26.setObjectName(_fromUtf8("label_26"))
        self.gridLayout_7.addWidget(self.label_26, 1, 0, 1, 1)
        self.pushButton_data = QtGui.QPushButton(self.layoutWidget_3)
        self.pushButton_data.setMinimumSize(QtCore.QSize(75, 0))
        self.pushButton_data.setObjectName(_fromUtf8("pushButton_data"))
        self.gridLayout_7.addWidget(self.pushButton_data, 3, 0, 1, 1, QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
        self.label_36 = QtGui.QLabel(self.layoutWidget_3)
        self.label_36.setObjectName(_fromUtf8("label_36"))
        self.gridLayout_7.addWidget(self.label_36, 0, 0, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem3, 2, 0, 1, 1)
        self.tabWidget.addTab(self.tab_2, _fromUtf8("About"))

        self.retranslateUi(SpatialPositionModel_dialog_base)
        self.tabWidget.setCurrentIndex(0)
        self.StewartComboBox_unknwPts.setCurrentIndex(-1)
        self.StewartComboBox_matdist.setCurrentIndex(-1)
        self.StewartComboBox_mask.setCurrentIndex(-1)
        QtCore.QObject.connect(self.button_box, QtCore.SIGNAL(_fromUtf8("accepted()")), SpatialPositionModel_dialog_base.accept)
        QtCore.QObject.connect(self.button_box, QtCore.SIGNAL(_fromUtf8("rejected()")), SpatialPositionModel_dialog_base.reject)
        QtCore.QMetaObject.connectSlotsByName(SpatialPositionModel_dialog_base)

    def retranslateUi(self, SpatialPositionModel_dialog_base):
        SpatialPositionModel_dialog_base.setWindowTitle(_translate("SpatialPositionModel_dialog_base", "Stewart Potentials", None))
        self.label_27.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Output</p></body></html>", None))
        self.radioButton_vector.setText(_translate("SpatialPositionModel_dialog_base", "Vector", None))
        self.radioButton_raster.setText(_translate("SpatialPositionModel_dialog_base", "Raster", None))
        self.label_29.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Custom distance matrix (opt.)</span></p></body></html>", None))
        self.label_30.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Location of unknown points (opt.)</span></p></body></html>", None))
        self.label_37.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Custom expression (opt.)</span></p></body></html>", None))
        self.label_6.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Beta</p></body></html>", None))
        self.label.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Point layer</p></body></html>", None))
        self.label_4.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Resolution <span style=\" font-size:10pt;\">(</span><span style=\" font-size:10pt; font-style:italic;\">km</span><span style=\" font-size:10pt;\">)</span></p></body></html>", None))
        self.label_5.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Span <span style=\" font-size:10pt;\">(</span><span style=\" font-size:10pt; font-style:italic;\">km</span><span style=\" font-size:10pt;\">)</span></p></body></html>", None))
        self.label_2.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Field name</p></body></html>", None))
        self.label_24.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Spatial interaction function</p></body></html>", None))
        self.StewartcomboBox_function.setItemText(0, _translate("SpatialPositionModel_dialog_base", "Exponential", None))
        self.StewartcomboBox_function.setItemText(1, _translate("SpatialPositionModel_dialog_base", "Pareto", None))
        self.label_7.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-size:12pt; font-weight:600;\">Stewart potentials</span></p></body></html>", None))
        self.label_38.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Break values (opt.)</span></p></body></html>", None))
        self.label_3.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-style:italic;\">Mask layer (optionnal)</span></p></body></html>", None))
        self.label_25.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\">Number of class</p></body></html>", None))
        self.StewartpushButton_clear.setText(_translate("SpatialPositionModel_dialog_base", "Reset field values", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("SpatialPositionModel_dialog_base", "Stewart potentials", None))
        self.textBrowser_2.setHtml(_translate("SpatialPositionModel_dialog_base", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:11pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600; text-decoration: underline;\">Usage</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Point layer</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> The set of known observations to estimate the potentials from.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Field name</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> Name of the variable from </span><span style=\" font-size:9pt; font-style:italic;\">Point layer</span><span style=\" font-size:9pt;\"> from which potentials are computed (quantitative variable without negative values). If no one is selected, 1 is used for each observation of </span><span style=\" font-size:9pt; font-style:italic;\">Point layer</span><span style=\" font-size:9pt;\">.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Mask layer</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> (</span><span style=\" font-size:9pt; font-style:italic;\">optional</span><span style=\" font-size:9pt;\">) The spatial extent of this object is used to create the regularly spaced grid of point and the output will be cropped by it shape.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Spatial interaction function</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> The spatial interaction function to use. Options are &quot;</span><span style=\" font-size:9pt; font-style:italic;\">pareto</span><span style=\" font-size:9pt;\">&quot; (means power law) or &quot;</span><span style=\" font-size:9pt; font-style:italic;\">exponential</span><span style=\" font-size:9pt;\">&quot;. If &quot;</span><span style=\" font-size:9pt; font-style:italic;\">pareto</span><span style=\" font-size:9pt;\">&quot; the interaction is defined as: </span><span style=\" font-family:\'sans-serif\'; font-size:9pt; font-weight:600;\">(1 + alpha * mDistance) ^ (-beta)</span><span style=\" font-size:9pt;\">. If &quot;</span><span style=\" font-size:9pt; font-style:italic;\">exponential</span><span style=\" font-size:9pt;\">&quot; the interaction is defined as: </span><span style=\" font-family:\'sans-serif\'; font-size:9pt; font-weight:600;\">exp(- alpha * mDistance ^ beta)</span><span style=\" font-size:9pt;\">. The alpha parameter is computed from parameters given by the user (</span><span style=\" font-size:9pt; font-style:italic;\">beta</span><span style=\" font-size:9pt;\"> and </span><span style=\" font-size:9pt; font-style:italic;\">span</span><span style=\" font-size:9pt;\">).</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Resolution</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> Resolution (</span><span style=\" font-size:9pt; font-style:italic;\">in map units</span><span style=\" font-size:9pt;\">) of the regular grid the create (and resolution of the output raster for Reilly and Huff algorithms - overriden by the </span><span style=\" font-size:9pt; font-style:italic;\">Location of unknown points</span><span style=\" font-size:9pt;\"> if provided for stewart\'s potentials).</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Span</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> Distance (in map unit, or in the unit of the </span><span style=\" font-size:9pt; font-style:italic;\">Distance matrix</span><span style=\" font-size:9pt;\"> if provided) where the density of probability of the spatial interaction function equals 0.5.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Beta</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> Impedance factor for the spatial interaction function.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Number of classes</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> Number of classes for the polygon shapefile to create.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic; color:#d1d1d1;\">Location of unknown points</span><span style=\" font-size:9pt; font-weight:600; color:#d1d1d1;\"> :</span><span style=\" font-size:9pt; color:#d1d1d1;\"> The set of unknown points for which the function computes the estimates. Using an unknown points layer will override the </span><span style=\" font-size:9pt; font-style:italic; color:#d1d1d1;\">Resolution</span><span style=\" font-size:9pt; color:#d1d1d1;\"> parameter.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic;\">Custom expression</span><span style=\" font-size:9pt; font-weight:600;\"> :</span><span style=\" font-size:9pt;\"> An expression to compute and combine potentials on many variables. Only name of numeric fields and basic operators (</span><span style=\" font-size:10pt; font-weight:600;\">*</span><span style=\" font-size:9pt;\">, </span><span style=\" font-size:10pt; font-weight:600;\">+</span><span style=\" font-size:9pt;\">, </span><span style=\" font-size:10pt; font-weight:600;\">/</span><span style=\" font-size:9pt;\"> and </span><span style=\" font-size:10pt; font-weight:600;\">-</span><span style=\" font-size:9pt;\">) are allowed (so the expression should be looking like </span><span style=\" font-family:\'sans-serif\'; font-size:8pt; font-weight:600;\">&quot;gdppps2008&quot; / &quot;pop2008&quot;</span><span style=\" font-family:\'sans-serif\'; font-size:8pt;\">)</span><span style=\" font-size:9pt;\">. Using a custom expression will override the </span><span style=\" font-size:9pt; font-style:italic;\">Field name</span><span style=\" font-size:9pt;\"> parameter.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:9pt; font-weight:600; font-style:italic; color:#d1d1d1;\">Distance matrix</span><span style=\" font-size:9pt; font-weight:600; color:#d1d1d1;\"> :</span><span style=\" font-size:9pt; color:#d1d1d1;\"> A distance matrix between the provided points of the layer given in </span><span style=\" font-size:9pt; font-style:italic; color:#d1d1d1;\">Point layer</span><span style=\" font-size:9pt; color:#d1d1d1;\"> and in </span><span style=\" font-size:9pt; font-style:italic; color:#d1d1d1;\">Location of unknown points</span><span style=\" font-size:9pt; color:#d1d1d1;\"> fields.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p align=\"center\" style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600; text-decoration: underline;\">References</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'sans-serif\'; font-size:8pt; color:#000000;\">- COMMENGES H., GIRAUD T. (2015) </span><a href=\"https://cran.r-project.org/web/packages/SpatialPosition/vignettes/SpatialPosition.html\"><span style=\" font-family:\'MS Shell Dlg 2\'; font-size:8pt; text-decoration: underline; color:#0000ff;\">Introduction to the SpatialPosition package</span></a></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'sans-serif\'; font-size:8pt; color:#000000;\">- HUFF D. (1964) Defining and Estimating a Trading Area. Journal of Marketing, 28: 34-38.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'sans-serif\'; font-size:8pt; color:#000000;\">- REILLY, W. J. (1931) The law of retail gravitation, W. J. Reilly, New York.</span></p>\n"
"<p style=\" margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'sans-serif\'; font-size:8pt; color:#000000;\">- STEWART J.Q. (1942) &quot;Measure of the influence of a population at a distance&quot;, Sociometry, 5(1): 63-71.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:7px; margin-bottom:7px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:\'sans-serif\'; font-size:8pt; color:#000000;\"><br /></p></body></html>", None))
        self.label_26.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-size:9pt;\">Partial python port of &quot;</span><span style=\" font-size:9pt; font-weight:600;\">SpatialPosition</span><span style=\" font-size:9pt;\">&quot; </span><span style=\" font-size:9pt; font-weight:600;\">R package </span><span style=\" font-size:9pt;\">(with the gracious consent of<br/>the authors H. Commenges &amp; T. Giraud)<br/><br/>-&gt; Original documentation and source code are available on the </span><a href=\"https://cran.r-project.org/package=SpatialPosition\"><span style=\" font-size:9pt; text-decoration: underline; color:#0000ff;\">CRAN</span></a><span style=\" font-size:9pt;\"> or on </span><a href=\"https://github.com/Groupe-ElementR/SpatialPosition\"><span style=\" font-size:9pt; text-decoration: underline; color:#0000ff;\">GitHub</span></a><span style=\" font-size:9pt;\">.<br/>-&gt; R vignettes explaining the main concepts and showing usecases:<br/></span><a href=\"https://cran.r-project.org/web/packages/SpatialPosition/vignettes/SpatialPosition.html\"><span style=\" font-size:9pt; text-decoration: underline; color:#0000ff;\">Introduction to the SpatialPosition package</span></a><span style=\" font-size:9pt;\"> &amp; </span><a href=\"https://cran.r-project.org/web/packages/SpatialPosition/vignettes/StewartExample.html\"><span style=\" font-size:9pt; text-decoration: underline; color:#0000ff;\">Stewart Potential : a Use Case</span></a></p><p align=\"center\"><span style=\" font-size:9pt;\">Plugin authors : M. Viry &amp; T. Giraud<br/></span></p></body></html>", None))
        self.pushButton_data.setText(_translate("SpatialPositionModel_dialog_base", "Load a sample of data", None))
        self.label_36.setText(_translate("SpatialPositionModel_dialog_base", "<html><head/><body><p align=\"center\"><span style=\" font-size:10pt; font-weight:600;\">About</span></p></body></html>", None))
