# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SpatialPositionModel
                                 A QGIS plugin
 Compute spatial interaction models
                             -------------------
        begin                : 2015-12-10
        copyright            : (C) 2015 by #H
        email                : mth@#!.org
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load SpatialPositionModel class from file SpatialPositionModel.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .SpatialPositionModel import SpatialPositionModel
    return SpatialPositionModel(iface)
