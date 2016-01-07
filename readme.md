SpatialPositionModel
--------------------

Partial python port of R SpatialPosition package.
Original package (documentation and R source code) are available on [GitHub](https://github.com/Groupe-ElementR/SpatialPosition) or on the [CRAN](https://cran.r-project.org/web/packages/SpatialPosition/).

Three models are implemented :

**Stewart Potentials**  
<img src="misc/stewart_screenshot.png" width="433" height="230">

**Reilly catchment areas**  
<img src="misc/reilly_screenshot.png" width="384" height="209">

**Huff probabilistic catchment areas**  
<img src="misc/huff_screenshot.png" width="411" height="207">

#### Installation :
This plugin hasn't been submitted on QGIS plugin repository yet.
It can be installed like this :
```
git clone https://github.com/mthh/SpatialPositionModel.git
cd SpatialPositionModel
make deploy
```

The plugin is now installed and will be available to load in the plugin window of QGIS (think also to allow experimental plugins) under the name **SpatialPositionModel**.
