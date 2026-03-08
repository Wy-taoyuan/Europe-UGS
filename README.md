# Urban Tree and Grassland Cover Mapping in European Cities

This repository contains the Google Earth Engine (GEE) code for our study: **A Novel Framework and Automated Algorithm for Subpixel Classification Mapping of Trees and Grasslands in European Cities from Sentinel-1/2 Satellite imagery**.

To support fine-scale urban green space assessment in Europe, we developed an automated sub-pixel classification framework for mapping **trees** and **grasslands** in complex urban environments. The method integrates Sentinel-1 and Sentinel-2 imagery, spectral unmixing, multi-source feature fusion, and Random Forest modeling on the GEE platform. Using this framework, we generated urban tree and grassland cover products for **103 cities across 36 European countries**, providing a refined and systematic perspective on urban green structure at the continental scale.

![](/figure/framework.png)

In this study, we first constructed a training dataset containing both **pure samples** and **non-pure samples** within urban areas. Pure samples were automatically extracted from existing land cover products, while non-pure samples were assigned soft labels through linear spectral mixture analysis. Based on these labeled samples, we built a **117-dimensional feature set** including optical features, radar features, spectral indices, and DEM-derived topographic variables, and then trained a Random Forest model for sub-pixel estimation of tree, grassland, and other land cover fractions.

The released GEE code follows this workflow and supports: urban boundary loading, automatic seed sample generation, optional merging of manually prepared samples, spectral unmixing-based sample refinement, Random Forest prediction, map visualization, and GeoTIFF export of the final products. The current implementation uses Sentinel-2 SR, Sentinel-1 GRD, ESA WorldCover, Dynamic World, and SRTM DEM datasets directly in GEE.

![](/figure/workflow_details.png)

## Results

### Product accuracy

Experimental results show that the proposed method achieves strong agreement with manually annotated reference data. Statistical assessment based on 300 urban blocks (100 m × 100 m) reported **R² = 0.892** and **RMSE = 0.078** for tree cover, and **R² = 0.883** and **RMSE = 0.072** for grassland cover. The generated products also outperform existing large-scale land cover datasets such as **ESA WorldCover**, **Dynamic World V1**, and **ESRI Land Cover** in urban environments.

![](/figure/accuracy.png)

### Urban tree and grassland cover in Europe

Based on the proposed framework, we generated tree and grassland cover products for **103 European cities**. The average urban tree and grassland coverages across Europe are **21.85%** and **26.55%**, respectively, while the population-weighted average coverages are **16.91%** and **19.09%**. Further analysis shows that **latitude** and **precipitation** are dominant drivers of urban tree distribution, whereas **population density** and **longitude** play important roles in explaining grassland patterns.

![](/figure/graphical_abstract.png)

## Package pre-requisites

The code runs on the **Google Earth Engine JavaScript API**. To use this repository, you need:

- A Google Earth Engine account
- Access to the GEE Code Editor
- A valid urban boundary asset
- Optional manually labeled sample assets for refinement

The workflow currently relies on the following datasets in GEE:

- `COPERNICUS/S2_SR`
- `COPERNICUS/S1_GRD`
- `ESA/WorldCover/v200/2021`
- `GOOGLE/DYNAMICWORLD/V1`
- `USGS/SRTMGL1_003`

## Running preparation

Before running the code, set your own urban boundary asset and optional manually prepared sample asset in the user configuration section:

```javascript
var CITY_ASSET = 'users/your_name/your_city_boundary';
var MANUAL_SHP_ASSET = 'users/your_name/your_manual_samples'; // optional
var DRIVE_FOLDER = 'GEE_exports';
var FILE_PREFIX = 'your_city_name';

var YEAR = 2024;
var MONTHS = [1, 12];
var CLOUD_PCT = 20;
var SCALE = 10;

var N_TREE = 1000;
var N_GRASS = 1000;
var N_OTHER = 1000;
var ERODE_PIX = 3;

var UNMIX_BANDS = ['B4', 'B8', 'B11'];
var LAMBDA = 1e-3;
var TRAIN_BUF = 10000;
var KEEP_FRAC = 0.5;
```

The main workflow includes:

1. generating high-purity tree, grass, and other seed samples;
2. merging optional manual samples;
3. refining endmembers through spectral unmixing;
4. training a Random Forest model with multi-source features;
5. predicting Tree / Grass / Other fractions for the target city;
6. exporting the final GeoTIFF products to Google Drive.

### Run in Google Earth Engine

```javascript
runAll();
```

### Output files

The script exports three GeoTIFF products:

- `*_Tree_rf_mprob_thr_f32`
- `*_Grass_rf_mprob_thr_f32`
- `*_Other_rf_mprob_thr_f32`

These files represent the thresholded prediction results for tree, grassland, and other land cover classes, respectively.

## Acknowledgement

This repository is built upon our study on automated sub-pixel mapping of urban trees and grasslands in Europe. We thank the providers of Sentinel-1/2, ESA WorldCover, Dynamic World, SRTM DEM, and related open geospatial datasets used in this work. We also acknowledge prior studies on spectral unmixing, urban vegetation mapping, and large-scale geospatial analysis on Google Earth Engine that inspired this implementation.

## Citation

```bibtex
@misc{GuoUrbanTreeGrassEurope,
  title  = {A Novel Framework and Automated Algorithm for Subpixel Classification Mapping of Trees and Grasslands in European Cities from Sentinel-1/2 Satellite imagery},
  author = {Jianhua Guo and Danfeng Hong and Zhiheng Liu and Xiao Xiang Zhu},
  note   = {Please update this entry with the final publication information},
}
```
