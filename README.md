# Urban Tree and Grassland Cover Mapping in European Cities

This repository contains the Google Earth Engine (GEE) code for our study: **A Novel Framework and Automated Algorithm for Subpixel Classification Mapping of Trees and Grasslands in European Cities from Sentinel-1/2 Satellite imagery**.

To support fine-scale urban green space assessment in Europe, we developed an automated sub-pixel classification framework for mapping trees and grasslands in complex urban environments. This method is based on Sentinel-1/2 images, and uses the random forest (RF) technique on the Google Earth Engine (GEE) platform to achieve sub-pixel mapping of urban trees and grasslands. Using this method, we created 10-m resolution tree and grassland cover products for 103 European cities and made them freely available to the community (https://doi.org/10.6084/m9.figshare.31566298).

![](/figure/framework.png)

In this study, we address training dataset construction, classification model training, the prediction of urban tree and grassland products, and the evaluation of vegetation coverage.

## Results

### UTC and UGC detection result

![](/figure/Munich_detection_result.png)


## Package pre-requisites

The code runs on the **Google Earth Engine JavaScript API**. To use this repository, you need:

- A Google Earth Engine account
- Access to the GEE Code Editor
- A valid urban boundary asset

The workflow currently relies on the following datasets in GEE:

- `COPERNICUS/S2_SR`
- `COPERNICUS/S1_GRD`
- `ESA/WorldCover/v200/2021`
- `GOOGLE/DYNAMICWORLD/V1`
- `USGS/SRTMGL1_003`

## Running preparation

Before running the code, set your own urban boundary asset in the user configuration section:

```javascript
var CITY_ASSET = 'users/your_name/your_city_boundary';
```

### Run in Google Earth Engine

```javascript
runAll();
```

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
