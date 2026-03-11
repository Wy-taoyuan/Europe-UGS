/********************* 用户配置区 *********************/
var MANUAL_SHP_ASSET = '';
var ALL_CITIES_ASSET = 'projects/earthengine-legacy/assets/users/wangwenna111/Europe_allcity_Merged_Finall2';

// 导出文件夹
var EXPORT_FOLDER = 'GEE_Unmix';

// 默认年份（UI 会覆盖）
var YEAR = 2024;
var MONTHS = [1, 12];
var CLOUD_PCT = 20;
var SCALE = 10;

// 自动纯样本数量
var N_TREE  = 500;
var N_GRASS = 500;
var N_IMP   = 500;
var N_SOIL  = 500;
var N_WATER = 200;   
var ERODE_PIX = 3;

// Fisher 相关
var FISHER_BANDS = ['B2','B3','B4','B8','B11','B12'];
var FSR_ALPHA = 0.3;   // FSR = 0.3 * sigma(EDA)
var KNN_K = 5;
var REPEAT_SCALE = 20;
var TRAIN_BUF = 10000;

// 基础特征
var S2_BANDS_ALL = [
  'B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
  'NDVI','EVI','GCI','SAVI','NDBI','NDWI'
];
var KEY_BAND = 'B4_mean';


/********************* 通用工具 *********************/
function s2Collection(region, year, months){
  var t0 = ee.Date.fromYMD(year, months[0], 1);
  var t1 = ee.Date.fromYMD(year, months[1], 1).advance(1, 'month');
  return ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(region)
    .filterDate(t0, t1)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PCT));
}

function addIndices(img){
  var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
  var evi = img.expression(
    '2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
    {NIR: img.select('B8'), RED: img.select('B4'), BLUE: img.select('B2')}
  ).rename('EVI');
  var gci = img.expression('(NIR/GREEN)-1', {
    NIR: img.select('B8'),
    GREEN: img.select('B3')
  }).rename('GCI');
  var savi = img.expression(
    '((NIR-RED)/(NIR+RED+0.5))*1.5',
    {NIR: img.select('B8'), RED: img.select('B4')}
  ).rename('SAVI');
  var ndbi = img.normalizedDifference(['B11','B8']).rename('NDBI');
  var ndwi = img.normalizedDifference(['B3','B8']).rename('NDWI');

  var base = img.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12']);
  return base.addBands([ndvi, evi, gci, savi, ndbi, ndwi]);
}

function maskCount(img, region, scale){
  var d = ee.Dictionary(
    img.mask().reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: region,
      scale: scale,
      maxPixels: 1e12,
      bestEffort: true
    })
  );
  var vals = d.values();
  var hasVal = ee.Number(vals.size()).gt(0);
  return ee.Number(ee.Algorithms.If(hasVal, vals.get(0), 0));
}

function stratifiedPointsSafe(seedImg, classValue, nPerClass, region, labelName){
  var cnt = maskCount(seedImg, region, SCALE);
  print('候选像元数量 - ' + labelName + ':', cnt);

  return ee.FeatureCollection(ee.Algorithms.If(
    cnt.gt(0),
    seedImg.updateMask(seedImg).rename('class').multiply(classValue).toInt()
      .stratifiedSample({
        numPoints: nPerClass,
        classBand: 'class',
        region: region,
        scale: SCALE,
        geometries: true,
        classValues: [classValue],
        classPoints: [nPerClass],
        seed: 42
      }),
    ee.FeatureCollection([])
  ));
}

function tryLoadManualSamples(path){
  if (!path || path === '') return {has: false, fc: ee.FeatureCollection([])};
  var has = false;
  var fc = ee.FeatureCollection([]);
  try {
    fc = ee.FeatureCollection(path);
    var c = fc.limit(1).size().getInfo();
    has = (c >= 0);
  } catch (err) {
    print('⚠️ 手动样本资产不可用/不存在：', path);
    has = false;
    fc = ee.FeatureCollection([]);
  }
  return {has: has, fc: fc};
}

function imageMeanDict(img, bands, region){
  return ee.Dictionary(img.select(bands).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: SCALE,
    maxPixels: 1e12,
    bestEffort: true
  }));
}

function zeroVector(n){
  return ee.Array(ee.List.repeat(0, n));
}

function zeroMatrix(n){
  var rows = ee.List.sequence(1, n).map(function(_){
    return ee.List.repeat(0, n);
  });
  return ee.Array(rows);
}

function identityMatrix(n){
  return ee.Array.identity(n);
}

function ensureTGISLabels(fc, cityName){
  return ee.FeatureCollection(fc).map(function(f){
    var cls = ee.Number(ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get('class'), null), -999, f.get('class')
    ));

    var t = ee.Number(ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get('TreeRatio'), null),
      ee.Algorithms.If(cls.eq(1), 1, 0),
      f.get('TreeRatio')
    ));

    var g = ee.Number(ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get('GrassRatio'), null),
      ee.Algorithms.If(cls.eq(2), 1, 0),
      f.get('GrassRatio')
    ));

    var i = ee.Number(ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get('ImpRatio'), null),
      ee.Algorithms.If(cls.eq(3), 1, 0),
      f.get('ImpRatio')
    ));

    var s = ee.Number(ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get('SoilRatio'), null),
      ee.Algorithms.If(cls.eq(4), 1, 0),
      f.get('SoilRatio')
    ));

    var sm = t.add(g).add(i).add(s);

    var tn = ee.Number(ee.Algorithms.If(sm.gt(0), t.divide(sm), 0));
    var gn = ee.Number(ee.Algorithms.If(sm.gt(0), g.divide(sm), 0));
    var inb = ee.Number(ee.Algorithms.If(sm.gt(0), i.divide(sm), 0));
    var sn = ee.Number(ee.Algorithms.If(sm.gt(0), s.divide(sm), 0));

    return f.set({
      TreeRatio: tn,
      GrassRatio: gn,
      ImpRatio: inb,
      SoilRatio: sn,
      city: cityName
    });
  });
}

function setTGISClassFromRatios(fc){
  return ee.FeatureCollection(fc).map(function(f){
    var t = ee.Number(f.get('TreeRatio'));
    var g = ee.Number(f.get('GrassRatio'));
    var i = ee.Number(f.get('ImpRatio'));
    var s = ee.Number(f.get('SoilRatio'));

    var arr = ee.List([t, g, i, s]);
    var mx = ee.Number(arr.reduce(ee.Reducer.max()));
    var idx = ee.Number(arr.indexOf(mx)).add(1); // 1..4
    return f.set('tgis', idx);
  });
}

function seedsFromWorldCoverTGIS(region){
  var wc = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(region);
  var s2Med = s2Collection(region, YEAR, MONTHS).map(addIndices).median();

  var ndvi = s2Med.select('NDVI');
  var ndwi = s2Med.select('NDWI');
  var ndbi = s2Med.select('NDBI');

  var tree  = wc.eq(10);  // Tree cover
  var grass = wc.eq(30);  // Grassland
  var imp   = wc.eq(50);  // Built-up
  var soil  = wc.eq(60);  // Bare / sparse vegetation
  var water = wc.eq(80);  // Water

  tree = tree.focal_min({radius: ERODE_PIX, units:'pixels'})
    .updateMask(ndvi.gte(0.55).and(ndwi.lt(0.15)));

  grass = grass.focal_min({radius: ERODE_PIX, units:'pixels'})
    .updateMask(ndvi.gte(0.30).and(ndvi.lt(0.65)).and(ndwi.lt(0.15)));

  imp = imp.focal_min({radius: ERODE_PIX, units:'pixels'})
    .updateMask(ndvi.lt(0.20).and(ndwi.lt(0.10)).and(ndbi.gt(0)));

  soil = soil.focal_min({radius: ERODE_PIX, units:'pixels'})
    .updateMask(ndvi.lt(0.25).and(ndwi.lt(0.10)).and(ndbi.lt(0.20)));

  water = water.focal_min({radius: 1, units:'pixels'})
    .updateMask(ndwi.gt(0.30));

  return {
    tree: tree.selfMask(),
    grass: grass.selfMask(),
    imp: imp.selfMask(),
    soil: soil.selfMask(),
    water: water.selfMask()
  };
}


/********************* Fisher 变换 *********************/
function fisherFitFromSamples(sampleFC, bandList, classProp, classValues, fallbackImg, fallbackRegion){
  bandList = ee.List(bandList);
  classValues = ee.List(classValues);

  var nb = ee.Number(bandList.size());
  var fallbackDict = imageMeanDict(fallbackImg, bandList, fallbackRegion);

  var classStats = classValues.map(function(cv){
    cv = ee.Number(cv);
    var sub = sampleFC.filter(ee.Filter.eq(classProp, cv));
    var n = ee.Number(sub.size());

    var meanArr = ee.Array(ee.List(bandList).map(function(b){
      b = ee.String(b);
      return ee.Number(ee.Algorithms.If(
        n.gt(0),
        sub.aggregate_mean(b),
        ee.Number(ee.Algorithms.If(fallbackDict.contains(b), fallbackDict.get(b), 0))
      ));
    }));

    var covArr = ee.Array(ee.Algorithms.If(
      n.gt(1),
      sub.reduceColumns({
        reducer: ee.Reducer.centeredCovariance(),
        selectors: bandList
      }).get('covariance'),
      identityMatrix(nb).multiply(1e-6)
    ));

    return ee.Dictionary({
      classValue: cv,
      n: n,
      mean: meanArr,
      cov: covArr
    });
  });

  var totalN = ee.Number(classStats.iterate(function(d, acc){
    d = ee.Dictionary(d);
    return ee.Number(acc).add(ee.Number(d.get('n')));
  }, 0)).max(1);

  var globalMean = ee.Array(classStats.iterate(function(d, acc){
    d = ee.Dictionary(d);
    acc = ee.Array(acc);
    var n = ee.Number(d.get('n'));
    var mu = ee.Array(d.get('mean'));
    return acc.add(mu.multiply(n));
  }, zeroVector(nb))).divide(totalN);

  var Sw = ee.Array(classStats.iterate(function(d, acc){
    d = ee.Dictionary(d);
    acc = ee.Array(acc);
    var n = ee.Number(d.get('n'));
    var cov = ee.Array(d.get('cov'));
    return acc.add(cov.multiply(n.subtract(1).max(1)));
  }, zeroMatrix(nb))).add(identityMatrix(nb).multiply(1e-6));

  var Sb = ee.Array(classStats.iterate(function(d, acc){
    d = ee.Dictionary(d);
    acc = ee.Array(acc);
    var n = ee.Number(d.get('n'));
    var mu = ee.Array(d.get('mean'));
    var diff = mu.subtract(globalMean).reshape([nb, 1]);
    return acc.add(diff.matrixMultiply(diff.matrixTranspose()).multiply(n));
  }, zeroMatrix(nb)));

  var mat = Sw.matrixInverse().matrixMultiply(Sb);
  var eig = mat.eigen();       // 每行：[eigenValue, eigenVector...]
  var eigVecs = eig.slice(1, 1);
  var outDim = ee.Number(classValues.size()).subtract(1); // 5类 -> 4维 Fisher 空间

  // 直接取前 outDim 个特征向量
  var W = eigVecs.slice(0, 0, outDim);

  return {
    W: W,
    outDim: outDim
  };
}

function projectImageToFisher(img, bandList, W, outNames){
  var x = img.select(bandList).toArray().toArray(1); // nb x 1
  var Wimg = ee.Image.constant(W);                   // d x nb
  return Wimg.matrixMultiply(x)
    .arrayProject([0])
    .arrayFlatten([outNames])
    .float();
}

function sampleFisherAndReflectance(fisherImg, coreImg, fc, props){
  return fisherImg.addBands(coreImg).sampleRegions({
    collection: fc,
    properties: props,
    scale: SCALE,
    geometries: true
  }).filter(ee.Filter.neq(FISHER_BANDS[0], -9999));
}


/********************* 端元提纯：KNN + FSR *********************/
function classCenterArray(fc, bandList, classProp, classValue, fallbackImg, fallbackRegion){
  var sub = fc.filter(ee.Filter.eq(classProp, classValue));
  var n = ee.Number(sub.size());
  var fbDict = imageMeanDict(fallbackImg, bandList, fallbackRegion);

  return ee.Array(ee.List(bandList).map(function(b){
    b = ee.String(b);
    return ee.Number(ee.Algorithms.If(
      n.gt(0),
      sub.aggregate_mean(b),
      ee.Number(ee.Algorithms.If(fbDict.contains(b), fbDict.get(b), 0))
    ));
  }));
}

function distanceImageToCenter(fisherImg, fisherBands, centerArr){
  var centerImg = ee.Image.constant(centerArr).rename(fisherBands);
  var diff2 = fisherImg.select(fisherBands).subtract(centerImg).pow(2);
  return diff2.reduce(ee.Reducer.sum()).sqrt().rename('fsDist');
}

function extractEndmemberByKnnAndFSR(fisherImg, coreImg, knnClassImg, classValue, pureTrainFC, fisherBands, region){
  var center = classCenterArray(pureTrainFC, fisherBands, 'tgis', classValue, fisherImg, region);
  var dImg = distanceImageToCenter(fisherImg, fisherBands, center);
  var candidateMask = knnClassImg.eq(classValue);

  var sigmaDict = ee.Dictionary(
    dImg.updateMask(candidateMask).reduceRegion({
      reducer: ee.Reducer.stdDev(),
      geometry: region,
      scale: SCALE,
      maxPixels: 1e12,
      bestEffort: true
    })
  );

  var sigma = ee.Number(ee.Algorithms.If(
    sigmaDict.contains('fsDist'),
    sigmaDict.get('fsDist'),
    0.05
  ));

  var fsr = sigma.multiply(FSR_ALPHA).max(0.01);
  var keptMask = candidateMask.and(dImg.lte(fsr));

  var fisherMeanDict = ee.Dictionary(
    fisherImg.select(fisherBands).updateMask(keptMask).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: SCALE,
      maxPixels: 1e12,
      bestEffort: true
    })
  );

  var reflMeanDict = ee.Dictionary(
    coreImg.select(FISHER_BANDS).updateMask(keptMask).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: SCALE,
      maxPixels: 1e12,
      bestEffort: true
    })
  );

  var fisherFallback = classCenterArray(pureTrainFC, fisherBands, 'tgis', classValue, fisherImg, region);
  var reflFallback   = classCenterArray(pureTrainFC, FISHER_BANDS, 'tgis', classValue, coreImg, region);

  var fisherMean = ee.Array(ee.List(fisherBands).map(function(b, idx){
    b = ee.String(b);
    idx = ee.Number(idx);
    return ee.Number(ee.Algorithms.If(
      fisherMeanDict.contains(b),
      fisherMeanDict.get(b),
      fisherFallback.get([idx])
    ));
  })).reshape([fisherBands.length, 1]);

  var reflMean = ee.Array(ee.List(FISHER_BANDS).map(function(b, idx){
    b = ee.String(b);
    idx = ee.Number(idx);
    return ee.Number(ee.Algorithms.If(
      reflMeanDict.contains(b),
      reflMeanDict.get(b),
      reflFallback.get([idx])
    ));
  })).reshape([FISHER_BANDS.length, 1]);

  return {
    center: center,
    fsr: fsr,
    fisher: fisherMean,
    refl: reflMean,
    mask: keptMask
  };
}


/********************* Fisher 空间约束解混 *********************/
// 用“sum-to-one 增广最小二乘 + 非负裁剪 + 归一化”近似 FCLS
function fisherConstrainedUnmix(fisherImg, fisherBands, E, outNames){
  var p = outNames.length;

  var oneRow = ee.Array([ee.List.repeat(1, p)]); // 1 x p
  var Eaug = ee.Array.cat([E, oneRow], 0);       // (d+1) x p

  var G = Eaug.matrixTranspose()
    .matrixMultiply(Eaug)
    .add(identityMatrix(p).multiply(1e-6));

  var pinv = G.matrixInverse().matrixMultiply(Eaug.matrixTranspose());

  var yaug = fisherImg.select(fisherBands)
    .addBands(ee.Image.constant(1).rename('const1'))
    .toArray()
    .toArray(1); // (d+1) x 1

  var f = ee.Image.constant(pinv).matrixMultiply(yaug)
    .arrayProject([0])
    .arrayFlatten([outNames]);

  var fPos = f.max(0);
  var sm = fPos.reduce(ee.Reducer.sum()).rename('sumF');

  return fPos.divide(sm).where(sm.eq(0), 0).float();
}

function reconRmseFromFisher(fisherImg, fisherBands, E, abundImg, outNames){
  var y = fisherImg.select(fisherBands).toArray().toArray(1);   // d x 1
  var f = abundImg.select(outNames).toArray().toArray(1);       // p x 1
  var recon = ee.Image.constant(E).matrixMultiply(f);           // d x 1
  var resid = y.subtract(recon);

  return resid.multiply(resid)
    .arrayReduce(ee.Reducer.mean(), [0])
    .arrayProject([0])
    .arrayFlatten([['UMX_RMSE']])
    .sqrt()
    .float();
}

function expandSoftLabelsToPseudoClasses(fc, repeatScale){
  var n = ee.Number(fc.size());
  var list = fc.toList(n);
  var empty = ee.FeatureCollection([]);

  return ee.FeatureCollection(ee.Algorithms.If(
    n.gt(0),
    ee.List.sequence(0, n.subtract(1)).iterate(function(idx, acc){
      idx = ee.Number(idx);
      var f = ee.Feature(list.get(idx));

      var probs = ee.List([
        ee.Number(f.get('TreeRatio')),
        ee.Number(f.get('GrassRatio')),
        ee.Number(f.get('ImpRatio')),
        ee.Number(f.get('SoilRatio'))
      ]);

      var mx = ee.Number(probs.reduce(ee.Reducer.max()));
      var argmax = ee.Number(probs.indexOf(mx));

      var fcOne = ee.FeatureCollection(
        ee.List.sequence(0, 3).map(function(c){
          c = ee.Number(c);
          var p = ee.Number(probs.get(c));

          var nRep = p.multiply(repeatScale).round().toInt();
          nRep = ee.Number(ee.Algorithms.If(c.eq(argmax), nRep.max(1), nRep));

          return ee.FeatureCollection(ee.Algorithms.If(
            nRep.gt(0),
            ee.FeatureCollection(ee.List.sequence(1, nRep).map(function(k){
              return ee.Feature(null, f.toDictionary()).set('cls', c);
            })),
            ee.FeatureCollection([])
          ));
        })
      ).flatten();

      return ee.FeatureCollection(acc).merge(fcOne);
    }, empty),
    empty
  ));
}


/********************* 主分析 *********************/
function runAnalysis(region, cityName, samplesAll){
  var trainAOI = samplesAll.geometry().buffer(TRAIN_BUF);
  var samplesLabeled = ensureTGISLabels(samplesAll, cityName);

  var pureSamples = ee.FeatureCollection(samplesLabeled.filter(
    ee.Filter.or(
      ee.Filter.eq('TreeRatio', 1),
      ee.Filter.eq('GrassRatio', 1),
      ee.Filter.eq('ImpRatio', 1),
      ee.Filter.eq('SoilRatio', 1)
    )
  ));

  var nonPureSamples = ee.FeatureCollection(samplesLabeled.filter(
    ee.Filter.and(
      ee.Filter.neq('TreeRatio', 1),
      ee.Filter.neq('GrassRatio', 1),
      ee.Filter.neq('ImpRatio', 1),
      ee.Filter.neq('SoilRatio', 1)
    )
  ));

  pureSamples = setTGISClassFromRatios(pureSamples);

  /**** 2) 影像准备 ****/
  var statsReducer = ee.Reducer.mean()
    .combine(ee.Reducer.min(), '', true)
    .combine(ee.Reducer.max(), '', true)
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.percentile([25,50,75], ['p25','p50','p75']), '', true);

  // 城市
  var s2_city_col = s2Collection(region, YEAR, MONTHS);
  var s2idx_city = s2_city_col.map(addIndices);
  var feature_city_base = s2idx_city.select(S2_BANDS_ALL).reduce(statsReducer).clip(region).unmask(-9999);
  var s2core_city = s2_city_col.median().select(FISHER_BANDS).clip(region).unmask(-9999);

  // 训练 AOI
  var s2_train_col = s2Collection(trainAOI, YEAR, MONTHS);
  var s2idx_train = s2_train_col.map(addIndices);
  var feature_train_base = s2idx_train.select(S2_BANDS_ALL).reduce(statsReducer).clip(trainAOI).unmask(-9999);
  var s2core_train = s2_train_col.median().select(FISHER_BANDS).clip(trainAOI).unmask(-9999);

  // S1
  var s1_city = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(region)
    .filterDate(ee.Date.fromYMD(YEAR, MONTHS[0], 1), ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1,'month'))
    .filter(ee.Filter.eq('instrumentMode','IW'))
    .filter(ee.Filter.eq('orbitProperties_pass','DESCENDING'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'));

  var vv_city = s1_city.select('VV').median().rename('VV_median');
  var vh_city = s1_city.select('VH').median().rename('VH_median');

  var s1_train = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(trainAOI)
    .filterDate(ee.Date.fromYMD(YEAR, MONTHS[0], 1), ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1,'month'))
    .filter(ee.Filter.eq('instrumentMode','IW'))
    .filter(ee.Filter.eq('orbitProperties_pass','DESCENDING'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'));

  var vv_train = s1_train.select('VV').median().rename('VV_median');
  var vh_train = s1_train.select('VH').median().rename('VH_median');

  // DEM
  var dem = ee.Image('USGS/SRTMGL1_003');
  var elevation = dem.rename('elevation');
  var slope = ee.Terrain.slope(dem).rename('slope');
  var aspect = ee.Terrain.aspect(dem).rename('aspect');

  var waterSeeds = seedsFromWorldCoverTGIS(region);
  var waterPts = stratifiedPointsSafe(waterSeeds.water, 5, N_WATER, region, 'Water(aux)')
    .map(function(f){ return f.set('tgis', 5); });

  /**** 4) Fisher 拟合 ****/
  var fisherTrainFC = sampleFisherAndReflectance(
    s2core_city,
    s2core_city,
    pureSamples.merge(waterPts),
    ['tgis']
  );

  var fisherFit = fisherFitFromSamples(
    fisherTrainFC,
    FISHER_BANDS,
    'tgis',
    [1,2,3,4,5],   // T,G,I,S,W -> 4个 Fisher 轴
    s2core_city,
    region
  );

  var fisherNames = ['F1','F2','F3','F4'];

  var fisher_city = projectImageToFisher(s2core_city, FISHER_BANDS, fisherFit.W, fisherNames);
  var fisher_train = projectImageToFisher(s2core_train, FISHER_BANDS, fisherFit.W, fisherNames);

  /**** 5) Fisher 空间 KNN 筛候选端元 ****/
  var pureCity = pureSamples.filterBounds(region);

  var pureFishCityTrain = sampleFisherAndReflectance(
    fisher_city,
    s2core_city,
    pureCity,
    ['tgis']
  );

  var knn = ee.Classifier.smileKNN(KNN_K).train({
    features: pureFishCityTrain,
    classProperty: 'tgis',
    inputProperties: fisherNames
  });

  var knnClassImg = fisher_city.classify(knn).rename('knnClass');

  /**** 6) 各类端元：KNN + FSR 提纯 ****/
  var emTree = extractEndmemberByKnnAndFSR(
    fisher_city, s2core_city, knnClassImg, 1, pureFishCityTrain, fisherNames, region
  );
  var emGrass = extractEndmemberByKnnAndFSR(
    fisher_city, s2core_city, knnClassImg, 2, pureFishCityTrain, fisherNames, region
  );
  var emImp = extractEndmemberByKnnAndFSR(
    fisher_city, s2core_city, knnClassImg, 3, pureFishCityTrain, fisherNames, region
  );
  var emSoil = extractEndmemberByKnnAndFSR(
    fisher_city, s2core_city, knnClassImg, 4, pureFishCityTrain, fisherNames, region
  );

  print('FSR(Tree/Grass/Imp/Soil):', emTree.fsr, emGrass.fsr, emImp.fsr, emSoil.fsr);

  // Fisher 端元矩阵：4 x 4，列顺序 T/G/I/S
  var E = ee.Array.cat([emTree.fisher, emGrass.fisher, emImp.fisher, emSoil.fisher], 1);

  /**** 7) Fisher 空间解混：得到软标签 ****/
  var abund_city = fisherConstrainedUnmix(
    fisher_city, fisherNames, E,
    ['UMX_Tree','UMX_Grass','UMX_Imp','UMX_Soil']
  );

  var abund_train = fisherConstrainedUnmix(
    fisher_train, fisherNames, E,
    ['UMX_Tree','UMX_Grass','UMX_Imp','UMX_Soil']
  );

  var rmse_city = reconRmseFromFisher(
    fisher_city, fisherNames, E, abund_city,
    ['UMX_Tree','UMX_Grass','UMX_Imp','UMX_Soil']
  );

  var rmse_train = reconRmseFromFisher(
    fisher_train, fisherNames, E, abund_train,
    ['UMX_Tree','UMX_Grass','UMX_Imp','UMX_Soil']
  );

  
  /**** 8) 非纯样本赋软标签 ****/
  var mixedSoft = ee.FeatureCollection(ee.Algorithms.If(
    nonPureSamples.size().gt(0),
    abund_city.sampleRegions({
      collection: nonPureSamples,
      scale: SCALE,
      geometries: true
    }).map(function(f){
      return f.set({
        TreeRatio:  f.get('UMX_Tree'),
        GrassRatio: f.get('UMX_Grass'),
        ImpRatio:   f.get('UMX_Imp'),
        SoilRatio:  f.get('UMX_Soil')
      });
    }),
    ee.FeatureCollection([])
  ));

  var labeledSamples = pureSamples.merge(mixedSoft).distinct(['.geo']);
  print('用于 RF 的总标签样本数:', labeledSamples.size());

  /**** 9) 构造 RF 特征 ****/
  var feature_city = feature_city_base
    .addBands([vh_city, vv_city, aspect, elevation, slope])
    .addBands(abund_city)
    .addBands(rmse_city)
    .clip(region);

  var feature_train = feature_train_base
    .addBands([vh_train, vv_train, aspect, elevation, slope])
    .addBands(abund_train)
    .addBands(rmse_train)
    .clip(trainAOI);

  var allBands = feature_city.bandNames();
  print('feature_city band count:', allBands.size());

  var trainBase = feature_train.sampleRegions({
    collection: labeledSamples,
    properties: ['TreeRatio','GrassRatio','ImpRatio','SoilRatio'],
    scale: SCALE,
    geometries: false
  })
  .filter(ee.Filter.neq(KEY_BAND, -9999))
  .filter(ee.Filter.notNull(['TreeRatio','GrassRatio','ImpRatio','SoilRatio']));
  var trainFC = expandSoftLabelsToPseudoClasses(trainBase, REPEAT_SCALE);
  var rfMulti = ee.Classifier.smileRandomForest(200)
    .setOutputMode('MULTIPROBABILITY')
    .train({
      features: trainFC,
      classProperty: 'cls',   // 0=Tree,1=Grass,2=Imp,3=Soil
      inputProperties: allBands
    });

  var probArr = feature_city.classify(rfMulti);

  var probs = probArr.arrayFlatten([['Tree','Grass','Impervious','Soil']]);
  var sumP = probs.reduce(ee.Reducer.sum()).rename('sumP');
  var probs1 = probs.divide(sumP).where(sumP.eq(0), 0).float();
  var rgbSoft = ee.Image.cat([
    abund_city.select('UMX_Imp'),
    abund_city.select('UMX_Tree'),
    abund_city.select('UMX_Grass')
  ]);
  Map.addLayer(rgbSoft, {min:[0,0,0], max:[1,1,1]}, cityName + ' RGB SoftLabel (I/T/G)');

  var rgbProb = ee.Image.cat([
    probs1.select('Impervious'),
    probs1.select('Tree'),
    probs1.select('Grass')
  ]);
  Map.addLayer(rgbProb, {min:[0,0,0], max:[1,1,1]}, cityName + ' RGB Prob (I/T/G)');

  /**** 11) 导出 ****/
  function exportToDrive(image, description, fileNamePrefix){
    Export.image.toDrive({
      image: image.clamp(0,1).toFloat(),
      description: description,
      folder: EXPORT_FOLDER,
      fileNamePrefix: fileNamePrefix,
      region: region,
      scale: SCALE,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF',
      formatOptions: {cloudOptimized: true}
    });
  }
  exportToDrive(probs1.select('Tree'),       cityName + '_Prob_Tree',       cityName + '_Prob_Tree');
  exportToDrive(probs1.select('Grass'),      cityName + '_Prob_Grass',      cityName + '_Prob_Grass');
}


/********************* UI *********************/
var allCitiesFC = ee.FeatureCollection(ALL_CITIES_ASSET);

var CURRENT_CITY_NAME = null;
var CURRENT_REGION = null;

var uiPanel = ui.Panel({
  style: {width: '430px', position: 'top-left', padding: '8px'}
});

uiPanel.add(ui.Label({
  value: '📌 城市整年分析（Fisher 软标签 + RF 多概率）',
  style: {fontWeight: 'bold', fontSize: '14px'}
}));

var citySelect = ui.Select({placeholder: '加载城市列表中...'});
allCitiesFC.aggregate_array('city').distinct().sort().evaluate(function(list){
  citySelect.items().reset(list);
  citySelect.setPlaceholder('选择城市');
});

var yearBox = ui.Textbox({
  placeholder: '年份（默认 ' + YEAR + '）',
  value: ''
});

var confirmBtn = ui.Button({
  label: '✅ 显示边界',
  style: {stretch: 'horizontal'}
});

var runBtn = ui.Button({
  label: '▶️ 运行整年分析',
  style: {stretch: 'horizontal'}
});

uiPanel.add(ui.Panel([ui.Label('城市'), citySelect], ui.Panel.Layout.flow('horizontal')));
uiPanel.add(ui.Panel([ui.Label('年份(可选)'), yearBox], ui.Panel.Layout.flow('horizontal')));
uiPanel.add(confirmBtn);
uiPanel.add(runBtn);

Map.add(uiPanel);

confirmBtn.onClick(function(){
  var cityName = citySelect.getValue();
  if (!cityName){
    print('⚠️ 请先选择城市');
    return;
  }

  var cityFC = allCitiesFC.filter(ee.Filter.eq('city', cityName));
  var region = cityFC.union(1).geometry();

  Map.clear();
  Map.add(uiPanel);
  Map.centerObject(region, 10);
  Map.addLayer(cityFC.style({color:'black', fillColor:'00000000', width:2}), {}, cityName + ' 边界');

  CURRENT_CITY_NAME = cityName;
  CURRENT_REGION = region;
});

runBtn.onClick(function(){
  var cityName = CURRENT_CITY_NAME || citySelect.getValue();
  if (!cityName){
    print('⚠️ 请先选择城市（可先点“显示边界”）');
    return;
  }

  var region = CURRENT_REGION;
  if (!region){
    var cityFC = allCitiesFC.filter(ee.Filter.eq('city', cityName));
    region = cityFC.union(1).geometry();
  }

  var yrTxt = yearBox.getValue();
  YEAR = (yrTxt && yrTxt.trim() !== '') ? parseInt(yrTxt, 10) : YEAR;
  MONTHS = [1, 12];

  Map.clear();
  Map.add(uiPanel);
  Map.centerObject(region, 10);
  Map.addLayer(
    ee.FeatureCollection(ALL_CITIES_ASSET)
      .filter(ee.Filter.eq('city', cityName))
      .style({color:'black', fillColor:'00000000', width:2}),
    {},
    cityName + ' 边界'
  );

  print('▶️ 开始分析：城市=' + cityName + '，YEAR=' + YEAR + '（整年 1–12 月）');
  runAllForUI(region, cityName);
});


/********************* 运行入口 *********************/
function runAllForUI(region, cityName){
  var seeds = seedsFromWorldCoverTGIS(region);

  var treePts  = stratifiedPointsSafe(seeds.tree,  1, N_TREE,  region, 'Tree');
  var grassPts = stratifiedPointsSafe(seeds.grass, 2, N_GRASS, region, 'Grass');
  var impPts   = stratifiedPointsSafe(seeds.imp,   3, N_IMP,   region, 'Impervious');
  var soilPts  = stratifiedPointsSafe(seeds.soil,  4, N_SOIL,  region, 'Soil');

  var autoPure = ee.FeatureCollection(treePts)
    .merge(grassPts)
    .merge(impPts)
    .merge(soilPts)
    .map(function(f){
      var c = ee.Number(f.get('class'));
      return f.set({
        TreeRatio:  ee.Number(ee.Algorithms.If(c.eq(1), 1, 0)),
        GrassRatio: ee.Number(ee.Algorithms.If(c.eq(2), 1, 0)),
        ImpRatio:   ee.Number(ee.Algorithms.If(c.eq(3), 1, 0)),
        SoilRatio:  ee.Number(ee.Algorithms.If(c.eq(4), 1, 0)),
        city: cityName
      });
    })
    .distinct(['.geo']);

  var manualLoad = tryLoadManualSamples(MANUAL_SHP_ASSET);

  var manualFiltered = ee.FeatureCollection(ee.Algorithms.If(
    manualLoad.has,
    ensureTGISLabels(ee.FeatureCollection(manualLoad.fc).filterBounds(region), cityName),
    ee.FeatureCollection([])
  ));

  var samplesAll = ee.FeatureCollection(ee.Algorithms.If(
    manualLoad.has,
    autoPure.merge(manualFiltered),
    autoPure
  )).distinct(['.geo']);

  var nTree  = samplesAll.filter(ee.Filter.eq('TreeRatio', 1)).size();
  var nGrass = samplesAll.filter(ee.Filter.eq('GrassRatio', 1)).size();
  var nImp   = samplesAll.filter(ee.Filter.eq('ImpRatio', 1)).size();
  var nSoil  = samplesAll.filter(ee.Filter.eq('SoilRatio', 1)).size();
  var nAll   = samplesAll.size();

  runAnalysis(region, cityName, samplesAll);
}
