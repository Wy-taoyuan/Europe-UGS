/********************* 用户配置区 *********************/
var MANUAL_SHP_ASSET = '';
var DRIVE_FOLDER = 'GEE_Unmix';

var YEAR = 2024;
var MONTHS = [1, 12];
var CLOUD_PCT = 20;
var SCALE = 10;

var N_TREE = 400;
var N_GRASS = 400;
var N_IMPERV = 400;
var N_SOIL = 400;
var ERODE_PIX = 3;

var UNMIX_BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
var LAMBDA = 1e-3;

var KNN_K = 60;
var FSR_PCTL = 85;

var TRAIN_BUF = 10000;
var KEY_BAND = 'B4_mean';

var ALL_CITIES_ASSET = '';
var allCitiesFC = ee.FeatureCollection(ALL_CITIES_ASSET);

/********************* 通用小工具 *********************/
function s2Collection(region, year, months) {
  var t0 = ee.Date.fromYMD(year, months[0], 1);
  var t1 = ee.Date.fromYMD(year, months[1], 1).advance(1, 'month');
  return ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(region)
    .filterDate(t0, t1)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PCT));
}

function addIndices(img) {
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = img.expression(
    '2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
    {NIR: img.select('B8'), RED: img.select('B4'), BLUE: img.select('B2')}
  ).rename('EVI');
  var gci = img.expression('(NIR/GREEN)-1', {
    NIR: img.select('B8'),
    GREEN: img.select('B3')
  }).rename('GCI');
  var savi = img.expression('((NIR-RED)/(NIR+RED+0.5))*1.5', {
    NIR: img.select('B8'),
    RED: img.select('B4')
  }).rename('SAVI');
  var ndbi = img.normalizedDifference(['B11', 'B8']).rename('NDBI');
  var ndwi = img.normalizedDifference(['B3', 'B8']).rename('NDWI');

  var base = img.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12']);
  return base.addBands([ndvi, evi, gci, savi, ndbi, ndwi]);
}

function maskCount(img, region, scale) {
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

function stratifiedPointsSafe(seedImg, classValue, nPerClass, region) {
  var cnt = maskCount(seedImg, region, SCALE);
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

function tryLoadManualSamples(path) {
  if (!path || path === '') return {has: false, fc: ee.FeatureCollection([])};
  var has = false;
  var fc = ee.FeatureCollection([]);
  try {
    fc = ee.FeatureCollection(path);
    var c = fc.limit(1).size().getInfo();
    has = (c >= 0);
  } catch (err) {
    has = false;
    fc = ee.FeatureCollection([]);
  }
  return {has: has, fc: fc};
}

function ensureRF3OneHot(fc) {
  return ee.FeatureCollection(fc).map(function(f) {
    var hasTree = f.get('TreeRatio');
    var hasGrass = f.get('GrassRatio');
    var hasOther = f.get('OtherRatio');

    var needFill = ee.Algorithms.If(
      ee.Algorithms.IsEqual(hasTree, null),
      true,
      ee.Algorithms.If(
        ee.Algorithms.IsEqual(hasGrass, null),
        true,
        ee.Algorithms.IsEqual(hasOther, null)
      )
    );

    var c = ee.Number(ee.Algorithms.If(f.get('class'), f.get('class'), -999));

    var lab = ee.Dictionary(ee.Algorithms.If(
      needFill,
      ee.Algorithms.If(
        c.eq(1), ee.Dictionary({'TreeRatio':1,'GrassRatio':0,'OtherRatio':0}),
        ee.Algorithms.If(
          c.eq(2), ee.Dictionary({'TreeRatio':0,'GrassRatio':1,'OtherRatio':0}),
          ee.Algorithms.If(
            c.eq(3), ee.Dictionary({'TreeRatio':0,'GrassRatio':0,'OtherRatio':1}),
            ee.Dictionary({'TreeRatio':0,'GrassRatio':0,'OtherRatio':0})
          )
        )
      ),
      ee.Dictionary({})
    ));

    return f.set(lab);
  }).filter(
    ee.Filter.or(
      ee.Filter.eq('TreeRatio', 1),
      ee.Filter.or(
        ee.Filter.eq('GrassRatio', 1),
        ee.Filter.eq('OtherRatio', 1)
      )
    )
  );
}

/********************* TGIS *********************/
function seedsFromWorldCoverTGIS(region) {
  var wc = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(region);
  var s2Med = s2Collection(region, YEAR, MONTHS).map(addIndices).median();

  var ndvi = s2Med.select('NDVI');
  var ndwi = s2Med.select('NDWI');

  var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterBounds(region)
    .filterDate(
      ee.Date.fromYMD(YEAR, MONTHS[0], 1),
      ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1, 'month')
    )
    .select('label')
    .mode()
    .clip(region);

  var built = dw.eq(6);
  var bare = dw.eq(7);
  var water = wc.eq(80).or(ndwi.gte(0.3));

  var tree = wc.eq(10)
    .focal_min({radius: ERODE_PIX, units: 'pixels'})
    .updateMask(ndvi.gte(0.55))
    .updateMask(water.not());

  var grass = wc.eq(30)
    .focal_min({radius: ERODE_PIX, units: 'pixels'})
    .updateMask(ndvi.gte(0.25).and(ndvi.lt(0.75)))
    .updateMask(water.not());

  var imperv = wc.eq(50).or(built)
    .focal_min({radius: ERODE_PIX, units: 'pixels'})
    .updateMask(ndvi.lt(0.25))
    .updateMask(ndwi.lt(0.15));

  var soil = wc.eq(60).or(bare)
    .focal_min({radius: ERODE_PIX, units: 'pixels'})
    .updateMask(ndvi.lt(0.20))
    .updateMask(ndwi.lt(0.15))
    .updateMask(imperv.not())
    .updateMask(water.not());

  return {
    tree: tree.selfMask(),
    grass: grass.selfMask(),
    imperv: imperv.selfMask(),
    soil: soil.selfMask()
  };
}

function buildAutoSamplesTGIS(region, cityName) {
  var seeds = seedsFromWorldCoverTGIS(region);

  var treePts = stratifiedPointsSafe(seeds.tree, 1, N_TREE, region);
  var grassPts = stratifiedPointsSafe(seeds.grass, 2, N_GRASS, region);
  var impervPts = stratifiedPointsSafe(seeds.imperv, 3, N_IMPERV, region);
  var soilPts = stratifiedPointsSafe(seeds.soil, 4, N_SOIL, region);

  var auto4 = ee.FeatureCollection(treePts)
    .merge(grassPts)
    .merge(impervPts)
    .merge(soilPts)
    .map(function(f) {
      var c = ee.Number(f.get('class'));

      var tgis = ee.Dictionary(ee.Algorithms.If(
        c.eq(1), ee.Dictionary({'TreeRatio':1,'GrassRatio':0,'ImpervRatio':0,'SoilRatio':0}),
        ee.Algorithms.If(
          c.eq(2), ee.Dictionary({'TreeRatio':0,'GrassRatio':1,'ImpervRatio':0,'SoilRatio':0}),
          ee.Algorithms.If(
            c.eq(3), ee.Dictionary({'TreeRatio':0,'GrassRatio':0,'ImpervRatio':1,'SoilRatio':0}),
            ee.Dictionary({'TreeRatio':0,'GrassRatio':0,'ImpervRatio':0,'SoilRatio':1})
          )
        )
      ));

      var tgo = ee.Dictionary(ee.Algorithms.If(
        c.eq(1), ee.Dictionary({'RF_TreeRatio':1,'RF_GrassRatio':0,'RF_OtherRatio':0}),
        ee.Algorithms.If(
          c.eq(2), ee.Dictionary({'RF_TreeRatio':0,'RF_GrassRatio':1,'RF_OtherRatio':0}),
          ee.Dictionary({'RF_TreeRatio':0,'RF_GrassRatio':0,'RF_OtherRatio':1})
        )
      ));

      return f.set(tgis).set(tgo).set({'city': cityName});
    })
    .distinct(['.geo']);

  return auto4;
}

/********************* Fisher-transformed LSMA 工具函数 *********************/
function listToColArray(list) {
  list = ee.List(list);
  return ee.Array(list.map(function(v) { return [v]; }));
}

function meanVectorFC(fc, bandList, fallbackImg, fbGeom) {
  fc = ee.FeatureCollection(fc);
  bandList = ee.List(bandList);
  var count = fc.size();

  var vals = bandList.map(function(b) {
    b = ee.String(b);
    return ee.Algorithms.If(
      count.gt(0),
      fc.aggregate_mean(b),
      fallbackImg.select([b]).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: fbGeom,
        scale: 10,
        maxPixels: 1e10,
        bestEffort: true
      }).get(b)
    );
  });
  return listToColArray(vals);
}

function safeCov(fc, bandList) {
  fc = ee.FeatureCollection(fc);
  bandList = ee.List(bandList);

  var n = ee.Number(fc.size());
  var p = bandList.size();
  var eye = ee.Array.identity(p);

  return ee.Array(ee.Algorithms.If(
    n.gt(1),
    ee.Array(
      bandList.map(function(bi) {
        bi = ee.String(bi);
        var mui = ee.Number(fc.aggregate_mean(bi));

        return bandList.map(function(bj) {
          bj = ee.String(bj);
          var muj = ee.Number(fc.aggregate_mean(bj));

          var withCovTerm = fc.map(function(f) {
            var xi = ee.Number(f.get(bi));
            var xj = ee.Number(f.get(bj));
            var term = xi.subtract(mui).multiply(xj.subtract(muj));
            return f.set('cov_term_tmp', term);
          });

          var s = ee.Number(withCovTerm.aggregate_sum('cov_term_tmp'));
          return s.divide(n.subtract(1));
        });
      })
    ),
    eye.multiply(1e-6)
  ));
}

function buildFisherTransform(pureSampledFC, bandList, fallbackImg, fbGeom) {
  var p = ee.List(bandList).size();
  var I = ee.Array.identity(p);

  var c0 = pureSampledFC.filter(ee.Filter.eq('cls', 0));
  var c1 = pureSampledFC.filter(ee.Filter.eq('cls', 1));
  var c2 = pureSampledFC.filter(ee.Filter.eq('cls', 2));
  var c3 = pureSampledFC.filter(ee.Filter.eq('cls', 3));

  var n0 = ee.Number(c0.size());
  var n1 = ee.Number(c1.size());
  var n2 = ee.Number(c2.size());
  var n3 = ee.Number(c3.size());

  var mu0 = meanVectorFC(c0, bandList, fallbackImg, fbGeom);
  var mu1 = meanVectorFC(c1, bandList, fallbackImg, fbGeom);
  var mu2 = meanVectorFC(c2, bandList, fallbackImg, fbGeom);
  var mu3 = meanVectorFC(c3, bandList, fallbackImg, fbGeom);

  var Sw = safeCov(c0, bandList).multiply(n0.subtract(1).max(0))
    .add(safeCov(c1, bandList).multiply(n1.subtract(1).max(0)))
    .add(safeCov(c2, bandList).multiply(n2.subtract(1).max(0)))
    .add(safeCov(c3, bandList).multiply(n3.subtract(1).max(0)))
    .add(I.multiply(LAMBDA));

  var invSw = Sw.matrixInverse();

  var w0 = invSw.matrixMultiply(mu0.subtract(mu3));
  var w1 = invSw.matrixMultiply(mu1.subtract(mu3));
  var w2 = invSw.matrixMultiply(mu2.subtract(mu3));

  var W = ee.Array.cat([
    w0.matrixTranspose(),
    w1.matrixTranspose(),
    w2.matrixTranspose()
  ], 0);

  return W;
}

function projectImageToFisher(img, W, bandList) {
  var arr = img.select(bandList).toArray().toArray(1);
  return ee.Image.constant(W).matrixMultiply(arr)
    .arrayProject([0])
    .arrayFlatten([['F1', 'F2', 'F3']]);
}

function projectFCToFisher(fc, bandList, W) {
  fc = ee.FeatureCollection(fc);
  bandList = ee.List(bandList);

  return fc.map(function(f) {
    var x = listToColArray(bandList.map(function(b) {
      return f.get(ee.String(b));
    }));
    var z = W.matrixMultiply(x);
    return f.set({
      F1: z.get([0, 0]),
      F2: z.get([1, 0]),
      F3: z.get([2, 0])
    });
  });
}

function projectVectorToFisher(vecCol, W) {
  return W.matrixMultiply(vecCol);
}

function colArrayToList3(colArr) {
  return ee.List([
    colArr.get([0, 0]),
    colArr.get([1, 0]),
    colArr.get([2, 0])
  ]);
}

function refineEndmemberByKNN_FSR(sampledFC, fallbackVec, bandList, W, fallbackImg, fbGeom) {
  sampledFC = ee.FeatureCollection(sampledFC);

  var proj = projectFCToFisher(sampledFC, bandList, W);

  var cF1 = ee.Number(proj.aggregate_mean('F1'));
  var cF2 = ee.Number(proj.aggregate_mean('F2'));
  var cF3 = ee.Number(proj.aggregate_mean('F3'));

  var withDist = proj.map(function(f) {
    var d = ee.Number(f.get('F1')).subtract(cF1).pow(2)
      .add(ee.Number(f.get('F2')).subtract(cF2).pow(2))
      .add(ee.Number(f.get('F3')).subtract(cF3).pow(2))
      .sqrt();
    return f.set('fDist', d);
  });

  var fsr = ee.Number(
    withDist.reduceColumns(ee.Reducer.percentile([FSR_PCTL]), ['fDist'])
      .get('p' + FSR_PCTL)
  );

  var kept = withDist
    .filter(ee.Filter.lte('fDist', fsr))
    .sort('fDist')
    .limit(KNN_K);

  return ee.Array(ee.Algorithms.If(
    kept.size().gt(0),
    meanVectorFC(kept, bandList, fallbackImg, fbGeom),
    fallbackVec
  ));
}

function normAbund4(img4) {
  var sum = img4.reduce(ee.Reducer.sum()).rename('UMX_sum');
  var nrm = img4.divide(sum).where(sum.eq(0), 0)
    .rename(['UMXn_Tree', 'UMXn_Grass', 'UMXn_Imperv', 'UMXn_Soil'])
    .float();
  return {nrm: nrm, sum: sum.float()};
}

function reconRmseFT(endmemberList, fisherImg, abundImg4) {
  var E = ee.Array(endmemberList).matrixTranspose();
  var Eimg = ee.Image.constant(E);

  var obs = fisherImg.toArray().toArray(1);
  var ab = abundImg4.toArray().toArray(1);

  var recon = Eimg.matrixMultiply(ab);
  var resid = obs.subtract(recon);

  return resid.multiply(resid)
    .arrayReduce(ee.Reducer.mean(), [0])
    .arrayProject([0])
    .arrayFlatten([['UMX_RMSE']])
    .sqrt()
    .float();
}


function collapseTGIS4toTGO3(abund4) {
  return ee.Image.cat([
    abund4.select('UMXn_Tree').rename('Soft_Tree'),
    abund4.select('UMXn_Grass').rename('Soft_Grass'),
    abund4.select('UMXn_Imperv').add(abund4.select('UMXn_Soil')).rename('Soft_Other')
  ]).float();
}

function toRF3Labels(fc) {
  return ee.FeatureCollection(fc).map(function(f) {
    var t = ee.Number(ee.Algorithms.If(f.get('TreeRatio'), f.get('TreeRatio'), 0));
    var g = ee.Number(ee.Algorithms.If(f.get('GrassRatio'), f.get('GrassRatio'), 0));
    var i = ee.Number(ee.Algorithms.If(f.get('ImpervRatio'), f.get('ImpervRatio'), 0));
    var s = ee.Number(ee.Algorithms.If(f.get('SoilRatio'), f.get('SoilRatio'), 0));

    var otherRaw = ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get('OtherRatio'), null),
      i.add(s),
      f.get('OtherRatio')
    );
    var o = ee.Number(otherRaw);

    return f.set({
      TreeRatio: t,
      GrassRatio: g,
      OtherRatio: o
    });
  });
}

function assignSoftLabelsToMixed3(samplesFC, abund4img) {
  var soft3img = collapseTGIS4toTGO3(abund4img);
  var fc3 = toRF3Labels(samplesFC);

  var tagged = fc3.map(function(f) {
    var t = ee.Number(ee.Algorithms.If(f.get('TreeRatio'), f.get('TreeRatio'), 0));
    var g = ee.Number(ee.Algorithms.If(f.get('GrassRatio'), f.get('GrassRatio'), 0));
    var o = ee.Number(ee.Algorithms.If(f.get('OtherRatio'), f.get('OtherRatio'), 0));
    var sum = t.add(g).add(o);

    var isPure = ee.Algorithms.If(
      t.eq(1),
      true,
      ee.Algorithms.If(g.eq(1), true, o.eq(1))
    );

    var isPure3 = ee.Algorithms.If(
      isPure,
      ee.Algorithms.If(sum.eq(1), 1, 0),
      0
    );

    return f.set('isPure3', isPure3);
  });

  var pure = tagged.filter(ee.Filter.eq('isPure3', 1));
  var mixed = tagged.filter(ee.Filter.eq('isPure3', 0));

  var mixedSoft = soft3img.sampleRegions({
    collection: mixed,
    scale: 10,
    geometries: true
  }).map(function(f) {
    return f.set({
      TreeRatio: f.get('Soft_Tree'),
      GrassRatio: f.get('Soft_Grass'),
      OtherRatio: f.get('Soft_Other')
    });
  });

  return pure.merge(mixedSoft);
}

function addHardClassFromRatios(fc) {
  return ee.FeatureCollection(fc).map(function(f) {
    var t = ee.Number(ee.Algorithms.If(f.get('TreeRatio'), f.get('TreeRatio'), 0));
    var g = ee.Number(ee.Algorithms.If(f.get('GrassRatio'), f.get('GrassRatio'), 0));
    var o = ee.Number(ee.Algorithms.If(f.get('OtherRatio'), f.get('OtherRatio'), 0));

    var maxV = ee.Number(ee.List([t, g, o]).reduce(ee.Reducer.max()));

    var cls = ee.Number(ee.Algorithms.If(
      t.eq(maxV), 0,
      ee.Algorithms.If(g.eq(maxV), 1, 2)
    ));

    return f.set({
      cls: cls,
      label_conf: maxV
    });
  });
}

/********************* 导出 *********************/
function exportToDrive(image, description, fileNamePrefix, region) {
  var out = image.clamp(0, 1).toFloat();
  Export.image.toDrive({
    image: out,
    description: description,
    folder: DRIVE_FOLDER,
    fileNamePrefix: fileNamePrefix,
    region: region,
    scale: 10,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}

/********************* 主流程 *********************/
function runAnalysis(region, cityName, endmemberSamples4, rfSamplesRaw3orMixed, trainBufferMeters) {
  var TRAIN_BUF_LOCAL = trainBufferMeters || TRAIN_BUF;
  var trainAOI = rfSamplesRaw3orMixed.geometry().buffer(TRAIN_BUF_LOCAL);

  var S2_BANDS_ALL = [
    'B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
    'NDVI','EVI','GCI','SAVI','NDBI','NDWI'
  ];

  var statsReducer = ee.Reducer.mean()
    .combine(ee.Reducer.min(), '', true)
    .combine(ee.Reducer.max(), '', true)
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.percentile([25,50,75], ['p25','p50','p75']), '', true);

  var s2_city = s2Collection(region, YEAR, MONTHS);
  var s2idx_city = s2_city.map(addIndices);
  var base_city = s2idx_city.select(S2_BANDS_ALL).reduce(statsReducer).clip(region).unmask(-9999);
  var s2core_city = s2_city.select(UNMIX_BANDS).median().clip(region).unmask(-9999);

  var s2_train = s2Collection(trainAOI, YEAR, MONTHS);
  var s2idx_train = s2_train.map(addIndices);
  var base_train = s2idx_train.select(S2_BANDS_ALL).reduce(statsReducer).clip(trainAOI).unmask(-9999);
  var s2core_train = s2_train.select(UNMIX_BANDS).median().clip(trainAOI).unmask(-9999);

  var s1_city = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(region)
    .filterDate(
      ee.Date.fromYMD(YEAR, MONTHS[0], 1),
      ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1, 'month')
    )
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
  var vv_city = s1_city.select('VV').median().rename('VV_median');
  var vh_city = s1_city.select('VH').median().rename('VH_median');

  var s1_train = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(trainAOI)
    .filterDate(
      ee.Date.fromYMD(YEAR, MONTHS[0], 1),
      ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1, 'month')
    )
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
  var vv_train = s1_train.select('VV').median().rename('VV_median');
  var vh_train = s1_train.select('VH').median().rename('VH_median');

  var dem = ee.Image('USGS/SRTMGL1_003').unmask(0);
  var elevation = dem.rename('elevation');
  var slope = ee.Terrain.slope(dem).rename('slope');
  var aspect = ee.Terrain.aspect(dem).rename('aspect');

  var samplesCity = endmemberSamples4.filterBounds(region);

  var treePureCity = samplesCity.filter(ee.Filter.eq('TreeRatio', 1));
  var grassPureCity = samplesCity.filter(ee.Filter.eq('GrassRatio', 1));
  var impervPureCity = samplesCity.filter(ee.Filter.eq('ImpervRatio', 1));
  var soilPureCity = samplesCity.filter(ee.Filter.eq('SoilRatio', 1));

  function sampleForEM(imgCore, fc, cls) {
    return imgCore.sampleRegions({
      collection: fc,
      scale: 10,
      geometries: false
    })
    .filter(ee.Filter.neq(UNMIX_BANDS[0], -9999))
    .map(function(f) { return f.set('cls', cls); });
  }

  var tS0 = sampleForEM(s2core_city, treePureCity, 0);
  var gS0 = sampleForEM(s2core_city, grassPureCity, 1);
  var iS0 = sampleForEM(s2core_city, impervPureCity, 2);
  var sS0 = sampleForEM(s2core_city, soilPureCity, 3);

  var pureSampled0 = tS0.merge(gS0).merge(iS0).merge(sS0);

  var tEnd0 = meanVectorFC(tS0, UNMIX_BANDS, s2core_city, region);
  var gEnd0 = meanVectorFC(gS0, UNMIX_BANDS, s2core_city, region);
  var iEnd0 = meanVectorFC(iS0, UNMIX_BANDS, s2core_city, region);
  var sEnd0 = meanVectorFC(sS0, UNMIX_BANDS, s2core_city, region);

  var Wf = buildFisherTransform(pureSampled0, UNMIX_BANDS, s2core_city, region);

  var tEnd = refineEndmemberByKNN_FSR(tS0, tEnd0, UNMIX_BANDS, Wf, s2core_city, region);
  var gEnd = refineEndmemberByKNN_FSR(gS0, gEnd0, UNMIX_BANDS, Wf, s2core_city, region);
  var iEnd = refineEndmemberByKNN_FSR(iS0, iEnd0, UNMIX_BANDS, Wf, s2core_city, region);
  var sEnd = refineEndmemberByKNN_FSR(sS0, sEnd0, UNMIX_BANDS, Wf, s2core_city, region);

  var tEndF = projectVectorToFisher(tEnd, Wf);
  var gEndF = projectVectorToFisher(gEnd, Wf);
  var iEndF = projectVectorToFisher(iEnd, Wf);
  var sEndF = projectVectorToFisher(sEnd, Wf);

  var endmembersFisher = ee.List([
    colArrayToList3(tEndF),
    colArrayToList3(gEndF),
    colArrayToList3(iEndF),
    colArrayToList3(sEndF)
  ]);

  var fisher_city = projectImageToFisher(s2core_city, Wf, UNMIX_BANDS);
  var fisher_train = projectImageToFisher(s2core_train, Wf, UNMIX_BANDS);

  var umx_city_img = fisher_city.unmix(endmembersFisher, true, true)
    .rename(['UMX_Tree','UMX_Grass','UMX_Imperv','UMX_Soil'])
    .clamp(0, 1);

  var umx_train_img = fisher_train.unmix(endmembersFisher, true, true)
    .rename(['UMX_Tree','UMX_Grass','UMX_Imperv','UMX_Soil'])
    .clamp(0, 1);

  var nCity = normAbund4(umx_city_img);
  var nTrain = normAbund4(umx_train_img);

  var rmse_city = reconRmseFT(endmembersFisher, fisher_city, umx_city_img);
  var rmse_train = reconRmseFT(endmembersFisher, fisher_train, umx_train_img);

  var soft3_city = collapseTGIS4toTGO3(nCity.nrm);
  var soft3_train = collapseTGIS4toTGO3(nTrain.nrm);

  var rfSamplesAll3 = assignSoftLabelsToMixed3(rfSamplesRaw3orMixed, nCity.nrm);

  var feature_city_base = base_city.addBands([vh_city, vv_city, aspect, elevation, slope]).clip(region);
  var feature_train_base = base_train.addBands([vh_train, vv_train, aspect, elevation, slope]).clip(trainAOI);

  var feature_city = feature_city_base
    .addBands([nCity.nrm, soft3_city, nCity.sum, rmse_city])
    .clip(region);

  var feature_train = feature_train_base
    .addBands([nTrain.nrm, soft3_train, nTrain.sum, rmse_train])
    .clip(trainAOI);

  var allBands = feature_city.bandNames();
  var classNames = ['Tree', 'Grass', 'Other'];

  var rfSamplesAll3Cls = addHardClassFromRatios(rfSamplesAll3);
  rfSamplesAll3Cls = rfSamplesAll3Cls.filter(ee.Filter.gte('label_conf', 0.50));

  var rfSamples = feature_train.sampleRegions({
      collection: rfSamplesAll3Cls,
      properties: ['cls', 'label_conf'],
      scale: 10,
      geometries: false
    })
    .filter(ee.Filter.neq(KEY_BAND, -9999))
    .filter(ee.Filter.gt('UMXn_Tree', -0.5));

  var rfMulti = ee.Classifier.smileRandomForest(200)
    .setOutputMode('MULTIPROBABILITY')
    .train({
      features: rfSamples,
      classProperty: 'cls',
      inputProperties: allBands
    });

  var probArr = feature_city.classify(rfMulti);
  var probs1 = probArr.arrayFlatten([classNames]).float();

  var sumP = probs1.reduce(ee.Reducer.sum()).rename('sumP');
  probs1 = probs1.divide(sumP).where(sumP.eq(0), 0).float();

  exportToDrive(probs1.select('Tree'), cityName + '_Tree', cityName + '_Tree', region);
  exportToDrive(probs1.select('Grass'), cityName + '_Grass', cityName + '_Grass', region);
}

/********************* UI控件 *********************/
var CURRENT_CITY_NAME = null;
var CURRENT_REGION = null;

var uiPanel = ui.Panel({
  style: {width: '430px', position: 'top-left', padding: '8px'}
});



var citySelect = ui.Select({placeholder: '加载城市列表中…'});
allCitiesFC.aggregate_array('city').distinct().sort().evaluate(function(list) {
  citySelect.items().reset(list);
  citySelect.setPlaceholder('选择城市');
});

var yearBox = ui.Textbox({
  placeholder: '年份（默认 ' + YEAR + '）',
  value: ''
});

var confirmBtn = ui.Button({
  label: '✅ 确认城市',
  style: {stretch: 'horizontal'}
});

var runBtn = ui.Button({
  label: '▶️ 运行整年分析',
  style: {stretch: 'horizontal'}
});

uiPanel.add(ui.Panel([ui.Label('城市'), citySelect], ui.Panel.Layout.flow('horizontal')));
uiPanel.add(ui.Panel([ui.Label('年份'), yearBox], ui.Panel.Layout.flow('horizontal')));
uiPanel.add(confirmBtn);
uiPanel.add(runBtn);
Map.add(uiPanel);

confirmBtn.onClick(function() {
  var cityName = citySelect.getValue();
  if (!cityName) return;

  var cityFC = allCitiesFC.filter(ee.Filter.eq('city', cityName));
  var region = cityFC.union(1).geometry();

  Map.clear();
  Map.add(uiPanel);
  Map.centerObject(region, 10);

  CURRENT_CITY_NAME = cityName;
  CURRENT_REGION = region;

  var yrTxt = yearBox.getValue();
  YEAR = (yrTxt && yrTxt.trim() !== '') ? parseInt(yrTxt, 10) : YEAR;
});

runBtn.onClick(function() {
  var cityName = CURRENT_CITY_NAME || citySelect.getValue();
  if (!cityName) return;

  var region = CURRENT_REGION;
  if (!region) {
    var cityFC = allCitiesFC.filter(ee.Filter.eq('city', cityName));
    region = cityFC.union(1).geometry();
  }

  var yrTxt = yearBox.getValue();
  YEAR = (yrTxt && yrTxt.trim() !== '') ? parseInt(yrTxt, 10) : YEAR;
  MONTHS = [1, 12];

  Map.clear();
  Map.add(uiPanel);
  Map.centerObject(region, 10);

  runAllForUI(region, cityName);
});

/********************* UI 封装运行 *********************/
function runAllForUI(region, cityName) {
  var auto4 = buildAutoSamplesTGIS(region, cityName);

  var autoRF3 = auto4.map(function(f) {
    return f.set({
      TreeRatio: ee.Number(f.get('RF_TreeRatio')),
      GrassRatio: ee.Number(f.get('RF_GrassRatio')),
      OtherRatio: ee.Number(f.get('RF_OtherRatio'))
    });
  });

  var manualLoad = tryLoadManualSamples(MANUAL_SHP_ASSET);
  var manualRF3 = ee.FeatureCollection([]);

  if (manualLoad.has) {
    manualRF3 = ensureRF3OneHot(manualLoad.fc).filterBounds(region);
  }

  var rfSamplesAll = ee.FeatureCollection(
    ee.Algorithms.If(
      manualLoad.has,
      autoRF3.merge(manualRF3),
      autoRF3
    )
  ).distinct(['.geo']);

  runAnalysis(region, cityName, auto4, rfSamplesAll, TRAIN_BUF);
}
