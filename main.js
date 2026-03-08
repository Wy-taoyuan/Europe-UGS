
/********************* 用户配置区 *********************/
// 城市边界（FeatureCollection 或单 Feature 的资产路径）
var CITY_ASSET = 'users/wangwenna111/Europe/Germany_Oberbayern';

// 可选：手动绘制样本（SHP转资产后的路径；字段与自动样本一致，或至少包含 class=1/2/3）
// 若没有手动样本，请留空字符串 ''
var MANUAL_SHP_ASSET = 'projects/total-contact-445408-m8/assets/Germany_Oberbayern_self'; // e.g. 'users/xxx/manual_samples_cityX'

// 文件导出设置（可选）
var DRIVE_FOLDER = 'GEE_exports';
var FILE_PREFIX  = 'Germany_Oberbayern';

// 时间与数据源参数
var YEAR = 2024;
var MONTHS = [1, 12];     // [起始月, 结束月]，闭区间
var CLOUD_PCT = 20;
var SCALE = 10;           // 采样/计算尺度（m）

// 自动样本参数
var N_TREE  = 1000;
var N_GRASS = 1000;
var N_OTHER = 1000;
var ERODE_PIX = 3;        // 腐蚀像素（抠掉边界混合像元）

// 解混 & RF 参数
var UNMIX_BANDS = ['B4','B8','B11'];
var LAMBDA = 1e-3;         // 岭回归系数（pinv 稳定）
var TRAIN_BUF = 10000;     // 训练 AOI = 全部样本的缓冲（m）
var KEEP_FRAC = 0.5;       // 端元提纯时保留比例（Top-K占比）
var KEY_BAND = 'B4_mean';  // reduce 后将产生 B4_mean，用于过滤无效

/********************* 通用小工具 *********************/
function s2Collection(region, year, months){
  var t0 = ee.Date.fromYMD(year, months[0], 1);
  var t1 = ee.Date.fromYMD(year, months[1], 1).advance(1, 'month');
  return ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(region)
    .filterDate(t0, t1)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PCT));
}

// 统一指数：兼容自动样本与特征工程
function addIndices(img){
  var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
  var evi  = img.expression('2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))', {
    NIR: img.select('B8'), RED: img.select('B4'), BLUE: img.select('B2')
  }).rename('EVI');
  var gci  = img.expression('(NIR/GREEN)-1', { NIR: img.select('B8'), GREEN: img.select('B3') }).rename('GCI');
  var savi = img.expression('((NIR-RED)/(NIR+RED+0.5))*1.5', { NIR: img.select('B8'), RED: img.select('B4') }).rename('SAVI');
  var ndbi = img.normalizedDifference(['B11','B8']).rename('NDBI');
  var ndwi = img.normalizedDifference(['B3','B8']).rename('NDWI');
  var base = img.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12']);
  return base.addBands([ndvi,evi,gci,savi,ndbi,ndwi]);
}

function getSafeRegion(cityAssetPath){
  var fc = ee.FeatureCollection(cityAssetPath);
  var first = ee.Feature(fc.first());
  var geom = first.geometry();
  var region = ee.Geometry(ee.Algorithms.If(geom, geom, fc.geometry().dissolve()));
  return {fc: fc, region: region};
}

function maskCount(img, region, scale){
  var d = ee.Dictionary(
    img.mask().reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: region, scale: scale,
      maxPixels: 1e12, bestEffort: true
    })
  );
  var vals = d.values();
  var hasVal = ee.Number(vals.size()).gt(0);
  return ee.Number(ee.Algorithms.If(hasVal, vals.get(0), 0));
}

// 若手动样本只有 class（1/2/3），补齐 one-hot 标签
function ensureOneHot(fc){
  return ee.FeatureCollection(fc).map(function(f){
    var hasTree   = f.get('TreeRatio');
    var hasGrass  = f.get('GrassRatio');
    var hasOther  = f.get('OtherRatio');
    var needFill  = ee.Algorithms.IsEqual(hasTree, null)
      .or(ee.Algorithms.IsEqual(hasGrass, null))
      .or(ee.Algorithms.IsEqual(hasOther, null));
    var c = ee.Number(f.get('class'));
    var lab = ee.Dictionary(ee.Algorithms.If(needFill,
      ee.Algorithms.If(
        c.eq(1), ee.Dictionary({'TreeRatio':1,'GrassRatio':0,'OtherRatio':0}),
        ee.Algorithms.If(
          c.eq(2), ee.Dictionary({'TreeRatio':0,'GrassRatio':1,'OtherRatio':0}),
          ee.Dictionary({'TreeRatio':0,'GrassRatio':0,'OtherRatio':1})
        )
      ),
      ee.Dictionary({})
    ));
    return f.set(lab);
  });
}

/********************* 自动样本：高纯净种子（Other = 非树非草） *********************/
function seedsFromWorldCoverRefined(region) {
  var wc = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(region);
  var s2Med = s2Collection(region, YEAR, MONTHS).map(addIndices).median();
  var ndvi = s2Med.select('NDVI');
  var ndwi = s2Med.select('NDWI');

  var tree = wc.eq(10);
  var grass = wc.eq(30);
  var other = wc.neq(10).and(wc.neq(30));  // 非树非草

  tree = tree.focal_min({radius: ERODE_PIX, units:'pixels'})
            .updateMask(ndvi.gte(0.5).and(ndwi.lt(0.2)));
  grass = grass.focal_min({radius: ERODE_PIX, units:'pixels'})
              .updateMask(ndvi.gte(0.3).and(ndvi.lt(0.6)).and(ndwi.lt(0.2)));
  other = other.focal_min({radius: ERODE_PIX, units:'pixels'})
              .updateMask(ndvi.lt(0.2));  // 排除潜在绿地

  // 动态地表：进一步过滤已建区/道路
  var dw = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
              .filterDate(ee.Date.fromYMD(YEAR, MONTHS[0], 1),
                          ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1,'month'))
              .select('label')
              .mode()
              .clip(region);
  var built = dw.eq(6).or(dw.eq(7));  // built + road

  tree  = tree.updateMask(built.not()).selfMask();
  grass = grass.updateMask(built.not()).selfMask();
  other = other.updateMask(ee.Image(1)).selfMask(); // 其他保留
  return {tree: tree, grass: grass, other: other};
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

/********************* 尝试加载手动样本（若存在） *********************/
function tryLoadManualSamples(path){
  if (!path || path === '') return {has: false, fc: ee.FeatureCollection([])};
  var has = false, fc = ee.FeatureCollection([]);
  try {
    fc = ee.FeatureCollection(path);
    var c = fc.limit(1).size().getInfo(); // 触发检查
    has = (c >= 0);
  } catch (err) {
    print('⚠️ 手动样本资产不可用/不存在：', path);
    has = false;
    fc = ee.FeatureCollection([]);
  }
  // 不再 ensureOneHot，直接返回
  return {has: has, fc: fc};
}


/********************* —— 解混 + RF 主流程（接收 samplesAll） *********************/
function runAnalysis(region, cityName, samplesAll, trainBufferMeters, keepFraction) {

  var TRAIN_BUF_LOCAL = trainBufferMeters || TRAIN_BUF;
  var KEEP_FRAC_LOCAL = ee.Number(keepFraction || KEEP_FRAC);

  function meanVectorOrFallback(sampled, bandList, fallbackImg, fbGeom){
    bandList = ee.List(bandList);
    var count = sampled.size();
    var vec = bandList.map(function(b){
      var v = sampled.aggregate_mean(ee.String(b));
      return ee.Algorithms.If(v, v, null);
    });
    var meanArr = ee.Array(vec);
    var fbArr = ee.Array(bandList.map(function(b){
      return fallbackImg.select([ee.String(b)]).reduceRegion({
        reducer: ee.Reducer.percentile([10,90]),
        geometry: fbGeom, scale: 10, maxPixels: 1e10
      }).values().reduce(ee.Reducer.mean());
    }));
    return ee.Array(ee.Algorithms.If(count.gt(0), meanArr, fbArr));
  }

  function selectTopFrac(fcSampled, band, keepFrac){
    var n = ee.Number(fcSampled.size());
    var k = n.multiply(keepFrac).ceil();
    return ee.FeatureCollection(ee.Algorithms.If(
      n.gt(0),
      fcSampled.sort(band, false).limit(k),
      fcSampled
    ));
  }

  function normAbund(img3) {
    var sum = img3.reduce(ee.Reducer.sum()).rename('UMX_sum');
    var nrm = img3.divide(sum).where(sum.eq(0), 0)
      .rename(['UMXn_Tree','UMXn_Grass','UMXn_Other']).float();
    return {nrm: nrm, sum: sum.float()};
  }

  function reconRmse(Eimg, observedColVec, abundImg3) {
    var abundCol = abundImg3.toArray().toArray(1); // 3x1
    var recon    = Eimg.matrixMultiply(abundCol);   // 3x1
    var resid    = observedColVec.subtract(recon);  // 3x1
    var rmse = resid.multiply(resid)
      .arrayReduce(ee.Reducer.mean(), [0])
      .arrayProject([0])
      .arrayFlatten([['UMX_RMSE']])
      .sqrt().float();
    return rmse;
  }

  /**** ---------- 样本：城内 for 端元；全域 for 训练 ---------- ****/
  var samplesCity = samplesAll.filterBounds(region);

  var treePureCity  = samplesCity.filter(ee.Filter.eq('TreeRatio', 1));
  var grassPureCity = samplesCity.filter(ee.Filter.eq('GrassRatio', 1));
  var otherPureCity = samplesCity.filter(ee.Filter.eq('OtherRatio', 1));

  var treePureAll  = samplesAll.filter(ee.Filter.eq('TreeRatio', 1))
                    .map(function(f){ return f.set({'TreeRatio':1, 'GrassRatio':0, 'OtherRatio':0}); });
  var grassPureAll = samplesAll.filter(ee.Filter.eq('GrassRatio', 1))
                    .map(function(f){ return f.set({'TreeRatio':0, 'GrassRatio':1, 'OtherRatio':0}); });
  var otherPureAll = samplesAll.filter(ee.Filter.eq('OtherRatio', 1))
                    .map(function(f){ return f.set({'TreeRatio':0, 'GrassRatio':0, 'OtherRatio':1}); });

  var pureAll = treePureAll.merge(grassPureAll).merge(otherPureAll);

  var cityPurePoints = ee.FeatureCollection([treePureCity, grassPureCity, otherPureCity])
    .flatten().map(function(f){ return f.setGeometry(f.geometry().centroid(10)); });

  /**** ---------- 训练AOI ---------- ****/
  var trainAOI = samplesAll.geometry().buffer(TRAIN_BUF_LOCAL);

  /**** ---------- 影像准备：城市区 & 训练AOI ---------- ****/
  var S2_BANDS_ALL = [
    'B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12',
    'NDVI','EVI','GCI','SAVI','NDBI','NDWI'
  ];

  var statsReducer = ee.Reducer.mean()
    .combine(ee.Reducer.min(),    '', true)
    .combine(ee.Reducer.max(),    '', true)
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.percentile([25, 50, 75], ['p25','p50','p75']), '', true);

  // 城市
  var s2_city = s2Collection(region, YEAR, MONTHS);
  var s2idx_city = s2_city.map(addIndices);
  var base_city = s2idx_city.select(S2_BANDS_ALL).reduce(statsReducer).clip(region).unmask(-9999);
  var s2core_city = s2_city.select(UNMIX_BANDS).median().clip(region).unmask(-9999);

  // 训练 AOI
  var s2_train = s2Collection(trainAOI, YEAR, MONTHS);
  var s2idx_train = s2_train.map(addIndices);
  var base_train = s2idx_train.select(S2_BANDS_ALL).reduce(statsReducer).clip(trainAOI).unmask(-9999);
  var s2core_train = s2_train.select(UNMIX_BANDS).median().clip(trainAOI).unmask(-9999);

  // S1
  var s1_city = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(region).filterDate(ee.Date.fromYMD(YEAR, MONTHS[0], 1),
                                    ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1,'month'))
    .filter(ee.Filter.eq('instrumentMode','IW'))
    .filter(ee.Filter.eq('orbitProperties_pass','DESCENDING'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'));
  var vv_city = s1_city.select('VV').median().rename('VV_median');
  var vh_city = s1_city.select('VH').median().rename('VH_median');

  var s1_train = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(trainAOI).filterDate(ee.Date.fromYMD(YEAR, MONTHS[0], 1),
                                      ee.Date.fromYMD(YEAR, MONTHS[1], 1).advance(1,'month'))
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

  /**** ---------- 端元：初解混 → Top-K 提纯 → 最终端元 ---------- ****/
  function sampleForEM(imgCore, fc){
    return imgCore.sampleRegions({collection: fc, scale:10, geometries:false})
                  .filter(ee.Filter.neq(UNMIX_BANDS[0], -9999));
  }
  var tS0 = sampleForEM(s2core_city, treePureCity);
  var gS0 = sampleForEM(s2core_city, grassPureCity);
  var oS0 = sampleForEM(s2core_city, otherPureCity);

  var tEnd0 = meanVectorOrFallback(tS0, UNMIX_BANDS, s2core_city, region);
  var gEnd0 = meanVectorOrFallback(gS0, UNMIX_BANDS, s2core_city, region);
  var oEnd0 = meanVectorOrFallback(oS0, UNMIX_BANDS, s2core_city, region);

  var E0img = ee.Image.constant(ee.Array.cat([tEnd0, gEnd0, oEnd0], 1));
  var ET0   = E0img.matrixTranspose();
  var A0    = ET0.matrixMultiply(E0img);
  var I3    = ee.Image.constant(ee.Array.identity(3));
  var pinv0 = A0.add(I3.multiply(LAMBDA)).matrixInverse().matrixMultiply(ET0);

  var observed_city0 = s2core_city.toArray().toArray(1);
  var umx_city_img0 = pinv0.matrixMultiply(observed_city0)
    .arrayProject([0]).arrayFlatten([['UMX_Tree','UMX_Grass','UMX_Other']]).clamp(0,1);

  var cityPurePoints = ee.FeatureCollection([treePureCity, grassPureCity, otherPureCity])
    .flatten().map(function(f){ return f.setGeometry(f.geometry().centroid(10)); });

  var pts_with_umx0 = umx_city_img0.sampleRegions({
    collection: cityPurePoints, scale: 10, geometries: true
  });

  var tPure_pts = selectTopFrac(pts_with_umx0, 'UMX_Tree',  KEEP_FRAC_LOCAL);
  var gPure_pts = selectTopFrac(pts_with_umx0, 'UMX_Grass', KEEP_FRAC_LOCAL);
  var oPure_pts = selectTopFrac(pts_with_umx0, 'UMX_Other', KEEP_FRAC_LOCAL);

  function keepStats(name, fcAll, fcKept){
    var nAll = ee.Number(fcAll.size());
    var nKept = ee.Number(fcKept.size());
    var frac = nKept.divide(nAll).multiply(100);
    print('🧼 提纯后端元点数（' + name + '）: kept=', nKept, ' / total=', nAll, ' (', frac.format('%.1f'), '%)');
  }
  keepStats('Tree',  pts_with_umx0, tPure_pts);
  keepStats('Grass', pts_with_umx0, gPure_pts);
  keepStats('Other', pts_with_umx0, oPure_pts);

  function endFromPts(imgCore, pts, bandList){
    var samp = imgCore.sampleRegions({collection: pts, scale: 10, geometries: false});
    return meanVectorOrFallback(samp, bandList, imgCore, region);
  }
  var tEnd = ee.Array(ee.Algorithms.If(tPure_pts.size().gt(0), endFromPts(s2core_city, tPure_pts, UNMIX_BANDS), tEnd0));
  var gEnd = ee.Array(ee.Algorithms.If(gPure_pts.size().gt(0), endFromPts(s2core_city, gPure_pts, UNMIX_BANDS), gEnd0));
  var oEnd = ee.Array(ee.Algorithms.If(oPure_pts.size().gt(0), endFromPts(s2core_city, oPure_pts, UNMIX_BANDS), oEnd0));

  var Eimg = ee.Image.constant(ee.Array.cat([tEnd, gEnd, oEnd], 1));
  var ET   = Eimg.matrixTranspose();
  var A    = ET.matrixMultiply(Eimg);
  var I    = ee.Image.constant(ee.Array.identity(3));
  var pinv = A.add(I.multiply(LAMBDA)).matrixInverse().matrixMultiply(ET);

  /**** ---------- UMX + 派生特征 ---------- ****/
  var observed_city = s2core_city.toArray().toArray(1);
  var umx_city_img = pinv.matrixMultiply(observed_city)
    .arrayProject([0]).arrayFlatten([['UMX_Tree','UMX_Grass','UMX_Other']]).clamp(0,1);
  var nCity = normAbund(umx_city_img);
  var rmse_city  = reconRmse(Eimg, observed_city, umx_city_img);

  var observed_train = s2core_train.toArray().toArray(1);
  var umx_train_img = pinv.matrixMultiply(observed_train)
    .arrayProject([0]).arrayFlatten([['UMX_Tree','UMX_Grass','UMX_Other']]).clamp(0,1);
  var nTrain = normAbund(umx_train_img);
  var rmse_train = reconRmse(Eimg, observed_train, umx_train_img);

  Map.addLayer(nCity.nrm.select('UMXn_Tree'),  {min:0, max:1, palette:['white','#FF0000']}, 'UMXn_Tree (city)');
  Map.addLayer(nCity.nrm.select('UMXn_Grass'), {min:0, max:1, palette:['white','#006837']}, 'UMXn_Grass (city)');

  var feature_city_base  = base_city.addBands([vh_city, vv_city, aspect, elevation, slope]).clip(region);
  var feature_train_base = base_train.addBands([vh_train, vv_train, aspect, elevation, slope]).clip(trainAOI);
  var feature_city  = feature_city_base.addBands([nCity.nrm,  nCity.sum,  rmse_city ]).clip(region);
  var feature_train = feature_train_base.addBands([nTrain.nrm, nTrain.sum, rmse_train]).clip(trainAOI);

  print('📑 feature_city band count:', feature_city.bandNames().size());
  //print('📑 feature_city bands (first 40):', feature_city.bandNames().slice(0,40));

  /**** ---------- RF（单模型多概率） ---------- ****/
  var allBands = feature_city.bandNames();

  var pureAll_cls = pureAll.map(function(f){
    var cls = ee.Algorithms.If(ee.Number(f.get('TreeRatio')).eq(1), 0,
              ee.Algorithms.If(ee.Number(f.get('GrassRatio')).eq(1), 1, 2));
    return f.set('cls', cls);
  });
  var classNames = ee.List(['Tree','Grass','Other']);

  var trainCls = feature_train.sampleRegions({
      collection: pureAll_cls, properties: ['cls'], scale: 10, geometries: false
    })
    .filter(ee.Filter.neq(KEY_BAND, -9999))
    .filter(ee.Filter.gt('UMXn_Tree', -0.5));

  //print('📦 RF(多概率) 训练样本数：', trainCls.size());

  var rfMulti = ee.Classifier.smileRandomForest(200)
    .setOutputMode('MULTIPROBABILITY')
    .train({ features: trainCls, classProperty: 'cls', inputProperties: allBands });

  var probArr = feature_city.classify(rfMulti);
  var probs   = probArr.arrayFlatten([classNames]);
  var sumP    = probs.reduce(ee.Reducer.sum()).rename('sumP');
  var probs1  = probs.divide(sumP).where(sumP.eq(0), 0).float();

  var visTree  = {min:0, max:1, palette:['white','#FF0000']};
  var visGrass = {min:0, max:1, palette:['white','#006837']};
  var visOther = {min:0, max:1, palette:['white','#000000']};
  Map.addLayer(probs1.select('Tree'),  visTree,  'Tree (RF multi-prob)');
  Map.addLayer(probs1.select('Grass'), visGrass, 'Grass (RF multi-prob)');
  Map.addLayer(probs1.select('Other'), visOther, 'Other (RF multi-prob)');

  // var sumCheck = probs1.reduce(ee.Reducer.sum()).rename('Sum(Tree+Grass+Other)');
  // print('🔎 Sum(Tree+Grass+Other) min/max：', sumCheck.reduceRegion({
  //   reducer: ee.Reducer.minMax(), geometry: region, scale: 30, maxPixels: 1e9
  // }));

  var rgb = ee.Image.cat([probs1.select('Other'), probs1.select('Tree'), probs1.select('Grass')]);
  Map.addLayer(rgb, {min:[0,0,0], max:[1,1,1]}, cityName + ' RGB (RF multi-prob)');

  // 可选阈值硬化与导出
  var HARDEN_LOW  = 0.10;
  var HARDEN_HIGH = 0.90;
  function harden01(img){
    img = img.float();
    var out = img.where(img.gte(HARDEN_HIGH), 1);
    out = out.where(img.lte(HARDEN_LOW), 0);
    return out;
  }
  var probs_thr = ee.Image.cat([
    harden01(probs1.select('Tree')).rename('Tree'),
    harden01(probs1.select('Grass')).rename('Grass'),
    harden01(probs1.select('Other')).rename('Other')
  ]).float();

  Map.addLayer(probs_thr.select('Tree'),  visTree,  'Tree (thr, RF multi-prob)');
  Map.addLayer(probs_thr.select('Grass'), visGrass, 'Grass (thr, RF multi-prob)');
  Map.addLayer(probs_thr.select('Other'), visOther, 'Other (thr, RF multi-prob)');

  function exportToDrive(image, description, fileNamePrefix) {
    Export.image.toDrive({
      image: image,
      description: description,
      folder: 'GEE_Unmix',
      fileNamePrefix: fileNamePrefix,
      region: region,
      scale: 10,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF',
      formatOptions: {cloudOptimized: true}
    });
  }
  exportToDrive(probs_thr.select('Tree'),  cityName + '_Tree_rf_mprob_thr_f32',  cityName + '_Tree_rf_mprob_thr_f32');
  exportToDrive(probs_thr.select('Grass'), cityName + '_Grass_rf_mprob_thr_f32', cityName + '_Grass_rf_mprob_thr_f32');
  exportToDrive(probs_thr.select('Other'), cityName + '_Other_rf_mprob_thr_f32', cityName + '_Other_rf_mprob_thr_f32');
}

/********************* 主控：自动 + 手动合并 → 解混 + RF *********************/
function runAll(){
  // 区域与命名
  var SAFE = getSafeRegion(CITY_ASSET);
  var cityFC = SAFE.fc;
  var region = SAFE.region;

  // cityName 用 FC 的 name 字段，否则用 FILE_PREFIX（需 client 字符串）
  var cityName = ee.String(ee.Algorithms.If(
    cityFC.first().get('name'),
    cityFC.first().get('name'),
    FILE_PREFIX
  )).getInfo();

  // 视图：用 cityFC 来 style（Geometry 不支持 .style()）
  Map.centerObject(cityFC, 10);
  Map.addLayer(
    cityFC.style({color: 'black', fillColor: '00000000', width: 2}),
    {},
    cityName + ' 边界'
  );

  // 自动样本生成
  var seeds = seedsFromWorldCoverRefined(region);
  var treePts  = stratifiedPointsSafe(seeds.tree,  1, N_TREE,  region, 'Tree');
  var grassPts = stratifiedPointsSafe(seeds.grass, 2, N_GRASS, region, 'Grass');
  var otherPts = stratifiedPointsSafe(seeds.other, 3, N_OTHER, region, 'Other');

  var autoAll = ee.FeatureCollection(treePts).merge(grassPts).merge(otherPts)
    .map(function(f){
      var c = ee.Number(f.get('class'));
      var lab = ee.Dictionary(ee.Algorithms.If(
        c.eq(1), ee.Dictionary({'TreeRatio':1,'GrassRatio':0,'OtherRatio':0}),
        ee.Algorithms.If(
          c.eq(2), ee.Dictionary({'TreeRatio':0,'GrassRatio':1,'OtherRatio':0}),
          ee.Dictionary({'TreeRatio':0,'GrassRatio':0,'OtherRatio':1})
        )
      ));
      return f.set(lab).set({'city': cityName});
    })
    .distinct(['.geo']);

  // 手动样本（如存在）合并
  var manualLoad = tryLoadManualSamples(MANUAL_SHP_ASSET);
  var samplesAll = (manualLoad.has ? autoAll.merge(manualLoad.fc) : autoAll).distinct(['.geo']);

  // 可视化样本
  Map.addLayer(samplesAll.filter(ee.Filter.eq('TreeRatio',1)),  {color:'#00b050'}, 'Tree (samples)');
  Map.addLayer(samplesAll.filter(ee.Filter.eq('GrassRatio',1)), {color:'#92d050'}, 'Grass (samples)');
  Map.addLayer(samplesAll.filter(ee.Filter.eq('OtherRatio',1)), {color:'#ff0000'}, 'Other (samples)');

  // 统计并打印
  var nTree  = samplesAll.filter(ee.Filter.eq('TreeRatio',1)).size();
  var nGrass = samplesAll.filter(ee.Filter.eq('GrassRatio',1)).size();
  var nOther = samplesAll.filter(ee.Filter.eq('OtherRatio',1)).size();
  var nAll   = samplesAll.size();
  print('样本点数量（Tree / Grass / Other / All）:', nTree, nGrass, nOther, nAll);
  if (manualLoad.has) print('✅ 已合并手动样本：', MANUAL_SHP_ASSET);
  else print('ℹ️ 未提供/未加载手动样本，仅使用自动样本。');

  // （可选）导出样本
  Export.table.toDrive({
    collection: ee.FeatureCollection(ee.Algorithms.If(nAll.gt(0), samplesAll, ee.FeatureCollection([]))),
    description: cityName + '_samples_auto_manual',
    folder: DRIVE_FOLDER,
    fileNamePrefix: FILE_PREFIX + '_samples',
    fileFormat: 'GeoJSON'
  });

  // 进入解混 + RF
  runAnalysis(region, cityName, samplesAll, TRAIN_BUF, KEEP_FRAC);
}

// 运行
runAll();

