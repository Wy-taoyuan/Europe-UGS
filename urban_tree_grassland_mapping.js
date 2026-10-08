
/***************************************************************
 CONFIGURATION AND CONSTANTS
 ***************************************************************/

var MANUAL_SHP_ASSET =
  ' ';

var ALL_CITIES_ASSET =
  ' ';

var YEAR = 2024;

var MONTHS = [1, 12];

var CLOUD_PCT = 10;

var SCALE = 10;

var N_TREE = 250;

var N_GRASS = 250;

var N_IMPERVIOUS = 250;

var N_SOIL = 250;

var N_MIXED = 500;

var ERODE_PIX = 3;


var FISHER_BANDS = [

  'B2',
  'B3',
  'B4',
  'B5',
  'B6',
  'B7',
  'B8',
  'B8A',
  'B11',
  'B12'

];


var FISHER_DIMS = 3;
var FISHER_REG = 1e-6;
var KNN_K = 8;
var KNN_MIN = 4;
var KNN_RADIUS = 3.0;
var PSEUDO_KEEP_PERCENTILE = 90;
var FCLS_RIDGE = 1e-8;
var RF_TREES = 200;
var RF_SEED = 2026;
var LOGRATIO_EPS = 0.001;
var TRAIN_BUF = 10000;
var DRIVE_FOLDER =
  'GEE_LocalFisher_CompositionalRF';
var KEY_BAND = 'B4_mean';

/***************************************************************
PREPROCESSING
 ***************************************************************/

function s2Collection(region, year, months) {

  var start =
    ee.Date.fromYMD(
      year,
      months[0],
      1
    );

  var end =
    ee.Date.fromYMD(
      year,
      months[1],
      1
    ).advance(
      1,
      'month'
    );

  return ee.ImageCollection(
      'COPERNICUS/S2_SR'
    )
    .filterBounds(region)
    .filterDate(
      start,
      end
    )
    .filter(
      ee.Filter.lt(
        'CLOUDY_PIXEL_PERCENTAGE',
        CLOUD_PCT
      )
    );
}

function addIndices(img) {

  var ndvi =
    img.normalizedDifference([
      'B8',
      'B4'
    ])
    .rename(
      'NDVI'
    );

  var evi =
    img.expression(
      '2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
      {

        NIR:
          img.select('B8'),

        RED:
          img.select('B4'),

        BLUE:
          img.select('B2')

      }
    )
    .rename(
      'EVI'
    );

  var gci =
    img.expression(
      '(NIR/GREEN)-1',
      {

        NIR:
          img.select('B8'),

        GREEN:
          img.select('B3')

      }
    )
    .rename(
      'GCI'
    );

  var savi =
    img.expression(
      '((NIR-RED)/(NIR+RED+0.5))*1.5',
      {

        NIR:
          img.select('B8'),

        RED:
          img.select('B4')

      }
    )
    .rename(
      'SAVI'
    );

  var ndbi =
    img.normalizedDifference([
      'B11',
      'B8'
    ])
    .rename(
      'NDBI'
    );

  var ndwi =
    img.normalizedDifference([
      'B8',
      'B11'
    ])
    .rename(
      'NDWI'
    );

  var base =
    img.select([

      'B2',
      'B3',
      'B4',
      'B5',
      'B6',
      'B7',
      'B8',
      'B8A',
      'B11',
      'B12'

    ]);

  return base.addBands([

    ndvi,
    evi,
    gci,
    savi,
    ndbi,
    ndwi

  ]);
}

/***************************************************************
 * MODULE 3: SAMPLE CONSTRUCTION
 ***************************************************************/

function maskCount(
  img,
  region,
  scale
) {

  var d =
    ee.Dictionary(
      img.mask().reduceRegion({

        reducer:
          ee.Reducer.sum(),

        geometry:
          region,

        scale:
          scale,

        maxPixels:
          1e12,

        bestEffort:
          true

      })
    );

  var vals =
    d.values();

  return ee.Number(
    ee.Algorithms.If(

      ee.Number(
        vals.size()
      ).gt(0),

      vals.get(0),

      0

    )
  );
}


function buildPureAndMixedMasks(region) {

  var wc =
    ee.Image(
      'ESA/WorldCover/v200/2021'
    )
    .select(
      'Map'
    )
    .clip(
      region
    );

  var s2Med =
    s2Collection(
      region,
      YEAR,
      MONTHS
    )
    .map(
      addIndices
    )
    .median()
    .clip(
      region
    );

  var ndvi =
    s2Med.select(
      'NDVI'
    );

  var ndwi =
    s2Med.select(
      'NDWI'
    );

  var ndbi =
    s2Med.select(
      'NDBI'
    );
  var dw =
    ee.ImageCollection(
      'GOOGLE/DYNAMICWORLD/V1'
    )
    .filterBounds(
      region
    )
    .filterDate(

      ee.Date.fromYMD(
        YEAR,
        1,
        1
      ),

      ee.Date.fromYMD(
        YEAR,
        12,
        1
      ).advance(
        1,
        'month'
      )

    )
    .select(
      'label'
    )
    .mode()
    .clip(
      region
    );

  var dwTrees =
    dw.eq(1);

  var dwGrass =
    dw.eq(2);

  var dwBuilt =
    dw.eq(6);

  var dwBare =
    dw.eq(7);

  var water =
    ndwi.gte(
      0.30
    );
  var tree =
    wc.eq(10)
      .and(
        dwTrees
      )

      .focal_min({

        radius:
          ERODE_PIX,

        units:
          'pixels'

      })

      .updateMask(

        ndvi.gte(0.50)

          .and(
            ndwi.lt(0.20)
          )

          .and(
            dwBuilt.not()
          )

          .and(
            dwBare.not()
          )

      )

      .selfMask();
  var grass =
    wc.eq(30)
      .and(
        dwGrass
      )

      .focal_min({

        radius:
          ERODE_PIX,

        units:
          'pixels'

      })

      .updateMask(

        ndvi.gte(0.25)

          .and(
            ndvi.lt(0.65)
          )

          .and(
            ndwi.lt(0.20)
          )

          .and(
            dwBuilt.not()
          )

      )

      .selfMask();
  var imperviousBase =
    wc.eq(50)
      .and(
        dwBuilt
      );

  var impervious =
    imperviousBase

      .focal_min({

        radius:
          ERODE_PIX,

        units:
          'pixels'

      })

      .updateMask(

        ndvi.lt(0.25)

          .and(
            ndwi.lt(0.20)
          )

          .and(
            ndbi.gt(-0.10)
          )

      )

      .selfMask();

  var soilBase =
    wc.eq(60)

      .and(
        dwBare
      )

      .and(
        imperviousBase.not()
      );

  var soil =
    soilBase

      .focal_min({

        radius:
          ERODE_PIX,

        units:
          'pixels'

      })

      .updateMask(

        ndvi.lt(0.30)

          .and(
            ndwi.lt(0.20)
          )

          .and(
            water.not()
          )

      )

      .selfMask();
  var pureUnion =
    tree.unmask(0)

      .or(
        grass.unmask(0)
      )

      .or(
        impervious.unmask(0)
      )

      .or(
        soil.unmask(0)
      );

  var valid =
    s2Med
      .select('B4')
      .mask()

      .and(
        water.not()
      );

  var mixed =
    valid

      .and(
        pureUnion.not()
      )

      .selfMask()

      .rename(
        'mixed'
      );

  return {

    tree:
      tree,

    grass:
      grass,

    impervious:
      impervious,

    soil:
      soil,

    mixed:
      mixed

  };
}
function stratifiedPointsSafe(
  mask,
  classValue,
  n,
  region,
  name
) {

  var count =
    maskCount(
      mask,
      region,
      SCALE
    );

  return ee.FeatureCollection(

    ee.Algorithms.If(

      count.gt(0),

      mask
        .rename(
          'class'
        )
        .multiply(
          classValue
        )
        .toInt()
        .stratifiedSample({

          numPoints:
            n,

          classBand:
            'class',

          region:
            region,

          scale:
            SCALE,

          classValues:
            [classValue],

          classPoints:
            [n],

          seed:
            42,

          geometries:
            true

        }),

      ee.FeatureCollection([])

    )
  );
}

function randomMixedPoints(
  mask,
  n,
  region
) {

  return mask.sample({

      region:
        region,

      scale:
        SCALE,

      numPixels:
        n,

      seed:
        2026,

      geometries:
        true,

      tileScale:
        4

    })

    .map(
      function(f) {

        return f.set({

          sample_type:
            'mixed_auto'

        });

      }
    );
}
function tryLoadManualSamples(path) {

  if (
    !path ||
    path === ''
  ) {

    return {

      has:
        false,

      fc:
        ee.FeatureCollection([])

    };

  }

  var has =
    false;

  var fc =
    ee.FeatureCollection([]);

  try {

    fc =
      ee.FeatureCollection(
        path
      );

    fc.limit(1)
      .size()
      .getInfo();

    has =
      true;

  }

  catch (err) {

    has =
      false;

  }

  return {

    has:
      has,

    fc:
      fc

  };
}

function normalizeManualLabels(fc) {

  return ee.FeatureCollection(fc)

    .map(
      function(f) {

        var names =
          f.propertyNames();

        var hasT =
          names.contains(
            'TreeRatio'
          );

        var hasG =
          names.contains(
            'GrassRatio'
          );

        var hasO =
          names.contains(
            'OtherRatio'
          );

        var hasClass =
          names.contains(
            'class'
          );

        var hasRatio =
          ee.Algorithms.If(

            hasT,

            ee.Algorithms.If(

              hasG,

              ee.Algorithms.If(

                hasO,

                true,

                false

              ),

              false

            ),

            false

          );

        var tRaw =
          ee.Number(
            ee.Algorithms.If(

              hasRatio,

              f.get(
                'TreeRatio'
              ),

              ee.Algorithms.If(

                hasClass,

                ee.Algorithms.If(

                  ee.Number(
                    f.get('class')
                  ).eq(1),

                  1,

                  0

                ),

                -1

              )

            )
          );

        var gRaw =
          ee.Number(
            ee.Algorithms.If(

              hasRatio,

              f.get(
                'GrassRatio'
              ),

              ee.Algorithms.If(

                hasClass,

                ee.Algorithms.If(

                  ee.Number(
                    f.get('class')
                  ).eq(2),

                  1,

                  0

                ),

                -1

              )

            )
          );

        var oRaw =
          ee.Number(
            ee.Algorithms.If(

              hasRatio,

              f.get(
                'OtherRatio'
              ),

              ee.Algorithms.If(

                hasClass,

                ee.Algorithms.If(

                  ee.Number(
                    f.get('class')
                  ).eq(3),

                  1,

                  0

                ),

                -1

              )

            )
          );

        var ok =
          ee.Number(
            ee.Algorithms.If(
              tRaw.gte(0),
              1,
              0
            )
          )

          .multiply(

            ee.Number(
              ee.Algorithms.If(
                gRaw.gte(0),
                1,
                0
              )
            )

          )

          .multiply(

            ee.Number(
              ee.Algorithms.If(
                oRaw.gte(0),
                1,
                0
              )
            )

          );

        var t =
          tRaw.max(0);

        var g =
          gRaw.max(0);

        var o =
          oRaw.max(0);

        var sum =
          t.add(g)
            .add(o);

        var tN =
          ee.Number(
            ee.Algorithms.If(

              sum.gt(0),

              t.divide(sum),

              0

            )
          );

        var gN =
          ee.Number(
            ee.Algorithms.If(

              sum.gt(0),

              g.divide(sum),

              0

            )
          );

        var oN =
          ee.Number(
            ee.Algorithms.If(

              sum.gt(0),

              o.divide(sum),

              0

            )
          );

        var finalOK =
          ok.multiply(

            ee.Number(
              ee.Algorithms.If(

                sum.gt(0),

                1,

                0

              )
            )

          );

        return f.set({

          TreeRatio:
            tN,

          GrassRatio:
            gN,

          OtherRatio:
            oN,

          sample_type:
            'manual',

          label_ok:
            finalOK

        });

      }
    )

    .filter(
      ee.Filter.eq(
        'label_ok',
        1
      )
    );
}

/***************************************************************
 * MODULE 4: FISHER PROJECTION UTILITIES
 ***************************************************************/

function featureMeanVector(
  fc,
  bands
) {

  bands =
    ee.List(
      bands
    );

  return ee.Array(

    bands.map(
      function(b) {

        return fc.aggregate_mean(
          ee.String(b)
        );

      }
    )

  );
}

function featureVector(
  feature,
  bands
) {

  feature =
    ee.Feature(
      feature
    );

  bands =
    ee.List(
      bands
    );

  return ee.Array(

    bands.map(
      function(b) {

        return ee.Number(
          feature.get(
            ee.String(b)
          )
        );

      }
    )

  );
}

function outerProduct(
  vector,
  p
) {

  vector =
    ee.Array(
      vector
    );

  var column =
    vector.reshape(

      ee.List([
        p,
        1
      ])

    );

  return column.matrixMultiply(
    column.matrixTranspose()
  );
}

function classScatter(
  fc,
  bands,
  meanVector
) {

  var p =
    ee.Number(
      ee.List(
        bands
      ).size()
    );

  var zero =
    ee.Array.identity(
      p
    )
    .multiply(0);

  return ee.Array(

    ee.FeatureCollection(fc)
      .iterate(

        function(item, accumulator) {

          var f =
            ee.Feature(
              item
            );

          var acc =
            ee.Array(
              accumulator
            );

          var x =
            featureVector(
              f,
              bands
            );

          var d =
            x.subtract(
              meanVector
            );

          return acc.add(

            outerProduct(
              d,
              p
            )

          );

        },

        zero

      )

  );
}
function buildFisherProjection(
  treeFC,
  grassFC,
  imperviousFC,
  soilFC,
  bands
) {

  bands =
    ee.List(
      bands
    );

  var p =
    ee.Number(
      bands.size()
    );

  var nT =
    ee.Number(
      treeFC.size()
    );

  var nG =
    ee.Number(
      grassFC.size()
    );

  var nI =
    ee.Number(
      imperviousFC.size()
    );

  var nS =
    ee.Number(
      soilFC.size()
    );

  var nAll =
    nT
      .add(nG)
      .add(nI)
      .add(nS);

  var mT =
    featureMeanVector(
      treeFC,
      bands
    );

  var mG =
    featureMeanVector(
      grassFC,
      bands
    );

  var mI =
    featureMeanVector(
      imperviousFC,
      bands
    );

  var mS =
    featureMeanVector(
      soilFC,
      bands
    );

  var meanAll =
    mT.multiply(nT)

      .add(
        mG.multiply(nG)
      )

      .add(
        mI.multiply(nI)
      )

      .add(
        mS.multiply(nS)
      )

      .divide(
        nAll
      );

  var Sw =
    classScatter(
      treeFC,
      bands,
      mT
    )

    .add(
      classScatter(
        grassFC,
        bands,
        mG
      )
    )

    .add(
      classScatter(
        imperviousFC,
        bands,
        mI
      )
    )

    .add(
      classScatter(
        soilFC,
        bands,
        mS
      )
    )

    .divide(
      nAll.max(1)
    );

  var dT =
    mT.subtract(
      meanAll
    );

  var dG =
    mG.subtract(
      meanAll
    );

  var dI =
    mI.subtract(
      meanAll
    );

  var dS =
    mS.subtract(
      meanAll
    );

  var Sb =
    outerProduct(
      dT,
      p
    )
    .multiply(
      nT
    )

    .add(

      outerProduct(
        dG,
        p
      )
      .multiply(
        nG
      )

    )

    .add(

      outerProduct(
        dI,
        p
      )
      .multiply(
        nI
      )

    )

    .add(

      outerProduct(
        dS,
        p
      )
      .multiply(
        nS
      )

    )

    .divide(
      nAll.max(1)
    );
  var identity =
    ee.Array.identity(
      p
    );

  var SwReg =
    Sw.add(

      identity.multiply(
        FISHER_REG
      )

    );
  var swEigen =
    SwReg.eigen();

  var swValues =
    swEigen.slice(
      1,
      0,
      1
    );

  var swVectors =
    swEigen.slice(
      1,
      1
    );

  var invSqrtValues =
    swValues

      .add(
        1e-12
      )

      .pow(
        -0.5
      );

  var DInvSqrt =
    invSqrtValues
      .matrixToDiag();

  var W =
    DInvSqrt
      .matrixMultiply(
        swVectors
      );

  var SbWhite =
    W

      .matrixMultiply(
        Sb
      )

      .matrixMultiply(
        W.matrixTranspose()
      );

  var fisherEigen =
    SbWhite.eigen();

  var fisherVectors =
    fisherEigen.slice(
      1,
      1
    );

  var topVectors =
    fisherVectors.slice(
      0,
      0,
      FISHER_DIMS
    );

  var F =
    topVectors.matrixMultiply(
      W
    );

  return {

    F:
      F,

    Sw:
      Sw,

    Sb:
      Sb

  };
}

function addFisherCoordinates(
  fc,
  F,
  bands
) {

  var p =
    ee.Number(
      ee.List(
        bands
      ).size()
    );

  return ee.FeatureCollection(fc)

    .map(
      function(f) {

        var x =
          featureVector(
            f,
            bands
          )
          .reshape(

            ee.List([
              p,
              1
            ])

          );

        var y =
          ee.Array(F)
            .matrixMultiply(
              x
            );

        return f.set({

          F1:
            y.get([
              0,
              0
            ]),

          F2:
            y.get([
              1,
              0
            ]),

          F3:
            y.get([
              2,
              0
            ])

        });

      }
    );
}

function fisherCoordinateStats(fc) {

  var f1 =
    ee.List(
      fc.aggregate_array(
        'F1'
      )
    );

  var f2 =
    ee.List(
      fc.aggregate_array(
        'F2'
      )
    );

  var f3 =
    ee.List(
      fc.aggregate_array(
        'F3'
      )
    );

  return ee.Dictionary({

    m1:
      ee.Number(
        f1.reduce(
          ee.Reducer.mean()
        )
      ),

    m2:
      ee.Number(
        f2.reduce(
          ee.Reducer.mean()
        )
      ),

    m3:
      ee.Number(
        f3.reduce(
          ee.Reducer.mean()
        )
      ),

    s1:
      ee.Number(
        f1.reduce(
          ee.Reducer.stdDev()
        )
      ).max(
        1e-9
      ),

    s2:
      ee.Number(
        f2.reduce(
          ee.Reducer.stdDev()
        )
      ).max(
        1e-9
      ),

    s3:
      ee.Number(
        f3.reduce(
          ee.Reducer.stdDev()
        )
      ).max(
        1e-9
      )

  });
}

function addStandardizedFisher(
  fc,
  stats
) {

  stats =
    ee.Dictionary(
      stats
    );

  var m1 =
    ee.Number(
      stats.get(
        'm1'
      )
    );

  var m2 =
    ee.Number(
      stats.get(
        'm2'
      )
    );

  var m3 =
    ee.Number(
      stats.get(
        'm3'
      )
    );

  var s1 =
    ee.Number(
      stats.get(
        's1'
      )
    );

  var s2 =
    ee.Number(
      stats.get(
        's2'
      )
    );

  var s3 =
    ee.Number(
      stats.get(
        's3'
      )
    );

  return ee.FeatureCollection(fc)

    .map(
      function(f) {

        return f.set({

          Z1:
            ee.Number(
              f.get(
                'F1'
              )
            )
            .subtract(m1)
            .divide(s1),

          Z2:
            ee.Number(
              f.get(
                'F2'
              )
            )
            .subtract(m2)
            .divide(s2),

          Z3:
            ee.Number(
              f.get(
                'F3'
              )
            )
            .subtract(m3)
            .divide(s3)

        });

      }
    );
}

/***************************************************************
 * MODULE 5: ENDMEMBER SELECTION AND UNCERTAINTY
 ***************************************************************/

function packLibrary(
  fc,
  bands
) {

  fc =
    ee.FeatureCollection(fc);

  bands =
    ee.List(
      bands
    );

  var size =
    fc.size();

  var z1 =
    fc.aggregate_array(
      'Z1'
    );

  var z2 =
    fc.aggregate_array(
      'Z2'
    );

  var z3 =
    fc.aggregate_array(
      'Z3'
    );

  var bandDict =
    ee.Dictionary(
      bands.iterate(

        function(b, accumulator) {

          accumulator =
            ee.Dictionary(
              accumulator
            );

          b =
            ee.String(
              b
            );

          return accumulator.set(

            b,

            fc.aggregate_array(
              b
            )

          );

        },

        ee.Dictionary({})

      )
    );

  return ee.Dictionary({

    size:
      size,

    z1:
      z1,

    z2:
      z2,

    z3:
      z3,

    bands:
      bandDict

  });
}

function packedKNN(
  mixedFeature,
  packed
) {

  mixedFeature =
    ee.Feature(
      mixedFeature
    );

  packed =
    ee.Dictionary(
      packed
    );

  var n =
    ee.Number(
      packed.get(
        'size'
      )
    );

  var z1List =
    ee.List(
      packed.get(
        'z1'
      )
    );

  var z2List =
    ee.List(
      packed.get(
        'z2'
      )
    );

  var z3List =
    ee.List(
      packed.get(
        'z3'
      )
    );

  var q1 =
    ee.Number(
      mixedFeature.get(
        'Z1'
      )
    );

  var q2 =
    ee.Number(
      mixedFeature.get(
        'Z2'
      )
    );

  var q3 =
    ee.Number(
      mixedFeature.get(
        'Z3'
      )
    );

  var indices =
    ee.List.sequence(
      0,
      n.subtract(1)
    );

  var distances =
    indices.map(
      function(index) {

        index =
          ee.Number(
            index
          );

        var d1 =
          ee.Number(
            z1List.get(index)
          )
          .subtract(q1);

        var d2 =
          ee.Number(
            z2List.get(index)
          )
          .subtract(q2);

        var d3 =
          ee.Number(
            z3List.get(index)
          )
          .subtract(q3);

        return d1.pow(2)

          .add(
            d2.pow(2)
          )

          .add(
            d3.pow(2)
          )

          .sqrt();

      }
    );

  var validIndices =
    indices

      .map(
        function(index) {

          index =
            ee.Number(index);

          var distance =
            ee.Number(
              distances.get(index)
            );

          return ee.Algorithms.If(

            distance.lte(
              KNN_RADIUS
            ),

            index,

            null

          );

        }
      )

      .removeAll([
        null
      ]);

  var validDistances =
    validIndices.map(
      function(index) {

        return distances.get(
          ee.Number(index)
        );

      }
    );

  var sorted =
    validIndices.sort(
      validDistances
    );

  var topK =
    sorted.slice(
      0,
      KNN_K
    );

  return ee.Dictionary({

    indices:
      topK,

    radiusCount:
      validIndices.size(),

    distances:
      distances

  });
}

function localEndmemberFromPacked(
  packed,
  indices,
  bands
) {

  packed =
    ee.Dictionary(
      packed
    );

  indices =
    ee.List(
      indices
    );

  bands =
    ee.List(
      bands
    );

  var bandDict =
    ee.Dictionary(
      packed.get(
        'bands'
      )
    );

  var vector =
    bands.map(
      function(b) {

        b =
          ee.String(
            b
          );

        var source =
          ee.List(
            bandDict.get(
              b
            )
          );

        var selected =
          indices.map(
            function(index) {

              return source.get(
                ee.Number(index)
              );

            }
          );

        return ee.Number(
          selected.reduce(
            ee.Reducer.mean()
          )
        );

      }
    );

  return ee.Array(
    vector
  );
}

function globalPackedMean(
  packed,
  bands
) {

  packed =
    ee.Dictionary(
      packed
    );

  bands =
    ee.List(
      bands
    );

  var bandDict =
    ee.Dictionary(
      packed.get(
        'bands'
      )
    );

  return ee.Array(

    bands.map(
      function(b) {

        return ee.Number(

          ee.List(
            bandDict.get(
              ee.String(b)
            )
          )
          .reduce(
            ee.Reducer.mean()
          )

        );

      }
    )

  );
}

function safeLocalEndmember(
  packed,
  knnResult,
  bands
) {

  packed =
    ee.Dictionary(
      packed
    );

  knnResult =
    ee.Dictionary(
      knnResult
    );

  var indices =
    ee.List(
      knnResult.get(
        'indices'
      )
    );

  var local =
    ee.Array(
      ee.Algorithms.If(

        indices.size().gt(0),

        localEndmemberFromPacked(
          packed,
          indices,
          bands
        ),

        globalPackedMean(
          packed,
          bands
        )

      )
    );

  return local;
}

function packedDispersion(
  packed,
  indices
) {

  packed =
    ee.Dictionary(
      packed
    );

  indices =
    ee.List(
      indices
    );

  var z1 =
    ee.List(
      packed.get(
        'z1'
      )
    );

  var z2 =
    ee.List(
      packed.get(
        'z2'
      )
    );

  var z3 =
    ee.List(
      packed.get(
        'z3'
      )
    );

  var selected1 =
    indices.map(
      function(i) {

        return z1.get(
          ee.Number(i)
        );

      }
    );

  var selected2 =
    indices.map(
      function(i) {

        return z2.get(
          ee.Number(i)
        );

      }
    );

  var selected3 =
    indices.map(
      function(i) {

        return z3.get(
          ee.Number(i)
        );

      }
    );

  var result =
    ee.Number(
      ee.Algorithms.If(

        indices.size().gt(1),

        (function() {

          var m1 =
            ee.Number(
              selected1.reduce(
                ee.Reducer.mean()
              )
            );

          var m2 =
            ee.Number(
              selected2.reduce(
                ee.Reducer.mean()
              )
            );

          var m3 =
            ee.Number(
              selected3.reduce(
                ee.Reducer.mean()
              )
            );

          var squared =
            ee.List.sequence(
              0,
              indices.size().subtract(1)
            )
            .map(
              function(k) {

                k =
                  ee.Number(
                    k
                  );

                var d1 =
                  ee.Number(
                    selected1.get(k)
                  )
                  .subtract(m1);

                var d2 =
                  ee.Number(
                    selected2.get(k)
                  )
                  .subtract(m2);

                var d3 =
                  ee.Number(
                    selected3.get(k)
                  )
                  .subtract(m3);

                return d1.pow(2)
                  .add(
                    d2.pow(2)
                  )
                  .add(
                    d3.pow(2)
                  );

              }
            );

          return ee.Number(
            squared.reduce(
              ee.Reducer.mean()
            )
          ).sqrt();

        })(),

        999

      )
    );

  return result;
}

/***************************************************************
 * MODULE 6: PSEUDO-LABEL GENERATION
 ***************************************************************/

var FCLS_SUBSETS = [

  [0],
  [1],
  [2],
  [3],

  [0,1],
  [0,2],
  [0,3],
  [1,2],
  [1,3],
  [2,3],

  [0,1,2],
  [0,1,3],
  [0,2,3],
  [1,2,3],

  [0,1,2,3]

];

function fclsCandidate(
  E,
  y,
  subset
) {

  var m =
    subset.length;

  var columns =
    [];

  for (
    var i = 0;
    i < m;
    i++
  ) {

    columns.push(

      E.slice(
        1,
        subset[i],
        subset[i] + 1
      )

    );

  }

  var Es =
    ee.Array.cat(
      columns,
      1
    );

  var ET =
    Es.matrixTranspose();

  var H =
    ET

      .matrixMultiply(
        Es
      )

      .multiply(2)

      .add(

        ee.Array.identity(
          m
        )
        .multiply(
          FCLS_RIDGE
        )

      );

  var ones =
    [];

  for (
    var j = 0;
    j < m;
    j++
  ) {

    ones.push([
      1
    ]);

  }

  var oneCol =
    ee.Array(
      ones
    );

  var oneRow =
    oneCol.matrixTranspose();

  var top =
    ee.Array.cat(
      [
        H,
        oneCol
      ],
      1
    );

  var bottom =
    ee.Array.cat(
      [

        oneRow,

        ee.Array([
          [0]
        ])

      ],
      1
    );

  var KKT =
    ee.Array.cat(
      [
        top,
        bottom
      ],
      0
    );

  var rhsTop =
    ET

      .matrixMultiply(
        y
      )

      .multiply(2);

  var rhs =
    ee.Array.cat(
      [

        rhsTop,

        ee.Array([
          [1]
        ])

      ],
      0
    );

  var solution =
    KKT

      .matrixInverse()

      .matrixMultiply(
        rhs
      );

  var fSmall =
    solution.slice(
      0,
      0,
      m
    );

  var fullValues =
    [];

  for (
    var k = 0;
    k < 4;
    k++
  ) {

    var position =
      subset.indexOf(
        k
      );

    if (
      position >= 0
    ) {

      fullValues.push(

        ee.Number(
          fSmall.get([
            position,
            0
          ])
        )

      );

    }

    else {

      fullValues.push(
        ee.Number(0)
      );

    }

  }

  var valid =
    ee.Number(1);

  for (
    var q = 0;
    q < 4;
    q++
  ) {

    valid =
      valid.multiply(

        ee.Number(
          ee.Algorithms.If(

            ee.Number(
              fullValues[q]
            )
            .gte(
              -1e-8
            ),

            1,

            0

          )
        )

      );

  }

  var positive =
    [];

  var sum =
    ee.Number(0);

  for (
    var h = 0;
    h < 4;
    h++
  ) {

    var value =
      ee.Number(
        fullValues[h]
      )
      .max(0);

    positive.push(
      value
    );

    sum =
      sum.add(
        value
      );

  }

  var rows =
    [];

  for (
    var z = 0;
    z < 4;
    z++
  ) {

    rows.push([

      ee.Number(
        positive[z]
      )
      .divide(
        sum.max(
          1e-12
        )
      )

    ]);

  }

  var fractions =
    ee.Array(
      rows
    );

  var residual =
    y.subtract(

      E.matrixMultiply(
        fractions
      )

    );

  var sse =
    ee.Number(

      residual

        .matrixTranspose()

        .matrixMultiply(
          residual
        )

        .get([
          0,
          0
        ])

    );

  return ee.Dictionary({

    valid:
      valid,

    sse:
      sse,

    fractions:
      fractions

  });
}

function fclsSolve(
  E,
  y
) {

  E =
    ee.Array(E);

  y =
    ee.Array(y);

  var candidates =
    [];

  for (
    var i = 0;
    i < FCLS_SUBSETS.length;
    i++
  ) {

    candidates.push(

      fclsCandidate(

        E,

        y,

        FCLS_SUBSETS[i]

      )

    );

  }

  var initial =
    ee.Dictionary({

      valid:
        0,

      sse:
        1e30,

      fractions:
        ee.Array([
          [0],
          [0],
          [0],
          [0]
        ])

    });

  return ee.Dictionary(

    ee.List(
      candidates
    )
    .iterate(

      function(item, previous) {

        var current =
          ee.Dictionary(
            item
          );

        var best =
          ee.Dictionary(
            previous
          );

        var valid =
          ee.Number(
            current.get(
              'valid'
            )
          );

        var better =
          ee.Number(
            current.get(
              'sse'
            )
          )
          .lt(
            ee.Number(
              best.get(
                'sse'
              )
            )
          );

        var choose =
          valid.multiply(

            ee.Number(
              ee.Algorithms.If(

                better,

                1,

                0

              )
            )

          );

        return ee.Algorithms.If(

          choose.eq(1),

          current,

          best

        );

      },

      initial

    )

  );
}

function localizedPseudoLabel(
  feature,
  treePack,
  grassPack,
  imperviousPack,
  soilPack,
  F
) {

  feature =
    ee.Feature(
      feature
    );

  var tKNN =
    packedKNN(
      feature,
      treePack
    );

  var gKNN =
    packedKNN(
      feature,
      grassPack
    );

  var iKNN =
    packedKNN(
      feature,
      imperviousPack
    );

  var sKNN =
    packedKNN(
      feature,
      soilPack
    );

  var nT =
    ee.Number(
      tKNN.get(
        'radiusCount'
      )
    );

  var nG =
    ee.Number(
      gKNN.get(
        'radiusCount'
      )
    );

  var nI =
    ee.Number(
      iKNN.get(
        'radiusCount'
      )
    );

  var nS =
    ee.Number(
      sKNN.get(
        'radiusCount'
      )
    );

  var valid =
    ee.Number(
      ee.Algorithms.If(
        nT.gte(
          KNN_MIN
        ),
        1,
        0
      )
    )

    .multiply(

      ee.Number(
        ee.Algorithms.If(
          nG.gte(
            KNN_MIN
          ),
          1,
          0
        )
      )

    )

    .multiply(

      ee.Number(
        ee.Algorithms.If(
          nI.gte(
            KNN_MIN
          ),
          1,
          0
        )
      )

    )

    .multiply(

      ee.Number(
        ee.Algorithms.If(
          nS.gte(
            KNN_MIN
          ),
          1,
          0
        )
      )

    );

  var p =
    ee.Number(
      ee.List(
        FISHER_BANDS
      ).size()
    );

  var tEnd =
    safeLocalEndmember(
      treePack,
      tKNN,
      FISHER_BANDS
    )
    .reshape(
      ee.List([
        p,
        1
      ])
    );

  var gEnd =
    safeLocalEndmember(
      grassPack,
      gKNN,
      FISHER_BANDS
    )
    .reshape(
      ee.List([
        p,
        1
      ])
    );

  var iEnd =
    safeLocalEndmember(
      imperviousPack,
      iKNN,
      FISHER_BANDS
    )
    .reshape(
      ee.List([
        p,
        1
      ])
    );

  var sEnd =
    safeLocalEndmember(
      soilPack,
      sKNN,
      FISHER_BANDS
    )
    .reshape(
      ee.List([
        p,
        1
      ])
    );

  var EOriginal =
    ee.Array.cat(
      [

        tEnd,
        gEnd,
        iEnd,
        sEnd

      ],
      1
    );

  var EFisher =
    ee.Array(F)
      .matrixMultiply(
        EOriginal
      );

  var y =
    ee.Array([

      [
        ee.Number(
          feature.get(
            'F1'
          )
        )
      ],

      [
        ee.Number(
          feature.get(
            'F2'
          )
        )
      ],

      [
        ee.Number(
          feature.get(
            'F3'
          )
        )
      ]

    ]);

  var solution =
    fclsSolve(
      EFisher,
      y
    );

  var fractions =
    ee.Array(
      solution.get(
        'fractions'
      )
    );

  var tree =
    ee.Number(
      fractions.get([
        0,
        0
      ])
    );

  var grass =
    ee.Number(
      fractions.get([
        1,
        0
      ])
    );

  var impervious =
    ee.Number(
      fractions.get([
        2,
        0
      ])
    );

  var soil =
    ee.Number(
      fractions.get([
        3,
        0
      ])
    );

  var other =
    impervious.add(
      soil
    );

  var rmse =
    ee.Number(
      solution.get(
        'sse'
      )
    )
    .divide(
      FISHER_DIMS
    )
    .sqrt();

  var uTree =
    packedDispersion(

      treePack,

      ee.Dictionary(
        tKNN
      ).get(
        'indices'
      )

    );

  var uGrass =
    packedDispersion(

      grassPack,

      ee.Dictionary(
        gKNN
      ).get(
        'indices'
      )

    );

  var uImp =
    packedDispersion(

      imperviousPack,

      ee.Dictionary(
        iKNN
      ).get(
        'indices'
      )

    );

  var uSoil =
    packedDispersion(

      soilPack,

      ee.Dictionary(
        sKNN
      ).get(
        'indices'
      )

    );

  var uncertainty =
    uTree
      .add(uGrass)
      .add(uImp)
      .add(uSoil)
      .divide(4);

  return feature.set({

    TreeRatio:
      tree,

    GrassRatio:
      grass,

    ImperviousRatio:
      impervious,

    SoilRatio:
      soil,

    OtherRatio:
      other,

    Fisher_RMSE:
      rmse,

    EndmemberUncertainty:
      uncertainty,

    nTreeNN:
      nT,

    nGrassNN:
      nG,

    nImperviousNN:
      nI,

    nSoilNN:
      nS,

    pseudo_valid:
      valid,

    sample_type:
      'mixed_local_fisher'

  });
}

/***************************************************************
 * MODULE 7: PREDICTOR CONSTRUCTION
 ***************************************************************/

function runAnalysis(
  region,
  cityName,
  pureTree,
  pureGrass,
  pureImpervious,
  pureSoil,
  mixedPoints,
  manualSamples
) {

  var allGeometry =
    pureTree
      .merge(
        pureGrass
      )
      .merge(
        pureImpervious
      )
      .merge(
        pureSoil
      )
      .merge(
        mixedPoints
      )
      .merge(
        manualSamples
      );

  var trainAOI =
    allGeometry
      .geometry()
      .buffer(
        TRAIN_BUF
      );

  var S2_BANDS_ALL = [

    'B2',
    'B3',
    'B4',
    'B5',
    'B6',
    'B7',
    'B8',
    'B8A',
    'B11',
    'B12',

    'NDVI',
    'EVI',
    'GCI',
    'SAVI',
    'NDBI',
    'NDWI'

  ];

  var reducer =
    ee.Reducer.mean()

      .combine(
        ee.Reducer.min(),
        '',
        true
      )

      .combine(
        ee.Reducer.max(),
        '',
        true
      )

      .combine(
        ee.Reducer.stdDev(),
        '',
        true
      )

      .combine(

        ee.Reducer.percentile(
          [25,50,75],
          [
            'p25',
            'p50',
            'p75'
          ]
        ),

        '',
        true

      );

  var s2City =
    s2Collection(
      region,
      YEAR,
      MONTHS
    );

  var baseCity =
    s2City
      .map(
        addIndices
      )
      .select(
        S2_BANDS_ALL
      )
      .reduce(
        reducer
      )
      .clip(
        region
      )
      .unmask(
        -9999
      );

  var fisherCity =
    s2City
      .select(
        FISHER_BANDS
      )
      .median()
      .multiply(
        0.0001
      )
      .clip(
        region
      );

  var s2Train =
    s2Collection(
      trainAOI,
      YEAR,
      MONTHS
    );

  var baseTrain =
    s2Train
      .map(
        addIndices
      )
      .select(
        S2_BANDS_ALL
      )
      .reduce(
        reducer
      )
      .clip(
        trainAOI
      )
      .unmask(
        -9999
      );

  var fisherTrain =
    s2Train
      .select(
        FISHER_BANDS
      )
      .median()
      .multiply(
        0.0001
      )
      .clip(
        trainAOI
      );

  var d0 =
    ee.Date.fromYMD(
      YEAR,
      1,
      1
    );

  var d1 =
    ee.Date.fromYMD(
      YEAR,
      12,
      1
    ).advance(
      1,
      'month'
    );

  function getS1(regionX) {

    return ee.ImageCollection(
        'COPERNICUS/S1_GRD'
      )

      .filterBounds(
        regionX
      )

      .filterDate(
        d0,
        d1
      )

      .filter(
        ee.Filter.eq(
          'instrumentMode',
          'IW'
        )
      )

      .filter(
        ee.Filter.eq(
          'orbitProperties_pass',
          'DESCENDING'
        )
      )

      .filter(
        ee.Filter.listContains(
          'transmitterReceiverPolarisation',
          'VV'
        )
      )

      .filter(
        ee.Filter.listContains(
          'transmitterReceiverPolarisation',
          'VH'
        )
      );
  }

  var s1City =
    getS1(
      region
    );

  var vvCity =
    s1City
      .select(
        'VV'
      )
      .median()
      .rename(
        'VV_median'
      )
      .unmask(
        -9999
      );

  var vhCity =
    s1City
      .select(
        'VH'
      )
      .median()
      .rename(
        'VH_median'
      )
      .unmask(
        -9999
      );

  var s1Train =
    getS1(
      trainAOI
    );

  var vvTrain =
    s1Train
      .select(
        'VV'
      )
      .median()
      .rename(
        'VV_median'
      )
      .unmask(
        -9999
      );

  var vhTrain =
    s1Train
      .select(
        'VH'
      )
      .median()
      .rename(
        'VH_median'
      )
      .unmask(
        -9999
      );

  var dem =
    ee.Image(
      'USGS/SRTMGL1_003'
    )
    .unmask(0);

  var elevation =
    dem.rename(
      'elevation'
    );

  var slope =
    ee.Terrain
      .slope(
        dem
      )
      .rename(
        'slope'
      );

  var aspect =
    ee.Terrain
      .aspect(
        dem
      )
      .rename(
        'aspect'
      );

  var featureCity =
    baseCity
      .addBands([
        vhCity,
        vvCity,
        elevation,
        slope,
        aspect
      ])
      .clip(
        region
      );

  var featureTrain =
    baseTrain
      .addBands([
        vhTrain,
        vvTrain,
        elevation,
        slope,
        aspect
      ])
      .clip(
        trainAOI
      );

  var allBands =
    featureCity
      .bandNames();

  var treeLib =
    fisherCity.sampleRegions({

      collection:
        pureTree,

      scale:
        SCALE,

      geometries:
        false

    });

  var grassLib =
    fisherCity.sampleRegions({

      collection:
        pureGrass,

      scale:
        SCALE,

      geometries:
        false

    });

  var impLib =
    fisherCity.sampleRegions({

      collection:
        pureImpervious,

      scale:
        SCALE,

      geometries:
        false

    });

  var soilLib =
    fisherCity.sampleRegions({

      collection:
        pureSoil,

      scale:
        SCALE,

      geometries:
        false

    });

  var fisherModel =
    buildFisherProjection(

      treeLib,

      grassLib,

      impLib,

      soilLib,

      FISHER_BANDS

    );

  var F =
    fisherModel.F;

  treeLib =
    addFisherCoordinates(
      treeLib,
      F,
      FISHER_BANDS
    );

  grassLib =
    addFisherCoordinates(
      grassLib,
      F,
      FISHER_BANDS
    );

  impLib =
    addFisherCoordinates(
      impLib,
      F,
      FISHER_BANDS
    );

  soilLib =
    addFisherCoordinates(
      soilLib,
      F,
      FISHER_BANDS
    );

  var allPure =
    treeLib
      .merge(
        grassLib
      )
      .merge(
        impLib
      )
      .merge(
        soilLib
      );

  var fisherStats =
    fisherCoordinateStats(
      allPure
    );

  treeLib =
    addStandardizedFisher(
      treeLib,
      fisherStats
    );

  grassLib =
    addStandardizedFisher(
      grassLib,
      fisherStats
    );

  impLib =
    addStandardizedFisher(
      impLib,
      fisherStats
    );

  soilLib =
    addStandardizedFisher(
      soilLib,
      fisherStats
    );

  var treePack =
    packLibrary(
      treeLib,
      FISHER_BANDS
    );

  var grassPack =
    packLibrary(
      grassLib,
      FISHER_BANDS
    );

  var impPack =
    packLibrary(
      impLib,
      FISHER_BANDS
    );

  var soilPack =
    packLibrary(
      soilLib,
      FISHER_BANDS
    );

  var mixedSpectral =
    fisherTrain
      .sampleRegions({

        collection:
          mixedPoints,

        scale:
          SCALE,

        geometries:
          true,

        tileScale:
          4

      });

  mixedSpectral =
    addFisherCoordinates(

      mixedSpectral,

      F,

      FISHER_BANDS

    );

  mixedSpectral =
    addStandardizedFisher(

      mixedSpectral,

      fisherStats

    );

  var pseudoAll =
    mixedSpectral.map(
      function(f) {

        return localizedPseudoLabel(

          f,

          treePack,

          grassPack,

          impPack,

          soilPack,

          F

        );

      }
    );

  var pseudoValid =
    pseudoAll.filter(
      ee.Filter.eq(
        'pseudo_valid',
        1
      )
    );

  var validCount =
    pseudoValid.size();

  var percentileName =
    'p' +
    PSEUDO_KEEP_PERCENTILE;

  var rmseThreshold =
    ee.Number(
      ee.Algorithms.If(

        validCount.gt(0),

        pseudoValid.reduceColumns({

          reducer:
            ee.Reducer.percentile([
              PSEUDO_KEEP_PERCENTILE
            ]),

          selectors:
            [
              'Fisher_RMSE'
            ]

        })
        .get(
          percentileName
        ),

        1e30

      )
    );

  var uncertaintyThreshold =
    ee.Number(
      ee.Algorithms.If(

        validCount.gt(0),

        pseudoValid.reduceColumns({

          reducer:
            ee.Reducer.percentile([
              PSEUDO_KEEP_PERCENTILE
            ]),

          selectors:
            [
              'EndmemberUncertainty'
            ]

        })
        .get(
          percentileName
        ),

        1e30

      )
    );

  var mixedSoft =
    pseudoValid

      .filter(
        ee.Filter.lte(
          'Fisher_RMSE',
          rmseThreshold
        )
      )

      .filter(
        ee.Filter.lte(
          'EndmemberUncertainty',
          uncertaintyThreshold
        )
      );

  var pureTreeLabels =
    pureTree.map(
      function(f) {

        return f.set({

          TreeRatio:
            1,

          GrassRatio:
            0,

          OtherRatio:
            0,

          sample_type:
            'pure_tree'

        });

      }
    );

  var pureGrassLabels =
    pureGrass.map(
      function(f) {

        return f.set({

          TreeRatio:
            0,

          GrassRatio:
            1,

          OtherRatio:
            0,

          sample_type:
            'pure_grass'

        });

      }
    );

  var pureImpLabels =
    pureImpervious.map(
      function(f) {

        return f.set({

          TreeRatio:
            0,

          GrassRatio:
            0,

          OtherRatio:
            1,

          sample_type:
            'pure_impervious'

        });

      }
    );

  var pureSoilLabels =
    pureSoil.map(
      function(f) {

        return f.set({

          TreeRatio:
            0,

          GrassRatio:
            0,

          OtherRatio:
            1,

          sample_type:
            'pure_soil'

        });

      }
    );

  var labeled =
    pureTreeLabels

      .merge(
        pureGrassLabels
      )

      .merge(
        pureImpLabels
      )

      .merge(
        pureSoilLabels
      )

      .merge(
        mixedSoft
      )

      .merge(
        manualSamples
      );

  var training =
    featureTrain
      .sampleRegions({

        collection:
          labeled,

        properties: [

          'TreeRatio',
          'GrassRatio',
          'OtherRatio',
          'sample_type'

        ],

        scale:
          SCALE,

        geometries:
          false,

        tileScale:
          4

      })

      .filter(
        ee.Filter.neq(
          KEY_BAND,
          -9999
        )
      );

  /***************************************************************
   * MODULE 8: RANDOM-FOREST REGRESSION
   ***************************************************************/

  var logTraining =
    training.map(
      function(f) {

        var tree =
          ee.Number(
            f.get(
              'TreeRatio'
            )
          )
          .max(0)
          .add(
            LOGRATIO_EPS
          );

        var grass =
          ee.Number(
            f.get(
              'GrassRatio'
            )
          )
          .max(0)
          .add(
            LOGRATIO_EPS
          );

        var other =
          ee.Number(
            f.get(
              'OtherRatio'
            )
          )
          .max(0)
          .add(
            LOGRATIO_EPS
          );

        return f.set({

          zTree:
            tree
              .divide(other)
              .log(),

          zGrass:
            grass
              .divide(other)
              .log()

        });

      }
    );

  var rfTree =
    ee.Classifier
      .smileRandomForest({

        numberOfTrees:
          RF_TREES,

        variablesPerSplit:
          null,

        minLeafPopulation:
          1,

        bagFraction:
          0.7,

        maxNodes:
          null,

        seed:
          RF_SEED

      })

      .setOutputMode(
        'REGRESSION'
      )

      .train({

        features:
          logTraining,

        classProperty:
          'zTree',

        inputProperties:
          allBands

      });

  var rfGrass =
    ee.Classifier
      .smileRandomForest({

        numberOfTrees:
          RF_TREES,

        variablesPerSplit:
          null,

        minLeafPopulation:
          1,

        bagFraction:
          0.7,

        maxNodes:
          null,

        seed:
          RF_SEED + 1

      })

      .setOutputMode(
        'REGRESSION'
      )

      .train({

        features:
          logTraining,

        classProperty:
          'zGrass',

        inputProperties:
          allBands

      });

  var zTree =
    featureCity

      .classify(
        rfTree
      )

      .rename(
        'zTree'
      )

      .clamp(
        -15,
        15
      );

  var zGrass =
    featureCity

      .classify(
        rfGrass
      )

      .rename(
        'zGrass'
      )

      .clamp(
        -15,
        15
      );

  var rTree =
    zTree.exp();

  var rGrass =
    zGrass.exp();

  var rOther =
    ee.Image.constant(1);

  var denominator =
    rTree

      .add(
        rGrass
      )

      .add(
        rOther
      );

  var treeFraction =
    rTree

      .divide(
        denominator
      )

      .rename(
        'Tree'
      )

      .float();

  var grassFraction =
    rGrass

      .divide(
        denominator
      )

      .rename(
        'Grass'
      )

      .float();

  var otherFraction =
    rOther

      .divide(
        denominator
      )

      .rename(
        'Other'
      )

      .float();

  var fractions =
    ee.Image.cat([

      treeFraction,

      grassFraction,

      otherFraction

    ])
    .clip(
      region
    );

  var fractionSum =
    fractions

      .reduce(
        ee.Reducer.sum()
      )

      .rename(
        'FractionSum'
      );

  /***************************************************************
   * MODULE 9: VISUALIZATION AND EXPORT
   ***************************************************************/

  Map.addLayer(

    mixedSoft,

    {
      color:
        'FFFF00'
    },

    'Reliable Local Fisher pseudo labels',

    false

  );

  Map.addLayer(

    treeFraction,

    {

      min:
        0,

      max:
        1,

      palette: [

        'ffffff',

        '006400'

      ]

    },

    'Tree fractional cover'

  );

  Map.addLayer(

    grassFraction,

    {

      min:
        0,

      max:
        1,

      palette: [

        'ffffff',

        '7CFC00'

      ]

    },

    'Grass fractional cover'

  );

  Map.addLayer(

    otherFraction,

    {

      min:
        0,

      max:
        1,

      palette: [

        'ffffff',

        '000000'

      ]

    },

    'Other fractional cover'

  );

  Map.addLayer(

    fractionSum,

    {

      min:
        0.999,

      max:
        1.001

    },

    'Tree + Grass + Other',

    false

  );

  var rgb =
    ee.Image.cat([

      otherFraction,

      treeFraction,

      grassFraction

    ]);

  Map.addLayer(

    rgb,

    {

      min: [
        0,
        0,
        0
      ],

      max: [
        1,
        1,
        1
      ]

    },

    cityName +
      ' fractional-cover RGB'

  );

  function exportToDrive(
    image,
    description,
    prefix
  ) {

    Export.image.toDrive({

      image:
        image
          .clamp(
            0,
            1
          )
          .toFloat(),

      description:
        description,

      folder:
        DRIVE_FOLDER,

      fileNamePrefix:
        prefix,

      region:
        region,

      scale:
        SCALE,

      maxPixels:
        1e13,

      fileFormat:
        'GeoTIFF',

      formatOptions: {

        cloudOptimized:
          true

      }

    });
  }

  exportToDrive(

    treeFraction,

    cityName +
      '_' +
      YEAR +
      '_Tree_LocalFisherRF',

    cityName +
      '_' +
      YEAR +
      '_Tree_LocalFisherRF'

  );

  exportToDrive(

    grassFraction,

    cityName +
      '_' +
      YEAR +
      '_Grass_LocalFisherRF',

    cityName +
      '_' +
      YEAR +
      '_Grass_LocalFisherRF'

  );

  exportToDrive(

    otherFraction,

    cityName +
      '_' +
      YEAR +
      '_Other_LocalFisherRF',

    cityName +
      '_' +
      YEAR +
      '_Other_LocalFisherRF'

  );

  exportToDrive(

    fractions,

    cityName +
      '_' +
      YEAR +
      '_TreeGrassOther_LocalFisherRF',

    cityName +
      '_' +
      YEAR +
      '_TreeGrassOther_LocalFisherRF'

  );

}

/***************************************************************
 * MODULE 10: USER INTERFACE AND EXECUTION
 ***************************************************************/

function runAllForUI(
  region,
  cityName
) {

  var masks =
    buildPureAndMixedMasks(
      region
    );

  var treePoints =
    stratifiedPointsSafe(

      masks.tree,

      1,

      N_TREE,

      region,

      'Tree'

    );

  var grassPoints =
    stratifiedPointsSafe(

      masks.grass,

      2,

      N_GRASS,

      region,

      'Grass'

    );

  var imperviousPoints =
    stratifiedPointsSafe(

      masks.impervious,

      3,

      N_IMPERVIOUS,

      region,

      'Impervious'

    );

  var soilPoints =
    stratifiedPointsSafe(

      masks.soil,

      4,

      N_SOIL,

      region,

      'Soil'

    );

  var mixedPoints =
    randomMixedPoints(

      masks.mixed,

      N_MIXED,

      region

    );

  var manualLoad =
    tryLoadManualSamples(
      MANUAL_SHP_ASSET
    );

  var manualSamples =
    ee.FeatureCollection(

      ee.Algorithms.If(

        manualLoad.has,

        normalizeManualLabels(

          manualLoad.fc
            .filterBounds(
              region
            )

        ),

        ee.FeatureCollection([])

      )
    );

  Map.addLayer(

    treePoints,

    {
      color:
        '006400'
    },

    'Pure Tree'

  );

  Map.addLayer(

    grassPoints,

    {
      color:
        '7CFC00'
    },

    'Pure Grass'

  );

  Map.addLayer(

    imperviousPoints,

    {
      color:
        '555555'
    },

    'Pure Impervious'

  );

  Map.addLayer(

    soilPoints,

    {
      color:
        '8B4513'
    },

    'Pure Soil'

  );

  Map.addLayer(

    mixedPoints,

    {
      color:
        'FFFF00'
    },

    'Mixed candidates'

  );

  Map.addLayer(

    manualSamples,

    {
      color:
        '00FFFF'
    },

    'Manual',

    false

  );

  runAnalysis(

    region,

    cityName,

    treePoints,

    grassPoints,

    imperviousPoints,

    soilPoints,

    mixedPoints,

    manualSamples

  );
}

var allCitiesFC =
  ee.FeatureCollection(
    ALL_CITIES_ASSET
  );

var CURRENT_CITY_NAME =
  null;

var CURRENT_REGION =
  null;

var uiPanel =
  ui.Panel({

    style: {

      width:
        '440px',

      position:
        'top-left',

      padding:
        '8px'

    }

  });

uiPanel.add(

  ui.Label({

    value:
      'Localized Fisher-FCLS + Compositional RF',

    style: {

      fontWeight:
        'bold',

      fontSize:
        '15px'

    }

  })

);

uiPanel.add(

  ui.Label({

    value:
      'Tree / Grass / Other fractional cover'

  })

);

var citySelect =
  ui.Select({

    placeholder:
      'Loading city list...'

  });

allCitiesFC
  .aggregate_array(
    'city'
  )
  .distinct()
  .sort()
  .evaluate(

    function(list) {

      citySelect
        .items()
        .reset(
          list
        );

      citySelect
        .setPlaceholder(
          'Select city'
        );

    }

  );

var yearBox =
  ui.Textbox({

    placeholder:
      'Year, default ' +
      YEAR,

    value:
      ''

  });

var confirmBtn =
  ui.Button({

    label:
      'Show boundary',

    style: {

      stretch:
        'horizontal'

    }

  });

var runBtn =
  ui.Button({

    label:
      'Run Local Fisher-FCLS + RF',

    style: {

      stretch:
        'horizontal'

    }

  });

uiPanel.add(

  ui.Panel(

    [

      ui.Label(
        'City'
      ),

      citySelect

    ],

    ui.Panel.Layout.flow(
      'horizontal'
    )

  )

);

uiPanel.add(

  ui.Panel(

    [

      ui.Label(
        'Year'
      ),

      yearBox

    ],

    ui.Panel.Layout.flow(
      'horizontal'
    )

  )

);

uiPanel.add(
  confirmBtn
);

uiPanel.add(
  runBtn
);

Map.add(
  uiPanel
);

confirmBtn.onClick(

  function() {

    var cityName =
      citySelect.getValue();

    if (
      !cityName
    ) {

      return;

    }

    var cityFC =
      allCitiesFC.filter(

        ee.Filter.eq(
          'city',
          cityName
        )

      );

    var region =
      cityFC
        .union(1)
        .geometry();

    Map.clear();

    Map.add(
      uiPanel
    );

    Map.centerObject(
      region,
      10
    );

    Map.addLayer(

      cityFC.style({

        color:
          'black',

        fillColor:
          '00000000',

        width:
          2

      }),

      {},

      cityName +
        ' boundary'

    );

    CURRENT_CITY_NAME =
      cityName;

    CURRENT_REGION =
      region;

  }

);

runBtn.onClick(

  function() {

    var cityName =
      CURRENT_CITY_NAME ||
      citySelect.getValue();

    if (
      !cityName
    ) {

      return;

    }

    var cityFC =
      allCitiesFC.filter(

        ee.Filter.eq(
          'city',
          cityName
        )

      );

    var region =
      CURRENT_REGION ||
      cityFC
        .union(1)
        .geometry();

    var yr =
      yearBox.getValue();

    if (
      yr &&
      yr.trim() !== ''
    ) {

      YEAR =
        parseInt(
          yr,
          10
        );

    }

    MONTHS =
      [
        1,
        12
      ];

    Map.clear();

    Map.add(
      uiPanel
    );

    Map.centerObject(
      region,
      10
    );

    Map.addLayer(

      cityFC.style({

        color:
          'black',

        fillColor:
          '00000000',

        width:
          2

      }),

      {},

      cityName +
        ' boundary'

    );

    runAllForUI(

      region,

      cityName

    );

  }

);
