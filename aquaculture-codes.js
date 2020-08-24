/*
  If you run the entire code, it will report an out-of-memory error.
  Therefore, intermediate results need to be exported to run the code in sections.
*/	
var params = {
    'cloudcover': 20,
    'startDate': '2018-01-01',
    'endDate': '2018-12-31',
    'geometry': ee.FeatureCollection("users/xiazilong123/important/aquaculture/Shanghai-border")
};
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
var s2 = ee.ImageCollection("COPERNICUS/S2"),
    s1 = ee.ImageCollection("COPERNICUS/S1_GRD")

function clipImage(image){
    image=image.clip(params.geometry);
    return image;
}

var getPerimeter = function(image){     
  var a1 = ee.Kernel.fixed(3, 3,[[0,-1,0],[0,1,0],[0,0,0]]);
  var a2 = ee.Kernel.fixed(3, 3,[[0,0,0],[-1,1,0],[0,0,0]]);
  var a3 = ee.Kernel.fixed(3, 3,[[0,0,0],[0,1,-1],[0,0,0]]);
  var a4 = ee.Kernel.fixed(3, 3,[[0,0,0],[0,1,0],[0,-1,0]]);
  var r1 = image.convolve(a1).where(image.convolve(a1).neq(0),1);
  var r2 = image.convolve(a2).where(image.convolve(a2).neq(0),1);
  var r3 = image.convolve(a3).where(image.convolve(a3).neq(0),1);
  var r4 = image.convolve(a4).where(image.convolve(a4).neq(0),1);
  var result=r1.add(r2).add(r3).add(r4);
  return result;
};
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
var SentinelFunctions = {
  applyNDWI: function(image) {
    var ndwi = image.normalizedDifference(['green','nir']);
    return ndwi.select([0], ['ndwi']);
  },
  applyNDVI: function(image) {
    var ndvi = image.normalizedDifference(['nir','red']);
    return ndvi.select([0], ['ndvi']);
  }
};

s2 = s2.filterDate(params.startDate,params.endDate)
                  .filterBounds(params.geometry)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', params.cloudcover))
                  .map(maskS2clouds)
                  .select(['B2','B3','B4','B8','B11','B12'],
                       [ 'blue', 'green', 'red', 'nir','swir1', 'swir2'])
                  .map(clipImage);
s1 = s1.filterDate(params.startDate,params.endDate)
                  .filterBounds(params.geometry)
                  .map(clipImage);

var NDWI = s2.map(SentinelFunctions.applyNDWI).reduce(ee.Reducer.intervalMean(75, 100), 10);
var NDWIstdDev = s2.map(SentinelFunctions.applyNDWI).reduce(ee.Reducer.stdDev(), 10);
var NDVI = s2.map(SentinelFunctions.applyNDVI).reduce(ee.Reducer.intervalMean(75, 100), 10);
var NDVIstdDev = s2.map(SentinelFunctions.applyNDVI).reduce(ee.Reducer.stdDev(), 10);
var SAR = s1.reduce(ee.Reducer.intervalMean(75, 100), 10).select('VH_mean')

//Generate object's features at different thresholds
function getIndex(threshold,image){ 
    var name=threshold.toString().replace(".","-");
    var water=image.gt(threshold);
    var bandName = water.bandNames().get(0)
    var area=water.addBands(water).reduceConnectedComponents(ee.Reducer.sum(),bandName,2000).unmask(-999).rename('area'+name);
    var perimeterPixels = getPerimeter(water).int().rename('perimeterPixels');
    var perimeter = perimeterPixels.addBands(water).reduceConnectedComponents(ee.Reducer.sum(),bandName,2000).unmask(-999).rename('perimeter'+name);
    var ndvi = NDVI.addBands(water).reduceConnectedComponents(ee.Reducer.mean(),bandName,2000).unmask(-999).rename('ndvi'+name);
    var ndvistddev = NDVIstdDev.addBands(water).reduceConnectedComponents(ee.Reducer.mean(),bandName,2000).unmask(-999).rename('ndvistddev'+name);
    var ndwi = NDWI.rename('ndwi').addBands(water).reduceConnectedComponents(ee.Reducer.mean(),bandName,2000).unmask(-999).rename('ndwi'+name);
    var ndwistddev = NDWIstdDev.addBands(water).reduceConnectedComponents(ee.Reducer.mean(),bandName,2000).unmask(-999).rename('ndwistddev'+name);
    var vhmean = SAR.addBands(water).reduceConnectedComponents(ee.Reducer.mean(),bandName,2000).unmask(-999).rename('vhmean'+name);
    var LSI = perimeter.multiply(0.25).divide(area.sqrt()).unmask(-999).rename('LSI'+name);  //形状指数
    var sizes = ee.Image.pixelLonLat().addBands(water)
        .reduceConnectedComponents(ee.Reducer.minMax(),bandName, 2000);
    var width = sizes.select('longitude_max').subtract(sizes.select('longitude_min')).rename('width');
    var height = sizes.select('latitude_max').subtract(sizes.select('latitude_min')).rename('height');
    var ratio = height.divide(width).unmask(-999).rename('ratio'+name);
    var compact = width.multiply(height).divide(area).unmask(-999).rename('compact'+name);
    return area.addBands(perimeter).addBands(ndvi).addBands(ndvistddev).addBands(ndwi).addBands(ndwistddev).addBands(vhmean)
                     .addBands(LSI).addBands(ratio).addBands(compact);
}

var objectPropertiesImage = ee.Image.cat([
      getIndex(0,NDWI),getIndex(0.05,NDWI),getIndex(0.1,NDWI),getIndex(0.15,NDWI),getIndex(0.2,NDWI),
      getIndex(0.25,NDWI),getIndex(0.3,NDWI),
]).float();
var training = ee.FeatureCollection("users/xiazilong123/important/aquaculture/training");
var bands = objectPropertiesImage.bandNames();
var classifier = ee.Classifier.randomForest(400).train(training, 'class', bands)
var classification = objectPropertiesImage.classify(classifier);  
classification=classification.mask(NDWI.gt(0.1)).unmask(0).clip(params.geometry);

/*  
   Post-classification
*/
/*  
   Although our results in the article only show three categories: 1.aquaculture ponds 2.Other water bodies
   0.Non-water. 
   Actually we have 6 types of water bodies in the final classification results(before integrated updating): 0.Non-water 1.aquaculture ponds
   2. wetland 3.lake 4.river 5.shadows from buildings (6.ocean) 
   However, because the focus of this article is not on the classification of water bodies, 
   we did not consider the typicality of non-aquaculture pond ground samples, so the accuracy of other classifications is not very high. 
   If you need water classification, please reset your own samples.
*/
var area_fish=classification.eq(1).connectedPixelCount(1024,false)
                .where(classification.neq(1),classification.eq(1));
var area_lake=classification.eq(3).connectedPixelCount(1024,false)
                .where(classification.neq(3),classification.eq(3));
var area_river=classification.eq(4).connectedPixelCount(1024,false)
                .where(classification.neq(4),classification.eq(4));   			
var ocean=classification.eq(2).addBands(classification.eq(2).rename('ocean')).reduceConnectedComponents(ee.Reducer.sum(),'ocean',2000).unmask(-999)
          .eq(-999).where(classification.neq(2),classification.eq(2)); //We just need to extract the approximate extent of ocean 
var river_to_lake= area_river.lt(12).and(area_river.gt(0)).add(area_lake.gt(45)).connectedPixelCount(1024,false)
                   .neq(area_river).and(area_river.lt(12).and(area_river.gt(0))).multiply(3);		  
var fish_to_lake= area_fish.lt(12).and(area_fish.gt(0)).add(area_lake.gt(45)).connectedPixelCount(1024,false)
                   .neq(area_fish).and(area_fish.lt(12).and(area_fish.gt(0))).multiply(3);
var fish_to_lake= area_fish.lt(12).and(area_fish.gt(0)).add(area_lake.gt(45)).connectedPixelCount(1024,false)
                   .neq(area_fish).and(area_fish.lt(12).and(area_fish.gt(0))).multiply(3);                 
var fish_to_river= area_fish.lt(12).and(area_fish.gt(0)).add(area_river.gt(20)).connectedPixelCount(1024,false)
                   .neq(area_fish).and(area_fish.lt(12).and(area_fish.gt(0))).multiply(4);                   
var lake_to_fish =  area_lake.lt(12).and(area_lake.gt(0)).add(area_fish.gt(20)).connectedPixelCount(1024,false)
                   .neq(area_lake).and(area_lake.lt(12).and(area_lake.gt(0))).multiply(1);
var river_to_fish =  area_river.lt(12).and(area_river.gt(0)).add(area_fish.gt(20)).connectedPixelCount(1024,false)
                   .neq(area_river).and(area_river.lt(12).and(area_river.gt(0))).multiply(1);
var fish_to_ocean =  area_fish.lt(20).and(area_fish.gt(0)).add(ocean).connectedPixelCount(1024,false)
                   .neq(area_fish).and(area_fish.lt(12).and(area_fish.gt(0))) .multiply(6) ;                                      
var result =  classification.where(fish_to_lake,fish_to_lake) .where(river_to_lake,river_to_lake)   
                 .where(lake_to_fish,lake_to_fish).where(river_to_fish,river_to_fish)   
                 .where(fish_to_ocean,fish_to_ocean) .where(ocean,ocean.multiply(6));  
/*
    Integrated updating

*/				 
function reclass(img){ 
	return img=img.remap([0,1,2,3,4,5,6],[0,1,2,2,2,2,0]);
}	
/*
  result~result4 are classification results of different years obtained in the previous step.
*/	
var reference =ee.ImageCollection([reclass(result),reclass(result2),reclass(result3),reclass(result4)]).mode();
function Eliminate(img){
	var area_fish=img.eq(1).connectedPixelCount(1024,false)
                .where(classification.neq(1),classification.eq(1));
	var area_nonfish=classification.eq(2).connectedPixelCount(1024,false)
                .where(classification.neq(2),classification.eq(2)); 
	var fish_to_nonfish= area_fish.lt(20).and(area_fish.gt(0)).add(area_nonfish.gt(45)).connectedPixelCount(1024,false)
                   .neq(area_fish).and(area_fish.lt(12).and(area_fish.gt(0))).multiply(2);
	var nonfish_to_fish= area_nonfish.lt(12).and(area_nonfish.gt(0)).add(area_fish.gt(20)).connectedPixelCount(1024,false)
                   .neq(area_nonfish).and(area_nonfish.lt(12).and(area_nonfish.gt(0))).multiply(1);
	return img.where(fish_to_nonfish,fish_to_nonfish).where(nonfish_to_fish,nonfish_to_fish); 
}
reference =Eliminate(reference)
var new_result = reclass(result).where(result.neq(0).and(reference.neq(0)),reference).where(result.eq(6),result.eq(6).multiply(2))
new_result =Eliminate(new_result)
			 
				 
				 
				 

