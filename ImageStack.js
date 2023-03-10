//Time Series EVI calculation//
// This code is a modified version of Pironkova et al. 2018//

//--------------------------------------------//
// 1. QA FUNCTIONS //
//--------------------------------------------//

//SATURATES MASK (L5/7)
var maskSaturate57 = function(image) {
    image = image.select('B3', 'B4', 'B1');
    var mask = image.lt(10000).and(image.gte(0));
    return image.updateMask(mask);
  }
  
  //SATURATES MASK (L8)
  var maskSaturate8 = function(image) {
    image = image.select('B4', 'B5', 'B2');
    var mask = image.lt(10000).and(image.gte(0));
    return image.updateMask(mask);
  }
  
  //PIXEL MASK (L8)
  var maskPixels8 = function(image) {
    var pixel_qa = image.select('pixel_qa');
    var mask = pixel_qa.eq(322);
    return image.updateMask(mask);
  }
  
  //PIXEL MASK (L5/L7)
  var maskPixels57 = function(image) {
    var pixel_qa = image.select('pixel_qa');
    var mask = pixel_qa.eq(66);
    return image.updateMask(mask);
  }
  
  //AEROSOL MASK (L8)
  var maskAerosol = function(image) {
    var aero_qa = image.select('sr_aerosol');
    var mask = aero_qa.neq(194).and(aero_qa.neq(224)).and(aero_qa.neq(160)).and(aero_qa.neq(130));
    return image.updateMask(mask);
  }
  
  //ATMOS OPACITY MASK (L5/L7) -to get rid of haze
  var maskHaze = function(image) {
    //Select band and multiply by scaling factor
    var atmos_qa = image.select('sr_atmos_opacity').multiply(0.0010);
    //Mask for non-hazy pixels and remove fill vals
    var mask = atmos_qa.lte(0.1).and(atmos_qa.gt(-9.9));
    return image.updateMask(mask);
  }
  
  //----------------------------------------------//
  // 2. EVI FUNCTIONS //
  //----------------------------------------------//
  
  //EVI FUNCTIONS FOR LANDSAT 5/7//
  var getEVI57 = function(image) {
    //Apply saturate function
    image = maskSaturate57(image);
    var nir = image.select('B4');
    var red = image.select('B3');
    var blue = image.select('B1');
    var evi = ((nir.subtract(red)).multiply(2.5)).divide((red.multiply(6)).add(nir).subtract(blue.multiply(7.5).add(1))).rename('nd');
    return (evi);
  }
  
  //EVI FUNCTIONS FOR LANDSAT 8//
  var getEVI8 = function(image) {
    //Apply saturate function
    image = maskSaturate8(image);
    var nir = image.select('B5');
    var red = image.select('B4');
    var blue = image.select('B2');
    var evi = ((nir.subtract(red)).multiply(2.5)).divide((red.multiply(6)).add(nir).subtract(blue.multiply(7.5).add(1))).rename('nd');
    return (evi);
  }
  //---------------------------------------------// 
  // 3. SENSOR CALLIBRATION FUNCTIONS // 
  //---------------------------------------------//
  var applyCoefficientL5 = function(image) {
    return image.select('nd').multiply(1.036);
  }
  
  var applyCoefficientL8 = function(image) {
    return image.select('nd').multiply(1.0863);
  }
  
  //---------------------------------------------// 
  // 4. MAIN FUNCTION // 
  //---------------------------------------------//
  var createEVIComposite = function(){ 
    //Set year range 
    // !!!! Adjust the years you want to cover here !!!! //
    var yearrangeStart = 2002; 
    var yearrangeStop = 2022; 
  
    var listfin = ee.Image([]);
    for(var loopYear = yearrangeStart; loopYear <= yearrangeStop; loopYear +=1){ 
      
      //Set year's date range //
      // !!!! Adjust the time period you want to cover within each year here !!!! //
  
      var start = ee.Date.fromYMD(loopYear, 6, 20); 
      var end = ee.Date.fromYMD(loopYear, 9, 25); 
    }
  }
  
  // Landsat 8
  // Loads the Landsat 8 surface reflection collection and applies all pre-processing
  // functions, including filtering by AOI, date, and cloud cover, applying pixel and
  // aerosol masks, calculating EVI, and applying sensor calibration coefficients.
  var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
    // Filter by AOI (imported at the top)
    .filterBounds(SBR)
    // Filter by date
    .filterDate(start, end)
    // Filter by CLOUD_COVER
    .filterMetadata('CLOUD_COVER', 'less_than', 80)
    // Apply pixel mask
    .map(maskPixels8)
    // Apply aerosol mask
    .map(maskAerosol)
    // Calculate EVI
    .map(getEVI8)
    // Apply sensor calibration coefficient
    .map(applyCoefficientL8);
  
  
  // Landsat 7
  // load Landsat 7 surface reflection collection and apply all
  // pre-processing functions
  // and calculate EVI
  
  var l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
    // Filter AOI (imported at the top)
    .filterBounds(SBR)
    // Filter by date
    .filterDate(start, end)
    // Filter CLOUD_COVER
    .filterMetadata('CLOUD_COVER', 'less_than', 100)
    // Apply Pixel Mask
    .map(maskPixels57)
    // Apply Haze Mask
    .map(maskHaze)
    // Calculate EVI
    .map(getEVI57);
  
  // Landsat 5
  // Loads the Landsat 5 surface reflection collection and applies all pre-processing
  // functions, including filtering by AOI, date, and cloud cover, calculating EVI, and
  // applying sensor calibration coefficients.
  var l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
    // Filter by AOI (imported at the top)
    .filterBounds(SBR)
    // Filter by date 
    .filterDate(start, end)
    // Filter by CLOUD_COVER
    .filterMetadata('CLOUD_COVER', 'less_than', 100)
    // Apply pixel mask
    .map(maskPixels57)
    // Apply haze mask
    .map(maskHaze)
    // Calculate EVI
    .map(getEVI57)
    // Apply sensor calibration coefficient
    .map(applyCoefficientL5);
  
  //Merge collections
  var mergedCollection = ee.ImageCollection(l8.merge(l7).merge(l5));
  
  //Create composite with median, clip to AOI and rename band based on year
  var finalOutput = mergedCollection.reduce(ee.Reducer.median()).clip(SBR).rename(loopYear.toString());
  
  //Generate filename for export
  var filename = ("EVI_Composite_").concat(loopYear.toString());
  
  // add mosaic of current year to image list
  listfin = ee.Image(listfin).addBands(finalOutput.rename(filename));
  
  //end of loop
  
  }
  // print out the stacked image list to check how many mosaics are contained
  print(listfin)
  
  // convert all datasets to the same data type (LS8 has another datatype thatn LS5 and 7)
  var listfin2 = listfin.toDouble()
  
  
  // Export the stack of EVI mosaics to Google Drive
  
  Export.image.toDrive({
    image: listfin2,
    description: 'EVI_Ts',
    //Landsat resolution is 30m
    scale: 500,
    region: SBR,
    //LCC_MNRF projection
    crs: 'epsg:4326',
    maxPixels: 800000000
  })