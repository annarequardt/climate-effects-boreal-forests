// !! important: make sure to import the ROI first (Krycklan-geometry) and rename it to geometry in GEE

var START = '2018-01-01'
var END = '2019-12-31'
var date_range = ee.Filter.date(START, END)

//THIS CODE WORKS WELL, SET CLOUD PERCENTAGE, TIME PERIOD AND IMPORT THE ROI OF THE STUDY AREA. 
//THEN RUN, AND YOU WILL ONLY HAVE TO CLICK DOZENS OF RUN IN THE TASKS
//Change Only herei
var S2_cloud_percent = 50

var threshold_cloudFraction = 20 //this is in percentage, and will be used for filtering after calculating cloud percentage in the ROI
var threshold_MaskedFraction_MOD = 50 //this is a percentage to exclude MODIS images with too many missing pixels


//Make sure to also import the polygon of the ROI and rename "table" to "geometry" 
//Change not needed for the rest!


//COMMON FUNCTIONS and VARIABLES
var batch = require('users/fitoprincipe/geetools:batch');
  
 var date_extractor = function(image) {
  //var date = ee.Date(image.get('system:time_start')).format("YYYY-MM-dd");
  var timeStamp = ee.Date(image.get('system:time_start'));
  var date = timeStamp.format("YYYY-MM-dd");
  timeStamp = timeStamp.format("YYYYMMddHHmmss");
  return image.set('timeStamp', timeStamp)
                .set('date', date);
  };
  
  var date_extractor_S2 = function(image) {
  //var date = ee.Date(image.get('system:time_start')).format("YYYY-MM-dd");
  var timeStamp = ee.Date(image.get('system:time_start'));
  var date = timeStamp.format("YYYY-MM-dd");
  timeStamp = timeStamp.format("YYYYMMddHHmmss");
  var Sun_Zenith = ee.String(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(10000).round().toInt());
  
  var Sun_Zenith_B2 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B2')).multiply(10000).round().toInt());
  var Sun_Zenith_B3 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B3')).multiply(10000).round().toInt());
  var Sun_Zenith_B4 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B4')).multiply(10000).round().toInt());
  var Sun_Zenith_B8 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8')).multiply(10000).round().toInt());
  
  var Sun_Azimuth = ee.String(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(10000).round().toInt());
  
  var Sun_Azimuth_B2 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B2')).multiply(10000).round().toInt());
  var Sun_Azimuth_B3 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B3')).multiply(10000).round().toInt());
  var Sun_Azimuth_B4 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B4')).multiply(10000).round().toInt());
  var Sun_Azimuth_B8 = ee.String(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8')).multiply(10000).round().toInt());
  return image.set('timeStamp', timeStamp)
                .set('date', date)
                .set('Sun_Zenith_Azimuth', Sun_Zenith
                      .cat('_').cat(Sun_Zenith_B2)
                      .cat('_').cat(Sun_Zenith_B3)
                      .cat('_').cat(Sun_Zenith_B4)
                      .cat('_').cat(Sun_Zenith_B8)
                      .cat('_').cat(Sun_Azimuth).cat('_')
                      .cat(Sun_Azimuth_B2).cat('_')
                      .cat(Sun_Azimuth_B3).cat('_')
                      .cat(Sun_Azimuth_B4).cat('_')
                      .cat(Sun_Azimuth_B8));
  };
  
  var renamer_MOD = function(image){
    var name = ee.String('MODIS_').cat(image.get('timeStamp'));
    return image.rename(name);
  };
  
  var id_resetter_MOD_ALBEDO = function(image){
    var name = ee.String('MODIS_ALBE_').cat(image.get('timeStamp'));
    return image.set('system:id', name);
  };
  
  var id_resetter_MOD_BRDF = function(image){
    var name = ee.String('MODIS_BRDF_').cat(image.get('timeStamp'));
    return image.set('system:id', name);
  };
  
  var id_resetter_S2 = function(image){
    var name = ee.String('S2_').cat(image.get('timeStamp')).cat('_').cat(image.get('Sun_Zenith_Azimuth'));
    return image.set('system:id', name)
                .set('id', name);
  };
  


/*var date_extractor_mod = function(image) {
  var index = image.getString('system:index');
  var date = index.slice(0,4).cat("-").cat(index.slice(5,7)).cat("-").cat(index.slice(8,10));
  date = ee.Date(date);
    return image.set('date', date);
  };
  
var date_extractor_s2 = function(image) {
  var index = image.getString('system:index');
  var date = index.slice(0,4).cat("-").cat(index.slice(4,6)).cat("-").cat(index.slice(6,8));
  date = ee.Date(date);
    return image.set('date', date);
  };*/
  


//Function for correcting the reflectance values with the scale factor for MODIS
function corrector_MOD(image) {
  return image.divide(1000)
              .set('timeStamp', image.get('timeStamp'))
              .set('date', image.get('date'));
}


//Function for correcting the reflectance values with the scale factor for S2
function corrector_S2(image) {
  return image.divide(10000)
              .set('timeStamp', image.get('timeStamp'))
              .set('date', image.get('date'))
              .set('Sun_Zenith_Azimuth', image.get('Sun_Zenith_Azimuth'));
}


/////
//Function for detecting which image has no masked pixel
function noMaskedCalc(image){
  var mask = image.select(['B2'])
                  .unmask(-99)
                  .eq(-99)
  
  var AllPixelsDict = mask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry,
    scale: 10,
    maxPixels: 3784216672400
  });
  
  var AllPixels = AllPixelsDict.get('B2');
  
  
  var MaskMasked = mask.selfMask();//Masks 0 pixels, i.e. non cloud, then we can count the rest, which are clouds
  
  var CloudPixelsDict = MaskMasked.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry,
    scale: 10,
    maxPixels: 3784216672400
  });
  
  var CloudPixels = CloudPixelsDict.get('B2')
  
  var fraction = (ee.Number(CloudPixels).divide(AllPixels))
  .multiply(100);
  return image.set('CloudFraction', fraction);
}

///
//Function for detecting which image has no masked pixel in MODIS Albedo images
function noMaskedCalcModAlbe(image){
  var mask = image.select(['Albedo_BSA_Band1'])
                  .unmask(-99)
                  .eq(-99)
  
  var AllPixelsDict = mask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry,
    scale: 500,
    maxPixels: 3784216672400
  });
  
  var AllPixels = AllPixelsDict.get('Albedo_BSA_Band1');
  
  
  var MaskMasked = mask.selfMask();//Masks 0 pixels, i.e. non cloud, then we can count the rest, which are clouds
  
  var MaskedPixelsDict = MaskMasked.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry,
    scale: 500,
    maxPixels: 3784216672400
  });
  
  var MaskedPixels = MaskedPixelsDict.get('Albedo_BSA_Band1')
  
  var fraction = (ee.Number(MaskedPixels).divide(AllPixels))
  .multiply(100);
  return image.set('MaskedFraction', fraction);
}



////
//Function for detecting which image has no masked pixel in MODIS BRDF images
function noMaskedCalcModBrdf(image){
  var mask = image.select(['BRDF_Albedo_Parameters_Band1_iso'])
                  .unmask(-99)
                  .eq(-99)
  
  var AllPixelsDict = mask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry,
    scale: 500,
    maxPixels: 3784216672400
  });
  
  var AllPixels = AllPixelsDict.get('BRDF_Albedo_Parameters_Band1_iso');
  
  
  var MaskMasked = mask.selfMask();//Masks 0 pixels, i.e. non cloud, then we can count the rest, which are clouds
  
  var MaskedPixelsDict = MaskMasked.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry,
    scale: 500,
    maxPixels: 3784216672400
  });
  
  var MaskedPixels = MaskedPixelsDict.get('BRDF_Albedo_Parameters_Band1_iso')
  
  var fraction = (ee.Number(MaskedPixels).divide(AllPixels))
  .multiply(100);
  return image.set('MaskedFraction', fraction);
}
///////////////////////////////////////////////////
//////////////////////////////////MOD ALBEDO images
///////////////////////////////////////////////////

var mod_albedo_dataset = ee.ImageCollection('MODIS/061/MCD43A3')
                  .filter(date_range)
                  .filterBounds(geometry)
                  .map(date_extractor);
                  
                  
/*print(mod_dataset);*/
var mod_albedo_dataset_sub = mod_albedo_dataset
                        .select([
                          'Albedo_BSA_Band1', 'Albedo_BSA_Band2', 'Albedo_BSA_Band3', 'Albedo_BSA_Band4',
                          'Albedo_WSA_Band1', 'Albedo_WSA_Band2', 'Albedo_WSA_Band3', 'Albedo_WSA_Band4'
                        ])
                        .map(corrector_MOD)
                        .map(id_resetter_MOD_ALBEDO)
                        .map(noMaskedCalcModAlbe)
                        .filter(ee.Filter.lt('MaskedFraction', threshold_MaskedFraction_MOD));


print('MODIS', mod_albedo_dataset_sub);
                        

///////////////////////////////////////////////////
//////////////////////////////////MOD BRDF images
///////////////////////////////////////////////////
var mod_brdf_dataset = ee.ImageCollection('MODIS/061/MCD43A1')
                  .filter(date_range)
                  .filterBounds(geometry)
                  .map(date_extractor);



var mod_brdf_dataset_sub = mod_brdf_dataset.select([
  'BRDF_Albedo_Parameters_Band1_iso', 'BRDF_Albedo_Parameters_Band1_geo', 'BRDF_Albedo_Parameters_Band1_vol', 
  'BRDF_Albedo_Parameters_Band2_iso', 'BRDF_Albedo_Parameters_Band2_geo', 'BRDF_Albedo_Parameters_Band2_vol',
  'BRDF_Albedo_Parameters_Band3_iso', 'BRDF_Albedo_Parameters_Band3_geo', 'BRDF_Albedo_Parameters_Band3_vol', 
  'BRDF_Albedo_Parameters_Band4_iso', 'BRDF_Albedo_Parameters_Band4_geo', 'BRDF_Albedo_Parameters_Band4_vol'
]).map(corrector_MOD)
  .map(id_resetter_MOD_BRDF)
  .map(noMaskedCalcModBrdf)
  .filter(ee.Filter.lt('MaskedFraction', threshold_MaskedFraction_MOD));
  
  
print('MODIS BRDF', mod_brdf_dataset_sub);
                        
///////////////////////////////////////////////////
//////////////////////////////////S2 images
///////////////////////////////////////////////////

/**
 * Function to mask clouds using the Sentinel-2 QA band
 * @param {ee.Image} image Sentinel-2 image
 * @return {ee.Image} cloud masked Sentinel-2 image
 */
/*function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask)
              .divide(10000)
              .set('timeStamp', image.get('timeStamp'))
              .set('date', image.get('date'));
}*/



//!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!
var AOI = geometry
var START_DATE = START
var END_DATE = END
var CLOUD_FILTER = 100
var CLD_PRB_THRESH = S2_cloud_percent
var NIR_DRK_THRESH = 0.15
var CLD_PRJ_DIST = 1
var BUFFER = 50


function get_s2_sr_cld_col(aoi, start_date, end_date) {
    // # Import and filter S2 SR.
    var s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER)))

    // # Import and filter s2cloudless.
    var s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    // # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals({
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))
}


function add_cloud_bands(img) {
    // # Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    // # Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    // # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))
    
}


function add_shadow_bands(img) {
    // # Identify water pixels from the SCL band.
    var not_water = img.select('SCL').neq(6)

    // # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    // # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    // # Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    // # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
}


function add_cld_shdw_mask(img) {
    // # Add cloud component bands.
    var img_cloud = add_cloud_bands(img)

    // # Add cloud shadow component bands.
    var img_cloud_shadow = add_shadow_bands(img_cloud)

    // # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    // # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    // # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject({'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))

    // # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)
}



function apply_cld_shdw_mask(img) {
    // # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    var not_cld_shdw = img.select('cloudmask').not()

    // # Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)
}

var s2_sr_cld_col = get_s2_sr_cld_col(AOI, START_DATE, END_DATE)
var s2_dataset = s2_sr_cld_col.map(add_cld_shdw_mask)
                             .map(apply_cld_shdw_mask)
                             .filterBounds(geometry)
                             .map(date_extractor_S2);

//!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!

/*var maskUnwantedPixels= function(img){
    var info = img.select('SCL');
    var mask = info.neq(1).or(info.neq(3).or(info.neq(8).or(info.neq(9).or(info.neq(10).or(info.neq(11))))));//read here for the meaning of each bit value https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED#bands
    img = img.updateMask(mask);
    return img;
}*/


/*var s2_dataset = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filter(date_range)
                  .filterBounds(geometry)
                  .map(date_extractor)
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',50))
                  .map(maskS2clouds);
                
*/
print('SentinelUnfiltered', s2_dataset);

var s2_dataset_sub = s2_dataset
                    .select(['B2','B3','B4','B8'])
                    .map(corrector_S2)
                    .map(id_resetter_S2)
                    .map(noMaskedCalc)
                    .filter(ee.Filter.lt('CloudFraction', threshold_cloudFraction));

print('Sentinel', s2_dataset_sub);


///////////////////////////////
///////////JOINING THE 2 COLLECTIONS
var filterTimeEq = ee.Filter.equals({
    leftField: 'date',
    rightField: 'date'
    });
//SENTINEL

var simpleJoin = ee.Join.simple()
var simpleJoinedCollection = simpleJoin.apply(s2_dataset_sub, mod_albedo_dataset_sub, filterTimeEq)
 

var collection_S2= ee.ImageCollection(simpleJoinedCollection)

print('selection', collection_S2)


//MODIS ALBEDO

var simpleJoin = ee.Join.simple()
var simpleJoinedCollection = simpleJoin.apply(mod_albedo_dataset_sub, collection_S2, filterTimeEq)
 

var collection_MOD_albedo= ee.ImageCollection(simpleJoinedCollection)

print('selection MODIS ALBEDO', collection_MOD_albedo)


//MODIS BRDF

var simpleJoin = ee.Join.simple()
var simpleJoinedCollection = simpleJoin.apply(mod_brdf_dataset_sub, collection_S2, filterTimeEq)
 

var collection_MOD_brdf= ee.ImageCollection(simpleJoinedCollection)

print('selection MODIS BRDF', collection_MOD_brdf)




//https://gis.stackexchange.com/questions/364758/renaming-images-of-imagecollection-to-add-band-name-at-end-on-google-earth-engin

/*batch.Download.ImageCollection.toDrive(collection_S2), '140binclip'+bands[i],   
                                         {scale: 10, region: geometry})*/
                                         
                                         
batch.Download.ImageCollection.toDrive(collection_S2, 'ANNA_S2_ALL', 
  {name: '{system:id}',
  scale: 10,
  maxPixels: 3784216672400,
  region: geometry,
  crs:'EPSG:3006'
  }
)




batch.Download.ImageCollection.toDrive(collection_MOD_brdf, 'ANNA_MODIS_BRDF_ALL', 
  {name: '{system:id}',
  scale: 500,
  maxPixels: 3784216672400,
  region: geometry,
  crs:'EPSG:3006'
  }
)