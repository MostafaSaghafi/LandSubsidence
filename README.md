# LandSubsidence
Negative Changes in SAR Amplitude using Sentinel-1 and sentinel-2


```javascript
Map.addLayer(shp, {},"My Polygon");
```
Adds the imported shapefile `shp` to the map as a layer. The `{}` signifies an empty dictionary for visualization parameters (using defaults), and "My Polygon" sets the layer name.  This `shp` defines your Area of Interest (AOI).


```javascript
var slope_threshold = .5; // unit: degree
var curv_threshold = -0.005; // unit: m/m^2
```
Sets the threshold values for slope (0.5 degrees) and curvature (-0.005 m/m²) which will be used for masking later.


```javascript
//define pre-event stack time period
var PreEventTime_1 = '2015-01-01T23:59'; // format: yyyy-mm-dd-HH:MM
var PreEventTime_2 = '2018-06-29T23:59'; // format: yyyy-mm-dd-HH:MM

//define post-event stack time period
var PostEventTime_1 = '2018-07-09T23:59'; // format: yyyy-mm-dd-HH:MM
var PostEventTime_2 = '2023-09-29T23:59'; // format: yyyy-mm-dd-HH:MM
```
Defines the start (`_1`) and end (`_2`) dates for the pre- and post-event periods.  These will be used to filter the satellite image collections.


```javascript
// cloud filter and mask for Sentinel-2 optical data
function maskS2clouds(image) {
  var qa = image.select('QA60');  // Select the QA60 band (quality assessment)
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;  // Left-shift creates a bitmask for the cloud bit (2^10)
  var cirrusBitMask = 1 << 11; // Bitmask for cirrus (2^11)
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0) // Check if the cloud bit is 0
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0)); // Check if the cirrus bit is 0
  return image.updateMask(mask).divide(10000); // Apply the mask and scale reflectance (0-1)
}
```
This function takes a Sentinel-2 image and masks out cloudy and cirrus pixels using the QA60 band. The `bitwiseAnd` operation checks the specific bits associated with clouds and cirrus.  Finally, the reflectance values are divided by 10000 to scale them (likely from 0-10000 to 0-1).


```javascript
// cloud filter and mask for Landsat 8 optical data (not used in the script)
function maskL8sr(image) {
  // ... (similar logic as maskS2clouds, but for Landsat 8)
}
```
This function is defined but not used in the script. It would mask clouds and cloud shadows in Landsat 8 data if it were called.


```javascript
// Load Sentinel-2 reflectance data.
var S2_PreEvent = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate(PreEventTime_1,PreEventTime_2)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) //only include images with less than 10% clouds
                  .filterBounds(shp)
                  .map(maskS2clouds);

var S2_PostEvent = ee.ImageCollection('COPERNICUS/S2') // Same as above, but for post-event
                  // ...
```
Loads Sentinel-2 image collections for pre- and post-event periods.  It filters the collections based on date, cloud cover percentage (less than 10%), and the area of interest (`shp`). The `maskS2clouds` function is applied to each image in the collection.



```javascript
// LOAD Shuttle Radar Topography Mission (SRTM) Digital Elevation Model (DEM)
var SRTM_dataset = ee.Image('USGS/SRTMGL1_003');
var elevation = SRTM_dataset.select('elevation');
var slope = ee.Terrain.slope(elevation); //slope in degrees
var mask_slope = slope.gte(slope_threshold); // slope mask: 1 where slope >= threshold, 0 elsewhere
var slope_mask = slope.updateMask(mask_slope); // Apply the slope mask to the slope image

// Calculate curvature
var smooth_curv = ee.Kernel.gaussian({ // Gaussian kernel for smoothing the DEM
  // ... parameters for the kernel
});
var xyDemGrad = elevation.convolve(smooth_curv).gradient(); // Gradient of smoothed elevation
var xGradient = xyDemGrad.select('x').gradient(); // Gradient of the x-component
var yGradient = xyDemGrad.select('y').gradient(); // Gradient of the y-component
var curvature = xGradient.select('x').add(yGradient.select('y')); // Calculate curvature
var mask_curvature = curvature.gte(curv_threshold); // Curvature mask
var curvature_mask = curvature.updateMask(mask_curvature); // Apply curvature mask
```
Loads the SRTM DEM, calculates slope and curvature, and creates masks based on the defined thresholds.  The curvature calculation involves smoothing the DEM with a Gaussian kernel before computing the second derivatives.


```javascript
// Define a Gaussian kernel to reduce noise in S1 scenes
var smooth_S1 = ee.Kernel.gaussian({
  // ... (parameters for S1 smoothing kernel)
});

// LOAD Sentinel-1 (S1) amplitude data VH polarization
var imgVH = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select('VH')
        .filterBounds(shp)
        .map(function(image) {
          var edge = image.lt(-30.0); // Identify low edge values (<-30dB)
          var maskedImage = image.mask().and(edge.not()); // Combine existing mask and edge removal
          return image.updateMask(maskedImage);
        });
```
Defines a smoothing kernel (`smooth_S1`) for Sentinel-1 data (though it's not used in the current calculations). Loads Sentinel-1 data, filtering by polarization ('VH'), instrument mode ('IW'), and the AOI. It then applies a function to mask out low edge values (less than -30dB) which are often noisy.


```javascript
var desc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')); // Filter for descending passes
var asc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));  // Filter for ascending passes

var PreEventPeriod = ee.Filter.date(PreEventTime_1,PreEventTime_2);  // Date filter for pre-event
var PostEventPeriod = ee.Filter.date(PostEventTime_1,PostEventTime_2); // Date filter for post-event

// calculate the median of the Pre-event S1 SAR Amplitude
var PreEventPeriod_asc = asc.filter(PreEventPeriod).median();  // Median of pre-event ascending
var PreEventPeriod_desc = desc.filter(PreEventPeriod).median(); // Median of pre-event descending

// calculate the median of the Post-event S1 SAR Amplitude
var PostEventPeriod_asc = asc.filter(PostEventPeriod).median(); // Median of post-event ascending
var PostEventPeriod_desc = desc.filter(PostEventPeriod).median();// Median of post-event descending


//Commented-out code (smoothing) - not active in this version of the script
// ...
```
Filters Sentinel-1 data into ascending (`asc`) and descending (`desc`) orbits. Creates date filters and calculates the median amplitude image for pre- and post-event periods for both ascending and descending orbits. Contains commented-out lines for optional smoothing, which are not active in the current code.


```javascript
// print out image information (number of images, image file name)
var num_asc_pre = asc.filter(PreEventPeriod).filterBounds(shp); // Pre-event ascending image collection
// ... (similar for other combinations of pre/post and asc/desc)

var count_asc_pre = num_asc_pre.sort('system:time_start').toList(5000,0).length(); // Count pre-event ascending images
// ... (similar counts for other combinations)

print("Number of pre-event ascending images: ", count_asc_pre);
// ... (print statements for other counts)

//Commented out print statements showing image collection info
// ...
```
Calculates and prints the number of images in each pre/post and asc/desc combination.  The `.toList(5000, 0)` method converts the collection to a list (limited to 5000 elements) to enable counting.  There are also some commented-out lines for printing more details about the image collections.

```javascript
// calculate the log ratio (using subtraction since data are in log scale) for Pre- and Post-event S1 SAR Amplitude
var A_ratio_desc = PreEventPeriod_desc.subtract(PostEventPeriod_desc); // Descending amplitude difference
var A_ratio_asc = PreEventPeriod_asc.subtract(PostEventPeriod_asc);   // Ascending amplitude difference
var A_ratio_avg_desc_asc = (A_ratio_asc.add(A_ratio_desc)).divide(2); // Average amplitude difference

////// Create a mask to retain only the negative values
var negativeMask = A_ratio_avg_desc_asc.lt(0); // Mask: 1 where average difference < 0, 0 elsewhere

//////// Visualize only the negative values
Map.addLayer(A_ratio_avg_desc_asc.updateMask(negativeMask).clip(shp), ColorScale, 'Negative Amplitude Change', true);

// define color palette
// ... (color palettes defined for various visualizations)
var ColorScale = {min: -10, max: 10, palette: ['0013ff','8178ff','ffffff','ff7e7e','ff0000']}; // for amplitude change (A_ratio)
var SlopeColorScale = {/* ... */}; // for slope
var AmplitudeColorScale = {/* ... */}; //for amplitude
var ColorCurv = {/* ... */}; // for curvature
var rgbVis = {/* ... */}; // for sentinel-2
```
Calculates the difference (log ratio) between pre- and post-event SAR amplitude for both ascending and descending orbits, then averages these differences. It creates a mask to highlight negative changes and adds a layer to the map visualizing these masked changes using `ColorScale`.  Color palettes are defined for visualizing various data.


```javascript
Map.centerObject(shp, 10); //zooms to center of AOI. The number determines the zoom level.


//// Sentinel-2 images pre_event
Map.addLayer(S2_PreEvent.median().clip(shp), rgbVis , 'S2 pre-event',true);
// ... (Export to Drive)

//// Sentinel-2 images post_event
Map.addLayer(S2_PostEvent.median().clip(shp), rgbVis, 'S2 post-event',true);
// ... (Export to Drive)

// ... (Map layers and exports for pre/post and asc/desc amplitude stacks)

Map.addLayer(A_ratio_asc.clip(shp), ColorScale, 'Amplitude ratio (asc)', true);
// ... (Export to Drive)

Map.addLayer(A_ratio_desc.clip(shp), ColorScale, 'Amplitude ratio (desc)', true);
// ... (Export to Drive)

Map.addLayer(curvature.clip(shp), ColorCurv,'SRTM curvature',false);
// ... (Export to Drive)


//Commented-out map layers - not active in this version of the script
// ...

Map.addLayer(A_ratio_avg_desc_asc.updateMask(mask_slope).updateMask(mask_curvature).clip(shp), ColorScale, 'S1 SAR amplitude change w/mask', true);


// Export final results to Drive
Export.image.toDrive({/* ... parameters for SAR_amplitude_change ... */});
Export.image.toDrive({/* ... parameters for SAR_amplitude_change_Threshold ... */});

```
Centers the map view on the AOI. Adds layers for visualization of pre- and post-event Sentinel-2, pre- and post-event ascending and descending Sentinel-1 amplitude, ascending and descending amplitude ratio, and curvature.  It exports these to Google Drive. Includes several commented-out sections that were present in earlier code versions but are not active now.  The active code includes a masked average amplitude change visualization and exports of the masked amplitude change and the mask itself.
