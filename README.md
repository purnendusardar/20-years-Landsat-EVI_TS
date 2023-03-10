# Pixel-Based Dense Time Series Analysis of LANDSAT EVI for 20-years
The time series stack of EVI was prepared using cloud-free Landsat 5 TM, Landsat 7 ETM+, and Landsat 8 OLI Surface Reflectance data from GEE Landsat Archive. The images were masked using the quality assessment function and then EVI for each month was stacked in a single raster image for further trend analysis of phenological stability

# Analysis Approach For The Final Map 

The analysis approach involves the following steps:

- Creating an image stack of the Landsat EVI data for each month using Earth Engine
- I used Landsat 5, 7 and Landsat 8 Surface Reflectance Data to generate the stack.
- I then used R to perform a pixel-based time series analysis using the *rtk* package in R
- The non-parametric Mann-Kendall test and Theil-Sen (TS) regression methods were used.
- The TS calculates the median of the paired slope throughout the time series. 
- I then visualized the results using ArcGIS Online maps.


# And ! Here is the final map! :metal:
![Map](/Result/Sample_map.png)

