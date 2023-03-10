#load packages 

require(raster)
require(rkt)
require(trend)
require(zoo)
require(tidyr)

rasterOptions(maxmemory=1e+06, chunksize=1e+07, progress = 'text')

# Set working directory (location of raster)
setwd("D:/Local Disk E/PhD_Data/Precipitation")
results <-"D:/Local Disk E/PhD_Data/Precipitation/results"

# stack all files with .tif extension in working directory

stacked_files <-stack(list.files(pattern = ".tif"))

#brick the stacked files

bricked_files <-brick(stacked_files)

#remove unnecessary files from memory
rm(stacked_files)
gc()

plot(bricked_files)

# Set year range for analysis
nlayers(bricked_files)
years =1:nlayers(bricked_files)
length(years)

# Analysis function
rktFun <-function(x) {
  if(all(is.na(x))){
    c(NA,NA,NA)
  } else {
    analysis <-rkt(years, x)
    a <-analysis$B # theil sen slope
    b <-analysis$sl # pvalue
    c <-analysis$tau # tau
    return(cbind(a, b, c))
  } }
rRaster <-calc(bricked_files, rktFun)


# Write to Results folder
writeRaster(rRaster[[1]], paste0(results,"/ts_slope.tif"), overwrite=T)
writeRaster(rRaster[[2]], paste0(results,"/mk_pvalue.tif"), overwrite=T)
writeRaster(rRaster[[3]], paste0(results,"/mk_tau.tif"), overwrite=T)


plot(rRaster)

# Load tau raster
tau <-raster(paste0(results, "/mk_tau.tif"))# Set p-values to create masks for
p_values <-c(0.01, 0.05, 0.1)
# Loop through p-values, producing a mask for each
for(i in 1:length(p_values)){
# Load the P-value raster
  p_value_raster <-raster(paste0(results, "/mk_pvalue.tif"))
  # Select current p-value
  p_val <-p_values[[i]]
  # Create string vers for filenaming
  p_val_str <-gsub("\\.", "", as.character(p_val))
  # Mask
  p_value_raster[p_value_raster > p_val] <-NA
  p_masked <-mask(tau, p_value_raster)
  # Write result
  writeRaster(p_masked, paste0(results, "/pvalue_mask", p_val_str, ".tif"),
              overwrite = TRUE)
  # Cleanup
}


plot(p_masked)

## 6. P-Value Masking with Significant Tau ##

# Isolate tau values > 0.4 and < -0.4
sigTau <-raster(paste0(results, "/mk_tau.tif"))
sigTau[sigTau>(-0.4) & sigTau<0.4] <-NA
# Loop through p-values, producing a mask for each
for(i in 1:length(p_values)){
  # Select current p-value
  p_val <-p_values[[i]]
  # Create string vers for filenaming
  p_val_str <-gsub("\\.", "", as.character(p_val))
  # Read current p-value raster
  p_value_raster <-raster(paste0(results, "/pvalue_mask", p_val_str,".tif"))
  # Mask significant tau with p-value raster
  tau_masked <-mask(sigTau, p_value_raster) # Write result
  writeRaster(tau_masked, paste0(results, "/tau_mask", p_val_str, ".tif"),
              overwrite = TRUE)
}

plot(tau_masked)


##Identify the turning points##
vals <- values(bricked_files)
head(vals)

# create empty matrix to store results
res <- matrix(nrow=nrow(vals), ncol=2)

# create empty matrix to store results
res <- matrix(nrow=nrow(vals), ncol=2)


# start looping through the pixels
for (i in 1:nrow(vals)){
  # get time series of first pixel
  x <- vals[i,]
  # check if there is data in the pixel
  if(all(is.na(x))){
    # if not, save NA, NA as result
    res[i,] <- c(NA,NA)
  } else {
    # if there is data, fill data gaps in the time series using simple interpolation
    x1 <- na.approx(na.approx(x))
    # apply the pettitt test
    analysis <- pettitt.test(x1)
    # extract the results
    a <-as.numeric(analysis$estimate)[1] # pettitt test (id at which time step the change occurred)
    b <-analysis$p.value
    # save the results
    res[i,] <- cbind(a, b)
  }
  # print current iteration
  print(i)
}

# get a single raster band of the time series (does not matter which one)

pettitt.p.val <- bricked_files[[1]]

# overwrite the values of the raster band with the p-values obtained for the pettitt-test
# this step only works because our results file has exactly the same number rows as the
# raster has pixels

values(pettitt.p.val) <- res[,2]

# plot the resulting raster

plot(pettitt.p.val)

# copy the results concerning turning points in a new variable

turnp <- res[,1]

# change each of the id values to the corresponding year

for (i in 1:length(years)){
  turnp[turnp==i] <- years[i]
  print(i)
}

# get a single raster band of the time series (does not matter which one)

pettitt.timep <- bricked_files[[1]]

# overwrite the pixel values with the results

values(pettitt.timep) <- turnp

# plot the turning point raster

plot(pettitt.timep)

mask <- pettitt.p.val < 00.5
plot(mask)

timep_fin <- mask(pettitt.timep, mask, maskvalue=0, updatevalue=NA)
plot(timep_fin)