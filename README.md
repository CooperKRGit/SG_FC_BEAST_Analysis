# SG_FC_BEAST_Analysis
Project files for the Spring Gully fractional cover analysis using the BEAST algorithm
---
title: "MajorProjectScript"
author: "CKR"
date: '2022-07-21'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## The following script was provided by Shawn Laffan and adapted

## Load necessary packages

```{r}

library("terra")
library("sf")
library("Rbeast")
library("stringr")
library("zoo")
library("lubridate")

```


```{r}
##  A set of rasters from 1987-2022 was obtained from the TERN portal below:
##  https://portal.tern.org.au/seasonal-fractional-cover-australia-coverage/22026

files = list.files(path="spring_gully", pattern=".tif$", full.names = TRUE)

target_band = 2  #  2 is green fraction

# bands were separated from each file and sorted by date
bands_by_date = list(
  `1` = list(),
  `2` = list(),
  `3` = list(),
  `4` = list()
)

for (f in files) {
  r = rast(f)
  # message(str_extract(f, "_m[0-9]{6}"))
  date = ymd(paste(str_extract(f, "[0-9]{6}"), "01"))
  # message(date)
  for (i in 1:4) {
    bands_by_date[[i]][[as.character(date)]] = r[[i]]
  }
}

rstacks = list()
for (i in 1:4) {
  rstacks[[i]] = rast(bands_by_date[[i]])
}

#  we processed these above in the loop but it is simpler to do them again
dates = ymd(paste(str_extract(files, "[0-9]{6}"), "01"))

#  Point locations - much easier with a point feature class loaded via sf
xlocs = c( 1751997, 1752065, 1752256, 1752186, 1752281, 1752102, 1751912, 1751138, 1751368, 1751878)
ylocs = c( -3852258, -3852333, -3852504, -3852484, -3852601, -3852470, -3852346, -3852129, -3852571, -3852677)
# The first 4 are in the swampland of Spring Gully itself, the next 3 are of the tree canopy surrounding the swampland of Spring Gully, the last three are from vegetated areas further from Spring Gully but still within the Spring Gully SFAZ region defined in the RNP Fire Management Strategy.

```


# Run the coordinates through the BEAST algorithm

```{r}

vals = terra::extract (rstacks[[target_band]], data.frame(x=xlocs, y=ylocs), list=FALSE)
#  make the ID val the row name
row.names(vals) = vals[,1]
vals = vals[,-1]  #  drop ID col
vals = t(vals)    #  transpose

#  now beastify each point
#  and collate each in lists
beast_results = list()
ts_objects    = list()
for (i in 1:ncol(vals)) {
  #  generic time series object
  ts = zoo(vals[,i], order.by = dates)
  plot(ts)
  ts_objects[[i]] = ts
  
  #  no season as these data have already been aggregated by season
  b = beast(vals[,i], time = dates, season='none')
  b$time = dates
  #Ensures dates are displayed along x-axis
  plot(
     b, 
     index = 1,
     vars  = c('st','s','scp','sorder','t','tcp','torder','slpsgn','o','ocp','error'),  
     col         = NULL, 
     main        = "BEAST decomposition and changepoint detection",
     xlab        = 'Time',
     ylab        = NULL,
     cex.main    = 1,
     cex.lab     = 1,  
     relative.heights = NULL,           
     interactive = FALSE,
     ncpStat     = c('median','mode','mean','pct90','max'),)

  beast_results[[i]] = b
}

```
