human_pop<-raster("0.Data/GHS_population/GHS_POP_E2015_GLOBE_R2019A_54009_1K_V1_0.tif")
settlements<-raster("0.Data/GHS_population/GHS_SMOD_POP2015_GLOBE_R2019A_54009_1K_V2_0.tif")

settlements_rural<-settlements %in% c(11, 12, 13)

human_rural<-human_pop*settlements_rural

writeRaster(human_rural, "0.Data/GHS_population/Calculated_human_rural_population.tif")
