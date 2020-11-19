# Quantifying Degree of Urbanisation

# Loading in packagees

library('easypackages')
library("devtools")
#devtools::install_github("16EAGLE/getSpatialData")
#devtools::install_github("vqv/ggbiplot")
#devtools::install_github("joshuakrobertson/brmsMethods")

libraries("brms", "rgdal", "gdalUtils", "RStoolbox", "getSpatialData", "rasterVis", "mapview", "loo",
"RColorBrewer", "grDevices", "gridExtra", "caret", "randomForest", "ranger", "MLmetrics", "nnet", "NeuralNetTools", "LiblineaR", "snow", "raster", "sf", "sp", "data.table", "tidyverse", "stringr", "doParallel", "parallel", "stars")

# Settomg directory for Sentinel-2 raster data

rasterOptions(tmpdir = "./rasterdump/cache")
login_CopHub(username = "joshuakrobertson", password = "BCCH2020!")

# Creating areas of interest and looping through Copernicus extraction

Erin = c(-80.129509 + 0.025, -80.129509 - 0.025, 43.730421 + 0.025, 43.730421 - 0.025); Guelph = c(-80.239789 + 0.025, -80.239789 - 0.025, 43.538196 + 0.025, 43.538196 - 0.025); Corwhin = c(-80.109827 + 0.025, -80.109827 - 0.025, 43.542421 + 0.025, 43.542421 - 0.025); Brantford = c(-80.302769 + 0.025, -80.302769 - 0.025, 43.133720 + 0.025, 43.133720 - 0.025);Cambridge = c(-80.291915 + 0.025, -80.291915 - 0.025, 43.381646 + 0.025, 43.381646 - 0.025); Ruthven = c(-79.872100 + 0.025, -79.872100 - 0.025, 42.968942 + 0.025, 42.968942 - 0.025)

Coordinates = list(Erin, Guelph, Corwhin, Brantford, Cambridge, Ruthven)
names(Coordinates) = c("Erin", "Guelph", "Corwhin", "Brantford", "Cambridge", "Ruthven")

# Creating box matrices in list form.

Coord_Matrices = vector("list", 6)

for (i in 1:length(Coord_Matrices)){
    Coord_Matrices[[i]] = matrix(data = c(Coordinates[[i]][1], Coordinates[[i]][3],
                       Coordinates[[i]][2], Coordinates[[i]][3],
                       Coordinates[[i]][2], Coordinates[[i]][4], 
                       Coordinates[[i]][1], Coordinates[[i]][4],
                       Coordinates[[i]][1], Coordinates[[i]][3]),
              ncol = 2, byrow = TRUE)
}

names(Coord_Matrices) = names(Coordinates)

# Looping through matrices and downloading corresponding Sentinel-2 rasters

for (i in 1:length(Coord_Matrices)){
    Dest = paste0("./rasterdump/Sentinel2/", names(Coord_Matrices[i]), "/")
    set_archive(Dest)

    set_aoi(Coord_Matrices[[i]])

    records = getSentinel_query(time_range = c("2018-04-01",
        as.character(Sys.Date())), platform = "Sentinel-2")

    records_filtered = records %>% 
    filter(level == "Level-2A") %>% 
    filter(cloudcov <= 1)

    datasets = getSentinel_data(records = records_filtered[1,])
    unzip(zipfile = list.files(paste0(Dest, "_datasets/Sentinel-2"), full.names = TRUE)[1], exdir = paste0(Dest, "_datasets/Sentinel-2"))
}

# Looping through collection of river, forest, building, street, and cropland data per region

Rivers = vector("list", 6); Roads = vector("list", 6); Forests = vector("list", 6); Buildings = vector("list", 6); Cropland = vector("list", 6) 

Source_List = list("/home/joshk/rasterdump/Rivers/ghy_000c11a_e.shp", "/home/joshk/rasterdump/Roads/lrnf000r16a_e.shp", "/home/joshk/rasterdump/Buildings/odb_ontario.shp", "/home/joshk/rasterdump/Forests/Wooded_Area.shp", "/home/joshk/rasterdump/Land_Use/IMG_AAFC_LANDUSE_Z17_2010.tif")

names(Source_List) = c("Rivers", "Roads", "Buildings", "Forests", "Cropland")

for (j in 1:length(Source_List)){
    for (i in 1:length(Coord_Matrices)){
        if (i == 1 & j %in% c(1,2,4)){
            Shape_Dat = st_read(Source_List[[j]])
        } else if (i == 1 & j == 3){
            Erin_Build = read.csv("/home/joshk/rasterdump/Buildings/Erin_Buildings_2.csv")
            Transitions = c(which(is.na(Erin_Build$X)))
            IDs = c()

            for (k in 1:length(Transitions)){
                if (k == 1){
                    IDs = c(IDs, rep(k, Transitions[k] - 1))
                } else if (k > 1){
                    IDs = c(IDs, rep(k, Transitions[k] - Transitions[k-1]))
                }
            }

            IDs = c(IDs, rep(length(Transitions) + 1, nrow(Erin_Build) - length(IDs)))
            Erin_Build = Erin_Build %>% 
            mutate("IDs" = IDs) %>% 
                dplyr::select(X, Polygon_0, IDs) %>% 
                mutate(X = as.numeric(as.character(X)), 
                Polygon_0 = as.numeric(as.character(Polygon_0)))

            Erin_Build_Simp = subset(Erin_Build, !is.na(X))
            EB_List = split(Erin_Build_Simp, f = Erin_Build_Simp$IDs)

            for (k in 1:length(EB_List)){
                EB_List[[k]] = EB_List[[k]]
                names(EB_List[[k]]) = c("long", "lat", "Polygon")
            } 

            EB_SP = vector("list", length(EB_List))

            for (k in 1:length(EB_List)){
                EB_SP[[k]] = SpatialPolygons(list(Polygons(list(Polygon(EB_List[[k]][,c(1,2)])),1)))
                proj4string(EB_SP[[k]]) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
            }

            EB_Collected = do.call(bind, EB_SP)
            Shape_Dat = st_as_sf(EB_Collected)
        } else if (i == 2 & j == 3){
            Shape_Dat = st_read(Source_List[[j]])
        } else if (i == 1 & j == 5){
            Shape_Dat = read_stars(Source_List[[j]])
        } 

        if (j %in% c(1:4)){

            Lin_Mat = c(Coord_Matrices[[i]][1,1], Coord_Matrices[[i]][2,1], Coord_Matrices[[i]][1,2], Coord_Matrices[[i]][3,2])
            xmin = Lin_Mat[2]; xmax = Lin_Mat[1]; ymin = Lin_Mat[4]; ymax = Lin_Mat[3]
    
            filter_func = sf::st_bbox(c(xmin = ifelse(is.na(xmin), -180, xmin), 
                ymin = ifelse(is.na(ymin),  -90,  ymin), 
                xmax = ifelse(is.na(xmax), +180, xmax), 
                ymax = ifelse(is.na(ymax),  +90, ymax)), 
                crs = st_crs(4269)) %>% 
            sf::st_as_sfc(.)

            if (j %in% c(2,3,4)){
                Shape_Dat = st_transform(Shape_Dat, 4269)
            }
            Trimmed_Data = sf::st_intersects(Shape_Dat, filter_func)
            To_Write = Shape_Dat[which(lengths(Trimmed_Data) != 0), ]

            if (j == 3 & i %in% c(2:6)){
                To_Write = st_zm(To_Write)    
            }

            st_write(To_Write, 
                paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(Coord_Matrices)[i], "_", names(Source_List)[j], ".shp"), driver = "ESRI Shapefile") 

            if (i == 6){
                rm(Shape_Dat)
            }
        } else if (j == 5){
                Lin_Mat = c(Coord_Matrices[[i]][1,1], Coord_Matrices[[i]][2,1], Coord_Matrices[[i]][1,2], Coord_Matrices[[i]][3,2])
                xmin = Lin_Mat[2]; xmax = Lin_Mat[1]; ymin = Lin_Mat[4]; ymax = Lin_Mat[3]
    
                filter_func = sf::st_bbox(c(xmin = ifelse(is.na(xmin), -180, xmin), 
                ymin = ifelse(is.na(ymin),  -90,  ymin), 
                xmax = ifelse(is.na(xmax), +180, xmax), 
                ymax = ifelse(is.na(ymax),  +90, ymax)), 
                crs = st_crs(4269)) %>% 
                sf::st_as_sfc(.)

                filter_func = st_transform(filter_func, 
                crs = '+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

                Trimmed_Data = st_crop(Shape_Dat, filter_func)
                Trimmed_Data = st_transform(Trimmed_Data, 
                    "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
                Shape_Data_Joined = st_as_sf(Trimmed_Data, as_points = FALSE, merge = FALSE)

                To_Write = Shape_Data_Joined %>% 
                    filter(IMG_AAFC_LANDUSE_Z17_2010.tif == 51) %>%
                    summarise(geometry = sf::st_union(geometry)) %>%
                    ungroup()

                st_write(To_Write, 
                paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(Coord_Matrices)[i], "_", names(Source_List)[j], ".shp"), driver = "ESRI Shapefile")         

                if (i == 6){
                    rm(Shape_Dat)
                }
        }
    }
}


# Reading in and compiling Sentinel-2 rasters

Base_Rasters = vector("list", 6)
Base_ggRasters = vector("list", 6)

for (j in 1:6){
    Base_fold = paste0("./rasterdump/Sentinel2/", names(Coord_Matrices[j]), "/_datasets/Sentinel-2/")
    jp2_folder = paste0(Base_fold, list.files(Base_fold, pattern = ".SAFE")[1], "/GRANULE/", list.files(paste0(Base_fold, "/", list.files(Base_fold, pattern = ".SAFE")[1], "/GRANULE"))[1], "/IMG_DATA/")

    jp2_10m = list.files(paste0(jp2_folder, "R10m"), pattern = ".*B.*jp2$", full.names = TRUE)
    jp2_20m = list.files(paste0(jp2_folder, "R20m"), pattern = ".*B.*[56781].*jp2$", full.names = TRUE)
    rst_lst = lapply(c(jp2_10m, jp2_20m), FUN = raster)

    extent_crop = Coord_Matrices[[j]] %>% 
        project(proj = proj4string(rst_lst[[1]])) %>% 
        apply(MARGIN = 2, FUN = range) %>%
        t %>%
        extent

    rst_crop_lst = lapply(rst_lst, FUN = raster::crop, y = extent_crop)
    names(rst_crop_lst) = str_extract(sapply(rst_crop_lst, names), "B.{2}")

    path_crop = paste0("/home/joshk/rasterdump/Sentinel2/", names(Coord_Matrices[j]), "/", names(rst_crop_lst), ".tif")

    for (i in 1:length(rst_crop_lst)){
    writeRaster(x = rst_crop_lst[[i]],
        filename = path_crop[i],
        format   = "GTiff",
        datatype = dataType(rst_crop_lst[[i]]),
        options  = "TFW=YES",
        overwrite = TRUE)
    }

    rst_for_prediction = vector(mode = "list", length = length(rst_crop_lst))
    names(rst_for_prediction) = names(rst_crop_lst)

    for (b in c("B05", "B06", "B07", "B8A", "B11", "B12")){
      rst_for_prediction[[b]] = raster::resample(x = rst_crop_lst[[b]],
                                                  y = rst_crop_lst$B02)
    }

    b_10m = c("B02", "B03", "B04", "B08")
    rst_for_prediction[b_10m] = rst_crop_lst[b_10m]
    brick_for_prediction = brick(rst_for_prediction)

    brick_for_prediction_norm = normImage(brick_for_prediction)
    names(brick_for_prediction_norm) = names(brick_for_prediction)

    Base_Rasters[[j]] = projectRaster(brick_for_prediction_norm, crs = '+proj=longlat +datum=WGS84')
    Base_ggRasters[[j]] = ggRGB(Base_Rasters[[j]])

    Base_Path = paste0("/home/joshk/rasterdump/Trimmed_Satellite/", names(Coord_Matrices[j]), "_Compiled.tif")

    Flat_Raster = merge(Base_Rasters[[j]], overlap = FALSE)

    writeRaster(x = Flat_Raster,
        filename = Base_Path,
        format   = "GTiff",
        datatype = dataType(Flat_Raster),
        options  = "TFW = YES",
        overwrite = TRUE)  
    }
}

# Compiling cropped shapefiles and writing to mosaic raster files 

Erin = vector("list", 5); Guelph = vector("list", 5); Corwhin = vector("list", 5); Brantford = vector("list", 5); Cambridge = vector("list", 5); Ruthven = vector("list", 5)

All_Shapes = list(Erin, Guelph, Corwhin, Brantford, Cambridge, Ruthven)
names(All_Shapes) = c("Erin", "Guelph", "Corwhin", "Brantford", "Cambridge", "Ruthven")

for (j in 1:6){
    for (i in 1:5){
        To_Read = paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(Coord_Matrices)[j], "_", names(Source_List)[i], ".shp")
        All_Shapes[[j]][[i]] = read_sf(To_Read)
    }
}

Sat_Images = vector('list', 6)

for (i in 1:6){
    Load_Sat = paste0("/home/joshk/rasterdump/Trimmed_Satellite/", names(Coord_Matrices[i]), "_Compiled.tif")
    Sat_Images[[i]] = raster(Load_Sat)
}

for (j in 1:6){

    # Laying in logical order

    C_base = All_Shapes[[j]][[5]] %>% 
    mutate("Type" = "Crop") %>%
    select(Type, geometry)
    C_base = rasterize(C_base, Sat_Images[[j]])
    C_base = projectRaster(C_base, crs = '+proj=longlat +datum=WGS84')
    C_vals = rep("1", length(C_base@data@values))
    C_vals[c(which(is.na(C_base@data@values)))] = "0"
    C_base@data@values = C_vals
    C_base@data@values = as.numeric(C_base@data@values)

    F_base = All_Shapes[[j]][[4]] %>% 
    mutate("Type" = "Forest") %>%
    select(Type, geometry)
    F_base = rasterize(F_base, Sat_Images[[j]])
    F_base = projectRaster(F_base, crs = '+proj=longlat +datum=WGS84')
    F_vals = rep("2", length(F_base@data@values))
    F_vals[c(which(is.na(F_base@data@values)))] = "0"
    F_base@data@values = as.numeric(F_vals)

    if (j %in% c(1:2,4:6)){
        R_base = All_Shapes[[j]][[1]] %>% 
        mutate("Type" = "Water") %>%
        select(Type, geometry)
        R_base = rasterize(R_base, Sat_Images[[j]])
        R_base = projectRaster(R_base, crs = '+proj=longlat +datum=WGS84')
        R_vals = rep("10", length(R_base@data@values))
        R_vals[c(which(is.na(R_base@data@values)))] = "0"
        R_base@data@values = as.numeric(R_vals)

        Mos1 = mosaic(C_base, F_base, R_base, fun = sum, na.rm = T)
        Mos1@data@values = as.integer(Mos1@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mos1@data@values)){
            if (Mos1@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mos1@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mos1@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 3){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 10){
                Mos_Vals[i] = 3
            } else if (Mos1@data@values[i] == 11){
                Mos_Vals[i] = 3
            } else if (Mos1@data@values[i] == 12){
                Mos_Vals[i] = 3
            }
        }
        Mos1@data@values = as.numeric(Mos_Vals)
    } else if (j == 3){
        Mos1 = mosaic(C_base, F_base, fun = sum, na.rm = T)
        Mos1@data@values = as.integer(Mos1@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mos1@data@values)){
            if (Mos1@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mos1@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mos1@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 3){
                Mos_Vals[i] = 2
            } 
        }
        Mos1@data@values = as.numeric(Mos_Vals)
    }

    S_base = All_Shapes[[j]][[2]] %>% 
    mutate("Type" = "Road") %>%
    select(Type, geometry)
    S_base = rasterize(S_base, Sat_Images[[j]])
    S_base = projectRaster(S_base, crs = '+proj=longlat +datum=WGS84')
    S_vals = rep("10", length(S_base@data@values))
    S_vals[c(which(is.na(S_base@data@values)))] = "0"
    S_base@data@values = as.numeric(S_vals)

    if (j %in% c(1,2,4,5)){
        B_base = All_Shapes[[j]][[3]] %>% 
        mutate("Type" = "Building") %>%
        select(Type, geometry)
        B_base = rasterize(B_base, Sat_Images[[j]])
        B_base = projectRaster(B_base, crs = '+proj=longlat +datum=WGS84')
        B_vals = rep("20", length(B_base@data@values))
        B_vals[c(which(is.na(B_base@data@values)))] = "0"
        B_base@data@values = as.numeric(B_vals)
  
        Mosaic_Full = mosaic(Mos1, S_base, B_base, fun=sum, na.rm=T)
        Mosaic_Full@data@values = as.integer(Mosaic_Full@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mosaic_Full@data@values)){
            if (Mosaic_Full@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mosaic_Full@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mosaic_Full@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mosaic_Full@data@values[i] == 3){
                Mos_Vals[i] = 3
            } else if (Mosaic_Full@data@values[i] == 10){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 11){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 12){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 13){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 20){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 21){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 22){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 23){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 30){
                Mos_Vals[i] = 5
            }
        }

        Mosaic_Full@data@values = as.numeric(Mos_Vals)

        writeRaster(x = Mosaic_Full,
        filename = paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[j]), ".tif"),
        format   = "GTiff",
        datatype = dataType(Mosaic_Full),
        options  = "TFW = YES",
        overwrite = TRUE)
        
    } else if (j %in% c(3,6)){
        Mosaic_Full = mosaic(Mos1, S_base, fun=sum, na.rm=T)
        Mosaic_Full@data@values = as.integer(Mosaic_Full@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mosaic_Full@data@values)){
            if (Mosaic_Full@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mosaic_Full@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mosaic_Full@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mosaic_Full@data@values[i] == 3){
                Mos_Vals[i] = 3
            } else if (Mosaic_Full@data@values[i] == 10){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 11){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 12){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 13){
                Mos_Vals[i] = 4
            } 
        }

        Mosaic_Full@data@values = as.numeric(Mos_Vals)

        writeRaster(x = Mosaic_Full,
        filename = paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[j]), ".tif"),
        format   = "GTiff",
        datatype = dataType(Mosaic_Full),
        options  = "TFW = YES",
        overwrite = TRUE)
    }
}

# Reading in written mosaics and plotting

Mosaics_List = vector("list", 6); Mosaics_List_Rasters = vector("list", 6); Mosaic_Plots = vector("list", 6)

for (i in 1:6){
    Mosaic_Path = paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[i]), ".tif")
    Mosaics_List_Rasters[[i]] = raster(paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[i]), ".tif"))
    Mosaics_List[[i]] = read_stars(Mosaic_Path)    
    Mosaics_List[[i]] = st_as_sf(Mosaics_List[[i]], as_points = FALSE, merge = TRUE)
    colnames(Mosaics_List[[i]])[1] = "Class"
    Mosaics_List[[i]]$Land_Use = "Unclassified"
    Mosaics_List[[i]]$Land_Use = as.character(Mosaics_List[[i]]$Land_Use)

    {
        Mosaics_List[[i]]$Land_Use[c(which(Mosaics_List[[i]]$Class == 0))] = "Unclassified"
        Mosaics_List[[i]]$Land_Use[c(which(Mosaics_List[[i]]$Class == 1))] = "Cropland"
        Mosaics_List[[i]]$Land_Use[c(which(Mosaics_List[[i]]$Class == 2))] = "Forest"
        Mosaics_List[[i]]$Land_Use[c(which(Mosaics_List[[i]]$Class == 3))] = "Water"
        Mosaics_List[[i]]$Land_Use[c(which(Mosaics_List[[i]]$Class == 4))] = "Street"
        Mosaics_List[[i]]$Land_Use[c(which(Mosaics_List[[i]]$Class == 5))] = "Building"
    }

    Mosaics_List[[i]]$Land_Use = factor(Mosaics_List[[i]]$Land_Use, levels = c("Unclassified","Building","Cropland","Forest","Street", "Water"))

    if (i %in% c(1,2,4,5)){
    Mosaic_Plots[[i]] = ggplot() + 
        geom_sf(data =  Mosaics_List[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("white","burlywood4","wheat","lightgreen","black","turquoise")) +
        theme_minimal()
    } else if (i == 3){
    Mosaic_Plots[[i]] = ggplot() + 
        geom_sf(data =  Mosaics_List[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("white","wheat","lightgreen","black")) +
        theme_minimal() 
    } else if (i == 6){
    Mosaic_Plots[[i]] = ggplot() + 
        geom_sf(data =  Mosaics_List[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("white","wheat","lightgreen","black","turquoise")) +
        theme_minimal()  
    }      
}

do.call("grid.arrange", c(Mosaic_Plots, ncol=3))

# Good. 

## Proceeding to preliminary data organisation for random forest algorithm. Note that this step is conducted to 
## classify buildings in remaining rural locations (here, Corwhin and Ruthven Park).

# First converting "polygon" rasters to points

Point_Frames = vector('list', 6)

# Option 1)

for (i in c(1,2,4,5)){
    Point_Rast = as.data.table(rasterToPoints(Mosaics_List_Rasters[[i]]))
    setnames(Point_Rast, old = names(Coord_Matrices)[i], new = "Class")

    Points_DF = SpatialPointsDataFrame(coords = Point_Rast[, .(x, y)],
    data = Point_Rast,
    proj4string = crs(Mosaics_List_Rasters[[i]]))

    # Pulling pixel values from Sentinel-2 images to overlay with "polygon" point values

    Class_Template = data.frame("Class" = c(0,1,2,3,4,5), 
    "Land_Use" = c("Unclassified", "Cropland", "Forest", "Water", "Street", "Building")
    )

    Proj_Rast = projectRaster(Base_Rasters[[i]], crs = crs(Points_DF))
    
    Point_Frames[[i]] = Proj_Rast %>% 
    raster::extract(y = Points_DF) %>%  
    as.data.frame %>% 
    mutate(Class = Points_DF@data$Class) %>% 
    left_join(Class_Template, by = c("Class")) %>% 
    mutate(Land_Use = factor(Land_Use)) %>% 
    select(-Class)

    # Quick histogram of bands

    # PFs %>% select(-Land_Use) %>%
    # melt(measure.vars = names(.)) %>% 
    # ggplot() +
    # geom_histogram(aes(value)) +
    # geom_vline(xintercept = 0, color = "gray70") +
    # facet_wrap(facets = vars(variable), ncol = 3)
}

# Option 2) Matching resolution of polygon to minimum resolution of Sentinel-2 image first

Point_Frames = vector('list', 6)

for (i in c(1,2,4,5)){
    Class_Template = data.frame("Class" = c(0,1,2,3,4,5), 
    "Land_Use" = c("Unclassified", "Cropland", "Forest", "Water", "Street", "Building")
    )

    Raster_Template = projectRaster(Base_Rasters[[i]], 
        crs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

    Raster_Template = raster(extent(Raster_Template$B02),
        resolution = 10,
        crs = projection(Raster_Template$B02))

    Raster_Template = projectRaster(Raster_Template, 
        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

    Poly_Matched = rasterize(Mosaics_List[[i]], Raster_Template, field = "Land_Use")
    Poly_DataTable = as.data.table(rasterToPoints(Poly_Matched))
    setnames(Poly_DataTable, old = "layer", new = "Land_Type")
    Poly_Points = SpatialPointsDataFrame(coords = Poly_DataTable[, .(x, y)],
        data = Poly_DataTable,
        proj4string = Poly_Matched@crs)

    Point_Frames[[i]] = Base_Rasters[[i]] %>% 
    raster::extract(y = Poly_Points) %>%  
    as.data.frame %>% 
    mutate(Land_Type = Poly_Points@data$Land_Type - 1) %>% 
    left_join(rename(Class_Template, "Land_Type" = Class), by = c("Land_Type")) %>% 
    mutate(Land_Use = factor(Land_Use)) %>% 
    select(-Land_Type)
}

# Creating training and test data with 75% applied to training
# First stripping out NAs. 

PFs = rbind(Point_Frames[[1]],Point_Frames[[2]],Point_Frames[[4]],Point_Frames[[5]])
PFs = na.omit(PFs)
PFs$Land_Use = factor(PFs$Land_Use) 

set.seed(101)

Train_Index = createDataPartition(PFs$Land_Use, p = 0.75, list = FALSE)
Train_Dat = PFs[Train_Index,]; Test_Dat = PFs[-Train_Index,]

table(Train_Dat$Land_Use); table(Test_Dat$Land_Use)

split_number = 10
set.seed(101)
Splits = createFolds(1:nrow(Train_Dat), k = split_number)
seeds = vector(mode = "list", length = split_number + 1)

for (i in 1:split_number){ 
    seeds[[i]] = sample.int(1000, split_number)
}
seeds[[split_number + 1]] = sample.int(1000, 1)

# Providing information for model control. Note that we are use cross-validation below (that is, "cv").

Mon_Control = trainControl(summaryFunction = multiClassSummary,
    method = "cv",
    number = split_number,
    search = "grid",
    classProbs = TRUE,
    savePredictions = TRUE,
    index = Splits,
    seeds = seeds)

# Trying neural network and SVM first 

nnet_grid = expand.grid(size = c(5, 10, 15),
                         decay = c(0.001, 0.01, 0.1))
svm_grid = expand.grid(cost = c(0.2, 0.5, 1),
                        Loss = c("L1", "L2"))

cluster_build = makeCluster(1/2 * detectCores()); registerDoParallel(cluster_build)

model_nnet = train(Land_Use ~ ., method = 'nnet', data = Train_Dat,
    importance = TRUE,
    maxit = 1000, 
    allowParallel = TRUE,
    tuneGrid = nnet_grid,
    trControl = Mon_Control)

model_svm = train(Land_Use ~ ., method = "svmLinear3", data = Train_Dat,
    allowParallel = TRUE,
    importance = TRUE,
    tuneGrid = svm_grid,
    trControl = Mon_Control)

stopCluster(cluster_build); rm(cluster_build)
registerDoSEQ()

# Summarising and saving model

model_svm$times$everything
saveRDS(model_svm, file = "/home/joshk/rasterdump/Random_Forest/SVM_Mod.rds")

# Constructing confusion matrix when contrasting with test data

Confused = confusionMatrix(data = predict(model_svm, newdata = Test_Dat),
    Test_Dat$Land_Use)

print(Confused)

# Very poor. Projecting onto map of Guelph

predict_svm = raster::predict(object = Base_Rasters[[2]], model = model_svm, type = 'raw')
viewRGB(brick(Base_Rasters[[2]]), r = 3, g = 2, b = 1) + mapView(predict_svm)

# Running random forests

#cluster_build = makeCluster(1/2 * detectCores()); registerDoParallel(cluster_build)

Random_Forest_Model = train(Land_Use ~ ., method = "rf", data = Train_Dat,
 #   allowParallel = TRUE,
    importance = TRUE,
    tuneGrid = data.frame(mtry = c(2,3,4,5,8)),
    trControl = Mon_Control)

stopCluster(cluster_build); rm(cluster_build)
registerDoSEQ()

# Viewing resultant model and viewing predictor importance for each land use type

Random_Forest_Model$finalModel
caret::varImp(Random_Forest_Model)$importance %>%
  as.matrix %>% 
  plot_ly(x = colnames(.), y = rownames(.), z = ., type = "heatmap",
          width = 350, height = 300)

###############################################################################################
# Simplifying polygon values to two levels (building and non-building) then re-running SVM.
# Note that Erin is excluded here temporarily then used to validate succes of model build with 
# urban data when predicting buildings in rural environments.
###############################################################################################

PF_Simple = rbind(Point_Frames[[2]],Point_Frames[[4]],Point_Frames[[5]])
PF_Simple = na.omit(PF_Simple)
PF_Simple$Land_Use = as.character(PF_Simple$Land_Use) 
PF_Simple$Land_Use[c(which(PF_Simple$Land_Use != "Building"))] = "Land"
PF_Simple$Land_Use = factor(PF_Simple$Land_Use) 

set.seed(101)

Train_Index = createDataPartition(PF_Simple$Land_Use, p = 0.75, list = FALSE)
Train_Dat = PF_Simple[Train_Index,]; Test_Dat = PF_Simple[-Train_Index,]

table(Train_Dat$Land_Use); table(Test_Dat$Land_Use)

split_number = 10
set.seed(101)
Splits = createFolds(1:nrow(Train_Dat), k = split_number)
seeds = vector(mode = "list", length = split_number + 1)

for (i in 1:split_number){ 
    seeds[[i]] = sample.int(1000, split_number)
}
seeds[[split_number + 1]] = sample.int(1000, 1)

Mon_Control = trainControl(summaryFunction = multiClassSummary,
    method = "cv",
    number = split_number,
    search = "grid",
    classProbs = TRUE,
    savePredictions = TRUE,
    index = Splits,
    seeds = seeds)

nnet_grid = expand.grid(size = c(5, 10, 15),
                         decay = c(0.001, 0.01, 0.1))

cluster_build = makeCluster(1/2 * detectCores()); registerDoParallel(cluster_build)

Simple_nnet = train(Land_Use ~ ., method = 'nnet', data = Train_Dat,
    importance = TRUE,
    maxit = 1000, 
    allowParallel = TRUE,
    tuneGrid = nnet_grid,
    trControl = Mon_Control)

stopCluster(cluster_build); rm(cluster_build)
registerDoSEQ()

Confused = confusionMatrix(data = predict(Simple_nnet, newdata = Test_Dat),
    Test_Dat$Land_Use)

print(Confused)

predict_nnet = raster::predict(object = Base_Rasters[[2]], model = Simple_nnet, type = 'raw')
viewRGB(brick(Base_Rasters[[2]]), r = 3, g = 2, b = 1) + mapView(predict_nnet)

# Much better! Saving model.

saveRDS(Simple_nnet, file = "/home/joshk/rasterdump/Random_Forest/Simple_NeuralNetwork_Model.rds")

# Collapsing point estimates into polygons and producing maps of true and estimated buildings for 
# visual comparisons.

Collapsed_Sample = rasterToPolygons(predict_nnet, dissolve = TRUE)
Collapsed_Sample = st_as_sf(Collapsed_Sample)
P1 = ggplot(data.frame("C" = c("Predicted", "Predicted"))) + 
    geom_sf(data = Collapsed_Sample, size = 0.5,
    aes(fill = factor(layer), colour = factor(layer))) + 
    scale_fill_manual(name = "Land\nCategorisation",
    labels = c("Building", "Other"), values = c("burlywood4", "white")) + 
    scale_colour_manual(name = "Land\nCategorisation", 
    labels = c("Building", "Other"), values = c("black", "white")) + 
    theme_bw() + facet_wrap(~C) + theme(panel.background = element_rect("white"),
    panel.grid.major = element_blank(), strip.text.x = element_text(size = 12)) + 
    xlab("Longitude") + ylab("Latitude")
    
Base_Sample = Mosaics_List[[2]] %>% 
    mutate(Class = as.character(Class), Land_Use = as.character(Land_Use))

Base_Sample$Class[c(which(Base_Sample$Land_Use == "Building"))] = "Building"
Base_Sample$Class[c(which(Base_Sample$Land_Use != "Building"))] = "Land"
Base_Sample = Base_Sample %>% mutate(Class = factor(Class)) %>% 
    rename("value" = Class) %>% 
    select(-Land_Use)
Base_Sample_Simple = as_Spatial(Base_Sample)
Base_Sample_Simple = aggregate(Base_Sample_Simple, by = list(Base_Sample_Simple$value), FUN="first") 
Base_Sample_Simple = st_as_sf(Base_Sample_Simple)

P2 = ggplot(data.frame("C" = c("True", "True"))) + 
    geom_sf(data = Base_Sample_Simple, size = 0.5, 
    aes(fill = factor(value), colour = factor(value))) + 
    scale_fill_manual(name = "Land\nCategorisation",
    labels = c("Building", "Other"), values = c("burlywood4", "white")) + 
    scale_colour_manual(name = "Land\nCategorisation",
    labels = c("Building", "Other"), values = c("black", "white")) + 
    theme_bw() + facet_wrap(~C) + theme(panel.background = element_rect("white"),
    panel.grid.major = element_blank(), strip.text.x = element_text(size = 12)) + 
    xlab("Longitude") + ylab("Latitude")

grid.arrange(P1, P2, ncol = 2)

# Removing left legend.

P1 = P1 + theme(legend.position = "none")

grid.arrange(P1, P2, ncol = 2)

#############################################
## Adding in Erin data and re-running model.
#############################################

PF_Expanded = rbind(Point_Frames[[3]],Point_Frames[[2]],Point_Frames[[4]],Point_Frames[[5]])
PF_Expanded = na.omit(PF_Expanded)
PF_Expanded$Land_Use = as.character(PF_Expanded$Land_Use) 
PF_Expanded$Land_Use[c(which(PF_Expanded$Land_Use != "Building"))] = "Land"
PF_Expanded$Land_Use = factor(PF_Expanded$Land_Use) 

set.seed(101)

Train_Index = createDataPartition(PF_Expanded$Land_Use, p = 0.75, list = FALSE)
Train_Dat = PF_Expanded[Train_Index,]; Test_Dat = PF_Expanded[-Train_Index,]

table(Train_Dat$Land_Use); table(Test_Dat$Land_Use)

split_number = 10
set.seed(101)
Splits = createFolds(1:nrow(Train_Dat), k = split_number)
seeds = vector(mode = "list", length = split_number + 1)

for (i in 1:split_number){ 
    seeds[[i]] = sample.int(1000, split_number)
}
seeds[[split_number + 1]] = sample.int(1000, 1)

nnet_grid = expand.grid(size = c(5, 10, 15),
                         decay = c(0.001, 0.01, 0.1))

Mon_Control = trainControl(summaryFunction = multiClassSummary,
    method = "cv",
    number = split_number,
    search = "grid",
    classProbs = TRUE,
    savePredictions = TRUE,
    index = Splits,
    seeds = seeds)

Expanded_nnet = train(Land_Use ~ ., method = 'nnet', data = Train_Dat,
    importance = TRUE,
    maxit = 1000, 
    allowParallel = TRUE,
    tuneGrid = nnet_grid,
    trControl = Mon_Control)

# Again, saving model

saveRDS(Expanded_nnet, file = "/home/joshk/rasterdump/Random_Forest/Expanded_NeuralNetwork_Model.rds")
#Expanded_nnet = readRDS(file = "/home/joshk/rasterdump/Random_Forest/Expanded_NeuralNetwork_Model.rds")
print(Expanded_nnet)

# Viewing accuracy and plotting

Confused = confusionMatrix(data = predict(Expanded_nnet, newdata = Test_Dat),
    Test_Dat$Land_Use)

print(Confused)

predict_nnet = raster::predict(object = Base_Rasters[[2]], model = Expanded_nnet, type = 'raw')
viewRGB(brick(Base_Rasters[[2]]), r = 3, g = 2, b = 1) + mapView(predict_nnet)

# Trying Corwhin

predict_nnet = raster::predict(object = Base_Rasters[[3]], model = Expanded_nnet, type = 'raw')
viewRGB(brick(Base_Rasters[[3]]), r = 3, g = 2, b = 1) + mapView(predict_nnet)

# Slight improvement. Grouping and plotting beside true values

Collapsed_Sample_Expanded = rasterToPolygons(predict_nnet, dissolve = TRUE)
Collapsed_Sample_Expanded = st_as_sf(Collapsed_Sample_Expanded)
P1 = ggplot(data.frame("C" = c("Predicted", "Predicted"))) + 
    geom_sf(data = Collapsed_Sample_Expanded, size = 0.5,
    aes(fill = factor(layer), colour = factor(layer))) + 
    scale_fill_manual(name = "Land\nCategorisation",
    labels = c("Building", "Other"), values = c("burlywood4", "beige")) + 
    scale_colour_manual(name = "Land\nCategorisation", 
    labels = c("Building", "Other"), values = c("black", "black")) + 
    theme_bw() + facet_wrap(~C) + theme(panel.background = element_rect("white"),
    panel.grid.major = element_blank(), strip.text.x = element_text(size = 12)) + 
    xlab("Longitude") + ylab("Latitude")
    
Base_Sample = Mosaics_List[[2]] %>% 
    mutate(Class = as.character(Class), Land_Use = as.character(Land_Use))

Base_Sample$Class[c(which(Base_Sample$Land_Use == "Building"))] = "Building"
Base_Sample$Class[c(which(Base_Sample$Land_Use != "Building"))] = "Land"
Base_Sample = Base_Sample %>% mutate(Class = factor(Class)) %>% 
    rename("value" = Class) %>% 
    select(-Land_Use)
Base_Sample_Simple = as_Spatial(Base_Sample)
Base_Sample_Simple = aggregate(Base_Sample_Simple, by = list(Base_Sample_Simple$value), FUN="first") 
Base_Sample_Simple = st_as_sf(Base_Sample_Simple)

P2 = ggplot(data.frame("C" = c("True", "True"))) + 
    geom_sf(data = Base_Sample_Simple, size = 0.5, 
    aes(fill = factor(value), colour = factor(value))) + 
    scale_fill_manual(name = "Land\nCategorisation",
    labels = c("Building", "Other"), values = c("burlywood4", "beige")) + 
    scale_colour_manual(name = "Land\nCategorisation",
    labels = c("Building", "Other"), values = c("black", "black")) + 
    theme_bw() + facet_wrap(~C) + theme(panel.background = element_rect("white"),
    panel.grid.major = element_blank(), strip.text.x = element_text(size = 12)) + 
    xlab("Longitude") + ylab("Latitude")

grid.arrange(P1, P2, ncol = 2)

# Still sparse, but acceptable if urban areas are scaled down. Replacing true building data with
# projected building data in complete mosaic and plotting. Note, however, that building layers cannot simply be "swapped" given that backgrounds may not be similar (that is, possible underlayers). Recompiling mosaics. Note that this requires emptying "All_Shapes" and "Sat_Images" lists. Furthermore, buildings are not permitted to exist on water.

Erin = vector("list", 5); Guelph = vector("list", 5); Corwhin = vector("list", 5); Brantford = vector("list", 5); Cambridge = vector("list", 5); Ruthven = vector("list", 5)

All_Shapes = list(Erin, Guelph, Corwhin, Brantford, Cambridge, Ruthven)
names(All_Shapes) = c("Erin", "Guelph", "Corwhin", "Brantford", "Cambridge", "Ruthven")

for (j in 1:6){
    for (i in 1:5){
        To_Read = paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(Coord_Matrices)[j], "_", names(Source_List)[i], ".shp")
        All_Shapes[[j]][[i]] = read_sf(To_Read)
    }
}

Sat_Images = vector('list', 6)

for (i in 1:6){
    Load_Sat = paste0("/home/joshk/rasterdump/Trimmed_Satellite/", names(Coord_Matrices[i]), "_Compiled.tif")
    Sat_Images[[i]] = raster(Load_Sat)
}

for (j in 1:6){

    C_base = All_Shapes[[j]][[5]] %>% 
    mutate("Type" = "Crop") %>%
    select(Type, geometry)
    C_base = rasterize(C_base, Sat_Images[[j]])
    C_base = projectRaster(C_base, crs = '+proj=longlat +datum=WGS84')
    C_vals = rep("1", length(C_base@data@values))
    C_vals[c(which(is.na(C_base@data@values)))] = "0"
    C_base@data@values = C_vals
    C_base@data@values = as.numeric(C_base@data@values)

    F_base = All_Shapes[[j]][[4]] %>% 
    mutate("Type" = "Forest") %>%
    select(Type, geometry)
    F_base = rasterize(F_base, Sat_Images[[j]])
    F_base = projectRaster(F_base, crs = '+proj=longlat +datum=WGS84')
    F_vals = rep("2", length(F_base@data@values))
    F_vals[c(which(is.na(F_base@data@values)))] = "0"
    F_base@data@values = as.numeric(F_vals)

    if (j %in% c(1:2,4:6)){
        R_base = All_Shapes[[j]][[1]] %>% 
        mutate("Type" = "Water") %>%
        select(Type, geometry)
        R_base = rasterize(R_base, Sat_Images[[j]])
        R_base = projectRaster(R_base, crs = '+proj=longlat +datum=WGS84')
        R_vals = rep("10", length(R_base@data@values))
        R_vals[c(which(is.na(R_base@data@values)))] = "0"
        R_base@data@values = as.numeric(R_vals)

        Mos1 = mosaic(C_base, F_base, R_base, fun = sum, na.rm = T)
        Mos1@data@values = as.integer(Mos1@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mos1@data@values)){
            if (Mos1@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mos1@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mos1@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 3){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 10){
                Mos_Vals[i] = 3
            } else if (Mos1@data@values[i] == 11){
                Mos_Vals[i] = 3
            } else if (Mos1@data@values[i] == 12){
                Mos_Vals[i] = 3
            }
        }
        Mos1@data@values = as.numeric(Mos_Vals)
    } else if (j == 3){
        Mos1 = mosaic(C_base, F_base, fun = sum, na.rm = T)
        Mos1@data@values = as.integer(Mos1@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mos1@data@values)){
            if (Mos1@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mos1@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mos1@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 3){
                Mos_Vals[i] = 2
            } 
        }
        Mos1@data@values = as.numeric(Mos_Vals)
    }

    # Now estimating building extents 

    Building_Estimates = raster::predict(object = Base_Rasters[[j]], model = Expanded_nnet, type = 'raw')
    Building_Estimates = projectRaster(Building_Estimates, crs = '+proj=longlat +datum=WGS84', method = "ngb")
    Building_Estimates@data@values = as.integer(Building_Estimates@data@values)
    B_vals = rep("10", length(Building_Estimates@data@values))
    B_vals[c(which(Building_Estimates@data@values == 2))] = "0"
    B_vals[c(which(is.na(Building_Estimates@data@values)))] = "0"
    Building_Estimates@data@values = as.integer(B_vals)

    S_base = All_Shapes[[j]][[2]] %>% 
    mutate("Type" = "Road") %>%
    select(Type, geometry)
    S_base = rasterize(S_base, Sat_Images[[j]])
    S_base = projectRaster(S_base, crs = '+proj=longlat +datum=WGS84')
    S_vals = rep("20", length(S_base@data@values))
    S_vals[c(which(is.na(S_base@data@values)))] = "0"
    S_base@data@values = as.numeric(S_vals)

    Mosaic_Full = mosaic(Mos1, Building_Estimates, S_base, fun=sum, na.rm=T)
    Mosaic_Full@data@values = as.integer(Mosaic_Full@data@values)
    Mos_Vals = c()
    for (i in 1:length(Mosaic_Full@data@values)){
        if (Mosaic_Full@data@values[i] == 0){
            Mos_Vals[i] = 0
        } else if (Mosaic_Full@data@values[i] == 1){
            Mos_Vals[i] = 1
        } else if (Mosaic_Full@data@values[i] == 2){
            Mos_Vals[i] = 2
        } else if (Mosaic_Full@data@values[i] == 3){
            Mos_Vals[i] = 3
        } else if (Mosaic_Full@data@values[i] == 10){
            Mos_Vals[i] = 4
        } else if (Mosaic_Full@data@values[i] == 11){
            Mos_Vals[i] = 4
        } else if (Mosaic_Full@data@values[i] == 12){
            Mos_Vals[i] = 4
        } else if (Mosaic_Full@data@values[i] == 13){
            Mos_Vals[i] = 3 # Note that water will take precedence over building here.
        } else if (Mosaic_Full@data@values[i] == 20){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 21){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 22){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 23){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 30){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 31){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 32){
            Mos_Vals[i] = 5
        } else if (Mosaic_Full@data@values[i] == 33){
            Mos_Vals[i] = 5
        }
    }

    Mosaic_Full@data@values = as.numeric(Mos_Vals)

    writeRaster(x = Mosaic_Full,
    filename = paste0("/home/joshk/rasterdump/Mosaics/Revised/", names(Coord_Matrices[j]), "_Projected_Buildings.tif"),
    format   = "GTiff",
    datatype = dataType(Mosaic_Full),
    options  = "TFW = YES",
    overwrite = TRUE)
}

# Reading in and plotting 

Mosaics_Est = vector("list", 6); Mosaics_Est_Rasters = vector("list", 6); Mosaics_Est_Plots = vector("list", 6)

for (i in 1:6){
    Mosaic_Path = paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[i]), "_Projected_Buildings.tif")
    Mosaics_Est_Rasters[[i]] = raster(paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[i]), "_Projected_Buildings.tif"))
    Mosaics_Est[[i]] = read_stars(Mosaic_Path)    
    Mosaics_Est[[i]] = st_as_sf(Mosaics_Est[[i]], as_points = FALSE, merge = TRUE)
    colnames(Mosaics_Est[[i]])[1] = "Class"
    Mosaics_Est[[i]]$Land_Use = "Unclassified"
    Mosaics_Est[[i]]$Land_Use = as.character(Mosaics_Est[[i]]$Land_Use)

    {
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 0))] = "Unclassified"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 1))] = "Cropland"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 2))] = "Forest"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 3))] = "Water"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 4))] = "Building"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 5))] = "Street"
    }

    Mosaics_Est[[i]]$Land_Use = factor(Mosaics_Est[[i]]$Land_Use, levels = c("Unclassified","Building","Cropland","Forest","Street", "Water"))

    if (i %in% c(1,2,4,5,6)){
    Mosaics_Est_Plots[[i]] = ggplot(data.frame("Y" = c(names(Coord_Matrices[i])))) + 
        geom_sf(data =  Mosaics_Est[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("white","burlywood4","wheat","lightgreen","black","turquoise")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12))
    } else if (i == 3){
    Mosaics_Est_Plots[[i]] = ggplot(data.frame("Y" = c(names(Coord_Matrices[i])))) + 
        geom_sf(data =  Mosaics_List[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("white","wheat","lightgreen","black")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12))
    }     
}

do.call("grid.arrange", c(Mosaics_Est_Plots, ncol=3))

# Note that building projections are failing for Corwhin. Running this site individually to assess methodological failings.

Corwhin_Crop = All_Shapes[[3]][[5]] %>% 
    mutate("Type" = "Crop") %>%
    select(Type, geometry)
Corwhin_Crop = rasterize(Corwhin_Crop, Sat_Images[[3]])
Corwhin_Crop = projectRaster(Corwhin_Crop, crs = '+proj=longlat +datum=WGS84')
C_vals = rep("1", length(Corwhin_Crop@data@values))
C_vals[c(which(is.na(Corwhin_Crop@data@values)))] = "0"
Corwhin_Crop@data@values = C_vals
Corwhin_Crop@data@values = as.numeric(Corwhin_Crop@data@values)
mapView(Corwhin_Crop)

Corwhin_Forest = All_Shapes[[3]][[4]] %>% 
    mutate("Type" = "Forest") %>%
    select(Type, geometry)
Corwhin_Forest = rasterize(Corwhin_Forest, Sat_Images[[3]])
Corwhin_Forest = projectRaster(Corwhin_Forest, crs = '+proj=longlat +datum=WGS84')
F_vals = rep("2", length(Corwhin_Forest@data@values))
F_vals[c(which(is.na(Corwhin_Forest@data@values)))] = "0"
Corwhin_Forest@data@values = as.numeric(F_vals)
mapView(Corwhin_Forest)

Corwhin_Buildings = raster::predict(object = Base_Rasters[[3]], model = Expanded_nnet, type = 'raw')
Corwhin_Buildings = projectRaster(Corwhin_Buildings, crs = '+proj=longlat +datum=WGS84', method = "ngb")
Corwhin_Buildings@data@values = as.integer(Corwhin_Buildings@data@values)
B_vals = rep("10", length(Corwhin_Buildings@data@values))
B_vals[c(which(Corwhin_Buildings@data@values == 2))] = "0"
B_vals[c(which(is.na(Corwhin_Buildings@data@values)))] = "0"
Corwhin_Buildings@data@values = as.integer(B_vals)
mapView(Corwhin_Buildings)

Corwhin_Streets = All_Shapes[[3]][[2]] %>% 
    mutate("Type" = "Road") %>%
    select(Type, geometry)
Corwhin_Streets = rasterize(Corwhin_Streets, Sat_Images[[3]])
Corwhin_Streets = projectRaster(Corwhin_Streets, crs = '+proj=longlat +datum=WGS84')
S_vals = rep("20", length(Corwhin_Streets@data@values))
S_vals[c(which(is.na(Corwhin_Streets@data@values)))] = "0"
Corwhin_Streets@data@values = as.numeric(S_vals)
mapView(Corwhin_Streets)

# All individual layes appear okay. Merging.

Corwhin_Mosaic = mosaic(Corwhin_Crop, Corwhin_Forest, Corwhin_Buildings, Corwhin_Streets, fun=sum, na.rm=T)
Corwhin_Mosaic@data@values = as.integer(Corwhin_Mosaic@data@values)
    Mos_Vals = c()
    for (i in 1:length(Corwhin_Mosaic@data@values)){
        if (Corwhin_Mosaic@data@values[i] == 0){
            Mos_Vals[i] = 0
        } else if (Corwhin_Mosaic@data@values[i] == 1){
            Mos_Vals[i] = 1
        } else if (Corwhin_Mosaic@data@values[i] == 2){
            Mos_Vals[i] = 2
        } else if (Corwhin_Mosaic@data@values[i] == 3){
            Mos_Vals[i] = 2
        } else if (Corwhin_Mosaic@data@values[i] == 10){
            Mos_Vals[i] = 4
        } else if (Corwhin_Mosaic@data@values[i] == 11){
            Mos_Vals[i] = 4
        } else if (Corwhin_Mosaic@data@values[i] == 12){
            Mos_Vals[i] = 4
        } else if (Corwhin_Mosaic@data@values[i] == 13){
            Mos_Vals[i] = 4
        } else if (Corwhin_Mosaic@data@values[i] == 20){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 21){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 22){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 23){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 30){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 31){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 32){
            Mos_Vals[i] = 5
        } else if (Corwhin_Mosaic@data@values[i] == 33){
            Mos_Vals[i] = 5
        }
    }

Corwhin_Mosaic@data@values = as.numeric(Mos_Vals)
mapView(Corwhin_Mosaic)

writeRaster(x = Corwhin_Mosaic,
    filename = paste0("/home/joshk/rasterdump/Mosaics/Revised/Corwhin_Projected_Buildings.tif"),
    format   = "GTiff",
    datatype = dataType(Corwhin_Mosaic),
    options  = "TFW = YES",
    overwrite = TRUE)

Mosaics_Est = vector("list", 6); Mosaics_Est_Rasters = vector("list", 6); Mosaics_Est_Plots = vector("list", 6)
Plot_Labels = c("Erin", "Guelph", "Corwhin", "Brantford", "Cambridge", "Ruthven Park")

for (i in 1:6){
    Mosaic_Path = paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[i]), "_Projected_Buildings.tif")
    Mosaics_Est_Rasters[[i]] = raster(paste0("/home/joshk/rasterdump/Mosaics/", names(Coord_Matrices[i]), "_Projected_Buildings.tif"))
    Mosaics_Est[[i]] = read_stars(Mosaic_Path)    
    Mosaics_Est[[i]] = st_as_sf(Mosaics_Est[[i]], as_points = FALSE, merge = TRUE)
    colnames(Mosaics_Est[[i]])[1] = "Class"
    Mosaics_Est[[i]]$Land_Use = "Bare Ground"
    Mosaics_Est[[i]]$Land_Use = as.character(Mosaics_Est[[i]]$Land_Use)

    {
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 0))] = "Bare Ground"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 1))] = "Cropland"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 2))] = "Forest"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 3))] = "Water"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 4))] = "Building"
        Mosaics_Est[[i]]$Land_Use[c(which(Mosaics_Est[[i]]$Class == 5))] = "Street"
    }

    if (i %in% c(1:2,4:6)){
        Mosaics_Est[[i]]$Land_Use = factor(Mosaics_Est[[i]]$Land_Use, levels = c("Bare Ground","Building","Cropland","Forest","Street", "Water"))
    } else if (i == 3){
        Mosaics_Est[[i]]$Land_Use = factor(Mosaics_Est[[i]]$Land_Use, levels = c("Bare Ground","Building","Cropland","Forest","Street"))
    }

    if (i %in% c(1,2,4,5,6)){
    Mosaics_Est_Plots[[i]] = ggplot(data.frame("Y" = c(Plot_Labels[i]))) + 
        geom_sf(data =  Mosaics_Est[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class",
        values = c("grey98","burlywood4","wheat","lightgreen","black","turquoise")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12)) + 
        theme(legend.position = "none", panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"))
    } else if (i == 3){
    Mosaics_Est_Plots[[i]] = ggplot(data.frame("Y" = c(Plot_Labels[i]))) + 
        geom_sf(data =  Mosaics_Est[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("grey98","burlywood4","wheat","lightgreen","black")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12)) + 
        theme(legend.position = "none", panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"))
    }     
}

do.call("grid.arrange", c(Mosaics_Est_Plots, ncol=3))

# Pulling out common legend

extractLegend = function(xplot) {
  grobs = ggplot_gtable(ggplot_build(xplot))
  g_title = which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[g_title]]
}

Legend_Plot = ggplot(data.frame("Y" = c(names(Coord_Matrices[1])))) + 
        geom_sf(data =  Mosaics_Est[[1]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = NULL, values = c("grey98","burlywood4","wheat","lightgreen","black","turquoise")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12)) + 
        theme(panel.grid.major = element_blank(), legend.position = "bottom")

Common_Legend = extractLegend(Legend_Plot)             

grid.arrange(arrangeGrob(Mosaics_Est_Plots[[1]], Mosaics_Est_Plots[[2]], Mosaics_Est_Plots[[3]],
    Mosaics_Est_Plots[[4]], Mosaics_Est_Plots[[5]], Mosaics_Est_Plots[[6]], nrow = 2),
    Common_Legend, nrow = 2, heights = c(15, 1))

To_Export = arrangeGrob(arrangeGrob(Mosaics_Est_Plots[[1]], Mosaics_Est_Plots[[2]], Mosaics_Est_Plots[[3]],
    Mosaics_Est_Plots[[4]], Mosaics_Est_Plots[[5]], Mosaics_Est_Plots[[6]], nrow = 2),
    Common_Legend, nrow = 2, heights = c(15, 1))

ggsave("/home/joshk/rasterdump/Mosaics/Panel_Plot.pdf", To_Export, height = 7.5, width = 10, dpi = 2000)
ggsave("/home/joshk/rasterdump/Mosaics/Panel_Plot.jpg", To_Export, height = 7.5, width = 10, dpi = 800)

# Lastly, re-arranging plots and labelling according to their ecotype

Rural_Label = ggplot(data.frame("X" = c(0,1), "Y" = c(0,1)), x = X, y = Y) + 
    annotate("text", x = 0.5, y = 0.5, label = expression(underline("Rural")), size = 4, colour = "black", angle = 90) + theme_void()

Urban_Label = ggplot(data.frame("X" = c(0,1), "Y" = c(0,1)), x = X, y = Y) + 
    annotate("text", x = 0.5, y = 0.5, label = expression(underline("Urban")), size = 4, colour = "black", angle = 90) + theme_void()

P1 = Mosaics_Est_Plots[[1]] + theme(plot.margin = unit(c(0.5,-0.5,0.5,-0.5), "cm")); P2 = Mosaics_Est_Plots[[2]] + theme(plot.margin = unit(c(0.5,-0.5,0.5,-0.5), "cm")); P3 = Mosaics_Est_Plots[[3]] + theme(plot.margin = unit(c(0.5,-0.5,0.5,-0.5), "cm")); P4 = Mosaics_Est_Plots[[4]] + theme(plot.margin = unit(c(0.5,-0.5,0.5,-0.5), "cm")); P5 = Mosaics_Est_Plots[[5]] + theme(plot.margin = unit(c(0.5,-0.5,0.5,-0.5), "cm")); P6 = Mosaics_Est_Plots[[6]] + theme(plot.margin = unit(c(0.5,-0.5,0.5,-0.5), "cm"))

grid.arrange(arrangeGrob(arrangeGrob(Urban_Label, Rural_Label, nrow = 2), 
    arrangeGrob(P2, P4, P5, P1, P3, P6, nrow = 2),
    widths = c(1,25)), Common_Legend, nrow = 2, heights = c(15, 1))

Labelled_Export = arrangeGrob(arrangeGrob(arrangeGrob(Urban_Label, Rural_Label, nrow = 2), 
    arrangeGrob(P2, P4, P5, P1, P3, P6, nrow = 2),
    widths = c(1,25)), Common_Legend, nrow = 2, heights = c(15, 1))

ggsave("/home/joshk/rasterdump/Mosaics/Panel_Plot_Labelled.pdf", Labelled_Export, height = 7.5, width = 10, dpi = 2000)
ggsave("/home/joshk/rasterdump/Mosaics/Panel_Plot_Labelled.jpg", Labelled_Export, height = 7.5, width = 10, dpi = 800)

# Good. Continuing on to quantify the degree of urbanisation at each capture local.
# Note that to do so, we will need to randomly sample more locations (that is, to collect sufficient 
# data for PCA construction). Below, data collection is conducted for 2 additional urban and 2 additional 
# rural locations.

# Urban locations

Hamilton = c(-79.869415 + 0.025, -79.869415 - 0.025, 43.227346 + 0.025, 43.227346 - 0.025); Kitchener = c(-80.532002 + 0.025, -80.532002 - 0.025, 43.460306 + 0.025, 43.460306 - 0.025); Gowanstown = c(-80.855071 + 0.025, -80.855071 - 0.025, 43.762024 + 0.025, 43.762024 - 0.025); Luther_Marsh = c(-80.419695 + 0.025, -80.419695 - 0.025, 43.938366 + 0.025, 43.938366 - 0.025)

R_Coordinates = list(Hamilton, Kitchener, Gowanstown, Luther_Marsh)
names(R_Coordinates) = c("Hamilton", "Kitchener", "Gowanstown", "Luther_Marsh")
R_Coord_Matrices = vector("list", 4)

for (i in 1:length(R_Coord_Matrices)){
    R_Coord_Matrices[[i]] = matrix(data = c(R_Coordinates[[i]][1], R_Coordinates[[i]][3],
                       R_Coordinates[[i]][2], R_Coordinates[[i]][3],
                       R_Coordinates[[i]][2], R_Coordinates[[i]][4], 
                       R_Coordinates[[i]][1], R_Coordinates[[i]][4],
                       R_Coordinates[[i]][1], R_Coordinates[[i]][3]),
              ncol = 2, byrow = TRUE)
}

names(R_Coord_Matrices) = names(R_Coordinates)

for (i in 1:length(R_Coord_Matrices)){
    Dest = paste0("./rasterdump/Sentinel2/", names(R_Coord_Matrices[i]), "/")
    set_archive(Dest)

    set_aoi(R_Coord_Matrices[[i]])

    records = getSentinel_query(time_range = c("2018-04-01",
        as.character(Sys.Date())), platform = "Sentinel-2")

    records_filtered = records %>% 
    filter(level == "Level-2A") %>% 
    filter(cloudcov <= 1)

    datasets = getSentinel_data(records = records_filtered[1,])
    unzip(zipfile = list.files(paste0(Dest, "_datasets/Sentinel-2"), full.names = TRUE)[1], exdir = paste0(Dest, "_datasets/Sentinel-2"))
}

# Again, trimming out street, river, forest and cropland data. Note that building data is excluded here.

Rivers = vector("list", 4); Roads = vector("list", 4); Forests = vector("list", 4); Cropland = vector("list", 4) 

R_Source_List = list("/home/joshk/rasterdump/Rivers/ghy_000c11a_e.shp", "/home/joshk/rasterdump/Roads/lrnf000r16a_e.shp", "/home/joshk/rasterdump/Forests/Wooded_Area.shp", "/home/joshk/rasterdump/Land_Use/IMG_AAFC_LANDUSE_Z17_2010.tif")

names(R_Source_List) = c("Rivers", "Roads", "Forests", "Cropland")

for (j in 1:length(R_Source_List)){
    for (i in 1:length(R_Coord_Matrices)){
        if (i == 1 & j %in% c(1:3)){
            Shape_Dat = st_read(R_Source_List[[j]])
        } else if (i == 1 & j == 4){
            Shape_Dat = read_stars(R_Source_List[[j]])
        } 

        if (j %in% c(1:3)){
            Lin_Mat = c(R_Coord_Matrices[[i]][1,1], R_Coord_Matrices[[i]][2,1], R_Coord_Matrices[[i]][1,2], R_Coord_Matrices[[i]][3,2])
            xmin = Lin_Mat[2]; xmax = Lin_Mat[1]; ymin = Lin_Mat[4]; ymax = Lin_Mat[3]
    
            filter_func = sf::st_bbox(c(xmin = ifelse(is.na(xmin), -180, xmin), 
                ymin = ifelse(is.na(ymin),  -90,  ymin), 
                xmax = ifelse(is.na(xmax), +180, xmax), 
                ymax = ifelse(is.na(ymax),  +90, ymax)), 
                crs = st_crs(4269)) %>% 
            sf::st_as_sfc(.)
            if (j %in% c(2,3)){
                Shape_Dat = st_transform(Shape_Dat, 4269)
            }

            Trimmed_Data = sf::st_intersects(Shape_Dat, filter_func)
            To_Write = Shape_Dat[which(lengths(Trimmed_Data) != 0), ]

            st_write(To_Write, 
                paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(R_Coord_Matrices)[i], "_", names(R_Source_List)[j], ".shp"), driver = "ESRI Shapefile") 

            if (i == 4){
                rm(Shape_Dat)
            }
        } else if (j == 4){
                Lin_Mat = c(R_Coord_Matrices[[i]][1,1], R_Coord_Matrices[[i]][2,1], R_Coord_Matrices[[i]][1,2], R_Coord_Matrices[[i]][3,2])
                xmin = Lin_Mat[2]; xmax = Lin_Mat[1]; ymin = Lin_Mat[4]; ymax = Lin_Mat[3]
    
                filter_func = sf::st_bbox(c(xmin = ifelse(is.na(xmin), -180, xmin), 
                ymin = ifelse(is.na(ymin),  -90,  ymin), 
                xmax = ifelse(is.na(xmax), +180, xmax), 
                ymax = ifelse(is.na(ymax),  +90, ymax)), 
                crs = st_crs(4269)) %>% 
                sf::st_as_sfc(.)

                filter_func = st_transform(filter_func, 
                crs = '+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

                Trimmed_Data = st_crop(Shape_Dat, filter_func)
                Trimmed_Data = st_transform(Trimmed_Data, 
                    "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
                Shape_Data_Joined = st_as_sf(Trimmed_Data, as_points = FALSE, merge = FALSE)

                To_Write = Shape_Data_Joined %>% 
                    filter(IMG_AAFC_LANDUSE_Z17_2010.tif == 51) %>%
                    summarise(geometry = sf::st_union(geometry)) %>%
                    ungroup()

                st_write(To_Write, 
                paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(R_Coord_Matrices)[i], "_", names(R_Source_List)[j], ".shp"), driver = "ESRI Shapefile")         

                if (i == 4){
                    rm(Shape_Dat)
                }
        }
    }
}

# Flattening Sentinel-2 images for use in projections.

R_Base_Rasters = vector("list", 4)

for (j in 1:4){
    Base_fold = paste0("./rasterdump/Sentinel2/", names(R_Coord_Matrices[j]), "/_datasets/Sentinel-2/")
    jp2_folder = paste0(Base_fold, list.files(Base_fold, pattern = ".SAFE")[1], "/GRANULE/", list.files(paste0(Base_fold, "/", list.files(Base_fold, pattern = ".SAFE")[1], "/GRANULE"))[1], "/IMG_DATA/")

    jp2_10m = list.files(paste0(jp2_folder, "R10m"), pattern = ".*B.*jp2$", full.names = TRUE)
    jp2_20m = list.files(paste0(jp2_folder, "R20m"), pattern = ".*B.*[56781].*jp2$", full.names = TRUE)
    rst_lst = lapply(c(jp2_10m, jp2_20m), FUN = raster)

    extent_crop = R_Coord_Matrices[[j]] %>% 
        project(proj = proj4string(rst_lst[[1]])) %>% 
        apply(MARGIN = 2, FUN = range) %>%
        t %>%
        extent

    rst_crop_lst = lapply(rst_lst, FUN = raster::crop, y = extent_crop)
    names(rst_crop_lst) = str_extract(sapply(rst_crop_lst, names), "B.{2}")

    path_crop = paste0("/home/joshk/rasterdump/Sentinel2/", names(R_Coord_Matrices[j]), "/", names(rst_crop_lst), ".tif")

    for (i in 1:length(rst_crop_lst)){
    writeRaster(x = rst_crop_lst[[i]],
        filename = path_crop[i],
        format   = "GTiff",
        datatype = dataType(rst_crop_lst[[i]]),
        options  = "TFW=YES",
        overwrite = TRUE)
    }

    rst_for_prediction = vector(mode = "list", length = length(rst_crop_lst))
    names(rst_for_prediction) = names(rst_crop_lst)

    for (b in c("B05", "B06", "B07", "B8A", "B11", "B12")){
      rst_for_prediction[[b]] = raster::resample(x = rst_crop_lst[[b]],
                                                  y = rst_crop_lst$B02)
    }

    b_10m = c("B02", "B03", "B04", "B08")
    rst_for_prediction[b_10m] = rst_crop_lst[b_10m]
    brick_for_prediction = brick(rst_for_prediction)

    brick_for_prediction_norm = normImage(brick_for_prediction)
    names(brick_for_prediction_norm) = names(brick_for_prediction)

    R_Base_Rasters[[j]] = projectRaster(brick_for_prediction_norm, crs = '+proj=longlat +datum=WGS84')
    Base_Path = paste0("/home/joshk/rasterdump/Trimmed_Satellite/", names(R_Coord_Matrices[j]), "_Compiled.tif")

    Flat_Raster = merge(R_Base_Rasters[[j]], overlap = FALSE)

    writeRaster(x = Flat_Raster,
        filename = Base_Path,
        format   = "GTiff",
        datatype = dataType(Flat_Raster),
        options  = "TFW = YES",
        overwrite = TRUE)  
    }
}

# Projecting buildings and combining layers into mosaics per region.

Hamilton = vector("list", 5); Kitchener = vector("list", 5); Gowanstown = vector("list", 5); Luther_Marsh = vector("list", 5)

R_All_Shapes = list(Hamilton, Kitchener, Gowanstown, Luther_Marsh)
names(R_All_Shapes) = c("Hamilton", "Kitchener", "Gowanstown", "Luther_Marsh")

for (j in 1:4){
    for (i in 1:4){
        To_Read = paste0("/home/joshk/rasterdump/Trimmed_Shapefiles/", names(R_Coord_Matrices)[j], "_", names(R_Source_List)[i], ".shp")
        R_All_Shapes[[j]][[i]] = read_sf(To_Read)
    }
}

R_Sat_Images = vector('list', 4)

for (i in 1:4){
    Load_Sat = paste0("/home/joshk/rasterdump/Trimmed_Satellite/", names(R_Coord_Matrices[i]), "_Compiled.tif")
    R_Sat_Images[[i]] = raster(Load_Sat)
}

for (j in 1:4){

    C_base = R_All_Shapes[[j]][[4]] %>% 
    mutate("Type" = "Crop") %>%
    select(Type, geometry)
    C_base = rasterize(C_base, R_Sat_Images[[j]])
    C_base = projectRaster(C_base, crs = '+proj=longlat +datum=WGS84')
    C_vals = rep("1", length(C_base@data@values))
    C_vals[c(which(is.na(C_base@data@values)))] = "0"
    C_base@data@values = C_vals
    C_base@data@values = as.numeric(C_base@data@values)

    F_base = R_All_Shapes[[j]][[3]] %>% 
    mutate("Type" = "Forest") %>%
    select(Type, geometry)
    F_base = rasterize(F_base, R_Sat_Images[[j]])
    F_base = projectRaster(F_base, crs = '+proj=longlat +datum=WGS84')
    F_vals = rep("2", length(F_base@data@values))
    F_vals[c(which(is.na(F_base@data@values)))] = "0"
    F_base@data@values = as.numeric(F_vals)

    if (j %in% c(2,4)){
        R_base = R_All_Shapes[[j]][[1]] %>% 
        mutate("Type" = "Water") %>%
        select(Type, geometry)
        R_base = rasterize(R_base, R_Sat_Images[[j]])
        R_base = projectRaster(R_base, crs = '+proj=longlat +datum=WGS84')
        R_vals = rep("10", length(R_base@data@values))
        R_vals[c(which(is.na(R_base@data@values)))] = "0"
        R_base@data@values = as.numeric(R_vals)

        Mos1 = mosaic(C_base, F_base, R_base, fun = sum, na.rm = T)
        Mos1@data@values = as.integer(Mos1@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mos1@data@values)){
            if (Mos1@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mos1@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mos1@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 3){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 10){
                Mos_Vals[i] = 3
            } else if (Mos1@data@values[i] == 11){
                Mos_Vals[i] = 3
            } else if (Mos1@data@values[i] == 12){
                Mos_Vals[i] = 3
            }
        }
        Mos1@data@values = as.numeric(Mos_Vals)
    } else if (j %in% c(1,3)){
        Mos1 = mosaic(C_base, F_base, fun = sum, na.rm = T)
        Mos1@data@values = as.integer(Mos1@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mos1@data@values)){
            if (Mos1@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mos1@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mos1@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mos1@data@values[i] == 3){
                Mos_Vals[i] = 2
            } 
        }
        Mos1@data@values = as.numeric(Mos_Vals)
    }

    Building_Estimates = raster::predict(object = R_Base_Rasters[[j]], model = Expanded_nnet, type = 'raw')
    Building_Estimates = projectRaster(Building_Estimates, crs = '+proj=longlat +datum=WGS84', method = "ngb")
    Building_Estimates@data@values = as.integer(Building_Estimates@data@values)
    B_vals = rep("10", length(Building_Estimates@data@values))
    B_vals[c(which(Building_Estimates@data@values == 2))] = "0"
    B_vals[c(which(is.na(Building_Estimates@data@values)))] = "0"
    Building_Estimates@data@values = as.integer(B_vals)

    S_base = R_All_Shapes[[j]][[2]] %>% 
    mutate("Type" = "Road") %>%
    select(Type, geometry)
    S_base = rasterize(S_base, R_Sat_Images[[j]])
    S_base = projectRaster(S_base, crs = '+proj=longlat +datum=WGS84')
    S_vals = rep("20", length(S_base@data@values))
    S_vals[c(which(is.na(S_base@data@values)))] = "0"
    S_base@data@values = as.numeric(S_vals)

    if (j %in% c(2,4)){
        Mosaic_Full = mosaic(Mos1, Building_Estimates, S_base, fun=sum, na.rm=T)
        Mosaic_Full@data@values = as.integer(Mosaic_Full@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mosaic_Full@data@values)){
            if (Mosaic_Full@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mosaic_Full@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mosaic_Full@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mosaic_Full@data@values[i] == 3){
                Mos_Vals[i] = 3
            } else if (Mosaic_Full@data@values[i] == 10){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 11){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 12){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 13){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 20){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 21){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 22){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 23){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 30){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 31){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 32){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 33){
                Mos_Vals[i] = 5
            }
        }

        Mosaic_Full@data@values = as.numeric(Mos_Vals)

        writeRaster(x = Mosaic_Full,
        filename = paste0("/home/joshk/rasterdump/Mosaics/", names(R_Coord_Matrices[j]), "_Projected_Buildings.tif"),
        format   = "GTiff",
        datatype = dataType(Mosaic_Full),
        options  = "TFW = YES",
        overwrite = TRUE)
    } else if (j %in% c(1,3)){
        Mosaic_Full = mosaic(Mos1, Building_Estimates, S_base, fun=sum, na.rm=T)
        Mosaic_Full@data@values = as.integer(Mosaic_Full@data@values)
        Mos_Vals = c()
        for (i in 1:length(Mosaic_Full@data@values)){
            if (Mosaic_Full@data@values[i] == 0){
                Mos_Vals[i] = 0
            } else if (Mosaic_Full@data@values[i] == 1){
                Mos_Vals[i] = 1
            } else if (Mosaic_Full@data@values[i] == 2){
                Mos_Vals[i] = 2
            } else if (Mosaic_Full@data@values[i] == 3){
                Mos_Vals[i] = 2
            } else if (Mosaic_Full@data@values[i] == 10){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 11){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 12){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 13){
                Mos_Vals[i] = 4
            } else if (Mosaic_Full@data@values[i] == 20){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 21){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 22){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 23){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 30){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 31){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 32){
                Mos_Vals[i] = 5
            } else if (Mosaic_Full@data@values[i] == 33){
                Mos_Vals[i] = 5
            }
        }

        Mosaic_Full@data@values = as.numeric(Mos_Vals)

        writeRaster(x = Mosaic_Full,
        filename = paste0("/home/joshk/rasterdump/Mosaics/", names(R_Coord_Matrices[j]), "_Projected_Buildings.tif"),
        format   = "GTiff",
        datatype = dataType(Mosaic_Full),
        options  = "TFW = YES",
        overwrite = TRUE)
    }
}

# Test plotting

R_Mosaics = vector("list", 4); R_Mosaics_Rasters = vector("list", 4); R_Mosaics_Plots = vector("list", 4)
Plot_Labels = c("Hamilton", "Kitchener", "Gowanstown", "Luther Marsh")

for (i in 1:4){
    Mosaic_Path = paste0("/home/joshk/rasterdump/Mosaics/", names(R_Coord_Matrices[i]), "_Projected_Buildings.tif")
    R_Mosaics_Rasters[[i]] = raster(paste0("/home/joshk/rasterdump/Mosaics/", names(R_Coord_Matrices[i]), "_Projected_Buildings.tif"))
    R_Mosaics[[i]] = read_stars(Mosaic_Path)    
    R_Mosaics[[i]] = st_as_sf(R_Mosaics[[i]], as_points = FALSE, merge = TRUE)
    colnames(R_Mosaics[[i]])[1] = "Class"
    R_Mosaics[[i]]$Land_Use = "Bare Ground"
    R_Mosaics[[i]]$Land_Use = as.character(R_Mosaics[[i]]$Land_Use)

    {
        R_Mosaics[[i]]$Land_Use[c(which(R_Mosaics[[i]]$Class == 0))] = "Bare Ground"
        R_Mosaics[[i]]$Land_Use[c(which(R_Mosaics[[i]]$Class == 1))] = "Cropland"
        R_Mosaics[[i]]$Land_Use[c(which(R_Mosaics[[i]]$Class == 2))] = "Forest"
        R_Mosaics[[i]]$Land_Use[c(which(R_Mosaics[[i]]$Class == 3))] = "Water"
        R_Mosaics[[i]]$Land_Use[c(which(R_Mosaics[[i]]$Class == 4))] = "Building"
        R_Mosaics[[i]]$Land_Use[c(which(R_Mosaics[[i]]$Class == 5))] = "Street"
    }

    if (i %in% c(2,4)){
        R_Mosaics[[i]]$Land_Use = factor(R_Mosaics[[i]]$Land_Use, levels = c("Bare Ground","Building","Cropland","Forest","Street", "Water"))
    } else if (i %in% c(1,3)){
        R_Mosaics[[i]]$Land_Use = factor(R_Mosaics[[i]]$Land_Use, levels = c("Bare Ground","Building","Cropland","Forest","Street"))
    }

    if (i %in% c(2,4)){
    R_Mosaics_Plots[[i]] = ggplot(data.frame("Y" = c(Plot_Labels[i]))) + 
        geom_sf(data =  R_Mosaics[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class",
        values = c("grey98","burlywood4","wheat","lightgreen","black","turquoise")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12)) + 
        theme(legend.position = "none", panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"))
    } else if (i %in% c(1,3)){
    R_Mosaics_Plots[[i]] = ggplot(data.frame("Y" = c(Plot_Labels[i]))) + 
        geom_sf(data =  R_Mosaics[[i]], aes(fill = Land_Use), colour = "black", size = 0.1) + 
        scale_fill_manual(name = "Land Class", values = c("grey98","burlywood4","wheat","lightgreen","black")) +
        theme_minimal() + 
        facet_wrap(~Y) + 
        theme(strip.text.x = element_text(size = 12)) + 
        theme(legend.position = "none", panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"))
    }     
}

do.call("grid.arrange", c(R_Mosaics_Plots, ncol=2))

# Acceptable, but may be worthwhile to train model on building data including hamilton and kitchener.
# First, however, extracting pointwise data and running PCA.

PCA_List = vector("list", 10)
for (i in 1:10){
    if (i %in% c(1:6)){
        Full_Set = as.data.frame(rasterToPoints(Mosaics_Est_Rasters[[i]]))
        colnames(Full_Set)[3] = "Land_Use"
        PCA_List[[i]] = Full_Set %>% group_by(Land_Use) %>% 
        summarise("Count" = n()) %>% mutate("Site" = names(Coord_Matrices)[i])
    } else if (i %in% c(7:10)){
        Full_Set = as.data.frame(rasterToPoints(R_Mosaics_Rasters[[i-6]]))
        colnames(Full_Set)[3] = "Land_Use"
        PCA_List[[i]] = Full_Set %>% group_by(Land_Use) %>% 
        summarise("Count" = n()) %>% mutate("Site" = names(R_Coord_Matrices)[i-6])
    }
}

PCA_Dat = bind_rows(PCA_List) %>% pivot_wider(names_from = "Land_Use", values_from = "Count") %>%
    rename("Unclassifed" = "0", "Cropland" = "1", "Forest" = "2", "Water" = "3", "Building" = "4", "Street" = "5") %>%
    replace_na(list(Water = 0)) %>% as.data.frame(.)

PCA_Out = prcomp(PCA_Dat[,-1], center = TRUE, scale = TRUE)    
plot(PCA_Out)
summary(PCA_Out)

# Low percent of variance explained (~59%) as expected with such low quantities of data. Producing biplot

require(ggbiplot)

biplot_urban = ggbiplot(PCA_Out, size = 1) + theme_bw() + theme(panel.grid.major = element_blank()) + 
    xlim(c(-2,2)) + ylim(c(-2, 1.5)) + 
    xlab("Standardised PC1") + ylab("Standardised PC2") 

ggsave("/home/joshk/Desktop/Urban_Biplot.jpg", biplot_urban, dpi = 800)
ggsave("/home/joshk/Desktop/Urban_Biplot.pdf", biplot_urban, dpi = 800)

# PC1 axis alligns with degree of urbanisation (streets, buildings and unclassified). Extracting values and comparing means for our sample.

require(ggsignif)

PCA_Dat$Urbanisation = PCA_Out$x[,1]
Our_Sample = subset(PCA_Dat, !(Site %in% c("Hamilton", "Kitchener", "Gowanstown","Luther_Marsh")))
Our_Sample$N_Urbanisation = Our_Sample$Urbanisation + 10
Our_Sample$N_Urbanisation = ((Our_Sample$N_Urbanisation - min(Our_Sample$N_Urbanisation))/(max(Our_Sample$N_Urbanisation) - min(Our_Sample$N_Urbanisation)))*100
Our_Sample$Ecotype = ifelse(Our_Sample$Site %in% c("Guelph", "Cambridge", "Brantford"), "Urban", "Rural") 

# Calculating Bayes factor for comparison

require("brmsMethods")

Wide_Sample = data.frame("Urban" = subset(Our_Sample, Ecotype == "Urban")$N_Urbanisation,
    "Rural" = subset(Our_Sample, Ecotype == "Rural")$N_Urbanisation)

DF_Wide = build_hdf(vars = list(Wide_Sample$Urban, Wide_Sample$Rural),
  priors = list(rnorm(nrow(Wide_Sample), 50, 2), rnorm(nrow(Wide_Sample), 50, 2)),
  names = c("Urban", "Rural"))

hypothesis_test = hypothesis_df("Urban > Rural", DF_Wide, class = "b", alpha = 0.05)

plot(hypothesis_test)
print(hypothesis_test)

Fig_Save = plot(hypothesis_test, plot = F, theme = theme_get())[[1]]

Fig_Save = Fig_Save + theme_bw() + ylab("Density") + xlab("Difference in Degree of Urbanisation\n(Urban - Rural)") + 
theme(strip.text.x = element_blank())

ggsave("/home/joshk/Desktop/Urban_Hypothesis_Test.jpg", Fig_Save, dpi = 800)
ggsave("/home/joshk/Desktop/Urban_Hypothesis_Test.pdf", Fig_Save, dpi = 800)

# Simply too small of a data-set to know, but support for hypothesis.

Urban_Plot = ggplot(Our_Sample, aes(x = Ecotype, y = N_Urbanisation, fill = Ecotype)) + 
    stat_summary(geom = "errorbar", fun.data = "mean_se", size = 1, width = 0.4, colour = "black") + 
    stat_summary(geom = "point", fun = "mean", size = 5, pch = 21, colour = "black") + 
    theme_bw() + theme(legend.position = "none") + xlab("Capture Ecotype") + 
    ylab("Degree of Urbanisation (PC1: 0-100)") + 
    scale_fill_manual(values = c(viridis::viridis(n = 10)[4], "black")) + 
    ylim(c(0,100)) +
    geom_signif(comparisons = list(c("Urban", "Rural")), y_position = 95, annotation=c("*"),
        textsize = 10) 

ggsave("/home/joshk/rasterdump/Ecotype_Comparison_Plot.jpg", Urban_Plot, width = 8, height = 8, dpi = 800)
ggsave("/home/joshk/rasterdump/Ecotype_Comparison_Plot.pdf", Urban_Plot, width = 8, height = 8, dpi = 800)

plot_test = ggplot(data = data.frame("X" = c(0,1), "Y" = c(0,1)), aes(x = X, y = Y)) + 
    annotate("text", label = "Stress Response", size = 6, x = 0.45, y = 0.8) + 
    annotate("text", label = "Homeostasis", size = 6, x = 0.2, y = 0.7) + 
    annotate("text", label = "Direct Coping Response", size = 6, x = 0.8, y = 0.7) + 
    annotate("text", label = "Cellular Survival", size = 6, x = 0.2, y = 0.6) + 
    annotate("text", label = "Organismal Survival", size = 6, x = 0.8, y = 0.6) + 
    xlim(c(0,1)) + 
    geom_segment(x = 0.43, y = 0.77, xend = 0.2, yend = 0.73, size = 1, colour = "black", arrow = arrow(length = unit(0.03, "npc"))) + 
    geom_segment(x = 0.47, y = 0.77, xend = 0.8, yend = 0.73, size = 1, colour = "black", arrow = arrow(length = unit(0.03, "npc"))) + 
    geom_segment(x = 0.2, y = 0.67, xend = 0.2, yend = 0.63, size = 1, colour = "black", arrow = arrow(length = unit(0.03, "npc"))) +         
    geom_segment(x = 0.8, y = 0.67, xend = 0.8, yend = 0.63, size = 1, colour = "black", arrow = arrow(length = unit(0.03, "npc"))) + 
    geom_segment(x = 0.3, y = 0.7, xend = 0.63, yend = 0.7, size = 1, colour = "firebrick", arrow = arrow(length = unit(0.03, "npc"))) + 
    geom_segment(x = 0.35, y = 0.6, xend = 0.65, yend = 0.6, size = 1, colour = "black", arrow = arrow(length = unit(0.03, "npc"))) + 
    geom_segment(x = 0.63, y = 0.7, xend = 0.3, yend = 0.7, size = 1, colour = "firebrick", arrow = arrow(length = unit(0.03, "npc"))) + 
    geom_segment(x = 0.65, y = 0.6, xend = 0.35, yend = 0.6, size = 1, colour = "black", arrow = arrow(length = unit(0.03, "npc")))

plot_test 

ggsave("/home/joshk/Pictures/GE_Diagram.jpg", plot_test, width = 8.75, height = 2.5, dpi = 800)