library(sdm)
library(usdm)
library(dismo)
library(rJava)
library(raster)
library(parallel)
library(doParallel)
library(snowfall)
library(tidyverse)


baseline_simulated <- function(d_rsl){
  
  ev7 = read_csv(paste0(d_rsl,"Alternaria blotch_model_evaluation_tss07.csv"))
  sdm_methods = c('brt','mars','maxent','rf','svm')
  sdmnum <- 1:length(sdm_methods)
  m <- readRDS(paste0(d_rsl,"Alternaria blotch_trained_sdms.rds"))
  
  ### variables names
  variable_names_sorts = c("Bio1","Bio2", "Bio3" ,"Bio4","Bio5","Bio6","Bio7", "Bio8", 
                           "Bio9","Bio10","Bio11","Bio12", "Bio13" ,"Bio14", "Bio15",
                           "Bio16","Bio17","Bio18","Bio19", 
                           "Aspect","Slope","Curvature", "Elevation")
  
  
  out_folder_path = file.path(d_rsl,"1970-2000")
  # check dir presence?
  if (!dir.exists(out_folder_path)) {
    # do not presence,create it
    dir.create(out_folder_path)
  }
  
  en_mth=mean(ev7$threshold)   ##threshold
  # get input envs
  env_bio = stack("./data/bioclimate_variables/1970-2000/wc2.1_2.5m_bio1970-2000.tif")
  names(env_bio) = variable_names_sorts
  
  
  #### ensemble model
  enf = ensemble(m,env_bio, paste0(out_folder_path,"/","1970-2020",'_ensemble.tif'),
                 setting=list(method='weighted',id=ev7$modelID,stat="TSS", opt=2),overwrite = T)
  
  #----------Binary map and map change
  pa= raster(enf)
  pa[] = ifelse(enf[]>=en_mth,1,0)
  writeRaster(pa,paste0(out_folder_path,"/","1970-2000",'_ensemble_Binary.tif'),overwrite=TRUE)
  
  
  baseline_diff_sdms_projecting = function(i){
    smethod <- sdm_methods[i]
    sdmev7 <- subset(ev7, methods == smethod)
    sdmmth <- mean(sdmev7$threshold)
    pre <- predict(m, newdata = env_bio, method = smethod, mean = T,
                   filename = paste0(out_folder_path, "/", "1970-2000", '_', smethod, '_difsdms.tif'),
                   overwrite = T,
                   parallelSetting = list(ncore = 5)) 
    # bin map
    paf <- raster(pre)
    paf[] <- ifelse(pre[] >= sdmmth, 1, 0)
    writeRaster(paf, paste0(out_folder_path, "/", "1970-2000", '_', smethod, '_difsdms_Binary.tif'), overwrite = TRUE)
  }
  
  
  # ## Init snowfall
  sfInit(parallel=TRUE, cpus=5)  ## select 5CPUs
  
  ## Export packages
  sfLibrary('sdm', character.only=TRUE)
  sfLibrary('raster', character.only=TRUE)
  sfLibrary('usdm', character.only=TRUE)
  sfLibrary('rJava', character.only=TRUE)
  sfLibrary('parallel', character.only=TRUE)
  sfLibrary('snowfall', character.only=TRUE)
  
  ## ## Export variables
  sfExportAll()
  mySFModelsOut = sfLapply(sdmnum, baseline_diff_sdms_projecting)
  ## ## stop snowfall
  sfStop(nostop=FALSE)
  
}





future_predicet <- function(d_rsl){
  
  d_rsl <- "./results/base/"
  ev7 = read_csv(paste0(d_rsl,"Alternaria blotch_model_evaluation_tss07.csv"))
  m <- readRDS(paste0(d_rsl,"Alternaria blotch_trained_sdms.rds"))
  pa <- raster(paste0(d_rsl,"1970-2000/1970-2000_ensemble_Binary.tif"))
  
  ssps <- c("126", "245", "370","585")
  times_jurs <- c("2021-2040","2041-2060", "2061-2080","2081-2100")
  
  # Variables names
  variable_names_sorts <- c("Bio1","Bio2", "Bio3" ,"Bio4","Bio5","Bio6","Bio7", "Bio8", 
                            "Bio9","Bio10","Bio11","Bio12", "Bio13" ,"Bio14", "Bio15",
                            "Bio16","Bio17","Bio18","Bio19", 
                            "Aspect","Slope","Curvature", "Elevation")
  
  
  future_periods_ensemble_predicting <- function(params){
    times_jur <- as.character(params$times_jur)
    ssp <- as.character(params$ssp)
    gcm <- "ensemble_gcms"
    
    out_folder_times <- file.path(d_rsl, times_jur)
    if (!dir.exists(out_folder_times)) {
      dir.create(out_folder_times)
    }
    out_folder_path <- file.path(out_folder_times, ssp)
    if (!dir.exists(out_folder_path)) {
      dir.create(out_folder_path)
    }
    
    folder_path_time_ssp <- file.path("./data/bioclimate_variables/", times_jur, ssp)
    files_bio <- list.files(folder_path_time_ssp, pattern = "wc2.1_2.5m_bio.tif$", full.names = TRUE)
    
    env_bio <- stack(files_bio)
    names(env_bio) <- variable_names_sorts
    
    en_mth <- mean(ev7$threshold)
    enf <- ensemble(m, env_bio, paste0(out_folder_path, "/", times_jur, ssp, 
                                       '_', gcm, '_ensemble.tif'),
                    setting = list(method = 'weighted', id = ev7$modelID, 
                                   stat = "TSS", opt = 2), overwrite = TRUE)
    
    paf <- raster(enf)
    paf[] <- ifelse(enf[] >= en_mth, 5, 3)
    writeRaster(paf, paste0(out_folder_path, "/", times_jur, ssp, '_', gcm, 
                            '_ensemble_Binary.tif'), overwrite = TRUE)
    pac <- paf - pa
    writeRaster(pac, paste0(out_folder_path, "/", times_jur, ssp, '_', gcm, 
                            '_ensemble_Change.tif'), overwrite = TRUE)
  }
  
  
  future_periods_predicting <- function(params) {
    
    
    times_jur <- as.character(params$times_jur)
    ssp <- as.character(params$ssp)
    gcm <- as.character(params$gcms)
    smethod <- as.character(params$sdm_methods)
    
    out_folder_times <- file.path(d_rsl, times_jur)
    if (!dir.exists(out_folder_times)) {
      dir.create(out_folder_times)
    }
    out_folder_path <- file.path(out_folder_times, ssp)
    if (!dir.exists(out_folder_path)) {
      dir.create(out_folder_path)
    }
    
    folder_path_time_ssp <- file.path("./data/bioclimate_variables/", times_jur, ssp)
    files_bio <- list.files(folder_path_time_ssp, pattern = paste0(gcm,".tif$"), full.names = TRUE)
    env_bio <- stack(files_bio)
    names(env_bio) <- variable_names_sorts
    
    sdmev7 <- subset(ev7, methods == smethod)
    sdmmth <- mean(sdmev7$threshold)
    
    pre <- predict(m, newdata = env_bio, method = smethod, mean = T,
                   filename = paste0(out_folder_path, "/", times_jur, ssp, 
                                     '_', gcm, '_', smethod, '_ssp_gcm_sdm.tif'),
                   overwrite = TRUE)
    
    paf <- raster(pre)
    paf[] <- ifelse(pre[] >= sdmmth, 5, 3)
    writeRaster(paf, paste0(out_folder_path, "/", times_jur, ssp, '_', gcm, 
                            '_', smethod, '_ssp_gcm_sdm_Binary.tif'), overwrite = TRUE)
  }
  
  
  parallel_future_ensemble_predicting(){
    
    ## Continue with existing snowfall code to manage times_jur and ssp combinations
    ## Initialize snowfall
    sfInit(parallel = TRUE, cpus = 4) # set 16 cpus
    sfLibrary("sdm", character.only = TRUE)
    sfLibrary("raster", character.only = TRUE)
    sfLibrary("usdm", character.only = TRUE)
    sfLibrary("rJava", character.only = TRUE)
    sfLibrary("parallel", character.only = TRUE)
    sfLibrary("snowfall", character.only = TRUE)
    # sfExportAll()
    
    # Combine time periods and SSPs into a list of parameters
    params_list <- expand.grid(times_jurs = times_jurs, ssps = ssps, KEEP.OUT.ATTRS = FALSE)
    
    ## Export all necessary variables and functions to the slave nodes
    sfExport("params_list", "future_periods_ensemble_predicting", "d_rsl","ev7",
             "m", "variable_names_sorts", "pa")
    
    ## Apply the function across all combinations using snowfall
    execution_time <- system.time({
      mySFModelsOut <- sfLapply(1:nrow(params_list), 
                                function(x) future_periods_ensemble_predicting(params_list[x, ]))
    })
    ## Stop snowfall
    sfStop(nostop = FALSE)
    print(execution_time/3600)
  }
  
  parallel_future_Predicting(){
    
    sdm_methods <- c('brt','mars','maxent','rf','svm')
    
    # do not save ensemble future climate 
    gcms <- c("BCC-CSM2-MR","CanESM5","CNRM-CM6-1","CNRM-ESM2-1","MIROC6")
    
    ## Continue with existing snowfall code to manage times_jur and ssp combinations
    ## Initialize snowfall
    sfInit(parallel = TRUE, cpus = 8) # set 8 cpus
    sfLibrary("sdm", character.only = TRUE)
    sfLibrary("raster", character.only = TRUE)
    sfLibrary("usdm", character.only = TRUE)
    sfLibrary("rJava", character.only = TRUE)
    sfLibrary("parallel", character.only = TRUE)
    sfLibrary("snowfall", character.only = TRUE)
    # sfExportAll()
    
    # Combine time periods and SSPs into a list of parameters
    params_list <- expand.grid(times_jurs = times_jurs, ssps = ssps,gcms = gcms,
                               sdm_methods = sdm_methods, KEEP.OUT.ATTRS = FALSE)
    
    ## Export all necessary variables and functions to the slave nodes
    sfExport("params_list", "future_periods_predicting", "d_rsl", "ev7", "m", "variable_names_sorts", "pa")
    
    ## Apply the function across all combinations using snowfall
    execution_time <- system.time({
      mySFModelsOut <- sfLapply(1:nrow(params_list), 
                                function(x) future_periods_predicting(params_list[x, ]))
    })
    ## Stop snowfall
    sfStop(nostop = FALSE)
    print(execution_time/3600)
    
  }
  
  
  parallel_future_ensemble_predicting()
  parallel_future_Predicting()
  
}


simulated_and_predicted <- function(d_rsl){
  baseline_simulated()
  rm(list=ls())
  future_predicet()
}


simulated_and_predicted("./results/base_cor_0.8/")



