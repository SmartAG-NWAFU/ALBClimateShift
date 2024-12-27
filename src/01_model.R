library(sdm)
library(usdm)
library(dismo)
library(rJava)
library(raster)
library(parallel)
library(doParallel)
library(snowfall)
library(tidyverse)

rm(list=ls())


created_presence_pseudo_absence_points <- function(){
  library(sf)
  # filename <- paste0("./data/species/",disease_name,"_thin1.csv")
  # species_data = read.csv(filename)
  set.seed(1024)
  species_data <- read.csv("./data/species/Alternaria blotch_thin1.csv")
  species_data$records <- 1
  npnt <- nrow(species_data)
  nback <- ifelse(npnt < 1000,npnt,1000)
  
  country <- read_sf("./data/bounds/2023年县级/县级.shp")
  
  apple_disease_area = read.csv("./data/species/statistics_diseases_mean_area_2011_2019.csv")
  
  apple_disease_area_df <- left_join(country, apple_disease_area[c("JGDWMC","Alternaria.blotch")], 
                                     by = c("地名" = "JGDWMC"))
  
  temp_df = subset(apple_disease_area_df, is.na(apple_disease_area_df$Alternaria.blotch))
  
  env_bio <- stack("./data/bioclimate_variables/1970-2000/wc2.1_2.5m_bio1970-2000.tif")
  
  temp = mask(crop(env_bio[[1]], temp_df), temp_df)
  background_points = sdm::background(temp, n = nback, method = 'gRandom')
  
  data_temp <- as.data.frame(background_points[, c(1, 2)])
  names(data_temp) <- c("lon", "lat")
  
  # Add record columns
  data_temp$records <- 0
  
  sp_data <- bind_rows(species_data[c("lon","lat","records")], data_temp)
  write.csv(sp_data, "./data/species/presence_pseudo_absence.csv")
  
}





species_distribution_model <- function(d_rsl, disease_name, species_data,sdm_methods, biom){

  sp = species_data[, c("lon","lat","records")] 
  spdata=sp
  coordinates(sp)= ~lon + lat
  
  sdm_data = sdmData(formula=records~., train=sp, predictors= biom,
              parallelSetting = list(ncore = 5, method = 'parallel'))
  
  m = sdm(records~., sdm_data, methods = sdm_methods, replication='sub',test.percent=30, n=10,
          parallelSetting = list(ncore = 5, method = 'parallel', fork = F))  
  
  ## save models
  saveRDS(m,paste0(d_rsl,disease_name,"_trained_sdms.rds"))
  
  # m <- readRDS(paste0(d_rsl,"Alternaria blotch_trained_sdms.rds"))
  
  ####Model Performance
  ev=getEvaluation(m,stat=c('AUC','TSS','COR','threshold','Kappa',opt=2))  #opt =2 means "max(se+sp)" =maxSSS ?getEvaluation
  ev$methods = rep(sdm_methods, each = 10)
  write.csv(ev, file =paste0(d_rsl,"/",disease_name,'_model_evaluation.csv'), row.names = F)
  ev7=subset(ev,TSS>=0.7)   #####Choose the best model
  write.csv(ev7,file =paste0(d_rsl,"/",disease_name,'_model_evaluation_tss07.csv'), row.names = F)
  
  
  # ####Variable Importance
  # imp = getVarImp(m, id="ensemble", setting = list(method = 'weighted', stat='TSS', id = ev7$modelID))
  # write.csv(imp@varImportance,file =paste0(d_rsl,"/",disease_name,'_variable_importance.csv'), row.names = F )
  
  # ####Variable Importance
  imp = getVarImp(m)
  write.csv(imp@varImportanceMean,file =paste0(d_rsl,"/",disease_name,'_variable_importance.csv'), row.names = F )
  
  
  
  ###########Obtain response curve value
  dt = as.data.frame(m@data, grp=c('train'))
  si = nrow(spdata)
  rc = new('.responseCurve')
  nf = m@setting@featureFrame@predictors
  dm = data.frame(matrix(nrow=si,ncol=length(nf)))
  colnames(dm) = nf
  for (n in nf) dm[,n] = mean(dt[,n],na.rm=TRUE)
  
  for (n in nf) {
    dv = dm
    if (n == "Bio9"){
      dv[,n] = seq(min(dt[,n],na.rm=TRUE), max(dt[,n],na.rm=TRUE)+10,length.out = si)
    }
    else{
      dv[,n] = seq(min(dt[,n],na.rm=TRUE), max(dt[,n],na.rm=TRUE),length.out = si)
    }
    p = predict(m, newdata=dv)
    rc@response[[n]] = data.frame(dv[,n],rowMeans(p[,1:ncol(p)]))
    colnames(rc@response[[n]]) = c(n,'Response')
  }
  
  res_result = do.call(cbind,rc@response)
  write.csv(res_result,file =paste0(d_rsl,"/", disease_name,'_ReponseCurve_Value.csv'), row.names = F)
  
}



fitted_sdm_model <- function(disease_name){
  
  # disease_name <- "Alternaria blotch"

  species_data = read.csv("./data/species/presence_pseudo_absence.csv")
  
  sdm_methods = c('brt','mars','maxent','rf','svm')
  
  ### variables names
  variable_names_sorts = c("Bio1","Bio2", "Bio3" ,"Bio4","Bio5","Bio6","Bio7", "Bio8", 
                           "Bio9","Bio10","Bio11","Bio12", "Bio13" ,"Bio14", "Bio15",
                           "Bio16","Bio17","Bio18","Bio19", 
                           "Aspect","Slope","Curvature", "Elevation")
  # baseline
  bio = stack("./data/bioclimate_variables/1970-2000/wc2.1_2.5m_bio1970-2000.tif")
  names(bio) = variable_names_sorts
  
  
  ## save correlation matrix
  as.data.frame(bio,xy=F,na.rm = TRUE) %>% 
    cor() %>% 
    write.csv("./data/bio_correlation_matrix.csv")
  
  #-------Filter variables
  
  ## cor 0.8
  # set model output dir
  d_rsl <- "./results/base_cor_0.8/"
  v = vifcor(as.data.frame(bio), th = 0.8)
  biom = exclude(bio,v)
  

  ## save correlation matrix
  as.data.frame(biom,xy=F,na.rm = TRUE) %>% 
    cor() %>% 
    write.csv(paste0(d_rsl,"filted_bio_to_biom_correlation_matrix.csv"))

  species_distribution_model(d_rsl, disease_name, species_data,sdm_methods, biom)
  
}

# created_presence_pseudo_absence_points()


fitted_sdm_model("Alternaria blotch")












