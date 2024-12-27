library(tidyverse)
library(raster)


# baseline
read_history_projected_results <- function() {
  
  allfiels_paths <- fs::dir_ls(d_rsl, recurse = T,
                               glob = ".*_ensemble.tif")
  
  filepaths = str_extract(allfiels_paths, ".*(1970-2000).*\\.tif") %>% 
    na.omit()
  rastertemp = stack(filepaths) 
  names(rastertemp) = c('Alternaria_blotch')
  return(rastertemp)
}

# future
read_future_projected_results <- function() {
  
  allfiels_paths <- fs::dir_ls(d_rsl, recurse = T,
                               glob = ".*gcms_ensemble.tif")
  rastertemp = stack(allfiels_paths) 
  
  names(rastertemp) = str_extract(names(rastertemp), "\\d{4}.\\d{7}")
  return(rastertemp)
}

# calculated pixel area
calculate_pixel_area <- function(resolution, extent) {
  # extract resolution
  res_x <- resolution[1]
  res_y <- resolution[2]
  
  # deg to km
  deg_to_km <- function(degree, lat) {
    lat_rad <- pi * lat / 180
    km_per_degree <- 111.32 * cos(lat_rad)
    return(degree * km_per_degree)
  }
  
  # mean latitude
  mean_lat <- (extent[3] + extent[4]) / 2
  
  # pixel_width_km
  pixel_width_km <- deg_to_km(res_x, mean_lat)
  pixel_height_km <- 111.32 * res_y  
  
  # pixel_area_km2
  pixel_area_km2 <- pixel_width_km * pixel_height_km
  return(pixel_area_km2)
}

get_break_points <- function(d_rsl){
  # read threshold values
  model_evaluation = read.csv(paste0(d_rsl,'ensemble_sdms_model_evaluation.csv'))
  initial_value = mean(model_evaluation$Threshold)
  b = round((1 - initial_value) / 3, 2)
  breaks_points = c(0, initial_value[1], initial_value[1] + b[1], initial_value[1] + 2*b[1], 1)
  return(breaks_points)
}


processing_models_evaluation <- function(){
  model_evaluation_data = read.csv(paste0(d_rsl,"Alternaria blotch_model_evaluation.csv"))
  
  temp = model_evaluation_data %>% 
    dplyr::select(c("AUC","TSS","COR","threshold")) %>% 
    mutate(methods = "Ensemble")
  model_evaluation_data2 = model_evaluation_data %>%
    dplyr::select(c("methods","AUC","TSS","COR","threshold")) %>% 
    rbind(temp)
  write.csv(model_evaluation_data2,paste0(d_rsl,"included_ensemble_model_evaluation.csv"))
  
  ## save models evaluation results with ensemble sdms 
  model_evaluation_data2 %>% 
    dplyr::select(c("methods","AUC","TSS","COR","threshold")) %>% 
    group_by(methods) %>%
    summarize(
      AUC = round(mean(AUC, na.rm = TRUE), 2),
      TSS = round(mean(TSS, na.rm = TRUE), 2),
      COR = round(mean(COR, na.rm = TRUE), 2),
      Threshold = round(mean(threshold, na.rm = TRUE), 2)
    ) %>% 
    write.csv(paste0(d_rsl,"ensemble_sdms_model_evaluation.csv"))
}
# processing_models_evaluation()


calculated_bioclimatic_variables_importance <- function(){
  data <- read_csv(paste0(d_rsl,"Alternaria blotch_variable_importance.csv"))
  tempt1 = data %>%
    dplyr::select(starts_with("corTest")) %>%
    pivot_longer(cols = everything(),
                 names_to = c("test_name",".value"),
                 names_pattern = "(corTest)\\.(.*)") %>%
    rename(values = corTest) %>%
    arrange(desc(values))
  #
  tempt2 = data %>%
    dplyr::select(starts_with("AUCTest")) %>%
    pivot_longer(cols = everything(),
                 names_to = c("test_name",".value"),
                 names_pattern = "(AUCtest)\\.(.*)")%>%
    rename(values = AUCtest)%>%
    arrange(desc(values))

  data2 = rbind(tempt1, tempt2) %>%
    mutate(values = values*100,lower = lower *100,upper = upper *100) # change unit to
  
  # data <- data %>% 
  #   mutate(AUCtest = AUCtest*100,corTest = corTest*100) 
  
  ## save environment importance as table
  data2 %>% write.csv(paste0(d_rsl,"bioclimatic_variables_importance.csv"))
}
# calculated_bioclimatic_variables_importance()


calculated_response_recurves <- function(){
  data = read.csv(paste0(d_rsl,"Alternaria blotch_ReponseCurve_Value.csv"))
  
  # pivot to long
  data_long = data %>% pivot_longer(cols = everything(), names_to = c("variable",".value"), 
                                    names_sep = "\\.") %>% 
    pivot_longer(-c("variable","Response"),values_to = "x_bio") %>% 
    dplyr::select(c("variable","x_bio","Response")) %>% 
    set_names(c("variables","x_bio","response"))
  
  # delete NA rows
  df_filtered = data_long[complete.cases(data_long$x_bio), ]
  
  # Define the desired order of variables
  # desired_order = c("Bio9", "Elevation", "Bio4", "Bio14")
  
  ## change unit of bio4
  df_filtered = df_filtered %>%
    mutate(x_bio = ifelse(variables == "Bio4", x_bio * 0.01, x_bio))
  
  df_filtered %>% write.csv(paste0(d_rsl,"response_recurves.csv"))
  
  # cross_pints = df_filtered %>%
  #   subset(variables %in% c("Bio9", "Elevation", "Bio4", "Bio14")) %>%
  #   mutate(variables = factor(variables, levels = desired_order,
  #                             labels=c("Bio9(°C)", "Elevation(m)", "Bio4(°C)", "Bio14(mm)")),
  #          response = round(response,2)) %>%
  #   subset(response == 0.46) %>% 
  #   mutate(x_bio = round(x_bio,2)) 
  # 
  # max_response_ponit_bio9 = df_filtered %>%
  #   filter(variables == "Bio8") %>%
  #   slice(which.max(response)) %>%
  #   pull(x_bio)
}
# calculated_response_recurves()


risk_area_and_lon_lat <- function(){
  # get break points
  breaks_alternaria_blotch = get_break_points(d_rsl)
  
  #read history and future data
  history_results <- read_history_projected_results()
  future_resutls <- read_future_projected_results()
  
  #calculated area of pixel
  # check resolution
  # calculated  pixel area (km2)
  pixel_area_km2 <- calculate_pixel_area(res(history_results), extent(history_results))
  
  data <- as.data.frame(history_results, xy = T, na.rm = T) %>%
    rename("lon" = x, "lat" = y) %>%
    mutate(values = cut(Alternaria_blotch, breaks = breaks_alternaria_blotch,
                        labels = c('US', 'LS', 'MS', 'HS')),
           periods_ssps = "Baseline") %>% 
    filter(values != "US") %>% 
    dplyr::select(c(lat, values)) 

  
  data$lat_group <- cut(data$lat, breaks = seq(21, ceiling(max(data$lat)) + 5, by = 5),
                        include.lowest = TRUE)
  
    df <- data %>%
    group_by(lat_group, values) %>%
    summarize(count = n(), .groups = 'drop') %>% 
    pivot_wider(names_from = values, values_from = count) %>% 
    mutate(Total = rowSums(across(c(LS, MS, HS)), na.rm = TRUE)) %>% 
    pivot_longer(-lat_group, names_to = "values", values_to = "count")
  
  
  
  ssps <- c("SSP126", "SSP245", "SSP370", "SSP585")

  expanded_df <- df[rep(1:nrow(df), each = 4), ]
  expanded_df$ssps <- rep(ssps, times = nrow(df))
  final_df <- expanded_df[, c("ssps", "values", "lat_group","count")]
  
  baseline_results <- final_df %>% 
    mutate(periods = "Baseline")
  
  future_data <- as.data.frame(future_resutls, xy = T, na.rm = T) %>%
    rename("lon" = x, "lat" = y) %>%
    mutate_at(vars(-lon, -lat), ~cut(., breaks = breaks_alternaria_blotch,
                                     labels = c('US','LS','MS','HS'))) %>%
    pivot_longer(-c(lon, lat), names_to = "times",values_to = "values",values_drop_na = TRUE) %>%
    mutate(periods = str_extract(times,"(?<=X).{9}") %>% str_replace("\\.", "-"),
           ssps = paste0("SSP", str_extract(times, "\\d{3}$")),.before = 3) %>% 
    dplyr::select(c("lon","lat","periods","ssps","values")) %>% 
    mutate(periods = case_when(
      periods == "2021-2040" ~ "2030s",
      periods == "2041-2060" ~ "2050s",
      periods == "2061-2080" ~ "2070s",
      periods == "2081-2100" ~ "2090s",
    )) %>% 
    filter(values != "US") %>% 
    dplyr::select(c(lat, periods, ssps, values))
   
  future_data$lat_group <- cut(future_data$lat,
                        breaks = seq(21, ceiling(max(future_data$lat)) + 5, by = 5),
                        include.lowest = TRUE)
  
  future_result <- future_data %>%
    group_by(periods, ssps, lat_group, values) %>%
    summarize(count = n(), .groups = 'drop') %>% 
    pivot_wider(names_from = values, values_from = count) %>% 
    mutate(Total = rowSums(across(c(LS, MS, HS)), na.rm = TRUE)) %>% 
    pivot_longer(-c(periods, ssps, lat_group), names_to = "values", values_to = "count")

  results <- bind_rows(future_result, baseline_results) %>% 
    filter(values == "Total") %>% 
    mutate(periods = factor(periods, 
                            levels = c("Baseline", "2030s", "2050s", "2070s", "2090s"))) %>%
   group_by(periods, lat_group) %>%
   summarize(mean_count = mean(count), se = sd(count), .groups = 'drop')

  p = ggplot(results, aes(x = lat_group, y = mean_count, fill = periods)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                  position = position_dodge(0.9), width = 0.25) +
    labs(x = "Lat ", y = "Total risk", title = "") +
    theme_test()
  
  ggsave("./fig/lat_risk_area.jpg", p, width = 6, height = 5, dpi = 300)

}



calculated_area_diff_periods_and_regions <- function(){
  ## all reagions
  caculate_all_regions_area <- function(){
    # get break points
    breaks_alternaria_blotch = get_break_points(d_rsl)
    
    #read history and future data
    history_results <- read_history_projected_results()
    future_resutls <- read_future_projected_results()
    
    #calculated area of pixel
    # check resolution
    # calculated  pixel area (km2)
    pixel_area_km2 <- calculate_pixel_area(res(history_results), extent(history_results))
    
    history_area <- as.data.frame(history_results, xy = T, na.rm = T) %>%
      rename("lon" = x, "lat" = y) %>%
      mutate(values = cut(Alternaria_blotch, breaks = breaks_alternaria_blotch,
                          labels = c('US', 'LS', 'MS', 'HS')),
             periods_ssps = "Baseline") %>%
      dplyr::select(c("periods_ssps","values")) %>% 
      # remove unsuitable lable
      # filter(values != "US") %>%
      group_by(periods_ssps, values) %>%
      summarise(areas = n()*pixel_area_km2, .groups = 'drop') %>% 
      mutate(areas = round((areas) *0.0001, 2))   # units are 10^4 km2 
    
    
    future_diff_area <- as.data.frame(future_resutls, xy = T, na.rm = T) %>%
      rename("lon" = x, "lat" = y) %>%
      mutate_at(vars(-lon, -lat), ~cut(., breaks = breaks_alternaria_blotch,
                                       labels = c('US','LS','MS','HS'))) %>%
      pivot_longer(-c(lon, lat), names_to = "times",values_to = "values",values_drop_na = TRUE) %>%
      mutate(periods = str_extract(times,"(?<=X).{9}") %>% str_replace("\\.", "-"),
             ssps = paste0("SSP", str_extract(times, "\\d{3}$")),.before = 3) %>% 
      dplyr::select(c("periods","ssps","values")) %>% 
      mutate(periods = case_when(
        periods == "2021-2040" ~ "2030s",
        periods == "2041-2060" ~ "2050s",
        periods == "2061-2080" ~ "2070s",
        periods == "2081-2100" ~ "2090s",
      )) %>% 
      group_by(periods, ssps, values) %>%
      summarise(areas = n()*pixel_area_km2, .groups = 'drop') %>% 
      mutate(areas = round((areas) *0.0001, 2)) %>%    # units are 10^4 km2 
      mutate(periods_ssps = paste(ssps, periods, sep = "_")) %>% 
      filter(periods_ssps %in% c("SSP126_2050s","SSP245_2050s","SSP370_2050s","SSP585_2050s",
                                 "SSP126_2090s","SSP245_2090s","SSP370_2090s","SSP585_2090s")) %>%
      dplyr::select(c("periods_ssps","values","areas"))
    
    
    sutiable_area <- bind_rows(history_area,future_diff_area)  
    write.csv(sutiable_area, paste0(d_rsl,"sutiable_area_2050s_2090s.csv"),row.names = FALSE)
    
    
    
    ## calculate table2 
    future_diff_area_table2 <- as.data.frame(future_resutls, xy = T, na.rm = T) %>%
      rename("lon" = x, "lat" = y) %>%
      mutate_at(vars(-lon, -lat), ~cut(., breaks = breaks_alternaria_blotch,
                                       labels = c('US','LS','MS','HS'))) %>%
      pivot_longer(-c(lon, lat), names_to = "times",values_to = "values",values_drop_na = TRUE) %>%
      mutate(periods = str_extract(times,"(?<=X).{9}") %>% str_replace("\\.", "-"),
             ssps = paste0("SSP", str_extract(times, "\\d{3}$")),.before = 3) %>% 
      dplyr::select(c("periods","ssps","values")) %>% 
      mutate(periods = case_when(
        periods == "2021-2040" ~ "2030s",
        periods == "2041-2060" ~ "2050s",
        periods == "2061-2080" ~ "2070s",
        periods == "2081-2100" ~ "2090s",
      )) %>% 
      group_by(periods, ssps, values) %>%
      summarise(areas = n()*pixel_area_km2, .groups = 'drop') %>% 
      mutate(areas = round((areas) *0.0001, 2)) %>%    # units are 10^4 km2 
      dplyr::select(c("periods","ssps","values","areas")) 
    
    future_diff_area_table2_results <- future_diff_area_table2 %>% 
      pivot_wider(id_cols = c(periods, ssps), names_from = values, values_from = areas) %>% 
      mutate(Total = LS + MS + HS) %>% 
      dplyr::select(c(periods, ssps, LS, MS, HS, Total))
      
    
    history_area_results <- history_area %>% 
      pivot_wider(id_cols = periods_ssps,names_from = values, values_from = areas) %>% 
      mutate(Total = LS + MS + HS) %>% 
      dplyr::select(c(periods_ssps,LS, MS, HS, Total)) %>% 
      rename(ssps = periods_ssps) %>% 
      mutate(periods = "1970-2000",.before = 2)
      
    sutiable_area <- bind_rows(history_area_results,future_diff_area_table2_results)  
    write.csv(sutiable_area, paste0(d_rsl,"table2_different_periods_suitable_area.csv"),
              row.names = FALSE)
  
  }
  
  ## five apple planting areas
  
  caculate_diff_planting_area <- function(){
    
    # get break points
    breaks_alternaria_blotch = get_break_points(d_rsl)
    
    #read history and future data
    history_results <- read_history_projected_results()
    future_resutls <- read_future_projected_results()
    
    apple_planting_zones = stack('./data/apple_planting_zones/five_zones/apple_planting_five_zones.tif') %>% 
      resample(history_results, method = "ngb")
    
    ## change values dtype to int
    values(apple_planting_zones) = as.integer(values(apple_planting_zones))
    
    baseline_data1 = stack(history_results, apple_planting_zones) %>% 
      as.data.frame(xy = FALSE) %>% 
      drop_na() %>% 
      rename("planting"="apple_planting_five_zones")
    
    # ## the values used for suitable area and unsuitable area
    
    baseline_data2 = baseline_data1 %>% 
      mutate(values = cut(Alternaria_blotch, breaks = breaks_alternaria_blotch,
                          labels = c('US', 'LS', 'MS', 'HS')),
             periods_ssps = "Baseline") %>%
      dplyr::select(c("planting","periods_ssps","values")) %>% 
      
      group_by(planting, periods_ssps, values) %>%
      summarise(areas = n()*pixel_area_km2, .groups = 'drop') %>% 
      mutate(areas = round((areas) *0.0001, 2))    # units are 10^4 km2 
    
    
    ### future
    future_data1 = stack(future_resutls, apple_planting_zones) %>% 
      as.data.frame(xy = T) %>% 
      rename("lon" = x, "lat" = y) %>%
      drop_na() %>% 
      rename("planting"="apple_planting_five_zones")
    
    # ## the values used for suitable area and unsuitable area
    
    future_data2 = future_data1 %>% 
      mutate_at(vars(-lon, -lat,-planting), ~cut(., breaks = breaks_alternaria_blotch,
                                                 labels = c('US','LS','MS','HS'))) %>%
      pivot_longer(-c(lon, lat, planting), names_to = "times",values_to = "values",values_drop_na = TRUE) %>%
      mutate(periods = str_extract(times,"(?<=X).{9}") %>% str_replace("\\.", "-"),
             ssps = paste0("SSP", str_extract(times, "\\d{3}$")),.before = 3) %>% 
      dplyr::select(c("planting","periods","ssps","values")) %>% 
      mutate(periods = case_when(
        periods == "2021-2040" ~ "2030s",
        periods == "2041-2060" ~ "2050s",
        periods == "2061-2080" ~ "2070s",
        periods == "2081-2100" ~ "2090s",
      )) %>% 
      group_by(planting,periods, ssps, values) %>%
      summarise(areas = n()*pixel_area_km2, .groups = 'drop') %>% 
      mutate(areas = round((areas) *0.0001, 2)) %>%    # units are 10^4 km2 
      mutate(periods_ssps = paste(ssps, periods, sep = "_")) %>% 
      filter(periods_ssps %in% c("SSP126_2050s","SSP245_2050s","SSP370_2050s","SSP585_2050s",
                                 "SSP126_2090s","SSP245_2090s","SSP370_2090s","SSP585_2090s")) %>%
      dplyr::select(c("planting","periods_ssps","values","areas"))
    
    
    sutiable_area_future <- bind_rows(baseline_data2, future_data2)  
    write.csv(sutiable_area_future, paste0(d_rsl,"diff_planting_sutiable_area_2050s_2090s.csv"),
              row.names = FALSE)
  }
  
  
  caculate_all_regions_area()
  each_region <- caculate_diff_planting_area()
  
  ## summary all results
  
  all_regions <- read_csv(paste0(d_rsl,"sutiable_area_2050s_2090s.csv")) %>% 
    filter(values != "US") %>% 
    mutate(planting = "All_regions")
  
  each_region <- read_csv(paste0(d_rsl,"diff_planting_sutiable_area_2050s_2090s.csv")) %>% 
    filter(values != "US") %>% 
    mutate(planting = as.character(planting))
  combined_data <- bind_rows(all_regions, each_region) %>% 
    pivot_wider(names_from = values,values_from = areas) 
  write.csv(combined_data, paste0(d_rsl, "diff_planting_sutiable_area_2050s_2090s_summary.csv"))
  
}
# calculated_area_diff_periods_and_regions()


calculated_diff_planting_prec_temperature_dem <- function(){
  china <- read_sf("./data/bounds/国界/国家矢量.shp")


  dem_df = raster("./data/dem_5km.tif") 
  
  
  apple_planting_zones <- stack('./data/apple_planting_zones/five_zones/apple_planting_five_zones.tif') %>% 
    resample(dem_df, method = "ngb")
  
  # change values dtype to int
  values(apple_planting_zones) <- as.integer(values(apple_planting_zones))
  
  results <- stack(dem_df, apple_planting_zones) %>% 
    as.data.frame(xy = FALSE) %>% 
    drop_na() %>% 
    rename("planting" = "apple_planting_five_zones") %>% 
    group_by(planting) %>% 
    summarise(mean_dem_5km = mean(dem_5km, na.rm = TRUE)) 

  ###prec
  
  prec <- raster("./data/china_mean_prec_tem/prec/annual_prec_sum.tif") %>% 
    mask(china)
  apple_planting_zones <- stack('./data/apple_planting_zones/five_zones/apple_planting_five_zones.tif') %>% 
    resample(prec, method = "ngb")
  
  # change values dtype to int
  values(apple_planting_zones) <- as.integer(values(apple_planting_zones))
  
  results <- stack(prec, apple_planting_zones) %>% 
    as.data.frame(xy = FALSE) %>% 
    drop_na() %>% 
    rename("planting" = "apple_planting_five_zones") %>% 
    group_by(planting) %>% 
    summarise(mean_annual_prec_sum = mean(annual_prec_sum , na.rm = TRUE)) 
  
  #### temper

  temper <- raster("./data/china_mean_prec_tem/temperature/mean_annual_temperature.tif") %>% 
    mask(china) 

  apple_planting_zones <- stack('./data/apple_planting_zones/five_zones/apple_planting_five_zones.tif') %>% 
    resample(temper, method = "ngb")
  
  # change values dtype to int
  values(apple_planting_zones) <- as.integer(values(apple_planting_zones))
  
  results <- stack(temper, apple_planting_zones) %>% 
    as.data.frame(xy = FALSE) %>% 
    drop_na() %>% 
    rename("planting" = "apple_planting_five_zones") %>% 
    group_by(planting) %>% 
    summarise(mean_mean_annual_temperature  = mean(mean_annual_temperature  , na.rm = TRUE))

}


calculated_diff_levels_mean_area_change_every_20_years <- function(){

  data <- read_csv(paste0(d_rsl,"table2_different_periods_suitable_area.csv"))
  mean_change_data <- data %>% 
    dplyr::select(c("ssps","periods","LS","MS","HS","Total")) %>% 
    filter(periods %in% c("1970-2000","2081-2100"))
  
  # Extract historical values
  history_values <- mean_change_data %>% filter(ssps == "history")
  
  # Calculate the changes
  mean_change_data_changes <- mean_change_data %>%
    filter(ssps != "history") %>%
    mutate(
      LS_change = (LS - history_values$LS)/5,
      MS_change = (MS - history_values$MS)/5,
      HS_change = (HS - history_values$HS)/5,
      Total_change = (Total - history_values$Total)/5
    )
  
  # Select relevant columns for the new data frame
  mean_change_data_changes <- mean_change_data_changes %>%
    dplyr::select(ssps, LS_change, MS_change, HS_change, Total_change)
  
  # Transform data to long format
  long_data <- pivot_longer(mean_change_data_changes, cols = -ssps, names_to = "Category",
                            values_to = "Value")
  long_data$Category <- factor(long_data$Category,
                               levels = c("Total_change","LS_change","MS_change","HS_change"))
  long_data %>% write_csv(paste0(d_rsl,"diff_levels_mean_area_change_every_20_years.csv"))
}


calculated_diff_levels_area_change_from_baseline_to_2090s <- function(){
  
  data <- read_csv(paste0(d_rsl,"table2_different_periods_suitable_area.csv"))
  change_data <- data %>% 
    dplyr::select(c("ssps","periods","LS","MS","HS","Total")) %>% 
    filter(periods %in% c("1970-2000","2090s"))
  
  # Extract historical values
  history_values <- change_data %>% filter(ssps == "Baseline")
  
  ## Calculate the changes rate
  change_data_changes <- change_data %>%
    filter(ssps != "Baseline") %>%
    mutate(
      LS_change = (LS - history_values$LS)/history_values$LS,
      MS_change = (MS - history_values$MS)/history_values$MS,
      HS_change = (HS - history_values$HS)/history_values$HS,
      Total_change = (Total - history_values$Total)/history_values$Total
    )
  
  # Select relevant columns for the new data frame
  change_data_changes <- change_data_changes %>%
    dplyr::select(ssps, LS_change, MS_change, HS_change, Total_change)
  
  name_csv <- "diff_levels_area_change_diff_levels_area_change_from_baseline_to_2090s_rate.csv"

  change_data_changes %>% 
  write_csv(paste0(d_rsl,name_csv))
  
  
  ## Calculate the changes
  change_data_changes <- change_data %>%
    filter(ssps != "Baseline") %>%
    mutate(
      LS_change = (LS - history_values$LS),
      MS_change = (MS - history_values$MS),
      HS_change = (HS - history_values$HS),
      Total_change = (Total - history_values$Total)
    )
  
  # Select relevant columns for the new data frame
  change_data_changes <- change_data_changes %>%
    dplyr::select(ssps, LS_change, MS_change, HS_change, Total_change)
  
  name_csv <- "diff_levels_area_change_diff_levels_area_change_from_baseline_to_2090s.csv"
  
  change_data_changes %>% 
    write_csv(paste0(d_rsl,name_csv))
  
  
  
  
}




zhongxin_change <- function(){
  library(terra)
  
  # Formula for calculating the center of gravity
  
  gravity_xy = function(p, x, y){
    gx = sum(p*x)/sum(p)
    gy = sum(p*y)/sum(p)
    gxy = tibble(gx, gy)
    gxy
  }
  
  allfiels_paths <- fs::dir_ls(d_rsl, recurse = T,
                               glob = ".*_ensemble.tif")
  rastertemp = stack(allfiels_paths) 
  # names(rastertemp) = c('Alternaria_blotch')
  
  diseases_gxy <- terra::as.data.frame(rastertemp, xy = T)%>% 
    pivot_longer(-c(x, y), 
                 names_to = "times",
                 values_to = "values",
                 values_drop_na = TRUE) %>%
    
    group_by(times) %>%   
    summarize(gxy = list(gravity_xy(values, x, y))) %>%  
    unnest(gxy) 
  
  baseline_gxy <- diseases_gxy %>% 
    mutate(periods = str_extract(times,"(?<=0.8.).{9}") %>% str_replace("\\.", "-"),
           .before = 1) %>% 
    filter(periods == "1970-2000") %>% 
    dplyr::select(periods, gx, gy) 
  
  ssps_values <- c("SSP126", "SSP245", "SSP370", "SSP585")
  
  baseline_gxy_expanded <- baseline_gxy %>%
    tidyr::expand(periods, ssps = ssps_values) %>%
    mutate(gx = 108, gy = 35.4)
  
  future_gxy <- diseases_gxy %>% 
    mutate(periods = str_extract(times,"(?<=0.8.).{9}") %>% str_replace("\\.", "-"),
           .before = 1) %>% 
    filter(periods != "1970-2000") %>% 
    mutate(ssps = paste0("SSP", str_extract(times, "\\d{3}(?=_ensemble)")),
           .before = 2) %>%
    dplyr::select(periods,ssps, gx, gy)
  
  combined_gxy <- bind_rows(baseline_gxy_expanded, future_gxy)
  write.csv(combined_gxy, file = paste0(d_rsl, "gxy_change.csv"),row.names = FALSE) 
  
  
  ### calculate remove the distance
  
  # read DEM data
  dem_df = raster("./data/dem_5km.tif") %>%
    projectRaster(crs=4326)
  
  
  # Extract DEM values based on gx and gy coordinates
  dem_values <- extract(dem_df, cbind(combined_gxy$gx, combined_gxy$gy))
  
  # Add the DEM values to the tibble
  combined_gxy <- combined_gxy %>%
    mutate(dem_value = dem_values)
  
  library(sf)

  combined_gxy_gxy_sf <- combined_gxy %>%
    st_as_sf(coords = c("gx", "gy"), crs = 4326)

  
  # Function to calculate distances between consecutive points in each group
  calculate_distances <- function(data) {
    if (nrow(data) <= 1) {
      return(data.frame(periods = data$periods, distance = NA_real_,
                        diff_dem = NA_real_))   # Return NA if there is only one record
    }
    
    # Calculate distances between each point and its predecessor
    distances <- st_distance(data$geometry)
    # Extract the distances between consecutive points from the matrix
    distance_vector <- distances[row(distances) == col(distances) - 1]
    
    
    # Calculate change of dem values
    dem_vector <- data$dem_value - lag(data$dem_value)
    
    # Append calculated distances to the dataframe
    data <- data %>%
      mutate(distance = c(NA, as.numeric(distance_vector)),
             diff_dem = c(as.numeric(dem_vector))) # Assign NA to the first period without a predecessor
    
    return(data)
  }
  
  # Process the entire dataset, calculate distances for each ssps group, and select relevant columns
  results <- combined_gxy_gxy_sf %>%
    group_by(ssps) %>%
    arrange(ssps, periods) %>%  # Ensure each group is sorted by time
    group_modify(~ calculate_distances(.x)) %>% 
    filter(periods!="1970-2000") %>% 
    dplyr::select(ssps, periods, distance,diff_dem)
  
  write.csv(results, file = paste0(d_rsl, "gxy_change_distance.csv"),row.names = FALSE) 
  
}


uncertain_analysis <- function(){
  ## baseline suitable area
  allfiels_paths_ana <- fs::dir_ls(d_rsl, recurse = T,
                                   glob = ".*_Binary.tif")
  baseline_raster = str_extract(allfiels_paths_ana, ".*(1970-2000).*\\.tif") %>% 
    na.omit() %>% 
    .[!str_detect(., "ensemble_Binary.tif")] %>% 
    stack() 

  baseline_raster_ana = baseline_raster %>% 
    as.data.frame(xy = FALSE) %>% 
    drop_na() 
  
  baseline_data_ana = baseline_raster_ana %>% 
    # group_by(apple_panlting_reagions6_5km) %>% 
    summarise(across(everything(),~sum(. == 1),.names = "{.col}__1")) %>% 
    pivot_longer(everything(), names_to = c("species_times",".value"),
                 names_sep = "__", values_drop_na = T) %>%
    mutate(sdms = str_extract(species_times,"(?<=\\d{4}_)[A-Za-z]+(?=_)"),.before = 2) %>% 
    dplyr::select(c("sdms","1")) %>%
    set_names("sdms", "baseline_area")

  
  ### future 
  allfiels_paths_sgd = fs::dir_ls(d_rsl, recurse = T,
                                  glob = ".*ssp_gcm_sdm_Binary.tif")
  sgd_raster = allfiels_paths_sgd %>% 
    na.omit() %>% 
    stack() 
  sgd_raster2 = sgd_raster %>% 
    as.data.frame(xy = F, na.rm = T)  
  
  data_sgd = sgd_raster2 %>% 
    summarise(across(everything(),~sum(. == 5),.names = "{.col}__5")) %>% 
    pivot_longer(everything(), names_to = c("species_times",".value"),
                 names_sep = "__", values_drop_na = T) %>%
    mutate(
      # times = str_extract(species_times,"\\d+\\.\\d{4}") %>% str_replace("\\.", "-"),
      times = str_extract(species_times,"(?<=0.8.)\\d+\\.\\d{4}") %>% str_replace("\\.", "-"),
      ssps = str_extract(species_times,"(?<=.\\d{4})\\d{3}"),
      sdms = str_extract(species_times,"[A-Za-z]+(?=_ssp_gcm_sdm)"),
      gcms = str_extract(species_times,"(?<=\\d{7}_)[A-Za-z0-9.]+(?=_)"),.before = 1) %>% 
    dplyr::select(c("times","ssps","sdms","gcms","5")) %>%
    set_names("times","ssps","sdms","gcms","area")  

  Three_ways_anova_analisis = function(data){
    data$sdms = factor(data$sdms)
    data$gcms = factor(data$gcms)
    data$ssps = factor(data$ssps)
    
    # three-way ANOVA
    model = lm(area ~ ssps * sdms * gcms, data=data)
    anova_result = anova(model)
    return(anova_result$`Sum Sq`)
  }
  result = data_sgd %>%
    group_nest(times) %>%
    mutate(anova_result = map(data, Three_ways_anova_analisis))
  column_names = c("Scens", "SDMs", "GCMs", "Scens×SDMs",
                   "Scens×GCMs", "SDMs×GCMs", "Scens×SDMs×GCMs", "Residuals")
                                                                                                                                                                                                                                                                                       
  finally_data = result %>%
    dplyr::select(times, anova_result) %>%
    unnest(anova_result) %>% 
    mutate(uncertainty_ele = rep(column_names, times = 4)) %>% 
    pivot_wider(names_from = uncertainty_ele,values_from =  anova_result) %>% 
    dplyr::select(c("times","Scens", "SDMs", "GCMs", "Scens×SDMs",
                    "Scens×GCMs", "SDMs×GCMs", "Scens×SDMs×GCMs")) %>% 
    mutate(total = rowSums(dplyr::select(., starts_with("Scens"), starts_with("SDMs"),
                                         starts_with("GCMs")))) %>% 
    # calculating probation
    mutate_at(vars( -times, -total), ~ ./total * 100) %>%
    dplyr::select(-total) %>% 
    pivot_longer(-c(times), names_to = "uncertainty", values_to = "contribution")
  finally_data$uncertainty = factor(finally_data$uncertainty,
                                    levels = c("SDMs", "Scens", "GCMs", "Scens×SDMs",
                                               "Scens×GCMs", "SDMs×GCMs", "Scens×SDMs×GCMs")) 
  finally_data %>% write_csv(paste0(d_rsl,"uncertain_analysis_results.csv"))
  
  
  # mean contribution
  average_contribution <- finally_data %>%
    group_by(uncertainty) %>%
    summarize(average_contribution = mean(contribution))
  
  average_contribution %>% 
    write_csv(paste0(d_rsl,"mean_uncertain_analysis_results.csv"))
}

results_analysis <- function(d_rsl){
  
  processing_models_evaluation()
  calculated_bioclimatic_variables_importance()
  calculated_response_recurves()
  caculated_diff_levels_area_change()
  calculated_diff_levels_mean_area_change_every_20_years()
  uncertain_analysis()
}


results_analysis("./results/base_cor_0.8/")


d_rsl <- "./results/base_cor_0.8/"








