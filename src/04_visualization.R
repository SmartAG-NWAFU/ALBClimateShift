rm(list = ls())

library(tidyverse)
library(sf)
library(raster)
library(terra)
library(ggspatial)
library(extrafont)

# global font size
family_bont = "Times New Roman"
base_size = 12


d_rsl <- "./results/base_cor_0.8/"


# font_import()  
loadfonts()

albers = "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

china_province = read_sf("./data/bounds//省级/省级.shp")
eng_name_province <- read_csv("./data/bounds/省级/provinces_name_english.csv")
china_province <- left_join(china_province,eng_name_province,by =  c("省" = "CNAME"))

nineline = read_sf("./data/bounds/九段线/九段线.shp")
china <- read_sf("./data/bounds/国界/国家矢量.shp")
apple_planting = read_sf('./data/apple_planting_zones/five_zones/apple_planting_area.shp')


get_break_points <- function(d_rsl){
  # read threshold values
  model_evaluation = read.csv(paste0(d_rsl,'ensemble_sdms_model_evaluation.csv'))
  initial_value = mean(model_evaluation$Threshold)
  b = round((1 - initial_value) / 3, 2)
  breaks_points = c(0, initial_value[1], initial_value[1] + b[1], initial_value[1] + 2*b[1], 1)
  return(breaks_points)
}


fig_ALB_records <- function(){
  
  library(RColorBrewer)
  library(ggpattern)
  # library(Cairo)
  library(viridis)
  library(cowplot)
  
  
  alternaria_blotch_records <- read_csv("./data/species/Alternaria blotch_thin1.csv")
  
  ## transform projection
  sp_datas = st_as_sf(alternaria_blotch_records, coords = c("lon", "lat"), crs = 4326) %>% 
    st_transform(crs = albers) 
  
  ## read dem data 
  dem_df = raster("./data/dem_5km.tif") %>% 
    projectRaster(crs=albers) %>%
    as.data.frame(xy=T)
  
  
  plot = ggplot() +
    # geom_raster(data = dem_df, aes(x = x, y = y, fill = dem_5km), interpolate = TRUE, alpha = 0.7) +
    # scale_fill_gradientn(colors = colors, na.value = "white") +
    # geom_sf(data = china_province, fill = NA, colour = "black", linewidth = 0.25) +
    geom_sf(data = nineline, colour = "black", linewidth = 0.25) +
    geom_sf(data = china, fill = "#E6E6E6", colour = "black", linewidth = 0.25) +
    geom_sf(data = apple_planting, aes(fill=planting), colour = "black", linewidth = 0.3) +
    geom_point(data = sp_datas, aes(x = st_coordinates(geometry)[, "X"],
                                    y = st_coordinates(geometry)[, "Y"],
                                    color = "Presence points (n=261)"), size = 2.5) +
    
    scale_color_manual(values = c("black" = "black", "Presence points (n=261)" = "red")) +
    scale_fill_manual(values = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6"),
                      labels = c("Sub-region I (the Loess Plateau)", 
                                 "Sub-region ll (the Bohai Bay)", 
                                 "Sub-region lll (the Southwest Cool Highland)", 
                                 "Sub-region lV (the old course of the Yellow River)", 
                                 "Sub-region V (the Specialty region)")) +
    geom_sf_text(data = apple_planting, aes(label = planting), size = 8) +
    coord_sf(crs = albers, 
             ylim = c(1500000, 6000000),
             xlim = c(-3100000, 2000000),
             expand = FALSE) +
    labs(x = NULL, y = NULL, fill = "", color = "") +
    theme_bw() +
    theme(axis.text = element_text("Times New Roman", color = "black",size=base_size),
          axis.title = element_text("Times New Roman", color = "black",size=base_size),
          strip.text = element_text("Times New Roman", color = "black",size=base_size),
          legend.text = element_text("Times New Roman", color = "black",size=base_size),
          legend.position = c(0.01, 0.39),
          legend.justification = c(0, 1),
          legend.box = "vertical",
          # legend.spacing = unit(0.05,"cm"),  # Adjust the space between legends
          # legend.margin = margin(t = 0.1, unit = "cm")  # Adjust the margin around the legend
          legend.background = element_rect(fill = alpha("white", 0))
    ) +  
    
    annotation_scale(location = "bl",pad_y = unit(0.35, "cm"),text_cex = 1.2, text_family = family_bont) + 
    annotation_north_arrow(height = unit(2.5, "cm"),width = unit(2.5, "cm"),
                           location = "tl", style = north_arrow_nautical(
                             fill = c("grey35", "white"), line_col = "grey15"))
  
  nanhai = ggplot() +
    geom_sf(data = nineline, colour = "black", linewidth = 0.25) +
    geom_sf(data = china, fill = "#E6E6E6", colour = "black", linewidth = 0.25) +
    coord_sf(ylim = c(-4028017,-1877844),xlim = c(117131.4,2115095),crs="+proj=laea +lat_0=40 +lon_0=104")+
    theme(legend.position = "",
          aspect.ratio = 1.25, 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="grey10",linetype=1,size=0.5))
  
  fig = ggdraw() +
    draw_plot(plot) +
    draw_plot(nanhai, x = 0.8, y = 0.005, width = 0.12, height = 0.3)
  
  # Save the plot as a JPG file
  ggsave("./fig/fig1_apple_disease_records.jpg",fig, width = 10, height = 8, dpi = 300)
  
  
  
  ## DEM
  colors = c("#33A02C", "#B2DF8A", "#FDBF6F", "#1F78B4", "#999999", "#E31A1C", "#E6E6E6", "#A6CEE3")
  dem_map <- ggplot() +
    geom_raster(data = dem_df, aes(x = x, y = y, fill = dem_5km), interpolate = TRUE) +
    scale_fill_gradientn(colors = colors, 
                         na.value = "white",
                         # limits = c(0, 8800),
                         name = "DEM(m)"
    ) +
    geom_sf(data = china, fill = NA, colour = "black", linewidth = 0.25)+
    coord_sf(crs = albers, 
             ylim = c(1500000, 6000000),
             xlim = c(-3100000, 2000000),
             expand = FALSE) +
    theme_test() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text("Times New Roman", color = "black",size=6),
      legend.text = element_text("Times New Roman", color = "black",size=5),
      legend.position = c(0.15, 0.2),
      legend.key.size = unit(0.2, "cm") 
    )
  
  # Save the plot as a JPG file
  ggsave("./fig/fig1_dem_map.jpg",dem_map , width = 2.5, height = 2, dpi = 300)
  
  
  ## mean prec
  
  library(raster)
  library(viridis)
  library(RColorBrewer)
  
  
  colors <- c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")
  
  prec <- raster("./data/china_mean_prec_tem/prec/annual_prec_sum.tif")
  prec_chian <- mask(prec, china) %>% projectRaster(crs=albers) %>%
    as.data.frame(xy=T)
  
  prec_map <- ggplot() +
    geom_raster(data = prec_chian, aes(x = x, y = y, fill = annual_prec_sum), interpolate = TRUE) +
    scale_fill_gradientn(colors = colors, 
                         na.value = "white",
                         limits = c(16, 3200),
                         name = "MAP(mm)") +
    geom_sf(data = china, fill = NA, colour = "black", linewidth = 0.25)+
    coord_sf(crs = albers, 
             ylim = c(1500000, 6000000),
             xlim = c(-3100000, 2000000),
             expand = FALSE) +
    theme_test() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text("Times New Roman", color = "black",size=6),
      legend.text = element_text("Times New Roman", color = "black",size=5),
      legend.position = c(0.15, 0.2),
      legend.key.size = unit(0.2, "cm") 
    )
  
  # Save the plot as a JPG file
  ggsave("./fig/fig1_prec_map.jpg",prec_map, width = 2.5, height = 2, dpi = 300)
  
  
  ## mean temperature
  
  temper <- raster("./data/china_mean_prec_tem/temperature/mean_annual_temperature.tif")
  temper_chian <- mask(temper, china) %>% projectRaster(crs=albers) %>%
    as.data.frame(xy=T)
  
  temper_map <- ggplot() +
    geom_raster(data = temper_chian, aes(x = x, y = y, fill = mean_annual_temperature), interpolate = TRUE) +
    scale_fill_gradientn(
      colors = rev(brewer.pal(11, "RdYlGn")), 
      na.value = "white", 
      limits = c(-20, 25),
      name = expression(MAT ~ (degree*C))
    ) +
    geom_sf(data = china, fill = NA, colour = "black", linewidth = 0.25) +
    coord_sf(crs = albers, 
             ylim = c(1500000, 6000000),
             xlim = c(-3100000, 2000000),
             expand = FALSE) +
    theme_test() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text("Times New Roman", color = "black",size=6),
      legend.text = element_text("Times New Roman", color = "black",size=5),
      legend.position = c(0.15, 0.2),
      legend.key.size = unit(0.2, "cm") 
    )
  
  # Save the plot as a JPG file
  ggsave("./fig/fig1_temper_map.jpg",temper_map, width = 2.5, height = 2, dpi = 300)
  
}



fig_baseline_map <- function(){
  
  library(dplyr)
  library(RColorBrewer)
  library(ggpattern)
  library(Cairo)
  library(viridis)
  library(cowplot)
  library(ggspatial)
  
  
  # read model simulation
  allfiels_paths <- fs::dir_ls(d_rsl, recurse = T,
                               glob = ".*_ensemble.tif")
  
  filepaths = str_extract(allfiels_paths, ".*(1970-2000).*\\.tif") %>% 
    na.omit()
  
  rastertemp = terra::rast(filepaths) %>% 
    project(albers, method = "bilinear")
  names(rastertemp) = c('Alternaria_blotch')
  
  
  breaks_Alternaria_blotch = get_break_points(d_rsl)
  
  # colnames(change_ensemble)
  rastertemp_discrese = as.data.frame(rastertemp, xy = T, na.rm = T) %>% 
    rename("lon" = x, "lat" = y)%>% 
    mutate(Alternaria_blotch = cut(Alternaria_blotch, breaks = breaks_Alternaria_blotch,
                                   labels = c('US', 'LS', 'MS', 'HS'))) %>% 
    pivot_longer(cols = -c(lon,lat), names_to = "species", values_to = "class")
  
  
  class_area <- rastertemp_discrese %>%
    count(class) %>%
    mutate(area = round((n *18.86954) *0.0001, 2)) %>% 
    subset(class %in% c("LS", "MS", "HS")) 
  
  
  # bar_statistic_areas
  bar_statistic_areas <- ggplot(data = class_area, aes(x = class, y = area, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
    scale_fill_manual(values = c("LS" = "#4daf4a", "MS" = "#377eb8", "HS" = "#e41a1c")) +
    labs(x = "", y = bquote(Area~(10^4~Km^2)), title = "") +
    theme_test() +
    theme(
      text = element_text(family = "Times New Roman", color = "black"),  
      axis.text = element_text(family = "Times New Roman", color = "black", size = 8),  
      axis.title = element_text(family = "Times New Roman", color = "black", size = 10),
      # axis.title.y = element_text(vjust = 1),
      legend.position = ""
    )
  # Save the plot as a JPG file
  ggsave("./fig/fig4_baseline_b.jpg", bar_statistic_areas, 
         width = 3, height = 2, dpi = 300)
  
  
  p = ggplot() +
    geom_raster(data = rastertemp_discrese, aes(x = lon, y = lat, fill = factor(class))) +
    scale_fill_manual(values = c("#b7b6b9","#4daf4a","#377eb8","#e41a1c")) +
    geom_sf(data = china_province, fill = NA, colour = "grey35", linewidth = 0.3) +
    geom_sf(data = nineline, colour = "black", linewidth = 0.25) +
    geom_sf(data = apple_planting, fill = NA, colour = "black", linewidth = 0.5) +  
    geom_sf(data = china, fill = NA, colour = "black", linewidth = 0.1) +
    geom_sf_text(data = china_province, aes(label = abb_name), size = 4) +
    # geom_sf_text(data = apple_planting, aes(label = planting), size = 8) +
    coord_sf(crs = albers, ylim = c(1500000,6000000), xlim = c(-3100000,2000000), expand = FALSE) +
    theme_test() + 
    labs(x = NULL, y = NULL,fill= "Suitable area") +
    theme(axis.text = element_text("Times New Roman", color = "black",size=base_size),
          axis.title = element_text("Times New Roman", color = "black",size=base_size),
          strip.text = element_text("Times New Roman", color = "black",size=base_size),
          legend.title = element_blank(),  
          legend.text = element_text(family = "Times New Roman", color = "black", size = base_size), 
          legend.position = c(0.52, 0.8),
          legend.direction = "vertical",
          # legend.spacing.y = unit(1, "lines"),
          legend.key.height = unit(1.5, "lines"),
          legend.key.width = unit(3, "lines"),
    ) +
    guides(fill = guide_legend(ncol = 2)) +
    annotation_scale(location = "bl",pad_y = unit(0.35, "cm"),text_cex = 1.2, text_family = family_bont) + 
    annotation_north_arrow(height = unit(2.5, "cm"),width = unit(2.5, "cm"),
                           location = "tl", style = north_arrow_nautical(
                             fill = c("grey35", "white"), line_col = "grey15"))
  # combined_plot <- ggdraw() +
  #   draw_plot(p) +
  #   draw_plot(bar_statistic_areas, x = 0.08, y = 0.08, width = 0.3, height = 0.25)
  
  # Save the plot as a JPG file
  ggsave("./fig/fig4_baseline_a.jpg", p, width = 10, height = 9, dpi = 300)
}




fig_climate_suitability_future_change_map <- function(){
  
  
  allfiels_paths <- fs::dir_ls(d_rsl, recurse = T,
                               glob = ".*gcms_ensemble_Change.tif")
  
  # filepaths2 = str_extract(allfiels_paths2, ".*(2041-2060245|2041-2060585|2081-2100245|2081-2100585).*\\.tif") %>% 
  #   na.omit()
  
  rastertemp = stack(allfiels_paths) %>% 
    projectRaster(crs = albers,method = "ngb")
  
  names(rastertemp) = str_extract(names(rastertemp), "\\d{4}.\\d{7}")
  
  rastertemp_discrese = as.data.frame(rastertemp, xy = T, na.rm = T) %>% 
    rename("lon" = x, "lat" = y)
  
  rastertemp_discrese = rastertemp_discrese %>%  
    pivot_longer(-c(lon, lat), names_to = "times",values_to = "values",values_drop_na = TRUE) %>%
    mutate(periods = str_extract(times,"(?<=X).{9}") %>% str_replace("\\.", "-"),
           ssps = paste0("SSP", str_extract(times, "\\d{3}$")),.before = 3) %>% 
    dplyr::select(c("lon","lat","periods","ssps","values")) 
  
  rastertemp_discrese <- rastertemp_discrese %>% 
    mutate(periods = case_when(
      periods == "2021-2040" ~ "2030s",
      periods == "2041-2060" ~ "2050s",
      periods == "2061-2080" ~ "2070s",
      periods == "2081-2100" ~ "2090s",
    ))
  
  
  p2 = ggplot() +
    geom_raster(data = rastertemp_discrese, aes(x = lon, y = lat, fill = factor(values)))+
    scale_fill_manual(values=c("#1a9641","#b7b6b6","#fdae61","#d7191c"),
                      labels = c("Contraction", "Still unsuitable", "Still suitable", "Expansions"))+
    geom_sf(data = china_province, fill = NA, colour = "grey35", linewidth = 0.1) +
    # geom_sf(data = nineline, colour = "black", linewidth = 0.25) +
    geom_sf(data = apple_planting, fill = NA, colour = "black", linewidth = 0.3) +  
    geom_sf(data = china, fill = NA, colour = "black", linewidth = 0.1) +
    facet_grid(ssps~periods)+
    coord_sf(crs = albers,
             ylim = c(1500000,6000000),
             xlim = c(-3100000,2000000),
             expand = FALSE)+
    theme_test() +
    labs(x = NULL,y=NULL)+
    theme(axis.text = element_text("Times New Roman", color = "black",size=16),
          axis.title = element_text("Times New Roman", color = "black",size=16),
          strip.text = element_text("Times New Roman", color = "black",size=16),
          legend.text = element_text(family = "Times New Roman", color = "black", size = 16), 
          legend.position = "bottom", 
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.key = element_rect(size = 4), # Increase the legend size
          legend.key.height = unit(2, "lines"),
          legend.key.width = unit(4, "lines"))
  
  # Save the plot as a JPG file
  ggsave("./fig/fig5_climate_changing_for_baseline.jpg", p2, width = 12, height = 12, dpi = 300)
}


fig_uncertain_analysis <- function(){
  data <- read_csv(paste0(d_rsl,"uncertain_analysis_results.csv")) %>% 
    mutate(times = case_when(
      times == "2021-2040" ~ "2030s",
      times == "2041-2060" ~ "2050s",
      times == "2061-2080" ~ "2070s",
      times == "2081-2100" ~ "2090s",
    )) %>% 
    mutate(uncertainty = 
             factor(uncertainty, 
                    levels = c("SDMs", "Scens", "GCMs", "Scens×SDMs", "Scens×GCMs",
                               "SDMs×GCMs", "Scens×SDMs×GCMs")))
  plot = data %>% 
    ggplot(aes(x = uncertainty, y = contribution, fill = times)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(),
             width = 0.7, size = 0.25) +
    scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
    theme_bw()+
    labs(x="", y="Contributions of uncertainties (%)",fill= "Periods")+
    theme(axis.text = element_text("Times New Roman", color = "black",size=7),
          axis.title = element_text("Times New Roman", color = "black",size=7),
          strip.text = element_text("Times New Roman", color = "black",size=7),
          axis.text.x = element_text(angle = 45, hjust = 1),# Setting font for facet titles
          legend.position = c(0.8,0.6),
          legend.title = element_blank(),
          legend.text = element_text("Times New Roman", color = "black",size=7),
          # legend.key = element_rect(size = 2.0), # Increase the legend size
          # legend.key.height = unit(1.5, "lines"),
          # legend.key.width = unit(1.5, "lines")
    ) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = c(seq(0, 100, by = 25)),
                       limits = c(0, 100)) 
  ggsave("./fig/fig6_uncertainty_analysis.jpg", plot, width = 2.75,
         height = 2.75, dpi = 300)
}

