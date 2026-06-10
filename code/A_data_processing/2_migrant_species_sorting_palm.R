# Script to find candidate species for GAM model fitting
# Adding code to extract elevation bands as well here

library(readr)
library(dplyr)
library(dggridR)
library(data.table)
library(lubridate)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
year <- as.numeric(unlist(args[1]))
#year <- 2004

setwd("Avian_Pheno_Database/")

source("code/0_phenometric_functions.R")

#Read in last three years of data, 

f1_name <- paste("data/1_year_files/",year-1,"_all.txt",sep="")
f2_name <- paste("data/1_year_files/",year,"_all.txt",sep="")
f3_name <- paste("data/1_year_files/",year+1,"_all.txt",sep="")

data_in_1 <- fread(f1_name,sep = "|",header = FALSE,na.strings = "") %>%
  dplyr::select(c(1,2,3,4,5,6,7,9,12)) %>%
  filter(substr(V6,6,7)%in%c("11","12"))
gc()

data_in_3 <- fread(f3_name,sep = "|",header = FALSE,na.strings = "") %>%
  dplyr::select(c(1,2,3,4,5,6,7,9,12)) %>%
  filter(substr(V6,6,7)%in%c("01","02"))
gc()

data_in_2 <- fread(f2_name,sep = "|",header = FALSE,na.strings = "") %>%
  dplyr::select(c(1,2,3,4,5,6,7,9,12))
gc()

data_in <- rbindlist(list(data_in_1, data_in_2,data_in_3), 
                     use.names = TRUE, fill = TRUE) %>% 
  as.data.frame()

data_in_1 <- data_in_2 <- data_in_3<- NULL

nrow(data_in)

colnames(data_in) <- c("CATEGORY","COMMON.NAME","OBSERVATION.COUNT",
                       "LATITUDE","LONGITUDE",'OBSERVATION.DATE',
                       "SAMPLING.EVENT.IDENTIFIER",'DURATION.MINUTES',
                       "GROUP.IDENTIFIER")

#This speeds up processing
data_in$OBSERVATION.DATE <- as.character(data_in$OBSERVATION.DATE)

#Set up grids
hexgrid6 <- dggridR::dgconstruct(res = 6) 
data_in$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                         in_lon_deg = data_in$LONGITUDE, 
                                         in_lat_deg = data_in$LATITUDE)[[1]]

#Read in unique checklists
grid_cell_checklist <- unique(dplyr::select(data_in,c("SAMPLING.EVENT.IDENTIFIER","cell","DURATION.MINUTES")))

#Change effort from minutes to hours
grid_cell_checklist$DURATION.MINUTES <- grid_cell_checklist$DURATION.MINUTES/60
colnames(grid_cell_checklist)[c(1,3)] <- c("checklist_num","EH")

#Get unique sampling events by removing duplicates
sampling_events <- dplyr::select(data_in,c("SAMPLING.EVENT.IDENTIFIER","GROUP.IDENTIFIER")) %>% as.data.frame()
#distinct(sampling_events,GROUP.IDENTIFIER)
sampling_events <- sampling_events[!duplicated(sampling_events[,"GROUP.IDENTIFIER"], incomparables = NA),]
sampling_events <- unique(dplyr::select(sampling_events,"SAMPLING.EVENT.IDENTIFIER"))

#Add julian days, lat long
julian_days <- unique(dplyr::select(data_in,c("SAMPLING.EVENT.IDENTIFIER","OBSERVATION.DATE","LATITUDE","LONGITUDE")))
julian_days$julian_day <- yday(julian_days$OBSERVATION.DATE)
sampling_events <- merge(sampling_events,julian_days,by="SAMPLING.EVENT.IDENTIFIER")

#Check for leap years
leap_yr_t_minus_one <- leap_year(year - 1)

if(leap_yr_t_minus_one){
  days_in_year_tm <- 366
} else {
  days_in_year_tm <- 365
}

leap_yr_t <- leap_year(year)
if(leap_yr_t){
  days_in_year <- 366
} else {
  days_in_year <- 365
}


#Adjust julian dates, so that 0 is Dec 31st, previous year
sampling_events$julian_day[year(sampling_events$OBSERVATION.DATE) == (year-1)] <- 
  sampling_events$julian_day[year(sampling_events$OBSERVATION.DATE) == (year-1)] - days_in_year_tm

#Adjust julian dates, so that Jan. 1 of year t+1 is either 366 or 367, depending on leap year
sampling_events$julian_day[year(sampling_events$OBSERVATION.DATE) == (year+1)] <- 
  sampling_events$julian_day[year(sampling_events$OBSERVATION.DATE) == (year+1)] + days_in_year

#Filter to checklists within specified time window, rename columns
sampling_events <- sampling_events %>%
  filter(julian_day >= -60 & julian_day <= (days_in_year+61)) %>%
  dplyr::select(c("SAMPLING.EVENT.IDENTIFIER","julian_day","LATITUDE","LONGITUDE"))
colnames(sampling_events) <- c("checklist_num","julian_day","LATITUDE","LONGITUDE") 

#Merge together with effort
spec_matrix <- merge(sampling_events,grid_cell_checklist,by="checklist_num")

lat_longs <- dplyr::select(spec_matrix,c("LONGITUDE","LATITUDE"))
# Read in elevation file 
# (can be downloaded here: 
# https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation
# )
elev_medians <- terra::rast(x = "data/input/elevation/GMTED2010_30ARC_SEC_MEANS.tif")
spec_matrix$elev <- terra::extract(elev_medians,as.data.frame(lat_longs))[,2]

#Clear memory
elev_medians <- lat_longs <- julian_days <- grid_cell_checklist <- NULL
gc()

# There may be some NAs in the elevation raster in super remote areas
# In these cases, we want to set to sea level.
if (sum(is.na(spec_matrix$elev)>0)){
  #Print out if there are NAs
  print(paste("NA cells:",sum(is.na(spec_matrix$elev)),sep = " "))
  spec_matrix$elev[which(is.na(spec_matrix$elev))] <- 0
}

#Changed 7/9/25 - add elevational factors here
spec_matrix$elev_factor <- floor(spec_matrix$elev/500)
spec_matrix$elev_factor[which(spec_matrix$elev_factor < 0)] <- 0 # Below sea level
spec_matrix$elev_factor[which(spec_matrix$elev_factor > 6)] <- 6 # Above 3000 meters

#Save as factor
spec_matrix$elev_factor <- as.factor(spec_matrix$elev_factor)

#Remove Lat-long, no longer needed
spec_matrix <- dplyr::select(spec_matrix,c("checklist_num","julian_day","cell","EH","elev","elev_factor"))

#Filter out areas with less than X checklists & < X unique days
cell_sampling_information <- spec_matrix %>% 
  group_by(cell) %>%
  summarise(n_checklists = n(),
            n_unq_days = length(unique(julian_day))) %>%
  filter(n_unq_days >= 100)

well_sampled_cells <- cell_sampling_information %>% select(cell) %>%
  unlist() %>%
  as.vector()

print("Cells with over 100 checklists = ")
print(length(well_sampled_cells))

# At this point, filter the observation data to remove observation counts that equal zero (this is possible!)
# Also filter out everything except ISSF, Species, Form records
# Select only Common Name and Sampling event identifier
data_in <- data_in %>% 
  filter(CATEGORY %in% c("species","issf","form") & OBSERVATION.COUNT != 0) %>%
  select(COMMON.NAME,SAMPLING.EVENT.IDENTIFIER,cell)


number_spec_tested <- length(unique(data_in$COMMON.NAME))
print(number_spec_tested)
#For each of those cells, read in data from that cell, figure out which species are candidates

spec_report_card <- data.frame(year=double(),each_cell=double(),each_spec=character(),
                               season_in=character(),arr_or_dep_in=character(),
                               elev_band = double(),
                               accepted=double(),
                               time_window_start_in=double(),time_window_end_in=double(),
                               total_detections=double(),elev_band_detections = double(),
                               elev_mode_all_checklists = double(),
                               perc_checklists_reported=double(),
                               num_detection_dates=double(),
                               biggest_gap=double(),
                               detections_before_60=double(),detections_after_150=double(),
                               detections_before_200=double(),detections_after_300=double())

for (each_cell in well_sampled_cells){
  gc()
  #Cell-year list for species with at least one accepted migration state
  cl_sp_year_dt <- c()
  
  print(each_cell)
  #each_cell <- 390
  spec_matrix_cell <- spec_matrix %>% filter(cell==each_cell)
  
  cell_spec_data <- data_in %>% filter(cell == each_cell)
  
  #Only calculate for species seen at least 30 times (impossible to pass filters otherwise)
  
  well_represented_species <- cell_spec_data %>% 
    count(COMMON.NAME) %>% 
    filter(n >= 30) %>% 
    select(COMMON.NAME) %>%
    unlist() %>%
    as.vector()
  
  #Only run if more than 1 species candidate
  if(length(well_represented_species)>0){
    #For each species, add columns, and P/NP information for each checklist
    ws_spec_data <- cell_spec_data %>%
      filter(COMMON.NAME %in% well_represented_species) %>%
      distinct(SAMPLING.EVENT.IDENTIFIER,COMMON.NAME) %>%
      mutate(p=1)
    
    wide_p_dt <- ws_spec_data %>% 
      pivot_wider(names_from = COMMON.NAME,values_from = p,values_fill = 0)
    
    colnames(wide_p_dt)[1] <- "checklist_num"
    
    spec_matrix_cell <- merge(spec_matrix_cell,wide_p_dt,by="checklist_num", all.x = TRUE)
    
    #Remaining NAs are from checklists with no species of interest
    spec_matrix_cell[(is.na(spec_matrix_cell))] <- 0
    
    # Now take this zero-filled list and run it through checks for each of the 4 conditions
    # Do this using a dataframe with pre-filled time windows
    for (nn in 1:4){
      season_in <- season_time_windows$season[nn]
      arr_or_dep_in <- season_time_windows$arr_or_dep[nn]
      time_window_start_in <- season_time_windows$time_window_start[nn]
      time_window_end_in <- season_time_windows$time_window_end[nn]
      
      # Group checklists by 2-week periods, sample evenly with replacement from those
      # The idea here is to mitigate the effects of cells with uneven sampling effort
      # over the course of a year. 
      possible_breaks <- floor((30 + time_window_start_in:time_window_end_in)/14) %>% 
        as.data.frame()
      colnames(possible_breaks) <- "breaks"
      
      freq_breaks <- possible_breaks %>% group_by(breaks) %>% summarise(count=n())
      
      spec_matrix_cell_filt <- spec_matrix_cell %>% 
        filter(julian_day >= time_window_start_in & julian_day <= time_window_end_in)
      
      #Across each elevation band, get the total number of detections of each species
      # And the total number of checklist, and the number of unique days
      elev_diagnostics <- spec_matrix_cell_filt %>%
        select(elev_factor, julian_day, 7:ncol(spec_matrix_cell_filt)) %>%
        group_by(elev_factor) %>%
        summarise(
          across(2:(ncol(spec_matrix_cell_filt)-5), ~ sum(.x, na.rm = TRUE), .names = "{.col}"), # Summing multiple columns
          count = n(),
          unq_dates = n_distinct(julian_day),
          .groups = "drop"
        ) %>% filter(
          unq_dates >= 100
        )
      
      # Only continue if there are more than zero suitable elevation bands
      if (nrow(elev_diagnostics) > 0){
        set.seed(each_cell)
        stratified_sample <- spec_matrix_cell_filt %>%
          mutate(breaks = floor((30 + julian_day)/14)) %>%
          left_join(freq_breaks,by="breaks") %>%
          group_by(breaks,elev_factor) %>%
          sample_n(300,weight = count, replace=TRUE)
        
        #Get the largest gap in the data
        biggest_gap_eb <- spec_matrix_cell_filt %>% 
          group_by(elev_factor) %>% 
          summarise(gap_size = max(diff(sort(unlist(unique(julian_day))))))
        
        elev_mode_all_checklists <- elev_diagnostics$elev_factor[which.max(elev_diagnostics$count)]
        
        for(each_spec in well_represented_species){
          
          spec_col_index <- which(colnames(spec_matrix_cell_filt) == each_spec)
          total_detections <- sum(spec_matrix_cell_filt[,spec_col_index])
          
          min_elev_band_count <- elev_diagnostics[,c(1,
                                                     which(colnames(elev_diagnostics) == each_spec),
                                                     (ncol(elev_diagnostics)-1):ncol(elev_diagnostics))]
          if (nrow(min_elev_band_count) > 0){
            
            colnames(min_elev_band_count)[2] <- "detections"
            
            for(mm in 1:nrow(min_elev_band_count)){
              
              elev_to_sample <- min_elev_band_count$elev_factor[mm]
              
              elev_band_detections <- min_elev_band_count$detections[mm]
              
              perc_checklists_reported <- min_elev_band_count$detections[mm]/min_elev_band_count$count[mm]
              
              num_detection_dates <- min_elev_band_count$unq_dates[mm]
              
              biggest_gap <- biggest_gap_eb$gap_size[which(biggest_gap_eb$elev_factor==elev_to_sample)]
              
              strat_sample_elev_band <- stratified_sample %>% filter(elev_factor == elev_to_sample)
              spec_col_index_strat <- which(colnames(strat_sample_elev_band) == each_spec)
              total_detections_strat <- sum(strat_sample_elev_band[,spec_col_index_strat])
              
              jd_p_sp <- strat_sample_elev_band$julian_day[which(strat_sample_elev_band[,spec_col_index_strat]==1)]
              
              # Very rarely, the stratified sample with have zero detections,
              # In this case, accepted == 0 and all the detection < or > filters go to 1
              if(length(jd_p_sp) == 0){
                detections_before_60 <- 1
                detections_after_150 <- 1
                detections_before_200 <- 1
                detections_after_300 <- 1
              } else {
                detections_before_60 <- sum(jd_p_sp < 60)/total_detections_strat
                detections_after_150 <- sum(jd_p_sp > 150)/total_detections_strat
                detections_before_200 <- sum(jd_p_sp < 200)/total_detections_strat
                detections_after_300 <- sum(jd_p_sp > 300)/total_detections_strat
              }
              
              if(perc_checklists_reported >= 0.005 & 
                 elev_band_detections >= 30 & 
                 biggest_gap < 14 & 
                 length(jd_p_sp) > 0 &
                 check_filters(arr_or_dep_in, season_in, detections_before_60,
                               detections_after_150,detections_before_200,
                               detections_after_300,num_detection_dates)){
                
                #Add species name to list
                cl_sp_year_dt <- c(cl_sp_year_dt,each_spec)
                
                accepted <- 1
                
              } else {
                accepted <- 0
              } 
              spec_report_card <- add_row(spec_report_card,
                                          year,each_cell,
                                          each_spec,season_in,
                                          arr_or_dep_in,
                                          elev_band = as.numeric(as.character(min_elev_band_count$elev_factor[mm])),
                                          accepted,time_window_start_in,time_window_end_in,
                                          total_detections,
                                          elev_band_detections = elev_band_detections,
                                          elev_mode_all_checklists=as.numeric(elev_mode_all_checklists),
                                          perc_checklists_reported,num_detection_dates,biggest_gap,
                                          detections_before_60,detections_after_150,detections_before_200,detections_after_300)
              
            }
            
          } else {
            # No suitable elevation bands!
            
          }
        }
      }
    }
    
    #Filter spec_matrix_cell_filt to only species that are migratory (as defined by the filters above)
    mig_spec_index <- colnames(spec_matrix_cell) %in% unique(cl_sp_year_dt)
    
    #Only save zero-filled file if migratory species exist
    if(sum(mig_spec_index)>0){
      saveRDS(spec_matrix_cell[,c(1:6,which(mig_spec_index))],paste("data/input/phenology/zf_files/",year,"/zf_dt_",year,"_",each_cell,".rds",sep=""))
    }
  }
}
spec_report_card <- as.data.frame(spec_report_card)
colnames(spec_report_card) <- c("year",
                                "cell",
                                "spec",
                                "season",
                                "arr_or_dep",
                                "elev_band",
                                "accepted",
                                "time_window_start_in",
                                'time_window_end_in',
                                "total_detections",
                                "elev_band_detections",
                                "elev_mode_all_checklists",
                                "perc_checklists_reported",
                                "num_detection_dates",
                                "biggest_gap",
                                "detections_before_60",
                                "detections_after_150",
                                "detections_before_200",
                                "detections_after_300")


saveRDS(spec_report_card,paste("data/input/phenology/report_cards/rc_",year,".rds",sep=""))

#Get number of unique species with at least one accepted record
mig_species_list <- spec_report_card %>% dplyr::select(spec,accepted) %>% filter(accepted == 1) %>% unique()

accepted_list <- spec_report_card %>% filter(accepted==1)

print("GAMs to fit:")
accepted_list %>% select(spec,cell) %>% unique() %>% count() %>% as.numeric()
print("Metrics to estimate:")
nrow(accepted_list)
print("Unique cell-species-season metrics")
accepted_list %>% select(spec,cell,arr_or_dep,season)%>% unique() %>% count() %>% as.numeric()
print("Unique Species")
length(unique(accepted_list$spec))
print("Unique Locations")
length(unique(accepted_list$cell))
print("Metric breakdown:")
table(accepted_list$season,accepted_list$arr_or_dep)
