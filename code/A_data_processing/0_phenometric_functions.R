### Functions and data for repeated code ###

season_time_windows <- data.frame(
  season = c("spring","spring","fall","fall"),
  arr_or_dep = c("arrival","departure","arrival","departure"),
  time_window_start = c(-30,-60,150,100),
  time_window_end = c(250,200,427,396)
)

#Simple function for setting up filters by season
check_filters <- function(arr_or_dep,season, detections_before_60,detections_after_150,
                          detections_before_200, detections_after_300,
                          num_detection_dates){
  if (arr_or_dep == "arrival" & season == "spring"){
    passed_filters <- detections_before_60 <= 0.15 & # Removes resident, mostly resident species
      detections_before_200 >= 0.10 & # Removes only fall migrants
      num_detection_dates >= 5
    return(passed_filters)
    
  } else if (arr_or_dep == "arrival" & season == "fall"){
    passed_filters <- detections_before_200 <= 0.10 & # Removes breeders, residents 
      detections_after_300 >= 0.25 &
      num_detection_dates >= 5
    
  } else if (arr_or_dep == "departure" & season == "spring"){
    passed_filters <- detections_before_60 >= 0.20 & 
      detections_after_150 <= 0.15 & 
      num_detection_dates >= 5
    
  } else if (arr_or_dep == "departure" & season == "fall"){
    passed_filters <- detections_after_300 <= 0.15 & 
      num_detection_dates >= 5
  } else {
    print("Failed")
  }
  return(passed_filters)
}
