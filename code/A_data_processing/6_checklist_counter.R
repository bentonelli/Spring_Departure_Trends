
total_checklists <- 0
for (each_year in 2002:2024){
  files_in <- list.files(paste0("~/Documents/Coding/R/Main/Avian_Pheno_Database/data/zf_files/",each_year),full.names = TRUE)
  
  for(each_file in files_in){
    fl_in <- readRDS(each_file)
    total_checklists <- total_checklists+nrow(fl_in)
  }
}

print(total_checklists)
#Total:
#90 million