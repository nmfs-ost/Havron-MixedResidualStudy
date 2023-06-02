library(httr)
library(here)

## Update results from github site
update_results <- function(update = FALSE){

  ## Define path to results files
  path <- paste0(here(), "/results/")

  if(update){
    req <- GET("https://api.github.com/repos/Cole-Monnahan-NOAA/mixed_resids/git/trees/main?recursive=1")
    stop_for_status(req)
    filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)
    res.list <- grep("results/", filelist, value = TRUE, fixed = TRUE)

    githubURL <- "https://github.com/Cole-Monnahan-NOAA/mixed_resids/tree/main/results"
    for(i in seq_along(res.list)){
      file <-  gsub("results/", "", res.list[[i]], perl = TRUE)
      download.file(githubURL, paste0(path, file))
    }
  }
}
