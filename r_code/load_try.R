################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load TRY data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_try <- function(){
  trySpecies <- fread("./data/try_impute/try_species.csv")
  try <- fread("./data/try_impute/try_nums_bhpmf_miss.csv")
  tryImputed <- apply_mice(try, 5)
  tryOut <- bind_cols(trySpecies, tryImputed) %>% as.data.frame
  return(tryOut)
}
