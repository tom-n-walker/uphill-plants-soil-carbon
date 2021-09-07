################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load glasshouse species level data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_gh_plants <- function(){
  # load
  plants <- fread("./data/glasshouse_plants.csv") %>%
    as.data.frame %>%
    # select traits of interest
    select(Pot:PFT, rAGB.mg:rTB.mg, SLA.cm2.g, Amax, gsmax) %>%
    # drop NAs
    drop_na
  # return
  return(plants)
}
