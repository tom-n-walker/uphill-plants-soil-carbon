################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load soil data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_soil <- function(){
  # load pools
  pools <- fread("./data/field_soil_data.csv") %>%
    # edit treatment and column names, remove species data
    mutate(Treatment = substr(Treatment, nchar(Treatment), nchar(Treatment))) %>%
    mutate(Site = tolower(Site)) %>%
    # add leading zero to block
    mutate(Block = ifelse(
      nchar(Block) == 3,
      Block,
      paste0(substr(Block, 1, 1), "0", substr(Block, 2, 2))
    )) %>%
    # sort column names and select variables
    dplyr::rename_with(.cols = Sample:Treatment, tolower) %>%
    select(-sample) %>%
    as.data.frame
  # load lavey extra data and combine
  lav <- fread("./data/field_deep_soil_data.csv") %>%
    as.data.frame %>%
    select(Block:Elevation, MR.ugC.g.h:CUE) %>%
    rename_with(.cols = Block:Elevation, tolower) %>%
    mutate(site = "lavey") %>%
    # add leading zero to block
    mutate(block = ifelse(
      nchar(block) == 3,
      block,
      paste0(substr(block, 1, 1), "0", substr(block, 2, 2))
    )) %>%
    mutate(treatment = substr(treatment, nchar(treatment), nchar(treatment)))
  # combine
  out <- pools %>%
    left_join(., lav, by = c("site", "block", "treatment", "elevation")) %>%
    # setup proper units
    rename(Rm = MR.ugC.g.h, 
           Gm = MG.ugC.g.h, 
           RmM = MRm.ugCmic.g.h,
           GmM = MGm.ugCmic.g.h) %>%
    mutate(Cmic = Cmic/1000,
           CUE = CUE*100,
           RmM = RmM * 1000,
           GmM = GmM * 1000)
  # return
  return(out)
}
