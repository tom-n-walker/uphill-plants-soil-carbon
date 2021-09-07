################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load flux data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_fluxes <- function(){
  # load fluxes
  calFlux <- fread("./data/field_calanda_gases.csv") %>%
    mutate(Treatment = substr(Treatment, nchar(Treatment), nchar(Treatment))) %>%
    mutate(site = "calanda") %>%
    # add leading zero to block
    mutate(Block = ifelse(
      nchar(Block) == 3,
      Block,
      paste0(substr(Block, 1, 1), "0", substr(Block, 2, 2))
    )) %>%
    select(site, Date, Block, Treatment, Temp.C:Soil.WVC, ER) %>%
    dplyr::rename_with(.cols = Date:Treatment, tolower)
  lavFlux <- fread("./data/field_lavey_gases.csv") %>%
    mutate(Treatment = substr(Treatment, nchar(Treatment), nchar(Treatment))) %>%
    mutate(site = "lavey") %>%
    # add leading zero to block
    mutate(Block = ifelse(
      nchar(Block) == 3,
      Block,
      paste0(substr(Block, 1, 1), "0", substr(Block, 2, 2))
    )) %>%
    select(site, Date, Block, Treatment, Temp.C:Soil.WVC, ER) %>%
    dplyr::rename_with(.cols = Date:Treatment, tolower)
  # collate & return
  fluxes <- bind_rows(calFlux, lavFlux) %>% as.data.frame
  return(fluxes)
}
