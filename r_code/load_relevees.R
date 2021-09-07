################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load relevee data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_relevees <- function(){
  # load quadrat information
  quad <- list(
    cal = fread("./data/relevee_received/quad_cal.csv", data.table = F),
    n01 = fread("./data/relevee_received/quad_nes_01.csv", data.table = F),
    n03 = fread("./data/relevee_received/quad_nes_03.csv", data.table = F),
    pra = fread("./data/relevee_received/quad_pra.csv", data.table = F),
    rio = fread("./data/relevee_received/quad_rionda.csv", data.table = F)
  )
  # load relevee data
  rele <- list(
    cal = fread("./data/relevee_received/rel_calanda.csv", data.table = F),
    n01 = fread("./data/relevee_received/rel_nesel_01.csv", data.table = F),
    n03 = fread("./data/relevee_received/rel_nesel_03.csv", data.table = F),
    pra = fread("./data/relevee_received/rel_pra.csv", data.table = F),
    rio = fread("./data/relevee_received/rel_rionda.csv", data.table = F)
  )
  # load collar information
  collars <- fread("./data/field_collars.csv", data.table = F) %>%
    select(grid_id, grid, marc_treatment)
  # load focal information
  focals <- fread(paste0(getwd(), "/data/field_focal_names.csv"), data.table = F)
  # collate data into list
  out <- list(
    quadrats = quad,
    relevees = rele,
    collars = collars,
    focals = focals
  )
  # return
  return(out)
}
