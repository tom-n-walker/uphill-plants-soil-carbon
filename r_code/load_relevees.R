################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load relevee data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_relevees <- function(){
  # specify directory
  dir_quad <- paste0("./data/relevee_received/quadrant_info/")
  dir_rele <- paste0("./data/relevee_received/relevees/")
  # load quadrat information
  quad <- list(
    cal = fread(paste0(dir_quad, "quad_cal.csv"), data.table = F),
    n01 = fread(paste0(dir_quad, "quad_nes_01.csv"), data.table = F),
    n03 = fread(paste0(dir_quad, "quad_nes_03.csv"), data.table = F),
    pra = fread(paste0(dir_quad, "quad_pra.csv"), data.table = F),
    rio = fread(paste0(dir_quad, "quad_rionda.csv"), data.table = F)
  )
  # load relevee data
  rele <- list(
    cal = fread(paste0(dir_rele, "rel_calanda.csv"), data.table = F),
    n01 = fread(paste0(dir_rele, "rel_nesel_01.csv"), data.table = F),
    n03 = fread(paste0(dir_rele, "rel_nesel_03.csv"), data.table = F),
    pra = fread(paste0(dir_rele, "rel_pra.csv"), data.table = F),
    rio = fread(paste0(dir_rele, "rel_rionda.csv"), data.table = F)
  )
  # load collar information
  collars <- fread("./data/collars.csv", data.table = F) %>%
    select(grid_id, grid, marc_treatment, tom_treatment = treatment)
  # load focal information
  focals <- fread(paste0(getwd(), "/data/focal_names.csv"), data.table = F)
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
