################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Load glasshouse pot data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

load_gh_soil <- function(){
  # load pot level data
  pots <- fread("./data/glasshouse_pot_data.csv") %>%
    as.data.frame %>%
    # drop mixed treatment
    filter(Treatment != "Mixed") %>%
    # rename treatments
    mutate(Treatment = substr(Treatment, 1, 1)) %>%
    rename_with(.cols = Pot:Block, tolower) %>%
    # create root to shoot ratio
    mutate(R2S = BGB.g/AGB.g)
  # load fluxes
  fluxes <- fread("./data/glasshouse_incubation_co2.csv") %>%
    as.data.frame %>%
    # remove mixed treatment
    filter(Treatment != "M") %>%
    # constrain to 6 weeks
    filter(Hours < (6 * 7 * 24 + 0.2)) %>%
    select(Day, Hours, Sample:R.ugC.g.h)
  # load soil model and add to pot data
  twoPool <- fread("./data/glasshouse_soil_pool_model.csv") %>%
    as.data.frame %>%
    filter(!is.na(p)) %>%
    # drop mixed treatment
    filter(treatment != "M") %>%
    # bind to pot data
    select(sample, p:b)
  # subset and collate data
  soil <- pots %>%
    select(pot:block, Cmic.ugC.g:baseMR.ugC.g.h, DOM.a350:DOM.C6) %>%
    left_join(., twoPool, by = c("pot" = "sample"))
  pots <- pots %>%
    select(pot:block, AGB.g:Tbio.g, R2S, cSLA.cm2.g:cGS.umol.m2.s) %>%
    filter(treatment != "B")
  # return
  out <- list(
    pot_plants = pots,
    pot_soil = soil,
    mic_resp = fluxes
  )
  return(out)
}
