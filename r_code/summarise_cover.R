################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Summarise cover at plot level
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

summarise_cover <- function(raw_relevees, named){
  # summarise cover by site, grid_id, species
  sum_cover <- named %>%
    group_by(site, grid_id, accepted_name) %>%
    summarize(total_cover = sum(total_cover)) %>%
    mutate(accepted_name = make.names(accepted_name)) %>%
    ungroup
  # group by site, nest, map for cover, treatments
  allCover <- sum_cover %>%
    group_by(site) %>%
    nest %>%
    # pivot data wider and replace missing with zeros
    mutate(all_cover = map(data, pivot_wider, "grid_id", "accepted_name", values_from = "total_cover")) %>%
    mutate(all_cover = map(all_cover, function(x) replace(x, is.na(x), 0))) %>%
    mutate(treatments = map(all_cover, ~join_treats(., raw_relevees$collars))) %>%
    mutate(all_cover = map(all_cover, ~select(., -grid_id)))
  # map for focal ID
  allCover$focals <- list(
    select_focals(raw_relevees$focals, "calanda"),
    select_focals(raw_relevees$focals, "lavey")
  )
  # map for focal/bkgnd cover, rel_abund all covers, biomass
  allCover <- allCover %>%
    # get focal cover and summarise at group level
    mutate(focal_cover = map2(all_cover, focals, subset_focals)) %>%
    mutate(group_cover = map2(all_cover, focal_cover, biomass)) %>%
    # remove unknown, bare ground and moss from background community
    mutate(all_cover = map(all_cover, function(x) select(x, -bare.ground, -mosses, -unknown))) %>%
    # subset for background cover
    mutate(bkgnd_cover = map2(all_cover, focals, subset_bckgnd)) %>%
    # recalculate as relative abundances
    mutate(focal_ra = map(focal_cover, ra_t_df)) %>%
    mutate(bkgnd_ra = map(bkgnd_cover, ra_t_df)) %>%
    mutate(all_ra = map(all_cover, ra_t_df))
  # subset for columns of interest
  out <- allCover %>%
    select(site, treatments, 
           all_cover, focal_cover, bkgnd_cover, group_cover,
           all_ra, focal_ra, bkgnd_ra)
  # return
  return(out)
}
