################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Combine field data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

combine_field_data <- function(cover_data, trait_data, raw_soil){
  ## CWM traits ----
  # make species names compatible
  trait_data$accepted_name <- make.names(trait_data$accepted_name)
  # gather trait data
  wide_traits <- list(leaf_area = sub_col_spread(trait_data, "leaf_area"),
                      leaf_C = sub_col_spread(trait_data, "leaf_C"),
                      leaf_N = sub_col_spread(trait_data, "leaf_N"),
                      plant_height = sub_col_spread(trait_data, "plant_height"),
                      seed_mass = sub_col_spread(trait_data, "seed_mass"),
                      SLA = sub_col_spread(trait_data, "SLA"))  
  cover_data$traits <- list(wide_traits, wide_traits)
  # generate data frames of CWM traits (see small functions)
  allCover <- cover_data %>%
    mutate(all_cwm = map2(all_cover, traits, add_cwm_traits),
           focal_cwm = map2(focal_cover, traits, add_cwm_traits),
           bkgnd_cwm = map2(bkgnd_cover, traits, add_cwm_traits)) %>%
    select(-traits)
  ## Add soil data ----
  # nest soil data
  nestSoil <- raw_soil %>%
    as_tibble %>%
    arrange(site) %>%
    group_by(site) %>%
    nest
  # bind to cover data
  out <- left_join(allCover, nestSoil, by = "site") %>%
    mutate(soil_pools = map2(treatments, data, ~match_soil(.x, .y, match = "treatment"))) %>%
    select(-data)
  # return
  return(out)
}
