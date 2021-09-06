################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Clean species names
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

clean_species_names <- function(formatted){
  # species names
  species <- unique(formatted$species)
  # call name resolver
  resolved <- gnr_resolve(species, 
                          best_match_only = T, 
                          with_canonical_ranks = T)
  # take resolved, select, join to data, merge gnr/original names, select
  cleaned <- resolved %>%
    select(user_supplied_name, accepted_name = matched_name2) %>%
    left_join(formatted, ., by = c("species" = "user_supplied_name")) %>%
    select(site, grid_id, species, accepted_name, total_cover)
  # deal with fuckers
  cleaned[cleaned$species == "bare ground", "accepted_name"] <- "bare ground"
  cleaned[cleaned$species == "Erigeron-Aster Group", "accepted_name"] <- "unknown"
  cleaned[cleaned$species == "Aster-Erigeron Group", "accepted_name"] <- "unknown"
  cleaned[cleaned$species == "Sebastian has photo_UKS", "accepted_name"] <- "unknown"
  cleaned[cleaned$species == "Grass Group", "accepted_name"] <- "unknown"
  cleaned[cleaned$species == "Lychen Group", "accepted_name"] <- "mosses"
  cleaned[cleaned$species == "Moos Group", "accepted_name"] <- "mosses"
  cleaned[cleaned$species == "CAREXLACCA", "accepted_name"] <- "Carex flacca"
  cleaned[cleaned$species == "Baby carlina_UKS", "accepted_name"] <- "Carlina"
  # return species table
  cleaned <- cleaned %>% select(-species)
  return(cleaned)
}
