################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Process TRY data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    27 May 2021
#### ---------------------------------------------------------------------------

process_try <- function(cover_data, full_try){
  ## Collate my species names ----
  # get all species
  sppList <- unique(do.call(c, lapply(cover_data$all_cover, colnames)))
  # format output
  allSpp <- sppList %>%
    # split string by period and make data frame
    str_split(., "\\.", simplify = T) %>%
    as.data.frame %>%
    # add original accepted name (changing period for space)
    mutate(accepted_name = str_replace_all(sppList, "\\.", " ")) %>%
    # rename columns
    select(accepted_name, genus = V1, species = V2, subspecies = V3) %>%
    # replace empty elements with NA
    replace(., . == "", NA)
  ## Subset and join TRY to my species ----
  # match to TRY 
  sppMatch <- !is.na(match(full_try$accepted_name, allSpp$accepted_name, incomparables = NA))
  genMatch <- !is.na(match(full_try$genus, allSpp$genus, incomparables = NA))
  # join to species and/or names
  sppTraits <- allSpp %>%
    left_join(., full_try[sppMatch, ], "accepted_name")
  genTraits <- allSpp %>%
    left_join(., full_try[sppMatch, ], c("genus" = "accepted_name"))
  # add genus-level traits where species traits absent
  sppMissing <- is.na(sppTraits$species)
  sppTraits[sppMissing, 12:19] <- genTraits[sppMissing, 12:19]
  ## Format selected trait data ----
  # impute, set up data frame and return
  imputed <- apply_mice(select(sppTraits, leaf_area:stem_density), 5)
  okTraits <- bind_cols(select(sppTraits, accepted_name), imputed)
  # identify as lowland or alpine
  lowlands <- map(cover_data$focal_cover, colnames) %>%
    do.call(c, .) %>%
    unique %>%
    str_replace_all(., "\\.", " ")
  out <- okTraits %>%
    mutate(is_focal = ifelse(accepted_name %in% lowlands, "yes", "no")) %>%
    select(accepted_name, is_focal, leaf_area:SLA)
  return(out)  
}
