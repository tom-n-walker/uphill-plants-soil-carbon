################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Small functions
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

# calculate CWMs for trait data
add_cwm_traits <- function(.x, .y){
  # apply a weighted average to each  bin
  cwm_list <- lapply(.y, function(y){
    match_trait <- match(colnames(.x), colnames(y))
    trait_value <- y[, match_trait] %>% as.matrix %>% as.vector
    cwm_trait <- .x %>%
      apply(1, function(.x){weighted.mean(w = .x, x = trait_value)})
    return(cwm_trait)
  })
  # bind into data frame
  cwm_bound <- do.call(cbind, cwm_list) %>% as.data.frame %>% as_tibble
  # return
  return(cwm_bound)
}

# calculate biomass from cover data
biomass <- function(x, y){
  # index for vegetation types
  moss_index <- colnames(x) %in% "mosses"
  bare_index <- colnames(x) %in% "bare.ground"
  vasc_index <- !(colnames(x) %in% c("mosses", "bare.ground"))
  # select vegetation types
  bio_out <- data.frame(moss_bio = x[, moss_index] %>% rowSums,
                        vasc_bio = x[, vasc_index] %>% rowSums,
                        bare_bio = x[, bare_index] %>% rowSums,
                        vege_bio = x[, !bare_index] %>% rowSums,
                        focal_bio = rowSums(y)) %>%
    as_tibble %>%
    mutate(bkgnd_bio = vege_bio - focal_bio)
  # return
  return(bio_out)
}

# format treatment information from collar data to grid-id
join_treats <- function(.x, .y){
  out <- left_join(.x, .y) %>%
    # format treatments correctly
    mutate(elevation = substr(marc_treatment, 1, 1),
           block = paste0(elevation, num_in_string(marc_treatment)),
           treatment = substr(marc_treatment, nchar(marc_treatment), nchar(marc_treatment))) %>%
    # add leading zero to block
    mutate(block = ifelse(
      nchar(block) == 3,
      block,
      paste0(substr(block, 1, 1), "0", substr(block, 2, 2))
    )) %>%
    # select columns
    select(grid_id, elevation:treatment)
  return(out)
}

# match soil data rows to plant data treatments
match_soil <- function(.x, .y, match){
  # matching variables
  match_by <- c("elevation","block", "treatment")
  names(match_by) <- c("elevation", "block", match)
  # join data and return
  out <- left_join(.x, .y, by = match_by) %>%
    select(Soil.temp:CUE)
  return(out)
}

# get number from string
num_in_string <- function(x){
  out <- str_extract_all(x, "[:digit:]", simplify = T) %>%
    apply(1, paste0, collapse = "")
  return(out)
}

#  take DF, calculate relative abundance, transpose, make DF
ra_t_df <- function(x){
  out <- x %>%
    apply(1, function(x) x/sum(x, na.rm = T)) %>%
    t
  out[is.na(out)] <- 0
  out <- as_tibble(out)
  return(out)
}

# subset focal data set by site-presences
select_focals <- function(focals, site){
  tf <- focals[, site] == "Y"
  out <- focals %>% filter(tf) %>%
    transmute(accepted_name = make.names(accepted_name))
  return(out)
}


# subset, make column names, spread ----
sub_col_spread <- function(traits, subset){
  # select subset
  now <- traits[, c("accepted_name", subset)]
  colnames(now) <- c("accepted_name", "value")
  out <- now %>%
    pivot_wider(names_from = accepted_name, values_from = value)
  return(out)
}

# subset cover data for background plants only
subset_bckgnd <- function(.x, .y){
  match_focals <- colnames(.x) %in% .y$accepted_name
  out <- .x[, !match_focals]
  return(out)
}

# subset cover data for focals only 
subset_focals <- function(.x, .y){
  match_focals <- colnames(.x) %in% .y$accepted_name
  out <- .x[, match_focals]
  return(out)
}
