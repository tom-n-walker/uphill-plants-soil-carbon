################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Make response data for exploring soil C loss
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    31 May 2021
#### ---------------------------------------------------------------------------

subset_field_plants <- function(field_data){
  ## Summarise soil data ----
  # unnest and delete high elevation control
  factorSoil <- field_data %>%
    # select nested variables
    select(site, treatments, soil_pools) %>%
    # unnest and ungroup
    unnest(cols = everything()) %>%
    ungroup %>%
    # select columns
    select(site, elevation:treatment, Csoil) %>%
    # make data frame
    as.data.frame
  # Calculate soil C loss by block
  soilCloss <- factorSoil %>%
    # remove high site control (not important here) and subset columns
    filter(treatment != "C") %>%
    # arrange by block and group
    arrange(site, block, treatment) %>%
    # group by site and block, calculate response ratio
    group_by(site, block) %>%
    summarise(soilCloss = Csoil[1]/Csoil[2]) %>%
    select(site, block, soilCloss) %>%
    as.data.frame
  
  ## Summarise basic PCC data ----
  # basic formatting
  basicPCC <- field_data %>%
    # do nmds on background relative abundance
    mutate(bkgnd_nmds = map(
      bkgnd_ra, 
      function(x){
        nmds <- metaMDS(x, k = 2, trymax = 1000)
        out <- as.data.frame(nmds$points)
        colnames(out) <- c("bgPCC1", "bgPCC2")
        return(out)
      }
    )) %>%
    # do PCA on CWM traits
    mutate(cwm_pca = map(
      all_cwm,
      function(x){
        # transform skewed variables
        xReady <- x %>%
          mutate(across(c(leaf_area, plant_height, seed_mass), log10))
        pca <- prcomp(x, scale = T, center = T)
        out <- as.data.frame(pca$x[, 1:3])
        colnames(out) <- c("cwmPC1", "cwmPC2", "cwmPC3")
        return(out)
      }
    )) %>%
    # select columns
    select(site, treatments, group_cover, bkgnd_nmds, focal_cover, all_cwm, cwm_pca)
  # Make long
  factorPCC <- basicPCC %>%
    # select nested variables
    select(site:group_cover, bkgnd_nmds, all_cwm, cwm_pca) %>%
    # unnest and ungroup
    unnest(cols = everything()) %>%
    ungroup %>%
    # add soil C and blocking effect
    mutate(Csoil = factorSoil$Csoil) %>%
    mutate(site_block = paste0(substr(site, 1, 1), block)) %>%
    # select columns and make data frame
    select(site, block, site_block, treatment, Csoil, moss_bio:cwmPC3) %>%
    as.data.frame
  # create focal biomass for adding to response data frame (below)
  focalBio <- factorPCC %>% 
    filter(treatment == "I") %>% 
    arrange(site, block) %>% 
    .$focal_bio
  # calculate response ratios
  pccResponse <- factorPCC %>%
    # remove high site control (not important here) and subset columns
    filter(treatment != "C") %>%
    # arrange by block and group
    arrange(site, site_block, block, treatment) %>%
    # group by site and block, calculate response ratio
    group_by(site, site_block, block) %>%
    summarise(across(c(-treatment, -focal_bio), function(x){x[1]/x[2]})) %>%
    # ungroup and add focal biomass
    ungroup %>%
    mutate(focal_bio = focalBio) %>%
    select(site:bkgnd_bio, focal_bio, bgPCC1:cwmPC3) %>%
    as.data.frame
  
  ## Build focal data ----
  # basic focal cover processing
  focals <- basicPCC %>%
    # select nested variables
    select(site, treatments, focal_cover) %>%
    # filter for WI treatment
    mutate(focals_wide = map2(
      treatments, focal_cover, 
      ~filter(.y, .x$treatment == "I") %>% as.data.frame
    )) %>%
    # add soil C data (simple for wide data)
    bind_cols(., soilCloss %>% group_by(site) %>% nest) %>%
    mutate(all_wide = map2(
      focals_wide, data, 
      ~bind_cols(.y, .x) %>%
        as.data.frame
    )) %>%
    # pivot longer to make long data
    mutate(all_long = map(
      all_wide, 
      ~pivot_longer(
        ., 
        c(-block, -soilCloss), 
        names_to = "species", 
        values_to = "cover") %>%
        as.data.frame
    )) %>%
    # subset columns
    select(site = `site...1`, all_wide, all_long)
  
  ## Output ----
  out <- list(
    plantsFull = factorPCC,
    plantsRR = pccResponse,
    focals = focals
  )
  return(out)
}
