################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Basic formatting of relevee data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

compute_cover <- function(raw_relevees){
  ## Filter for 2017 year only ----
  quad_2017 <- lapply(raw_relevees$quadrats, function(x) x %>% filter(year == 2017))
  rele_2017 <- lapply(raw_relevees$relevees, function(x) x %>% filter(year == 2017))
  ## Match quadrants and relevees to collar data ----
  # quadrants
  quad_matched <- lapply(
    quad_2017, 
    semi_join, 
    y = raw_relevees$collars, 
    by = "grid_id"
  )
  # relevees
  rele_matched <- lapply(
    rele_2017,
    semi_join,
    y = raw_relevees$collars, 
    by = "grid_id"
  )
  ## Quantify bare ground and species cover ----
  # bare ground
  quad_bare <- lapply(
    quad_matched, 
    function(x){
      x %>%
        filter(quadrant_fail != "M") %>%
        group_by(grid_id) %>%
        summarize(total_cover = sum(cover_cm2_bare_soil, na.rm = T)) %>%
        ungroup %>% 
        as.data.frame %>%
        mutate(species = "bare ground")
    }
  )
  # species covers
  rele_sum <- lapply(
    rele_matched, 
    function(x){
      x %>%
        group_by(grid_id, species) %>%
        summarize(total_cover = sum(cover_cm2, na.rm = T)) %>%
        ungroup %>% 
        as.data.frame
    }
  )
  ## Format output ----
  # unlist to data frames
  bare_long <- do.call(rbind, quad_bare)
  rele_long <- do.call(rbind, rele_sum)
  all_long <- rbind(bare_long, rele_long)
  # identify site information
  hilo_id <- substr(all_long$grid_id, 1, 3)
  calanda <- hilo_id %in% c("Cal", "Nes")
  # combine, add site_id, rearrange
  all_out <- all_long %>%
    mutate(site = if_else(calanda, "calanda", "lavey")) %>%
    arrange(grid_id, species) %>%
    remove_rownames
  # return
  return(all_out)
}
