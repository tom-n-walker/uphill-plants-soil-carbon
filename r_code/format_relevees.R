################################################################################
#### Project: Lowland plant migrations alpine soil C loss
#### Title:   Function | Compile plant community data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 May 2021
#### ---------------------------------------------------------------------------

format_relevees <- function(raw_relevees){
  formatted <- compute_cover(raw_relevees)
  named <- clean_species_names(formatted)
  summarised <- summarise_cover(raw_relevees, named)
}
