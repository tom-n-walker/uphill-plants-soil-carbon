################################################################################
#### Project: SNF field experiment
#### Title:   Drake plan
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 March 2021
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())

# configure default R session options (no factors, bias against scientific #s)
options(stringsAsFactors = F,
        scipen = 6)

## Libraries ----
source("packages.r")

## Code ----
sapply(list.files("./r_code",full.names = T), source)


#### PLANS ---------------------------------------------------------------------

## Load data ----
load_plan <- drake_plan(
  raw_relevees = load_relevees(),
  raw_soil = load_soil(),
  flux_data = load_fluxes(),
  full_try = load_try(),
  gh_plants = load_gh_plants(),
  gh_soil = load_gh_soil()
)

## Format field data ----
field_plan <- drake_plan(
  cover_data = format_relevees(raw_relevees = raw_relevees),
  trait_data = process_try(cover_data = cover_data, full_try = full_try),
  field_data = combine_field_data(cover_data = cover_data, trait_data = trait_data, raw_soil = raw_soil),
  field_data_subset = subset_field_plants(field_data = field_data)
)

# ## Write all data to file ----
# export_plan <- drake_plan(
#   localExport = target(
#     command = {
#       save(
#         # Original species information
#         edAllSpp,
#         # Species information
#         pfSpecies, edSpecies,
#         # Phylogenies
#         pfPhylogeny, edPhylogeny,
#         # Trait data
#         pfTraits, edTraits,
#         pfTraitsPCA, edTraitsPCA,
#         # Metabolite data
#         pfMtbs, edMtbs,
#         pfMtbsDiv, edMtbsDiv,
#         pfMtbsPCoA, edMtbsPCoA,
#         # Biogeography data
#         pfBioGeo, edBioGeo,
#         pfClimPCA, edClimPCA,
#         # Output
#         file = file_out("./data/exported/all_data.RData")
#       )
#       save(
#         # Full biogeographic data for plotting
#         pfBioGeoAll, edBioGeoAll,
#         # Output
#         file = file_out("./data/exported/geo_data.RData")
#       )
#       
#     }
#   )
# )


#### MAKE ----------------------------------------------------------------------

## Collate plans ----
weAreGo <- bind_rows(
  load_plan,
  field_plan
)

## Make ----
make(weAreGo)
