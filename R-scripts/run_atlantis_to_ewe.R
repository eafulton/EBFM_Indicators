## Stand alone script to convert Atlantis output into the EwE style files used by
## EBFM_indicator_calculation_template.qmd
##
## Fill in the settings below and run (e.g. source("R-scripts/run_atlantis_to_ewe.R") from the
## EBFM_Indicators directory). Then set up the qmd as for an EwE model, using the directories, file
## names and area_model printed at the end of this script.
##
## See the notes at the top of R-scripts/atlantis_to_ewe_funcs.R for how each EwE quantity is
## derived from the Atlantis output.

## Clean up the workspace
rm(list=ls())

source("R-scripts/atlantis_to_ewe_funcs.R")

########################### SETTINGS ###########################

### The Atlantis model files ----
atlantis_model_dir <- "/Users/ful083/runs/run_area/SEAP_runs/SEAP_wind/SEAP_Wind_Runs"
groups_file <- file.path(atlantis_model_dir, "NESP_groups.csv")
bgm_file <- file.path(atlantis_model_dir, "SEAP_NESP_final_xy_fixed.bgm")
bio_prm_file <- file.path(atlantis_model_dir, "NESP_bio.prm")
fisheries_file <- file.path(atlantis_model_dir, "CEPFisheries_FishedCEP.csv")  # Set to NULL if no fisheries
init_nc_file <- file.path(atlantis_model_dir, "NESP_init.nc")   # Used to get Linf for the Species_ID template

### The Atlantis output of the scenario of interest (e.g. historical fishing) ----
scenario_out_dir <- file.path(atlantis_model_dir, "NESPout")
scenario_prefix <- "NESPout"      # The output file name given to Atlantis with -o (without the .nc)

### The Atlantis output of the unfished run (set to NULL if using a burn-in period in the scenario run instead) ----
nofish_out_dir <- NULL
nofish_prefix <- "NESPout"

### Where to put the EwE style files - use the same values as in the qmd set.file.paths chunk ----
mainDir <- "/Users/ful083/Work/EBFM_Atlantis/SEAP_wind"   # Where the qmd, Species_ID.csv, Fleet_ID.csv and matrices are kept
InputDir <- file.path(mainDir, "Atlantis_as_EwE_output")  # Equivalent of the EwE output directory
InputsubDir <- "Wind_scenario"                            # Sub directory for the scenario of interest
BosubDir <- "NoFishing"                                   # Note the qmd reads the no fishing files from file.path(mainDir, BosubDir)

### File names - use the same values as in the qmd file.names chunk ----
species_ID_filename <- "Species_ID.csv"
fleet_ID_filename <- "Fleet_ID.csv"
consumption_filename <- "Atlantis_Consumption.csv"
diet_filename <- "Atlantis_Diet_composition.csv"

### Conversion parameters ----
model_name <- "SEAP NESP"
start_year <- 2000       # Calendar year of the first year of the Atlantis run
area_km2 <- NULL         # Model area in km2 - NULL means use the area of the non-boundary boxes in the bgm file
ecopath_year_row <- 1    # Year (row) of the scenario run used for the Ecopath style consumption and diet matrices
pisc_threshold <- 0.5    # Proportion of fish (and squid) in the diet above which a fish is classed Piscivorous
# Optionally assign groups as Piscivorous or Planktivorous yourself (by Atlantis code or name) - these override
# the diet based classification above, e.g. pisczoopl_groups <- c(FTB = "Piscivorous", Sardine = "Planktivorous")
pisczoopl_groups <- NULL
include_m0 <- TRUE       # Include non-predation natural mortality (from SpecificMort.txt) in the total mortality
overwrite_id_files <- FALSE  # If FALSE existing Species_ID / Fleet_ID files are kept and templates written alongside
catch_source <- "auto"   # Fisheries output to use: "auto", "txt" (CatchPerFishery/DiscardPerFishery.txt), "nc" (CATCH.nc) or "total" (Catch.txt/Discard.txt as one fleet)

########################### CONVERSION ###########################

### Scenario of interest ----
atl_scenario <- convert_atlantis_to_ewe(
  atlantis_out_dir = scenario_out_dir, output_prefix = scenario_prefix,
  groups_file = groups_file, bgm_file = bgm_file, bio_prm_file = bio_prm_file,
  ewe_out_dir = file.path(InputDir, InputsubDir),
  fisheries_file = fisheries_file, init_nc_file = init_nc_file,
  network_dir = mainDir, model_name = model_name, scenario_name = InputsubDir,
  start_year = start_year, area_km2 = area_km2, ecopath_year_row = ecopath_year_row,
  consumption_filename = consumption_filename, diet_filename = diet_filename,
  species_ID_filename = species_ID_filename, fleet_ID_filename = fleet_ID_filename,
  overwrite_id_files = overwrite_id_files, pisc_threshold = pisc_threshold, pisczoopl_groups = pisczoopl_groups,
  include_m0 = include_m0, catch_source = catch_source)

### Unfished run ----
if (!is.null(nofish_out_dir)) {
  atl_nofish <- convert_atlantis_to_ewe(
    atlantis_out_dir = nofish_out_dir, output_prefix = nofish_prefix,
    groups_file = groups_file, bgm_file = bgm_file, bio_prm_file = bio_prm_file,
    ewe_out_dir = file.path(mainDir, BosubDir),
    fisheries_file = fisheries_file, init_nc_file = init_nc_file,
    network_dir = NULL, model_name = model_name, scenario_name = BosubDir,
    start_year = start_year, area_km2 = atl_scenario$area_model, include_m0 = include_m0,
    catch_source = catch_source)
}

########################### SETTINGS FOR THE QMD ###########################

cat("\nSettings to use in EBFM_indicator_calculation_template.qmd:\n")
cat(paste0('mainDir <- "', mainDir, '"\n'))
cat(paste0('InputDir <- "', InputDir, '"\n'))
cat(paste0('InputsubDir <- "', InputsubDir, '"\n'))
cat(paste0('BosubDir <- "', BosubDir, '"\n'))
cat(paste0('species_ID_filename <- "', basename(atl_scenario$files$species_ID), '"\n'))
cat(paste0('fleet_ID_filename <- "', basename(atl_scenario$files$fleet_ID), '"\n'))
cat(paste0('consumption_filename <- "', consumption_filename, '"\n'))
cat(paste0('diet_filename <- "', diet_filename, '"\n'))
cat(paste0("bo_in_same_file <- ", ifelse(is.null(nofish_out_dir), 1, 0), "\n"))
cat(paste0("area_model <- ", round(atl_scenario$area_model, 3), "\n"))
cat(paste0("# Years written: ", min(atl_scenario$years), " - ", max(atl_scenario$years),
           " (", atl_scenario$nyears, " rows) - choose bo_yr_row and yr_ignore accordingly\n"))
