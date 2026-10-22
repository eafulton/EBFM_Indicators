## Functions to convert Atlantis output into the EwE (Ecopath with Ecosim) style csv
## files read by EBFM_indicator_calculation_template.qmd (and the routines in funcs.R)
##
## Files created (in the formats exported by EwE 6.7):
##   Ecosim style (one directory per scenario):
##     biomass_annual.csv              t/km2 per group per year
##     catch_annual.csv                total catch (landings + discards) t/km2/yr per group
##     consumption-biomass_annual.csv  Q/B per group per year (yr-1)
##     mortality_annual.csv            total mortality Z (= M0 + M2 + F) per group per year (yr-1)
##     tl_annual.csv                   trophic level per group per year (calculated from the realised diet)
##     landings_annual.csv             landings per fleet per group (long format) t/km2/yr
##     catch-fleet-group_annual.csv    catch (landings + discards) per fleet per group (long format) t/km2/yr
##     mort-fleet-group_annual.csv     fishing mortality per fleet per group (long format) yr-1
##   Ecopath style (for the network/hub species analysis):
##     <consumption_filename>          consumption matrix (prey rows x predator columns) t/km2/yr
##     <diet_filename>                 diet composition matrix (prey rows x predator columns) proportions
##   Helper files (templates - CHECK AND EDIT THEM):
##     <species_ID_filename>           Species_ID.csv as described in "Calculating Hub index and ETI from Ecopath with Ecosim.docx"
##     <fleet_ID_filename>             Fleet_ID.csv
##
## How the EwE quantities are derived from Atlantis:
##   Biomass     - <prefix>BiomIndx.txt (tonnes, whole model domain excluding boundary boxes). The annual
##                 value is the mean of the snapshots from the start to the end of each model year.
##   Consumption - <prefix>PROD.nc "<Name><cohort>_Eat" (age structured groups) and "<Name>Grazing" (biomass
##                 pools) treated as mg N m-3 d-1 and multiplied by box volume (the same convention as
##                 atlantistools::calculate_consumed_biomass). Boundary boxes are excluded. The annual value is
##                 the mean daily rate across the PROD.nc records in that year * 365.
##                 This total is split across prey using the diet proportions in <prefix>DietCheck.txt
##                 (per predator cohort), giving an annual prey x predator consumption matrix.
##   Q/B         - annual consumption of the predator / its annual biomass
##   Z           - M0 + M2 + F, where
##                   M2 = annual consumption of the group by all predators / its biomass (from the matrix above)
##                   F  = annual catch (landings + discards) / biomass
##                   M0 = non-predation natural mortality from the M1 columns of <prefix>SpecificMort.txt,
##                        averaged across cohorts and converted to an instantaneous rate -log(1 - M1)
##   TL          - solved each year from the realised diet (TL = 1 + sum(diet * TL prey)), with primary
##                 producers and detritus set to exactly 1 (funcs.R identifies producers by TL == 1)
##   Catch       - <prefix>CatchPerFishery.txt and <prefix>DiscardPerFishery.txt (tonnes per fleet per year),
##                 or if they are not present <prefix>CATCH.nc ("<Code>_Catch_FC<n>", "<Code>_Discard_FC<n>"),
##                 or <prefix>Catch.txt / <prefix>Discard.txt (one aggregate fleet). If there is no fisheries
##                 output (e.g. a no fishing run) a single dummy fleet with zero catch is written.
## All the amounts are divided by the model area (km2) so they are per km2 like EwE output. The qmd multiplies
## them back up by area_model, so set area_model in the qmd to the area returned by convert_atlantis_to_ewe().
##
## Atlantis outputs required: BiomIndx.txt, DietCheck.txt, SpecificMort.txt, PROD.nc and the main output .nc
## (for box volumes). Fishery outputs are optional (see above).
##
## Stanzas (convert_atlantis_to_ewe(..., stanzas = TRUE)):
##   By default each Atlantis group is a single EwE group. With stanzas = TRUE the age structured groups (those with
##   "<Name><cohort>_Nums" in the main output .nc) are instead split into stanzas, like EwE multi-stanza groups. Each
##   stanza is written as its own group ("<Code>_<label>", "<Name>_<label>") in place of the Atlantis group.
##   Stanza split - by default the Atlantis juvenile / adult split (cohorts before <Code>_age_mat are "juv", the
##                  rest "adult"). stanza_breaks can give other splits per group (the first cohort of each stanza
##                  after the first, e.g. list(FTS = c(3, 6)) gives cohorts 1-2, 3-5 and 6+) and stanza_labels the
##                  labels to use (e.g. list(FTS = c("larval", "juv", "adult"))). stanza_groups limits which groups are split.
##   Biomass     - Nums * (StructN + ResN) * X_CN * k_wetdry / 1e9 per cohort from the main output .nc (boundary
##                 boxes excluded, annual mean of the snapshots as for BiomIndx.txt). These cohort biomasses are
##                 used to share the BiomIndx.txt biomass of the group across its stanzas.
##   Predators   - each cohort's consumption (PROD.nc) and diet (DietCheck.txt) is given to the stanza it is in.
##   Prey        - DietCheck.txt only gives the prey group, so the consumption of a split prey group by a predator
##                 cohort is shared across the prey cohorts following Avail_Fish() in atvertprocesses.c - in
##                 proportion to pPREY<prey stage><Pred><pred stage> (from the biology prm file) x the size based
##                 feeding window (KLP_, KUP_, Kmax_coefft_, li_a_, li_b_ with UseHardFeedingWindow and
##                 UseBiLogisticFeedingWindow) x the prey cohort biomass. This is done box by box, weighted by where the
##                 predator cohort fed (its Eat per box). Sizes are the annual mean individual StructN and ResN of
##                 each cohort. The refuge, stock, turbidity and pH scalars in Avail_Fish() are not included.
##   Catch       - each group's catch per fleet is shared across its stanzas using the catch and discards per cohort
##                 in CATCH.nc ("<Name><cohort>_Catch", "<Name><cohort>_Discards", numbers x mean individual weight).
##                 If those are not in the output the catch is shared in proportion to the stanza biomass.
##   M0          - biomass weighted mean across the cohorts in the stanza of the SpecificMort.txt M1 values.

########################### LIBRARY FUNCTIONS ###########################

# Load (installing if needed) the libraries used by the Atlantis conversion
##' @title Load libraries used by the Atlantis to EwE conversion
##' @author Beth Fulton
atl_load_libraries <- function() {
  for (this_package in c("ncdf4", "stringr", "stringdist")) {
    if (!requireNamespace(this_package, quietly = TRUE)) {
      install.packages(this_package, dependencies = TRUE)
    }
  }
  library(ncdf4)
}

########################### ATLANTIS INPUT FILE FUNCTIONS ###########################

# Find a column in a dataframe ignoring case (Atlantis group file headers vary between versions)
##' @title Case insensitive column lookup
##' @author Beth Fulton
atl_col <- function(df, col_name, default = NA) {
  this_col <- which(tolower(names(df)) == tolower(col_name))
  if (length(this_col) == 0) {
    return(rep(default, nrow(df)))
  }
  return(df[[this_col[1]]])
}

# Read the Atlantis functional groups csv file - keeps only the groups that are turned on
# The ID given to each group (1..n in the order of the groups file) is the EwE style group ID
##' @title Read Atlantis functional groups file
##' @return Dataframe of active groups
##' @author Beth Fulton
atl_read_groups <- function(groups_file) {
  fgs <- read.csv(groups_file, header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)

  group_type <- atl_col(fgs, "GroupType")
  if (all(is.na(group_type))) {
    group_type <- atl_col(fgs, "InvertType")  # Older Atlantis group files
  }

  groups <- data.frame(
    Code = atl_col(fgs, "Code"),
    Name = atl_col(fgs, "Name"),
    IsTurnedOn = atl_col(fgs, "IsTurnedOn", 1),
    NumCohorts = atl_col(fgs, "NumCohorts", 1),
    NumGeneTypes = atl_col(fgs, "NumGeneTypes", 1),
    NumAgeClassSize = atl_col(fgs, "NumAgeClassSize", 1),
    IsFished = atl_col(fgs, "IsFished", 0),
    IsImpacted = atl_col(fgs, "IsImpacted", 0),
    isTAC = atl_col(fgs, "isTAC", 0),
    GroupType = toupper(group_type),
    stringsAsFactors = FALSE
  )
  groups$FileRow <- seq_len(nrow(groups))  # Position in the groups file (and so in the pPREY vectors)
  groups <- groups[groups$IsTurnedOn > 0, ]
  groups$NumCohorts <- groups$NumCohorts * groups$NumGeneTypes
  groups$ID <- seq_len(nrow(groups))

  # Broad group categories used in the conversion
  groups$IsDetritus <- as.integer(groups$GroupType %in% c("LAB_DET", "REF_DET", "CARRION") |
                                    grepl("_DET$", groups$GroupType))
  groups$IsProducer <- as.integer(groups$GroupType %in% c("LG_PHY", "SM_PHY", "DINOFLAG", "MICROPHTYBENTHOS",
                                                          "PHYTOBEN", "SEAGRASS", "TURF", "ICE_DIATOMS",
                                                          "MANGROVES", "SEAWEED", "PHYTO"))
  # Atlantis treats FISH_INVERT groups as fish (so as vertebrates)
  groups$IsVertebrate <- as.integer(groups$GroupType %in% c("FISH", "FISH_INVERT", "SHARK", "MAMMAL", "BIRD"))

  # The Atlantis group each output group comes from and the cohorts it holds - the same as the group itself
  # unless the groups are split into stanzas (see atl_build_stanzas())
  groups$Parent <- seq_len(nrow(groups))
  groups$AtlCode <- groups$Code
  groups$AtlName <- groups$Name
  groups$IsStanza <- 0
  groups$FirstCohort <- 1
  groups$LastCohort <- groups$NumCohorts

  rownames(groups) <- groups$Code
  return(groups)
}

# Read a value from an Atlantis prm file (the first number after the parameter name)
##' @title Read parameter value from Atlantis prm file
##' @author Beth Fulton
atl_prm_value <- function(prm_lines, param_name, default = NA) {
  this_line <- grep(paste0("^", param_name, "[[:space:]]"), prm_lines, value = TRUE)
  if (length(this_line) == 0) {
    return(default)
  }
  this_value <- suppressWarnings(as.numeric(strsplit(trimws(this_line[1]), "[[:space:]]+")[[1]][2]))
  if (is.na(this_value)) {
    return(default)
  }
  return(this_value)
}

# Read the vectors of the parameters matching a pattern from an Atlantis prm file (the values are on the
# first non comment line after the parameter name, e.g. the pPREY availabilities)
##' @title Read vector parameters from Atlantis prm file
##' @return Named list of numeric vectors
##' @author Beth Fulton
atl_prm_vectors <- function(prm_lines, name_pattern) {
  these_lines <- grep(paste0("^", name_pattern), prm_lines)
  vals <- lapply(these_lines, function(i) {
    for (j in (i + 1):length(prm_lines)) {
      this_text <- trimws(sub("#.*$", "", prm_lines[j]))
      if (nchar(this_text) > 0) break
    }
    suppressWarnings(as.numeric(strsplit(this_text, "[[:space:]]+")[[1]]))
  })
  names(vals) <- sub("[[:space:]].*$", "", trimws(prm_lines[these_lines]))
  return(vals)
}

# Point in polygon (ray casting) - used to identify the boundary boxes the same way Atlantis does
##' @title Point in polygon test
##' @author Beth Fulton
atl_point_in_poly <- function(x, y, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- FALSE
  j <- n
  for (i in 1:n) {
    if (((poly_y[i] > y) != (poly_y[j] > y)) &&
        (x < (poly_x[j] - poly_x[i]) * (y - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])) {
      inside <- !inside
    }
    j <- i
  }
  return(inside)
}

# Read the box information from the bgm file. Boundary boxes are identified as Atlantis does (box
# inside point outside of the bnd_vert polygon, or botz >= 0)
##' @title Read Atlantis bgm file
##' @return Dataframe of box number, area (m2), depth and whether it is a boundary box
##' @author Beth Fulton
atl_read_bgm <- function(bgm_file) {
  bgm <- readLines(bgm_file, warn = FALSE)

  get_box_field <- function(field) {
    these_lines <- grep(paste0("^box[0-9]+\\.", field, "[[:space:]]"), bgm, value = TRUE)
    box_id <- as.integer(sub("^box([0-9]+)\\..*", "\\1", these_lines))
    vals <- strsplit(trimws(sub(paste0("^box[0-9]+\\.", field), "", these_lines)), "[[:space:]]+")
    list(box_id = box_id, vals = vals)
  }

  area <- get_box_field("area")
  botz <- get_box_field("botz")
  inside <- get_box_field("inside")

  boxes <- data.frame(box = area$box_id, area = as.numeric(sapply(area$vals, `[`, 1)))
  boxes$botz <- as.numeric(sapply(botz$vals, `[`, 1))[match(boxes$box, botz$box_id)]
  boxes$inside_x <- as.numeric(sapply(inside$vals, `[`, 1))[match(boxes$box, inside$box_id)]
  boxes$inside_y <- as.numeric(sapply(inside$vals, `[`, 2))[match(boxes$box, inside$box_id)]
  boxes <- boxes[order(boxes$box), ]

  # Boundary polygon
  bnd_lines <- grep("^bnd_vert", bgm, value = TRUE)
  boxes$is_boundary <- boxes$botz >= 0
  if (length(bnd_lines) > 2) {
    bnd <- do.call(rbind, lapply(strsplit(trimws(bnd_lines), "[[:space:]]+"), function(v) as.numeric(v[2:3])))
    for (i in seq_len(nrow(boxes))) {
      if (!atl_point_in_poly(boxes$inside_x[i], boxes$inside_y[i], bnd[, 1], bnd[, 2])) {
        boxes$is_boundary[i] <- TRUE
      }
    }
  }

  return(boxes)
}

########################### TIME HANDLING FUNCTIONS ###########################

# Atlantis Time columns are days since the start of the run (a model year is 365 days).
# Outputs that report something accumulated over the period ending at Time t (consumption, diet,
# mortality, catch) are allocated to model year ceiling(t/365) - 1 (so t = 365 is year 0)
##' @title Convert Atlantis time (days) of period outputs to model year index (0 based)
##' @author Beth Fulton
atl_period_year <- function(time_days) {
  yr <- ceiling(time_days / 365 - 1e-6) - 1
  yr[time_days <= 0] <- -1
  return(yr)
}

# Keep only the first row for each key (Atlantis writes some files twice at year ends and at the
# end of the run - the first row is the one with the full period accumulated)
##' @title Remove duplicated rows written at the same time
##' @author Beth Fulton
atl_first_per_key <- function(df, key_cols) {
  df[!duplicated(df[, key_cols, drop = FALSE]), , drop = FALSE]
}

########################### ATLANTIS OUTPUT LOAD FUNCTIONS ###########################

# Annual biomass (tonnes) - mean of the BiomIndx snapshots from the start to the end of each model year
##' @title Load annual biomass from BiomIndx.txt
##' @return Matrix years x groups (tonnes)
##' @author Beth Fulton
atl_annual_biomass <- function(biomindx_file, groups, nyears) {
  bi <- read.table(biomindx_file, header = TRUE, check.names = FALSE)
  bi <- atl_first_per_key(bi, "Time")

  missing_groups <- setdiff(groups$Code, names(bi))
  if (length(missing_groups) > 0) {
    warning(paste("Groups missing from BiomIndx.txt (biomass set to 0):", paste(missing_groups, collapse = " ")))
  }

  bio <- matrix(0, nrow = nyears, ncol = nrow(groups), dimnames = list(NULL, groups$Code))
  for (k in 1:nyears) {
    these_rows <- which(bi$Time >= (k - 1) * 365 - 1e-6 & bi$Time <= k * 365 + 1e-6)
    for (g in intersect(groups$Code, names(bi))) {
      bio[k, g] <- mean(bi[these_rows, g])
    }
  }

  return(bio)
}

# Box volumes (m3) through time from the main output netcdf file, with boundary boxes set to zero
##' @title Load box volumes from Atlantis output netcdf file
##' @return List of volume array (layer x box x time) and the time (days) of each record
##' @author Beth Fulton
atl_box_volumes <- function(out_nc_file, boxes) {
  nc <- nc_open(out_nc_file)
  on.exit(nc_close(nc))
  vol <- ncvar_get(nc, "volume", collapse_degen = FALSE)  # z, b, t
  vol_time <- ncvar_get(nc, "t") / 86400

  boundary_boxes <- boxes$box[boxes$is_boundary] + 1   # Atlantis box numbers start at 0
  boundary_boxes <- boundary_boxes[boundary_boxes <= dim(vol)[2]]
  vol[, boundary_boxes, ] <- 0
  vol[is.na(vol)] <- 0

  return(list(vol = vol, vol_time = vol_time))
}

# Biomass and individual size of each cohort of the age structured groups from the main output netcdf file.
# Biomass (tonnes) = Nums * (StructN + ResN) * X_CN * k_wetdry / 1e9 with boundary boxes excluded. The annual
# values are the mean of the snapshots from the start to the end of each model year (as for BiomIndx.txt)
##' @title Load annual biomass and individual size per cohort
##' @return List of arrays years x group x cohort of biomass (tonnes) and individual StructN and ResN (mg N,
##'   numbers weighted mean), the biomass per box (years x group x cohort x box) and whether each group has cohort tracers
##' @author Beth Fulton
atl_cohort_info <- function(out_nc_file, groups, boxes, nyears, mgN_to_tonnes) {
  nc <- nc_open(out_nc_file)
  on.exit(nc_close(nc))
  nc_time <- ncvar_get(nc, "t") / 86400
  nbox <- nc$dim$b$len
  boundary_boxes <- boxes$box[boxes$is_boundary] + 1   # Atlantis box numbers start at 0
  boundary_boxes <- boundary_boxes[boundary_boxes <= nbox]
  year_rec <- lapply(1:nyears, function(k) which(nc_time >= (k - 1) * 365 - 1e-6 & nc_time <= k * 365 + 1e-6))

  max_cohort <- max(groups$NumCohorts)
  bio <- array(0, dim = c(nyears, nrow(groups), max_cohort), dimnames = list(NULL, groups$Code, NULL))
  sn <- bio
  rn <- bio
  bio_box <- array(0, dim = c(nyears, nrow(groups), max_cohort, nbox), dimnames = list(NULL, groups$Code, NULL, NULL))
  has_tracers <- rep(FALSE, nrow(groups))

  for (g in which(groups$NumCohorts > 1)) {
    var_stem <- paste0(groups$Name[g], seq_len(groups$NumCohorts[g]))
    if (!all(paste0(var_stem, c("_Nums", "_StructN", "_ResN")) %in% names(nc$var))) next
    has_tracers[g] <- TRUE

    for (cohort in seq_along(var_stem)) {
      get_var <- function(suffix) {
        x <- ncvar_get(nc, paste0(var_stem[cohort], suffix), collapse_degen = FALSE)  # layer x box x time
        x[is.na(x)] <- 0
        x
      }
      nums <- get_var("_Nums")
      nums[, boundary_boxes, ] <- 0
      sn_ind <- get_var("_StructN")
      rn_ind <- get_var("_ResN")
      ntime <- dim(nums)[3]
      layer_sum <- function(x) matrix(colSums(matrix(x, nrow = dim(nums)[1])), nrow = nbox, ncol = ntime)  # box x time
      bio_bt <- layer_sum(nums * (sn_ind + rn_ind)) * mgN_to_tonnes
      nums_t <- colSums(layer_sum(nums))
      sn_t <- colSums(layer_sum(nums * sn_ind))
      rn_t <- colSums(layer_sum(nums * rn_ind))

      for (k in 1:nyears) {
        these_rec <- year_rec[[k]]
        if (length(these_rec) == 0) next
        bio_box[k, g, cohort, ] <- rowMeans(bio_bt[, these_rec, drop = FALSE])
        bio[k, g, cohort] <- sum(bio_box[k, g, cohort, ])
        if (sum(nums_t[these_rec]) > 0) {
          sn[k, g, cohort] <- sum(sn_t[these_rec]) / sum(nums_t[these_rec])
          rn[k, g, cohort] <- sum(rn_t[these_rec]) / sum(nums_t[these_rec])
        }
      }
    }
  }

  return(list(bio = bio, bio_box = bio_box, sn = sn, rn = rn, has_tracers = has_tracers))
}

# Total consumption (tonnes per year) of each predator cohort from the Eat and Grazing variables in PROD.nc
##' @title Load annual consumption per predator cohort
##' @param by_box If TRUE also return the consumption per box (used to share the consumption of split prey groups
##'   across their stanzas)
##' @return Array years x predator x cohort (tonnes per year), or if by_box a list of that (total) and the
##'   array years x predator x cohort x box (box)
##' @author Beth Fulton
atl_annual_consumption <- function(prod_nc_file, vol_info, groups, nyears, mgN_to_tonnes, by_box = FALSE) {
  nc <- nc_open(prod_nc_file)
  on.exit(nc_close(nc))

  prod_time <- ncvar_get(nc, "t") / 86400
  prod_year <- atl_period_year(prod_time)

  # Match each PROD.nc record to the nearest record in the main output file (for the volumes)
  vol_index <- sapply(prod_time, function(t) which.min(abs(vol_info$vol_time - t)))
  vol_wc <- apply(vol_info$vol, c(2, 3), sum)  # box x time
  nlayer <- dim(vol_info$vol)[1]
  nbox <- dim(vol_info$vol)[2]

  max_cohort <- max(groups$NumCohorts)
  cons <- array(0, dim = c(nyears, nrow(groups), max_cohort),
                dimnames = list(NULL, groups$Code, NULL))
  cons_box <- NULL
  if (by_box) {
    cons_box <- array(0, dim = c(nyears, nrow(groups), max_cohort, nbox), dimnames = list(NULL, groups$Code, NULL, NULL))
  }
  has_cons_data <- rep(FALSE, nrow(groups))

  for (g in seq_len(nrow(groups))) {
    if (groups$IsDetritus[g] > 0 || groups$IsProducer[g] > 0) next

    # Age structured groups have one Eat variable per cohort, biomass pools have a Grazing variable
    var_names <- paste0(groups$Name[g], seq_len(groups$NumCohorts[g]), "_Eat")
    if (!all(var_names %in% names(nc$var))) {
      var_names <- paste0(groups$Name[g], "Grazing")
      if (!(var_names %in% names(nc$var))) next
    }
    has_cons_data[g] <- TRUE

    for (cohort in seq_along(var_names)) {
      eat <- ncvar_get(nc, var_names[cohort], collapse_degen = FALSE)
      eat[is.na(eat)] <- 0
      if (length(dim(eat)) == 3) {
        # layer x box x time - use the volume of each layer
        eat_mass <- eat * vol_info$vol[, , vol_index, drop = FALSE]
        daily_box <- matrix(colSums(matrix(eat_mass, nrow = nlayer)), nrow = nbox)
      } else {
        # box x time - use the volume of the whole box
        daily_box <- eat * vol_wc[, vol_index, drop = FALSE]
      }
      daily <- colSums(daily_box)
      for (k in 1:nyears) {
        these_rec <- which(prod_year == (k - 1))
        if (length(these_rec) > 0) {
          cons[k, g, cohort] <- mean(daily[these_rec]) * 365 * mgN_to_tonnes
          if (by_box) {
            cons_box[k, g, cohort, ] <- rowMeans(daily_box[, these_rec, drop = FALSE]) * 365 * mgN_to_tonnes
          }
        }
      }
    }
  }

  no_data <- groups$Code[!has_cons_data & groups$IsDetritus < 1 & groups$IsProducer < 1]
  if (length(no_data) > 0) {
    message(paste("Consumers with no Eat or Grazing variable in PROD.nc (consumption set to 0):",
                  paste(no_data, collapse = " ")))
  }

  if (by_box) {
    return(list(total = cons, box = cons_box))
  }
  return(cons)
}

# Annual diet proportions per predator cohort from DietCheck.txt
##' @title Load annual diet composition per predator cohort
##' @return Array years x predator x cohort x prey (proportions, NA where no diet information)
##' @author Beth Fulton
atl_annual_diet <- function(dietcheck_file, groups, nyears) {
  dc <- read.table(dietcheck_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  id_cols <- intersect(c("Time", "Predator", "Cohort", "Stock", "Updated"), names(dc))
  prey_cols <- intersect(groups$Code, setdiff(names(dc), id_cols))
  if (!("Cohort" %in% names(dc))) dc$Cohort <- 0
  if (!("Stock" %in% names(dc))) dc$Stock <- 0
  dc <- atl_first_per_key(dc, c("Time", "Predator", "Cohort", "Stock"))
  dc$Year <- atl_period_year(dc$Time)
  dc <- dc[dc$Year >= 0 & dc$Year < nyears & dc$Predator %in% groups$Code, ]

  # Drop rows with no diet recorded (predator absent)
  dc_prop <- as.matrix(dc[, prey_cols])
  row_tot <- rowSums(dc_prop)
  keep <- row_tot > 0
  dc <- dc[keep, ]
  dc_prop <- dc_prop[keep, , drop = FALSE] / row_tot[keep]

  max_cohort <- max(groups$NumCohorts)
  diet <- array(NA, dim = c(nyears, nrow(groups), max_cohort, nrow(groups)),
                dimnames = list(NULL, groups$Code, NULL, groups$Code))

  key <- paste(dc$Year, dc$Predator, dc$Cohort)
  for (this_key in unique(key)) {
    these_rows <- which(key == this_key)
    k <- dc$Year[these_rows[1]] + 1
    pred <- dc$Predator[these_rows[1]]
    cohort <- dc$Cohort[these_rows[1]] + 1
    if (cohort > max_cohort) next
    diet[k, pred, cohort, prey_cols] <- colMeans(dc_prop[these_rows, , drop = FALSE])  # Average across stocks and times in the year
  }
  diet[, , , setdiff(groups$Code, prey_cols)][is.na(diet[, , , setdiff(groups$Code, prey_cols)])] <- 0

  return(diet)
}

# Combine consumption per predator cohort and diet proportions into annual prey x predator consumption
##' @title Calculate annual consumption matrix
##' @param stanza Stanza information from atl_setup_stanzas() - if given the matrix is between the stanza (output)
##'   groups rather than the Atlantis groups
##' @return Array years x prey x predator (tonnes per year)
##' @author Beth Fulton
atl_consumption_matrix <- function(cons, diet, groups, bio_t, stanza = NULL) {
  nyears <- dim(cons)[1]
  nsp <- nrow(groups)
  out_codes <- if (is.null(stanza)) groups$Code else stanza$groups$Code
  Q <- array(0, dim = c(nyears, length(out_codes), length(out_codes)), dimnames = list(NULL, out_codes, out_codes))
  not_allocated <- c()

  # Bacteria are not in DietCheck.txt - they take up labile and refractory detritus, so if there is
  # no diet information for them split their uptake across those in proportion to the detritus biomass
  bact <- which(groups$GroupType %in% c("PL_BACT", "SED_BACT", "ICE_BACT"))
  det <- which(groups$GroupType %in% c("LAB_DET", "REF_DET"))
  for (j in bact) {
    for (k in 1:nyears) {
      for (cohort in 1:dim(diet)[3]) {
        if (all(is.na(diet[k, j, cohort, ])) && length(det) > 0) {
          diet[k, j, cohort, ] <- 0
          diet[k, j, cohort, det] <- pmax(bio_t[k, det], 1e-20) / sum(pmax(bio_t[k, det], 1e-20))
        }
      }
    }
  }

  for (j in seq_len(nsp)) {
    if (sum(cons[, j, ]) <= 0) next

    # Fall back diets - predator average for that year, then predator average over all years
    all_year_diet <- apply(diet[, j, , , drop = FALSE], 4, function(x) mean(x, na.rm = TRUE))

    for (k in 1:nyears) {
      year_diet <- apply(diet[k, j, , , drop = FALSE], 4, function(x) mean(x, na.rm = TRUE))
      for (cohort in 1:dim(cons)[3]) {
        this_cons <- cons[k, j, cohort]
        if (this_cons <= 0) next
        this_diet <- diet[k, j, cohort, ]
        if (any(is.na(this_diet)) || sum(this_diet) <= 0) this_diet <- year_diet
        if (any(is.na(this_diet)) || sum(this_diet) <= 0) this_diet <- all_year_diet
        if (any(is.na(this_diet)) || sum(this_diet) <= 0) {
          not_allocated <- c(not_allocated, groups$Code[j])
          next
        }
        prey_cons <- this_cons * this_diet / sum(this_diet)
        if (is.null(stanza)) {
          Q[k, , j] <- Q[k, , j] + prey_cons
        } else {
          pred_out <- stanza$cohort_map[j, cohort]
          Q[k, , pred_out] <- Q[k, , pred_out] + atl_split_prey_cons(prey_cons, stanza, k, j, cohort)
        }
      }
    }
  }

  if (length(not_allocated) > 0) {
    warning(paste("No DietCheck information so consumption could not be allocated to prey for:",
                  paste(unique(not_allocated), collapse = " ")))
  }

  return(Q)
}

# Non-predation natural mortality (M0) from the M1 columns of SpecificMort.txt
##' @title Load annual M0 per group
##' @param cohort_bio Biomass per cohort (from atl_cohort_info()) used to weight the cohorts within a stanza
##' @return Matrix years x groups (yr-1)
##' @author Beth Fulton
atl_annual_m0 <- function(specificmort_file, groups, nyears, cohort_bio = NULL) {
  m0 <- matrix(0, nrow = nyears, ncol = nrow(groups), dimnames = list(NULL, groups$Code))
  if (is.null(specificmort_file) || !file.exists(specificmort_file)) {
    warning("No SpecificMort.txt file so M0 set to 0")
    return(m0)
  }

  sm <- read.table(specificmort_file, header = TRUE, check.names = FALSE)
  sm <- atl_first_per_key(sm, "Time")
  sm$Year <- atl_period_year(sm$Time)

  m1_cols <- grep("-M1$", names(sm), value = TRUE)
  m1_code <- sub("-.*$", "", m1_cols)
  m1_cohort <- suppressWarnings(as.integer(sub("^[^-]+-([0-9]+)-.*$", "\\1", m1_cols))) + 1  # Columns are <Code>-<cohort>-S<stock>-M1

  for (k in 1:nyears) {
    this_row <- which(sm$Year == (k - 1))
    if (length(this_row) == 0) next
    for (o in which(groups$AtlCode %in% m1_code)) {
      these_cols <- m1_code == groups$AtlCode[o]
      if (groups$IsStanza[o] > 0) {
        these_cols <- these_cols & m1_cohort %in% (groups$FirstCohort[o]:groups$LastCohort[o])
      }
      m1 <- unlist(sm[this_row[1], m1_cols[these_cols]])
      m1 <- pmin(pmax(m1, 0), 0.99)
      if (groups$IsStanza[o] > 0 && !is.null(cohort_bio)) {
        # Biomass weighted across the cohorts in the stanza (and averaged across stocks)
        wgt <- cohort_bio[k, groups$Parent[o], m1_cohort[these_cols]]
        if (sum(wgt) > 0) {
          m0[k, o] <- sum(wgt * -log(1 - m1)) / sum(wgt)
          next
        }
      }
      m0[k, o] <- mean(-log(1 - m1))
    }
  }

  return(m0)
}

# Fleet list - from the fisheries csv file if given
##' @title Read Atlantis fisheries file
##' @return Dataframe of fleets
##' @author Beth Fulton
atl_read_fleets <- function(fisheries_file) {
  if (is.null(fisheries_file) || !file.exists(fisheries_file)) {
    return(NULL)
  }
  fsh <- read.csv(fisheries_file, header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  fleets <- data.frame(Code = atl_col(fsh, "Code"), Name = atl_col(fsh, "Name"), stringsAsFactors = FALSE)
  fleets$ID <- seq_len(nrow(fleets))
  return(fleets)
}

# Read one of the per fishery text files (CatchPerFishery.txt or DiscardPerFishery.txt)
##' @title Load per fishery catch or discards text file
##' @return Array years x fleet x group (tonnes per year)
##' @author Beth Fulton
atl_read_per_fishery <- function(this_file, fleets, groups, nyears) {
  out <- array(0, dim = c(nyears, nrow(fleets), nrow(groups)),
               dimnames = list(NULL, fleets$Code, groups$Code))
  if (!file.exists(this_file)) {
    return(out)
  }
  df <- read.table(this_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  df <- atl_first_per_key(df, c("Time", "Fishery"))
  df$Year <- atl_period_year(df$Time)
  df <- df[df$Year >= 0 & df$Year < nyears, ]

  # Fishery column can be the fleet code or a (0 based) fleet number
  if (all(df$Fishery %in% fleets$Code)) {
    df$fleet_index <- match(df$Fishery, fleets$Code)
  } else {
    df$fleet_index <- as.integer(df$Fishery) + 1
  }

  these_groups <- intersect(groups$Code, names(df))
  for (i in seq_len(nrow(df))) {
    out[df$Year[i] + 1, df$fleet_index[i], these_groups] <- out[df$Year[i] + 1, df$fleet_index[i], these_groups] +
      unlist(df[i, these_groups])
  }
  return(out)
}

# Annual landings and discards per fleet per group
##' @title Load annual Atlantis fisheries output
##' @return List of fleets and arrays (years x fleet x group, tonnes per year) of landings and discards
##' @author Beth Fulton
atl_annual_catch <- function(atlantis_out_dir, output_prefix, fleets, groups, boxes, nyears, catch_source = "auto") {
  use_source <- function(this_source) catch_source %in% c("auto", this_source)
  out_file <- function(suffix) file.path(atlantis_out_dir, paste0(output_prefix, suffix))

  if (is.null(fleets) && use_source("txt") && file.exists(out_file("CatchPerFishery.txt"))) {
    df <- read.table(out_file("CatchPerFishery.txt"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    fleets <- data.frame(Code = unique(df$Fishery), Name = unique(df$Fishery), stringsAsFactors = FALSE)
    fleets$ID <- seq_len(nrow(fleets))
  }

  if (!is.null(fleets) && use_source("txt") && file.exists(out_file("CatchPerFishery.txt"))) {
    message("Fisheries data taken from CatchPerFishery.txt and DiscardPerFishery.txt")
    # CatchPerFishery.txt is the total catch (landings + discards) - checked against Catch.txt (landings)
    # and Discard.txt - so the landings are the catch less the discards. Groups that are not isFished
    # can still be discarded but are left out of CatchPerFishery.txt, hence the floor at zero
    catch_pf <- atl_read_per_fishery(out_file("CatchPerFishery.txt"), fleets, groups, nyears)
    discards <- atl_read_per_fishery(out_file("DiscardPerFishery.txt"), fleets, groups, nyears)
    landings <- pmax(catch_pf - discards, 0)
    return(list(fleets = fleets, landings = landings, discards = discards))
  }

  if (!is.null(fleets) && use_source("nc") && file.exists(out_file("CATCH.nc"))) {
    message("Fisheries data taken from CATCH.nc")
    landings <- array(0, dim = c(nyears, nrow(fleets), nrow(groups)), dimnames = list(NULL, fleets$Code, groups$Code))
    discards <- landings
    nc <- nc_open(out_file("CATCH.nc"))
    on.exit(nc_close(nc))
    catch_year <- atl_period_year(ncvar_get(nc, "t") / 86400)
    internal_boxes <- boxes$box[!boxes$is_boundary] + 1
    for (g in seq_len(nrow(groups))) {
      for (f in seq_len(nrow(fleets))) {
        for (this_type in c("Catch", "Discard")) {
          var_name <- paste0(groups$Code[g], "_", this_type, "_FC", f)
          if (!(var_name %in% names(nc$var))) next
          val <- ncvar_get(nc, var_name, collapse_degen = FALSE)   # box x time
          val[is.na(val)] <- 0
          per_rec <- colSums(val[internal_boxes[internal_boxes <= nrow(val)], , drop = FALSE])
          for (k in 1:nyears) {
            tot <- sum(per_rec[catch_year == (k - 1)])
            if (this_type == "Catch") landings[k, f, g] <- tot else discards[k, f, g] <- tot
          }
        }
      }
    }
    return(list(fleets = fleets, landings = landings, discards = discards))
  }

  if (use_source("total") && file.exists(out_file("Catch.txt"))) {
    message("Fisheries data taken from Catch.txt and Discard.txt (aggregated into a single fleet)")
    fleets <- data.frame(Code = "All_fleets", Name = "All_fleets", ID = 1, stringsAsFactors = FALSE)
    landings <- array(0, dim = c(nyears, 1, nrow(groups)), dimnames = list(NULL, fleets$Code, groups$Code))
    discards <- landings
    for (this_type in c("Catch", "Discard")) {
      if (!file.exists(out_file(paste0(this_type, ".txt")))) next
      df <- read.table(out_file(paste0(this_type, ".txt")), header = TRUE, check.names = FALSE)
      df <- atl_first_per_key(df, "Time")
      df$Year <- atl_period_year(df$Time)
      these_groups <- intersect(groups$Code, names(df))
      for (k in 1:nyears) {
        this_row <- which(df$Year == (k - 1))
        if (length(this_row) == 0) next
        if (this_type == "Catch") {
          landings[k, 1, these_groups] <- unlist(df[this_row[1], these_groups])
        } else {
          discards[k, 1, these_groups] <- unlist(df[this_row[1], these_groups])
        }
      }
    }
    return(list(fleets = fleets, landings = landings, discards = discards))
  }

  message("No fisheries output found - writing a single dummy fleet with zero catch")
  fleets <- data.frame(Code = "No_fishing", Name = "No_fishing", ID = 1, stringsAsFactors = FALSE)
  landings <- array(0, dim = c(nyears, 1, nrow(groups)), dimnames = list(NULL, fleets$Code, groups$Code))
  return(list(fleets = fleets, landings = landings, discards = landings))
}

########################### STANZA FUNCTIONS ###########################

# Split the age structured groups into stanzas (like EwE multi-stanza groups). Each stanza is an output group
# holding a contiguous range of cohorts. By default the split is the Atlantis juvenile / adult split (cohorts
# before <Code>_age_mat are juveniles). Groups that are not split stay as a single output group
##' @title Define the stanza (output) groups
##' @param has_tracers Whether each group has cohort tracers in the main output .nc (only those can be split)
##' @param stanza_groups Codes of the groups to split (NULL for all the groups with cohort tracers)
##' @param stanza_breaks Optional named list (by group code) of the first cohort of each stanza after the first
##' @param stanza_labels Optional named list (by group code) of the labels of the stanzas
##' @return List of the output groups (a groups dataframe with one row per stanza), the output group of each
##'   group x cohort (cohort_map), whether each group is split and the Atlantis life stage (0 juvenile, 1 adult) of each cohort
##' @author Beth Fulton
atl_build_stanzas <- function(groups, bio_prm_lines, has_tracers, stanza_groups = NULL, stanza_breaks = NULL,
                              stanza_labels = NULL) {
  nsp <- nrow(groups)
  max_cohort <- max(groups$NumCohorts)

  # Atlantis life stage of each cohort (as set in atbiolsetup.c) - used for the pPREY availabilities whatever the
  # stanza split. Age class = cohort / NumGeneTypes and a cohort is an adult once its age class reaches age_mat
  age_mat <- sapply(groups$Code, function(cd) atl_prm_value(bio_prm_lines, paste0(cd, "_age_mat"), 0))
  life_stage <- matrix(NA, nsp, max_cohort, dimnames = list(groups$Code, NULL))
  for (g in seq_len(nsp)) {
    age_class <- floor((seq_len(groups$NumCohorts[g]) - 1) / groups$NumGeneTypes[g])
    life_stage[g, seq_len(groups$NumCohorts[g])] <- as.integer(groups$NumCohorts[g] > 1 & age_class >= age_mat[g])
  }

  if (is.null(stanza_groups)) {
    stanza_groups <- groups$Code[has_tracers]
  }
  not_found <- setdiff(stanza_groups, groups$Code[has_tracers])
  if (length(not_found) > 0) {
    warning(paste("These groups can not be split into stanzas (not active or no cohort tracers in the output .nc):",
                  paste(not_found, collapse = " ")))
  }
  not_split <- setdiff(names(stanza_breaks), stanza_groups)
  if (length(not_split) > 0) {
    warning(paste("stanza_breaks given for groups that are not being split (ignored):", paste(not_split, collapse = " ")))
  }

  out_rows <- list()
  cohort_map <- matrix(NA, nsp, max_cohort, dimnames = list(groups$Code, NULL))
  is_split <- rep(FALSE, nsp)
  for (g in seq_len(nsp)) {
    cd <- groups$Code[g]
    cohorts <- seq_len(groups$NumCohorts[g])
    stanza_of_cohort <- rep(1, length(cohorts))
    labels <- NULL
    if (cd %in% stanza_groups && has_tracers[g]) {
      if (!is.null(stanza_breaks[[cd]])) {
        stanza_of_cohort <- findInterval(cohorts, sort(unique(c(1, stanza_breaks[[cd]]))))
      } else {
        stanza_of_cohort <- life_stage[g, cohorts] + 1
        labels <- c("juv", "adult")[sort(unique(stanza_of_cohort))]
      }
      stanza_of_cohort <- match(stanza_of_cohort, sort(unique(stanza_of_cohort)))  # Drop any empty stanzas
    }
    nstanza <- max(stanza_of_cohort)

    if (nstanza < 2) {
      if (cd %in% stanza_groups && has_tracers[g]) {
        message(paste(cd, "has only one stanza with this split so it is not split"))
      }
      out_rows[[length(out_rows) + 1]] <- groups[g, ]
      cohort_map[g, cohorts] <- length(out_rows)
      next
    }

    is_split[g] <- TRUE
    if (!is.null(stanza_labels[[cd]])) {
      if (length(stanza_labels[[cd]]) != nstanza) {
        stop(paste("stanza_labels for", cd, "has", length(stanza_labels[[cd]]), "labels but there are", nstanza, "stanzas"))
      }
      labels <- stanza_labels[[cd]]
    }
    if (is.null(labels)) {
      labels <- paste0("s", seq_len(nstanza))
    }
    for (s in seq_len(nstanza)) {
      this_row <- groups[g, ]
      this_row$Code <- paste0(cd, "_", labels[s])
      this_row$Name <- paste0(groups$Name[g], "_", labels[s])
      this_row$IsStanza <- 1
      this_row$FirstCohort <- min(cohorts[stanza_of_cohort == s])
      this_row$LastCohort <- max(cohorts[stanza_of_cohort == s])
      out_rows[[length(out_rows) + 1]] <- this_row
      cohort_map[g, cohorts[stanza_of_cohort == s]] <- length(out_rows)
    }
  }

  out_groups <- do.call(rbind, out_rows)
  out_groups$ID <- seq_len(nrow(out_groups))
  rownames(out_groups) <- out_groups$Code

  return(list(groups = out_groups, cohort_map = cohort_map, is_split = is_split, life_stage = life_stage))
}

# The pPREY availabilities (pSPVERTeat in Atlantis) of each prey group and stage to each predator group and stage.
# Age structured predators have pPREY<prey stage><Pred><pred stage> parameters, the others pPREY<Pred>
##' @title Read the pPREY availabilities from the biology prm file
##' @return Array predator x predator stage (1 juvenile, 2 adult) x prey x prey stage
##' @author Beth Fulton
atl_read_pprey <- function(bio_prm_lines, groups) {
  nsp <- nrow(groups)
  pprey <- atl_prm_vectors(bio_prm_lines, "pPREY")
  avail <- array(NA, dim = c(nsp, 2, nsp, 2), dimnames = list(groups$Code, NULL, groups$Code, NULL))
  for (j in seq_len(nsp)) {
    for (pred_stage in 1:2) {
      for (prey_stage in 1:2) {
        vals <- pprey[[paste0("pPREY", prey_stage, groups$Code[j], pred_stage)]]
        if (is.null(vals)) vals <- pprey[[paste0("pPREY", groups$Code[j])]]
        if (!is.null(vals) && length(vals) >= max(groups$FileRow)) {
          avail[j, pred_stage, , prey_stage] <- vals[groups$FileRow]
        }
      }
    }
    # Groups with a single life stage only have the stage 1 values
    for (prey_stage in 1:2) {
      if (all(is.na(avail[j, 2, , prey_stage]))) avail[j, 2, , prey_stage] <- avail[j, 1, , prey_stage]
    }
  }
  return(avail)
}

# Parameters of the size based feeding window of each group (see Avail_Fish() in atvertprocesses.c)
##' @title Read the feeding window parameters from the biology prm file
##' @return List of the window flags and the per group parameters
##' @author Beth Fulton
atl_read_gape <- function(bio_prm_lines, groups) {
  get_par <- function(prefix) {
    sapply(groups$Code, function(cd) atl_prm_value(bio_prm_lines, paste0(prefix, cd)))
  }
  list(hard = atl_prm_value(bio_prm_lines, "UseHardFeedingWindow", 0) > 0,
       bilogistic = atl_prm_value(bio_prm_lines, "UseBiLogisticFeedingWindow", 0) > 0,
       X_CN = atl_prm_value(bio_prm_lines, "X_CN", 5.7),
       k_wetdry = atl_prm_value(bio_prm_lines, "k_wetdry", 20),
       KLP = get_par("KLP_"), KUP = get_par("KUP_"), Kmax = get_par("Kmax_coefft_"),
       li_a = get_par("li_a_"), li_b = get_par("li_b_"))
}

# Size based availability of the prey cohorts to a predator cohort - the sizeScalar of Avail_Fish() in
# atvertprocesses.c. If the predator has no cohort sizes (e.g. invertebrates) or no window parameters it is 1
##' @title Size based feeding window
##' @param pred_sn,pred_rn Individual structural and reserve N of the predator cohort (mg N)
##' @param prey_sn,prey_rn Individual structural and reserve N of each prey cohort (mg N)
##' @return Vector of the availability scalar (0 - 1) of each prey cohort
##' @author Beth Fulton
atl_gape_scalar <- function(gape, pred_code, prey_code, pred_sn, pred_rn, prey_sn, prey_rn) {
  klp <- gape$KLP[[pred_code]]
  kup <- gape$KUP[[pred_code]]
  kmax <- gape$Kmax[[pred_code]]
  if (pred_sn <= 0 || is.na(klp) || is.na(kup)) {
    return(rep(1, length(prey_sn)))
  }

  if (gape$hard) {
    return(as.numeric(prey_sn >= klp * pred_sn & prey_sn <= kup * pred_sn))
  }
  if (is.na(kmax)) {
    return(rep(1, length(prey_sn)))
  }

  if (!gape$bilogistic) {
    # Humped (potentially skewed) curve on structural N
    rel_size <- prey_sn / (kup * pred_sn)
    scalar <- rel_size * exp(kmax * (1 - rel_size))
  } else {
    # Bi-logistic curve on length (cm) from the length-weight relationship
    li_a <- c(gape$li_a[[pred_code]], gape$li_a[[prey_code]])
    li_b <- c(gape$li_b[[pred_code]], gape$li_b[[prey_code]])
    if (any(is.na(c(li_a, li_b)))) {
      return(rep(1, length(prey_sn)))
    }
    get_len <- function(wgt_mgN, a, b) ((wgt_mgN * gape$k_wetdry * gape$X_CN / 1000) / a)^(1 / b)
    pred_len <- get_len(pred_sn + pred_rn, li_a[1], li_b[1])
    prey_len <- get_len(prey_sn + prey_rn, li_a[2], li_b[2])
    maxavail <- klp + (kup - klp) * 0.5
    small_side <- prey_len <= maxavail * pred_len
    xmid <- ifelse(small_side, (klp + (maxavail - klp) * 0.5) * pred_len, (kup - (kup - maxavail) * 0.5) * pred_len)
    scalar <- ifelse(small_side, 1 / (1 + exp(-kmax * (prey_len - xmid))), 1 / (1 + exp(kmax * (prey_len - xmid))))
  }

  scalar[is.na(scalar)] <- 0
  return(pmin(pmax(scalar, 0), 1))
}

# Proportion of the consumption of a prey group by one predator cohort taken from each prey cohort. In each box
# the prey cohorts are eaten in proportion to pPREY x size window x prey cohort biomass (Avail_Fish() without the
# refuge, stock, turbidity and pH scalars), and the boxes are weighted by the consumption of the predator cohort in them
##' @title Share of a predator cohort's consumption of a prey group from each prey cohort
##' @return Vector of proportions per prey cohort (sums to 1)
##' @author Beth Fulton
atl_prey_cohort_fractions <- function(stanza, k, j, cohort, i) {
  ag <- stanza$atl_groups
  ci <- stanza$cohort
  prey_cohorts <- seq_len(ag$NumCohorts[i])

  avail <- stanza$avail[j, stanza$life_stage[j, cohort] + 1, i, stanza$life_stage[i, prey_cohorts] + 1]
  avail[is.na(avail)] <- 1
  if (ag$IsVertebrate[i] > 0) {
    avail <- avail * atl_gape_scalar(stanza$gape, ag$Code[j], ag$Code[i], ci$sn[k, j, cohort], ci$rn[k, j, cohort],
                                     ci$sn[k, i, prey_cohorts], ci$rn[k, i, prey_cohorts])
  }

  prey_bio <- matrix(ci$bio_box[k, i, prey_cohorts, ], nrow = length(prey_cohorts))  # cohort x box
  # Only the relative availability matters - rescaling avoids underflow when all the prey are far outside the window
  if (max(avail) > 0) avail <- avail / max(avail)
  avail_box <- avail * prey_bio
  domain <- rowSums(avail_box)
  if (sum(domain) <= 0) {
    # Nothing available with the stage and size limits at the annual mean sizes - share by the prey biomass instead
    stanza$log$no_window <- c(stanza$log$no_window, paste0(ag$Code[j], cohort, "-", ag$Code[i]))
    domain <- rowSums(prey_bio)
    if (sum(domain) <= 0) domain <- rep(1, length(prey_cohorts))
    return(domain / sum(domain))
  }

  pred_box <- stanza$cons_box[k, j, cohort, ]
  if (sum(pred_box) <= 0) {
    return(domain / sum(domain))
  }
  # Consumption in boxes where none of the prey is available (e.g. the prey moved out during the year) is
  # shared using the whole of model domain availability
  box_tot <- colSums(avail_box)
  in_box <- box_tot > 0
  frac <- as.vector(avail_box[, in_box, drop = FALSE] %*% (pred_box[in_box] / box_tot[in_box])) +
    domain / sum(domain) * sum(pred_box[!in_box])
  return(frac / sum(frac))
}

# Map the consumption of one predator cohort per Atlantis prey group onto the stanza (output) prey groups
##' @title Share prey consumption across the prey stanzas
##' @return Vector of consumption per output group
##' @author Beth Fulton
atl_split_prey_cons <- function(prey_cons, stanza, k, j, cohort) {
  out <- numeric(nrow(stanza$groups))
  unsplit <- which(!stanza$is_split)
  out[stanza$cohort_map[unsplit, 1]] <- prey_cons[unsplit]
  for (i in which(stanza$is_split & prey_cons > 0)) {
    frac <- atl_prey_cohort_fractions(stanza, k, j, cohort, i)
    stanza_frac <- rowsum(frac, stanza$cohort_map[i, seq_along(frac)])
    these_out <- as.integer(rownames(stanza_frac))
    out[these_out] <- out[these_out] + prey_cons[i] * stanza_frac[, 1]
  }
  return(out)
}

# Biomass of each stanza (output) group - the BiomIndx.txt biomass of the group shared across its stanzas in
# proportion to the cohort biomasses from the main output .nc
##' @title Annual biomass per stanza
##' @return Matrix years x output groups (tonnes)
##' @author Beth Fulton
atl_stanza_biomass <- function(bio_t, stanza) {
  og <- stanza$groups
  ci <- stanza$cohort
  bio <- matrix(0, nrow = nrow(bio_t), ncol = nrow(og), dimnames = list(NULL, og$Code))
  for (o in seq_len(nrow(og))) {
    g <- og$Parent[o]
    if (og$IsStanza[o] < 1) {
      bio[, o] <- bio_t[, g]
      next
    }
    cohort_tot <- apply(ci$bio[, g, , drop = FALSE], 1, sum)
    stanza_tot <- apply(ci$bio[, g, og$FirstCohort[o]:og$LastCohort[o], drop = FALSE], 1, sum)
    bio[, o] <- bio_t[, g] * ifelse(cohort_tot > 0, stanza_tot / cohort_tot, 0)
  }

  # Check the cohort biomasses agree with BiomIndx.txt
  split_groups <- which(stanza$is_split)
  ratio <- sapply(split_groups, function(g) {
    cohort_tot <- apply(ci$bio[, g, , drop = FALSE], 1, sum)
    these <- bio_t[, g] > 0
    if (!any(these)) return(c(NA, NA))
    range(cohort_tot[these] / bio_t[these, g])
  })
  if (length(split_groups) > 0 && any(is.finite(ratio))) {
    message(paste0("Annual cohort biomass from the output .nc / BiomIndx.txt biomass of the split groups ranges from ",
                   signif(min(ratio, na.rm = TRUE), 4), " to ", signif(max(ratio, na.rm = TRUE), 4),
                   " (they differ where BiomIndx.txt is written less often than the output .nc and as BiomIndx.txt",
                   " includes migrants outside the model) - the BiomIndx.txt biomass is shared across the stanzas"))
  }
  return(bio)
}

# Share each group's catch per fleet across its stanzas using the catch (and discards) per cohort in CATCH.nc
# ("<Name><cohort>_Catch" and "<Name><cohort>_Discards", numbers) times the mean individual weight. If those are
# not in the output the catch is shared in proportion to the stanza biomass (the same share for every fleet)
##' @title Share the catch across the stanzas
##' @return The fishery list with the landings and discards arrays (years x fleet x output group)
##' @author Beth Fulton
atl_stanza_catch <- function(fishery, catch_nc_file, stanza, nyears, mgN_to_tonnes) {
  ag <- stanza$atl_groups
  og <- stanza$groups
  ci <- stanza$cohort
  cohort_catch <- list(Catch = array(0, dim = dim(ci$bio)), Discards = array(0, dim = dim(ci$bio)))

  if (file.exists(catch_nc_file)) {
    nc <- nc_open(catch_nc_file)
    on.exit(nc_close(nc))
    catch_year <- atl_period_year(ncvar_get(nc, "t") / 86400)
    for (g in which(stanza$is_split)) {
      for (this_type in names(cohort_catch)) {
        for (cohort in seq_len(ag$NumCohorts[g])) {
          var_name <- paste0(ag$Name[g], cohort, "_", this_type)
          if (!(var_name %in% names(nc$var))) next
          val <- ncvar_get(nc, var_name, collapse_degen = FALSE)
          val[is.na(val)] <- 0
          per_rec <- if (length(dim(val)) > 1) colSums(matrix(val, ncol = dim(val)[length(dim(val))])) else val
          ind_wgt <- (ci$sn[, g, cohort] + ci$rn[, g, cohort]) * mgN_to_tonnes
          for (k in 1:nyears) {
            cohort_catch[[this_type]][k, g, cohort] <- sum(per_rec[catch_year == (k - 1)]) * ind_wgt[k]
          }
        }
      }
    }
  }

  # Share of each stanza - from the cohort catch if there is any that year, otherwise the stanza biomass
  used_biomass <- c()
  get_share <- function(cc, cc_backup, total_catch) {
    share <- matrix(1, nrow = nyears, ncol = nrow(og))
    for (o in which(og$IsStanza > 0)) {
      g <- og$Parent[o]
      all_cohorts <- seq_len(ag$NumCohorts[g])
      these <- og$FirstCohort[o]:og$LastCohort[o]
      for (k in 1:nyears) {
        if (sum(cc[k, g, all_cohorts]) > 0) {
          share[k, o] <- sum(cc[k, g, these]) / sum(cc[k, g, all_cohorts])
        } else if (sum(cc_backup[k, g, all_cohorts]) > 0) {
          share[k, o] <- sum(cc_backup[k, g, these]) / sum(cc_backup[k, g, all_cohorts])
        } else {
          if (sum(total_catch[k, , g]) > 0) used_biomass <<- c(used_biomass, ag$Code[g])
          bio_tot <- sum(ci$bio[k, g, all_cohorts])
          share[k, o] <- if (bio_tot > 0) sum(ci$bio[k, g, these]) / bio_tot else 1 / sum(og$Parent == g)
        }
      }
    }
    share
  }
  expand <- function(x, share) {
    out <- array(0, dim = c(dim(x)[1:2], nrow(og)), dimnames = list(NULL, dimnames(x)[[2]], og$Code))
    for (o in seq_len(nrow(og))) out[, , o] <- x[, , og$Parent[o]] * share[, o]
    out
  }

  fishery$landings <- expand(fishery$landings, get_share(cohort_catch$Catch, cohort_catch$Discards, fishery$landings))
  fishery$discards <- expand(fishery$discards, get_share(cohort_catch$Discards, cohort_catch$Catch, fishery$discards))
  if (length(used_biomass) > 0) {
    message(paste("No catch per cohort in", basename(catch_nc_file), "for these groups so their catch is shared",
                  "across the stanzas by biomass:", paste(unique(used_biomass), collapse = " ")))
  }
  return(fishery)
}

# Set up everything needed to split the age structured groups into stanzas
##' @title Set up the stanza split
##' @return List used by atl_consumption_matrix() and the other stanza functions, including the output groups
##'   (groups) and their biomass (biomass, years x output groups in tonnes)
##' @author Beth Fulton
atl_setup_stanzas <- function(out_nc_file, groups, boxes, bio_prm_lines, nyears, mgN_to_tonnes, cons_box, bio_t,
                              stanza_groups = NULL, stanza_breaks = NULL, stanza_labels = NULL) {
  message("Reading the cohort biomass and sizes from the output .nc for the stanza split")
  cohort_info <- atl_cohort_info(out_nc_file, groups, boxes, nyears, mgN_to_tonnes)
  stanza <- atl_build_stanzas(groups, bio_prm_lines, cohort_info$has_tracers, stanza_groups, stanza_breaks, stanza_labels)
  stanza$atl_groups <- groups
  stanza$cohort <- cohort_info
  stanza$cons_box <- cons_box
  stanza$avail <- atl_read_pprey(bio_prm_lines, groups)
  stanza$gape <- atl_read_gape(bio_prm_lines, groups)
  stanza$log <- new.env()
  stanza$log$no_window <- character(0)
  stanza$biomass <- atl_stanza_biomass(bio_t, stanza)

  if (atl_prm_value(bio_prm_lines, "flag_olddiet", 0) > 0) {
    warning("flag_olddiet is on - the prey stanza split assumes the new diet calculations (pPREY used as pSPVERTeat)")
  }
  message(paste(sum(stanza$is_split), "groups split into", sum(stanza$groups$IsStanza), "stanzas -",
                nrow(stanza$groups), "groups in the EwE style output"))
  return(stanza)
}

########################### DERIVED QUANTITIES ###########################

# Trophic level per group per year from the realised diet. Primary producers and detritus are TL 1
##' @title Calculate trophic levels
##' @return Matrix years x groups
##' @author Beth Fulton
atl_trophic_levels <- function(Q, groups) {
  nyears <- dim(Q)[1]
  nsp <- nrow(groups)
  tl <- matrix(1, nrow = nyears, ncol = nsp, dimnames = list(NULL, groups$Code))
  base_groups <- which(groups$IsProducer > 0 | groups$IsDetritus > 0)

  # Diet summed over all years - used if a consumer did not feed in a year
  Q_all <- apply(Q, c(2, 3), sum)

  # Consumers that never fed in the run are given a trophic level of 2
  no_diet <- setdiff(which(colSums(Q_all) <= 0), base_groups)
  tl[, no_diet] <- 2
  if (length(no_diet) > 0) {
    warning(paste("No consumption recorded for these consumers so their trophic level is set to 2:",
                  paste(groups$Code[no_diet], collapse = " ")))
  }
  consumers <- setdiff(seq_len(nsp), c(base_groups, no_diet))
  if (length(consumers) == 0) return(tl)

  for (k in 1:nyears) {
    Qk <- Q[k, , ]
    for (j in consumers) {
      if (sum(Qk[, j]) <= 0) Qk[, j] <- Q_all[, j]
    }

    # Solve TL = 1 + sum(diet * TL_prey) for the consumers, the TL of all other groups is fixed
    DC <- sweep(Qk, 2, pmax(colSums(Qk), 1e-30), "/")  # prey x predator proportions
    A <- t(DC[consumers, consumers, drop = FALSE])
    rhs <- 1 + as.vector(t(DC[-consumers, consumers, drop = FALSE]) %*% tl[k, -consumers])
    this_tl <- tryCatch(solve(diag(length(consumers)) - A, rhs), error = function(e) NULL)
    if (is.null(this_tl)) {
      # Singular system (e.g. closed cannibalistic loop) so iterate instead
      this_tl <- rep(2, length(consumers))
      for (iter in 1:200) this_tl <- rhs + as.vector(A %*% this_tl)
    }
    tl[k, consumers] <- this_tl
  }

  return(tl)
}

########################### EWE FORMAT WRITING FUNCTIONS ###########################

# Format numbers as EwE does in its csv files
##' @title Format numbers for EwE style output
##' @author Beth Fulton
ewe_num <- function(x) {
  x[is.na(x) | !is.finite(x)] <- 0
  out <- formatC(x, digits = 7, format = "g")
  out <- trimws(out)
  return(out)
}

# The header block EwE puts at the top of every Ecosim csv file (the qmd skips these lines)
##' @title Create EwE style csv header
##' @author Beth Fulton
ewe_header <- function(data_label, model_name, model_source, scenario_name, start_year) {
  c('"<HEADER software/>"',
    'EwEVersion,"Atlantis output converted by atlantis_to_ewe_funcs.R"',
    paste0('Date,"', format(Sys.time(), "%d/%m/%Y %H:%M:%S"), '"'),
    '"<HEADER ecopath/>"',
    paste0('ModelName,"', model_name, '"'),
    paste0("ModelSource,", model_source),
    '"<HEADER ecosim/>"',
    paste0('EcosimScenario,"', scenario_name, '"'),
    "TimeSeries,-",
    paste0("StartYear,", start_year),
    '"<HEADER end/>"',
    "",
    paste0("Data,", data_label),
    "")
}

# Write a wide format Ecosim annual file (year rows, group columns)
##' @title Write EwE wide format annual csv
##' @author Beth Fulton
ewe_write_wide <- function(filename, data_label, values, years, header_info) {
  header <- do.call(ewe_header, c(list(data_label = data_label), header_info))
  col_header <- paste0("year\\group,", paste(seq_len(ncol(values)), collapse = ","))
  body <- sapply(seq_len(nrow(values)), function(k) paste(c(years[k], ewe_num(values[k, ])), collapse = ","))
  writeLines(c(header, col_header, body), filename)
}

# Write a long format Ecosim annual file (year, fleet, group, value)
##' @title Write EwE long format annual csv
##' @author Beth Fulton
ewe_write_long <- function(filename, data_label, values, years, header_info) {
  header <- do.call(ewe_header, c(list(data_label = data_label), header_info))
  nfleet <- dim(values)[2]
  nsp <- dim(values)[3]
  body <- character(0)
  for (k in seq_along(years)) {
    for (f in 1:nfleet) {
      body <- c(body, paste(years[k], f, 1:nsp, ewe_num(values[k, f, ]), sep = ","))
    }
  }
  writeLines(c(header, "year,fleet,group,value", body), filename)
}

# Write an Ecopath style matrix (prey rows x predator columns) as saved by the "Save to csv" button
##' @title Write EwE Ecopath style consumption or diet matrix
##' @author Beth Fulton
ewe_write_matrix <- function(filename, mat, group_names, is_diet = FALSE) {
  nsp <- length(group_names)
  fmt <- function(x) {
    out <- formatC(x, digits = 10, format = "f")
    out[x == 0] <- ""   # EwE leaves empty cells where there is no flow
    out
  }

  # EwE leaves out the columns of groups that do not consume anything. clean_up_network_files() in
  # funcs.R expects some columns to be missing and appends them back on at the end (and do_Gao()
  # then uses the column position as the group ID), so only the trailing non-consumer columns
  # (typically the detritus groups) are left out here
  ncol_keep <- nsp
  while (ncol_keep > 1 && sum(mat[, ncol_keep]) <= 0) ncol_keep <- ncol_keep - 1
  if (ncol_keep == nsp) {
    warning(paste("The last group in the matrix is a consumer, so no columns can be left out of", filename,
                  "- clean_up_network_files() in funcs.R needs at least one trailing non-consumer column"))
  }
  mat <- mat[, 1:ncol_keep, drop = FALSE]

  lines <- paste0(',"Prey \\ predator",', paste(1:ncol_keep, collapse = ","))
  for (i in 1:nsp) {
    lines <- c(lines, paste0(i, ',"', group_names[i], '",', paste(fmt(mat[i, ]), collapse = ",")))
  }
  import_row <- formatC(rep(0, ncol_keep), digits = 10, format = "f")
  sum_row <- formatC(colSums(mat), digits = 10, format = "f")
  if (is_diet) {
    lines <- c(lines,
               paste0(',Import,', paste(import_row, collapse = ",")),
               paste0(',Sum,', paste(sum_row, collapse = ",")),
               paste0(',"(1 - Sum)",', paste(formatC(ifelse(colSums(mat) > 0, 1 - colSums(mat), 0), digits = 10, format = "f"), collapse = ",")))
  } else {
    lines <- c(lines,
               paste0(nsp + 1, ',Import,', paste(import_row, collapse = ",")),
               paste0(nsp + 2, ',Sum,', paste(sum_row, collapse = ",")))
  }
  writeLines(lines, filename)
}

########################### HELPER (ID) FILE FUNCTIONS ###########################

# Create a template Species_ID.csv from the Atlantis group definitions. The values are first guesses
# and need to be checked by the user (see "Calculating Hub index and ETI from Ecopath with Ecosim.docx")
##' @title Create template species ID dataframe
##' @return Dataframe in the Species_ID.csv format
##' @author Beth Fulton
atl_make_species_id <- function(groups, bio_prm_lines, init_nc_file, diet_ref, pisc_threshold = 0.5,
                                pisczoopl_groups = NULL) {
  nsp <- nrow(groups)
  gt <- groups$GroupType
  is_tep <- gt %in% c("MAMMAL", "BIRD", "REPTILE", "SHARK")
  is_fish <- gt %in% c("FISH")
  is_jelly <- grepl("jell|gelat|medus|salp", groups$Name, ignore.case = TRUE)
  is_zoo <- gt %in% c("SM_ZOO", "MED_ZOO", "LG_ZOO")
  is_bact <- gt %in% c("PL_BACT", "SED_BACT", "ICE_BACT")
  is_phyto <- groups$IsProducer > 0 & !(gt %in% c("PHYTOBEN", "SEAGRASS", "MANGROVES", "SEAWEED", "TURF"))
  is_macro <- groups$IsProducer > 0 & !is_phyto
  is_habitat <- gt %in% c("SED_EP_FF", "CORAL", "SPONGE") | is_macro

  # Pelagic vs Demersal for fish from flagdem in the biology prm file
  flagdem <- sapply(groups$AtlCode, function(cd) atl_prm_value(bio_prm_lines, paste0("flagdem", cd), 1))

  pel_dem <- rep("Invert", nsp)
  pel_dem[is_fish] <- ifelse(flagdem[is_fish] > 0, "Demersal", "Pelagic")
  pel_dem[is_tep] <- "TEP"
  pel_dem[is_zoo | is_bact] <- "Micro"
  pel_dem[is_jelly] <- "Jelly"
  pel_dem[is_phyto] <- "Phytoplankton"
  pel_dem[is_macro] <- "Macroalgae"
  pel_dem[groups$IsDetritus > 0] <- "Detritus"

  # Piscivorous vs Planktivorous for fish - based on the proportion of fish and squid in the diet
  is_fish_prey <- gt %in% c("FISH", "SHARK", "CEP") |
    (gt == "FISH_INVERT" & grepl("squid|ceph|octo|cuttle", groups$Name, ignore.case = TRUE))
  diet_prop <- sweep(diet_ref, 2, pmax(colSums(diet_ref), 1e-30), "/")
  fish_in_diet <- colSums(diet_prop[is_fish_prey, , drop = FALSE])
  pisc_zoo <- pel_dem
  pisc_zoo[is_fish] <- ifelse(fish_in_diet[is_fish] >= pisc_threshold, "Piscivorous", "Planktivorous")
  # plot_biomass_ratios() in funcs.R fails if there are no Piscivorous (or no Planktivorous) groups,
  # so if the threshold puts all the fish in one class split them at the median instead
  if (sum(is_fish) > 1 && length(unique(pisc_zoo[is_fish])) < 2) {
    pisc_zoo[is_fish] <- ifelse(fish_in_diet[is_fish] > median(fish_in_diet[is_fish]), "Piscivorous", "Planktivorous")
    warning(paste("No fish had a proportion of fish and squid in their diet on the other side of pisc_threshold =",
                  pisc_threshold, "so fish have been split into Piscivorous and Planktivorous at the median proportion",
                  "- please check PISCZOOPL in the species ID file"))
  }

  # User assigned Piscivorous / Planktivorous groups (by Atlantis code or name) override the diet based values
  if (length(pisczoopl_groups) > 0) {
    allowed <- c("Piscivorous", "Planktivorous")
    user_class <- allowed[match(tolower(pisczoopl_groups), tolower(allowed))]
    if (any(is.na(user_class))) {
      stop(paste("pisczoopl_groups values must be Piscivorous or Planktivorous, not:",
                 paste(unique(pisczoopl_groups[is.na(user_class)]), collapse = " ")))
    }
    # Names can be the output group (e.g. a stanza) or the Atlantis group (which then covers all its stanzas)
    user_index <- lapply(names(pisczoopl_groups), function(nm) {
      idx <- which(groups$Code == nm | groups$Name == nm)
      if (length(idx) == 0) idx <- which(groups$AtlCode == nm | groups$AtlName == nm)
      idx
    })
    if (any(lengths(user_index) == 0)) {
      stop(paste("pisczoopl_groups names not found as group codes or names:",
                 paste(names(pisczoopl_groups)[lengths(user_index) == 0], collapse = " ")))
    }
    user_class <- rep(user_class, lengths(user_index))
    user_index <- unlist(user_index)
    not_fish <- user_index[!is_fish[user_index]]
    if (length(not_fish) > 0) {
      message(paste("Note pisczoopl_groups includes groups that are not fish:", paste(groups$Code[not_fish], collapse = " "),
                    "- they will be included in the fish based indicators (e.g. LFI, piscivore:planktivore ratio)"))
    }
    pisc_zoo[user_index] <- user_class
  }
  if (!all(c("Piscivorous", "Planktivorous") %in% pisc_zoo)) {
    warning("There are no Piscivorous or no Planktivorous groups - plot_biomass_ratios() in funcs.R will fail")
  }

  # Maximum age - number of cohorts x years per cohort for age structured groups (the longevity of the species,
  # so the same for all the stanzas of a group)
  max_age <- groups$NumCohorts * groups$NumAgeClassSize
  default_age <- ifelse(is_zoo, 0.5, ifelse(is_bact, 0.1, ifelse(is_phyto, 0.3, ifelse(is_macro, 7,
                        ifelse(groups$IsDetritus > 0, 1, 5)))))
  max_age[groups$NumCohorts <= 1] <- default_age[groups$NumCohorts <= 1]

  # Maximum length - from the weight of the oldest cohort (of the stanza if split into stanzas) in the initial
  # conditions and the length-weight relationship (li_a, li_b) in the biology prm file
  default_len <- ifelse(is_zoo, 1, ifelse(is_bact, 0.001, ifelse(is_phyto, 0.01, ifelse(is_macro, 500,
                        ifelse(is_jelly, 10, ifelse(groups$IsDetritus > 0, 1, 10))))))
  linf <- default_len
  if (!is.null(init_nc_file) && file.exists(init_nc_file)) {
    X_CN <- atl_prm_value(bio_prm_lines, "X_CN", 5.7)
    k_wetdry <- atl_prm_value(bio_prm_lines, "k_wetdry", 20)
    nc <- nc_open(init_nc_file)
    for (g in which(groups$NumCohorts > 1)) {
      sn_name <- paste0(groups$AtlName[g], groups$LastCohort[g], "_StructN")
      rn_name <- paste0(groups$AtlName[g], groups$LastCohort[g], "_ResN")
      li_a <- atl_prm_value(bio_prm_lines, paste0("li_a_", groups$AtlCode[g]))
      li_b <- atl_prm_value(bio_prm_lines, paste0("li_b_", groups$AtlCode[g]))
      if (!(sn_name %in% names(nc$var)) || is.na(li_a) || is.na(li_b)) next
      get_val <- function(v) {
        x <- suppressWarnings(max(ncvar_get(nc, v), na.rm = TRUE))
        if (!is.finite(x) || x <= 0) x <- ncatt_get(nc, v, "_FillValue")$value
        x
      }
      wgt_g <- (get_val(sn_name) + get_val(rn_name)) * X_CN * k_wetdry / 1000
      if (is.finite(wgt_g) && wgt_g > 0) {
        linf[g] <- round((wgt_g / li_a)^(1 / li_b), 1)
      }
    }
    nc_close(nc)
  }

  # Classification - the fishing flags (IsFished, isTAC) are only used for fish, sharks and the
  # larger invertebrates, plankton and bacteria are always Robust and habitat formers always Habitat
  is_fishable <- is_fish | gt == "SHARK" |
    !(is_tep | is_zoo | is_bact | is_jelly | groups$IsProducer > 0 | groups$IsDetritus > 0 | gt %in% c("SM_INF", "LG_INF"))
  fished_class <- ifelse(groups$isTAC > 0, "Target", ifelse(groups$IsFished > 0, "Byproduct", NA))
  classification <- rep("Robust", nsp)
  classification[is_fish] <- "Bycatch"
  classification[is_tep] <- "Vulnerable"
  use_fished <- is_fishable & !is.na(fished_class)
  classification[use_fished] <- fished_class[use_fished]
  classification[is_habitat] <- "Habitat"
  classification[groups$IsDetritus > 0] <- "Detritus"

  df_species <- data.frame(ID = groups$ID,
                           Group = groups$Name,
                           PelDemID = pel_dem,
                           PISCZOOPL = pisc_zoo,
                           Linf_cm = linf,
                           MaxAge_yr = max_age,
                           IsDetritus = groups$IsDetritus,
                           Classification = classification,
                           stringsAsFactors = FALSE)
  return(df_species)
}

# Check the group names will be matched correctly when clean_up_network_files() in funcs.R reads the
# consumption and diet matrices (it cleans the names and then fuzzy matches them to Species_ID.csv)
##' @title Check group names survive the name matching in funcs.R
##' @author Beth Fulton
atl_check_name_matching <- function(group_names) {
  clean_str_local <- function(x) {
    x <- tolower(x)
    x <- gsub("[<>]", "", x)
    x <- gsub("[^a-z0-9\\s]", "", x, perl = TRUE)
    x <- gsub("[.\\s]+", "_", x, perl = TRUE)
    stringr::str_trim(x)
  }
  species_clean <- tolower(stringr::str_trim(gsub("[.\\s]+", " ", group_names, perl = TRUE)))
  matched <- sapply(clean_str_local(group_names), function(nm) {
    exact <- which(species_clean == nm)
    if (length(exact) > 0) return(group_names[exact[1]])
    group_names[which.min(stringdist::stringdist(nm, species_clean, method = "jw"))]
  })
  bad <- which(matched != group_names)
  if (length(bad) > 0) {
    warning(paste("These group names will be mismatched by clean_up_network_files() - rename them in",
                  "Species_ID.csv and the consumption and diet files:",
                  paste(group_names[bad], "->", matched[bad], collapse = "; ")))
  }
  return(length(bad) == 0)
}

########################### MAIN CONVERSION FUNCTION ###########################

# Convert one Atlantis run into the EwE style files needed by EBFM_indicator_calculation_template.qmd
##' @title Convert Atlantis output to EwE style output
##' @param atlantis_out_dir Directory holding the Atlantis output files
##' @param output_prefix The output file name prefix used in the Atlantis run (e.g. "NESPout")
##' @param groups_file Atlantis functional groups csv file
##' @param bgm_file Atlantis bgm (box geometry) file
##' @param bio_prm_file Atlantis biology prm file
##' @param ewe_out_dir Where to write the Ecosim style files (i.e. file.path(InputDir, InputsubDir) in the qmd)
##' @param fisheries_file Atlantis fisheries csv file (optional)
##' @param init_nc_file Atlantis initial conditions netcdf file (optional - used to get Linf for Species_ID.csv)
##' @param network_dir Where to write the Ecopath style consumption and diet matrices and the ID files (NULL to skip)
##' @param start_year Calendar year of the first model year
##' @param area_km2 Model area (km2) - if NULL it is the area of the non-boundary boxes in the bgm file
##' @param ecopath_year_row Year (row) of the run used for the Ecopath style consumption and diet matrices
##' @param pisc_threshold Proportion of fish (and squid) in the diet at or above which a fish is classed Piscivorous
##' @param pisczoopl_groups Optional named vector assigning groups (names are Atlantis codes or names) as "Piscivorous"
##'   or "Planktivorous" in the species ID file, overriding the diet based classification, e.g.
##'   c(FVT = "Piscivorous", Mackerel = "Planktivorous")
##' @param catch_source Which fisheries output to use: "auto" (first found of the order below), "txt" (CatchPerFishery.txt
##'   and DiscardPerFishery.txt), "nc" (CATCH.nc) or "total" (Catch.txt and Discard.txt as a single fleet)
##' @param stanzas If TRUE split the age structured groups into stanzas (see the notes at the top of this file).
##'   Use the same stanza settings for the scenario and unfished runs so they have the same groups
##' @param stanza_groups Codes of the groups to split into stanzas (NULL for all the age structured groups)
##' @param stanza_breaks Optional named list (by group code) of the first cohort of each stanza after the first,
##'   e.g. list(FTS = c(3, 6)). Groups not in the list are split into juveniles and adults (using <Code>_age_mat)
##' @param stanza_labels Optional named list (by group code) of the stanza labels, e.g. list(FTS = c("larval", "juv", "adult"))
##' @return List with the model area, years written and the file names
##' @author Beth Fulton
convert_atlantis_to_ewe <- function(atlantis_out_dir, output_prefix, groups_file, bgm_file, bio_prm_file,
                                    ewe_out_dir, fisheries_file = NULL, init_nc_file = NULL,
                                    network_dir = NULL, model_name = "Atlantis model",
                                    scenario_name = output_prefix, start_year = 1, area_km2 = NULL,
                                    ecopath_year_row = 1,
                                    consumption_filename = "Atlantis_Consumption.csv",
                                    diet_filename = "Atlantis_Diet_composition.csv",
                                    species_ID_filename = "Species_ID.csv",
                                    fleet_ID_filename = "Fleet_ID.csv",
                                    overwrite_id_files = FALSE, pisc_threshold = 0.5, pisczoopl_groups = NULL,
                                    include_m0 = TRUE,
                                    catch_source = "auto",
                                    stanzas = FALSE, stanza_groups = NULL, stanza_breaks = NULL,
                                    stanza_labels = NULL) {

  atl_load_libraries()
  out_file <- function(suffix) file.path(atlantis_out_dir, paste0(output_prefix, suffix))
  for (needed in c("BiomIndx.txt", "DietCheck.txt", "PROD.nc", ".nc")) {
    if (!file.exists(out_file(needed))) stop(paste("Missing Atlantis output file:", out_file(needed)))
  }

  ## Model set up ----
  groups <- atl_read_groups(groups_file)
  boxes <- atl_read_bgm(bgm_file)
  bio_prm_lines <- readLines(bio_prm_file, warn = FALSE)
  X_CN <- atl_prm_value(bio_prm_lines, "X_CN", 5.7)
  k_wetdry <- atl_prm_value(bio_prm_lines, "k_wetdry", 20)
  mgN_to_tonnes <- X_CN * k_wetdry / 1e9  # mg N -> mg dry -> mg wet -> tonnes wet

  if (is.null(area_km2)) {
    area_km2 <- sum(boxes$area[!boxes$is_boundary]) / 1e6
  }
  message(paste("Model area (non-boundary boxes):", round(area_km2, 1), "km2 - use this for area_model in the qmd"))

  # Number of complete model years in the run
  bi_time <- read.table(out_file("BiomIndx.txt"), header = TRUE, check.names = FALSE)$Time
  nyears <- floor(max(bi_time) / 365 + 1e-6)
  if (nyears < 1) stop("Less than a full year of Atlantis output")
  years <- start_year + 0:(nyears - 1)
  message(paste("Converting", nyears, "complete model years:", years[1], "-", years[nyears]))

  ## Load the Atlantis output ----
  bio_t <- atl_annual_biomass(out_file("BiomIndx.txt"), groups, nyears)
  vol_info <- atl_box_volumes(out_file(".nc"), boxes)
  cons <- atl_annual_consumption(out_file("PROD.nc"), vol_info, groups, nyears, mgN_to_tonnes, by_box = stanzas)
  diet <- atl_annual_diet(out_file("DietCheck.txt"), groups, nyears)
  fleets <- atl_read_fleets(fisheries_file)
  fishery <- atl_annual_catch(atlantis_out_dir, output_prefix, fleets, groups, boxes, nyears, catch_source)

  stanza <- NULL
  if (stanzas) {
    stanza <- atl_setup_stanzas(out_file(".nc"), groups, boxes, bio_prm_lines, nyears, mgN_to_tonnes, cons$box, bio_t,
                                stanza_groups, stanza_breaks, stanza_labels)
    cons <- cons$total
  }
  Q <- atl_consumption_matrix(cons, diet, groups, bio_t, stanza)

  if (stanzas) {
    if (length(stanza$log$no_window) > 0) {
      message(paste("For", length(unique(stanza$log$no_window)), "predator cohort - prey pairs none of the prey was in the",
                    "feeding window at the annual mean sizes, so their consumption was shared across the prey stanzas by biomass"))
    }
    # From here on the groups are the stanza (output) groups
    fishery <- atl_stanza_catch(fishery, out_file("CATCH.nc"), stanza, nyears, mgN_to_tonnes)
    bio_t <- stanza$biomass
    groups <- stanza$groups
  }
  m0 <- atl_annual_m0(if (include_m0) out_file("SpecificMort.txt") else NULL, groups, nyears, stanza$cohort$bio)

  ## Derived quantities ----
  bio_safe <- pmax(bio_t, 1e-20 * area_km2)
  catch_fg <- fishery$landings + fishery$discards           # years x fleet x group (tonnes)
  catch_t <- apply(catch_fg, c(1, 3), sum)                  # years x group (tonnes)
  if (is.null(dim(catch_t))) catch_t <- matrix(catch_t, nrow = nyears)
  colnames(catch_t) <- groups$Code

  cons_by_pred <- apply(Q, c(1, 3), sum)                    # years x group (tonnes)
  eaten <- apply(Q, c(1, 2), sum)                           # years x group (tonnes)
  QB <- cons_by_pred / bio_safe
  M2 <- eaten / bio_safe
  F_tot <- catch_t / bio_safe
  Z <- m0 + M2 + F_tot
  Z <- pmax(Z, 1e-6)  # keep production > 0 so the Green Band log regression works

  F_fleet <- catch_fg
  for (f in seq_len(dim(catch_fg)[2])) F_fleet[, f, ] <- catch_fg[, f, ] / bio_safe

  tl <- atl_trophic_levels(Q, groups)

  ## Write the Ecosim style files ----
  dir.create(ewe_out_dir, recursive = TRUE, showWarnings = FALSE)
  header_info <- list(model_name = model_name, model_source = normalizePath(atlantis_out_dir),
                      scenario_name = scenario_name, start_year = start_year)

  ewe_write_wide(file.path(ewe_out_dir, "biomass_annual.csv"), "Biomass", pmax(bio_t / area_km2, 1e-20), years, header_info)
  ewe_write_wide(file.path(ewe_out_dir, "catch_annual.csv"), "Catch", catch_t / area_km2, years, header_info)
  ewe_write_wide(file.path(ewe_out_dir, "consumption-biomass_annual.csv"), "ConsumptionBiomass", QB, years, header_info)
  ewe_write_wide(file.path(ewe_out_dir, "mortality_annual.csv"), "Mortality", Z, years, header_info)
  ewe_write_wide(file.path(ewe_out_dir, "tl_annual.csv"), "TL", tl, years, header_info)
  ewe_write_long(file.path(ewe_out_dir, "landings_annual.csv"), "Landings", fishery$landings / area_km2, years, header_info)
  ewe_write_long(file.path(ewe_out_dir, "catch-fleet-group_annual.csv"), "CatchFleetGroup", catch_fg / area_km2, years, header_info)
  ewe_write_long(file.path(ewe_out_dir, "mort-fleet-group_annual.csv"), "MortFleetGroup", F_fleet, years, header_info)
  message(paste("Ecosim style files written to", ewe_out_dir))

  ## Write the Ecopath style matrices and the helper ID files ----
  files_written <- list(ewe_out_dir = ewe_out_dir)
  if (!is.null(network_dir)) {
    dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)
    Q_ref <- Q[ecopath_year_row, , ]
    # Use the all years average for any predator that did not feed in the reference year
    Q_all <- apply(Q, c(2, 3), mean)
    no_feed <- which(colSums(Q_ref) <= 0 & colSums(Q_all) > 0)
    Q_ref[, no_feed] <- Q_all[, no_feed]
    diet_ref <- sweep(Q_ref, 2, pmax(colSums(Q_ref), 1e-30), "/")

    ewe_write_matrix(file.path(network_dir, consumption_filename), Q_ref / area_km2, groups$Name)
    ewe_write_matrix(file.path(network_dir, diet_filename), diet_ref, groups$Name, is_diet = TRUE)

    # Template ID files - not overwritten unless asked as the user will have edited them
    species_file <- file.path(network_dir, species_ID_filename)
    fleet_file <- file.path(network_dir, fleet_ID_filename)
    if (file.exists(species_file) && !overwrite_id_files) {
      species_file <- file.path(network_dir, paste0("atlantis_template_", species_ID_filename))
    }
    if (file.exists(fleet_file) && !overwrite_id_files) {
      fleet_file <- file.path(network_dir, paste0("atlantis_template_", fleet_ID_filename))
    }
    # The diet over the whole run is used to class fish as Piscivorous or Planktivorous
    df_species <- atl_make_species_id(groups, bio_prm_lines, init_nc_file, apply(Q, c(2, 3), sum), pisc_threshold,
                                      pisczoopl_groups)
    write.csv(df_species, species_file, row.names = FALSE)
    df_fleet_out <- data.frame(ID = fishery$fleets$ID, Fleet.name = fishery$fleets$Code)
    names(df_fleet_out)[2] <- "Fleet name"
    write.csv(df_fleet_out, fleet_file, row.names = FALSE)
    atl_check_name_matching(groups$Name)

    message(paste("Consumption and diet matrices (from year", years[ecopath_year_row], ") written to", network_dir))
    message(paste("Template species ID file written to", species_file, "- please check the PelDemID, PISCZOOPL, Linf_cm, MaxAge_yr and Classification entries"))
    message(paste("Fleet ID file written to", fleet_file))
    files_written <- c(files_written, list(consumption = file.path(network_dir, consumption_filename),
                                           diet = file.path(network_dir, diet_filename),
                                           species_ID = species_file, fleet_ID = fleet_file))
  }

  return(invisible(list(area_model = area_km2, years = years, nyears = nyears, groups = groups, stanzas = stanzas,
                        fleets = fishery$fleets, files = files_written,
                        biomass_t = bio_t, consumption_t = Q, catch_t = catch_t, QB = QB, Z = Z, M0 = m0,
                        M2 = M2, F = F_tot, TL = tl)))
}
