## Functions used when running the EBFM_indicator_calculation.qmd R code

########################### LIBRARY FUNCTIONS ###########################
# Helpful function for library install and load
# first, checks if libraries are installed
# if not, the library is installed prior to loading
##' @title Install-load function. This checks if a library is installed before attempting to load it
##' @return Installed and loaded libraries
##' @author Beth Fulton
install_load <- function (this_package)  {
  package <- c(this_package)
  if(package %in% rownames(installed.packages()))
    do.call('library', list(package))
  ## if package is not installed locally, download, then load
  else {
    install.packages(package, dependencies=TRUE)
    do.call("library", list(package))
  }
}

# Function to load all the needed libraries 
# This calls install_load to check if a library is installed before attempting to load it
##' @title Load libraries function. 
##' @author Beth Fulton
load_libraries <- function() {
  # reticulate::py_install("kaleido", method = "auto", pip = TRUE)
  
  ## Make sure quarto is in place
  install_load("quarto")
  
  ## String handling
  install_load("stringr")
  install_load("stringdist")
  
  ### Data handling and table reshaping ----
  install_load("tidyverse")
  install_load("reshape2")
  install_load("devtools")
  install_load("MASS")
  install_load("data.table")
  install_load("dplyr")
  
  ### Plotting ----
  install_load("ggplot2")
  install_load("RColorBrewer")
  install_load("plot3D")
  install_load("plotly")
  install_load("webshot")
  install_load("htmlwidgets")
  
  ### PCA and Clustering ----
  install_load("factoextra")
  install_load("FactoMineR")
  install_load("corrplot")
  install_load("ape")
  install_load("plot3D")
  install_load("heatmaply")
  install_load("data.table")
  
  ### Timeseries analysis ----
  install_load("zoo")
  
  # Other ----
  install_load("Rcpp")
  install_load("inline")
  install_load("vegan")
  install_load("patchwork")
  install_load("cowplot")
  install_load("ggpubr")
  install_load("matrixStats")
  install_load("drc")
  install_load("viridis")
  install_load("gridExtra")
  install_load("cowplot")
  
  # Network related ----
  install_load("data.table")
  install_load("magrittr")
  install_load("ggraph")
  install_load("qgraph")
  install_load("network")
  install_load("magrittr")
  install_load("stringr")
  install_load("visNetwork")
  install_load("intergraph")
  install_load("tidygraph")
  install_load("igraph") 

}


# Clean strings to remove special characters - need to do this for the network analysis
##' @title Clean strings to remove special characters
##' @author Beth Fulton
clean_str <- function(x) {
  x <- tolower(x)
  x <- gsub("[<>]", "", x)
  x <- gsub("[^a-z0-9\\s]", "", x, perl = TRUE)  # keep only letters, numbers, spaces
  x <- gsub("[.\\s]+", "_", x, perl = TRUE)      # replace whitespace with _
  x <- str_trim(x)
  x
}

# Function: clean the consumption and diet matrix names, match against cleaned master,
# but return the original unmodified master name from species_ID.csv
##' @title Finding closest match of column or row name from consumption and diet data with master list in species_ID.csv
##' @author Beth Fulton
match_to_master <- function(names_to_match, species_names, species_clean) {
  names_clean <- clean_str(names_to_match)
  
  matched <- sapply(names_clean, function(nm) {
    # Try exact match first
    exact <- which(species_clean == nm)
    if (length(exact) > 0) return(species_names[exact[1]])
    
    # Fuzzy fallback
    dists <- stringdist(nm, species_clean, method = "jw")
    species_names[which.min(dists)]
  })
  
  return(matched)
}

########################### DIRECTORY FUNCTIONS ###########################

# Create output directory if required
##' @title define directories - checks to see if all output directories exists and creates them if not
##' @author Beth Fulton
define_directories <-function(mainDir, OutsubDirStep1, OutsubDirStep2) {
  OutDir <- paste(mainDir,"/",OutsubDirStep1,"/", OutsubDirStep2,"/", sep="")
  
  ifelse(!dir.exists(file.path(mainDir, OutsubDirStep1)), dir.create(file.path(mainDir, OutsubDirStep1)), FALSE)
  ifelse(!dir.exists(file.path(mainDir, OutsubDirStep2)), dir.create(file.path(mainDir, OutsubDirStep2)), FALSE)
  if(!dir.exists(file.path(OutDir))) dir.create(file.path(OutDir))  
  
  return(OutDir)
}

########################### HUB SPECIES RELATED FUNCTIONS ###########################

clean_up_network_files <- function(consumption_filename, diet_filename, species_ID_filename) {
  
  ## Load the data ----
  # consumption data ----
  df_consumption <- read.table(consumption_filename,
                         header = TRUE, 
                         sep = ",")
  
  # diet composition data ----
  df_diet <- read.table(diet_filename,
                         header = TRUE, 
                         sep = ",")
  
  ## functional groups in numeric order with a number assigned
  df_species <- read.table(species_ID_filename,
                           header = TRUE, 
                           sep=",")
  
  ## Now clean up the names so they work in the network analysis
  # Master species names — left completely untouched from what is in species_ID.csv
  species_names <- df_species$Group
  
  # Cleaned version of master names just for comparison (so can get best match possible)
  species_clean <- tolower(str_trim(gsub("[.\\s]+", " ", species_names, perl = TRUE)))
  
  #### Clean up the consumption dataframe
  # Remove unneeded first column
  df_consumption <- df_consumption %>% select(-X)
  # Rename the column of species names to Group
  names(df_consumption)[names(df_consumption) == "Prey...predator"] <- "Group"
  # Clean up the names of the species to remove any special characters
  df_consumption$Group <- clean_str(df_consumption$Group)
  # Remove unneeded rows
  do_not_keep <- c("import","sum")
  df_consumption <- df_consumption %>% filter(!Group %in% do_not_keep)
  # Replace the empty cells
  df_consumption[is.na(df_consumption)] <- 0.0
  # Replace the group names with the names in Species_ID.csv
  df_consumption$Group <- match_to_master(df_consumption$Group, species_names, species_clean)
  # Find missing columns (not a square matrix yet so find what need to do to fix it)
  col_nums_present <- as.integer(sub("X", "", grep("^X\\d+$", names(df_consumption), value = TRUE)))
  all_nums <- 1:max(as.integer(rownames(df_consumption)))
  missing_nums <- setdiff(all_nums, col_nums_present)
  # Add zero columns for missing species
  missing_cols <- setNames(
    as.data.frame(matrix(0, nrow = nrow(df_consumption), ncol = length(missing_nums))),
    paste0("X", missing_nums)
  )
  df_consumption <- cbind(df_consumption, missing_cols)
  # Rename all X columns using row number to match Group
  x_cols <- grep("^X\\d+$", names(df_consumption))
  col_nums <- as.integer(sub("X", "", names(df_consumption)[x_cols]))
  names(df_consumption)[x_cols] <- df_consumption$Group[match(col_nums, as.integer(rownames(df_consumption)))]
  rownames(df_consumption) <- df_consumption$Group
  
  #### Clean up the diet dataframe
  # Remove unneeded first column
  df_diet <- df_diet %>% select(-X)
  # Rename the column of species names to Group
  names(df_diet)[names(df_diet) == "Prey...predator"] <- "Group"
  # Clean up the names of the species to remove any special characters
  df_diet$Group <- clean_str(df_diet$Group)
  # Remove unneeded rows
  do_not_keep <- c("import","sum","1_sum")
  df_diet <- df_diet %>% filter(!Group %in% do_not_keep)
  # Replace the empty cells
  df_diet[is.na(df_diet)] <- 0.0
  # Replace the group names with the names in Species_ID.csv
  df_diet$Group <- match_to_master(df_diet$Group, species_names, species_clean)
  # Find missing columns (not a square matrix yet so find what need to do to fix it)
  col_nums_present <- as.integer(sub("X", "", grep("^X\\d+$", names(df_diet), value = TRUE)))
  all_nums <- 1:max(as.integer(rownames(df_diet)))
  missing_nums <- setdiff(all_nums, col_nums_present)
  # Add zero columns for missing species
  missing_cols <- setNames(
    as.data.frame(matrix(0, nrow = nrow(df_diet), ncol = length(missing_nums))),
    paste0("X", missing_nums)
  )
  df_diet <- cbind(df_diet, missing_cols)
  # Rename all X columns using row number to match Group
  x_cols <- grep("^X\\d+$", names(df_diet))
  col_nums <- as.integer(sub("X", "", names(df_diet)[x_cols]))
  names(df_diet)[x_cols] <- df_diet$Group[match(col_nums, as.integer(rownames(df_diet)))]
  rownames(df_diet) <- df_diet$Group
  
  return(list(df_consumption = df_consumption, df_diet = df_diet, df_species = df_species))
}
  
# All the steps need to calculate the Hub species and update the Species_ID.csv accordingly
##' @title Hub species function - undertakes a networks analysis to identify hub species
##' @param file Uses consumption and diet matrices from an Ecopath model to calculate the hub_species
##' @author Beth Fulton
hub_species <- function(df_consumption, df_diet, df_species, species_ID_filename, OutDir, use_weights, detailed_check = 0, use_top_5_pct_only = 1) {

  ## Load the data ----
  # consumption data ----
  consData <- df_consumption %>% select(-Group)

  # diet composition data ----
  dietData <- df_diet %>% select(-Group)
  
  # trim the species info
  cols_to_remove <- c("PelDemID", "PISCZOOPL", "Linf_cm", "MaxAge_yr", "IsDetritus", "Classification")
  localIDS <- df_species %>% select(-all_of(cols_to_remove))
  
  ###### CREATE NETWORK STRUCTURE ----
  ## Create edges df ----
  ## Calculate Degree in and Degree out ----
  step1 <- reshape2::melt(as.matrix(consData))
  step1a <- reshape2::melt(as.matrix(dietData))
  names(step1)[names(step1) == "Var1"] <- "preyName"
  names(step1)[names(step1) == "Var2"] <- "predName"
  names(step1)[names(step1) == "value"] <- "propDiet"
  step1$predName <- gsub(".", " ", step1$predName, fixed=TRUE)  #Sort the issue with spaces in names
  
  # Replace names with IDs
  step2a <- step1 %>% dplyr::inner_join(localIDS,by=c("preyName" = "Group"))
  names(step2a)[names(step2a) == "ID"] <- "prey_ID"
  
  step2 <- step2a %>% dplyr::inner_join(localIDS,by=c("predName" = "Group"))
  names(step2)[names(step2) == "ID"] <- "pred_ID"
  
  nDiet <- step2[step2$propDiet != 0, ]
  
  # Identify all the links between the functional groups
  if (use_weights > 0) {
    nl <- subset(nDiet, select = -c(preyName, predName))
    names(nl)[names(nl) == "prey_ID"] <- "V1"
    names(nl)[names(nl) == "pred_ID"] <- "V2"
    nl$na <- "FALSE"
    nl$V1 <- as.numeric(nl$V1)
    nl$V2 <- as.numeric(nl$V2)
    nl$na <- as.logical.factor(nl$na)
    nl$na <- FALSE
    nl$weight <- nl$propDiet
    nl <- subset(nl, select = -c(propDiet))
  } else {
    nl <- subset(nDiet, select = -c(preyName, predName, propDiet))
    names(nl)[names(nl) == "prey_ID"] <- "V1"
    names(nl)[names(nl) == "pred_ID"] <- "V2"
    nl$na <- "FALSE"
    nl$V1 <- as.numeric(nl$V1)
    nl$V2 <- as.numeric(nl$V2)
    nl$na <- as.logical.factor(nl$na)
    nl$na <- FALSE
  }
  
  nlv <- localIDS
  nlv$na <- "FALSE"
  nlv$na <- as.logical.factor(nlv$na)
  nlv$na <- FALSE
  nlv$vertex.names <- nlv$Group
  names(nlv)[names(nlv) == "species_ID"] <- "intergraph_id"
  nlv <- subset(nlv, select = -c(Group))
  
  g <- asIgraph(nl, vertices=nlv, directed=TRUE)
  set_vertex_attr(g, "label", value=nlv[,2])
  
  
  # IGRAPH WORK UP OF THE NETWORK INFORMATION COLLATED ABOVE ----
  edges <- as_edgelist(g)   # Can also be extracted using E(g)
  nodes <- V(g)$vertex.names    # Name of vertices
  E(g)$weight  # Use this to see weight per edge (NULL here as incidence only)
  
  # View Adjacency matrix - usually used if debugging things so commented out for now
  # g[]
  # To see only the first line use g[1,]
  
  ## Plot the network ---
  OutFilename <- paste(OutDir,"/Network_plot.png",sep="")
  png(OutFilename, 1200, 800)
  plot(g)
  dev.off()
  
  #### OPTIONAL NETWORK ANALYSIS - WILL BE REPEATED AND SAVED BELOW
  if (detailed_check > 0) {
  
    ## Use tidygraph instead - create the network (graph) - can also use as_tbl_graph() to convert an igraph or netwprk 
    # library network to a tidygraph graph
    food_tidy <- as_tbl_graph(g, directed = TRUE)
    
    ## Now do real network analysis
    centRes <- centrality(food_tidy)
    dnonorm <- degree(g, mode = "all", normalized = FALSE)
    
    # Node strength (degree):
    centRes$OutDegree # Or InDegree, it's the same in unweighted networks
    
    # Closeness:
    centRes$Closeness
    
    # Betweenness:
    centRes$Betweenness
    
    # Plotting centrality measures
    centralityPlot(food_tidy)
    
    ## Other statistical approaches to calculate the indicators ----
    # Number of nodes
    length(V(g))
    
    # Find the standalone components - subwebs
    clusters(g)   # Output shows its one interconnected cluster (community)
    # If it wasn't all one community but you wanted to focus on one community you can suck out just that sub-web using
    # subweb<-induced.subgraph(food_tidy, which(clusters(newnet)$membership == x)  where x is the mmerbship number of the subweb you want
    
    # Average path length - can also be done using mean_distance(g)
    average.path.length(g)
    
    # Clustering coefficient
    transitivity(g)
    
    # Centralization degree score
    centralization.degree(g)$centralization
  }
  
  #### CORE CALCULATIONS FOR HUB SPECIES
  
  # iGraph degree routine Degree centrality is simplest of the methods, 
  # it measures the number of connections (links) between a node and all other nodes. 
  deg_all <- degree(g, mode = "all", normalized = FALSE) # All links in our out
  deg_out <- degree(g, mode = "out", normalized = FALSE) # Only links where the group is prey
  
  # Page rank score - originally defined by google to identify high traffic 
  # flow web pages in their network which needed maintaining, would be high search hits
  pg_rank <- page_rank(g, damping = 0.85) 
  
  # Create a new dataframe from these network indices
  df_network <- data.frame(deg_all, deg_out, pg_rank$vector)
  names(df_network)[names(df_network) == "pg_rank.vector"] <- "PageRank"
  
  if (detailed_check > 0) {
    deg_in <- degree(g, mode = "in", normalized = FALSE) # Only links where the group is predator
    
    # Degrees but without normalising first
    centrality_all <- degree(g, mode = "all", normalized = TRUE)
    centrality_in <- degree(g, mode = "in", normalized = TRUE)
    centrality_out <- degree(g, mode = "out", normalized = TRUE)
  
    # Closeness centrality is an evaluation of the proximity of a node to all other
    # nodes in a network, not only the nodes to which it is directly connected
    closeness_in <- closeness(g, mode="in", weights=NA, normalized=T) 
    closeness_out <- closeness(g, mode="out", weights=NA, normalized=T)
    closeness_all <- closeness(g, mode="all", weights=NA, normalized=T)
    
    # Betweenness centrality looks for chokepoints in the network
    betwn <- betweenness(g, directed=F, weights=NA, normalized = T)
    
    # SURF----
    dims <- dim(dietData)
    preydim <- dims[[1]]
    preddim <- dims[[2]]
    dD <- as.matrix(dietData)
    pstep <- matrix(0, preydim, preddim)
    
    # Assumes each row is a prey and columns is predators
    # The columns are diet vector (proportion of diet for predator j from prey i)
    L <- 0
    for(i in 1:preydim) {
      for(j in 1:preddim) {
        pstep[i,j] <- dD[i,j] * dD[i,j] 
        if (dD[i,j] > 0) {
          L <- L + 1
        }
      }
    }
    
    df <- as.data.frame(pstep)
    dfrowsum <- rowSums(df)
    SURF <- dfrowsum / L

    # SUPPLY CHAIN CRITICALITY ----
    dims <- dim(consData)
    cpreydim <- dims[[1]]
    cpreddim <- dims[[2]]
    cD <- as.matrix(consData)
    dfcolsum <- colSums(consData)
    step1 <- matrix(0, cpreydim, cpreddim)
    step2 <- rep(0, cpreddim)
    step3 <- matrix(0, cpreydim, cpreddim)
    
    # Assumes each row is a prey and columns is predators
    # The columns are consumption vector (proportion of consumption for predator j from prey i)
    for(i in 1:cpreydim) {
      for(j in 1:cpreddim) {
        if(dfcolsum[j] > 0) {
          step1[i,j] <- cD[i,j] / dfcolsum[j]
        }
      }
    }
    
    totcons <- sum(dfcolsum, na.rm = TRUE)
    step2 <- dfcolsum / totcons
    
    for(i in 1:cpreydim) {
      for(j in 1:cpreddim) {
        step3[i,j] <- step1[i,j] * step2[j] * step2[j]
      }
    }
    
    step3df <- as.data.frame(step3)
    CRITpred <- colSums(step3df)
    CRITprey <- rowSums(step3df)
    
    # In the case where all the indices are desired then add those to the datframe too
    new_cols <- data.frame(deg_in, centrality_all, centrality_in, centrality_out, 
                           closeness_in, closeness_out, closeness_all,
                           betwn, SURF, CRITpred, CRITprey)
    df_network <- cbind(df_network, new_cols)
  }
  
  ## SAVE THE RESULTS
  # Update the species_ID.csv file with hub information
  df_species <- Update_Species_ID_File(df_species, df_network, species_ID_filename, use_top_5_pct_only)
  
  # Store Ooutput of the netowkr anaalyses
  Network_index_filename <- paste(OutDir,"/","NetworkAnalysis_Outputs.csv",sep="")
  write.csv(df_network, file = Network_index_filename, row.names = FALSE)
  
  return(df_species)
}

# Updates the Species_ID.csv to include the information on hub species 
# and other classifications needed to calculate tte Green Band and ETI indicators
##' @param file Uses network indices calculated in hub_species function
##' @author Beth Fulton

Update_Species_ID_File <- function(df_species, df_network, species_ID_filename, use_top_5_pct_only) {

  ## FOCUS ON LIVING GROUPS
  df_combined <- cbind(df_species, df_network)
  df_living <- df_combined %>% filter(IsDetritus != 1) # Assumes 0 if not Detritus, 1 if is Detritus
  
  ## CALCULATE THE RANKS FOR deg_all, deg_out and pg_rank
  df_living <- df_living %>% 
    mutate(rank_Degree = rank(-deg_all, ties.method = "min"),
           rank_DegreeOut = rank(-deg_out, ties.method = "min"),
           rank_PageRank = rank(-PageRank, ties.method = "min")) 
  
  # Now per row (group) get the minimum rank across degree, degree out and page rank
  df_living <- df_living %>% mutate(
    hub_score = pmin(rank_Degree,rank_DegreeOut,rank_PageRank))
    
  # Rank teh hub scores and take the top 5% (10% for moderately small models as otehrwise too few to be useful)
  threshold <- ceiling(0.05 * nrow(df_living))
  if(use_top_5_pct_only < 1) {
    if((threshold) < 4 && (nrow(df_living) > 20)) {
      # Use top 10 for moderately small models
      threshold <- ceiling(0.1 * nrow(df_living))
    }
  }  
  top5pc <- df_living %>%
    filter(hub_score <= threshold) %>%
    ungroup() %>%
    pull(Group)
  
  # Update the Species_ID dataframe and file
  df_species$isHub <- 0
  
  # Add new column identifying hub species based on the raw calculation
  df_species <- df_species %>%
    mutate( isHub = if_else( Group %in% top5pc, 1, isHub))
  
  # Filter out phytoplankton groups as hub species
  do_not_keep <- c("Micro","Phytoplankton")
  df_species <- df_species %>%
    mutate( isHub = if_else( PelDemID %in% do_not_keep, 0, isHub))
  
  # Update the Classification field accordingly
  df_species <- df_species %>%
    mutate( Classification = if_else(
        isHub == 1, "Hub", Classification)
    )
  
  # Create the Class_code column of the Dataframe
  df_species <- df_species %>%
    mutate(
      Class_Code = case_when(
        Classification == "Detritus"   ~ 0L, 
        Classification == "Vulnerable" ~ 1L, 
        Classification == "Habitat"    ~ 2L, 
        Classification == "Byproduct"  ~ 3L, 
        Classification == "Bycatch"    ~ 4L, 
        Classification == "Target"     ~ 5L, 
        Classification == "Robust"     ~ 6L, 
        Classification == "Hub"        ~ 7L,
        TRUE ~ NA_integer_  
      )
    )
  
  ## OUTPUT FINAL FILE
  updated_species_filename <- paste(OutDir,"/","updated_Species_ID.csv", sep="")
  write.csv(df_species, file = updated_species_filename, row.names = FALSE)
  
  return(df_species)
  
}

########################### EWE DATA LOAD FUNCTIONS ###########################

# Helpful function for loading EwE files
##' @title Skip line function. This automatically determines the number of rows based on the ewe version
##' @param file EwE .csv output file
##' @return Number of lines to skip
##' @author Javier Porobic
skip_lines <- function(file){
  strs  <- readLines(file)
  lines <- grep('year\\\\group,', strs) -1
  return(lines)
}

# Function to load in the Ecosim outputs as the dataframe for calculating indicators
##' @title Create . This automatically determines the number of rows based on the ewe version
##' @param file EwE .csv output files
##' @return List of the dataframes used to calculate the EBFM indicators
##' @author Beth Fulton
create_dataframe <- function(InputDir, InputsubDir, df_species, fleet_ID_filename, 
                             bo_in_same_file = 1, bo_yr_row = 1, bo_conservative = 1,
                             chl_from_file = 0, yr_ignore = 1900) {
  # Read in biomass annual data ----
  BioFile <- paste(InputDir,"/",InputsubDir,"/","biomass_annual.csv",sep="")
  d <- read.csv(BioFile, header = FALSE)
  s_year = d[c(1:7), ]
  
  # Read in consumption data ----
  ConsumFile <- paste(InputDir,"/",InputsubDir,"/","consumption-biomass_annual.csv",sep="")
  skip_this <- skip_lines(ConsumFile)
  consumption <- read.csv(
    ConsumFile,
    header = TRUE,
    skip = skip_this,
    check.names = FALSE
  )
  names(consumption)[1] <- "Year"
  
  # Read in catch data ----
  CatchFile <- paste(InputDir,"/",InputsubDir,"/","catch_annual.csv",sep="")
  skip_this <- skip_lines(CatchFile)
  catch <- read.csv(CatchFile, header = T, skip = skip_this, check.names = FALSE) # because csv file has headers - will need to skip the first 9 rows; check names makes sure an X isn't added before the number for the species
  names(catch)[1] <- 'Year'   # changes the name of the first column to year instead of year//group
  
  # Read aggregate landings data ----
  LandFile <- paste(InputDir,"/",InputsubDir,"/","landings_annual.csv",sep="")
  landings <- read.csv(LandFile, header = T, skip = skip_this, check.names = FALSE) # because csv file has headers - will need to skip the first 9 rows; check names makes sure an X isn't added before the number for the species
  names(landings)[1] <- 'Year'   # changes the name of the first column to year instead of year//group
  names(landings)[3] <- 'species_id'   # changes the name of the first column to year instead of year//group
  
  # Read biomass data ----
  biomass <- read.csv(BioFile, header = T, skip = skip_this, check.names = FALSE)
  names(biomass)[1] <- 'Year'
  
  # Read total catch and landings data per fleet ----
  RemFile <- paste(InputDir,"/",InputsubDir,"/","catch-fleet-group_annual.csv",sep="")
  removals_fleet <- read.csv(RemFile, header = T, skip = skip_this, check.names = FALSE)
  names(removals_fleet)[1] <- 'Year'
  LFFile <- paste(InputDir,"/",InputsubDir,"/","landings_annual.csv",sep="")
  landings_fleet <- read.csv(LFFile, header = T, skip = skip_this, check.names = FALSE)
  names(landings_fleet)[1] <- 'Year'
  
  # Read species id - put names ---- used to read file now uses df_species
  id <- df_species
  
  # Read fleet id - put names ----
  df_fleet <- read.csv(fleet_ID_filename) 
  names(df_fleet)[names(df_fleet) == "Fleet.name"] <- "FleetName"
  
  # Read mortality data ----
  MortFile <- paste(InputDir,"/",InputsubDir,"/","mort-fleet-group_annual.csv",sep="")
  mortality_fleet <- read.csv(MortFile, header = T, skip = skip_this, check.names = FALSE)
  names(mortality_fleet)[1] <- 'Year'
  MortFile <- paste(InputDir,"/",InputsubDir,"/","mortality_annual.csv",sep="")
  mortality <- read.csv(MortFile, header = T, skip = skip_this, check.names = FALSE)
  names(mortality)[1] <- 'Year'
  
  # Read trophic level data ----
  TLFile <- paste(InputDir,"/",InputsubDir,"/","tl_annual.csv",sep="")
  tl_annual <- read.csv(TLFile, header = T, skip = skip_this, check.names = FALSE)
  names(tl_annual)[1] <- 'Year'   # changes the name of the first column to year instead of year//group
  
  # Read chl file if present ----
  if (chl_from_file > 0) {
    ChlFile <- paste(InputDir,"/",InputsubDir,"/","chl.csv",sep="")   # Assumes a format of Year chl
    dfChl <- read.csv(ChlFile, header = T, check.names = FALSE)
    
    # make sure has columns Year and Chl
  }
  
  ## Format data
  # Convert data into a data frame and pivot tables to tidy ----
  bio <- as.data.frame(biomass) %>% 
    pivot_longer(-Year, names_to = "species_id", values_to = "biomass_tonnes")
  bio$biomass_tonnes <- area_model * bio$biomass_tonnes
  
  cat <- as.data.frame(catch) %>% 
    pivot_longer(-Year, names_to = "species_id", values_to = "catch_tonnes")
  cat$catch_tonnes <- area_model * cat$catch_tonnes
  
  cons <- as.data.frame(consumption) %>% 
    pivot_longer(-Year, names_to = "species_id", values_to = "consumption")
  
  aggland <- landings %>%
    group_by(Year, species_id) %>%
    dplyr::summarise(landings_tonnes = sum(value))
  aggland$landings_tonnes <- area_model * aggland$landings_tonnes
  aggland$species_id <- as.character(aggland$species_id)
  
  tl <- as.data.frame(tl_annual) %>% 
    pivot_longer(-Year, names_to = "species_id", values_to = "trophic_level")
  
  land <- as.data.frame(landings_fleet) %>% 
    dplyr::rename(species_id = group) %>% # rename group to species_id
    dplyr::rename(fleet_landing_tonnes = value) %>% # rename value to landing_tonnes
    mutate(species_id = as.character(species_id)) %>% 
    unique()
  land$fleet_landing_tonnes <- area_model * land$fleet_landing_tonnes
  
  removals <- as.data.frame(removals_fleet) %>% 
    dplyr::rename(species_id = group) %>% # rename group to species_id
    dplyr::rename(fleet_removals_tonnes = value) %>% # rename value to landing_tonnes
    mutate(species_id = as.character(species_id)) %>% 
    unique()
  removals$fleet_removals_tonnes <- area_model * removals$fleet_removals_tonnes
  
  # Add species names to the data frame ----
  ## Change the COL_ID to species_id so can join
  id2 <- id %>% 
    dplyr::rename(species_id = ID) %>%  # specify that the dplyr rename is the function you want & rename ID to species_id
    #    select(-COL_ID) %>% # removes the COL_ID column
    mutate(species_id = as.character(species_id)) %>% # change to character so same as other data frame
    unique() # make sure only unique values are used.
  
  # Merge the catch and biomass into one data frame ----
  df1tmp <- full_join(bio, cat) # have two columns the same so they will automatically join on these
  df1 <- full_join(df1tmp, aggland)
  df1[is.na(df1)] <- 0
  df <- full_join(df1, id2) # data frame with biomass, catch and species names
  #df <- full_join(df1tmp, id2) # data frame with biomass, catch and species names
  df <- full_join(df, cons)
  
  # New df3 which has landings and total catch per fleet in one file
  df3tmp <- full_join(land, removals)
  
  # Join landings and species names ----
  df3tmpA <- full_join(df3tmp, id2) 
  df3 <- merge(df3tmpA, df_fleet, by.x = "fleet", by.y = "ID") 
  
  # add species names to trophic level data
  tl_species <- full_join(tl, id2)
  
  # Get reference values ----
  # Depending on the value of bo_conservative set Bo
  # bo_in_same_file = 1 then take Bs from row bo_yr_row
  # bo_in_same_file = 0 and bo_conservative = 0 then take Bs from the Bo file as dicatted by bo_yr_row
  # bo_in_same_file = 0 and bo_conservative = 1 then take max(first row of biofile, row of Bo file)
  if(bo_in_same_file < 1){
    # Load Bo and Mo from files
    # bo_skip_this <- bo_header_row_ID - 1
    BoBioFile <- paste(mainDir,"/",BosubDir,"/","biomass_annual.csv",sep="")
    bo_skip_this <- skip_lines(BoBioFile)
    BoBio <- read.csv(BoBioFile, header = T, skip = bo_skip_this, check.names = FALSE)
    MoMortFile <- paste(mainDir,"/",BosubDir,"/","mortality_annual.csv",sep="")
    MoMort <- read.csv(MoMortFile, header = T, skip = bo_skip_this, check.names = FALSE)
    names(BoBio)[1] <- 'Year'
    names(MoMort)[1] <- 'Year' 
    if(bo_yr_row < 1) {
      bo_yr_row <- length(BoBio[,1])
      RefB <- BoBio[bo_yr_row,]
      RefM <- MoMort[bo_yr_row,]
      
    } else {
      RefB <- BoBio[bo_yr_row,]
      RefM <- MoMort[bo_yr_row,]
    }
    
    if(bo_conservative > 0) {
      altB1 <- as.data.frame(RefB) %>% 
        pivot_longer(-Year, names_to = "species_id", values_to = "RefB")
      altBtmp <- biomass[1,]
      altB2 <- as.data.frame(altBtmp) %>% 
        pivot_longer(-Year, names_to = "species_id", values_to = "RefB")
      altB1$Year <- RefM[1,1]
      altB2$Year <- RefM[1,1]
      dftmp <- merge(altB1, altB2, by=c("Year","species_id"))
      dftmp$RefB <- ifelse(dftmp$RefB.x > dftmp$RefB.y, dftmp$RefB.x, dftmp$RefB.y)
      dfRefB <- dftmp
      dfRefB <- subset (dfRefB, select = -c(Year, RefB.x, RefB.y))
    } else {
      dfRefB <- as.data.frame(RefB) %>% 
        pivot_longer(-Year, names_to = "species_id", values_to = "RefB")
      dfRefB <- subset (dfRefB, select = -c(Year))
    }
    dfRefM <- as.data.frame(RefM) %>% 
      pivot_longer(-Year, names_to = "species_id", values_to = "RefM")
    dfRefM <- subset (dfRefM, select = -c(Year))
    
  } else {
    # Assumes has sensible reference year defined by user
    RefB <- biomass[bo_yr_row,]
    RefM <- mortality[bo_yr_row,]
    
    dfRefB <- as.data.frame(RefB) %>% 
      pivot_longer(-Year, names_to = "species_id", values_to = "RefB")
    dfRefB <- subset (dfRefB, select = -c(Year))
    dfRefM <- as.data.frame(RefM) %>% 
      pivot_longer(-Year, names_to = "species_id", values_to = "RefM")
    dfRefM <- subset (dfRefM, select = -c(Year))
    
  }
  dfRef <- merge(dfRefM, dfRefB, by=c("species_id")) # Intentional reuse
  dfRef$RefBo <- area_model * dfRef$RefB
  
  # Filter out burn-in period ----
  aggland <- dplyr::filter(aggland, Year > yr_ignore)
  bio <- dplyr::filter(bio, Year > yr_ignore)
  biomass <- dplyr::filter(biomass, Year > yr_ignore)
  cat <- dplyr::filter(cat, Year > yr_ignore)
  catch <- dplyr::filter(catch, Year > yr_ignore)
  cons <- dplyr::filter(cons, Year > yr_ignore)
  consumption <- dplyr::filter(consumption, Year > yr_ignore)
  df <- dplyr::filter(df, Year > yr_ignore)
  df1 <- dplyr::filter(df1, Year > yr_ignore)
  df3 <- dplyr::filter(df3, Year > yr_ignore)
  land <- dplyr::filter(land, Year > yr_ignore)
  landings <- dplyr::filter(landings, Year > yr_ignore)
  landings_fleet <- dplyr::filter(landings_fleet, Year > yr_ignore)
  mortality <- dplyr::filter(mortality, Year > yr_ignore)
  mortality_fleet <- dplyr::filter(mortality_fleet, Year > yr_ignore)
  removals <- dplyr::filter(removals, Year > yr_ignore)
  removals_fleet <- dplyr::filter(removals_fleet, Year > yr_ignore)
  tl <- dplyr::filter(tl, Year > yr_ignore)
  tl_annual <- dplyr::filter(tl_annual, Year > yr_ignore)
  tl_species <- dplyr::filter(tl_species, Year > yr_ignore)
  
  return(list(aggland = aggland, bio = bio, biomass = biomass, cat = cat, 
              catch = catch, cons = cons, consumption = consumption, id2 = id2,
              df = df, df1 = df1, df3 = df3, dfRef = dfRef, df_fleet = df_fleet,
              land = land, landings = landings, landings_fleet = landings_fleet, 
              mortality = mortality, mortality_fleet = mortality_fleet, 
              removals = removals, removals_fleet = removals_fleet, 
              tl = tl, tl_annual = tl_annual, tl_species = tl_species))
}
  
########################### DIAGNOSTIC PLOT FUNCTIONS ###########################

# Function to call all the time series plotting functions
##' @title Create time series plots
##' @param file dataframes loaded from EwE output files
##' @return Updated df3 dataframe
##' @author Beth Fulton
print_diagnostic_timeseries <- function(df, df3, dfRef, df_fleet, OutDir){
  
  # Biomass time series
  create_bio_ts_plots(df, dfRef, OutDir)
  
  # Catch time series
  create_catch_ts_plots(df3, df_fleet, OutDir)
  
  # Discard time series
  df3 <- create_discard_ts_plots(df3, df_fleet, OutDir)
  
  return(df3)
  
}

# Function to plot biomass time series - so you can check the time series 
# (and ths the rest of the indicators) makes sense
##' @title Create biomass time series plots per model group. This automatically determines the number of rows based on the ewe version
##' @param file dataframe df - contains the biological time series loaded from EwE output files
##' @author Beth Fulton
create_bio_ts_plots <- function(df, dfRef, OutDir) {
  this_nSPName <-length(unique(df$Group))
  this_SpeciesNames <- unique(df$Group)
  dftmp <- merge(df, dfRef, by = c("species_id"))
  
  # plot out catch for each species for each fleet 
  nLoop <- ceiling (this_nSPName / nPanel)
  counter1 <- 1
  counter2 <- nPanel
  for (i in 1:nLoop) {
    check_name <- this_SpeciesNames[counter1:counter2]  
    extractdf <- dplyr::filter(dftmp, Group %in% check_name)
    outfilepng <- paste("biomass",i,"-page",i,".png",sep="")
    
    p <- ggplot(extractdf, aes(x = Year, y = biomass_tonnes)) +
      geom_line(color = tscolor, linewidth = 1.2) +
      geom_line(aes(x = Year, y = RefBo), 
                color = "black", linewidth = 0.8, linetype="twodash") +
      theme(legend.position = "none") +
      expand_limits(y=0) +
      labs(x = "Year",
           y = "Biomass (tonnes)") +
      facet_wrap(~ Group, scales = 'free') +
      theme(axis.text.y=element_text(size=8),
            axis.text.x=element_text(size=8),
            axis.title=element_text(size=12,face="bold"), 
            strip.text = element_text(face="bold", size=10))
    
    print(p)
    
    ggsave(filename = outfilepng, 
           plot = p, 
           path = OutDir)
    counter1 <- counter1 + nPanel
    counter2 <- counter2 + nPanel
  }
}

# Function to plot catch time series per fleet - so you can check the time series 
# (and ths the rest of the indicators) makes sense
##' @title Create catch time series plots per model group per model fleet. 
##' @param file dataframe df - contains the biological time series loaded from EwE output files
##' @author Beth Fulton
create_catch_ts_plots <- function(df3, df_fleet, OutDir) {

  # subset out data by fleet
  # to see how many fleets
  df4 <- na.omit(df3) 
  #df4 <- df3  # Beth feels the NA omission is important, others disagree
  tmp <- unique(df4[c("fleet")])
  nFleet <- max(tmp[,1])
  
  for (i in 1:nFleet) {
    txttitle <- paste0("Catch time series by species by fleet ", i)
    catch_fi <- df3 %>% dplyr::filter(fleet == i)
    this_fleet <- df_fleet[i,2]
    this_nSPName <- length(unique(catch_fi$Group))
    this_SpeciesNames <- unique(catch_fi$Group)
    
    # plot out catch for each species for each fleet 
    nLoop <- ceiling (this_nSPName / nPanel)
    counter1 <- 1
    counter2 <- nPanel
    counter3 <- 1
    for (iLoop in 1:nLoop) {
      
      check_name <- this_SpeciesNames[counter1:counter2]  
      extractC <- dplyr::filter(catch_fi, Group %in% check_name)
      lenC <- dim(extractC)
      l1 <- lenC[1]
      if(l1 > 0) {
        outfilepng <- paste0("catch-fleet",i,"-page",counter3,".png")
        catch_pi <- ggplot(extractC, 
                           aes(x = Year, y = fleet_landing_tonnes)) +
          geom_line(color = tscolor, linewidth = 1.2) +
          theme(legend.position = "none") +
          expand_limits(y=0) +
          labs(
            title = txttitle,
            subtitle = this_fleet,
            x = "Year",
            y = "Amount of catch (tonnes)") +
          facet_wrap(~ Group, scales = 'free') +
          theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=8), axis.title=element_text(size=12,face="bold"), strip.text = element_text(face="bold", size=10))
        
        print(catch_pi)
        
        # save out plots
        ggsave(filename = outfilepng, plot = catch_pi, path = OutDir)
        counter3 <- counter3 + 1
      }
      
      counter1 <- counter1 + nPanel
      counter2 <- counter2 + nPanel
    }
  }
}

# Function to plot discard time series per fleet - so you can check the time series 
# (and ths the rest of the indicators) makes sense
##' @title Create discard time series plots per model group per model fleet. Calculated as discard = catch-landings
##' @param file dataframe df - contains the biological time series loaded from EwE output files
##' @return Updated df3 dataframe
##' @author Beth Fulton
create_discard_ts_plots <- function(df3, df_fleet, OutDir) {
  # Number of groups to print
  this_nSPName <-length(unique(df3$Group))
  this_SpeciesNames <- unique(df3$Group)
  
  this_nFName <- length(unique(df_fleet$FleetName))
  this_FleetNames <- df_fleet$FleetName
  
  ## First discard time series as output by Ecosim
  # check to see if there are any columns in df3 that have NA values
  list_na <- colnames(df3)[apply(df3, 2,anyNA)]
  
  # calculate discards - put this into df3 dataframe not new discards dataframe - adjusted script below to match
  df3$discard_tonnes <- df3$fleet_removals_tonnes - df3$fleet_landing_tonnes
  
  # Create new dataframes that is total discards across all fleets (per species and year) and across all species (per fleet)
  discardsBySpecies <- df3 %>%
    group_by(Year, species_id, Group) %>%
    dplyr::summarise(TotalDiscards = sum(discard_tonnes))
  
  discardsByFleet <- df3 %>%
    group_by(Year, FleetName) %>%
    dplyr::summarise(TotalDiscards = sum(discard_tonnes))
  
  # Convert any NA with zeros
  discardsBySpecies[is.na(discardsBySpecies)] <- 0
  discardsByFleet[is.na(discardsByFleet)] <- 0
  
  # Filter out any years that don't make sense
  discardsBySpecies <- discardsBySpecies %>% dplyr::filter(Year > 0)
  discardsByFleet <- discardsByFleet %>% dplyr::filter(Year > 0)
  
  # check the maximum discard  
  ymax <- 1.05 * max(discardsBySpecies$TotalDiscards, na.rm = TRUE)
  ymin <- 0.95 * min(discardsBySpecies$TotalDiscards, na.rm = TRUE)
  if (ymin < 0) { ymin <- 0 }
  
  # plot the discards for each species
  nLoop <- ceiling (this_nSPName / nPanel)
  counter1 <- 1
  counter2 <- nPanel
  for (iLoop in 1:nLoop) {
    outPlotName <- paste(OutDir,"/Total_Discards_per_species-page",iLoop,".png",sep="")
    check_name <- this_SpeciesNames[counter1:counter2]  
    extractC <- dplyr::filter(discardsBySpecies, Group %in% check_name)
    check_dim <- dim(extractC)
    check_this <- check_dim[[1]]
    if (check_this > 0) {
      p_dis <- ggplot(extractC, aes(x = Year,  y = TotalDiscards))+
        geom_line(color = tscolor, linewidth = 1.2) +
        #geom_point() +
        #ylim(ymin, ymax)+ 
        theme(legend.position = "none")+
        labs(
          title = "Discard tonnes of species by species",
          x = "Year",
          y = "Amount of catch discarded (tonnes)"
        ) +
        facet_wrap(~ Group, scales = 'free_y', ncol = 6) +
        theme(axis.text.y = element_text(size=8),
              axis.text.x = element_text(size=8),
              axis.title = element_text(size=12,face="bold"), 
              strip.text = element_text(face="bold", size=10))
      
      print(p_dis)
      
      ggsave(p_dis, file=outPlotName)
    }
    
    counter1 <- counter1 + nPanel
    counter2 <- counter2 + nPanel
  }
  
  ## Second maximum discards and total discards
  # check the maximum discard  
  ymax <- 1.05 * max(discardsByFleet$TotalDiscards, na.rm = TRUE)
  ymin <- 0.95 * min(discardsByFleet$TotalDiscards, na.rm = TRUE)
  if (ymin < 0) { ymin <- 0 }
  
  # have a look at each species by fleet
  nLoop <- ceiling (this_nFName / nPanel)
  counter1 <- 1
  counter2 <- nPanel
  for (iLoop in 1:nLoop) {
    outPlotName <- paste(OutDir,"/Total_Discards_per_fleet-page",iLoop,".png",sep="")
    check_name <- this_FleetNames[counter1:counter2]  
    extractC <- dplyr::filter(discardsByFleet, FleetName %in% check_name)
    p_dis2 <- ggplot(extractC, 
               aes(x = Year,  y = TotalDiscards))+
      geom_line(color = tscolor, linewidth = 1.2) +
      # geom_point() +
      # ylim(ymin, ymax)+ 
      theme(legend.position = "none") +
      labs(
        title = "Discarded tonnes of catch by fleet",
        x = "Year",
        y = "Amount of catch discarded (tonnes)"
      ) +
      facet_wrap(~ FleetName, scales = 'free_y', ncol = 6) +
      theme(axis.text.y = element_text(size=8),
            axis.text.x = element_text(size=8),
            axis.title = element_text(size=12,face="bold"), 
            strip.text = element_text(face="bold", size=10))
    
    print(p_dis2)
    
    ggsave(p_dis2, file=outPlotName)
    
    counter1 <- counter1 + nPanel
    counter2 <- counter2 + nPanel
  }
  
  ## Finally discard rates
  # calculate discard rate
  df3$discard_rate <- df3$discard_tonnes / df3$fleet_removals_tonnes
  
  discard_stats <- df3 %>% 
    dplyr::group_by(FleetName, Year) %>% 
    dplyr::summarise(discard_mean = mean(discard_rate, na.rm = TRUE), 
                     discard_mean_sd = sd(discard_rate, na.rm = TRUE),
                     discard_n = n()) %>% 
    dplyr::filter(!is.na(FleetName))
  
  # Reset any NAs (std deviations when n is 1)
  discard_stats[is.na(discard_stats)] <- 0
  discard_stats$max_CI <- discard_stats$discard_mean + (1.96 * discard_stats$discard_mean_sd / sqrt(discard_stats$discard_n))
  discard_stats$min_CI <- discard_stats$discard_mean - (1.96 * discard_stats$discard_mean_sd / sqrt(discard_stats$discard_n))
  
  # maximum and minimun discard rates
  ymax <- 1.05 * max(discard_stats$max_CI, na.rm = TRUE)
  ymin <- 0.95 * min(discard_stats$min_CI, na.rm = TRUE)
  if (ymin < 0) { ymin <- 0 }
  
  # plot the discard rate per fleeet
  nLoop <- ceiling (this_nFName / nPanel)
  counter1 <- 1
  counter2 <- nPanel
  for (iLoop in 1:nLoop) {
    outPlotName <- paste(OutDir,"/Discard_rates_per_fleet-page",iLoop,".png",sep="")
    check_name <- this_FleetNames[counter1:counter2]  
    extractC <- dplyr::filter(discard_stats, FleetName %in% check_name)
    p_dr <- ggplot(extractC, aes(x = Year,  y = discard_mean))+
      geom_line(color = tscolor, linewidth = 1.2) +
      #geom_point() +
      #ylim(ymin, ymax)+ 
      geom_line(aes(y = max_CI, color=tscolor), linetype="twodash") +
      geom_line(aes(y = min_CI, color=tscolor), linetype="twodash") +
      theme(legend.position = "none")+
      labs(
        title = "Discard rates by fleet",
        x = "Year",
        y = "Discard rates (proportion)"
      ) +
      facet_wrap(~ FleetName, scales = 'free_y', ncol = 6) +
      theme(axis.text.y = element_text(size=8),
            axis.text.x = element_text(size=8),
            axis.title = element_text(size=12,face="bold"), 
            strip.text = element_text(face="bold", size=10))
    
    print(p_dr)
    
    ggsave(p_dr, file=outPlotName)
    
    counter1 <- counter1 + nPanel
    counter2 <- counter2 + nPanel
  }
  
  return(df3)
  
}

# Function to calculate PCA of catch composition - using spectral decomposition 
# approach via the princomp approach
##' @title Principle component analysis
##' @return List of PCA results
##' @author Beth Fulton
PCA_through_time <- function(df, df3) {
  
  # Recast so the species are the columns - starting by getting catch per year from old df3
  # CdataS <- df3 %>%
  #  group_by(Year, Group) %>%
  #  dplyr::summarise(TotalYield = mean(catch_tonnes))
  # Now just take it straight from df
  CdataS <- df
  
  # Replace spaces in names with "_"
  CdataS$Group <- gsub("\\(", "", CdataS$Group)
  CdataS$Group <- gsub("\\/", "", CdataS$Group)
  CdataS$Group <- gsub(")", "", CdataS$Group)
  CdataS$Group <- gsub("&", "", CdataS$Group)
  CdataS$Group <- gsub("-", "", CdataS$Group)
  CdataS$Group <- gsub(" ", "", CdataS$Group)
  # CdataS$Group <- gsub("  ", " ", CdataS$Group)
  # CdataS$Group <- gsub(" ", "_", CdataS$Group)
  # SP_as_colA <- reshape2::dcast(CdataS, Year ~ Group, value.var = "TotalYield")
  
  SP_as_colA <- reshape2::dcast(CdataS, Year ~ Group, 
                                value.var = "catch_tonnes")
  
  # Replace NA with zeros
  SP_as_colA[is.na(SP_as_colA)] <- 0
  
  # Strip out columns of all zeros
  SP_as_col <- SP_as_colA[, colSums(SP_as_colA != 0) > 0]
  
  # Check dimensions
  dimC <- dim(SP_as_col)
  this_nSP <- dimC[[2]]
  this_nYr <- dimC[[1]]
  
  can_plot <- 1
  if(this_nSP > this_nYr) {
    can_plot <- 0
  }
  
  if (can_plot > 0) { # Can do a spectral decomposition as more years than groups
    # Get the column headers
    pc.f <- formula(paste("~", paste(names(SP_as_col)[2:dimC[2]], collapse = "+")))
    
    # PCA calculations - using spectral decomposition approach via the princomp approach
    pl.pca <- princomp(pc.f, cor=TRUE, data=SP_as_col)
    
    # Put on row named
    row.names(pl.pca$scores) <- SP_as_col$Year
  } else {
    pl.pca <- 0
  }
  
  return(list(pca_data = SP_as_col, pca_results = pl.pca, can_plot = can_plot))
}


# Function to calculate PCA of catch composition
##' @title Principle component analysis
##' @return PCA biplot
##' @author Beth Fulton
plot_biplot <- function(pl.pca, OutDir) {
  outBiplot <- paste(OutDir,"/PCA_Biplot_Thru_Time.png",sep="")   
  dfPCA <- data.frame(comp1=pl.pca$scores[,1],
                      comp2=pl.pca$scores[,2])
  
  p_bi <- ggplot(data = dfPCA, aes(x=comp1, y=comp2, group=1)) +
    geom_point(size=5, aes(colour=rownames(pl.pca$scores))) +
    geom_path(size = 0.2) +
    geom_text(label=rownames(pl.pca$scores)) + 
    theme(legend.position="none")
  
  ggsave(p_bi, file=outBiplot)
  
  pca_p <- ggplot(data = dfPCA, aes(x=comp1, y=comp2, group=1)) +
    geom_point(size=5, aes(colour=rownames(pl.pca$scores))) +
    geom_path(size = 0.2) +
    geom_text(label=rownames(pl.pca$scores)) + 
    theme(legend.position="none")
  
  return(pca_p)
  
}

# Function to calculate PCA of catch composition - using spectral decomposition 
# approach via the princomp approach. This includes plotting new PCA results
##' @title Principle component analysis - singular value decomposition
##' @return List of plots of the PCA results
##' @author Beth Fulton

##' @title Principle component analysis - spectral decomposition 
##' @return List of PCA results
##' @author Beth Fulton
PCA_of_composition <- function(SP_as_col) {
  
  ## Use prcomp() - this uses singular value decomposition 
  # This can handle cases where there is more variables than observations
  res.pca <- prcomp(SP_as_col, scale = TRUE)
  
  # Scree plot
  res_plot1 <- fviz_eig(res.pca)
  
  # Biplot
  res_plot2 <- fviz_pca_ind(res.pca,
               col.ind = "contrib", # Color by congtribution
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,     # Avoid text overlapping
               #label=SP_as_col$Year
  ) +
    labs(title ="PCA", x = "PC1", y = "PC2")
  
  # Weighted variables on biplots
  res_plot3 <- fviz_pca_var(res.pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = PCA_color_gradient,
               repel = TRUE     # Avoid text overlapping
  )
  
  # Compute hierarchical clustering on principal components
  res.pca4 <- PCA(SP_as_col, ncp = 3, graph = FALSE)
  res.hcpc <- HCPC(res.pca4, graph = FALSE)
  
  # Dendrogram of PCA scores
  res_plot4 <- fviz_dend(res.hcpc, 
            cex = 0.7,                     # Label size
            palette = "jco",               # Color palette see ?ggpubr::ggpar
            rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
            rect_border = "jco",           # Rectangle color
            labels_track_height = 0.8      # Augment the room for labels
  )
  
  # Biplot clusters of the PCA factors
  res_plot5 <- fviz_cluster(res.hcpc,
               repel = TRUE,            # Avoid label overlapping
               show.clust.cent = TRUE, # Show cluster centers
               palette = "jco",         # Color palette see ?ggpubr::ggpar
               ggtheme = theme_bw(),
               main = "Factor map"
  )
  
  # Save all plots
  outScree <- paste(OutDir,"/PCA_SVD_Scree.png",sep="")   
  ggsave(plot = res_plot1, file=outScree)
  outBiPlot <- paste(OutDir,"/PCA_SVD_bioplot.png",sep="")   
  ggsave(plot = res_plot2, file=outBiPlot)
  outVar <- paste(OutDir,"/PCA_SVD_weighted_var.png",sep="")  
  ggsave(plot = res_plot3, file=outVar)
  outDendro <- paste(OutDir,"/PCA_SVD_dendrogram.png",sep="")  
  ggsave(plot = res_plot4, file=outDendro)
  outClust <- paste(OutDir,"/PCA_SVD_cluster.png",sep="")  
  ggsave(plot = res_plot5, file=outClust)
  
  pca_res <- list(prcomp_pca = res.pca, prcomp_hcpc = res.hcpc, 
       res_plot1 = res_plot1, res_plot2 = res_plot2, res_plot3 = res_plot3,
       res_plot4 = res_plot4, res_plot5 = res_plot5)
  
  
  return(pca_res)
  
}

# Plot a heatmap of the contribution of each group to the catch through time. Ranks contribution through time
# of species independently and then clusters the ranking to plot species with similar patterns adjactent.
# Annotate this later to mark points where technology, fishing bhaviour, management requirements or poilicy changed.
# Also make on significant environmental events. Helps with interpretation
##' @title Calculate heatmap.
##' @return List of plots of the indicators
##' @author Beth Fulton
do_heatmap <- function(df, OutDir){
  CdataS <- df
  
  # Replace spaces in names with "_"
  CdataS$Group <- gsub("\\(", "", CdataS$Group)
  CdataS$Group <- gsub("\\/", "", CdataS$Group)
  CdataS$Group <- gsub(")", "", CdataS$Group)
  CdataS$Group <- gsub("&", "", CdataS$Group)
  CdataS$Group <- gsub("-", "", CdataS$Group)
  CdataS$Group <- gsub(" ", "", CdataS$Group)
  SP_as_colA <- reshape2::dcast(CdataS, Year ~ Group, value.var = "landings_tonnes")
  
  # Replace NA with zeros
  SP_as_colA[is.na(SP_as_colA)] <- 0
  
  # Strip out columns of all zeros
  SP_as_col <- SP_as_colA[, colSums(SP_as_colA != 0) > 0]
  
  # Create row names and trim off first column
  rownames(SP_as_col) <- SP_as_col[,1]
  SP_as_col <- SP_as_col %>% dplyr::select(-Year)
  
  # Transpose SP_as_col so years are column headers
  # Here the va;ues are proportion of the catch but could be absolute catch
  Yr_as_col <- t(SP_as_col)
  
  #Replace NA with 0
  Yr_as_col[is.na(Yr_as_col)] <- 0
  Prop_Catch <- apply(Yr_as_col, 2, function(i) i/sum(i))
  
  nSP <- length(Yr_as_col[,1])
  nYr <- length(Yr_as_col[1,])
  
  # Create new df to populate
  drRanked <- data.frame(matrix(ncol = nYr, nrow = nSP))
  colnames(drRanked)[1:nYr] = as.character(colnames(Yr_as_col)[1:nYr])
  rownames(drRanked)[1:nSP] = as.character(rownames(Yr_as_col)[1:nSP])
  
  # o ranking (using order() routine as want largest value to be rated 1)
  for (i in 1:nSP) {
    drRanked[i,] <- frank(Yr_as_col[i,], ties.method ="dense")
  }
  
  # Create the heatmap
  HeatmapFilename <- paste0(OutDir, "/Landings heatmap.pdf") 
  HeatmapTitle <- paste0("Catch composition heatmap - catch data") 
  p_heat <-heatmaply(drRanked,
                xlab = "Year",
                ylab = "Species",
                main = HeatmapTitle,
                dendrogram = "row",
                fontsize_col = 6,
                fontsize_row = 8
                #file = HeatmapFilename
  )
  
  # ggsave(p, file = HeatmapFilename)
  
  htmlwidgets::saveWidget(p_heat, "temp.html")
  webshot::webshot("temp.html", paste0(OutDir, "Landings heatmap.png"), 
                   delay = 3, vwidth = 1200, vheight = 800)
  #file.remove("temp.html")
  
  # No longer putting out proportional heatmaps as not helpful
  #HeatmapFilename <- paste("Proportion of landings heatmap.pdf",sep="") 
  #HeatmapTitle <- paste("Catch composition heatmap - proportion of catchh",sep="") 
  #p_propheat <- heatmaply(Prop_Catch,
  #               xlab = "Year",
  #               ylab = "Species",
  #               main = HeatmapTitle,
  #               dendrogram = "row",
  #               fontsize_col = 6,
  #               fontsize_row = 8#,
  #                #file = HeatmapFilename
  #)
  ## ggsave(p_propheat, file = HeatmapFilename)
  
  #htmlwidgets::saveWidget(p_propheat, "temp.html")
  #webshot::webshot("temp.html", 
  #                 paste0(OutDir, "Proportion of landings heatmap.png"), 
  #                 delay = 3, vwidth = 1200, vheight = 800)
  #file.remove("temp.html")
  
  
}

########################### OTHER INDICATOR FUNCTIONS ###########################


# Function to calculate the indicators recommended by INDISEAS and subsequent work by Jason Link
##' @title Calculate EBFM indicators.
##' @return List of plots of the indicators
##' @author Beth Fulton
calculate_EBFM_indices <-function(df, df3, df_fleet, mortality, biomass, id2, OutDir, effort_filename, chl_from_file, yr_ignore) {
  
  # Explotation rate
  plot_exploitation_rate(df, OutDir)
  
 # plot 1/ (biomass/landings)
  plot_biomass_landings(df, OutDir)
  
  # Plot biomass ratios
  plot_biomass_ratios(df, OutDir)
  
  # Plot of CV ratios
  plot_CV(df, OutDir)
  
  # Plot INDISEAS indicators
  plot_INDISEAS_indices(df, OutDir)
  
  # Plot system level productivity and exploitation indices recommended by Link and Watson 2019
  plot_Link_Watson(df, mortality, biomass, id2, chl_from_file, OutDir)
  
  # Cumulative biomass-TL plot
  cumB_TL_plot <- plot_CumB_TL(df, tl)
  
  return(cumB_TL_plot)
  
}

# Function to calculate the indicators recommended by INDISEAS and subsequent work by Jason Link
##' @title Calculate and plot exploitation rate
##' @author Beth Fulton
plot_exploitation_rate <- function(df, OutDir) {
  # Local number of groups to plot
  this_nSPName <- length(unique(df$Group))
  this_SpeciesNames <- unique(df$Group)
  
  # calculate exploitation rate
  df$exploit_rate <- 100 * (df$catch_tonnes / df$biomass_tonnes)
  
  # check for N values
  list_na <- colnames(df)[apply(df, 2,anyNA)]
  
  # check the maximum exploitation rate to make sure not over 100%
  ymax <- 1.05 * max(df$exploit_rate)
  subset(df, exploit_rate == max(exploit_rate)) # gives the row/s with the max exploitation rate - might need to check original data!
  
  # plot trends in exploitation rate
  nLoop <- ceiling (this_nSPName / nPanel)
  counter1 <- 1
  counter2 <- nPanel
  for (iLoop in 1:nLoop) {
    outPlotName <- paste(OutDir,"/Exploitation_rates_per_species-page",iLoop,".png",sep="")
    check_name <- this_SpeciesNames[counter1:counter2]  
    extractC <- dplyr::filter(df, Group %in% check_name)
    p_ex <- ggplot(data = extractC, mapping = aes(x = Year, y = exploit_rate)) +
      geom_line(color = tscolor, linewidth = 1.2) +
      expand_limits(y=0) +
      theme(legend.position = "none")+
      labs(
        title = "Expolitation rate by species",
        x = "Year",
        y = "Percentage of stock exploited") +
      facet_wrap(~ Group, scales = 'free') +
      theme(axis.text.y = element_text(size=8),
            axis.text.x = element_text(size=8),
            axis.title = element_text(size=12,face="bold"), 
            strip.text = element_text(face="bold", size=10))
    
    print(p_ex) 
    
    ggsave(p_ex, file=outPlotName)
    
    counter1 <- counter1 + nPanel
    counter2 <- counter2 + nPanel
  }
  
}


# Function to plot 1 / (biomass / landings) at the system level
##' @title Calculate and plot 1 / (biomass / landings) at the system level
##' @author Beth Fulton
plot_biomass_landings <- function(df, OutDir) {
  
  # biomass time series using 1/(biomass/landings)
  # Sum biomass of those species with landings and the calculate the index
  biomass_calc_tmp <-df %>% dplyr::filter(landings_tonnes > 0)
  biomass_calc <- biomass_calc_tmp %>%
    group_by(Year) %>%
    dplyr::summarise(TotAvailBio = sum(biomass_tonnes), TotLandings = sum(landings_tonnes))
  biomass_calc$bioindx  <- 1.0 / (biomass_calc$TotAvailBio / biomass_calc$TotLandings) 
  
  # plot the index
  plot_biomass_indx <- paste(OutDir,"/one_on_bio_on_landings_indx.png",sep="")
  
  p_bl <- ggplot(biomass_calc, aes (x = Year, y = bioindx)) +
    geom_line(color = aggtscolor, linewidth = 1.5) +
    labs(
      title = "1/(Biomass/Landings) time series",
      x = "Year",
      y = "1/(Biomass/Landings)"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(p_bl)
  
  ggsave(p_bl, file=plot_biomass_indx)
  
}


# Function to plot biomass ratios Biomass ratios piscivorous:zooplanktivorous fish
# and pelagic:demersal)
# Reads in classification from species_ID.csv file 
# Then calculate the ratios
##' @title Calculate and plot biomass ratios
##' @author Beth Fulton
plot_biomass_ratios <- function(df, OutDir) {
  
  # get aggregate values - Pelagic:Demersal
  PDtmp <- df %>%
    group_by(Year, PelDemID) %>%
    dplyr::summarise(BioPerClass = sum(biomass_tonnes), CatPerClass = sum(biomass_tonnes))
  PDindB <- reshape2::dcast(PDtmp, Year ~ PelDemID, value.var = "BioPerClass")
  PDindB[is.na(PDindB)] <- 0 # Replace NA with zeros
  if (length(PDindB$Demersal) > 0) {
    PDindB$PDb <- PDindB$Pelagic / PDindB$Demersal
  }else {
    PDindB$PDb <- 0  # Not valid
  }
  
  PDindC <- reshape2::dcast(PDtmp, Year ~ PelDemID, value.var = "CatPerClass")
  PDindC[is.na(PDindC)] <- 0 # Replace NA with zeros
  if (length(PDindC$Demersal) > 0) {
    PDindC$PDc <- PDindC$Pelagic / PDindC$Demersal
  } else {
    PDindC$PDc <- 0 # Not valid
  }
  
  PDind <- merge(PDindB, PDindC, by = c("Year"))
  
  # plot the index
  plot_pd_indx <- paste(OutDir,"/pelagic_demersal_biomass_ratio_indx.png",sep="")
  p_bd <- ggplot(PDind, aes (x = Year, y = PDb), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    geom_line(aes (x = Year, y = PDc), color = aggtscolor2, linewidth = 1.2) +
    theme(legend.position = "right")+
    expand_limits(y=0.5) +
    labs(
      title = "Pelagic:Demersal time series",
      subtitle = "Dark line is biomass, light line is catch",
      x = "Year",
      y = "Pelagic:Demersal"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(p_bd)
  
  ggsave(p_bd, file=plot_pd_indx)
  
  
  # get aggregate values - Piscivorous:Zooplanktivorous
  PZtmp <- df %>%
    group_by(Year, PISCZOOPL) %>%
    dplyr::summarise(BioPerClass = sum(biomass_tonnes), CatPerClass = sum(biomass_tonnes))
  PZindB <- reshape2::dcast(PZtmp, Year ~ PISCZOOPL, value.var = "BioPerClass")
  PZindB[is.na(PZindB)] <- 0 # Replace NA with zeros
  if(length(PZindB$Planktivorous) > 0) {
    PZindB$PZb <- PZindB$Piscivorous / PZindB$Planktivorous
  } else {
    PZindB$PZb <- 0 # Not valid
  }
  
  PZindC <- reshape2::dcast(PZtmp, Year ~ PISCZOOPL, value.var = "CatPerClass")
  PZindC[is.na(PZindC)] <- 0 # Replace NA with zeros
  if(length(PZindC$Planktivorous) > 0){
    PZindC$PZc <- PZindC$Piscivorous / PZindC$Planktivorous
  } else {
    PZindC$PZc <- 0
  }
  
  PZind <- merge(PZindB, PZindC, by = c("Year"))
  
  # plot the index
  plot_pd_indx <- paste(OutDir,"/pisciv_zoopl_biomass_ratio_indx.png",sep="")
  p_PZ <- ggplot(PZind, aes (x = Year, y = PZb), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    geom_line(aes (x = Year, y = PZc), color = aggtscolor2, linewidth = 1.2) +
    theme(legend.position = "right")+
    expand_limits(y=0.5) +
    labs(
      title = "Piscivore:Zooplanktivore time series",
      subtitle = "Dark line is biomass, light line is catch",
      x = "Year",
      y = "Piscivore:Zooplanktivore"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(p_PZ)
  
  ggsave(p_PZ, file=plot_pd_indx)
  
}


# Function to plot  CVi vs CVt -coefficient of variation of ecosystem (t) versus individual (i)
##' @title Calculate and plot CV ratio
##' @author Beth Fulton
plot_CV <- function(df, OutDir) {
 
  ## Make sure have SP_as_col
  CdataS <- df
  
  # Replace spaces in names with "_"
  CdataS$Group <- gsub("\\(", "", CdataS$Group)
  CdataS$Group <- gsub("\\/", "", CdataS$Group)
  CdataS$Group <- gsub(")", "", CdataS$Group)
  CdataS$Group <- gsub("&", "", CdataS$Group)
  CdataS$Group <- gsub("-", "", CdataS$Group)
  CdataS$Group <- gsub(" ", "", CdataS$Group)
  #CdataS$Group <- gsub("  ", " ", CdataS$Group)
  #CdataS$Group <- gsub(" ", "_", CdataS$Group)
  #SP_as_colA <- reshape2::dcast(CdataS, Year ~ Group, value.var = "TotalYield")
  
  SP_as_colA <- reshape2::dcast(CdataS, Year ~ Group, value.var = "landings_tonnes")
  
  # Replace NA with zeros
  SP_as_colA[is.na(SP_as_colA)] <- 0
  
  # Strip out columns of all zeros
  SP_as_col <- SP_as_colA[, colSums(SP_as_colA != 0) > 0]
  
  ## Melted version that have to cast
  outNameA <- paste(OutDir,"/CV_3year_period.csv",sep="")
  outNameB <- paste(OutDir,"/CV_5year_period.csv",sep="")
  plot_indx <- paste(OutDir,"/CV_indx.png",sep="")
  
  nCol <- dim(SP_as_col)[2]
  NameEx <- colnames(SP_as_col)[2:nCol]
  nSPNameEx <- nCol - 1
  #NameEx[nSPNameEx]
  NameExT <- append(NameEx,"Total")
  SP_as_col$Total <- rowSums(SP_as_col[ , c(2:nCol)], na.rm=TRUE)
  
  windowA <- 3
  windowB <- 5
  
  funMeanA <- function(x) zoo::rollmean(x, windowA, na.pad=TRUE, align="right")
  funMeanB <- function(x) zoo::rollmean(x, windowB, na.pad=TRUE, align="right")
  funSDA <- function(x) zoo::rollapply(x, width = windowA, FUN = sd, fill=NA, align='right') 
  funSDB <- function(x) zoo::rollapply(x, width = windowB, FUN = sd, fill=NA, align='right') 
  
  dfoutSDA <- SP_as_col %>% 
    mutate_at(NameExT, funSDA)
  dfoutSDB <- SP_as_col %>% 
    mutate_at(NameExT, funSDB)
  dfoutMA <- SP_as_col %>% 
    mutate_at(NameExT, funMeanA)
  dfoutMB <- SP_as_col %>% 
    mutate_at(NameExT, funMeanB)
  
  NLose <- windowA - 1
  dfCVA <- 100.0 * dfoutSDA / dfoutMA
  dfCVA[is.na(dfCVA)] <- 0 # Replace NA with zeros
  dfCVtest <- subset (dfCVA, select = -c(Year, Total))
  dfCVA$SppAboveCVt <- rowSums(dfCVtest > dfCVA$Total)
  dfCVA$PropAboveCVt <- dfCVA$SppAboveCVt / nSPNameEx
  dfCVA$Year <- dfoutSDA$Year
  dfCVA <- dfCVA[-c(1:NLose),] # Drop rows of NAs from the start
  
  NLose <- windowB - 1
  dfCVB <- 100.0 * dfoutSDB / dfoutMB
  dfCVB[is.na(dfCVB)] <- 0 # Replace NA with zeros
  dfCVtest <- subset (dfCVB, select = -c(Year, Total))
  dfCVB$SppAboveCVt <- rowSums(dfCVtest > dfCVB$Total)
  dfCVB$PropAboveCVt <- dfCVB$SppAboveCVt / nSPNameEx
  dfCVB$Year <- dfoutSDB$Year
  dfCVB <- dfCVB[-c(1:NLose),] # Drop rows of NAs from the start
  
  # save results file
  write.csv(dfCVA, file = outNameA, row.names = FALSE)
  write.csv(dfCVB, file = outNameB, row.names = FALSE)
  
  # Plot results
  p_CV <- ggplot(data = NULL) +
    geom_line(data = dfCVA, aes (x = Year, y = PropAboveCVt), color = aggtscolor, linewidth = 1.5) +
    geom_line(data = dfCVB, aes (x = Year, y = PropAboveCVt), color = aggtscolor2, linewidth = 1.2, linetype="twodash") +
    theme(legend.position = "right")+
    expand_limits(y=0) +
    labs(
      title = "Proportion of species with CV > CVtotal",
      x = "Year",
      y = "Proportion"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(p_CV)
  
  ggsave(p_CV, file=plot_indx)
  
}

# Function to plot a number of system level indices recommended by INDISEAS
# Mean length in a community or catch - Assumes Linf is in Species_ID.csv file
# Mean trophic level of community, catch - Trophic level calculated by EwE
# model based on diet (not biomass).
# Mean longevity - Assumes MaxAge is in Species_ID.csv file
# Proporation of large fish (LFI) - based on Linf in Species_ID.csv file
##' @title Calculate system level indicators
##' @author Beth Fulton
plot_INDISEAS_indices <- function(df, OutDir) {
  
  # This assumes that size classes have been loaded as part of the Species_ID file
  KeepThese <- c("TEP","Piscivorous","Planktivorous","Invert")
  KeepFish <- c("Piscivorous","Planktivorous")
  minYear <- min(df$Year)
  maxYear <- max(df$Year)
  nYear <- maxYear - minYear + 1
  YearNum <- as.data.frame(seq(minYear, maxYear, by = 1))
  
  dfCalc <- data.frame(matrix(0, ncol = 9, nrow = nYear))
  rownames(dfCalc)[1:nYear] <- YearNum[,1]
  colnames(dfCalc)[1:9] <- c("Year","meanLengthBio","meanLengthCatch","meanTLBio","meanTLCatch","meanAgeBio","meanAgeCatch","LFIbio","LFIcatch")
  
  # Combine dataframes
  dftmp <- merge(df,tl_species, by = c("Year", "Group"))
  dftmp <- subset (dftmp, select = -c(species_id.y,PelDemID.y:MaxAge_yr.y))
  colnames(dftmp)[7:10] <- c("PelDemID","PISCZOOPL","Linf_cm","MaxAge_yr")
  colnames(dftmp)[3] <- c("species_id")
  
  counter <- 0
  for (iYear in minYear:maxYear) {
    counter <- counter + 1
    
    # Do biomass calculations first
    dfExtract <- dplyr::filter(dftmp, Year == iYear & PISCZOOPL %in% KeepThese)
    dfExtract$SizeCalcB <- dfExtract$biomass_tonnes * dfExtract$Linf_cm
    dfExtract$TLB <- dfExtract$biomass_tonnes * dfExtract$trophic_level
    dfExtract$AgeB <- dfExtract$biomass_tonnes * dfExtract$MaxAge_yr
    meanL <- sum(dfExtract$SizeCalcB) / sum(dfExtract$biomass_tonnes)
    meanTL <- sum(dfExtract$TLB) / sum(dfExtract$biomass_tonnes)
    meanAge <- sum(dfExtract$AgeB) / sum(dfExtract$biomass_tonnes)
    dfCalc$Year[counter] <- iYear
    dfCalc$meanLengthBio[counter] <- meanL
    dfCalc$meanTLBio[counter] <- meanTL 
    dfCalc$meanAgeBio[counter] <- meanAge 
    
    # Do LFI for biomass
    dfExtract <- dplyr::filter(dfExtract, PISCZOOPL %in% KeepFish)
    totbio <- sum(dfExtract$biomass_tonnes)
    dfExtract <- dplyr::filter(dfExtract, Linf_cm > 40.0)
    totLFI <- sum(dfExtract$biomass_tonnes)
    LFI <- totLFI / totbio
    dfCalc$LFIbio[counter] <- LFI 
    
    # Now do it based on landings by only including those species that are fished
    dfExtract <- dplyr::filter(dftmp, Year == iYear)
    ndim <- dim(dfExtract)
    num_row <- ndim[1]
    dfExtract$Include <- 0
    for (i in 1:num_row ) {
      if(dfExtract$landings_tonnes[i] > 0) {
        dfExtract$Include[i] <- 1
      } else {
        dfExtract$Include[i] <- 0
      }
    }
    dfExtract <- dplyr::filter(dfExtract, Include > 0)
    dfExtract$SizeCalcL <- dfExtract$landings_tonnes * dfExtract$Linf_cm
    dfExtract$TLCatch <- dfExtract$landings_tonnes * dfExtract$trophic_level
    dfExtract$AgeC <- dfExtract$landings_tonnes * dfExtract$MaxAge_yr
    meanL <- sum(dfExtract$SizeCalcL) / sum(dfExtract$landings_tonnes)
    meanTL <- sum(dfExtract$TLCatch) / sum(dfExtract$landings_tonnes)
    meanAge <- sum(dfExtract$AgeC) / sum(dfExtract$landings_tonnes)
    dfCalc$meanLengthCatch[counter] <- meanL
    dfCalc$meanTLCatch[counter] <- meanTL 
    dfCalc$meanAgeCatch[counter] <- meanAge 
    
    # Do LFI for catch
    dfExtract <- dplyr::filter(dfExtract, PISCZOOPL %in% KeepFish)
    totbio <- sum(dfExtract$landings_tonnes)
    dfExtract <- dplyr::filter(dfExtract, Linf_cm > 40.0)
    totLFI <- sum(dfExtract$landings_tonnes)
    LFI <- totLFI / totbio
    dfCalc$LFIcatch[counter] <- LFI 
    
  }
  
  # save results file
  outName <- paste(OutDir,"/Aggregate_indices.csv",sep="")
  write.csv(dfCalc, file = outName, row.names = FALSE)
  
  # Plot the results
  plot_indx <- paste(OutDir,"/mean_size_bio_indx.png",sep="")
  pLFI <- ggplot(dfCalc, aes (x = Year, y = meanLengthBio), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "Mean size time series - based on macrofauna biomass",
      x = "Year",
      y = "Size (cm)"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pLFI)
  
  ggsave(pLFI, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_TL_bio_indx.png",sep="")
  pTL <- ggplot(dfCalc, aes (x = Year, y = meanTLBio), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "Mean trophic level series - based on macrofauna biomass",
      x = "Year",
      y = "Trophic level"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pTL)
  
  ggsave(pTL, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_Age_bio_indx.png",sep="")
  pAB <- ggplot(dfCalc, aes (x = Year, y = meanAgeBio), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "Mean age time series - based on macrofauna biomass",
      x = "Year",
      y = "Longevity (years)"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pAB)
  
  ggsave(pAB, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_LFI_bio_indx.png",sep="")
  pLFIb <- ggplot(dfCalc, aes (x = Year, y = LFIbio), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "LFI time series - based on macrofauna biomass",
      x = "Year",
      y = "LFI"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pLFIb)
  
  ggsave(pLFIb, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_size_landings_indx.png",sep="")
  pSL <- ggplot(dfCalc, aes (x = Year, y = meanLengthCatch), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "Mean size time series - based on landings",
      x = "Year",
      y = "Size (cm)"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pSL)
  
  ggsave(pSL, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_TL_landings_indx.png",sep="")
  pmTL <- ggplot(dfCalc, aes (x = Year, y = meanTLCatch), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "Mean trophic level series - based on landings",
      x = "Year",
      y = "Trophic level"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pmTL)
  
  ggsave(pmTL, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_Age_landings_indx.png",sep="")
  pAL <- ggplot(dfCalc, aes (x = Year, y = meanAgeCatch), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "Mean age time series - based on landings",
      x = "Year",
      y = "Longevity (years)"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pAL)
  
  ggsave(pAL, file=plot_indx)
  
  plot_indx <- paste(OutDir,"/mean_LFI_landings_indx.png",sep="")
  pLFIL <- ggplot(dfCalc, aes (x = Year, y = LFIcatch), ) +
    geom_line(color = aggtscolor, linewidth = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "LFI time series - based on landings",
      x = "Year",
      y = "LFI"
    ) + 
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=12,face="bold"), 
          strip.text = element_text(face="bold", size=10))
  
  print(pLFIL)
  
  ggsave(pLFIL, file=plot_indx)
  
}

# Function to plot Link and Watson 2019 system level indicators
# Ryther - catch per unit area (divide total catch by area of the model to get this value)
# Fogarty - total catch / total primary production (units are per mil - o/oo where 
# 0.001 is 1 o/oo)
# Friedland - use catch: Chlorophyll a (basically doing Fogarty but use Chl a 
# instead if primary production is not available. For EwE calcaulte off primary 
# production but also plot it up based on assuptions of 
# biomass:C -> c: Chla -> final calculations
##' @title Calculate system level indicators defined by Link and Watson 2019
##' @author Beth Fulton
plot_Link_Watson <- function(df, mortality, biomass, id2, chl_from_file, OutDir) {
  
  # Ryther 
  dfcat <- df %>%
    group_by(Year) %>%
    dplyr::summarise(tot_catch = sum(landings_tonnes))
  dfcat$Ryther <- dfcat$tot_catch / area_model
  
  # Fogarty 
  # Start by getting the production values
  production <- mortality * biomass
  production$Year <- biomass$Year
  dfProd <- as.data.frame(production) %>% 
    pivot_longer(-Year, names_to = "species_id", values_to = "production")
  dftmp <- full_join(dfProd, id2) # Reuse old temporary df
  df1tmp <- full_join(dftmp, tl) # Reuse old temporary df
  dfProd <- dplyr::filter(df1tmp, trophic_level == 1.0) # Intentional reuse
  dfProd <- dplyr::filter(dfProd, IsDetritus == 0)
  
  dftmp <- dfProd %>%
    group_by(Year) %>%
    dplyr::summarise(tot_prod = sum(production))  # Already per km2 as never multiplied production by area (and mortality and biomass is per km2)
  dfLW <- full_join(dfcat, dftmp)
  dfLW$Fogarty <- 1000.0 * dfLW$Ryther / dfLW$tot_prod
  
  # Friedland 
  # Calculate Chl a
  if (chl_from_file < 1) {
    X_CN <- 5.7 # Redfield assumption
    X_CHLN <- 7.0 # Assumed ratio of N:Chl a
    dftmp <- as.data.frame(biomass) %>% 
      pivot_longer(-Year, names_to = "species_id", values_to = "biomass_per_area")
    dftmp1 <- full_join(dftmp, id2) # Reuse old temporary df
    dftmp <- full_join(dftmp1, tl) # Reuse old temporary df
    dfProd <- dplyr::filter(dftmp, trophic_level == 1.0) # Intentional reuse
    dfProd <- dplyr::filter(dfProd, IsDetritus == 0)
    dfProd$Chl <- 1000000000.0 * ((dfProd$biomass_per_area / (X_CN * 20.0)) / X_CHLN)    # From wet weight to N to Chl a
    thislinetype <- "dashed"
  } else {
    dfProd <- dfChl 
    thislinetype <- "solid"
  }
  
  dftmp <- dfProd %>%
    group_by(Year) %>%
    dplyr::summarise(tot_chl = sum(Chl))  # Already per km2 
  dfLW <- full_join(dfLW, dftmp)
  dfLW$Friedland <- 1000.0 * dfLW$Ryther / dfLW$tot_prod
  
  # Plots
  outaplot <- paste(OutDir,"/Ryther_indx.png",sep="")
  check_name <- c("Ryther")
  dfLWplot <- dfLW %>% pivot_longer(-Year, names_to = "Index")
  dfLW <- dplyr::filter(dfLWplot, Index %in% check_name) 
  pRy <- ggplot(dfLW, aes(x=Year, y=value)) +
    geom_line(aes(linetype = Index), linewidth = 2, colour="black") + 
    geom_hline(yintercept = 1, colour = "gray47", linetype="dashed", linewidth=0.5) +
    geom_hline(yintercept = 3, colour = "orange", linetype="dashed", linewidth=0.5) +
    geom_hline(yintercept = 5, colour = "red", linetype="dashed", linewidth=0.5) +
    theme_bw() + theme(axis.title=element_text(size=14,face="bold")) +
    theme(legend.position = "none") + ggtitle("Ryther index") +
    labs(x="Year", y="tonnes km-2 year-1") +
    scale_linetype_manual(values=c("solid", thislinetype, "solid"))
  
  print(pRy)
  
  ggsave(pRy, file=outaplot)
  
  outaplot <- paste(OutDir,"/Friedland_ratio.png",sep="")
  check_name <- c("Friedland")
  dfLW <- dplyr::filter(dfLWplot, Index %in% check_name) 
  pFR <- ggplot(dfLW, aes(x=Year, y=value)) +
    geom_line(aes(linetype = Index), linewidth = 2, colour="darkcyan") + 
    geom_hline(yintercept = 1, colour = "gray47", linetype="dashed", linewidth=0.5) +
    theme_bw() + theme(axis.title=element_text(size=14,face="bold")) +
    theme(legend.position = "none") + ggtitle("Friedland Ratio") +
    labs(x="Year", y="0/00") +
    scale_linetype_manual(values=c("solid", thislinetype, "solid"))
  
  print(pFR)
  
  ggsave(pFR, file=outaplot)
  
  outaplot <- paste(OutDir,"/Fogarty_ratio.png",sep="")
  check_name <- c("Fogarty")
  dfLW <- dplyr::filter(dfLWplot, Index %in% check_name) 
  pFog <- ggplot(dfLW, aes(x=Year, y=value)) +
    geom_line(aes(linetype = Index), linewidth = 2, colour="darkgoldenrod1") + 
    geom_hline(yintercept = 1, colour = "gray47", linetype="dashed", linewidth=0.5) +
    theme_bw() + theme(axis.title=element_text(size=14,face="bold")) +
    theme(legend.position = "none") + ggtitle("Fogarty Ratio") +
    labs(x="Year", y="") +
    scale_linetype_manual(values=c("solid", thislinetype, "solid"))
  
  print(pFog)
  
  ggsave(pFog, file=outaplot)
  
}

# Cumulative biomass (survey or catch) vs Trophic level – steepness of 
# the curve and height of the asymptote are indicators of overall system state
##' @title Calculate and plot cumaultive biomass vs TL plots
##' @author Beth Fulton
plot_CumB_TL <- function(df, tl) {
  dfcat <- full_join(df, tl) # Reuse old temporary df
  
  # Get cumulative biomass per TL increment - Using Biomass
  #minYear <- min(dfcat$Year)
  minYear <- yr_ignore + 1
  maxYear <- max(dfcat$Year)
  nYear <- maxYear - minYear + 1
  YN <- seq(minYear,maxYear,1)
  YearNum <- as.data.frame(YN)
  
  cumRes <- data.frame(matrix(0, ncol = 8, nrow = nYear))
  rownames(cumRes)[1:nYear] <- YearNum[,1]
  colnames(cumRes)[1:8] <- c("Year","infTL","inflB","steepness","b1","b2","c","d")
  
  dfcreated <- 0
  
  for (iYear in minYear:maxYear) {
    thisYear <- iYear - minYear + 1
    
    # Filter by year
    dfCD <- filter(dfcat, Year == iYear)
    
    # sort by TL ascending
    dfCDsorted <- dfCD[order(dfCD$trophic_level), ]
    
    # calculate cumulative catch 
    dfCDsorted$cumB <- cumsum(dfCDsorted$landings_tonnes)
    
    # TotalBiomass
    TotB <- max(dfCDsorted$cumB)
    
    # Now make column PropCumB
    dfCDsorted$PropCumB <- dfCDsorted$cumB / TotB
    
    # Fit five-parameter regression - using drc
    fm1 <- drm(PropCumB~trophic_level, data=dfCDsorted, fct=baro5())
    dfCDsorted$curvefit = predict(fm1)
    
    # Store results for sorted materials
    if (dfcreated == 0) {
      dfBioResults <- dfCDsorted
      dfcreated <- 1
    } else {
      dfBioResults <-rbind(dfBioResults, dfCDsorted)
    }
    
    #plot(fm1)
    #ggplot(dfCDsorted, aes(x = TL)) +
    #  geom_point(aes(y = cumB), color ="deepskyblue3") +
    #  geom_point(aes(y = cumC), color ="darkorange") +
    #  geom_line(aes(y = curvefit), size = 1)
    
    # Calculate - B inflection and TL inflection (only do it for catch if non-zero)
    # This involves finding when the second derivative of the curve = 0
    # However the "e" coeffient given by fm1 is the inflection TL so simply read off the 
    # inflection point B using that TL value
    infTL <- fm1$fit$par[5]
    b1 <- fm1$fit$par[1]
    b2 <- fm1$fit$par[2]
    c <- fm1$fit$par[3]
    d <- fm1$fit$par[4]
    e <- infTL
    
    baro5eq = expression(c + ((d - c)/(1 + (1/(1 + exp((2 * b1 * b2/(b1 + b2)) * (log(x) - log(e))))) * (exp(b1 * (log(x) - log(e)))) + (1 - (1/(1 + exp((2 * b1 * b2/(b1 + b2)) * (log(x) - log(e)))))) * (exp(b2 * (log(x) - log(e)))))))
    x <- infTL
    infB <- eval(baro5eq)
    
    # Calculate - Steepness (the first derivative at the inflection point)
    dx2x <- D(baro5eq,"x")
    steep <- eval(dx2x)
    
    # Store indicator results
    cumRes[thisYear,1] <- iYear
    cumRes[thisYear,2] <- infTL
    cumRes[thisYear,3] <- infB
    cumRes[thisYear,4] <- steep
    cumRes[thisYear,5] <- b1
    cumRes[thisYear,6] <- b2
    cumRes[thisYear,7] <- c
    cumRes[thisYear,8] <- d
    
  }
  
  
  pCB0 <- ggplot(dfBioResults, aes(x = trophic_level, y = curvefit)) +
    geom_point(aes(colour = factor(Year))) +
    geom_line(aes(colour = factor(Year))) +
    scale_color_viridis(option = "D", discrete = TRUE) +
    theme(text = element_text(size = textsize)) 
  
  # Plot 6 panels - the 3 indicators through time and 
  # Reproduce plots from Libralato et al 2019 
  # Steep vs TL
  # B vs TL
  # Steep vs B
  
  # 1. infTL through time
  pCB1 <- ggplot(cumRes, aes(x = Year, y = infTL)) +
    geom_line(color ="deepskyblue3", linewidth = linesize) +
    theme(text = element_text(size = textsize)) 
  
  # 2. infB through time
  pCB2 <- ggplot(cumRes, aes(x = Year, y = inflB)) +
    geom_line(color ="darkorange", linewidth = linesize) +
    theme(text = element_text(size = textsize)) 
  
  # 3. steep through time
  pCB3 <- ggplot(cumRes, aes(x = Year, y = steepness)) +
    geom_line(color ="grey44", linewidth = linesize) +
    theme(text = element_text(size = textsize)) 
  
  # 4. Steep vs TL
  pCB4 <- ggplot(cumRes, aes(x = infTL, y = steepness)) +
    geom_point() +
    theme(text = element_text(size = textsize)) 
  
  # 5. B vs TL 
  pCB5 <- ggplot(cumRes, aes(x = infTL, y = inflB)) +
    geom_point() +
    theme(text = element_text(size = textsize)) 
  
  # 6. Steep vs B 
  pCB6 <- ggplot(cumRes, aes(x = inflB, y = steepness)) +
    geom_point() +
    theme(text = element_text(size = textsize)) 
  
  # Plot Catch through time too
  nSP <- length(unique(dfcat$Group))
  colPalette <- get_palette(c("#00AFBB", "#E7B800", "#FC4E07"), nSP)
  
  pCBC <- ggplot(dfcat, aes(x = Year, y = landings_tonnes, fill = Group, order = Group)) +
    geom_bar(colour = "black", stat="identity", size = 0.1) +
    scale_fill_manual(values = colPalette) +
    labs(x="Years", y = "Catch") +  theme(legend.position = "none", text = element_text(size = 12))
  
  # Plot Biomass
  pCBB <- ggplot(dfcat, aes(x = Year, y = biomass_tonnes, fill = Group, order = Group)) +
    geom_bar(colour = "black", stat="identity", size = 0.1) +
    scale_fill_manual(values = colPalette) +
    labs(x="Years", y = "Biomass") +  theme(legend.position = "none", text = element_text(size = 12))
  
  # Final panel plot
  outPlotName <- paste(OutDir,"/Cumulative-plots-Catch.png",sep="")
  
  # Arrange plots using arrangeGrob
  # returns a gtable (gt)
  #gt <- arrangeGrob(p0,   # Big # cumBio vs TL
  #                  pC, pB, p1, p2, p3, p4, p5, p6, # All other plots
  #                  ncol = 2, nrow = 6, 
  #                  layout_matrix = rbind(c(1,1), c(2,2), c(3,3), c(4,5), c(6,7), c(8,9)))
  # Add labels to the arranged plots
  gt <- arrangeGrob(pCB0,   # Big # cumBio vs TL
                    pCBC, pCB1, pCB2, pCB3, pCB4, pCB5, pCB6, # All other plots
                    ncol = 2, nrow = 6, 
                    layout_matrix = rbind(c(1,1), c(2,2), c(3,4), c(5,6), c(7,8)))
  # Add labels to the arranged plots
  pCBfinal <- as_ggplot(gt)                                 # transform to a ggplot
  #+ draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H")) # Add labels
  
  ggsave(pCBfinal, file=outPlotName)
  
  return(pCBfinal)
  
  
}

########################### NETWORK INDICATOR FUNCTIONS ###########################

# Calculating NOAA's FSSU. Look at Blim (overfished), Btarg and F vs Ftarg (overfishing) 
# Also Count of species in each state
# NOAA Index Scoring Methodology (they do it each quarter, here it is annual)
# NOAA Fisheries calculates an Fish Stock Sustainability Index score, 
# incorporating information from new stock assessments and stock status 
# determinations. The index is calculated on a 1,000 point scale using 
# the following methodology:
#
#  Step 1: Assign weighted criteria points for each stock based on the following:
#  Criteria Criteria Points:
#    1. "Overfished" status is known. 0.5 
#    2. "Overfishing" status is known. 0.5 
#    3. Overfishing is not occurring (for stocks with known "overfishing" status). 1.0 
#    4. Stock biomass is above the "overfished" level defined for the stock. 1.0 
#    5. Stock biomass is at or above 80% of the biomass that produces Btarget 1.0 
#
#  Step 2: Calculate the sum of criteria points for all index stocks. 
#
#  Step 3: Calculate maximum criteria points possible: multiply number of 
# index stocks x maximum criteria points per stock (4 points). 
#
#  Step 4: Calculate a raw total point score: divide sum of criteria 
#  points / maximum criteria points possible.
#
#  Step 5: Convert raw total point score to a 1,000 point scale: total raw point score*1,000.
##' @title Calculate NOAA FSSI
##' @return dfFSSI - the FSSI values per species
##' @author Beth Fulton
calculation_FSSI <- function(df, dfRef) {
  # Reload reference biomasses - just in case
  dftmp <- merge(df, dfRef, by=c("species_id")) # Intentional reuse
  dfFSSI <- dftmp %>% dplyr::filter(IsDetritus == 0)
  colnames(dfFSSI)[2] <- "Year"
  dfFSSI$Rel_B <- dfFSSI$biomass_tonnes / dfFSSI$RefBo
  dfFSSI$U<- dfFSSI$catch_tonnes / dfFSSI$biomass_tonnes
  dfFSSI$Utarg <- 0.5 * dfFSSI$RefM
  dfFSSI$Blim <- ifelse(dfFSSI$Class_Code < 7, 0.2, 0.4)
  dfFSSI$Btarg <- ifelse(dfFSSI$Class_Code == 1, 0.4, 
                         (ifelse(dfFSSI$Class_Code == 2, 0.6, 
                                 (ifelse(dfFSSI$Class_Code == 3, 0.4, 
                                         (ifelse(dfFSSI$Class_Code == 4, 0.4, 
                                                 (ifelse(dfFSSI$Class_Code == 5, 0.48, 
                                                         (ifelse(dfFSSI$Class_Code == 6, 0.5, 0.7)))))))))))
  dfFSSI$BtargB <- dfFSSI$Btarg * 0.8
  
  dfFSSI$B_Score1 <- 0.5  # Needs to be replaced with meaningful value when used with observational data - right now assume status vs Blim, Btarg above is sufficient
  dfFSSI$B_Score2 <- 0.5 # Needs to be replaced with meaningful value when used with observational data -- right now assume status vs U = 0.5M is sufficient test
  dfFSSI$B_Score3 <- ifelse(dfFSSI$U < dfFSSI$Utarg, 1, 0)
  dfFSSI$B_Score4 <- ifelse(dfFSSI$Rel_B > dfFSSI$Btarg, 1, 0)
  dfFSSI$B_Score5 <- ifelse(dfFSSI$Rel_B > dfFSSI$BtargB, 1, 0)
  dfFSSI$ScoreStatus <- ifelse(dfFSSI$Rel_B > dfFSSI$Btarg, "L", (ifelse(dfFSSI$Rel_B <= dfFSSI$Blim, "F", "A")))
  
  
  # Assumes Thresholds
  # 
  #  CLass    ID    Blim  Btarg
  #
  # Vulnerable  1   0.2   0.4  
  # Habitat     2   0.2   0.6
  # Byproduct   3   0.2   0.4
  # Bycatch     4   0.2   0.4
  # Target      5   0.2   0.48
  # Robust      6   0.2   0.5
  # Hub         7   0.4   0.7
  
  nSP <- length(unique(dfFSSI$Group))
  MaxScore <- 4 * nSP
  dfAggFSSI <- dfFSSI%>%
    group_by(Year) %>%
    dplyr::summarise(FSSI1 = sum(B_Score1), FSSI2 = sum(B_Score2), FSSI3 = sum(B_Score3), FSSI4 = sum(B_Score4), FSSI5 = sum(B_Score5))
  dfAggFSSI$Total_FSSI <- dfAggFSSI$FSSI1 + dfAggFSSI$FSSI2 + dfAggFSSI$FSSI3 + dfAggFSSI$FSSI4 + dfAggFSSI$FSSI5
  dfAggFSSI$Rel_FSSI <- dfAggFSSI$Total_FSSI / MaxScore
  minFSSI <- round(min(dfAggFSSI$Rel_FSSI), digits = 2)
  maxFSSI <- round(max(dfAggFSSI$Rel_FSSI), digits = 2)
  notemsg <- paste("Minimum: ", minFSSI, " Maximum: ", maxFSSI, sep="")
  
  plot_indx <- paste(OutDir,"/FSSI_indx.png",sep="")
  pFS <- ggplot(dfAggFSSI, aes (x = Year, y = Rel_FSSI), ) +
    geom_line(color = aggtscolor, size = 2) +
    theme(legend.position = "none")+
    expand_limits(y=0) +
    labs(
      title = "FFSI index - proportion of max possible score",
      subtitle = notemsg,
      x = "Year",
      y = "Relative FSSI"
    ) + 
    theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=10), axis.title=element_text(size=12,face="bold"), strip.text = element_text(face="bold", size=10))
  
  print(pFS)
  
  ggsave(pFS, file=plot_indx)
  
  return(dfFSSI)
}

# Calculation of structural indcators - creating plots for each one
##' @title Calculate NOAA FSSI
##' @return dfFSSI - the FSSI values per species
##' @author Beth Fulton
calculate_struct_indices <- function(df, dfRef, mortality, biomass, df_consumption, df_species, id2) {
  
  # Calculate the Greenband
  dfScored <- do_greenband(dfRef, mortality, biomass, id2)
  
  # Calculate Gao resilience
  Gaoindx <- do_Gao(df_consumption, df_species, df)
  
  return(list(dfScored = dfScored, Gaoindx = Gaoindx))
  
}


# Calcualtion of the green band Follows the approach of Heath - where the allowed
# catch is proportional to production. This is determined by using the predation profile
# of an unfished system to identify catch levels that should not put disproportionate pressure
# on the foodweb
##' @title Calculate Greenband
##' @author Beth Fulton
do_greenband <-function(dfRef, mortality, biomass, id2) {
  
  RefVal <- dfRef
  RefVal$RefP <- RefVal$RefB * RefVal$RefM
  idRef <- merge(id2, RefVal, by=c("species_id"))
  idRef <- idRef %>% dplyr::filter(IsDetritus < 1)
  
  # Estimate the initial slope
  idRef$logP <- log10(idRef$RefP)
  idRef$logB <- log10(idRef$RefB)
  linreg <- lm(formula = logB ~ logP, data = idRef)
  matrix_coef <- summary(linreg)$coefficients 
  reg_slope <- round(matrix_coef[2,1],3)
  reg_const <- round(matrix_coef[1,1],3)
  band_slope <- reg_slope + 1.0   # This is the slope of the grend band for the rest of the work
  lm_equation <- paste("y = ", reg_const, " + ", reg_slope, " x",sep="")
  this_r2 <- round(summary(linreg)$r.squared,3)
  lm_R2 <- paste("R2 = ",this_r2,sep="")
  
  ymin <- 1e-8
  xmin <- 1e-6
  ymax <- 1e+5
  xmax <- 1e+5
  lymin <- log10(1e-8)
  lxmin <- log10(1e-6)
  lymax <- log10(1e+5)
  lxmax <- log10(1e+5)
  
  outaplot <- paste(OutDir,"/Green_Band_Reference_Plot.png",sep="")
  
  p_ref <- ggplot(data = idRef, aes(x = RefP, y = RefB)) + 
    geom_point() +
    geom_smooth(method='lm') +
    scale_x_log10(limits=c(xmin, xmax)) + scale_y_log10(limits=c(ymin, ymax)) +
    labs(x="Production", y = "Biomass") +
    annotate("text", x=1e-4, y=3e+4, label= lm_equation) + 
    annotate("text", x = 1e-4, y=1e+4, label = lm_R2) +
    theme(axis.text=element_text(size=textsize,face="bold"), 
          axis.title=element_text(size=axistitle,face="bold"))
  
  print(p_ref)
  
  ggsave(p_ref, file=outaplot)
  
  # Get production once fished
  production <- mortality * biomass
  production$Year <- biomass$Year
  dfProd <- as.data.frame(production) %>% 
    pivot_longer(-Year, names_to = "species_id", values_to = "Production")
  dftmp <- full_join(dfProd, id2) # Reuse old temporary df
  dfBand <- full_join(dftmp, df) 
  dfBand <- dfBand %>% dplyr::filter(Class_Code > 0)
  
  # Rate in comparison to green band
  # Equations for key lines
  # Max_allowable  C = 0.5 * P
  # Max_green_band  C = P ^ band_slope
  # Min_green_band  C = C_Max_green_band * 0.01
  # Realised Max_green_band  C = min(Max_allowable, Max_green_band)
  # Realised Min_green_band  C = min(Max_allowable, Min_green_band)
  # Light 1 - L
  # Acceptable 2 - A
  # Fail 3 - F
  
  dfBand$MaxAllow <- 0.5 * dfBand$Production
  dfBand$MaxGreenBand <- dfBand$Production ^ band_slope
  dfBand$MinGreenBand <- dfBand$MaxGreenBand * 0.01
  dfBand$RealMax <- ifelse(dfBand$MaxAllow < dfBand$MaxGreenBand, dfBand$MaxAllow , dfBand$MaxGreenBand)
  dfBand$RealMin <- ifelse(dfBand$MaxAllow < dfBand$MinGreenBand, dfBand$MaxAllow , dfBand$MinGreenBand)
  dfBand$catch <- dfBand$catch_tonnes / area_model
  #dfBand$DistortScore <- ifelse(dfBand$catch_tonnes < dfBand$RealMin, "L", (ifelse(dfBand$catch_tonnes >= dfBand$RealMax, "F", "A")))
  dfBand$DistortScore <- ifelse(dfBand$catch < dfBand$RealMin, "L", (ifelse(dfBand$catch >= dfBand$RealMax, "F", "A")))
  
  Ofname <- paste(OutDir,"/dfBand.csv",sep="")
  write.csv(dfBand, file = Ofname, row.names = FALSE)
  
  dfScored <- merge(dfBand, dfRef, by=c("species_id"))
  names(dfScored)[names(dfScored) == "species_id.x"] <- "species_id"
  names(dfScored)[names(dfScored) == "Classification.x"] <- "Classification"
  names(dfScored)[names(dfScored) == "Class_Code.x"] <- "Class_Code"
  
  dfScored$Rel_B <- dfScored$biomass_tonnes / dfScored$RefBo
  dfScored$Min_B <- ifelse(dfScored$Class_Code == 1, 0.5, 
                           (ifelse(dfScored$Class_Code == 2, 0.3, 
                                   (ifelse(dfScored$Class_Code == 3, 0.35, 
                                           (ifelse(dfScored$Class_Code == 4, 0.2, 
                                                   (ifelse(dfScored$Class_Code == 5, 0.4, 
                                                           (ifelse(dfScored$Class_Code == 6, 0.4, 0.6)))))))))))
  dfScored$Max_B <- ifelse(dfScored$Class_Code == 1, 0.7, 
                           (ifelse(dfScored$Class_Code == 2, 0.6, 
                                   (ifelse(dfScored$Class_Code == 3, 0.6, 
                                           (ifelse(dfScored$Class_Code == 4, 0.6, 
                                                   (ifelse(dfScored$Class_Code == 5, 0.5, 
                                                           (ifelse(dfScored$Class_Code == 6, 0.5, 0.7)))))))))))
  dfScored$B_Score <- ifelse(dfScored$Rel_B > dfScored$Max_B, "L", (ifelse(dfScored$Rel_B <= dfScored$Min_B, "F", "A")))
  
  # Assumes Thresholds
  # 
  #  CLass    ID    Min   Max
  #
  # Vulnerable  1   0.5   0.7  
  # Habitat     2   0.3   0.6
  # Byproduct   3   0.35  0.6
  # Bycatch     4   0.2   0.6
  # Target      5   0.4   0.5
  # Robust      6   0.4   0.5
  # Hub         7   0.6   0.7
  
  # Aggregate Results
  dfAggScore_Distort <- dfScored %>% dplyr::count(Year, DistortScore)
  dfAggScore_B <- dfScored %>% dplyr::count(Year, B_Score)
  
  dfAggScore_Distort <- dfAggScore_Distort %>% dplyr::filter(!is.na(DistortScore))
  dfAggScore_B <- dfAggScore_B %>% dplyr::filter(!is.na(B_Score))
  
  # Plot Scores through time - as bar plot
  outPlotName <- paste(OutDir,"/GreenBand-aggregate-score.png",sep="")
  pA1 <- ggplot(data = dfAggScore_Distort, aes(x = Year, y = n, fill = DistortScore)) + geom_bar(colour = "black", stat="identity", size = 0.1) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme_bw() +
    labs(x="Years", y = "Ratings") +
    theme(axis.text=element_text(size=textsize,face="bold"), axis.title=element_text(size=axistitle,face="bold"))
  print(pA1)
  ggsave(pA1, file=outPlotName) 
  
  # Plot Scores through time - as bar plot
  outPlotName <- paste(OutDir,"/GreenBand-aggregate-scoreB.png",sep="")
  pA2 <- ggplot(data = dfAggScore_B, aes(x = Year, y = n, fill = B_Score)) + geom_bar(colour = "black", stat="identity", size = 0.1) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme_bw() +
    labs(x="Years", y = "Ratings") +
    theme(axis.text=element_text(size=textsize,face="bold"), axis.title=element_text(size=axistitle,face="bold"))
  print(pA2)
  ggsave(pA2, file=outPlotName) 
  
  Bounds <- data.frame(matrix(0, ncol = 6, nrow = 10))
  colnames(Bounds)[1:6] <- c("P","MaxAllow","MaxGreenBand","MinGreenBand","MinC","MaxC")
  listP <- c(1.00E-19, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)
  Bounds$P <- listP
  Bounds$MaxAllow <- 0.5 * Bounds$P
  Bounds$MaxGreenBand <- Bounds$P ^ band_slope
  Bounds$MinGreenBand <- Bounds$MaxGreenBand * 0.01
  Bounds$MaxC <- ifelse(Bounds$MaxAllow < Bounds$MaxGreenBand, Bounds$MaxAllow, Bounds$MaxGreenBand)
  Bounds$MinC <- ifelse(Bounds$MaxAllow < Bounds$MinGreenBand, Bounds$MaxAllow, Bounds$MinGreenBand)
  BoundsMelt <- reshape2::melt(Bounds,"P") 
  names(BoundsMelt)[2] <- 'Variables'
  names(BoundsMelt)[3] <- 'CValue'
  
  write.csv(BoundsMelt, 
            file = paste0(OutDir, "/Bounds.csv"), 
            row.names = FALSE)
  write.csv(dfScored, 
            file = paste0(OutDir, "/GreenBand.csv"), 
            row.names = FALSE)
  # dfScored <- read.csv(
  #   paste0(OutDir, "/GreenBand.csv"), 
  #   header = TRUE, check.names = FALSE)
  
  
  # Per km2 results
  nSpCode <- max(dfScored$Class_Code)
  for (i in 1:nSpCode) {
    dfThisCode <- dfScored %>% dplyr::filter(Class_Code == i)
    outPlotName <- paste(OutDir,"/GreenBand-",dfThisCode$Classification[1],"_Per_Area.png",sep="")
    
    dfMin <- BoundsMelt %>% dplyr::filter(Variables == "MinC")
    dfMax <- BoundsMelt %>% dplyr::filter(Variables == "MaxC")
    
    #ymin <- min(dfThisCode$C) * 0.1
    #xmin <- min(dfThisCode$P) * 0.1
    #ymax <- max(dfThisCode$C) * 10.0
    #xmax <- max(dfThisCode$P) * 10.0
    
    pA3 <- ggplot(data = dfThisCode, aes(x = Production, y = catch)) + 
      geom_point(data = dfThisCode, aes(size = Year, color = Group)) +
      geom_line(data = dfThisCode, aes(color = Group)) +
      geom_line(data = dfMin, aes(x = P, y = CValue), linetype = "dashed", color = "black") +
      geom_line(data = dfMax, aes(x = P, y = CValue), linetype = "dashed", color = "black") +
      scale_size_continuous(range = c(1, 3)) +
      scale_x_log10(limits=c(xmin, xmax)) + scale_y_log10(limits=c(ymin, ymax)) +
      labs(x="Production", y = "Catch") +
      theme(axis.text=element_text(size=textsize,face="bold"), axis.title=element_text(size=axistitle,face="bold"))
    
    print(pA3)
    
    ggsave(pA3, file=outPlotName)
    
  }
  
  ## Total (so area corrected) results
  ymin <- 1e-5 
  xmin <- 1e-5
  ymax <- 1e+10
  xmax <- 1e+10
  
  for (i in 1:nSpCode) {
    dfThisCode <- dfScored %>% dplyr::filter(Class_Code == i)
    outPlotName <- paste(OutDir,"/GreenBand-",dfThisCode$Classification[1],".png",sep="")
    
    dfThisCode$P<- dfThisCode$Production * area_model
    
    dfMin <- BoundsMelt %>% dplyr::filter(Variables == "MinC")
    dfMax <- BoundsMelt %>% dplyr::filter(Variables == "MaxC")
    
    dfMin$Pa <- dfMin$P * area_model
    dfMin$Ca <- dfMin$CValue * area_model
    dfMax$Pa <- dfMax$P * area_model
    dfMax$Ca <- dfMax$CValue * area_model
    
    #ymin <- min(dfThisCode$C) * 0.1
    #xmin <- min(dfThisCode$P) * 0.1
    #ymax <- max(dfThisCode$C) * 10.0
    #xmax <- max(dfThisCode$P) * 10.0
    
    pA4 <- ggplot(data = dfThisCode, aes(x = P, y = catch_tonnes)) + 
      geom_point(data = dfThisCode, aes(size = Year, color = Group)) +
      geom_line(data = dfThisCode, aes(color = Group)) +
      geom_line(data = dfMin, aes(x = Pa, y = Ca), linetype = "dashed", color = "black") +
      geom_line(data = dfMax, aes(x = Pa, y = Ca), linetype = "dashed", color = "black") +
      scale_size_continuous(range = c(0.5, 2)) +
      scale_x_log10(limits=c(xmin, xmax)) + scale_y_log10(limits=c(ymin, ymax)) +
      labs(x="Production", y = "Catch") +
      guides(col = guide_legend(ncol = 2)) +
      theme(axis.text=element_text(size=textsize,face="bold"), axis.title=element_text(size=axistitle,face="bold"))
    
    print(pA4)
    
    ggsave(pA4, file=outPlotName)
    
  }
  
  outPlotName <- paste(OutDir,"/GreenBand-All.png",sep="")
  # dfThisCode <- read.csv("GBsnapshot.csv", header = T)
  dfThisCode <- dfThisCode %>% dplyr::filter(catch_tonnes > 0)
  dfMin <- BoundsMelt %>% dplyr::filter(Variables == "MinC")
  dfMax <- BoundsMelt %>% dplyr::filter(Variables == "MaxC")
  
  dfMin$Pa <- dfMin$P * area_model
  dfMin$Ca <- dfMin$CValue * area_model
  dfMax$Pa <- dfMax$P * area_model
  dfMax$Ca <- dfMax$CValue * area_model
  
  pA5 <- ggplot(data = dfThisCode, aes(x = Production, y = catch_tonnes)) + 
    geom_point(data = dfThisCode, aes(color = Group)) +
    geom_line(data = dfMax, aes(x = Pa, y = Ca), linetype = "dashed", color = "black") +
    geom_line(data = dfMin, aes(x = Pa, y = Ca), linetype = "dashed", color = "black") +
    scale_size_continuous(range = c(0.5, 2)) +
    scale_x_log10(limits=c(xmin, xmax)) + scale_y_log10(limits=c(ymin, ymax)) +
    labs(x="Production", y = "Catch") +
    guides(col = guide_legend(ncol = 2)) +
    theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=textsize,face="bold"))
  
  print(pA5)
  
  ggsave(pA5, file=outPlotName)
  
  return(dfScored)
  
}

# Calcualtion of Gao Reilsience - which uses the density of connections and degree of heterogeneity of
# flow through sub-webs to abalytically predict resilience.
##' @title Calculate Gao Resilience
##' @author Beth Fulton
do_Gao <- function(df_consumption, df_species, df) {
  
  # To calculate the degree in and out through time
  # Load base consumption matrix
  localIDs <- subset(df_species, select = -c(PelDemID, PISCZOOPL, Linf_cm, MaxAge_yr, IsDetritus, Classification, Class_Code, isHub))
  
  # Calculate Gao index for that value of <s> so can check on whether above or below resilience horizon
  c_coefft <- 0.019
  b_coefft <- 0.8931
  a_coefft <- 5.3384
  H_ref <- 5.32
  s_ref <- 6.97
  
  # Get network indices through time - assumes that ref is for year 1
  minYear <- min(df$Year)
  maxYear <- max(df$Year)
  nYear <- maxYear - minYear + 1
  nSP <- length(unique(df$species_id))
  Gaoindx <- data.frame(matrix(0, ncol = 7, nrow = nYear))
  colnames(Gaoindx)[1:7] <- c("Year","Gao_density","Gao_s","Gao_H","Gao_indx","resil_H","GaoScore")
  rownames(Gaoindx)[1:nYear] <- seq(minYear, maxYear, by=1)
  dfDeg <- data.frame(matrix(0, ncol = 3, nrow = nSP))
  colnames(dfDeg)[1:3] <- c("species_id","Din","Dout")
  
  ref_yr <- df %>% dplyr::filter(Year == minYear)
  names(ref_yr)[names(ref_yr) == "consumption"] <- "ref_consum"
  ref_yr$ref_consum <- ref_yr$ref_consum + 0.0000001
  
  startScore <- -1
  startScoreH <- 0
  for (iYear in minYear:maxYear) {
    thisYear <- iYear - minYear + 1
    extractdf <- df %>% dplyr::filter(Year == iYear)
    dftmp <- merge(extractdf, dfRef, by=c("species_id"))
    dftmp$rel_B <-dftmp$biomass_tonnes / dftmp$RefBo
    dftmp$rel_B[dftmp$rel_B < 0.00001] <- 0  # So ave an effect on flows of being highly depleted
    newQ <- df_consumption %>%
      left_join(dftmp %>% select(Group, rel_B), by = "Group") %>%
      mutate(across(-c(Group, rel_B), ~ . * rel_B)) %>%
      select(-rel_B)
    
    rownames(newQ) <- newQ$Group
    newQ <- newQ %>% select(-Group)
    dietData <- newQ * 0  ## empty matrix
    for(iSP in 1:nSP){  # Treat iSP as col_ID doing calculation on
      dfPred <- dftmp %>% dplyr::filter(species_id == iSP)
      bioPred <- dfPred$rel_B
      dietData[,iSP] <- newQ[ ,iSP] * bioPred
    }
    
    # Calculate Degree in and Degree out
    step1 <- reshape2::melt(as.matrix(dietData))
    names(step1)[names(step1) == "Var1"] <- "preyName"
    names(step1)[names(step1) == "Var2"] <- "predName"
    names(step1)[names(step1) == "value"] <- "propDiet"
    step1$predName <- gsub(".", " ", step1$predName, fixed=TRUE)  #Sort the issue with spaces in names
    
    # Replace names with IDs
    step2 <- step1 %>% dplyr::inner_join(localIDs,by=c("preyName" = "Group"))
    names(step2)[names(step2) == "ID"] <- "prey_ID"
    step2 <- step2 %>% dplyr::inner_join(localIDs,by=c("predName" = "Group"))
    names(step2)[names(step2) == "ID"] <- "pred_ID"
    nDiet <- step2[step2$propDiet != 0, ]
    
    # Create edges df
    if (use_weights > 0) {
      nl <- subset(nDiet, select = -c(preyName, predName))
      names(nl)[names(nl) == "prey_ID"] <- "V1"
      names(nl)[names(nl) == "pred_ID"] <- "V2"
      nl$na <- "FALSE"
      nl$V1 <- as.numeric(nl$V1)
      nl$V2 <- as.numeric(nl$V2)
      nl$na <- as.logical.factor(nl$na)
      nl$na <- FALSE
      nl$weight <- nl$propDiet
      nl <- subset(nl, select = -c(propDiet))
    } else {
      nl <- subset(nDiet, select = -c(preyName, predName, propDiet))
      names(nl)[names(nl) == "prey_ID"] <- "V1"
      names(nl)[names(nl) == "pred_ID"] <- "V2"
      nl$na <- "FALSE"
      nl$V1 <- as.numeric(nl$V1)
      nl$V2 <- as.numeric(nl$V2)
      nl$na <- as.logical.factor(nl$na)
      nl$na <- FALSE
    }
    
    nlv <- localIDs
    nlv$na <- "FALSE"
    nlv$na <- as.logical.factor(nlv$na)
    nlv$na <- FALSE
    nlv$vertex.names <- nlv$Group
    names(nlv)[names(nlv) == "species_ID"] <- "intergraph_id"
    nlv <- subset(nlv, select = -c(Group))
    
    g <- asIgraph(nl, vertices=nlv, directed=TRUE)
    set_vertex_attr(g, "label", value=nlv[,2])
    #edges <- get.edgelist(g)   # Can also be extracted using E(g)
    #nodes <- V(g)$vertex.names    # Name of vertices
    
    food_tidy <- as_tbl_graph(g, directed = TRUE)
    centRes   <- centrality(food_tidy)
    deg <- degree(g, mode = "all", normalized = FALSE)
    d_in  <- degree(g, mode = "in", normalized = FALSE)
    d_out <- degree(g, mode = "out", normalized = FALSE)
    
    # Store value
    dfDeg$species_id <- seq(1, nSP, by=1)
    dfDeg$Din <- d_in
    dfDeg$Dout <- d_out
    dfDeg$density <- deg
    
    # Calculate H
    GaoStep1 <- merge(dftmp, dfDeg, by=c("species_id")) # intentionally reused
    #GaoStep1$WgtDin <- GaoStep1$Din * GaoStep1$biomass_tonnes / area_model
    #GaoStep1$WgtDout <- GaoStep1$Dout * GaoStep1$biomass_tonnes / area_model
    GaoStep1$WgtDin <- GaoStep1$Din * GaoStep1$consumption
    GaoStep1$WgtDout <- GaoStep1$Dout * GaoStep1$consumption
    GaoStep1$WDIWDO <- GaoStep1$WgtDin * GaoStep1$WgtDout
    GaoStep3 <- merge(GaoStep1, ref_yr, by=c("species_id"))
    
    step1 <- as.data.frame(GaoStep1$WgtDin)
    colnames(step1)[1] <- "x"
    step2 <- as.data.frame(GaoStep1$WgtDout)
    colnames(step2)[1] <- "x"
    step1 <- rbind(step1, step2)
    #Gao_density <- mean(step1$x)
    step3 <- GaoStep3$density * GaoStep3$consumption / GaoStep3$ref_consum
    Gao_density <- mean(step3)
    Gao_meanWDIWDO <- mean(GaoStep1$WDIWDO)
    Gao_meanWgtDin <- mean(GaoStep1$WgtDin)
    Gao_meanWgtDout <- mean(GaoStep1$WgtDout)
    Gao_sdWgtDin <- sd(GaoStep1$WgtDin)
    Gao_sdWgtDout <- sd(GaoStep1$WgtDout)
    Gao_s <- abs(Gao_meanWDIWDO-(Gao_meanWgtDin*Gao_meanWgtDout))/(Gao_sdWgtDin*Gao_sdWgtDout)
    Gao_H <- 0.01*(Gao_sdWgtDout*Gao_sdWgtDin)/Gao_density
    Gindx <- Gao_density+Gao_H*Gao_s
    
    Gaoindx[thisYear,1] <- iYear
    Gaoindx[thisYear,2] <- Gao_density
    Gaoindx[thisYear,3] <- Gao_s
    Gaoindx[thisYear,4] <- Gao_H
    Gaoindx[thisYear,5] <- Gindx
    
    #x <- Gao_s
    x <- Gao_density
    H <- x*x*c_coefft-b_coefft*x+a_coefft
    
    Gaoindx[thisYear,6] <- H
    # Do scoring of result with respect to resilience rating  
    #Gstep1 <- ifelse(Gao_s > s_ref, 1, 0)
    Gstep1 <- ifelse(Gao_density > s_ref, 1, 0)
    Gstep2 <- ifelse(Gao_H > H_ref, 1, 0)
    Gstep3 <- Gstep1 + Gstep2
    if(Gstep3 > 1) {
      GaoScore <- "A"
    }
    if(Gstep3 == 1) {
      GaoScore <- "M"
    }
    if(Gstep3 < 1) {
      GaoScore <- "M"
      if(Gao_H < H){
        GaoScore <- "F"
      }
      
    }
    #GaoScore <- ifelse (Gao_s > s_ref, "A",(ifelse(Gao_H < H, "F", (ifelse(Gao_H < H_ref, "M", "A")))))
    
    Gaoindx[thisYear,7] <- GaoScore
  }
  
  # Get refrence line for the plot
  x <- seq(0.01, 6.95, by=0.01) 
  n <- length(x)
  H <- x*x*c_coefft-b_coefft*x+a_coefft
  gao_ref_pts <- data.frame(matrix(0, ncol = 2, nrow = n))
  colnames(gao_ref_pts)[1:2] <- c("s","H")
  gao_ref_pts$s <- x
  gao_ref_pts$H <- H
  
  write.csv(Gaoindx,  file = paste0(OutDir, "/Gaoindx.csv"), row.names = FALSE)
  
  ymin <- 1e-1
  xmin <- 1e-2
  ymax <- 1e+3
  xmax <- 1e+3
  
  Dmax <- max(Gaoindx$Gao_density)
  Hmax <- max(Gaoindx$Gao_H)
  
  # Check if need to zoom in
  if (Dmax < 1e2) {
    if (Hmax < 1e2) {
      ymax <- 1e+2
      xmax <- 1e+2
    }
  }
  
  # Plot the results as a scatter and line
  outPlotName <- paste(OutDir,"/Gao_Index.png",sep="")
  dfGaoPlot <- filter(Gaoindx, Year > yr_ignore)
  pG1 <- ggplot(data = dfGaoPlot, aes(x = Gao_density, y = Gao_H)) +
    geom_point(data = dfGaoPlot, size=1, aes(colour=rownames(dfGaoPlot))) +
    #geom_jitter(data = Gaoindx, size=1, aes(colour=rownames(Gaoindx)), width = 0.05, height = 0.25) + # Colour based on Year?
    geom_path(size = 0.2) + #so connected up through time (colour based on Year)
    geom_line(data = gao_ref_pts, aes(x = s, y = H)) + # refrence line details
    scale_x_log10(limits=c(xmin, xmax)) + scale_y_log10(limits=c(ymin, ymax)) +
    theme_bw() +
    labs(title = "Gao index", x = "<s>", y = "H", color='Year') +
    theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=axistitle,face="bold"),title=element_text(size=20,face="bold"))
  
  print(pG1)
  
  ggsave(pG1, file=outPlotName)
  
  # Plot the score through time
  GaoindxPlot <- Gaoindx
  GaoindxPlot <- subset (GaoindxPlot, select = -c(Gao_density, Gao_s, Gao_H, Gao_indx, resil_H))
  GaoindxPlot$n <- 1
  outPlotName <- paste(OutDir,"/Gao_Rating.png",sep="")
  pG2 <- ggplot(data = GaoindxPlot, aes(x = Year, y = n, fill = GaoScore)) + geom_bar(colour = "black", stat="identity", size = 0.1) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme_bw() +
    labs(x="Years", y = "Gao Rating") +
    theme(axis.text=element_text(size=textsize,face="bold"), axis.title=element_text(size=axistitle,face="bold"))
  
  print(pG2)
  
  ggsave(pG2, file=outPlotName)
  
  return(Gaoindx)
  
}


# Calcualtion ETI (a measure of ecosystem health based on network theory adn ecology)
##' @title Calculate ETI
##' @author Beth Fulton
do_ETI <- function(dfFSSI, dfScored, GaoIndx, id2){
  
  #Possible results - A = Acceptable, F = Fail, L = Light
  Score_Types <- c("F-F", "F-A", "A-F", "F-L", "L-F", "A-A", "A-L", "L-A", "L-L")
  nScoreT <- length(Score_Types)
  dfScore <- data.frame(matrix(0, ncol = 4, nrow = nScoreT))
  colnames(dfScore)[1:4] <- c("CombinedScore","BaseScore","FiveCase","SixCase","SevenCase")
  dfScore$CombinedScore <- Score_Types
  dfScore$BaseScore <-c(0, 0.25, 0.25, 0.5, 0.5, 1.0, 2.0, 2.0, 2.5)
  dfScore$FiveCase <-c(0, 0.675, 1.0125, 1.35, 2.025, 2.7, 4.05, 5.4, 6.75)
  dfScore$SixCase <- c(0, 0.925, 1.3875, 1.85, 2.775, 3.7, 5.55, 7.4, 9.25)
  #dfScore$SevenCase <- c(0, 1.175, 1.7625, 2.35, 3.525, 4.7, 7.05, 9.4, 11.75)
  dfScore$SevenCase <- c(2.5, 3.7, 4.6, 5.5, 6.4, 7.3, 8.2, 9.1, 10.0)
  
  # Class weighting
  ClassTypes <- c("Vulnerable","Habitat","Byproduct","Bycatch","Target","Robust","Hub")
  nClassT <- length(ClassTypes)
  dfClassScore <- data.frame(matrix(0, ncol = 2, nrow = nClassT))
  colnames(dfClassScore)[1:2] <- c("Classification","ClassScoreVal")
  dfClassScore$Classification <- ClassTypes
  dfClassScore$ClassScoreVal <- c(1, 1, 0.75, 0.5, 1, 0.5, 1)
  
  # Construct dataframes
  dfFSSIreuse <- dfFSSI
  dfFSSIreuse <- subset (dfFSSIreuse, select = c(Year, species_id, Group, IsDetritus, Classification, ScoreStatus))
  dfFSSIreuse <- dfFSSIreuse %>% dplyr::filter(IsDetritus < 1)
  
  dfScoredreuse <- dfScored
  dfScoredreuse <- subset (dfScoredreuse, select = c(Year, species_id, Group, Classification, DistortScore))
  
  # Scores
  dftmp <- merge(dfScoredreuse, dfFSSIreuse, by=c("Year","species_id","Group","Classification"))  # Intentional reuse
  nThisClass <- length(unique(dftmp$Classification))
  dfComp <- merge(dftmp, dfClassScore, by=c("Classification"))
  dfComp$CombinedScore <- paste(dfComp$ScoreStatus,"-",dfComp$DistortScore, sep="")
  
  # Classes
  ClassList <- unique(dfComp$Classification)
  nPerClass <- id2 %>% dplyr::count(Classification)
  
  # Find the response bin
  CHeckValID <- ifelse(nClassT == 7, 1, ifelse(nClassT == 6, 2, 3))
  CHeckVal <- dfScore$FiveCase  # Set as default
  if (CHeckValID ==1) {CHeckVal <- dfScore$SevenCase}
  if (CHeckValID ==2) {CHeckVal <- dfScore$SixCase}
  
  # Final score
  minYear <- min(df$Year)
  maxYear <- max(df$Year)
  nYear <- maxYear - minYear + 1
  dfCIndx <- data.frame(matrix(0, ncol = 3, nrow = nYear))
  colnames(dfCIndx)[1:3] <- c("Year","BaseETIIndx","ETIIndx")
  
  for (iYear in minYear:maxYear) {
    thisYear <- iYear - minYear + 1
    extractdf <- dfComp %>% dplyr::filter(Year == iYear)
    dftmp <- merge(extractdf, dfScore, by=c("CombinedScore"))  # Intentional reuse
    dftmp$FinalScore <- paste(dftmp$Classification,"-",dftmp$CombinedScore, sep="")
    
    #Oname <- paste(OutDir,"/dftmp_",iYear,".csv",sep="")
    #write.csv(dftmp, file = Oname, row.names = FALSE)
    
    nClassScores <- dftmp %>% 
      group_by(CombinedScore,Classification) %>%
      dplyr::summarise(ScoreStep1 = sum(BaseScore))
    dftmp1 <- merge(nClassScores, dfClassScore, by=c("Classification")) # Intentional reuse
    dfCalc <- merge(dftmp1, nPerClass,  by=c("Classification"))
    dfCalc$WgtVal <- dfCalc$ScoreStep1 * dfCalc$ClassScoreVal / dfCalc$n
    
    #Oname <- paste(OutDir,"/dfcalc_",iYear,".csv",sep="")
    #write.csv(dftmp, file =Oname, row.names = FALSE)
    
    dfGao <- Gaoindx %>% dplyr::filter(Year == iYear)
    GaoState <- dfGao$GaoScore
    #resilWgt <- ifelse(GaoState == "A", 1.0, (ifelse(GaoState == "M", 0.8, 0.5)))
    resilWgt <- 0.8
    ansScore <- sum(dfCalc$WgtVal) * resilWgt
    ansScore <- round(ansScore, 1)
    
    dfCIndx[thisYear,1] <- iYear
    dfCIndx[thisYear,2] <- ansScore
    
    if ( ansScore <= CHeckVal[1]) {
      ansClass <- 1
    } else if ( ansScore <= CHeckVal[2]) {
      ansClass <- 2
    } else if ( ansScore <= CHeckVal[3]) {
      ansClass <- 3
    } else if ( ansScore <= CHeckVal[4]) {
      ansClass <- 4
    } else if ( ansScore <= CHeckVal[5]) {
      ansClass <- 5
    } else if ( ansScore <= CHeckVal[6]) {
      ansClass <- 6
    } else if ( ansScore <= CHeckVal[7]) {
      ansClass <- 7
    } else if ( ansScore <= CHeckVal[8]) {
      ansClass <- 8
    } else if ( ansScore <= CHeckVal[9]) {
      ansClass <- 9
    } else {
      ansClass <- 10
    }
    
    dfCIndx[thisYear,3] <- ansClass
  }  
  
  Ofname <- paste(OutDir,"/ETIout.csv",sep="")
  write.csv(dfCIndx, file = Ofname, row.names = FALSE)
  
  
  # Create a column of colours
  
  # Score       Rating
  #   <2.5      Collapse
  #   2.5-3.7   Shocks Likely
  #   3.7-4.6   Shocks Possible
  #   4.6-5.5   Low Integrity
  #   5.5-6.4   Med Integrity
  #   6.4-7.3   Med-High Integrity
  #   7.3-8.2   High Integrity"
  #   8.2-9.1   Very Robust Integrity
  #   9.1-10.0  Close to Pristine
  #   >10.0     Pristine
  
  dfCIndx$ETIIndx2 <- as.integer(dfCIndx$ETIIndx)
  dfCIndx <- dfCIndx %>% mutate(ETIIndx2 = round(ETIIndx, 0)) %>% merge(colors, by.x='ETIIndx2', by.y = 'val')
  myColors <- brewer.pal(11, "RdYlGn")
  
  # Plot results
  plot_indx <- paste(OutDir,"/ETI_continuous.png",sep="")
  pE1 <- ggplot(dfCIndx, aes (x = Year, y = BaseETIIndx, colour=factor(ETIIndx))) +
    geom_point(size = 2) + 
    expand_limits(y=0) + expand_limits(y=10) +
    scale_colour_manual(breaks=colors$val, values=colors$col, labels=colors$labels) +
    theme_dark() +
    guides(color = guide_legend(reverse=TRUE)) +
    labs(
      #title = "Basic Ecosystem Traits & Health Index (BI)",
      title = "Ecosystem Traits Index (ETI)",
      x = "Year",
      y = "ETI",
      colour = "ETI rating"
    ) + 
    theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=10), axis.title=element_text(size=12,face="bold"), strip.text = element_text(face="bold", size=10)) 
  
  print(pE1)
  
  ggsave(pE1, file=plot_indx)
  
  # Plot results
  plot_indx <- paste(OutDir,"/ETI_block_band.png",sep="")
  pE2 <- ggplot(dfCIndx, aes (x = Year, y = ETIIndx, colour=factor(ETIIndx))) +
    geom_point(size = 2) + 
    expand_limits(y=0) + expand_limits(y=10) +
    scale_colour_manual(breaks=colors$val, values=colors$col, labels=colors$labels) +
    theme_dark() +
    guides(color = guide_legend(reverse=TRUE)) +
    labs(
      title = "Ecosystem Traits Index (ETI)",
      x = "Year",
      y = "ETI",
      colour = "ETI rating"
    ) + 
    theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=10), axis.title=element_text(size=12,face="bold"), strip.text = element_text(face="bold", size=10))

  print(pE2)
  
  ggsave(pE2, file=plot_indx)
  
}

