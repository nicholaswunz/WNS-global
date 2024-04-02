# Load library
pacman::p_load(Matrix, tidyverse, colorspace, cowplot, nlme, writexl,
               brms, rstan, ape, phytools, ggtree, ggnewscale, rasterSp, raster)

#devtools::install_github("alrobles/geotax")
#library(geotax)

# Functions 
mytheme <- function() {
  theme_bw() +
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(0.2, 0.2, 0.2, 0.2), units = , "cm"),
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
    )
} # set up plot theme

stack <- function(species_names = NA, path = getwd(), filename = NA, ...){
  if(!is.na(filename)){
    if(file.exists(filename)){
      r_species <- raster::raster(filename) # raster the file
    }
  } else if(!exists("r_species")){
    # Check which files are already there
    available_files <- list.files(path)
    
    # which filenames and species names overlap
    if(!anyNA(species_names)){
      available_names <- sapply(available_files, FUN=function(x){
        paste(strsplit(as.character(x), split="_")[[1]][1],strsplit(as.character(x), split="_")[[1]][2])})
      available_files <- available_files[which(available_names %in% species_names)]
      rm(available_names)
    }
    
    r_species <- raster::raster(paste0(path, available_files[1]))
    for(i in 2:length(available_files)){
      r_species <- raster::stack(r_species, raster::raster(paste0(path, available_files[i])))
    }
    # Set all 0s to NA
    r_species[raster::getValues(r_species) == 0] <- NA
    
    # Save r_species to file
    if(!is.na(filename)){
      raster::writeRaster(r_species, filename=filename, ...)
    }
  }
  return(r_species)
}

# Set directory and file locations
setwd('/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/') # Macbook
data_path    <- file.path("Manuscript - WNS comparative/data")
phylo_path   <- file.path("Methodology/Phylogeny")
spatial_path <- file.path("Methodology/Spatial map")

## RAW DATA ## --------------------------
bat_dat <- read.csv(file.path(data_path, "bat_comp_data.csv")) %>%
  dplyr::select(suborder:log10_Pd_load) %>%
  dplyr::mutate(log10_Pd_load = ifelse(Pd_load_ng_mm2_forearm_length > 0, log10(Pd_load_ng_mm2_forearm_length), NA),
                detected      = factor(detected),
                detected_bin  = ifelse(detected == "yes", 1, 0),
                Pd            = factor(Pd),
                WNS           = factor(WNS),
                WNS_decline   = factor(WNS_declining),
                suborder      = factor(suborder),
                family        = factor(family),
                roost_pref    = factor(roost_pref, levels = c("trees", "broad", "caves")),
                roost_n       = as.numeric(roost_pref),
                risk          = factor(risk, levels = c("Least Concern", "Near Threatened", "Vulnerable", "Endangered",
                                                           "Critically Endangered", "Data Deficient")),
                hemisphere    = factor(hemisphere),
                continent     = factor(continent),
                continent_2   = dplyr::recode(continent,
                                              "North America" = "North America/South America",
                                              "South America" = "North America/South America",
                                              "Asia"          = "Asia/Europe",
                                              "Europe"        = "Asia/Europe",
                                              "Asia"          = "Asia/Africa",
                                              "Africa"        = "Asia/Africa",
                                              "Australia"     = "Australia/Oceania",
                                              "Oceania"       = "Australia/Oceania")) %>%
  dplyr::filter(risk != "" & hemisphere != "" &
                  species_iucn != "Murina tenebrosa" & # only known by holotype in Japan
                  species_iucn != "Mystacina robusta") # no confirmed sightings of this species since 1965


## Merge with Tanalgo data
#bat_dat <- read.csv(file.path(data_path, "bat_comp_data.csv"))

#tanalgo_dat <- read.csv(file.path(data_path, "Tanalgo_appendix_1.csv")) %>% # from Tanalgo et al 2022
  #dplyr::select(Suborder:Mass..g.) %>%
  #dplyr::mutate(species_iucn = Species.Name,
                #species_iucn = gsub('<ca>','', species_iucn),
                #species_iucn = str_trim(species_iucn, "right"))

#test <- bat_dat %>% 
  #dplyr::left_join(tanalgo_dat, by = "species_iucn")

#test <- test[!duplicated(test), ] # remove rows with duplicated species

# setdiff(tanalgo_dat$species_iucn, bat_dat$species_iucn) check differences in species list

#library("writexl")
#write_xlsx(test, "Manuscript - WNS comparative/bat_comp_data.xlsx")

## PHYLOGENETIC TREE ## ------------------------------------------------------------------------------------------------------------
# Load tree
phylo_tree <- ape::read.nexus(file.path(phylo_path, "bat_tree_2015_S17613.nex")) # time-calibrated tree
#plot(phylo_tree, cex = 0.3)
#axisPhylo()

# Pruning data and phylogeny
tree_tip_label    <- phylo_tree$tip.label # extract tree tip names
species_list      <- bat_dat$species_tree # extract species name from bat data
pruned_tree       <- ape::drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_list)) # prune phylo_tree to keep species from the raw data
pruned_tree       <- ape::drop.tip(pruned_tree, c("Mus_musculus", "Sorex_araneus", "Canis_lupus")) %>%
  phytools::force.ultrametric(pruned_tree, method = "extend") # remove outgroups and ultrametricize the tree
phylo_dat         <- subset(bat_dat, species_tree %in% pruned_tree$tip.label) 
phylo_dat$species_tree <- factor(phylo_dat$species_tree, levels = c(tree_tip_label))

ape::is.ultrametric(pruned_tree)

# check if lengths match for both data and tree
length(unique(tree_tip_label)) # 815
length(unique(pruned_tree$tip.label)) # 782
length(unique(phylo_dat$species_tree)) # 782

# Plot phylogeny with traits
rownames(phylo_dat) <- NULL
trait_dat           <- phylo_dat %>%
  tibble::column_to_rownames(var = 'species_tree') %>% 
  droplevels(phylo_dat$hemisphere)

fam_cls <- list(Pteropodidae     = as.vector(filter(phylo_dat, family == "Pteropodidae")$species_tree),
                Phyllostomidae   = as.vector(filter(phylo_dat, family == "Phyllostomidae")$species_tree),
                Vespertilionidae = as.vector(filter(phylo_dat, family == "Vespertilionidae")$species_tree),
                Hipposideridae   = as.vector(filter(phylo_dat, family == "Hipposideridae")$species_tree),
                Emballonuridae   = as.vector(filter(phylo_dat, family == "Emballonuridae")$species_tree),
                Megadermatidae   = as.vector(filter(phylo_dat, family == "Megadermatidae")$species_tree),
                Molossidae       = as.vector(filter(phylo_dat, family == "Molossidae")$species_tree),
                Natalidae        = as.vector(filter(phylo_dat, family == "Natalidae")$species_tree),
                Cistugidae       = as.vector(filter(phylo_dat, family == "Cistugidae")$species_tree),
                Craseonycteridae = as.vector(filter(phylo_dat, family == "Craseonycteridae")$species_tree),
                Furipteridae     = as.vector(filter(phylo_dat, family == "Furipteridae")$species_tree),
                Miniopteridae    = as.vector(filter(phylo_dat, family == "Miniopteridae")$species_tree),
                Mormoopidae      = as.vector(filter(phylo_dat, family == "Mormoopidae")$species_tree),
                Mystacinidae     = as.vector(filter(phylo_dat, family == "Mystacinidae")$species_tree),
                Myzopodidae      = as.vector(filter(phylo_dat, family == "Myzopodidae")$species_tree),
                Noctilionidae    = as.vector(filter(phylo_dat, family == "Noctilionidae")$species_tree),
                Nycteridae       = as.vector(filter(phylo_dat, family == "Nycteridae")$species_tree),
                Rhinolophidae    = as.vector(filter(phylo_dat, family == "Rhinolophidae")$species_tree),
                Rhinopomatidae   = as.vector(filter(phylo_dat, family == "Rhinopomatidae")$species_tree),
                Thyropteridae    = as.vector(filter(phylo_dat, family == "Thyropteridae")$species_tree))

ord_cls <- list(Yinpterochiroptera = as.vector(filter(phylo_dat, suborder == "Yinpterochiroptera")$species_tree),
                Yangochiroptera    = as.vector(filter(phylo_dat, suborder == "Yangochiroptera")$species_tree))

pruned_tree <- ggtree::groupOTU(pruned_tree, fam_cls, group_name = "family")
pruned_tree <- ggtree::groupOTU(pruned_tree, ord_cls, group_name = "suborder")

# Outline "Pteropodidae", "Hipposideridae", "Rhinolophidae", "Phyllostomidae", "Molossidae", "Vespertilionidae"

# Fig 1.
tree_plot <- ggtree::ggtree(pruned_tree, aes(color = family, linetype = suborder), layout = "circular", size = 0.5) +
  scale_linetype_manual(values = c("solid", "solid", "dotted")) +
  theme(legend.position = "bottom")

hemi_phylo <- ggtree::gheatmap(tree_plot + ggnewscale::new_scale_fill(), trait_dat %>% dplyr::select(hemisphere), color = NA, width = 0.08, colnames_angle = 90) +
  scale_fill_manual(values = c("#C7E5BE", "#64A79A", "#005D67"), name = "Hemisphere") + # North, North/South, South
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

cave_phylo <- ggtree::gheatmap(hemi_phylo + ggnewscale::new_scale_fill(), trait_dat %>% dplyr::select(roost_pref), color = NA, offset = 5, width = 0.08, colnames_angle = 90) +
  scale_fill_manual(values = c("white", "#200128"), name = "Cave use") + # broad and caves
  theme(legend.position = "bottom")

hib_phylo <- ggtree::gheatmap(cave_phylo + ggnewscale::new_scale_fill(), trait_dat %>% dplyr::select(hib_study), color = NA, offset = 10, width = 0.08, colnames_angle = 90) +
  scale_fill_manual(values = c("#1D396A"), name = "Hibernator", na.value = NA) + # yes and blank
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

ggtree::gheatmap(hib_phylo + ggnewscale::new_scale_fill(), trait_dat %>% dplyr::select(detected), offset = 15, width = 0.08, colnames_angle = 90, color = NA) +
  scale_fill_manual(values = c("#808080", "#DD3C26"), name = NULL, na.value = NA) +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

## PD BAT INCIDENCE MATRIX ## ----------------------------------
pd_dat  <- phylo_dat %>% dplyr::filter(detected != "") %>% dplyr::select(detected, species_tree)
pd_tree <- ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, pd_dat$species_tree))

wns_dat  <- phylo_dat %>% dplyr::filter(WNS != "") %>% dplyr::select(WNS, species_tree)
wns_tree <- ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, wns_dat$species_tree))

# check if lengths match for both data and tree
length(unique(pd_tree$tip.label)) # 65
length(unique(pd_dat$species_tree)) # 65
length(unique(wns_tree$tip.label)) # 52
length(unique(wns_dat$species_tree)) # 52

phylo_dat %>%
  dplyr::filter(detected != "") %>%
  group_by(genus) %>%
  summarise(n = length(genus)) %>%
  dplyr::arrange(-n)

phylo_dat %>%
  dplyr::filter(WNS_decline != "") %>%
  group_by(genus) %>%
  summarise(n = length(genus)) %>%
  dplyr::arrange(-n)

## For Pd-exposure risk
# Calculate distance matrix
pd_d <- ape::cophenetic.phylo(pd_tree) # computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
pd_d <- log10(pd_d + 1)

# Generate an incidence matrix from interaction data base
bat_i_matrix_pd <- geotax::get_incidence_matrix(pd_dat)[, colnames(pd_d)]

# Logistic regression coefficients for yes Pd and no Pd (only use yes)
pd_coef <- sapply(1:nrow(bat_i_matrix_pd), function(x) geotax::log_reg_boostrap(bat_i_matrix_pd[x, , drop = F], pd_d, 1000) ) %>% t

PD      <- seq(0, 60, 1) #x axis values across range of phylogenetic distances
pd_ps1 <- apply(pd_coef[ ,c(1,7)], 1, function(x) geotax::prob_logit(x, log10(PD + 1))) %>%
  data.frame() %>%
  dplyr::rename(no = X1, yes = X2)

pd_plot <- pd_ps1 %>%
  dplyr::select(yes) %>%
  tibble::add_column(pd = PD) %>%
  ggplot(aes(x = pd, y = yes)) +
  ylim(0, 1) +
  geom_line(size = 1) +
  labs(y = "Probability of Pd detection",
       x = "Phylogenetic distance from source to target host (my)") +
  mytheme()

## For WNS
# Calculate distance matrix
wns_d <- ape::cophenetic.phylo(wns_tree) # computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
wns_d <- log10(wns_d + 1)

# Generate an incidence matrix from interaction data base
bat_i_matrix_wns <- geotax::get_incidence_matrix(wns_dat)[, colnames(wns_d)]

# Logistic regression coefficients for yes Pd and no Pd (only use yes)
wns_coef <- sapply(1:nrow(bat_i_matrix_wns), function(x) geotax::log_reg_boostrap(bat_i_matrix_wns[x, , drop = F], wns_d, 1000) ) %>% t

# probability phylogenetic distance
PD      <- seq(0, 60, 1) #x axis values across range of phylogenetic distances
wns_ps1 <- apply(wns_coef[ ,c(1,7)], 1, function(x) geotax::prob_logit(x, log10(PD + 1))) %>%
  data.frame() %>%
  dplyr::rename(no = X1, yes = X2)

wns_plot <- wns_ps1 %>%
  dplyr::select(yes) %>%
  tibble::add_column(pd = PD) %>%
  ggplot(aes(x = pd, y = yes)) +
  ylim(0, 1) +
  geom_line(size = 1) +
  labs(y = "Probability of developing WNS",
       x = "Phylogenetic distance from source to target host (my)") +
  mytheme()

cowplot::plot_grid(pd_plot, wns_plot, ncol = 2)

levels(phylo_dat$continent_2)

# Sensitivity analysis - Eurasia only
asia_dat  <- phylo_dat %>% dplyr::filter(WNS != "" & continent_2 == "Asia/Europe") %>% dplyr::select(WNS, species_tree)
asia_tree <- ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, asia_dat$species_tree))

# Calculate distance matrix
asia_d <- ape::cophenetic.phylo(asia_tree) # computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
asia_d <- log10(asia_d + 1)

# Generate an incidence matrix from interaction data base
bat_i_matrix_asia <- geotax::get_incidence_matrix(asia_dat)[, colnames(asia_d)]

# Logistic regression coefficients for yes Pd and no Pd (only use yes)
asia_coef <- sapply(1:nrow(bat_i_matrix_asia), function(x) geotax::log_reg_boostrap(bat_i_matrix_asia[x, , drop = F], asia_d, 1000) ) %>% t

# probability phylogenetic distance
PD      <- seq(0, 60, 1) #x axis values across range of phylogenetic distances
asia_ps1 <- apply(asia_coef[ ,c(1,7)], 1, function(x) geotax::prob_logit(x, log10(PD + 1))) %>%
  data.frame() %>%
  dplyr::rename(no = X1, yes = X2)

asia_plot <- asia_ps1 %>%
  dplyr::select(yes) %>%
  tibble::add_column(pd = PD) %>%
  ggplot(aes(x = pd, y = yes)) +
  ylim(0, 1) +
  geom_line(size = 1) +
  labs(y = "Probability of developing WNS - Eurasia",
       x = "Phylogenetic distance from source to target host (my)") +
  mytheme()

cowplot::plot_grid(wns_plot, asia_plot, ncol = 2)


# Plot tree
tree_plot  <- ggtree(pd_tree, layout = "rectangular")
cont_phylo <- ggtree::gheatmap(tree_plot, trait_dat %>% dplyr::select(continent_2), offset = 0, width = 0.05, colnames_angle = 90, color = NA) +
  scale_fill_viridis_d(option = "viridis", name = NULL) +
  theme(legend.position = "bottom")

detect_phylo <- ggtree::gheatmap(cont_phylo + new_scale_fill(), trait_dat %>% dplyr::select(detected), offset = 3, width = 0.05, colnames_angle = 90, color = NA) +
  scale_fill_manual(values = c("#808080", "#e87674"), name = NULL, na.value = "#808080") +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

wns_phylo <- ggtree::gheatmap(detect_phylo + new_scale_fill(), trait_dat %>% dplyr::select(WNS), offset = 6, width = 0.05, colnames_angle = 90, color = NA) +
  scale_fill_manual(values = c("#808080", "#de5957"), name = NULL, na.value = "#808080") +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

declin_phylo <- ggtree::gheatmap(wns_phylo + new_scale_fill(), trait_dat %>% dplyr::select(WNS_decline), offset = 9, width = 0.05, colnames_angle = 90, color = NA) +
  scale_fill_manual(values = c("#808080", "#D6302D"), name = NULL, na.value = "#808080") +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

fam_phylo <- ggtree::gheatmap(declin_phylo + new_scale_fill(), trait_dat %>% dplyr::select(family), offset = 13, width = 0.01, colnames_angle = 90, color = NA) +
  scale_fill_viridis_d(option = "magma", name = NULL) +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

d2 <- data.frame(id = wns_tree$tip.label, log10_Pd_load = phylo_dat$log10_Pd_load[match(wns_tree$tip.label, phylo_dat$species_tree)])
trait_dat$species_tree <- phylo_dat$species_tree[match(rownames(trait_dat), phylo_dat$species_tree)]

p2 <- facet_plot(fam_phylo, panel = "dot", data = d2, geom = geom_point, aes(x = log10_Pd_load, colour = log10_Pd_load)) +
  colorspace::scale_color_continuous_diverging(palette = "Blue-Red", mid = -6, name = NULL) +
  theme_tree2() +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom")

# load the packaged
#library(grid)
#library(gtable)

gt <- ggplot_gtable(ggplot_build(p2))
gtable_show_layout(gt) # will show you the layout - very handy function
#gt # see plot layout in table format
gt$widths[7] = 0.3 * gt$widths[7] # in this case it was column 7 - reduce the width by a half
grid::grid.draw(gt) # plot with grid draw

# Phylogenetic signal - Pd load
pd_load <- setNames(trait_dat$log10_Pd_load, rownames(trait_dat))
phytools::phylosig(pd_tree, pd_load, method = "lambda", test = TRUE, nsim = 1000)

test <- trait_dat %>% dplyr::filter(continent_2 == "Asia/Europe")
test_load <- setNames(test$log10_Pd_load, rownames(test))
phytools::phylosig(asia_tree, test_load, method = "lambda", test = TRUE, nsim = 1000)

## BAT SPECIES RICHNESS ## -----------------------------------------------
filedir <- "/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/Manuscript - WNS comparative/"
#filedir <- "C:/Users/30061860/OneDrive - Western Sydney University/WNS bat project - WSU/Manuscript - Ecological naivety/"

# Convert shape files into rasters and save to file (downloaded shape files 08/03/2022)
# only extracted distribution of native range
#mammals_rast <- rasterSp::rasterizeIUCN(dsn = paste0(filedir, "/IUCN/MAMMALS_TERRESTRIAL_ONLY.shp"), resolution = 0.25, seasonal = c(1,2), origin = 1, presence = c(1,2), save = TRUE, path = paste0(filedir, "/SpeciesData_25/"))

# For AusBat data
#Moo <- raster::shapefile("/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/Manuscript - WNS comparative/AusBatMap/Miniopterus_orianae_oceanensis/Miniopterus_orianae_oceanensis.shp")
#Mob <- raster::shapefile("/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/Manuscript - WNS comparative/AusBatMap/Miniopterus_orianae_bassanii/Miniopterus_orianae_bassanii.shp")
#Moor <- raster::shapefile("/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/Manuscript - WNS comparative/AusBatMap/Miniopterus_orianae_orianae/Miniopterus_orianae_orianae.shp")
#Rr  <- raster::shapefile("C:/Users/30061860/OneDrive - Western Sydney University/WNS bat project - WSU/Manuscript - Ecological naivety/AusBatMap/Rhinolophus_robertsi/Rhinolophus_robertsi.shp")
#Mo <- rbind(Moo, Mob, Moor)
#Mo <- sp::spTransform(Mo, raster::crs("+proj=longlat +datum=WGS84 +no_defs"))

#pgeo <- sp::spTransform(Mo, raster::crs('+proj=longlat +datum=WGS84'))
#ext  <- raster::extent(bat_cave_sr)
#rr   <- raster::rasterize(pgeo, raster::raster(ext, res = 0.25), field = 1)
#crs(rr) <- "+proj=longlat +datum=WGS84 +no_defs"

#raster::writeRaster(rr, filename = "Miniopterus_orianae_0.25.tif", format = "GTiff", overwrite = TRUE)

bat_cave_sr <- rasterSp::calcSR(species_names = as.matrix(bat_dat %>% dplyr::filter(roost_pref == "caves") %>%
                                                       dplyr::select(species_iucn)), path = paste0(filedir, "/SpeciesData_25/"))
#saveRDS(bat_cave_sr, "Manuscript - WNS comparative/bat_cave_sr.rds")
bat_cave_sr <- readRDS("Manuscript - WNS comparative/bat_cave_sr.rds")

# Convert raster to matrix then to data frame
bat_cave_df <- raster::as.data.frame(raster::rasterToPoints(bat_cave_sr)) %>% rename(cave_sp_n = layer)

# Figure 1
world <- map_data("world")
richness_plot <- ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(long, lat, map_id = region),
                              fill = "dark grey", colour = NA) +
  geom_raster(data = bat_cave_df, aes(y = y, x = x, fill = log(cave_sp_n))) +
  mytheme() + ylab(NULL) + xlab(NULL) +
  scale_fill_viridis_c(option = "inferno", name = "Species richness", end = 0.96,
                       breaks = fivenum(log(bat_cave_df$cave_sp_n)),
                       labels = fivenum(bat_cave_df$cave_sp_n)) +
  labs(title = "Bat species richness known to roost in caves", x = NULL, y = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5)) +
  coord_fixed()

## RISK MAP ## -------------------------------------------
# Extract BIOCLIM data
bio_dat <- getData(name = "worldclim", var = "bio", res = 5)

bio_dat_resampled <- projectRaster(bio_dat, bat_cave_sr, method = 'ngb')
tmean <- raster::as.data.frame(bio_dat_resampled$bio1 / 10, xy = TRUE) %>% # divide by 10 to get Degrees Celsius
  mutate(dist_m = ifelse(bio1 == 'NA', NA, 50), # 50 m cave distance
         dist_100 = ifelse(bio1 == 'NA', NA, 100),) # 100 m cave distance

# Load cave temperature data and scale MAST and distance for model
cave_raw <- read.csv("/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/Methodology/Climate loggers/cave_temp.csv") %>%
  dplyr::select(!c("notes", "title", "ref", "X")) %>%
  dplyr::mutate(MAST.s = as.numeric(scale(MAST)),
                dist.s = as.numeric(scale(dist_m))) %>%
  dplyr::filter(bats != "yes" & season == "Winter" & MAST != "NA") # remove cave temps next to roosting

#Fit model with with MAST and distance for caves only (with exponential covariance structure)
temp_mod <- nlme::lme(cave_temp ~ MAST.s * dist.s, 
                      random = ~ 1 | site, 
                      data = cave_raw)

temp_mod2 <- nlme::lme(cave_temp ~ MAST.s, # without dist.s
                      random = ~ 1 | site, 
                      data = cave_raw)
summary(temp_mod)
r2(temp_mod, n = NULL)
r2(temp_mod2, n = NULL) # without cave depth
sjPlot::plot_model(temp_mod, type = "pred", terms = c("MAST.s", "dist.s"))

cave_raw$cave_temp_pred <- predict(temp_mod, newdata = cave_raw, type = "response")

valid_plot <- ggplot(cave_raw) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(aes(x = cave_temp, y = cave_temp_pred, colour = log(dist_m + 1))) +
  xlab("Observed cave temperature (°C)") + ylab("Predicted cave temperature (°C)") +
  scale_color_viridis_c(option = "mako", name = "Cave depth (m)", end = 0.96,
                        breaks = fivenum(log(cave_raw$dist_m + 1)),
                        labels = round(fivenum(cave_raw$dist_m))) +
  mytheme() + theme(legend.position = "bottom")

pred_plot <- ggplot(cave_raw) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(aes(x = MAST2, y = cave_temp_pred, colour = log(dist_m + 1))) +
  xlab("Surface temperature (°C)") + ylab("Predicted cave temperature (°C)") +
  scale_color_viridis_c(option = "mako", name = "Cave depth (m)", end = 0.96,
                        breaks = fivenum(log(cave_raw$dist_m + 1)),
                        labels = round(fivenum(cave_raw$dist_m))) +
  mytheme() + theme(legend.position = "bottom")

# Fig S1
cowplot::plot_grid(valid_plot, pred_plot, ncol = 2)

# Scale tmean to the mean/sd of the original dataset and predict internal cave temperatures using temp_mod
tcave <- tmean %>%
  filter(!is.na(bio1)) %>%
  mutate(MAST.s = as.numeric(scale(bio1, center = mean(cave_raw$MAST), scale = sd(cave_raw$MAST))),
         dist.s = as.numeric(scale(dist_m, center = mean(cave_raw$dist_m), scale = sd(cave_raw$dist_m))), 
         dist100.s = as.numeric(scale(dist_100, center = mean(cave_raw$dist_m), scale = sd(cave_raw$dist_m))),
         cave_temp = 11.093605 + (MAST.s * 6.884297) + (dist.s * 3.362448) - (MAST.s * dist.s * 0.918320),
         cave_temp100 = 11.093605 + (MAST.s * 6.884297) + (dist100.s * 3.362448) - (MAST.s * dist100.s * 0.918320)
  )

tcave %>%
  ggplot() + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(aes(x = bio1, y = cave_temp)) +
  geom_point(aes(x = bio1, y = cave_temp100, colour = "red")) +
  theme_classic()

tmean_df <- tcave %>% filter(y > -56 & bio1 > 0 & bio1 < 19.8)
tcave_df <- tcave %>% filter(y > -56 & cave_temp100 > 0 & cave_temp < 19.8)

tmean_bat_df <- merge(tmean_df, bat_cave_df, c("x" ,"y"))
tcave_bat_df <- merge(tcave_df, bat_cave_df, c("x" ,"y"))

# Fig 2
tmean_plot <- ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(map_id = region), fill = "light grey", colour = NA) +
  geom_raster(data = tmean_bat_df, aes(y = y, x = x, fill = log(cave_sp_n))) +
  mytheme() + ylab(NULL) + xlab(NULL) +
  scale_fill_viridis_c(option = "inferno", name = "Species richness", end = 0.96,
                       breaks = fivenum(log(tmean_bat_df$cave_sp_n)),
                       labels = fivenum(tmean_bat_df$cave_sp_n)) +
  #labs(title = "Cave roosting species richness within Pd growth limit \nbased on MAST", x = "Longitude", y = 'Latitude') +
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_x_continuous(expand = c(0, 0)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5)) +
  coord_fixed()

tcave_plot <- ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(map_id = region), fill = "light grey", colour = NA) +
  geom_raster(data = tcave_bat_df, aes(y = y, x = x, fill = log(cave_sp_n))) +
  mytheme() + ylab(NULL) + xlab(NULL) +
  scale_fill_viridis_c(option = "inferno", name = "Species richness", end = 0.96,
                       breaks = fivenum(log(tcave_bat_df$cave_sp_n)),
                       labels = fivenum(tcave_bat_df$cave_sp_n)) +
  #labs(title = "Cave roosting species richness within Pd growth limit\nbased on cave-corrected MAST", x = "Longitude", y = 'Latitude') +
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_x_continuous(expand = c(0, 0)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5)) +
  coord_fixed()

cowplot::plot_grid(tmean_plot, tcave_plot, ncol = 1, labels = c('a', 'b'))

## SPECIES RISK ## ---------------------------------
# Select only cave-roosting species from continents not in Asia, Europe, and North America
cave_species <- bat_dat %>% dplyr::filter(roost_pref == "caves",
                                          continent %in% c("South America", "Oceania", "Africa", "North America/South America",
                                                           "Australia", "Asia/Australia"),
                                          hemisphere != "Northern",
                                          family != "Pteropodidae", # not known to hibernate
                                          family != "Mormoopidae")  # not known to hibernate

iucn_list <- list.files(path = paste0("/Users/nicholaswu/Library/CloudStorage/OneDrive-WesternSydneyUniversity/WNS bat project - WSU/Manuscript - WNS comparative/SpeciesData_25/"), 
                             pattern = ".tif", full.names = TRUE)

sp_shp <- raster::stack(iucn_list)
#saveRDS(sp_shp, "Manuscript - WNS comparative/sp_shp.rds")
sp_shp <- readRDS("Manuscript - WNS comparative/sp_shp.rds")
names(sp_shp) <- sub("_0.25.*", "", names(sp_shp)) # remove _0.25 from species name

sp_names <- as.data.frame(names(sp_shp))

cave_sp_names <- as.matrix(cave_species %>% 
                             dplyr::select(species_iucn) %>%
                             mutate(species_iucn = sub(" ", "_", species_iucn))) 

cave_sp_shp <- raster::subset(sp_shp, cave_sp_names, value = T)

# Merge all raster layers and convert to dataframe
cave_sp_df <- raster::as.data.frame(cave_sp_shp, xy = T)
colnames(cave_sp_df) <- sub("_", " ", colnames(cave_sp_df)) # remove _0.25 from species name

# Extract MAST and filter layer by species distribution
bio_dat           <- raster::getData(name = "worldclim", var = "bio", res = 5)
bio_dat_resampled <- raster::projectRaster(bio_dat, cave_sp_shp, method = 'ngb') # reproject MAST raster to fit bat raster

tmean_df     <- tcave %>% dplyr::filter(y > -56 & !is.na(bio1))
tmean_pd_df  <- tcave %>% dplyr::filter(y > -56 & bio1 > 0 & bio1 < 19.8 & !is.na(bio1))
tmean_cave_pd_df <- tcave %>% dplyr::filter(y > -56 & cave_temp100 > 0 & cave_temp < 19.8 & !is.na(bio1))

# Merge MAST with species map
temp_bat_all <- merge(tmean_df, cave_sp_df, c("x" ,"y"))
temp_bat_pd  <- merge(tmean_pd_df, cave_sp_df, c("x" ,"y"))
temp_bat_pd_cave <- merge(tmean_cave_pd_df, cave_sp_df, c("x" ,"y"))

# calculate number of spatial grid cells occupied by each species
all_sum <- temp_bat_all %>% summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  dplyr::select(-c(1:10)) %>%
  tidyr::pivot_longer(everything(), names_to = "species", values_to = "total_sum")

pd_sum <- temp_bat_pd %>% summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  dplyr::select(-c(1:10)) %>%
  tidyr::pivot_longer(everything(), names_to = "species", values_to = "pd_sum") %>%
  merge(all_sum) %>%
  data.frame() %>%
  dplyr::mutate(percent = round(pd_sum / total_sum * 100, 2))

# calculate number of spatial grid cells (within Pd limit) occupied by each species, merge, and calculate difference between total and Pd limit
pd_sum_cave <- temp_bat_pd_cave %>% summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  dplyr::select(-c(1:10)) %>%
  tidyr::pivot_longer(everything(), names_to = "species", values_to = "pd_sum2") %>%
  merge(all_sum) %>%
  data.frame() %>%
  dplyr::mutate(species      = sub("_", " ", species),
                family       = cave_species$family[match(species, cave_species$species_iucn)],
                percent_cave = round(pd_sum2 / total_sum * 100, 2),
                continent    = cave_species$continent[match(species, cave_species$species_iucn)],
                hibernator   = cave_species$hib_study[match(species, cave_species$species_iucn)],
                roost_pref   = cave_species$roost_pref[match(species, cave_species$species_iucn)],
                risk         = cave_species$risk[match(species, cave_species$species_iucn)]) %>%
  dplyr::arrange(continent, -percent_cave) %>%
  dplyr::filter(percent_cave > 5)

final_table <- merge(pd_sum_cave, pd_sum, by = "species") %>%
  dplyr::arrange(continent, -percent)

write_xlsx(final_table,"overlap_table.xlsx")

## KARST MAP ## --------------------------------------
karst_poly  <- rgdal::readOGR("Methodology/Spatial map/WHYMAP_WOKAM/shp", layer = "whymap_karst__v1_poly") # import WOKAM shape file layers
bio_karst <- mask(bio_dat, karst_poly)

bio_karst_resampled <- projectRaster(bio_karst, bat_cave_sr, method = 'ngb')
tmean_karst <- raster::as.data.frame(bio_karst_resampled$bio1 / 10, xy = TRUE) %>% # divide by 10 to get Degrees Celsius
  mutate(dist_m = ifelse(bio1 == 'NA', NA, 50), # 50 m cave distance
         dist_100 = ifelse(bio1 == 'NA', NA, 100),) # 100 m cave distance

tmean_karst <- tmean_karst %>%
  filter(!is.na(bio1)) %>%
  mutate(MAST.s = as.numeric(scale(bio1, center = mean(cave_raw$MAST), scale = sd(cave_raw$MAST))),
         dist.s = as.numeric(scale(dist_m, center = mean(cave_raw$dist_m), scale = sd(cave_raw$dist_m))), 
         dist100.s = as.numeric(scale(dist_100, center = mean(cave_raw$dist_m), scale = sd(cave_raw$dist_m))),
         cave_temp = 11.093605 + (MAST.s * 6.884297) + (dist.s * 3.362448) - (MAST.s * dist.s * 0.918320),
         cave_temp100 = 11.093605 + (MAST.s * 6.884297) + (dist100.s * 3.362448) - (MAST.s * dist100.s * 0.918320)
  )

#tmean_karst_df <- tmean_karst %>% dplyr::filter(y > -56 & bio1 > 0 & bio1 < 19.8 & !is.na(bio1)) # filter out mean temp above 19.8 C and below 0 C
tcave_karst_df <- tmean_karst %>% dplyr::filter(y > -56 & cave_temp100 > 0 & cave_temp < 19.8 & !is.na(bio1)) # filter out mean temp above 19.8 C and below 0 C

#tmean_bat_df <- merge(tmean_karst_df, bat_cave_df, c("x" ,"y"))
tcave_bat_df <- merge(tcave_karst_df, bat_cave_df, c("x" ,"y"))

# Fig S2
ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(long, lat, map_id = region), fill = "light grey", colour = NA) +
  geom_polygon(data = karst_poly, aes(x = long, y = lat, group = id), fill = "black") +
  geom_raster(data = tcave_bat_df, aes(y = y, x = x), fill = "#CC1C2F") +
  mytheme() + ylab(NULL) + xlab(NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom") +
  coord_fixed()

## HIBERNATION STUDIES ## --------------------------------------
hib_dat <- read.csv(file.path(data_path, "hibernation_study.csv"))

# Fig S4
hib_dat %>%
  group_by(year, hemisphere, temp_sen) %>%
  summarise(n = length(unique(study_ID))) %>%
  ggplot(aes(x = year, y = n, fill = temp_sen)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("grey", "red")) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  facet_grid(hemisphere ~ .) +
  mytheme()