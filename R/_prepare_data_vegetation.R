# Beta diversity on dike grasslands
# Prepare species, sites, and traits data ####

# Markus Bauer
# 2024-02-28



# Content #####################################################################
# A Load data +++++++++++++++++++++++++++++++++++++
## 1 Sites
## 2 Species
## 3 Traits
## 4 Temperature and precipitation
## 5 Check data frames
# B Create variables ++++++++++++++++++++++++++++++
## 1 Create simple variables
## 2 Coverages
## 3 Alpha diversity
## 4 Biotope types
## 5 LCBD: Local contributions to beta diversity
## 6 Synchrony
## 7 PCAs: Principle component analyses
## 8 ESy: Calculate EUNIS habitat types
## 9 dbMEM: Distance-based Moran's eigenvector maps
## 10 TBI: Temporal beta diversity index
## 11 Finalization of data frames
# C Save processed data +++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(ellipsis)
suppressPackageStartupMessages(library(installr))
suppressPackageStartupMessages(library(lubridate))
library(tidyverse)
library(naniar)
suppressPackageStartupMessages(library(TNRS))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(FD))
suppressPackageStartupMessages(library(adespatial))
library(remotes)
#remotes::install_github(file.path("larsito", "tempo"))
library(tempo)


### Start ###
rm(list = ls())
#installr::updateR(browse_news = FALSE, install_R = TRUE, copy_packages = TRUE, copy_Rprofile.site = TRUE, keep_old_packages = TRUE, update_packages = TRUE, start_new_R = FALSE, quit_R = TRUE, print_R_versions = TRUE, GUI = TRUE)
#renv::status()



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Sites of dikes ###########################################################


sites_dikes <- read_csv(
  here("data", "raw", "data_raw_sites.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"), col_types = cols(
      .default = "?",
      id = "f",
      block = "f",
      survey_date.2017 = col_date(format = "%Y-%m-%d"),
      survey_date.2018 = col_date(format = "%Y-%m-%d"),
      survey_date.2019 = col_date(format = "%Y-%m-%d"),
      survey_date.2021 = col_date(format = "%Y-%m-%d")
    )
  ) %>%
  mutate(across(
    starts_with("vegetation_cover") |
      starts_with("vegetation_height") |
      starts_with("opensoil_cover") |
      starts_with("moss_cover") |
      starts_with("litter_cover") |
      starts_with("botanist"),
    ~ as.character(.x)
  )) %>%
  pivot_longer(
    starts_with("vegetation_cover.") |
      starts_with("vegetation_height.") |
      starts_with("opensoil_cover.") |
      starts_with("moss_cover.") |
      starts_with("litter_cover.") |
      starts_with("botanist."),
    names_to = c("x", "survey_year"),
    names_sep = "\\.",
    values_to = "n"
    ) %>%
  pivot_wider(names_from = x, values_from = n) %>%
  mutate(across(
    c(survey_year, vegetation_cover, vegetation_height,
      opensoil_cover, moss_cover, litter_cover),
    ~ as.numeric(.x)
  )) %>%
  mutate(
    survey_year_factor = factor(survey_year),
    survey_year_factor_minus = factor(survey_year - 1),
    construction_year_factor = factor(construction_year),
    construction_year_factor_plus = factor(construction_year + 1)
  ) %>%
  mutate(
    id = str_c(id, survey_year, sep = "_"),
    .keep = "all",
    id = paste0("X", id),
    plot = str_sub(id, start = 2, end = 3),
    position = str_sub(id, start = 5, end = 5),
    location_abb = str_sub(location, 1, 3),
    location_abb = str_to_upper(location_abb),
    location_abb = factor(location_abb,
      levels = unique(location_abb[order(construction_year)])
    ),
    location_construction_year = str_c(location_abb, construction_year,
                                       sep = "-")
  ) %>%
  select(
    -position, -starts_with("survey_date_"), -starts_with("topsoil_depth_"),
    -cn_level, -ends_with("_class")
  ) %>%
  relocate(
    plot, .after = id
  ) %>%
  relocate(
    c("location_abb", "location_construction_year"), .after = location
  ) %>%
  relocate(
    c("survey_year", "survey_year_factor", "survey_year_factor_minus"),
    .after = river_km
  ) %>%
  relocate(
    c("construction_year_factor", "construction_year_factor_plus"),
    .after = construction_year
  )

### Sites with spatial calculations ###
sites_with_spatial_data <- read_csv(
  here("data", "processed", "spatial", "sites_with_spatial_data.csv"),
  col_names = TRUE, na = c("", "NA", "na"),
  col_types = cols(.default = "?")
  ) %>%
  select(plot, longitude_center, latitude_center,
         river_distance, biotope_area, biotope_distance)
sites_dikes <- sites_dikes %>%
  left_join(sites_with_spatial_data, by = "plot")



## 2 Species ##################################################################


species_dikes <- data.table::fread(
  here("data", "raw", "data_raw_species.csv"),
  sep = ",", dec = ".", skip = 0, header = TRUE, na.strings = c("", "NA", "na"),
  colClasses = list(
    character = "name"
    )
  ) %>%
  ### Check that each species occurs at least one time ###
  group_by(name) %>%
  arrange(name) %>%
  select(name, all_of(sites_dikes$id)) %>%
  mutate(
    total = sum(c_across(starts_with("X")), na.rm = TRUE),
    presence = if_else(total > 0, 1, 0),
    name = factor(name)
  ) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))



### Create list with species names and their frequency ###
specieslist <- species_dikes %>%
  mutate_if(is.numeric, ~ 1 * (. != 0)) %>%
  mutate(
    sum = rowSums(across(where(is.numeric)), na.rm = TRUE), .keep = "unused"
  ) %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
#write_csv(specieslist, "specieslist_2022xxxx.csv")



## 3 Traits ###################################################################


traits <- read_csv(
  here("data", "raw", "data_raw_traits.csv"),
  col_names = TRUE, na = c("", "NA", "na"),
  col_types =
    cols(
      .default = "f",
      name = "c",
      sociology = "d",
      l = "d",
      t = "d",
      k = "d",
      f = "d",
      r = "d",
      n = "d"
    )) %>%
  separate(
    name, c("genus", "species", "ssp", "subspecies"), "_",
    remove = FALSE, extra = "drop", fill = "right"
  ) %>%
  mutate(
    genus = str_sub(genus, 1, 4),
    species = str_sub(species, 1, 4),
    subspecies = str_sub(subspecies, 1, 4),
    name = factor(name)
  ) %>%
  unite(abb, genus, species, subspecies, sep = "") %>%
  mutate(
    abb = str_replace(abb, "NA", ""),
    abb = as_factor(abb)
  ) %>%
  select(-ssp) %>%
  arrange(name)

### Check species names ###
anti_join(traits, species_dikes, by = "name") %>% select(name)
anti_join(species_dikes, traits, by = "name") %>% select(name)
data <- traits %>%
  mutate(name = str_replace_all(name, "_", " ")) %>%
  select(abb, name)
TNRS::TNRS(
    taxonomic_names = data,
    sources = c("wfo", "tropicos", "wcvp"),
    classification = "wfo",
    mode = "resolve"
    ) %>%
  dplyr::select(
    ID, Name_submitted, Overall_score, Source, Accepted_name,
    Accepted_name_rank, Accepted_name_author, Accepted_family,
    Accepted_name_url, Taxonomic_status
  ) %>%
  filter(
    Overall_score < 1 | Source != "wfo" | Taxonomic_status != "Accepted"
    ) %>%
  arrange(Source)


### Combine with species_dikes table ###
traits <- traits %>%
  semi_join(species_dikes, by = "name")



## 4 Temperature and precipitation  ###########################################


sites_temperature <- read_csv(
  here("data", "raw", "temperature", "data", "data_OBS_DEU_P1M_T2M.csv"),
  col_names = TRUE, na = c("", "NA", "na"), col_types =
    cols(
      .default = "?"
    )
) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  select(site, date, value)

sites_precipitation <-  read_csv(
  here("data", "raw", "precipitation", "data", "data_OBS_DEU_P1M_RR.csv"),
  col_names = TRUE, na = c("", "NA", "na"), col_types =
    cols(
      .default = "?"
      )
) %>%
  rename(date = Zeitstempel, value = Wert, site = SDO_ID) %>%
  select(site, date, value)



## 5 Check data frames ########################################################


### Check typos ###
sites_dikes %>%
  filter(!str_detect(id, "_seeded$")) %>%
  janitor::tabyl(vegetation_cover)
#sites_dikes %>% filter(vegetation_cover == 17)
species_dikes %>%
  select(-name) %>%
  unlist() %>%
  janitor::tabyl()
species_dikes %>% # Check special typos
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  filter(value == 8)

### Compare vegetation_cover and accumulated_cover ###
species_dikes %>%
  summarise(across(where(is.double), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "id", values_to = "value") %>%
  mutate(id = factor(id)) %>%
  full_join(sites_dikes, by = "id") %>%
  mutate(diff = (value - vegetation_cover)) %>%
  select(id, survey_year, value, vegetation_cover, diff) %>%
  filter(diff > 50 | diff < -30) %>%
  arrange(diff) %>%
  print(n = 100)

### Check plots over time ###
species_dikes %>%
  select(name, starts_with("X01")) %>%
  filter(if_any(starts_with("X"), ~ . > 0)) %>%
  print(n = 70)

### Check missing data ###
miss_var_summary(sites_dikes, order = TRUE)
vis_miss(sites_dikes, cluster = FALSE, sort_miss = TRUE)
vis_miss(traits, cluster = FALSE, sort_miss = TRUE)


rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "sites_temperature", "sites_precipitation",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ###########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Create simple variables ###################################################


traits <- traits %>%
  mutate(
    lean_indicator = if_else(
      !(is.na(table30)) | !(is.na(table33)) | !(is.na(table34)), "yes", "no"
    ),
    target = if_else(
      target_herb == "yes" | target_grass == "yes", "yes", "no"
    ),
    ruderal = if_else(
      sociology >= 3300 & sociology < 3700, "yes", "no"
    ),
    target_ellenberg = if_else(
      sociology >= 5300 & sociology < 5400, "dry_grassland", if_else(
        sociology >= 5400 & sociology < 6000, "hay_meadow", if_else(
          sociology >= 5100 & sociology < 5120, "nardus_grassland", if_else(
            sociology >= 5200 & sociology < 5300, "sand_and_rock_vegetation",
            if_else(
              sociology >= 6100 & sociology < 6200, "fringe_vegetation", "no"
            )
          )
        )
      )
    )
  )

sites_dikes <- sites_dikes %>%
  mutate(
    n_total_concentration = finematerial_depth * finematerial_density * 10 *
      n_total_ratio / 100,
    plot_age = survey_year - construction_year
  )



## 2 Coverages #################################################################


cover <- species_dikes %>%
  left_join(traits, by = "name") %>%
  select(
    name, family, target, target_herb, target_arrhenatherion, lean_indicator,
    nitrogen_indicator, ruderal, table33, target_ellenberg,
    starts_with("X")
  ) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  group_by(id)

#### Graminoid, herb, and total coverage ###
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = TRUE), .groups = "keep") %>%
  mutate(
    type = if_else(
      family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae",
      "graminoid_cover", "herb_cover"
    )
  ) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = TRUE), .groups = "keep") %>%
  spread(type, total) %>%
  mutate(accumulated_cover = graminoid_cover + herb_cover) %>%
  ungroup()

#### German biotope mapping: Target species' coverage ###
cover_target <- cover %>%
  filter(target == "yes") %>%
  summarise(target_cover = sum(n, na.rm = TRUE)) %>%
  mutate(target_cover = round(target_cover, 1)) %>%
  ungroup()

#### German biotope mapping: Target herb species' coverage ###
cover_target_herb <- cover %>%
  filter(target_herb == "yes") %>%
  summarise(target_herb_cover = sum(n, na.rm = TRUE)) %>%
  mutate(target_herb_cover = round(target_herb_cover, 1)) %>%
  ungroup()

#### German biotope mapping: Arrhenatherum species' cover ratio ###
cover_target_arrhenatherion <- cover %>%
  filter(target_arrhenatherion == "yes") %>%
  summarise(arrh_cover = sum(n, na.rm = TRUE)) %>%
  mutate(arrh_cover = round(arrh_cover, 1)) %>%
  ungroup()

#### German biotope mapping: Lean indicator's coverage ###
cover_lean_indicator <- cover %>%
  filter(lean_indicator == "yes") %>%
  summarise(lean_cover = sum(n, na.rm = TRUE)) %>%
  mutate(lean_cover = round(lean_cover, 1)) %>%
  ungroup()

#### German biotope mapping: Nitrogen indicator's coverage ###
cover_nitrogen_indicator <- cover %>%
  filter(nitrogen_indicator == "yes") %>%
  summarise(nitrogen_cover = sum(n, na.rm = TRUE)) %>%
  mutate(nitrogen_cover = round(nitrogen_cover, 1)) %>%
  ungroup()

#### Ruderal's coverage ###
cover_ruderal <- cover %>%
  filter(ruderal == "yes") %>%
  summarise(ruderal_cover = sum(n, na.rm = TRUE)) %>%
  mutate(ruderal_cover = round(ruderal_cover, 1)) %>%
  ungroup()

#### German biotope mapping: Table 33 species' coverage ###
cover_table33 <- cover %>%
  mutate(table33 = if_else(table33 == "4" |
    table33 == "3" |
    table33 == "2",
  "table33_cover", "other"
  )) %>%
  filter(table33 == "table33_cover") %>%
  summarise(table33_cover = sum(n, na.rm = TRUE)) %>%
  mutate(table33_cover = round(table33_cover, 1)) %>%
  ungroup()

#### Ellenberg target species' coverage ###
cover_ellenberg <- cover %>%
  filter(target_ellenberg != "no") %>%
  summarise(ellenberg_cover = sum(n, na.rm = TRUE)) %>%
  mutate(ellenberg_cover = round(ellenberg_cover, 1)) %>%
  ungroup()

sites_dikes <- sites_dikes %>%
  right_join(cover_total_and_graminoid, by = "id") %>%
  right_join(cover_target, by = "id") %>%
  right_join(cover_target_herb, by = "id") %>%
  right_join(cover_target_arrhenatherion, by = "id") %>%
  right_join(cover_lean_indicator, by = "id") %>%
  right_join(cover_nitrogen_indicator, by = "id") %>%
  right_join(cover_ruderal, by = "id") %>%
  right_join(cover_table33, by = "id") %>%
  right_join(cover_ellenberg, by = "id") %>%
  mutate(
    ellenberg_cover_ratio = ellenberg_cover / accumulated_cover,
    graminoid_cover_ratio = graminoid_cover / accumulated_cover
  )

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "sites_temperature", "sites_precipitation",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 3 Alpha diversity ##########################################################


### a Species richness -------------------------------------------------------

richness <- left_join(species_dikes, traits, by = "name") %>%
  select(
    name, rlg, rlb, target, target_herb, target_arrhenatherion,
    target_ellenberg, ffh6510, ffh6210, nitrogen_indicator, lean_indicator,
    table33, table34, starts_with("X")
  ) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("X")) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  group_by(id)

#### Total species richness ###
richness_total <- richness %>%
  summarise(species_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Red list Germany (species richness) ###
richness_rlg <- richness %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  summarise(rlg_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Red list Bavaria (species richness) ###
richness_rlb <- richness %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  summarise(rlb_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### German biotope mapping: target species (species richness) ###
richness_target <- richness %>%
  filter(target == "yes") %>%
  summarise(target_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### German biotope mapping: target herb species (species richness) ###
richness_target_herb <- richness %>%
  filter(target_herb == "yes") %>%
  summarise(target_herb_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### German biotope mapping: Arrhenatherion species (species richness) ###
richness_arrh <- richness %>%
  filter(target_arrhenatherion == "yes") %>%
  summarise(arrh_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### FFH6510 species (species richness) ###
richness_ffh6510 <- richness %>%
  filter(ffh6510 == "yes") %>%
  summarise(ffh6510_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### FFH6210 species (species richness) ###
richness_ffh6210 <- richness %>%
  filter(ffh6210 == "yes") %>%
  summarise(ffh6210_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### German biotope mapping: lean_indicator species (species richness) ###
richness_lean_indicator <- richness %>%
  filter(lean_indicator == "yes") %>%
  summarise(lean_indicator_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### German biotope mapping: table33 (species richness) ###
richness_table33_2 <- richness %>%
  filter(table33 == "2") %>%
  summarise(table33_2_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()
richness_table33_3 <- richness %>%
  filter(table33 == "3") %>%
  summarise(table33_3_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()
richness_table33_4 <- richness %>%
  filter(table33 == "4") %>%
  summarise(table33_4_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### German biotope mapping: table34 (species richness) ###
richness_table34_2 <- richness %>%
  filter(table34 == "2") %>%
  summarise(table34_2_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()
richness_table34_3 <- richness %>%
  filter(table34 == "3") %>%
  summarise(table34_3_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Ellenberg target species (species richness) ###
richness_ellenberg <- richness %>%
  filter(target_ellenberg != "no") %>%
  summarise(ellenberg_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

sites_dikes <- sites_dikes %>%
  right_join(richness_total, by = "id") %>%
  right_join(richness_rlg, by = "id") %>%
  right_join(richness_rlb, by = "id") %>%
  right_join(richness_target, by = "id") %>%
  right_join(richness_target_herb, by = "id") %>%
  right_join(richness_arrh, by = "id") %>%
  right_join(richness_ffh6510, by = "id") %>%
  right_join(richness_ffh6210, by = "id") %>%
  right_join(richness_lean_indicator, by = "id") %>%
  right_join(richness_table33_2, by = "id") %>%
  right_join(richness_table33_3, by = "id") %>%
  right_join(richness_table33_4, by = "id") %>%
  right_join(richness_table34_2, by = "id") %>%
  right_join(richness_table34_3, by = "id") %>%
  right_join(richness_ellenberg, by = "id") %>%
  mutate(
    ellenberg_richness_ratio = ellenberg_richness / species_richness
    )


### b Species eveness and shannon ---------------------------------------------

data <- species_dikes %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id") %>%
  diversity(index = "shannon") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") %>%
  mutate(id = factor(id)) %>%
  rename(shannon = value)
sites_dikes <- sites_dikes %>%
  left_join(data, by = "id") %>%
  mutate(eveness = shannon / log(species_richness))

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "sites_temperature", "sites_precipitation",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 4 Biotope types ############################################################


### a Calculate types --------------------------------------------------------

data <- sites_dikes %>%
  select(
    id, table33_2_richness, table33_3_richness, table33_4_richness,
    table33_cover, table34_2_richness, table34_3_richness, target_richness,
    target_herb_richness, arrh_richness, target_cover, lean_cover, arrh_cover,
    target_herb_cover, nitrogen_cover
  ) %>%
  mutate(
    table33Rich_proof = if_else(
      table33_2_richness >= 2 | table33_2_richness + table33_3_richness >= 3 |
        table33_2_richness + table33_3_richness + table33_4_richness >= 4,
      "yes", "no"
    ),
    table33_cover_proof = if_else(
      table33_cover >= 25,
      "yes", "no"
    ),
    table33_proof = if_else(
      table33Rich_proof == "yes" & table33_cover_proof == "yes",
      "yes", "no"
    ),
    table34Rich_proof = if_else(
      table34_2_richness >= 2 | table34_2_richness + table34_3_richness >= 3,
      "yes", "no"
    ),
    G312_GT6210_type = if_else(
      table33_proof == "yes" & table34Rich_proof == "yes",
      "yes", "no"
    ),
    GE_proof = if_else(
      target_herb_cover >= 12.5 & target_richness >= 20 & nitrogen_cover < 25,
      "yes", "no"
    ),
    G214_GE6510_type = if_else(
      GE_proof == "yes" & arrh_richness >= 1 & arrh_cover > 0.5 &
        lean_cover >= 25 & target_herb_cover >= 12.5,
      "yes", "no"
    ),
    G212_LR6510_type = if_else(
      GE_proof == "yes" & arrh_richness >= 1 & arrh_cover > 0.5 &
        lean_cover < 25,
      "yes", "no"
    ),
    G214_GE00BK_type = if_else(
      GE_proof == "yes" & arrh_richness == 0 & lean_cover >= 25 &
        target_herb_cover >= 12.5,
      "yes", if_else(
        GE_proof == "yes" & arrh_cover >= 0.5 & lean_cover >= 25 &
          target_herb_cover >= 12.5,
        "yes", "no"
      )
    ),
    G213_GE00BK_type = if_else(
      GE_proof == "yes" & arrh_richness == 0 & lean_cover >= 25 &
        target_herb_cover < 12.5,
      "yes", if_else(
        GE_proof == "yes" & arrh_cover >= 0.5 & lean_cover >= 25 &
          target_herb_cover < 12.5,
        "yes", "no"
      )
    ),
    G213_type = if_else(
      lean_cover >= 25, "yes", "no"
    ),
    G212_type = if_else(
      lean_cover >= 1 & lean_cover < 25 & target_herb_cover >= 12.5 &
        target_herb_richness >= 10,
      "yes", "no"
    ),
    G211_type = if_else(
      lean_cover >= 1 & lean_cover < 25 & target_herb_cover >= 1 &
        target_herb_richness >= 5,
      "yes", "no"
    ),
    biotope_type = if_else(
      G312_GT6210_type == "yes", "G312-GT6210", if_else(
        G214_GE6510_type == "yes", "G214-GE6510", if_else(
          G212_LR6510_type == "yes", "G212-LR6510", if_else(
            G214_GE00BK_type == "yes", "G214-GE00BK", if_else(
              G213_GE00BK_type == "yes", "G213-GE00BK", if_else(
                G213_type == "yes", "G213", if_else(
                  G212_type == "yes", "G212", if_else(
                    G211_type == "yes",
                    "G211", "other"
                  )
                )
              )
            )
          )
        )
      )
    ),
    biotope_type = as_factor(biotope_type)
  ) %>%
  select(id, biotope_type, -ends_with("proof"), -starts_with("table")) %>%
  mutate(
    ffh6510 = str_match(biotope_type, "6510"),
    ffh6210 = str_match(biotope_type, "6210"),
    baykompv = as_factor(str_sub(biotope_type, start = 1, end = 4))
  ) %>%
  unite(ffh, ffh6510, ffh6210, sep = "") %>%
  mutate(
    ffh = str_replace(ffh, "NA", ""),
    ffh = as_factor(str_replace(ffh, "NA", "non-FFH")),
    biotope_type = as_factor(biotope_type),
    biotope_points = if_else(
      biotope_type == "G312-GT6210", 13, if_else(
        biotope_type == "G214-GE6510", 12, if_else(
          biotope_type == "G212-LR6510", 9, if_else(
            biotope_type == "G214-GE00BK", 12, if_else(
              biotope_type == "G213-GE00BK", 9, if_else(
                biotope_type == "G213", 8, if_else(
                  biotope_type == "G212", 8, if_else(
                    biotope_type == "G211",
                    6, 0
                  )
                )
              )
            )
          )
        )
      )
    )
  )

sites_dikes <- sites_dikes %>%
  left_join(data, by = "id") %>%
  select(
    -target_herb_cover, -arrh_cover, -lean_cover, -nitrogen_cover,
    -table33_cover,
    -target_herb_richness, -arrh_richness, -lean_indicator_richness,
    -ffh6510_richness, -ffh6210_richness, -table33_2_richness,
    -table33_3_richness, -table33_4_richness, -table34_2_richness,
    -table34_3_richness
  )

traits <- traits %>%
  select(
    -target_arrhenatherion, -table30, -table33, -table34,
    -nitrogen_indicator, -nitrogen_indicator2, -lean_indicator,
    -grazing_indicator, -rlb, -target_herb, -target_grass,
    -sociology, -legal
  )


### b Calculate constancy of ffh evaluation -----------------------------------

data <- sites_dikes %>%
  select(id, plot, survey_year, ffh) %>%
  group_by(plot) %>%
  mutate(count = n()) %>%
  filter(count == max(count)) %>%
  pivot_wider(
    id_cols = -id, names_from = "survey_year", values_from = "ffh"
  ) %>%
  group_by(plot) %>%
  rename(x17 = "2017", x18 = "2018", x19 = "2019") %>%
  mutate(
    change_type = ifelse(
      (x17 == "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH"),
      "better", ifelse(
        (x17 != "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") |
          (x17 == "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH") |
          (x17 == "non-FFH" & x18 == "non-FFH" & x19 != "non-FFH") |
          (x17 != "non-FFH" & x18 != "non-FFH" & x19 == "non-FFH"),
        "change", ifelse(
          (x17 == "non-FFH" & x18 == "non-FFH" & x19 == "non-FFH"),
          "non-FFH", ifelse(
            x17 == "6510" & x18 == "6510" & x19 == "6510",
            "FFH6510", ifelse(
              x17 == "6210" & x18 == "6210" & x19 == "6210",
              "FFH6210", ifelse(
                x17 != "non-FFH" & x18 != "non-FFH" & x19 != "non-FFH",
                "any-FFH", "worse"
                )
              )
            )
          )
        )
      )
    ) %>%
  select(plot, change_type)

#sites_dikes <- left_join(sites_dikes, data, by = "plot")
# --> Data not used ###

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "sites_temperature", "sites_precipitation",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 5 LCBD: Local contributions to beta diversity ##############################


data_sites <- sites_dikes %>%
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n))
data_species <- species_dikes %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))

#### * 2017 ###
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2017")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2017 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none

#### * 2018 ###
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2018")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2018 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none

#### * 2019 ###
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2019")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2019 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none

#### * 2021 ###
data_species_lcbd <- data_species %>%
  filter(str_detect(id, "2021")) %>%
  column_to_rownames(var = "id")
data <- beta.div(data_species_lcbd, method = "hellinger", nperm = 9999)
data_2021 <- data$LCBD
row.names(data_species_lcbd[which(p.adjust(data$p.LCBD, "holm") <= 0.05), ])
#--> none


data <- c(data_2017, data_2018, data_2019, data_2021)
data <- data_sites %>%
  mutate(lcbd = data) %>%
  select(id, lcbd)
#sites_dikes <- left_join(sites_dikes, data, by = "id")
# --> Data not used ###

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "sites_temperature", "sites_precipitation",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 6 PCAs: Principle component analyses #######################################


### a Soil PCA  ---------------------------------------------------------------

data <- sites_dikes %>%
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n) & survey_year == 2017) %>%
  select(
    plot, ph,
    calciumcarbonat, humus, n_total_ratio, n_total_concentration, cn_ratio,
    sand, silt, clay,
    phosphorus, potassium, magnesium,
    topsoil_depth
  ) %>%
  column_to_rownames(var = "plot")

#### * PCA ####
pca <- rda(data, scale = TRUE)

summary(pca, axes = 0)
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))
biplot(pca, display = "species", scaling = 2)


summary_table_part_1 <- pca %>%
  summary() %>%
  magrittr::extract2("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC4, variables)
summary_table_part_2 <- pca %>%
  eigenvals() %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC4, variables)
pca_soil <- bind_rows(summary_table_part_1, summary_table_part_2) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))

data <- pca %>%
  summary() %>%
  magrittr::extract2("sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot") %>%
  tibble() %>%
  select(plot, PC1:PC4) %>%
  rename(
    pc1_soil = PC1,
    pc2_soil = PC2,
    pc3_soil = PC3,
    pc4_soil = PC4
  )
sites_dikes <- left_join(sites_dikes, data, by = "plot") %>%
  select(
    -c550, -calciumcarbonat, -c_organic, -humus, -humus_level,
    -n_total_ratio, -cn_ratio, -ph, -sand, -silt, -clay, -phosphorus,
    -potassium, -magnesium, -topsoil_depth, -finematerial_depth,
    -finematerial_density, -sceleton_ratio_volume, -sceleton_ratio_mass,
    -sceleton_level, -ufc_percent, -ufc_mm, -n_total_concentration
  )

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "sites_temperature", "sites_precipitation",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)


### b Climate PCA - Preparation  ----------------------------------------------

#### * Temperature ####

data <- sites_temperature %>%
  filter(date >= "2002-03-01") %>%
  mutate(
    site = factor(site),
    season = floor_date(date, "season"),
    year = year(season),
    season = month(season),
    season = factor(season),
    season = fct_recode(
      season,
      "spring" = "3",
      "summer" = "6",
      "autumn" = "9",
      "winter" = "12"
    )
  ) %>%
  group_by(year) %>%
  mutate(
    year_mean = mean(value),
    year_mean = round(year_mean, digits = 1),
    currentYear = if_else(season == "spring", 0, 1),
    currentYear = year + currentYear
  ) %>%
  # current year of 2021 is e.g. from summer 2020 to spring 2021 =
  # climate for survey_year
  group_by(currentYear) %>%
  mutate(
    current_mean = mean(value),
    current_mean = round(current_mean, digits = 1),
    current_mean = if_else(season == "spring", current_mean, NA_real_)
  ) %>%
  group_by(year, season, year_mean, currentYear, current_mean) %>%
  summarise(season_mean = round(mean(value), digits = 1), .groups = "keep") %>%
  pivot_wider(
    id_cols = c(year, year_mean, current_mean),
    names_from = season, values_from = season_mean
  ) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~ max(., na.rm = TRUE)),
            .groups = "keep") %>%
  # warnings because of lates year (summer, autumn, winter), can be ignored
  mutate(year = factor(year))

sites_dikes <- sites_dikes %>%
  left_join(data %>% select(-current_mean),
    by = c("construction_year_factor" = "year")
  ) %>%
  rename(
    "temp_spring_construction_year" = "spring",
    "temp_summer_construction_year" = "summer",
    "temp_autumn_construction_year" = "autumn",
    "temp_winter_construction_year" = "winter",
    "temp_mean_construction_year" = "year_mean"
  ) %>%
  left_join(data %>% select(-current_mean),
    by = c("construction_year_factor_plus" = "year")
  ) %>%
  rename(
    "temp_spring_construction_year_plus" = "spring",
    "temp_summer_construction_year_plus" = "summer",
    "temp_autumn_construction_year_plus" = "autumn",
    "temp_winter_construction_year_plus" = "winter",
    "temp_mean_construction_year_plus" = "year_mean"
  ) %>%
  left_join(data %>% select(year, current_mean, spring),
    by = c("survey_year_factor" = "year")
  ) %>%
  rename(
    "temp_mean_survey_year" = "current_mean",
    "temp_spring_survey_year" = "spring"
  ) %>%
  left_join(data %>% select(year, summer, autumn, winter),
    by = c("survey_year_factor_minus" = "year")
  ) %>%
  rename(
    "temp_summer_survey_year" = "summer",
    "temp_autumn_survey_year" = "autumn",
    "temp_winter_survey_year" = "winter"
  )

#### * Precipitation ####

data <- sites_precipitation %>%
  filter(date >= "2002-03-01") %>%
  mutate(
    site = factor(site),
    season = floor_date(date, "season"),
    year = year(season),
    season = month(season),
    season = factor(season),
    season = fct_recode(
      season,
      "spring" = "3",
      "summer" = "6",
      "autumn" = "9",
      "winter" = "12"
    )
  ) %>%
  group_by(site, year) %>%
  mutate(
    yearSum = round(sum(value), digits = 0),
    currentYear = if_else(season == "spring", 0, 1),
    currentYear = year + currentYear
  ) %>%
  # current year of 2021 is e.g. from summer 2020 to spring 2021 =
  # climate for survey_year
  group_by(site, currentYear) %>%
  mutate(
    currentSum = sum(value),
    currentSum = round(currentSum, 0)
  ) %>%
  group_by(site, season, year, currentYear, yearSum, currentSum) %>%
  summarise(
    seasonSum = round(sum(value), digits = 0),
    .groups = "keep"
  ) %>%
  group_by(year, season, currentYear) %>%
  summarise(
    season_mean = round(mean(seasonSum), digits = 0),
    year_mean = round(mean(yearSum), digits = 0),
    current_mean = round(mean(currentSum), digits = 0),
    .groups = "keep"
  ) %>%
  mutate(current_mean = if_else(season == "spring", current_mean, NA_real_)) %>%
  pivot_wider(id_cols = c(year, year_mean, current_mean),
              names_from = season, values_from = season_mean) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), ~ max(., na.rm = TRUE)),
    .groups = "keep"
  ) %>%
  # warnings because of lates year (summer, autumn, winter), can be ignored
  mutate(year = factor(year))

sites_dikes <- sites_dikes %>%
  left_join(data %>% select(-current_mean),
    by = c("construction_year_factor" = "year")
  ) %>%
  rename(
    "prec_spring_construction_year" = "spring",
    "prec_summer_construction_year" = "summer",
    "prec_autumn_construction_year" = "autumn",
    "prec_winter_construction_year" = "winter",
    "prec_mean_construction_year" = "year_mean"
  ) %>%
  left_join(data %>% select(-current_mean),
    by = c("construction_year_factor_plus" = "year")
  ) %>%
  rename(
    "prec_spring_construction_year_plus" = "spring",
    "prec_summer_construction_year_plus" = "summer",
    "prec_autumn_construction_year_plus" = "autumn",
    "prec_winter_construction_year_plus" = "winter",
    "prec_mean_construction_year_plus" = "year_mean"
  ) %>%
  left_join(data %>% select(year, current_mean, spring),
    by = c("survey_year_factor" = "year")
  ) %>%
  rename(
    "prec_mean_survey_year" = "current_mean",
    "prec_spring_survey_year" = "spring"
  ) %>%
  left_join(data %>% select(year, summer, autumn, winter),
    by = c("survey_year_factor_minus" = "year")
  ) %>%
  rename(
    "prec_summer_survey_year" = "summer",
    "prec_autumn_survey_year" = "autumn",
    "prec_winter_survey_year" = "winter"
  )


### c Climate PCA - Survey year  ----------------------------------------------

data <- sites_dikes %>%
  arrange(id) %>%
  column_to_rownames("id") %>%
  select(
    temp_mean_survey_year, temp_spring_survey_year, temp_summer_survey_year,
    temp_autumn_survey_year, temp_winter_survey_year,
    prec_mean_survey_year, prec_spring_survey_year, prec_summer_survey_year,
    prec_autumn_survey_year, prec_winter_survey_year
  )

#### * PCA ####
pca <- rda(X = decostand(data, method = "standardize"), scale = TRUE)


summary(pca, axes = 0)
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))
biplot(pca, display = "species", scaling = 2)


summary_table_part_1 <- pca %>%
  summary() %>%
  magrittr::extract2("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
summary_table_part_2 <- pca %>%
  eigenvals() %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
pca_survey_year <- bind_rows(summary_table_part_1, summary_table_part_2) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))


data <- pca %>%
  summary() %>%
  magrittr::extract2("sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot") %>%
  tibble() %>%
  select(plot, PC1:PC3) %>%
  rename(
    pc1_survey_year = PC1,
    pc2_survey_year = PC2,
    pc3_survey_year = PC3
  )
#sites_dikes <- left_join(sites_dikes, data, by = "plot")
# Data not used


### d Climate PCA - Construction year  ----------------------------------------

data <- sites_dikes %>%
  arrange(id) %>%
  select(
    plot,
    temp_mean_construction_year, temp_spring_construction_year,
    temp_summer_construction_year, temp_autumn_construction_year,
    temp_winter_construction_year,
    temp_mean_construction_year_plus, temp_spring_construction_year_plus,
    temp_summer_construction_year_plus, temp_autumn_construction_year_plus,
    temp_winter_construction_year_plus,
    prec_mean_construction_year, prec_spring_construction_year,
    prec_summer_construction_year, prec_autumn_construction_year,
    prec_winter_construction_year,
    prec_mean_construction_year_plus, prec_spring_construction_year_plus,
    prec_summer_construction_year_plus, prec_autumn_construction_year_plus,
    prec_winter_construction_year_plus
  ) %>%
  group_by(plot) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
  column_to_rownames("plot")

#### * PCA ####
pca <- rda(X = decostand(data, method = "standardize"), scale = TRUE)


summary(pca, axes = 0)
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))
biplot(pca, display = "species", scaling = 2)


summary_table_part_1 <- pca %>%
  summary() %>%
  magrittr::extract2("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
summary_table_part_2 <- pca %>%
  eigenvals() %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "variables") %>%
  tibble() %>%
  select(PC1:PC3, variables)
pca_construction_year <- bind_rows(
  summary_table_part_1, summary_table_part_2
  ) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3)))


data <- pca %>%
  summary() %>%
  magrittr::extract2("sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot") %>%
  tibble() %>%
  select(plot, PC1:PC3) %>%
  rename(
    pc1_construction_year = PC1,
    pc2_construction_year = PC2,
    pc3_construction_year = PC3
  )
sites_dikes <- left_join(sites_dikes, data, by = "plot") %>%
  select(-starts_with("temp"), -starts_with("prec"))

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 7 ESy: EUNIS expert vegetation classification system #######################


#### Start ###
### Bruelheide et al. 2021 Appl Veg Sci
### https://doi.org/10.1111/avsc.12562

expertfile <- "EUNIS-ESy-2020-06-08.txt" ### file of 2021 is not working

obs <- species_dikes %>%
  pivot_longer(cols = -name,
               names_to = "RELEVE_NR",
               values_to = "Cover_Perc") %>%
  rename(TaxonName = "name") %>%
  mutate(
    TaxonName = str_replace_all(TaxonName, "_", " "),
    TaxonName = str_replace_all(TaxonName, "ssp", "subsp."),
    TaxonName = as.factor(TaxonName),
    TaxonName = fct_recode(
      TaxonName,
      "Carex praecox" = "Carex praecox subsp. curvata",
      "Cerastium fontanum" = "Cerastium fontanum subsp. vulgare",
      "Clinopodium acinos" = "Acinos arvensis",
      "Ranunculus polyanthemos" = "Ranunculus serpens subsp. nemorosus",
      "Silene latifolia" = "Silene latifolia subsp. alba",
      "Vicia villosa" = "Vicia villosa subsp. varia"
      )
    ) %>%
  data.table::as.data.table()

header <- sites_dikes %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 31468) %>%
  sf::st_transform(4326) %>%
  rename(
    RELEVE_NR = id
    ) %>%
  mutate(
    "Altitude (m)" = 313,
    Latitude = sf::st_coordinates(.)[, 2],
    Longitude = sf::st_coordinates(.)[, 1],
    Country = "Germany",
    Coast_EEA = "N_COAST",
    Dunes_Bohn = "N_DUNES",
    Ecoreg = 686,
    dataset = "Danube_dikes"
    ) %>%
  select(RELEVE_NR, "Altitude (m)", Latitude, Longitude, Country,
         Coast_EEA, Dunes_Bohn, Ecoreg, dataset) %>%
  sf::st_drop_geometry()

setwd(here("R", "esy"))
source(here("R", "esy", "code", "prep.R"))

#### Step 1 and 2: Load and parse the expert file ###
source(here("R", "esy", "code", "step1and2_load-and-parse-the-expert-file.R"))

#### Step 3: Create a numerical plot x membership condition matrix  ###
plot.cond <- array(
  0,
  c(length(unique(obs$RELEVE_NR)), length(conditions)),
  dimnames = list(
    as.character(unique(obs$RELEVE_NR)),
    conditions
    )
  )

### Step 4: Aggregate taxon levels ###
source(here("R", "esy", "code", "step4_aggregate-taxon-levels.R"))

(data <- obs %>%
  group_by(TaxonName) %>%
  slice(1) %>%
  anti_join(AGG, by = c("TaxonName" = "ind")))

#### Step 5: Solve the membership conditions ###
mc <- 1
source(
  here(
    "R", "esy", "code", "step3and5_extract-and-solve-membership-conditions.R"
    )
  )

table(result.classification)
eval.EUNIS(which(result.classification == "V39")[1], "V39")

sites_dikes <- sites_dikes %>%
  mutate(
    esy = result.classification,
    esy = if_else(id == "X05_m_2021", "R1A", esy),
    esy = if_else(id == "X62_m_2019", "R", esy),
    esy = if_else(id == "X67_o_2021", "R", esy)
    )
table(sites_dikes$esy)

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 8 dbMEM: Distance-based Moran's eigenvector maps (41 plots) ################


# Borcard et al. (2018) Numerical Ecology https://doi.org/10.1007/978-1-4419-7976-6
# function 'quickMEM' from Numerical Ecology textbook:
source("https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R")
data_sites <- sites_dikes %>%
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n))
data_species <- species_dikes %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(id, names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))


### a 2017 ---------------------------------------------------------------------

data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2017) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = TRUE) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")


m <- quickMEM(
  data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  rangexy = TRUE,
  perm.max = 999
)


# --> undetrended data, R2adj of minimum (final) model = 0.05
m$RDA_test # p: 0.001
m$RDA_axes_test # RDA1: sig. axis
m$RDA # RDA1: 0.05

data <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1) %>%
  rename(mem1_2017 = MEM1)
sites_dikes <- left_join(sites_dikes, data, by = "id")


### b 2018 ---------------------------------------------------------------------

data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2018) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = TRUE) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")


### Try several times! ###
m <- quickMEM(
  data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  rangexy = TRUE,
  detrend = FALSE,
  perm.max = 999
)


# --> undetrended data, R2adj of minimum (final) model = 0.04
m$RDA_test # p: 0.006
m$RDA_axes_test #  RDA2: sig. axis
m$RDA # RDA2: 0.03

data <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM2) %>%
  rename(mem2_2018 = MEM2)
sites_dikes <- left_join(sites_dikes, data, by = "id")


### c 2019 ---------------------------------------------------------------------

data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2019) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = TRUE) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")


m <- quickMEM(
  data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  detrend = FALSE,
  rangexy = TRUE,
  perm.max = 999
)


# --> not detrended, R2adj of minimum (final) model = 0.056
m$RDA_test # p: 0.001
m$RDA_axes_test # RDA1: sig axis
m$RDA # RDA1: 0.04

data <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1) %>%
  rename(mem1_2019 = MEM1)
sites_dikes <- left_join(sites_dikes, data, by = "id")


### d 2021 ---------------------------------------------------------------------

data_sites_dbmem <- data_sites %>%
  filter(survey_year == 2021) %>%
  select(id, longitude, latitude)
data_species_dbmem <- data_species %>%
  semi_join(data_sites_dbmem, by = "id") %>%
  column_to_rownames(var = "id") %>%
  decostand("hellinger")
i <- (colSums(data_species_dbmem, na.rm = TRUE) != 0)
data_species_dbmem <- data_species_dbmem[, i]
data_sites_dbmem <- data_sites_dbmem %>%
  column_to_rownames("id")


m <- quickMEM(
  data_species_dbmem, data_sites_dbmem,
  alpha = 0.05,
  method = "fwd",
  detrend = FALSE,
  rangexy = TRUE,
  perm.max = 999
)


# --> undetrended data, R2adj of minimum (final) model: 0.064
m$RDA_test # p: 0.001
m$RDA_axes_test # RDA1, RDA2: sig axes
m$RDA # RDA1: 0.05, RDA2: 0.04

data <- m$dbMEM %>%
  rownames_to_column(var = "id") %>%
  select(id, MEM1, MEM2) %>%
  rename(mem1_2021 = MEM1, mem2_2021 = MEM2)
sites_dikes <- left_join(sites_dikes, data, by = "id")

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits",
      "pca_soil", "pca_construction_year", "pca_survey_year"
    )
  )
)



## 9 TBI: Temporal Beta diversity Index ######################################


### Preparation ###
data_sites <- sites_dikes %>%
  ### Choose only plots which were surveyed in each year
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(id, plot)
data_species <- species_dikes %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  left_join(traits %>% select(name, target_ellenberg), by = "name")
data_species_all <- data_species %>%
  pivot_longer(-c("name", "target_ellenberg"),
               names_to = "id", values_to = "value") %>%
  pivot_wider(id, names_from = "name", values_from = "value") %>%
  mutate(
    year = str_match(id, "\\d{4}"),
    plot = str_match(id, "\\d{2}"),
    year = factor(year),
    plot = factor(plot)
  ) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  select(plot, year, tidyselect::peek_vars(), -id)
data_species_target <- data_species %>%
  filter(target_ellenberg != "no") %>%
  pivot_longer(-c("name", "target_ellenberg"),
               names_to = "id",
               values_to = "value") %>%
  pivot_wider(id, names_from = "name", values_from = "value") %>%
  mutate(
    year = str_match(id, "\\d{4}"),
    plot = str_match(id, "\\d{2}"),
    year = factor(year),
    plot = factor(plot)
  ) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  select(plot, year, tidyselect::peek_vars(), -id)

for (i in unique(data_species_all$year)) {
  nam <- paste("data_species_all_", i, sep = "")

  assign(nam, data_species_all %>%
    filter(year == i) %>%
    select(-year) %>%
    column_to_rownames(var = "plot"))
}

for (i in unique(data_species_target$year)) {
  nam <- paste("data_species_target_", i, sep = "")

  assign(nam, data_species_target %>%
           filter(year == i) %>%
           select(-year) %>%
           column_to_rownames(var = "plot"))
}


### a Calculate TBI Presence - All species -------------------------------------

#### * 2017 vs. 2018 ###
res1718 <- TBI(
  data_species_all_2017, data_species_all_2018,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1718$BCD.summary # B: 0.223, C: 0.155, D: 0.378 (58.9% vs. 41.0%)
res1718$t.test_B.C # p.perm: 0.0058
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ###
res1819 <- TBI(
  data_species_all_2018, data_species_all_2019,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1819$BCD.summary # B: 0.118, C: 0.214, D: 0.332 (35.6% vs. 64.3%)
res1819$t.test_B.C # p.perm: 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819")
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ###
res1921 <- TBI(
  data_species_all_2019, data_species_all_2021,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1921$BCD.summary # B: 0.249, C: 0.140, D: 0.390 (63.8% vs. 36.1%)
res1921$t.test_B.C # p.perm: 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ###
res1719 <- TBI(
  data_species_all_2017, data_species_all_2019,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1719$BCD.summary # B: 0.186, C: 0.212, D: 0.399 (46.7% vs. 53.2%)
res1719$t.test_B.C # p.perm: 0.273
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ###
res1721 <- TBI(
  data_species_all_2017, data_species_all_2021,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1721$BCD.summary # B: 0.184, C: 0.450, D: 0.590 (59.0% vs. 40.9%)
res1721$t.test_B.C # p.perm: 0.0021
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
plot(res1721, type = "BC")


data_presence_all <- bind_rows(tbi1718, tbi1819, tbi1921, tbi1719, tbi1721) %>%
  mutate(presabu = "presence",
         pool = "all")


### b Calculate TBI Abundance - All species ------------------------------------

#### * 2017 vs. 2018 ###
res1718 <- TBI(
  data_species_all_2017, data_species_all_2018,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1718$BCD.summary # B: 0.213, C: 0.260, D: 0.473 (45.0% vs. 54.9%)
res1718$t.test_B.C # p.perm: 0.1756
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ###
res1819 <- TBI(
  data_species_all_2018, data_species_all_2019,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1819$BCD.summary # B: 0.167, C: 0.302, D: 0.470 (35.7% vs. 64.2%)
res1819$t.test_B.C # p.perm: 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819")
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ###
res1921 <- TBI(
  data_species_all_2019, data_species_all_2021,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1921$BCD.summary # B: 0.331, C: 0.168, D: 0.499 (66.3% vs. 33.6%)
res1921$t.test_B.C # p.perm: 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ###
res1719 <- TBI(
  data_species_all_2017, data_species_all_2019,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1719$BCD.summary # B: 0.210, C: 0.390, D: 0.601 (35.0% vs. 64.9%)
res1719$t.test_B.C # p.perm: 1e-04
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ###
res1721 <- TBI(
  data_species_all_2017, data_species_all_2021,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
  )
res1721$BCD.summary # B: 0.301, C: 0.319, D: 0.620 (48.5% vs. 51.4%)
res1721$t.test_B.C # p.perm: 0.598
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
plot(res1721, type = "BC")


data_abundance_all <- bind_rows(tbi1718, tbi1819, tbi1921, tbi1719, tbi1721) %>%
  mutate(presabu = "abundance",
         pool = "all")


### c Calculate TBI Presence - Target species ----------------------------------

#### * 2017 vs. 2018 ###
res1718 <- TBI(
  data_species_target_2017, data_species_target_2018,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1718$BCD.summary # B: 0.223, C: 0.155, D: 0.378 (58.9% vs. 41.0%)
res1718$t.test_B.C # p.perm: 0.0058
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ###
res1819 <- TBI(
  data_species_target_2018, data_species_target_2019,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1819$BCD.summary # B: 0.118, C: 0.214, D: 0.332 (35.6% vs. 64.3%)
res1819$t.test_B.C # p.perm: 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819")
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ###
res1921 <- TBI(
  data_species_target_2019, data_species_target_2021,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1921$BCD.summary # B: 0.249, C: 0.140, D: 0.390 (63.8% vs. 36.1%)
res1921$t.test_B.C # p.perm: 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ###
res1719 <- TBI(
  data_species_target_2017, data_species_target_2019,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1719$BCD.summary # B: 0.186, C: 0.212, D: 0.399 (46.7% vs. 53.2%)
res1719$t.test_B.C # p.perm: 0.273
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ###
res1721 <- TBI(
  data_species_target_2017, data_species_target_2021,
  method = "sorensen",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1721$BCD.summary # B: 0.184, C: 0.450, D: 0.590 (59.0% vs. 40.9%)
res1721$t.test_B.C # p.perm: 0.0021
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
plot(res1721, type = "BC")


data_presence_target <- bind_rows(
  tbi1718, tbi1819, tbi1921, tbi1719, tbi1721
  ) %>%
  mutate(presabu = "presence",
         pool = "target")


### d Calculate TBI Abundance - Target species --------------------------------

#### * 2017 vs. 2018 ###
res1718 <- TBI(
  data_species_target_2017, data_species_target_2018,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1718$BCD.summary # B: 0.213, C: 0.260, D: 0.473 (45.0% vs. 54.9%)
res1718$t.test_B.C # p.perm: 0.1756
tbi1718 <- as_tibble(res1718$BCD.mat) %>%
  mutate(comparison = "1718")
plot(res1718, type = "BC")

#### * 2018 vs. 2019 ###
res1819 <- TBI(
  data_species_target_2018, data_species_target_2019,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1819$BCD.summary # B: 0.167, C: 0.302, D: 0.470 (35.7% vs. 64.2%)
res1819$t.test_B.C # p.perm: 1e-04
tbi1819 <- as_tibble(res1819$BCD.mat) %>%
  mutate(comparison = "1819")
plot(res1819, type = "BC")

#### * 2019 vs. 2021 ###
res1921 <- TBI(
  data_species_target_2019, data_species_target_2021,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1921$BCD.summary # B: 0.331, C: 0.168, D: 0.499 (66.3% vs. 33.6%)
res1921$t.test_B.C # p.perm: 1e-04
tbi1921 <- as_tibble(res1921$BCD.mat) %>%
  mutate(comparison = "1921")
plot(res1921, type = "BC")

#### * 2017 vs. 2019 ###
res1719 <- TBI(
  data_species_target_2017, data_species_target_2019,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1719$BCD.summary # B: 0.210, C: 0.390, D: 0.601 (35.0% vs. 64.9%)
res1719$t.test_B.C # p.perm: 1e-04
tbi1719 <- as_tibble(res1719$BCD.mat) %>%
  mutate(comparison = "1719")
plot(res1719, type = "BC")

#### * 2017 vs. 2021 ###
res1721 <- TBI(
  data_species_target_2017, data_species_target_2021,
  method = "%diff",
  nperm = 9999, test.t.perm = TRUE, clock = TRUE
)
res1721$BCD.summary # B: 0.301, C: 0.319, D: 0.620 (48.5% vs. 51.4%)
res1721$t.test_B.C # p.perm: 0.598
tbi1721 <- as_tibble(res1721$BCD.mat) %>%
  mutate(comparison = "1721")
plot(res1721, type = "BC")


data_abundance_target <- bind_rows(
  tbi1718, tbi1819, tbi1921, tbi1719, tbi1721
  ) %>%
  mutate(presabu = "abundance",
         pool = "target")


### e Combine all data ---------------------------------------------------------

plot <- data_sites %>%
  filter(str_detect(id, "2017")) %>%
  pull(plot)
data <- add_row(data_presence_target, data_abundance_target) %>%
  add_row(data_presence_all) %>%
  add_row(data_abundance_all) %>%
  mutate(plot = rep(
    plot,
    length(data_abundance_all$comparison) *
      length(unique(sites_dikes$survey_year)) /
      length(data_species_all_2017[, 1])
    ))
sites_temporal <- sites_dikes %>%
  filter(survey_year_factor == "2017") %>%
  left_join(data, by = "plot") %>%
  rename(
    b = "B/(2A+B+C)", c = "C/(2A+B+C)", d = "D=(B+C)/(2A+B+C)",
  ) %>%
  mutate(
    plot = as_factor(plot)
  )

rm(
  list = setdiff(
    ls(), c(
      "sites_dikes", "species_dikes", "traits", "pca_soil",
      "pca_construction_year", "pca_survey_year"
    )
  )
)



## 10 Finalization ############################################################


### a Rounding ----------------------------------------------------------------

sites_dikes <- sites_dikes %>%
  mutate(across(
    c(
      pc1_soil, pc2_soil, pc3_soil,
      pc1_construction_year, pc2_construction_year, pc3_construction_year
    ),
    ~ round(.x, digits = 4)
  )) %>%
  mutate(across(
    c(
      ellenberg_cover_ratio, graminoid_cover_ratio,
      ellenberg_richness_ratio, shannon, eveness
    ),
    ~ round(.x, digits = 3)
  )) %>%
  mutate(across(
    c(river_distance, accumulated_cover),
    ~ round(.x, digits = 1)
  ))

sites_temporal <- sites_temporal %>%
  mutate(across(
    c(pc1_soil, pc2_soil, pc3_soil, river_distance, b, c, d),
    ~ round(.x, digits = 4)
  ))


### b Final selection of variables --------------------------------------------

sites_spatial <- sites_dikes %>%
  select(
    id, plot, block, latitude, longitude, survey_year, botanist, esy,
    # local
    exposition, orientation, pc1_soil, pc2_soil, pc3_soil,
    # space
    location, location_abb, location_construction_year, latitude, longitude,
    river_km, river_distance, biotope_distance, biotope_area,
    mem1_2017, mem2_2018, mem1_2019, mem1_2021, mem2_2021,
    # history
    construction_year, plot_age,
    pc1_construction_year, pc2_construction_year, pc3_construction_year,
    # response variables
    accumulated_cover, species_richness, eveness, shannon,
    graminoid_cover_ratio, ruderal_cover, ellenberg_cover,
    ellenberg_richness, ellenberg_richness_ratio, ellenberg_cover_ratio
  )

sites_temporal <- sites_temporal %>%
  select(
    plot, block, comparison, pool, presabu,
    # local
    exposition, orientation, pc1_soil, pc2_soil, pc3_soil,
    # space
    location, location_abb, location_construction_year,
    river_km, river_distance, biotope_distance, biotope_area,
    # history
    construction_year,
    # temporal beta-diversity indices (TBI)
    b, c, d
  )

sites_restoration <- sites_dikes %>%
  select(
    id, plot, block, latitude, longitude, botanist, survey_year, esy,
    # local
    exposition, orientation, pc1_soil, pc2_soil, pc3_soil,
    # space
    location, location_abb, location_construction_year,
    river_km, river_distance, biotope_distance, biotope_area,
    # history
    construction_year, plot_age,
    # response variables
    species_richness, eveness,
    accumulated_cover, graminoid_cover_ratio, ruderal_cover,
    ellenberg_richness, ellenberg_richness_ratio,
    rlg_richness, ellenberg_cover_ratio, ellenberg_cover,
    # legal evaluation
    biotope_type, ffh, baykompv, biotope_points
  )


### c Final selection of plots ------------------------------------------------

sites_spatial <- sites_spatial %>%
  ### Choose only plots which were surveyed in each year
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(-n) %>%
  mutate(plot = factor(plot)) # check number of plots

species_dikes <- species_dikes %>%
  ### remove species with no occurence
  group_by(name) %>%
  arrange(name) %>%
  select(name, all_of(sites_spatial$id)) %>%
  mutate(
    total = sum(c_across(starts_with("X")), na.rm = TRUE),
    presence = if_else(total > 0, 1, 0),
    name = factor(name)
  ) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence)


sites_temporal <- sites_temporal %>%
  semi_join(sites_spatial, by = "plot") %>%
  mutate(plot = factor(plot)) # check number of plots



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data #######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### * Sites data ####
write_csv(
  sites_spatial,
  here("data", "processed", "data_processed_sites_spatial.csv")
)
write_csv(
  sites_temporal,
  here("data", "processed", "data_processed_sites_temporal.csv")
)
write_csv(
  sites_restoration,
  here("data", "processed", "data_processed_sites_for_practioners.csv")
)

### * Species and traits data ####
write_csv(
  species_dikes,
  here("data", "processed", "data_processed_species.csv")
)
write_csv(
  traits,
  here("data", "processed", "data_processed_traits.csv")
)

### * PCA tables ####
write_csv(
  pca_soil,
  here("outputs", "statistics", "pca_soil.csv")
)
write_csv(
  pca_construction_year,
  here("outputs", "statistics", "pca_construction_year.csv")
)
