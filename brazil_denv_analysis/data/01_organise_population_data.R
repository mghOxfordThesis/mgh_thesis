#################################################################################################
#
# We read in and organise the population data. Corrections are based on the naming conventions
# detailed in https://www.citypopulation.de
#
# # There are two outputs: "bra2_population_data.rds" and "bra2_population_data_2000_2014.rds", 
# which are data frames with the yearly population size in each municipality for 2001-2021 and
# the yearly interpolated population size in each municipality for 2000-2014 respectively. 
#
#################################################################################################

#################################################################################################
# Definition of required packages and directories

# required packages
library(tidyverse)
library(readxl)
library(magrittr)

# important directories
POPULATION_DATA_FOLDER <- "bra2_population/"
CASE_DATA_FOLDER <- "bra2_cases/"

#################################################################################################
# Organize the population data

# read in population data
bra2_population_data <- read_excel(path=paste0(POPULATION_DATA_FOLDER, "bra2_population_data.xlsx"), 
                                   sheet="Tabela", na="...")
data_colnames <- as.numeric(bra2_population_data[3, -1])
data_colnames[1:2] <- c("city_id", "region_name")
bra2_population_data <- bra2_population_data[4:(nrow(bra2_population_data)-1), -1]
colnames(bra2_population_data) <- data_colnames

# split region name into state and municipality
bra2_population_data %<>%
  mutate(state_abbr=gsub(".*[(]([^)]+)[)].*", "\\1", region_name)) %>%
  mutate(city=gsub("\\s*\\([^\\)]+\\)\\s*$","", region_name)) %>%
  select(-region_name) %>%
  select(city_id, state_abbr, city, dplyr::everything())

# change state abbreviations to state names
state_abbr_mapping <- read_excel(path=paste0(CASE_DATA_FOLDER, "state_abbr_mapping.xlsx"))
bra2_population_data %<>%
  left_join(state_abbr_mapping, by=c("state_abbr"="code")) %>%
  select(-state_abbr) %>%
  select(city_id, state, city, dplyr::everything())
  
# correct city names according to naming conventions in citypopulation.de
fix_region_name <- function(state, city, city_id) {
  if (state=="Piauí" && city=="Pau D'Arco do Piauí") city <- "Pau d'Arco do Piauí"
  else if (state=="Rio Grande do Norte" && city=="Olho d'Água do Borges") city <- "Olho-d'Água do Borges"
  else if (state=="Sergipe" && city=="Amparo do São Francisco") city <- "Amparo de São Francisco"
  else if (state=="Bahia" && city=="Iuiu") city <- "Iuiú"
  else if (state=="Minas Gerais" && city=="Pingo d'Água") {city<-"Pingo-d'Água"}
  else if (state=="Santa Catarina" && city=="Grão-Pará") {city<-"Grão Pará"}
  else if (state=="Mato Grosso" && city=="Conquista D'Oeste") {city<-"Conquista d'Oeste"}
  else if (state=="Mato Grosso" && city=="Lambari D'Oeste") {city<-"Lambari d'Oeste"}
  return(c(state, city, city_id))
}
corrected_names <- t(apply(select(bra2_population_data, state, city, city_id), 1, 
                           function(x) fix_region_name(state=x[1], city=x[2], city_id=x[3])))
bra2_population_data$state <- corrected_names[, 1]
bra2_population_data$city <- corrected_names[, 2]
bra2_population_data$city_id <- as.numeric(corrected_names[, 3])

# reorganize the dataframe in long format
bra2_population_data %<>% 
  gather(key=year, value=n, `2001`:`2021`) %>%
  mutate(year=as.numeric(year), n=as.numeric(n))

#################################################################################################
# here we perform some checks on non-interpolated population data

# check that all city ids are unique (i.e. there are no duplicates)
count1 <- 0
duplicates <- bra2_population_data %>%
  select(state, city_id, city, year) %>%
  duplicated() %>%
  bra2_population_data[., ]
if (nrow(duplicates)!=0) {
  count1 <- 1
  print("There are duplicate entries in the population data.")
  print(duplicates)
}

# check that population data is available for all regions in the DENV case data
bra_denv_case_data <- readRDS(file.path(CASE_DATA_FOLDER, "bra_denv_case_data.rds"))
count2 <- 0
for (ii in 2000:2016) {
  region_names <- bra_denv_case_data %>%
    filter(year==ii) %>%
    select(state, city) %>%
    unique()
  region_names <- region_names[which(!grepl("Município ignorado", region_names$city, fixed=TRUE)), ]
  no_match_rows <- region_names %>%
    anti_join(unique(select(bra2_population_data, city_id, state, city)), by=c("state", "city"))
  if (nrow(no_match_rows)!=0) {
    count2 <- count2 + 1
    cat("The following regions do not have population data: \n")
    print(no_match_rows)
  }
}
if (count2==0) print("Population data available for all regions in the DENV case data.")

if((count1+count2)==0) {
  print("All tests passed!")
  print("Saving the population data...")
  saveRDS(bra2_population_data, file.path(POPULATION_DATA_FOLDER, "bra2_population_data.rds"))
}

#################################################################################################
# here we linearly interpolate the data for years which do not have population estimates

get_interpolated_pop_data <- function(this_state, this_city, this_city_id) {
  city_data <- bra2_population_data %>%
    filter(city_id==this_city_id, state==this_state, city==this_city) %>% 
    arrange(-desc(year))
  years <- city_data$year
  n <- city_data$n
  
  interpolated_data <- approx(x=years, y=n, xout=2000:2021, rule=2)
  interpolated_years <- interpolated_data$x
  interpolated_n <- interpolated_data$y
  adj_city_data <- data.frame(city_id=this_city_id, state=this_state, 
                              city=this_city, year=interpolated_years, 
                              n=interpolated_n)
  return(adj_city_data)
}
adj_bra2_population_data <- bra2_population_data %>%
  arrange(desc(state), desc(city), desc(-year)) %>% 
  select(state, city, city_id) %>%
  unique() %>% 
  apply(1, function(x) get_interpolated_pop_data(this_state=x[1], 
                                                 this_city=x[2], 
                                                 this_city_id=x[3])) %>%
  do.call(rbind, .)

#################################################################################################
# saving the population data

if((count1+count2)==0) {
  print("All tests passed!")
  print("Saving the interpolated population data...")
  adj_bra2_population_data %<>% select(-city_id)
  bra2_population_data_2000_2014 <- adj_bra2_population_data %>% filter(year<=2014)
  #saveRDS(adj_bra2_population_data, file.path(POPULATION_DATA_FOLDER, "adj_bra2_population_data.rds"))
  saveRDS(bra2_population_data_2000_2014, file.path(POPULATION_DATA_FOLDER, "bra2_population_data_2000_2014.rds"))
}

