#################################################################################################
#
# We organise the polygon data collected from GADM (https://gadm.org/) for Brazil at the 
# municipality level. The municipality names are changed to conform with the DENV case data which 
# we consider the standard. This data set is outdated with administrative boundaries defined when 
# the country had 5504 municipalities (approx. 2004). Now, it has 5570 municipalities with slightly 
# different boundaries. We manually correct the administrative boundaries by adding new municipalities
# and changing the boundaries of existing municipalities. 
# 
# There are two outputs: "simplified_bra2_sf_data.rds" and "bra2_sf_data.rds", which are simple
# feature objects of the administrative boundaries of the 5570 municipalities across Brazil. 
# The simplified version is mostly used in later analyses to decrease computational time. 
#
#################################################################################################

#################################################################################################
# Definition of required packages, directories and local data

# libraries
library(pacman)
p_load("tidyverse", "sf", "osmdata", "rgdal")

# local data:
bra2 <- readRDS("bra2_sf/gadm36_BRA_2_sf.rds") %>%
  transmute(state=NAME_1, city=NAME_2)

#################################################################################################
# let's define some helper functions needed to correct the boundaries of the administrative units

# remove_region removes a region (state, city) from the dataset (data)
remove_region <- function(data, city, state) {
  city_pos <- which(data$state==state & data$city==city)
  data[city_pos, ][, c(1, 2)] <- c(NA, NA)
  return(data)
}

# crop_region crops a region (to_crop_state, to_crop_city) according to the geographical boundary 
# of another region (crop_from_state, crop_from_city)
crop_region <- function(data, crop_from_city, crop_from_state, to_crop_city, to_crop_state) {
  crop_extent <- data[which(data$state==crop_from_state & data$city==crop_from_city), ]$geometry
  cropped_geometry <- st_difference(data[which(data$state==to_crop_state & data$city==to_crop_city), ]$geometry, 
                                    crop_extent)
  if (length(cropped_geometry)>1) {
    for (ii in seq_along(cropped_geometry)) {
      if (class(cropped_geometry[ii])[1] %in% c("sfc_MULTIPOLYGON", "sfc_POLYGON")) {
        cropped_geometry <- cropped_geometry[ii]
        break
      }
    }
  }
  data[which(data$state==to_crop_state & data$city==to_crop_city), ]$geometry <- cropped_geometry
  return(data)
}

# merge_duplicates merges duplicate entries for a given region (state, city)
merge_duplicates <- function(data, city, state) {
  city_pos <- which(data$state==state & data$city==city)
  union_geometry <- st_union(data[city_pos, ]$geometry)
  if (length(union_geometry)>1) union_geometry <- st_combine(data[city_pos, ]$geometry)
  data[city_pos, ][1, ]$geometry <- union_geometry
  data[city_pos, ][2:length(city_pos), c(1,2)] <- c(NA, NA)
  return(data)
}

# split_and_add takes the datset (data), the new city name (new_state, new_city), and the old city name 
# (old_state, old_city) that contains the new city and returns an update spatial data set
split_and_add <- function(data, desc, new_city, old_city=NULL, new_state, old_state=NULL) {
  # "Santa Cruz, Paraíba" is an exception which must be included manually
  if (desc=="Santa Cruz, Paraíba") {
    lon <- c(-39, -37, -37, -39)
    lat <- c(-7, -7, -6, -6)
    xym <- cbind(lon, lat) %>% as.data.frame()
    bound_poly <- xym %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_bbox() %>%
      st_as_sfc()
    new_city_polygon <- st_intersection(data[which(data$state==old_state & data$city==old_city), ]$geometry, 
                                        bound_poly) 
    new_city_polygon <- new_city_polygon %>%
      st_sf %>%
      st_cast("MULTIPOLYGON") %>%
      transmute(state=new_state, city=new_city)
    
    data[which(data$state==old_state & data$city==old_city), ]$geometry <- 
      st_difference(data[which(data$state==old_state & data$city==old_city), ]$geometry, 
                    new_city_polygon$geometry)[2] %>%
      st_cast("MULTIPOLYGON")
    data <- rbind(data, new_city_polygon)
    return(data)
  }
  
  # split the old city into the new city and the remainder of the old city
  new_city_polygon <- opq(bbox=getbb(desc)) %>%
    add_osm_feature(key = "boundary", value = "administrative") %>%
    add_osm_feature(key="admin_level", value=8) %>%
    osmdata_sf()
  new_city_polygon <- new_city_polygon$osm_multipolygons
  new_city_polygon <- new_city_polygon[which(new_city_polygon$name==new_city), ] %>%
    transmute(state=new_state, city=name)
  if (!is.null(old_city)) {
    new_geometry <- st_difference(data[which(data$state==old_state & data$city==old_city), ]$geometry, 
                                  new_city_polygon$geometry)
    if (length(new_geometry)>1) {
      for (ii in seq_along(new_geometry)) {
        if (class(new_geometry[ii])[1] %in% c("sfc_MULTIPOLYGON", "sfc_POLYGON")) {
          new_geometry <- new_geometry[ii]
          break
        }
      }
    }
    data[which(data$state==old_state & data$city==old_city), ]$geometry <- new_geometry
  }
  data <- rbind(data, new_city_polygon)
  return(data)
}

#################################################################################################
# we start by manually correcting the polygons of municipalities whose administrative boundaries 
# were outdated or incorrectly defined. 

bra2 <- split_and_add(bra2, desc="Governador Lindenberg", new_city="Governador Lindenberg", old_city="Colatina", new_state="Espírito Santo", old_state="Espírito Santo")
bra2 <- split_and_add(bra2, desc="Pau d'Arco do Piauí", new_city="Pau d'Arco do Piauí", old_city="Altos", new_state="Piauí", old_state="Piauí")
bra2 <- split_and_add(bra2, desc="Itapiranga, Amazonas", new_city="Itapiranga", new_state="Amazonas", old_state="Amazonas") %>%
  split_and_add(desc="Manaus, Amazonas", new_city="Manaus", old_city="Maués", new_state="Amazonas", old_state="Amazonas")
bra2 <- split_and_add(bra2, desc="Dois Irmãos do Tocantins, Tocantins", new_city="Dois Irmãos do Tocantins", old_city="Divinópolis do Tocantins", new_state="Tocantins", old_state="Tocantins") %>%
  split_and_add(desc="Guaraí, Tocantins", new_city="Guaraí", old_city="Tupirama", new_state="Tocantins", old_state="Tocantins")
bra2 <- split_and_add(bra2, desc="Aroeiras do Itaim, Piauí", new_city="Aroeiras do Itaim", old_city="Picos", new_state="Piauí", old_state="Piauí") %>%
  split_and_add(desc="Nazária, Piauí", new_city="Nazária", old_city="Teresina", new_state="Piauí", old_state="Piauí")
bra2 <- split_and_add(bra2, desc="Várzea, Paraíba", new_city="Várzea", new_state="Paraíba", old_city="Várzea", old_state="Rio Grande do Norte") %>%
  split_and_add(desc="Jundiá, Rio Grande do Norte", new_city="Jundiá", new_state="Rio Grande do Norte", old_city="Várzea", old_state="Rio Grande do Norte") %>%
  split_and_add(desc="Santa Cruz, Paraíba", new_city="Santa Cruz", new_state="Paraíba", old_city="Santa Cruz", old_state="Rio Grande do Norte")
bra2 <- split_and_add(bra2, desc="Orocó, Pernambuco", new_city="Orocó", new_state="Pernambuco", old_city="Orocó", old_state="Pernambuco") %>%
  split_and_add(desc="São José do Egito, Pernambuco", new_city="São José do Egito", new_state="Pernambuco", old_city="São João", old_state="Pernambuco")
bra2 <- split_and_add(bra2, desc="Jequiá da Praia, Alagoas", new_city="Jequiá da Praia", new_state="Alagoas", old_city="São Miguel dos Campos", old_state="Alagoas")
bra2 <- split_and_add(bra2, desc="Barrocas, Bahia", new_city="Barrocas", new_state="Bahia", old_city="Serrinha", old_state="Bahia") %>%
  split_and_add(desc="Luís Eduardo Magalhães, Bahia", new_city="Luís Eduardo Magalhães", new_state="Bahia", old_city="Barreiras", old_state="Bahia") %>%
  split_and_add(desc="Wanderley, Bahia", new_city="Wanderley", new_state="Bahia", old_city="Cotegipe", old_state="Bahia")
bra2 <- split_and_add(bra2, desc="Itaipé, Minas Gerais", new_city="Itaipé", new_state="Minas Gerais", old_city="Itajubá", old_state="Minas Gerais") %>%
  split_and_add(desc="Itaguara, Minas Gerais", new_city="Itaguara", new_state="Minas Gerais", old_city="Itaipé", old_state="Minas Gerais") %>%
  split_and_add(desc="Itacambira, Minas Gerais", new_city="Itacambira", new_state="Minas Gerais", old_city="Itaguara", old_state="Minas Gerais")
bra2 <- split_and_add(bra2, desc="Rio Piracicaba, Minas Gerais", new_city="Rio Piracicaba", new_state="Minas Gerais", old_city="Itamogi", old_state="Minas Gerais") %>%
  split_and_add(desc="Itamogi, Minas Gerais", new_city="Itamogi", new_state="Minas Gerais", old_city="Itamonte", old_state="Minas Gerais") %>%
  split_and_add(desc="Itamonte, Minas Gerais", new_city="Itamonte", new_state="Minas Gerais", old_city="Itanhandu", old_state="Minas Gerais") %>%
  split_and_add(desc="Itanhandu, Minas Gerais", new_city="Itanhandu", new_state="Minas Gerais", old_city="Itanhomi", old_state="Minas Gerais")
bra2 <- split_and_add(bra2, desc="Vazante, Minas Gerais", new_city="Vazante", new_state="Minas Gerais", old_city="Lagoa Grande", old_state="Minas Gerais")
bra2 <- split_and_add(bra2, desc="Itabirito, Minas Gerais", new_city="Itabirito", new_state="Minas Gerais", old_city="Itacarambi", old_state="Minas Gerais") %>%
  merge_duplicates(city="Itabirito", state="Minas Gerais") %>%
  split_and_add(desc="Itabira, Minas Gerais", new_city="Itabira", new_state="Minas Gerais", old_city="Itabirito", old_state="Minas Gerais") %>%
  split_and_add(desc="Santa Vitória, Minas Gerais", new_city="Santa Vitória", new_state="Minas Gerais", old_city="Itabira", old_state="Minas Gerais") %>%
  merge_duplicates(city="Santa Vitória", state="Minas Gerais")
bra2 <- split_and_add(bra2, desc="Itatiaiuçu, Minas Gerais", new_city="Itatiaiuçu", new_state="Minas Gerais", old_city="Itaú de Minas", old_state="Minas Gerais") %>%
  merge_duplicates(city="Itatiaiuçu", state="Minas Gerais") %>%
  remove_region(city="Itaú de Minas", state="Minas Gerais") %>%
  split_and_add(desc="Itaú de Minas, Minas Gerais", new_city="Itaú de Minas", new_state="Minas Gerais", old_city="Itabira", old_state="Minas Gerais") %>%
  merge_duplicates(city="Itabira", state="Minas Gerais")

bra2 <- split_and_add(bra2, desc="Itabirinha, Minas Gerais", new_city="Itabirinha", new_state="Minas Gerais", old_city="Itabirito", old_state="Minas Gerais") %>%
  split_and_add(desc="Itabirinha, Minas Gerais", new_city="Itabirinha", new_state="Minas Gerais") %>%
  merge_duplicates(city="Itabirinha", state="Minas Gerais")
bra2 <- split_and_add(bra2, desc="Mesquita, Rio de Janeiro", new_city="Mesquita", new_state="Rio de Janeiro", old_city="Nova Iguaçu", old_state="Rio de Janeiro")
bra2 <- split_and_add(bra2, desc="Ilha Comprida, São Paulo", new_city="Ilha Comprida", new_state="São Paulo", old_city="Cananéia", old_state="São Paulo")
bra2 <- split_and_add(bra2, desc="Balneário Rincão, Santa Catarina", new_city="Balneário Rincão", new_state="Santa Catarina", old_city="Içara", old_state="Santa Catarina")
bra2 <- split_and_add(bra2, desc="Pescaria Brava, Santa Catarina", new_city="Pescaria Brava", new_state="Santa Catarina", old_city="Laguna", old_state="Santa Catarina")
bra2 <- split_and_add(bra2, desc="Lagoa Santa, Goiás", new_city="Lagoa Santa", new_state="Goiás", old_city="Itajá", old_state="Goiás")
bra2 <- split_and_add(bra2, desc="Ipiranga de Goiás, Goiás", new_city="Ipiranga de Goiás", new_state="Goiás", old_city="Ceres", old_state="Goiás")
bra2 <- split_and_add(bra2, desc="Gameleira de Goiás, Goiás", new_city="Gameleira de Goiás", new_state="Goiás", old_city="Silvânia", old_state="Goiás")
bra2 <- split_and_add(bra2, desc="Damolândia, Goiás", new_city="Damolândia", new_state="Goiás", old_city="Inhumas", old_state="Goiás")
bra2 <- split_and_add(bra2, desc="Campo Limpo de Goiás, Goiás", new_city="Campo Limpo de Goiás", new_state="Goiás", old_city="Anápolis", old_state="Goiás")
bra2 <- split_and_add(bra2, desc="Britânia, Goiás", new_city="Britânia", new_state="Goiás", old_city="Aruanã", old_state="Goiás")
bra2 <- split_and_add(bra2, desc="Vale de São Domingos, Mato Grosso", new_city="Vale de São Domingos", new_state="Mato Grosso", old_city="Pontes e Lacerda", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Sorriso, Mato Grosso", new_city="Sorriso", new_state="Mato Grosso", old_city="Tapurah", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Santo Antônio do Leste, Mato Grosso", new_city="Santo Antônio do Leste", new_state="Mato Grosso", old_city="Novo São Joaquim", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Santa Rita do Trivelato, Mato Grosso", new_city="Santa Rita do Trivelato", new_state="Mato Grosso", old_city="Nova Mutum", old_state="Mato Grosso") %>%
  split_and_add(desc="Santa Rita do Trivelato, Mato Grosso", new_city="Santa Rita do Trivelato", new_state="Mato Grosso", old_city="Rosário Oeste", old_state="Mato Grosso") %>%
  merge_duplicates(city="Santa Rita do Trivelato", state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Santa Cruz do Xingu, Mato Grosso", new_city="Santa Cruz do Xingu", new_state="Mato Grosso", old_city="São José do Xingu", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Rondolândia, Mato Grosso", new_city="Rondolândia", new_state="Mato Grosso", old_city="Aripuanã", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Bom Jesus do Araguaia, Mato Grosso", new_city="Bom Jesus do Araguaia", new_state="Mato Grosso", old_city="Alto Boa Vista", old_state="Mato Grosso") %>%
  split_and_add(desc="Bom Jesus do Araguaia, Mato Grosso", new_city="Bom Jesus do Araguaia", new_state="Mato Grosso", old_city="Ribeirão Cascalheira", old_state="Mato Grosso") %>%
  merge_duplicates(city="Bom Jesus do Araguaia", state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Novo Santo Antônio, Mato Grosso", new_city="Novo Santo Antônio", new_state="Mato Grosso", old_city="Cocalinho", old_state="Mato Grosso") %>%
  split_and_add(desc="Novo Santo Antônio, Mato Grosso", new_city="Novo Santo Antônio", new_state="Mato Grosso", old_city="São Félix Xingu", old_state="Mato Grosso") %>%
  merge_duplicates(city="Novo Santo Antônio", state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Serra Nova Dourada, Mato Grosso", new_city="Serra Nova Dourada", new_state="Mato Grosso", old_city="Alto Boa Vista", old_state="Mato Grosso") %>%
  split_and_add(desc="Serra Nova Dourada, Mato Grosso", new_city="Serra Nova Dourada", new_state="Mato Grosso", old_city="São Félix Xingu", old_state="Mato Grosso") %>%
  merge_duplicates(city="Serra Nova Dourada", state="Mato Grosso") %>%
  split_and_add(desc="Serra Nova Dourada, Mato Grosso", new_city="Serra Nova Dourada", new_state="Mato Grosso", old_city="Bom Jesus do Araguaia", old_state="Mato Grosso") %>%
  merge_duplicates(city="Serra Nova Dourada", state="Mato Grosso")

bra2 <- split_and_add(bra2, desc="Nova Santa Helena, Mato Grosso", new_city="Nova Santa Helena", new_state="Mato Grosso", old_city="Itaúba", old_state="Mato Grosso") %>%
  crop_region(crop_from_city="Nova Santa Helena", crop_from_state="Mato Grosso", to_crop_city="Cláudia", to_crop_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Nova Nazaré, Mato Grosso", new_city="Nova Nazaré", new_state="Mato Grosso", old_city="Água Boa", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Itanhangá, Mato Grosso", new_city="Itanhangá", new_state="Mato Grosso", old_city="Tapurah", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Ipiranga do Norte, Mato Grosso", new_city="Ipiranga do Norte", new_state="Mato Grosso", old_city="Tapurah", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Colniza, Mato Grosso", new_city="Colniza", new_state="Mato Grosso", old_city="Aripuanã", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Conquista d'Oeste, Mato Grosso", new_city="Conquista d'Oeste", new_state="Mato Grosso", old_city="Pontes e Lacerda", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Lambari d'Oeste, Mato Grosso", new_city="Lambari d'Oeste", new_state="Mato Grosso", old_city="Cáceres", old_state="Mato Grosso") %>%
  merge_duplicates(city="Lambari d'Oeste", state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Curvelândia, Mato Grosso", new_city="Curvelândia", new_state="Mato Grosso", old_city="Cáceres", old_state="Mato Grosso")
bra2 <- split_and_add(bra2, desc="Ladário, Mato Grosso do Sul", new_city="Ladário", new_state="Mato Grosso do Sul", old_city="Corumbá", old_state="Mato Grosso do Sul")
bra2 <- split_and_add(bra2, desc="Figueirão, Mato Grosso do Sul", new_city="Figueirão", new_state="Mato Grosso do Sul", old_city="Camapuã", old_state="Mato Grosso do Sul")
bra2 <- split_and_add(bra2, desc="Westfália, Rio Grande do Sul", new_city="Westfália", new_state="Rio Grande do Sul", old_city="Imigrante", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Westfália", crop_from_state="Rio Grande do Sul", to_crop_city="Teutônia", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Tio Hugo, Rio Grande do Sul", new_city="Tio Hugo", new_state="Rio Grande do Sul", old_city="Ibirapuitã", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Tio Hugo", crop_from_state="Rio Grande do Sul", to_crop_city="Ernestina", to_crop_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Tio Hugo", crop_from_state="Rio Grande do Sul", to_crop_city="Victor Graeff", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="São Vicente do Sul, Rio Grande do Sul", new_city="São Vicente do Sul", new_state="Rio Grande do Sul", old_city="São Valentim do Sul", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="São Vicente do Sul", crop_from_state="Rio Grande do Sul", to_crop_city="São Valentim do Sul", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="São Pedro das Missões, Rio Grande do Sul", new_city="São Pedro das Missões", new_state="Rio Grande do Sul", old_city="Palmeira das Missões", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="São José do Sul, Rio Grande do Sul", new_city="São José do Sul", new_state="Rio Grande do Sul", old_city="Salvador do Sul", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="São José do Sul", crop_from_state="Rio Grande do Sul", to_crop_city="Montenegro", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="São José do Norte, Rio Grande do Sul", new_city="São José do Norte", new_state="Rio Grande do Sul", old_city="Porto Alegre", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Santa Margarida do Sul, Rio Grande do Sul", new_city="Santa Margarida do Sul", new_state="Rio Grande do Sul", old_city="São Gabriel", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Santa Cecília do Sul, Rio Grande do Sul", new_city="Santa Cecília do Sul", new_state="Rio Grande do Sul", old_city="Água Santa", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Santa Cecília do Sul", crop_from_state="Rio Grande do Sul", to_crop_city="Tapejara", to_crop_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Santa Cecília do Sul", crop_from_state="Rio Grande do Sul", to_crop_city="Ibiaçá", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Rolador, Rio Grande do Sul", new_city="Rolador", new_state="Rio Grande do Sul", old_city="São Luiz Gonzaga", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Quatro Irmãos, Rio Grande do Sul", new_city="Quatro Irmãos", new_state="Rio Grande do Sul", old_city="Erechim", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Quatro Irmãos", crop_from_state="Rio Grande do Sul", to_crop_city="Jacutinga", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Paulo Bento, Rio Grande do Sul", new_city="Paulo Bento", new_state="Rio Grande do Sul", old_city="Erechim", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Paulo Bento", crop_from_state="Rio Grande do Sul", to_crop_city="Barão de Cotegipe", to_crop_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Paulo Bento", crop_from_state="Rio Grande do Sul", to_crop_city="Ponte Preta", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Pinto Bandeira, Rio Grande do Sul", new_city="Pinto Bandeira", new_state="Rio Grande do Sul", old_city="Bento Gonçalves", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Pinhal da Serra, Rio Grande do Sul", new_city="Pinhal da Serra", new_state="Rio Grande do Sul", old_city="Esmeralda", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Pedras Altas, Rio Grande do Sul", new_city="Pedras Altas", new_state="Rio Grande do Sul", old_city="Herval", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Novo Xingu, Rio Grande do Sul", new_city="Novo Xingu", new_state="Rio Grande do Sul", old_city="Constantina", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Mato Queimado, Rio Grande do Sul", new_city="Mato Queimado", new_state="Rio Grande do Sul", old_city="Caibaté", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Lagoa Bonita do Sul, Rio Grande do Sul", new_city="Lagoa Bonita do Sul", new_state="Rio Grande do Sul", old_city="Sobradinho", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Jacuizinho, Rio Grande do Sul", new_city="Jacuizinho", new_state="Rio Grande do Sul", old_city="Espumoso", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Jacuizinho", crop_from_state="Rio Grande do Sul", to_crop_city="Campos Borges", to_crop_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Jacuizinho", crop_from_state="Rio Grande do Sul", to_crop_city="Salto do Jacuí", to_crop_state="Rio Grande do Sul")

bra2 <- split_and_add(bra2, desc="Itati, Rio Grande do Sul", new_city="Itati", new_state="Rio Grande do Sul", old_city="Terra de Areia", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Lajeado, Rio Grande do Sul", new_city="Lajeado", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, state="Rio Grande do Sul", city="Lajeado do Bugre")
bra2 <- split_and_add(bra2, desc="Lajeado do Bugre, Rio Grande do Sul", new_city="Lajeado do Bugre", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, state="Rio Grande do Sul", city="Lajedão")
bra2 <- split_and_add(bra2, desc="Forquetinha, Rio Grande do Sul", new_city="Forquetinha", new_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Cruzaltense, Rio Grande do Sul", new_city="Cruzaltense", new_state="Rio Grande do Sul", old_city="Campinas do Sul", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Coronel Pilar, Rio Grande do Sul", new_city="Coronel Pilar", new_state="Rio Grande do Sul", old_city="Garibaldi", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Coqueiro Baixo, Rio Grande do Sul", new_city="Coqueiro Baixo", new_state="Rio Grande do Sul", old_city="Nova Bréscia", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Capão do Cipó, Rio Grande do Sul", new_city="Capão do Cipó", new_state="Rio Grande do Sul", old_city="Santiago", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Capão do Cipó", crop_from_state="Rio Grande do Sul", to_crop_city="São Miguel das Missões", to_crop_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Capão do Cipó", crop_from_state="Rio Grande do Sul", to_crop_city="Bossoroca", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Capão Bonito do Sul, Rio Grande do Sul", new_city="Capão Bonito do Sul", new_state="Rio Grande do Sul", old_city="Lagoa Vermelha", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Canudos do Vale, Rio Grande do Sul", new_city="Canudos do Vale", new_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Bozano, Rio Grande do Sul", new_city="Bozano", new_state="Rio Grande do Sul", old_city="Ijuí", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Boa Vista do Incra, Rio Grande do Sul", new_city="Boa Vista do Incra", new_state="Rio Grande do Sul", old_city="Cruz Alta", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Boa Vista do Incra", crop_from_state="Rio Grande do Sul", to_crop_city="Fortaleza dos Valos", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Boa Vista do Cadeado, Rio Grande do Sul", new_city="Boa Vista do Cadeado", new_state="Rio Grande do Sul", old_city="Cruz Alta", old_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Boa Vista do Cadeado", crop_from_state="Rio Grande do Sul", to_crop_city="Augusto Pestana", to_crop_state="Rio Grande do Sul") %>%
  crop_region(crop_from_city="Boa Vista do Cadeado", crop_from_state="Rio Grande do Sul", to_crop_city="Ijuí", to_crop_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Arroio do Padre, Rio Grande do Sul", new_city="Arroio do Padre", new_state="Rio Grande do Sul", old_city="Pelotas", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Almirante Tamandaré do Sul, Rio Grande do Sul", new_city="Almirante Tamandaré do Sul", new_state="Rio Grande do Sul", old_city="Carazinho", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Aceguá, Rio Grande do Sul", new_city="Aceguá", new_state="Rio Grande do Sul", old_city="Bagé", old_state="Rio Grande do Sul")
bra2 <- split_and_add(bra2, desc="Três Barras, Santa Catarina", new_city="Três Barras", new_state="Santa Catarina", old_city="São Mateus do Sul", old_state="Paraná")
bra2 <- split_and_add(bra2, desc="Lapa, Paraná", new_city="Lapa", new_state="Paraná", old_city="Balsa Nova", old_state="Paraná")
bra2 <- split_and_add(bra2, desc="Santo Antônio do Jardim, São Paulo", new_city="Santo Antônio do Jardim", new_state="São Paulo", old_city="Balsa Nova", old_state="São Paulo")
bra2 <- split_and_add(bra2, desc="Piraju, São Paulo", new_city="Piraju", new_state="São Paulo", old_city="Pirajui", old_state="São Paulo")

bra2 <- remove_region(bra2, city="Aratuípe", state="Bahia") %>%
  remove_region(city="Aratuipe", state="Bahia") %>%
  split_and_add(desc="Aratuípe, Bahia", new_city="Aratuípe", new_state="Bahia", old_city="Andarai", old_state="Bahia")
bra2 <- split_and_add(bra2, desc="Arco-Íris, São Paulo", new_city="Arco-Íris", new_state="São Paulo") %>%
  remove_region(city="Arco-íris", state="São Paulo") 
bra2 <- split_and_add(bra2, desc="Altamira do Paraná, Paraná", new_city="Altamira do Paraná", new_state="Paraná") 
bra2 <- remove_region(bra2, city="Melgaço", state="Pará") %>%
  remove_region(city="Melgaco", state="Pará") %>%
  split_and_add(desc="Melgaço, Pará", new_city="Melgaço", new_state="Pará") 
bra2 <- remove_region(bra2, city="Andirá", state="Paraná") %>%
  remove_region(city="Andira", state="Paraná") %>%
  split_and_add(desc="Andirá, Paraná", new_city="Andirá", new_state="Paraná") 
bra2 <- remove_region(bra2, city="Tutóia", state="Maranhão") %>%
  remove_region(city="Tutoia", state="Maranhão") %>%
  split_and_add(desc="Tutóia, Maranhão", new_city="Tutóia", new_state="Maranhão") 
bra2 <- remove_region(bra2, city="Iporã", state="Paraná") %>%
  remove_region(city="Iporá", state="Paraná") %>%
  split_and_add(desc="Iporã, Paraná", new_city="Iporã", new_state="Paraná") 
bra2 <- remove_region(bra2, city="Boa Vista das Missões", state="Rio Grande do Sul") %>%
  remove_region(city="Boa Vista das Misses", state="Rio Grande do Sul") %>%
  split_and_add(desc="Boa Vista das Missões, Rio Grande do Sul", new_city="Boa Vista das Missões", new_state="Rio Grande do Sul") 
bra2 <- merge_duplicates(bra2, city="Aceguá", state="Rio Grande do Sul") 
bra2 <- merge_duplicates(bra2, city="Itaguara", state="Minas Gerais") 
bra2 <- remove_region(bra2, city="Itaipé", state="Minas Gerais") %>%
  remove_region(city="Itaguara", state="Minas Gerais") %>%
  remove_region(city="Itatiaiuçu", state="Minas Gerais") %>%
  remove_region(city="Itapeva", state="Minas Gerais") %>%
  split_and_add(desc="Itaipé, Minas Gerais", new_city="Itaipé", new_state="Minas Gerais") %>%
  split_and_add(desc="Itaguara, Minas Gerais", new_city="Itaguara", new_state="Minas Gerais") %>%
  split_and_add(desc="Itatiaiuçu, Minas Gerais", new_city="Itatiaiuçu", new_state="Minas Gerais") %>%
  split_and_add(desc="Itapeva, Minas Gerais", new_city="Itapeva", new_state="Minas Gerais") 
bra2 <- remove_region(bra2, city="Itamogi", state="Minas Gerais") %>%
  remove_region(city="Itambé do Mato Dentro", state="Minas Gerais") %>%
  remove_region(city="Itabira", state="Minas Gerais") %>%
  split_and_add(desc="Itamogi, Minas Gerais", new_city="Itamogi", new_state="Minas Gerais") %>%
  split_and_add(desc="Itambé do Mato Dentro, Minas Gerais", new_city="Itambé do Mato Dentro", new_state="Minas Gerais") %>%
  split_and_add(desc="Itabira, Minas Gerais", new_city="Itabira", new_state="Minas Gerais") 
bra2 <- remove_region(bra2, city="Itamonte", state="Minas Gerais") %>%
  remove_region(city="Itanhandu", state="Minas Gerais") %>%
  split_and_add(desc="Itamonte, Minas Gerais", new_city="Itamonte", new_state="Minas Gerais") %>%
  split_and_add(desc="Itanhandu, Minas Gerais", new_city="Itanhandu", new_state="Minas Gerais") 
bra2 <- remove_region(bra2, city="Pirajuí", state="São Paulo") %>% 
  remove_region(city="Pirajui", state="São Paulo") %>%
  split_and_add(desc="Pirajuí, São Paulo", new_city="Pirajuí", new_state="São Paulo") 
bra2 <- remove_region(bra2, city="Pinhal", state="São Paulo") 
bra2 <- remove_region(bra2, city="Conceição do Almeida", state="Bahia") %>%
  remove_region(city="Conceicao do Almeida", state="Bahia") %>% 
  split_and_add(desc="Conceição do Almeida, Bahia", new_city="Conceição do Almeida", new_state="Bahia") 
bra2 <- remove_region(bra2, city="Itacarambi", state="Minas Gerais") %>%
  remove_region(city="Itacarambira", state="Minas Gerais") %>%
  split_and_add(desc="Itacarambi, Minas Gerais", new_city="Itacarambi", new_state="Minas Gerais") %>%
  split_and_add(desc="Itacambira, Minas Gerais", new_city="Itacambira", new_state="Minas Gerais") 
bra2 <- remove_region(bra2, city="Altamira do Paran", state="Paraná") 
bra2 <- remove_region(bra2, city="Serra Nova Dourada", state="Mato Grosso") %>%
  split_and_add(desc="Serra Nova Dourada, Mato Grosso", new_city="Serra Nova Dourada", new_state="Mato Grosso") 
bra2 <- remove_region(bra2, city="Jiquiriçá", state="Bahia") %>%
  remove_region(city="Jeremoabo", state="Bahia") %>%
  split_and_add(desc="Jiquiriçá, Bahia", new_city="Jiquiriçá", new_state="Bahia") %>% 
  split_and_add(desc="Jeremoabo, Bahia", new_city="Jeremoabo", new_state="Bahia") 
bra2 <- remove_region(bra2, city="Brasópolis", state="Minas Gerais") %>%
  remove_region(city="Brasilândia de Minas", state="Minas Gerais") %>%
  remove_region(city="Braúnas", state="Minas Gerais") %>%
  split_and_add(desc="Brazópolis, Minas Gerais", new_city="Brazópolis", new_state="Minas Gerais") %>% 
  split_and_add(desc="Brasilândia de Minas, Minas Gerais", new_city="Brasilândia de Minas", new_state="Minas Gerais")  %>%
  split_and_add(desc="Braúnas, Minas Gerais", new_city="Braúnas", new_state="Minas Gerais") 

bra2 <- remove_region(bra2, city="Capitão Andrade", state="Minas Gerais") %>%
  remove_region(city="Capitão Enéas", state="Minas Gerais") %>%
  split_and_add(desc="Capitão Andrade, Minas Gerais", new_city="Capitão Andrade", new_state="Minas Gerais") %>% 
  split_and_add(desc="Capitão Enéas, Minas Gerais", new_city="Capitão Enéas", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Conceição das Pedras", state="Minas Gerais") %>%
  remove_region(city="Conceição das Alagoas", state="Minas Gerais") %>%
  split_and_add(desc="Conceição das Pedras, Minas Gerais", new_city="Conceição das Pedras", new_state="Minas Gerais") %>% 
  split_and_add(desc="Conceição das Alagoas, Minas Gerais", new_city="Conceição das Alagoas", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Coração de Jesus", state="Minas Gerais") %>%
  remove_region(city="Corinto", state="Minas Gerais") %>%
  remove_region(city="Cordisburgo", state="Minas Gerais") %>%
  remove_region(city="Coroaci", state="Minas Gerais") %>% 
  remove_region(city="Conceição das Alagoas", state="Minas Gerais") %>%
  remove_region(city="Cordislândia", state="Minas Gerais") %>%
  split_and_add(desc="Coração de Jesus, Minas Gerais", new_city="Coração de Jesus", new_state="Minas Gerais") %>% 
  split_and_add(desc="Corinto, Minas Gerais", new_city="Corinto", new_state="Minas Gerais") %>%
  split_and_add(desc="Cordisburgo, Minas Gerais", new_city="Cordisburgo", new_state="Minas Gerais") %>%
  split_and_add(desc="Coroaci, Minas Gerais", new_city="Coroaci", new_state="Minas Gerais") %>%
  split_and_add(desc="Cordislândia, Minas Gerais", new_city="Cordislândia", new_state="Minas Gerais") %>%
  split_and_add(desc="Conceição das Alagoas, Minas Gerais", new_city="Conceição das Alagoas", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Francisco Sá", state="Minas Gerais") %>%
  remove_region(city="Francisco Dumon", state="Minas Gerais") %>%
  remove_region(city="Franciscópolis", state="Minas Gerais") %>%
  remove_region(city="Francisco Badaró", state="Minas Gerais") %>% 
  split_and_add(desc="Francisco Sá, Minas Gerais", new_city="Francisco Sá", new_state="Minas Gerais") %>% 
  split_and_add(desc="Francisco Dumont, Minas Gerais", new_city="Francisco Dumont", new_state="Minas Gerais") %>%
  split_and_add(desc="Franciscópolis, Minas Gerais", new_city="Franciscópolis", new_state="Minas Gerais") %>%
  split_and_add(desc="Francisco Badaró, Minas Gerais", new_city="Francisco Badaró", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Itambacuri", state="Minas Gerais") %>%
  remove_region(city="Itamarati de Minas", state="Minas Gerais") %>%
  remove_region(city="Itamarandiba", state="Minas Gerais") %>%
  split_and_add(desc="Itambacuri, Minas Gerais", new_city="Itambacuri", new_state="Minas Gerais") %>% 
  split_and_add(desc="Itamarati de Minas, Minas Gerais", new_city="Itamarati de Minas", new_state="Minas Gerais") %>%
  split_and_add(desc="Itamarandiba, Minas Gerais", new_city="Itamarandiba", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Itapecerica", state="Minas Gerais") %>%
  remove_region(city="Itapagipe", state="Minas Gerais") %>%
  remove_region(city="Itapeva", state="Minas Gerais") %>%
  remove_region(city="Itaobim", state="Minas Gerais") %>% 
  remove_region(city="Itanhomi", state="Minas Gerais") %>%
  split_and_add(desc="Itapecerica, Minas Gerais", new_city="Itapecerica", new_state="Minas Gerais") %>% 
  split_and_add(desc="Itapagipe, Minas Gerais", new_city="Itapagipe", new_state="Minas Gerais") %>%
  split_and_add(desc="Itapeva, Minas Gerais", new_city="Itapeva", new_state="Minas Gerais") %>%
  split_and_add(desc="Itaobim, Minas Gerais", new_city="Itaobim", new_state="Minas Gerais") %>%
  split_and_add(desc="Itanhomi, Minas Gerais", new_city="Itanhomi", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Palma", state="Minas Gerais") %>%
  remove_region(city="Palmópolis", state="Minas Gerais") %>%
  split_and_add(desc="Palma, Minas Gerais", new_city="Palma", new_state="Minas Gerais") %>% 
  split_and_add(desc="Palmópolis, Minas Gerais", new_city="Palmópolis", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Santa Bárbara", state="Minas Gerais") %>%
  remove_region(city="Santa Bárbara do Monte Verde", state="Minas Gerais") %>%
  remove_region(city="Santa Bárbara do Leste", state="Minas Gerais") %>%
  remove_region(city="Santa Bárbara do Tugúrio", state="Minas Gerais") %>% 
  split_and_add(desc="Santa Bárbara, Minas Gerais", new_city="Santa Bárbara", new_state="Minas Gerais") %>% 
  split_and_add(desc="Santa Bárbara do Monte Verde, Minas Gerais", new_city="Santa Bárbara do Monte Verde", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Bárbara do Leste, Minas Gerais", new_city="Santa Bárbara do Leste", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Bárbara do Tugúrio, Minas Gerais", new_city="Santa Bárbara do Tugúrio", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Santa Rita de Caldas", state="Minas Gerais") %>%
  remove_region(city="Santa Rita do Ibitipoca", state="Minas Gerais") %>%
  remove_region(city="Santa Rita de Jacutinga", state="Minas Gerais") %>%
  remove_region(city="Santa Rita de Minas", state="Minas Gerais") %>% 
  remove_region(city="Santa Rita Itueto", state="Minas Gerais") %>%
  remove_region(city="Santa Rita do Sapucaí", state="Minas Gerais") %>%
  remove_region(city="Santa Rosa da Serra", state="Minas Gerais") %>%
  split_and_add(desc="Santa Rita de Caldas, Minas Gerais", new_city="Santa Rita de Caldas", new_state="Minas Gerais") %>% 
  split_and_add(desc="Santa Rita de Ibitipoca, Minas Gerais", new_city="Santa Rita de Ibitipoca", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Rita de Jacutinga, Minas Gerais", new_city="Santa Rita de Jacutinga", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Rita de Minas, Minas Gerais", new_city="Santa Rita de Minas", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Rita do Itueto, Minas Gerais", new_city="Santa Rita do Itueto", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Rita do Sapucaí, Minas Gerais", new_city="Santa Rita do Sapucaí", new_state="Minas Gerais") %>%
  split_and_add(desc="Santa Rosa da Serra, Minas Gerais", new_city="Santa Rosa da Serra", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Santo Antônio do Amparo", state="Minas Gerais") %>%
  remove_region(city="Santo Antônio do Aventureiro", state="Minas Gerais") %>%
  remove_region(city="Santo Antônio do Grama", state="Minas Gerais") %>%
  remove_region(city="San Antonio do Itambe", state="Minas Gerais") %>% 
  remove_region(city="Santo Antônio do Jacinto", state="Minas Gerais") %>%
  remove_region(city="Santo Antônio do Monte", state="Minas Gerais") %>%
  remove_region(city="Santo Antônio do Retiro", state="Minas Gerais") %>%
  remove_region(city="San Antonio do Rio Abai", state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Amparo, Minas Gerais", new_city="Santo Antônio do Amparo", new_state="Minas Gerais") %>% 
  split_and_add(desc="Santo Antônio do Aventureiro, Minas Gerais", new_city="Santo Antônio do Aventureiro", new_state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Grama, Minas Gerais", new_city="Santo Antônio do Grama", new_state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Itambé, Minas Gerais", new_city="Santo Antônio do Itambé", new_state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Jacinto, Minas Gerais", new_city="Santo Antônio do Jacinto", new_state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Monte, Minas Gerais", new_city="Santo Antônio do Monte", new_state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Retiro, Minas Gerais", new_city="Santo Antônio do Retiro", new_state="Minas Gerais") %>%
  split_and_add(desc="Santo Antônio do Rio Abaixo, Minas Gerais", new_city="Santo Antônio do Rio Abaixo", new_state="Minas Gerais")
bra2 <- split_and_add(bra2, desc="Mojuí dos Campos, Pará", new_city="Mojuí dos Campos", new_state="Pará", old_city="Santarém", old_state="Pará")
bra2 <- remove_region(bra2, city="Poço de José de Moura", state="Paraíba") %>%
  remove_region(city="Pombal", state="Paraíba") %>%
  split_and_add(desc="Poço de José de Moura, Paraíba", new_city="Poço de José de Moura", new_state="Paraíba") %>% 
  split_and_add(desc="Pombal, Paraíba", new_city="Pombal", new_state="Paraíba")

bra2 <- remove_region(bra2, city="Ivaté", state="Paraná") %>%
  remove_region(city="Ivatuva", state="Paraná") %>%
  split_and_add(desc="Ivaté, Paraná", new_city="Ivaté", new_state="Paraná") %>% 
  split_and_add(desc="Ivatuba, Paraná", new_city="Ivatuba", new_state="Paraná")
bra2 <- remove_region(bra2, city="Santa Maria do Oeste", state="Paraná") %>%
  remove_region(city="Santa Mônica", state="Paraná") %>%
  remove_region(city="Santa Mariana", state="Paraná") %>%
  split_and_add(desc="Santa Maria do Oeste, Paraná", new_city="Santa Maria do Oeste", new_state="Paraná") %>% 
  split_and_add(desc="Santa Mônica, Paraná", new_city="Santa Mônica", new_state="Paraná") %>% 
  split_and_add(desc="Santa Mariana, Paraná", new_city="Santa Mariana", new_state="Paraná")
bra2 <- remove_region(bra2, city="Orobó", state="Pernambuco") %>%
  remove_region(city="Orocó", state="Pernambuco") %>%
  split_and_add(desc="Orobó, Pernambuco", new_city="Orobó", new_state="Pernambuco") %>% 
  split_and_add(desc="Orocó, Pernambuco", new_city="Orocó", new_state="Pernambuco")
bra2 <- remove_region(bra2, city="Várzea", state="Rio Grande do Norte") %>%
  split_and_add(desc="Várzea, Rio Grande do Norte", new_city="Várzea", new_state="Rio Grande do Norte")
bra2 <- remove_region(bra2, city="Baro", state="Rio Grande do Sul") %>%
  remove_region(city="Barão do Triunfo", state="Rio Grande do Sul") %>%
  split_and_add(desc="Barão, Rio Grande do Sul", new_city="Barão", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Barão do Triunfo, Rio Grande do Sul", new_city="Barão do Triunfo", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, city="Capão da Canoa", state="Rio Grande do Sul") %>%
  remove_region(city="Capivari do Sul", state="Rio Grande do Sul") %>%
  remove_region(city="Capão do Leão", state="Rio Grande do Sul") %>%
  remove_region(city="Capela de Santana", state="Rio Grande do Sul") %>%
  remove_region(city="Capitão", state="Rio Grande do Sul") %>%
  split_and_add(desc="Capão da Canoa, Rio Grande do Sul", new_city="Capão da Canoa", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Capivari do Sul, Rio Grande do Sul", new_city="Capivari do Sul", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Capão do Leão, Rio Grande do Sul", new_city="Capão do Leão", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Capela de Santana, Rio Grande do Sul", new_city="Capela de Santana", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Capitão, Rio Grande do Sul", new_city="Capitão", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, city="Nova Palma", state="Rio Grande do Sul") %>%
  remove_region(city="Nova Petrópolis", state="Rio Grande do Sul") %>%
  remove_region(city="Nova Bassano", state="Rio Grande do Sul") %>%
  remove_region(city="Nova Prata", state="Rio Grande do Sul") %>%
  remove_region(city="Nova Pádua", state="Rio Grande do Sul") %>%
  split_and_add(desc="Nova Palma, Rio Grande do Sul", new_city="Nova Palma", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Nova Petrópolis, Rio Grande do Sul", new_city="Nova Petrópolis", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Nova Bassano, Rio Grande do Sul", new_city="Nova Bassano", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Nova Prata, Rio Grande do Sul", new_city="Nova Prata", new_state="Rio Grande do Sul") %>% 
  split_and_add(desc="Nova Pádua, Rio Grande do Sul", new_city="Nova Pádua", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, city="São Valentim do Sul", state="Rio Grande do Sul") %>% 
  split_and_add(desc="São Valentim do Sul, Rio Grande do Sul", new_city="São Valentim do Sul", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, city="São Miguel do Guaporé", state="Rondônia") %>%
  remove_region(city="São Felipe d'Oeste", state="Rondônia") %>%
  remove_region(city="Seringueiras", state="Rondônia") %>%
  remove_region(city="São Francisco do Guaporé", state="Rondônia") %>%
  remove_region(city="Costa Marques", state="Rondônia") %>%
  split_and_add(desc="São Miguel do Guaporé, Rondônia", new_city="São Miguel do Guaporé", new_state="Rondônia") %>% 
  split_and_add(desc="São Felipe d'Oeste, Rondônia", new_city="São Felipe d'Oeste", new_state="Rondônia") %>% 
  split_and_add(desc="Seringueiras, Rondônia", new_city="Seringueiras", new_state="Rondônia") %>% 
  split_and_add(desc="São Francisco do Guaporé, Rondônia", new_city="São Francisco do Guaporé", new_state="Rondônia") %>% 
  split_and_add(desc="Costa Marques, Rondônia", new_city="Costa Marques", new_state="Rondônia")
bra2 <- remove_region(bra2, city="Amparo", state="São Paulo") %>%
  remove_region(city="Américo de Campos", state="São Paulo") %>%
  split_and_add(desc="Amparo, São Paulo", new_city="Amparo", new_state="São Paulo") %>% 
  split_and_add(desc="Américo de Campos, São Paulo", new_city="Américo de Campos", new_state="São Paulo")
bra2 <- remove_region(bra2, city="Guararapes", state="São Paulo") %>%
  remove_region(city="Guaratinguetá", state="São Paulo") %>%
  remove_region(city="Guararema", state="São Paulo") %>%
  split_and_add(desc="Guararapes, São Paulo", new_city="Guararapes", new_state="São Paulo") %>% 
  split_and_add(desc="Guaratinguetá, São Paulo", new_city="Guaratinguetá", new_state="São Paulo") %>% 
  split_and_add(desc="Guararema, São Paulo", new_city="Guararema", new_state="São Paulo")
bra2 <- remove_region(bra2, city="Suzano", state="São Paulo") %>%
  remove_region(city="Suzanápolis", state="São Paulo") %>%
  split_and_add(desc="Suzano, São Paulo", new_city="Suzano", new_state="São Paulo") %>% 
  split_and_add(desc="Suzanápolis, São Paulo", new_city="Suzanápolis", new_state="São Paulo")
bra2 <- remove_region(bra2, city="Rio Paranaíba", state="Minas Gerais") %>%
  remove_region(city="Rio Paranaiba", state="Minas Gerais") %>%
  remove_region(city="Rio Pardo de Minas", state="Minas Gerais") %>%
  remove_region(city="Rio Piracicaba", state="Minas Gerais") %>%
  split_and_add(desc="Rio Paranaíba, Minas Gerais", new_city="Rio Paranaíba", new_state="Minas Gerais") %>% 
  split_and_add(desc="Rio Pardo de Minas, Minas Gerais", new_city="Rio Pardo de Minas", new_state="Minas Gerais") %>% 
  split_and_add(desc="Rio Piracicaba, Minas Gerais", new_city="Rio Piracicaba", new_state="Minas Gerais")
bra2 <- merge_duplicates(bra2, city="Itacambira", state="Minas Gerais")

bra2 <- remove_region(bra2, city="Tabatinga", state="Amazonas") %>%
  remove_region(city="São Sebastião do Uatumã", state="Amazonas") %>%
  split_and_add(desc="São Sebastião do Uatumã, Amazonas", new_city="São Sebastião do Uatumã", new_state="Amazonas") %>% 
  split_and_add(desc="Tabatinga, Amazonas", new_city="Tabatinga", new_state="Amazonas")
bra2 <- remove_region(bra2, city="Alexania", state="Goiás") %>%
  remove_region(city="Alexânia", state="Goiás") %>% 
  split_and_add(desc="Alexânia, Goiás", new_city="Alexânia", new_state="Goiás")
bra2 <- remove_region(bra2, city="Bonfinópolis", state="Goiás") %>%
  remove_region(city="Bonópolis", state="Goiás") %>% 
  split_and_add(desc="Bonfinópolis, Goiás", new_city="Bonfinópolis", new_state="Goiás") %>%
  split_and_add(desc="Bonópolis, Goiás", new_city="Bonópolis", new_state="Goiás")
bra2 <- remove_region(bra2, city="Igarapé do Meio", state="Maranhão") %>%
  remove_region(city="Igarapé Grande", state="Maranhão") %>% 
  split_and_add(desc="Igarapé do Meio, Maranhão", new_city="Igarapé do Meio", new_state="Maranhão") %>%
  split_and_add(desc="Igarapé Grande, Maranhão", new_city="Igarapé Grande", new_state="Maranhão")
bra2 <- remove_region(bra2, city="Água Clara", state="Mato Grosso do Sul") %>% 
  remove_region(city="Costa Rica", state="Mato Grosso do Sul") %>%
  remove_region(city="Chapadão do Sul", state="Mato Grosso do Sul") %>%
  split_and_add(desc="Chapadão do Sul, Mato Grosso do Sul", new_city="Chapadão do Sul", new_state="Mato Grosso do Sul") %>%
  split_and_add(desc="Costa Rica, Mato Grosso do Sul", new_city="Costa Rica", new_state="Mato Grosso do Sul") %>%
  split_and_add(desc="Água Clara, Mato Grosso do Sul", new_city="Água Clara", new_state="Mato Grosso do Sul") %>%
  split_and_add(desc="Paraíso das Águas, Mato Grosso do Sul", new_city="Paraíso das Águas", new_state="Mato Grosso do Sul")
bra2 <- remove_region(bra2, city="Serranópolis de Minas", state="Minas Gerais") %>%
  remove_region(city="Serrania", state="Minas Gerais") %>% 
  split_and_add(desc="Serranópolis de Minas, Minas Gerais", new_city="Serranópolis de Minas", new_state="Minas Gerais") %>%
  split_and_add(desc="Serrania, Minas Gerais", new_city="Serrania", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Divisópolis", state="Minas Gerais") %>%
  remove_region(city="Divisa Nova", state="Minas Gerais") %>% 
  remove_region(city="Divisa Alegre", state="Minas Gerais") %>% 
  split_and_add(desc="Divisópolis, Minas Gerais", new_city="Divisópolis", new_state="Minas Gerais") %>%
  split_and_add(desc="Divisa Nova, Minas Gerais", new_city="Divisa Nova", new_state="Minas Gerais") %>%
  split_and_add(desc="Divisa Alegre, Minas Gerais", new_city="Divisa Alegre", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="José Raydan", state="Minas Gerais") %>%
  remove_region(city="Josenópolis", state="Minas Gerais") %>% 
  remove_region(city="José Gonçalves de Minas", state="Minas Gerais") %>% 
  split_and_add(desc="José Raydan, Minas Gerais", new_city="José Raydan", new_state="Minas Gerais") %>%
  split_and_add(desc="Josenópolis, Minas Gerais", new_city="Josenópolis", new_state="Minas Gerais") %>%
  split_and_add(desc="José Gonçalves de Minas, Minas Gerais", new_city="José Gonçalves de Minas", new_state="Minas Gerais")

bra2 <- remove_region(bra2, city="Concórdia do Pará", state="Pará") %>%
  remove_region(city="Conceição do Araguaia", state="Pará") %>% 
  split_and_add(desc="Concórdia do Pará, Pará", new_city="Concórdia do Pará", new_state="Pará") %>%
  split_and_add(desc="Conceição do Araguaia, Pará", new_city="Conceição do Araguaia", new_state="Pará")
bra2 <- remove_region(bra2, city="Cuitegi", state="Paraíba") %>%
  remove_region(city="Cuité de Mamanguape", state="Paraíba") %>% 
  split_and_add(desc="Cuitegi, Paraíba", new_city="Cuitegi", new_state="Paraíba") %>%
  split_and_add(desc="Cuité de Mamanguape, Paraíba", new_city="Cuité de Mamanguape", new_state="Paraíba")
bra2 <- remove_region(bra2, city="Santo Antônio de Lisboa", state="Piauí") %>%
  remove_region(city="Santo Antônio dos Milagres", state="Piauí") %>% 
  split_and_add(desc="Santo Antônio de Lisboa, Piauí", new_city="Santo Antônio de Lisboa", new_state="Piauí") %>%
  split_and_add(desc="Santo Antônio dos Milagres, Piauí", new_city="Santo Antônio dos Milagres", new_state="Piauí")
bra2 <- remove_region(bra2, city="Caicó", state="Rio Grande do Norte") %>%
  remove_region(city="Caiçara do Norte", state="Rio Grande do Norte") %>% 
  remove_region(city="Caiçara do Rio do Vento", state="Rio Grande do Norte") %>% 
  split_and_add(desc="Caicó, Rio Grande do Norte", new_city="Caicó", new_state="Rio Grande do Norte") %>%
  split_and_add(desc="Caiçara do Norte, Rio Grande do Norte", new_city="Caiçara do Norte", new_state="Rio Grande do Norte") %>%
  split_and_add(desc="Caiçara do Rio do Vento, Rio Grande do Norte", new_city="Caiçara do Rio do Vento", new_state="Rio Grande do Norte")
bra2 <- remove_region(bra2, city="Balneário Arroio do Silva", state="Santa Catarina") %>%
  remove_region(city="Balneário Camboriú", state="Santa Catarina") %>% 
  remove_region(city="Balneário Barra do Sul", state="Santa Catarina") %>% 
  split_and_add(desc="Balneário Arroio do Silva, Santa Catarina", new_city="Balneário Arroio do Silva", new_state="Santa Catarina") %>%
  split_and_add(desc="Balneário Camboriú, Santa Catarina", new_city="Balneário Camboriú", new_state="Santa Catarina") %>%
  split_and_add(desc="Balneário Barra do Sul, Santa Catarina", new_city="Balneário Barra do Sul", new_state="Santa Catarina")
bra2 <- remove_region(bra2, city="Itapiranga", state="Santa Catarina") %>%
  split_and_add(desc="Itapiranga, Santa Catarina", new_city="Itapiranga", new_state="Santa Catarina")
bra2 <- remove_region(bra2, city="Figueirópolis", state="Tocantins") %>%
  remove_region(city="Filadélfia", state="Tocantins") %>% 
  remove_region(city="Fortaleza do Tabocão", state="Tocantins") %>% 
  remove_region(city="Fátima", state="Tocantins") %>% 
  split_and_add(desc="Figueirópolis, Tocantins", new_city="Figueirópolis", new_state="Tocantins") %>%
  split_and_add(desc="Filadélfia, Tocantins", new_city="Filadélfia", new_state="Tocantins") %>%
  split_and_add(desc="Tabocão, Tocantins", new_city="Tabocão", new_state="Tocantins") %>%
  split_and_add(desc="Fátima, Tocantins", new_city="Fátima", new_state="Tocantins")
bra2 <- remove_region(bra2, city="Estação", state="Rio Grande do Sul") %>%
  remove_region(city="Estância Velha", state="Rio Grande do Sul") %>% 
  split_and_add(desc="Estação, Rio Grande do Sul", new_city="Estação", new_state="Rio Grande do Sul") %>%
  split_and_add(desc="Estância Velha, Rio Grande do Sul", new_city="Estância Velha", new_state="Rio Grande do Sul")
bra2 <- remove_region(bra2, city="Formoso do Araguaia", state="Tocantins") %>%
  split_and_add(desc="Formoso do Araguaia, Tocantins", new_city="Formoso do Araguaia", new_state="Tocantins")
bra2 <- remove_region(bra2, city="Nova Brasilândia d'Oeste", state="Rondônia") %>%
  split_and_add(desc="Nova Brasilândia d'Oeste, Rondônia", new_city="Nova Brasilândia d'Oeste", new_state="Rondônia")
bra2 <- remove_region(bra2, city="Cajari", state="Maranhão") %>%
  split_and_add(desc="Cajari, Maranhão", new_city="Cajari", new_state="Maranhão")
bra2 <- remove_region(bra2, city="Canudos", state="Bahia") %>%
  remove_region(city="Coronel João Sá", state="Bahia") %>%
  remove_region(city="Pedro Alexandre", state="Bahia") %>%
  remove_region(city="Santa Brígida", state="Bahia") %>%
  remove_region(city="Paripiranga", state="Bahia") %>%
  split_and_add(desc="Canudos, Bahia", new_city="Canudos", new_state="Bahia") %>%
  split_and_add(desc="Coronel João Sá, Bahia", new_city="Coronel João Sá", new_state="Bahia") %>%
  split_and_add(desc="Santa Brígida, Bahia", new_city="Santa Brígida", new_state="Bahia") %>%
  split_and_add(desc="Paripiranga, Bahia", new_city="Paripiranga", new_state="Bahia") %>%
  split_and_add(desc="Pedro Alexandre, Bahia", new_city="Pedro Alexandre", new_state="Bahia")

bra2 <- remove_region(bra2, city="Jaíba", state="Minas Gerais") %>%
  remove_region(city="Verdelândia", state="Minas Gerais") %>%
  split_and_add(desc="Jaíba, Minas Gerais", new_city="Jaíba", new_state="Minas Gerais") %>%
  split_and_add(desc="Verdelândia, Minas Gerais", new_city="Verdelândia", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Almenara", state="Minas Gerais") %>%
  split_and_add(desc="Almenara, Minas Gerais", new_city="Almenara", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Bom Jardim de Minas", state="Minas Gerais") %>%
  split_and_add(desc="Bom Jardim de Minas, Minas Gerais", new_city="Bom Jardim de Minas", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Bom Jardim de Minas", state="Minas Gerais") %>%
  split_and_add(desc="Bom Jardim de Minas, Minas Gerais", new_city="Bom Jardim de Minas", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Itajubá", state="Minas Gerais") %>%
  split_and_add(desc="Itajubá, Minas Gerais", new_city="Itajubá", new_state="Minas Gerais")
bra2 <- remove_region(bra2, city="Araucária", state="Paraná") %>%
  split_and_add(desc="Araucária, Paraná", new_city="Araucária", new_state="Paraná")

bra2 <- remove_region(bra2, city="Corumbá", state="Mato Grosso do Sul") %>%
  split_and_add(desc="Corumbá, Mato Grosso do Sul", new_city="Corumbá", new_state="Mato Grosso do Sul")

#################################################################################################
# here we fix the names of regions in the GADM data set so that they are consistent with the DENV
# case data. We also remove some islands for which there is no case data

# fix_region_name takes a string for the state and a string for the city
# and returns a string vector of the corrected names. 
fix_region_name <- function(state, city) {
  if (is.na(state) && is.na(city)) return(c(state, city))
  
  # fix name
  if (state=="Acre" && city=="Assis Brazil") city<-"Assis Brasil"
  else if (state=="Alagoas" && city=="Ataléia") city<-"Atalaia"
  else if (state=="Alagoas" && city=="Maceió (capital)") city<-"Maceió"
  else if (state=="Alagoas" && city=="Marechal deodoro") city<-"Marechal Deodoro"
  else if (state=="Alagoas" && city=="Palmeira dos índios") city<-"Palmeira dos Índios"
  else if (state=="Amapá" && city=="Macapa") city<-"Macapá"
  else if (state=="Amapá" && city=="Pedra Branca do Amaparí") city<-"Pedra Branca do Amapari"
  else if (state=="Amazonas" && city=="Manicore") city<-"Manicoré"
  else if (state=="Amazonas" && city=="São Gabriel de Cahoeira") city<-"São Gabriel da Cachoeira"
  else if (state=="Bahia" && city=="Alcobaca") city<-"Alcobaça"
  else if (state=="Bahia" && city=="América dourada") city<-"América Dourada"
  else if (state=="Bahia" && city=="Andarai") city<-"Andaraí"
  else if (state=="Bahia" && city=="Barra da Choça") city<-"Barra do Choça"
  else if (state=="Bahia" && city=="Dias d'vila") city<-"Dias d'Ávila"
  else if (state=="Bahia" && city=="Lajedao") city<-"Lajedão"
  else if (state=="Bahia" && city=="Livramento do Brumado") city<-"Livramento de Nossa Senhora"
  else if (state=="Bahia" && city=="Madre de deus") city<-"Madre de Deus"
  else if (state=="Bahia" && city=="Oliveria dos Brejinhos") city<-"Oliveira dos Brejinhos"
  else if (state=="Bahia" && city=="Pau Brazil") city<-"Pau Brasil"
  else if (state=="Bahia" && city=="Riachao do Jacuipe") city<-"Riachão do Jacuípe"
  else if (state=="Bahia" && city=="Serra dourada") city<-"Serra Dourada"
  else if (state=="Bahia" && city=="Araças") {city<-"Araçás"}
  else if (state=="Bahia" && city=="Muquém de São Francisco") {city<-"Muquém do São Francisco"}
  else if (state=="Bahia" && city=="Santa Teresinha") {city<-"Santa Terezinha"}
  else if (state=="Ceará" && city=="Acarapé") city<-"Acarape"
  else if (state=="Ceará" && city=="Araçoiaba") city<-"Aracoiaba"
  else if (state=="Ceará" && city=="Caririaçú") city<-"Caririaçu"
  else if (state=="Ceará" && city=="Ipú") city<-"Ipu"
  else if (state=="Ceará" && city=="Itarumã") city<-"Itarema"
  else if (state=="Ceará" && city=="Misso Velha") city<-"Missão Velha"
  else if (state=="Ceará" && city=="Pacajús") city<-"Pacajus"
  else if (state=="Ceará" && city=="São João do Belm") city<-"Mauriti"
  else if (state=="Ceará" && city=="São Luiz do Curu") city<-"São Luís do Curu"
  else if (state=="Ceará" && city=="Quixada") city<-"Quixadá"
  else if (state=="Espírito Santo" && city=="Guia Branca") city<-"Águia Branca"
  else if (state=="Espírito Santo" && city=="Vitoria") city<-"Vitória"
  else if (state=="Espírito Santo" && city=="Atilio Vivacqua") {city<-"Atílio Vivácqua"}
  else if (state=="Goiás" && city=="Americano do Brazil") city<-"Americano do Brasil"
  else if (state=="Goiás" && city=="Aparecida do Rio doce") city<-"Aparecida do Rio Doce"
  else if (state=="Goiás" && city=="Araporã") state="Minas Gerais"
  else if (state=="Goiás" && city=="Brasabrantes") city<-"Brazabrantes"
  else if (state=="Goiás" && city=="Cachoeira de Goias") city<-"Cachoeira de Goiás"
  else if (state=="Goiás" && city=="Cachoeira dourada") city<-"Cachoeira Dourada"
  else if (state=="Goiás" && city=="Chapadinha") city<-"Chapadão do Céu"
  else if (state=="Goiás" && city=="Goiania") city<-"Goiânia"
  else if (state=="Goiás" && city=="Itaruma") city<-"Itarumã"
  else if (state=="Goiás" && city=="Mateira") city<-"Paranaiguara"
  else if (state=="Goiás" && city=="Novo Brazil") city<-"Novo Brasil"
  else if (state=="Goiás" && city=="Porteiro") city<-"Porteirão"
  else if (state=="Goiás" && city=="Santa Filomena do Maranhão") city<-"Santa Fé de Goiás"
  else if (state=="Goiás" && city=="Santa Rita de Araguaia") city<-"Santa Rita do Araguaia"
  else if (state=="Goiás" && city=="Santa Rita do Novo destino") city<-"Santa Rita do Novo Destino"
  else if (state=="Goiás" && city=="São Francisco de Goias") city<-"São Francisco de Goiás"
  else if (state=="Goiás" && city=="Varjao") city<-"Varjão"
  else if (state=="Maranhão" && city=="Alto Alegre do Maranho") city<-"Alto Alegre do Maranhão"
  else if (state=="Maranhão" && city=="Alto Parnaiba") city<-"Alto Parnaíba"
  else if (state=="Maranhão" && city=="Amapá do Maranho") city<-"Amapá do Maranhão"
  else if (state=="Maranhão" && city=="Anapuros") city<-"Anapurus"
  else if (state=="Maranhão" && city=="Bom Jardin") city<-"Bom Jardim"
  else if (state=="Maranhão" && city=="Chapadão do Céu") city<-"Chapadinha"
  else if (state=="Maranhão" && city=="Humberto Campos") city<-"Humberto de Campos"
  else if (state=="Maranhão" && city=="Mates do Norte") city<-"Matões do Norte"
  else if (state=="Maranhão" && city=="Peri-Mirim") city<-"Peri Mirim"
  else if (state=="Maranhão" && city=="Santa Fé de Goiás") city<-"Santa Filomena do Maranhão"
  else if (state=="Maranhão" && city=="São Luis") city<-"São Luís"
  else if (state=="Maranhão" && city=="São Luis Gonzaga do Maranhao") city<-"São Luís Gonzaga do Maranhão"
  else if (state=="Maranhão" && city=="Victorino Freire") city<-"Vitorino Freire"
  else if (state=="Mato Grosso" && city=="Barra dos Bugre") city<-"Barra do Bugres"
  else if (state=="Mato Grosso" && city=="Cabixi") state="Rondônia"
  else if (state=="Mato Grosso" && city=="CanaBrava do Norte") city<-"Canabrava do Norte"
  else if (state=="Mato Grosso" && city=="Cuiaba") city<-"Cuiabá"
  else if (state=="Mato Grosso" && city=="Figueirópolis d'Oeste") city<-"Figueirópolis D'Oeste"
  else if (state=="Mato Grosso" && city=="Glória d'Oeste") city<-"Glória D'Oeste"
  else if (state=="Mato Grosso" && city=="Luciára") city<-"Luciara"
  else if (state=="Mato Grosso" && city=="Mato Grosso") city<-"Vila Bela da Santíssima Trindade"
  else if (state=="Mato Grosso" && city=="São Félix Xingu") city<-"São Félix do Araguaia"
  else if (state=="Mato Grosso" && city=="Poxoréo") {city<-"Poxoréu"}
  else if (state=="Mato Grosso do Sul" && city=="Bataiporã") city<-"Batayporã"
  else if (state=="Mato Grosso do Sul" && city=="Fatima do Sul") city<-"Fátima do Sul"
  else if (state=="Minas Gerais" && city=="Alto Rio doce") city<-"Alto Rio Doce"
  else if (state=="Minas Gerais" && city=="Antonio Prado de Minas") city<-"Antônio Prado de Minas"
  else if (state=="Minas Gerais" && city=="Aracai") city<-"Araçaí"
  else if (state=="Minas Gerais" && city=="Araujos") city<-"Araújos"
  else if (state=="Minas Gerais" && city=="Ataleia") city<-"Ataléia"
  else if (state=="Minas Gerais" && city=="Bandiera do Sul") city<-"Bandeira do Sul"
  else if (state=="Minas Gerais" && city=="Bom despacho") city<-"Bom Despacho"
  else if (state=="Minas Gerais" && city=="Bras Pires") city<-"Brás Pires"
  else if (state=="Minas Gerais" && city=="Cachoeira de Pajes") city<-"Cachoeira de Pajeú"
  else if (state=="Minas Gerais" && city=="Cachoeira dourada") city<-"Cachoeira Dourada"
  else if (state=="Minas Gerais" && city=="Caconde") state="São Paulo"
  else if (state=="Minas Gerais" && city=="Campos Verdes de Goiás") city<-"Cana Verde"
  else if (state=="Minas Gerais" && city=="Caravalhopolis") city<-"Carvalhópolis"
  else if (state=="Minas Gerais" && city=="Carmo do Paranaiba") city<-"Carmo do Paranaíba"
  else if (state=="Minas Gerais" && city=="Cassiterita") city<-"Conceição da Barra de Minas"
  else if (state=="Minas Gerais" && city=="Chale") city<-"Chalé"
  else if (state=="Minas Gerais" && city=="Conceição do Para") city<-"Conceição do Pará"
  else if (state=="Minas Gerais" && city=="Couto de Magalhães") city<-"Couto de Magalhães de Minas"
  else if (state=="Minas Gerais" && city=="Curral de dentro") city<-"Curral de Dentro"
  else if (state=="Minas Gerais" && city=="Divinolandia de Minas") city<-"Divinolândia de Minas"
  else if (state=="Minas Gerais" && city=="Ecoporanga") state="Espírito Santo"
  else if (state=="Minas Gerais" && city=="Estrela dalva") city<-"Estrela Dalva"
  else if (state=="Minas Gerais" && city=="Felisberto Caldeira") city<-"São Gonçalo do Rio Preto"
  else if (state=="Minas Gerais" && city=="Gouvea") city<-"Gouveia"
  else if (state=="Minas Gerais" && city=="Guaranesia") city<-"Guaranésia"
  else if (state=="Minas Gerais" && city=="Guimarania") city<-"Guimarânia"
  else if (state=="Minas Gerais" && city=="Iaçu") city<-"Iapu"
  else if (state=="Minas Gerais" && city=="Lagoa dourada") city<-"Lagoa Dourada"
  else if (state=="Minas Gerais" && city=="Madre de deus de Minas") city<-"Madre de Deus de Minas"
  else if (state=="Minas Gerais" && city=="Pedra dourada") city<-"Pedra Dourada"
  else if (state=="Minas Gerais" && city=="Piedade do Ponte Nova") city<-"Piedade de Ponte Nova"
  else if (state=="Minas Gerais" && city=="Pingo d'Água") city<-"Pingo-d'Água"
  else if (state=="Minas Gerais" && city=="Pirauba") city<-"Piraúba"
  else if (state=="Minas Gerais" && city=="Piui") city<-"Piumhi"
  else if (state=="Minas Gerais" && city=="Queluzita") city<-"Queluzito"
  else if (state=="Minas Gerais" && city=="Rio doce") city<-"Rio Doce"
  else if (state=="Minas Gerais" && city=="Santana do deserto") city<-"Santana do Deserto"
  else if (state=="Minas Gerais" && city=="São Francisco de Oliveira") city<-"São Francisco de Paula"
  else if (state=="Minas Gerais" && city=="São Sebastio da Vargem Alegre") city<-"São Sebastião da Vargem Alegre"
  else if (state=="Minas Gerais" && city=="Serra da Saudad") city<-"Serra da Saudade"
  else if (state=="Minas Gerais" && city=="Turvolandia") city<-"Turvolândia"
  else if (state=="Pará" && city=="Almerim") city<-"Almeirim"
  else if (state=="Pará" && city=="Anajas") city<-"Anajás"
  else if (state=="Pará" && city=="Bagé") city<-"Bagre"
  else if (state=="Pará" && city=="Braganga") city<-"Bragança"
  else if (state=="Pará" && city=="Brazil Novo") city<-"Brasil Novo"
  else if (state=="Pará" && city=="Me do Rio") city<-"Mãe do Rio"
  else if (state=="Pará" && city=="Melgaco") city<-"Melgaço"
  else if (state=="Pará" && city=="Pau d'Arco") city<-"Pau D'Arco"
  else if (state=="Pará" && city=="Peixe Boi") city<-"Peixe-Boi"
  else if (state=="Pará" && city=="Eldorado dos Carajás") {city<-"Eldorado do Carajás"}
  else if (state=="Pará" && city=="Santa Isabel do Pará") {city<-"Santa Izabel do Pará"}
  else if (state=="Paraíba" && city=="Santarém") city<-"Joca Claudino"
  else if (state=="Paraíba" && city=="Aracagi") city<-"Araçagi"
  else if (state=="Paraíba" && city=="Boqueirao dos Cochos") city<-"Igaracy"
  else if (state=="Paraíba" && city=="Cachoeira dos índios") city<-"Cachoeira dos Índios"
  else if (state=="Paraíba" && city=="Cacimba de dentro") city<-"Cacimba de Dentro"
  else if (state=="Paraíba" && city=="Desterro de Malta") city<-"Vista Serrana"
  else if (state=="Paraíba" && city=="Lagoa de dentro") city<-"Lagoa de Dentro"
  else if (state=="Paraíba" && city=="Mongeiro") city<-"Mogeiro"
  else if (state=="Paraíba" && city=="Passabém") city<-"Passagem"
  else if (state=="Paraíba" && city=="Pedra Lavadra") city<-"Pedra Lavrada"
  else if (state=="Paraíba" && city=="Riacho") city<-"Riachão"
  else if (state=="Paraíba" && city=="Ricaho dos Cavalos") city<-"Riacho dos Cavalos"
  else if (state=="Paraíba" && city=="São Domingos de Pombal") city<-"São Domingos"
  else if (state=="Paraíba" && city=="São José do Belmonte") city<-"São José do Brejo do Cruz"
  else if (state=="Paraíba" && city=="São Miguel Taipu") city<-"São Miguel de Taipu"
  else if (state=="Paraíba" && city=="Seridó") {city<-"São Vicente do Seridó"}
  else if (state=="Paraíba" && city=="Várzea Branca") state="Piauí"
  else if (state=="Paraíba" && city=="Quixabá") city="Quixaba"
  else if (state=="Paraná" && city=="Amapora") city<-"Amaporã"
  else if (state=="Paraná" && city=="Antonio Olinto") city<-"Antônio Olinto"
  else if (state=="Paraná" && city=="Arapu") city<-"Arapuã"
  else if (state=="Paraná" && city=="Assis Chateaubri") city<-"Assis Chateaubriand"
  else if (state=="Paraná" && city=="Ataléia") city<-"Atalaia"
  else if (state=="Paraná" && city=="Campo") city<-"Campo Largo"
  else if (state=="Paraná" && city=="Cêrro Azul") city<-"Cerro Azul"
  else if (state=="Paraná" && city=="Conselheiro Mayrinck") city<-"Conselheiro Mairinck"
  else if (state=="Paraná" && city=="Coronel domingos Soares") city<-"Coronel Domingos Soares"
  else if (state=="Paraná" && city=="Diamante d'Oeste") city<-"Diamante D'Oeste"
  else if (state=="Paraná" && city=="Itambaraca") city<-"Itambaracá"
  else if (state=="Paraná" && city=="Jabuti") city<-"Jaboti"
  else if (state=="Paraná" && city=="Luiziânia") city<-"Luiziana"
  else if (state=="Paraná" && city=="Rancho Alegre d'Oeste") city<-"Rancho Alegre D'Oeste"
  else if (state=="Paraná" && city=="Salto do Londra") city<-"Salto do Lontra"
  else if (state=="Paraná" && city=="Santa Cruz de Monte Caste") city<-"Santa Cruz de Monte Castelo"
  else if (state=="Paraná" && city=="Santa Isabel do Oeste") city<-"Santa Izabel do Oeste"
  else if (state=="Paraná" && city=="Santo Antonio da Platina") city<-"Santo Antônio da Platina"
  else if (state=="Paraná" && city=="Santo Antonio do Caiuá") city<-"Santo Antônio do Caiuá"
  else if (state=="Paraná" && city=="Santo Antonio do Paraíso") city<-"Santo Antônio do Paraíso"
  else if (state=="Paraná" && city=="São Antonio de Sudoeste") city<-"Santo Antônio do Sudoeste"
  else if (state=="Paraná" && city=="Texeira Soares") city<-"Teixeira Soares"
  else if (state=="Paraná" && city=="Tibaji") city<-"Tibagi"
  else if (state=="Paraná" && city=="Venceslau Bras") city<-"Wenceslau Braz"
  else if (state=="Paraná" && city=="Vila Alta") city<-"Alto Paraíso"
  else if (state=="Pernambuco" && city=="Barra de Guabira") city<-"Barra de Guabiraba"
  else if (state=="Pernambuco" && city=="Brejo da Madre de deus") city<-"Brejo da Madre de Deus"
  else if (state=="Pernambuco" && city=="Belém de São Francisco") {city<-"Belém do São Francisco"}
  else if (state=="Pernambuco" && city=="Cabo") city<-"Cabo de Santo Agostinho"
  else if (state=="Pernambuco" && city=="Cachoerinha") city<-"Cachoeirinha"
  else if (state=="Pernambuco" && city=="Cortes") city<-"Cortês"
  else if (state=="Pernambuco" && city=="Goianá") city<-"Goiana"
  else if (state=="Pernambuco" && city=="Iguaraci") city<-"Iguaracy"
  else if (state=="Pernambuco" && city=="Igaracu") city<-"Igarassu"
  else if (state=="Pernambuco" && city=="Itambaracá") city<-"Ilha de Itamaracá"
  else if (state=="Pernambuco" && city=="Jupiá") city<-"Jupi"
  else if (state=="Pernambuco" && city=="Quixabá") city<-"Quixaba"
  else if (state=="Pernambuco" && city=="Salidao") city<-"Solidão"
  else if (state=="Pernambuco" && city=="São João do Belmonte") city<-"São José do Belmonte"
  else if (state=="Pernambuco" && city=="São Joaquin do Monte") city<-"São Joaquim do Monte"
  else if (state=="Pernambuco" && city=="Sitio dos Moreiras") city<-"Moreilândia"
  else if (state=="Pernambuco" && city=="Tambe") city<-"Itambé"
  else if (state=="Pernambuco" && city=="São Vicente Ferrer") {city<-"São Vicente Férrer"}
  else if (state=="Pernambuco" && city=="Lagoa do Itaenga") {city<-"Lagoa de Itaenga"}
  else if (state=="Piauí" && city=="Barra d'Alcântara") city<-"Barra D'Alcântara"
  else if (state=="Piauí" && city=="Brazileira") city<-"Brasileira"
  else if (state=="Piauí" && city=="Francisco Macêdo") city<-"Francisco Macedo"
  else if (state=="Piauí" && city=="Morro Cabeça No Tempo") city<-"Morro Cabeça no Tempo"
  else if (state=="Piauí" && city=="Olho d'água do Piauí") city<-"Olho D'Água do Piauí"
  else if (state=="Piauí" && city=="Pau D'Arco do Piauí") {city<-"São Valério"}
  else if (state=="Piauí" && city=="Pedro Li") city<-"Pedro II"
  else if (state=="Piauí" && city=="Santa Cruz do Piaui") city<-"Santa Cruz do Piauí"
  else if (state=="Piauí" && city=="São João Piaui") city<-"São João do Piauí"
  else if (state=="Piauí" && city=="São Juliao") city<-"São Julião"
  else if (state=="Piauí" && city=="São Miguel Tapuio") city<-"São Miguel do Tapuio"
  else if (state=="Rio de Janeiro" && city=="Campos") city<-"Campos dos Goytacazes"
  else if (state=="Rio de Janeiro" && city=="Carepebus") city<-"Carapebus"
  else if (state=="Rio de Janeiro" && city=="Conceicao Macabu") city<-"Conceição de Macabu"
  else if (state=="Rio de Janeiro" && city=="Engenheiro Paulo de Front") city<-"Engenheiro Paulo de Frontin"
  else if (state=="Rio de Janeiro" && city=="Parati") city<-"Paraty"
  else if (state=="Rio de Janeiro" && city=="Trajano de Morais") city<-"Trajano de Moraes"
  else if (state=="Rio de Janeiro" && city=="Valencia") city<-"Valença"
  else if (state=="Rio Grande do Norte" && city=="Fernando de Noronha") state<-"Pernambuco"
  else if (state=="Rio Grande do Norte" && city=="Governador Dix-Sept Rosad") city<-"Governador Dix-Sept Rosado"
  else if (state=="Rio Grande do Norte" && city=="Groaíras") city<-"Grossos"
  else if (state=="Rio Grande do Norte" && city=="Jardim-Piranhas") city<-"Jardim de Piranhas"
  else if (state=="Rio Grande do Norte" && city=="Junco") city<-"Messias Targino"
  else if (state=="Rio Grande do Norte" && city=="Lagoa de Anta") city<-"Lagoa d'Anta"
  else if (state=="Rio Grande do Norte" && city=="Lagoas de Velhos") city<-"Lagoa de Velhos"
  else if (state=="Rio Grande do Norte" && city=="Passabém") city<-"Passagem"
  else if (state=="Rio Grande do Norte" && city=="Poço Dantas") state="Paraíba"
  else if (state=="Rio Grande do Norte" && city=="Santana") city<-"Santana do Seridó"
  else if (state=="Rio Grande do Norte" && city=="São Miguel de Touros") city<-"São Miguel do Gostoso"
  else if (state=="Rio Grande do Norte" && city=="Presidente Juscelino") city <- "Serra Caiada"
  else if (state=="Rio Grande do Sul" && city=="Baje") city<-"Bagé"
  else if (state=="Rio Grande do Sul" && city=="Barao de Cotegipe") city<-"Barão de Cotegipe"
  else if (state=="Rio Grande do Sul" && city=="Cacique doble") city<-"Cacique Doble"
  else if (state=="Rio Grande do Sul" && city=="Camagua") city<-"Camaquã"
  else if (state=="Rio Grande do Sul" && city=="Campo Real") city<-"Não-Me-Toque"
  else if (state=="Rio Grande do Sul" && city=="Chiapeta") city<-"Chiapetta"
  else if (state=="Rio Grande do Sul" && city=="Dilermano de Aguiar") city<-"Dilermando de Aguiar"
  else if (state=="Rio Grande do Sul" && city=="Erval") city<-"Herval"
  else if (state=="Rio Grande do Sul" && city=="Inhacor") city<-"Inhacorá"
  else if (state=="Rio Grande do Sul" && city=="Maçambara") city<-"Maçambará"
  else if (state=="Rio Grande do Sul" && city=="Maraú") city<-"Marau"
  else if (state=="Rio Grande do Sul" && city=="Marcionilio Dias") city<-"Marcelino Ramos"
  else if (state=="Rio Grande do Sul" && city=="Maximiliano de Almaeida") city<-"Maximiliano de Almeida"
  else if (state=="Rio Grande do Sul" && city=="Palmitinhos") city<-"Palmitinho"
  else if (state=="Rio Grande do Sul" && city=="Portao") city<-"Portão"
  else if (state=="Rio Grande do Sul" && city=="Rio dos índios") city<-"Rio dos Índios"
  else if (state=="Rio Grande do Sul" && city=="Santana do Livramento") city<-"Sant'Ana do Livramento"
  else if (state=="Rio Grande do Sul" && city=="Santo Ángelo") city<-"Santo Ângelo"
  else if (state=="Rio Grande do Sul" && city=="São Miguel das Misses") city<-"São Miguel das Missões"
  else if (state=="Rio Grande do Sul" && city=="urea") city<-"Áurea"
  else if (state=="Rio Grande do Sul" && city=="Vitória das Misses") city<-"Vitória das Missões"
  else if (state=="Rio Grande do Sul" && city=="Vespasiano Correa") {city<-"Vespasiano Corrêa"}
  else if (state=="Rio Grande do Sul" && city=="Restinga Seca") {city<-"Restinga Sêca"}
  else if (state=="Rondônia" && city=="Alta Floresta d'Oeste") city<-"Alta Floresta D'Oeste"
  else if (state=="Rondônia" && city=="Alvorada d'Oeste") city<-"Alvorada D'Oeste"
  else if (state=="Rondônia" && city=="Espigão d'Oeste") city<-"Espigão D'Oeste"
  else if (state=="Rondônia" && city=="Machadinho") city<-"Machadinho D'Oeste"
  else if (state=="Rondônia" && city=="Alta Floresta d'Oeste") city<-"Alta Floresta D'Oeste"
  else if (state=="Rondônia" && city=="Nova Brasilândia d'Oeste") city<-"Nova Brasilândia D'Oeste"
  else if (state=="Rondônia" && city=="Santa Luzia d'Oeste") city<-"Santa Luzia D'Oeste"
  else if (state=="Rondônia" && city=="São Felipe d'Oeste") city<-"São Felipe D'Oeste"
  else if (state=="Santa Catarina" && city=="Florianopolis") city<-"Florianópolis"
  else if (state=="Santa Catarina" && city=="Gravataí") city<-"Gravatal"
  else if (state=="Santa Catarina" && city=="Ipirá") city<-"Ipira"
  else if (state=="Santa Catarina" && city=="Joinvile") city<-"Joinville"
  else if (state=="Santa Catarina" && city=="Mirim doce") city<-"Mirim Doce"
  else if (state=="Santa Catarina" && city=="Orleaes") city<-"Orleans"
  else if (state=="Santa Catarina" && city=="Paulo Lopez") city<-"Paulo Lopes"
  else if (state=="Santa Catarina" && city=="Piçarras") city<-"Balneário Piçarras"
  else if (state=="Santa Catarina" && city=="Ponta Alta") city<-"Ponte Alta"
  else if (state=="Santa Catarina" && city=="Presidente Castelo Branco") city<-"Presidente Castello Branco"
  else if (state=="Santa Catarina" && city=="Lauro Muller") {city<-"Lauro Müller"}
  else if (state=="Santa Catarina" && city=="Sul Brazil") city<-"Sul Brasil"
  else if (state=="São Paulo" && city=="Aguai") city<-"Aguaí"
  else if (state=="São Paulo" && city=="Alfredo Marconde") city<-"Alfredo Marcondes"
  else if (state=="São Paulo" && city=="Analandia") city<-"Analândia"
  else if (state=="São Paulo" && city=="Aparecida doeste") city<-"Aparecida d'Oeste"
  else if (state=="São Paulo" && city=="Aruja") city<-"Arujá"
  else if (state=="São Paulo" && city=="Avare") city<-"Avaré"
  else if (state=="São Paulo" && city=="Boa Esperanca do Sul") city<-"Boa Esperança do Sul"
  else if (state=="São Paulo" && city=="Bon Jesus dos Perdoes") city<-"Bom Jesus dos Perdões"
  else if (state=="São Paulo" && city=="Brauna") city<-"Braúna"
  else if (state=="São Paulo" && city=="Brodosqui") city<-"Brodowski"
  else if (state=="São Paulo" && city=="Catigua") city<-"Catiguá"
  else if (state=="São Paulo" && city=="Dulcinopolis") city<-"Dolcinópolis"
  else if (state=="São Paulo" && city=="Estrela do Oeste") city<-"Estrela d'Oeste"
  else if (state=="São Paulo" && city=="Embu") {city<-"Embu das Artes"}
  else if (state=="São Paulo" && city=="Ferno") city<-"Fernão"
  else if (state=="São Paulo" && city=="Ferraz de Vascon") city<-"Ferraz de Vasconcelos"
  else if (state=="São Paulo" && city=="Florínia") {city<-"Florínea"}
  else if (state=="São Paulo" && city=="Guarani do Oeste") city<-"Guarani d'Oeste"
  else if (state=="São Paulo" && city=="Guzolandia") city<-"Guzolândia"
  else if (state=="São Paulo" && city=="Ipaucu") city<-"Ipaussu"
  else if (state=="São Paulo" && city=="Jabuticabal") city<-"Jaboticabal"
  else if (state=="São Paulo" && city=="Luisiania") city<-"Luiziânia"
  else if (state=="São Paulo" && city=="Lupercio") city<-"Lupércio"
  else if (state=="São Paulo" && city=="Macedonia") city<-"Macedônia"
  else if (state=="São Paulo" && city=="Itaóca") {city<-"Itaoca"}
  else if (state=="São Paulo" && city=="Mombaça") city<-"Mombuca"
  else if (state=="São Paulo" && city=="Orlandia") city<-"Orlândia"
  else if (state=="São Paulo" && city=="Palmeira do Oeste") city<-"Palmeira d'Oeste"
  else if (state=="São Paulo" && city=="Paranaparema") city<-"Paranapanema"
  else if (state=="São Paulo" && city=="Piracununga") city<-"Pirassununga"
  else if (state=="São Paulo" && city=="Pontes Gestral") city<-"Pontes Gestal"
  else if (state=="São Paulo" && city=="Quitana") city<-"Quintana"
  else if (state=="São Paulo" && city=="Ribeirão dos índios") city<-"Ribeirão dos Índios"
  else if (state=="São Paulo" && city=="Ribeirao Preto") city<-"Ribeirão Preto"
  else if (state=="São Paulo" && city=="Salto do Pirapora") city<-"Salto de Pirapora"
  else if (state=="São Paulo" && city=="Santa Clara do Oeste") city<-"Santa Clara d'Oeste"
  else if (state=="São Paulo" && city=="Santa Lucia") city<-"Santa Lúcia"
  else if (state=="São Paulo" && city=="Santa Rita do Oeste") city<-"Santa Rita d'Oeste"
  else if (state=="São Paulo" && city=="São João das Duas Ponte") city<-"São João das Duas Pontes"
  else if (state=="São Paulo" && city=="Sertaozinho") city<-"Sertãozinho"
  else if (state=="São Paulo" && city=="Tejupa") city<-"Tejupá"
  else if (state=="Sergipe" && city=="Buquim") city<-"Boquim"
  else if (state=="Sergipe" && city=="Itaporanga dajuda") city<-"Itaporanga d'Ajuda"
  else if (state=="Sergipe" && city=="Nossa Senhora Aprecido") city<-"Nossa Senhora Aparecida"
  else if (state=="Sergipe" && city=="Maçambara") city<-"Macambira"
  else if (state=="Sergipe" && city=="Riachao do dantas") city<-"Riachão do Dantas"
  else if (state=="Sergipe" && city=="Umbauba") city<-"Umbaúba"
  else if (state=="Tocantins" && city=="Couto Magalhaes") city<-"Couto Magalhães"
  else if (state=="Tocantins" && city=="Dianopolis") city<-"Dianópolis"
  else if (state=="Tocantins" && city=="Lajedão") city<-"Lajeado"
  else if (state=="Tocantins" && city=="Tabocão") city<-"Fortaleza do Tabocão"
  else if (state=="Tocantins" && city=="Mosquito") city<-"Palmeiras do Tocantins"
  else if (state=="Tocantins" && city=="Paraná") city<-"Paranã"
  else if (state=="Tocantins" && city=="Pau d'Arco") city<-"Pau D'Arco"
  else if (state=="Tocantins" && city=="Ponte Alta do Norte") city<-"Ponte Alta do Tocantins"
  else if (state=="Tocantins" && city=="São Valério da Natividade") {city<-"São Valério"}
  
  # remove regions that were merged into other regions
  if ((state=="Mato Grosso do Sul" && city=="Sapezal") || 
      (state=="Minas Gerais" && city=="Itabirinha de Mantena") || 
      (state=="Paraná" && city=="Araújos") || 
      (state=="Santa Catarina" && city=="Floriniapolis") ||
      (state=="Minas Gerais" && city=="Chaveslandia")) {state=NA; city=NA; return(c(state, city))}
  
  # remove islands for which there is no epidemiological data
  if ((state=="Espírito Santo" && city=="Ilha Trindade") ||
      (state=="Espírito Santo" && city=="Ilhas de Martim Vaz") ||
      (state=="Rio Grande do Sul" && city=="Lagoa Mirim") ||
      (state=="São Paulo" && city=="Mangaratiba")) {state=NA; city=NA; return(c(state, city))}
  
  return(c(state, city))
}
corrected_regions <- bra2 %>%
  st_drop_geometry() %>%
  apply(1, function(x) fix_region_name(state=x[1], city=x[2])) %>%
  t() %>% 
  as.data.frame() %>%
  setNames(c("state", "city"))
bra2 <- bra2 %>%
  mutate(state=corrected_regions$state, city=corrected_regions$city) %>%
  na.omit()

#################################################################################################
# here we perform some tests on the modified GADM data set 

# check that spatial data is available for all regions in the DENV case data
bra_denv_case_data <- readRDS(file.path("bra2_cases/", "bra_denv_case_data.rds"))
count1 <- 0
for (ii in 2000:2016) {
  region_names <- bra_denv_case_data %>%
    filter(year==ii) %>%
    select(state, city) %>%
    unique()
  region_names <- region_names[which(!grepl("Município ignorado", region_names$city, fixed=TRUE)), ]
  no_match_rows <- region_names %>%
    anti_join(st_drop_geometry(bra2), by=c("state", "city"))
  if (nrow(no_match_rows)!=0) {
    count1 <- count1 + 1
    cat("The following regions do not have spatial data: \n")
    print(no_match_rows)
  }
}
if (count1==0) print("Spatial data available for all regions in the DENV case data.")

# check that there are no duplicate regions in the spatial data
count2 <- 0 
duplicate_rows <- bra2 %>%
  st_drop_geometry() %>%
  duplicated()
if (sum(duplicate_rows)!=0) {
  print("The following regions have duplicate entries: ")
  print(bra2[duplicate_rows, ])
  count2 <- count2 + 1
} else {
  print("Spatial data has no duplicate entries for a given region")
}

# check that all regions in the spatial data are represented in the DENV case data
count3 <- 0
for (ii in 2000:2016) {
  region_names <- bra_denv_case_data %>%
    filter(year==ii) %>%
    select(state, city) %>%
    unique()
  region_names <- region_names[which(!grepl("Município ignorado", region_names$city, fixed=TRUE)), ]
  no_match_rows <- st_drop_geometry(bra2) %>% 
    anti_join(region_names , by=c("state", "city"))
  if (nrow(no_match_rows)!=0) {
    count3 <- count3 + 1
    cat("The following regions are not in the DENV case data: \n")
    print(no_match_rows)
  }
}
if (count3==0) print("DENV case data contains entries (NA for no information) for all regions in spatial data.")

# check that no spatial polygon/multipolygon cover multiple regions
count4 <- 0 
# for (ii in 1:nrow(bra2)) {
#   region <- bra2[ii, ]
#   if ((region %>% st_drop_geometry() %>% as.character())[1] %in% checked) next
#   region_geometry <- st_make_valid(region$geometry)
#   region_centroid <- st_centroid(region_geometry)
#   if (!st_within(region_centroid, region_geometry, sparse=FALSE)) {
#     region_name <- region %>% st_drop_geometry() %>% as.character()
#     cat(paste0("Problem with the polygon/multipolygons of ", region_name[2], ", ", region_name[1], "\n"))
#     count4 <- count4 + 1
#   }
# }
# if (count4==0) print("No spatial polygons/multipolygons cover multiple regions in the spatial data set")

#################################################################################################
# here we compute some spatial statistics (e.g. area)
bra2 <- bra2 %>% 
  mutate(area_m2=as.numeric(st_area(st_make_valid(geometry)))) %>%
  select(state, city, area_m2)

#################################################################################################
# finally we save data
if ((count1+count2+count3+count4)==0) {
  print("All tests passed!")
  print("Saving the spatial polygon data...")
  saveRDS(bra2, "bra2_sf/bra2_sf_data.rds")
  simplified_bra2 <- st_simplify(st_make_valid(bra2), dTolerance=1000)
  saveRDS(simplified_bra2, "bra2_sf/simplified_bra2_sf_data.rds")
}
