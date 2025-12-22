library(data.table)
sediment <- fread("data/sediment_ages_v1.0_varves.csv")
colnames(sediment)
library(dplyr)
sediment %>%
  filter(sitename == "LakeAllie") %>%
  str()
sediment[!startsWith(id, "bchron_draw")][]
library(dplyr)

sediment %>%
  filter(!grepl("^(bchron_draw|bacon_draw)", id)) %>%
  group_by(sitename) %>%
  summarise(
    n_ids = n_distinct(id),
    ids = paste(unique(id), collapse = ", ")
  ) %>%
  arrange(desc(n_ids))

sediment %>%  pull(sitename) %>% unique()
keep_lakes <- c(
  "HostageLake",
  "FerryLake",
  "LilyLake",             # or "LilyLake(US:Minnesota)"
  "LoneStarLake",
  "HellHoleLake",
  "OzawindibLake",
  "MudLake" ,
  "PetersonLake",
  "RubyLake",
  "DarkLake" #or MudLake(US:Minnesota:Hubbard)
)

sediment %>%
  filter(sitename %in% keep_lakes) %>% 
  filter(!grepl("^(bchron_draw|bacon_draw)", id)) %>%
  distinct(sitename, id, lat, long)



library(dplyr)

lake_coords <- sediment %>%
  filter(sitename %in% keep_lakes & 
           !grepl("^(bchron_draw|bacon_draw)", id)) %>%
  select(sitename, lat, long, state, altitude, id) %>%
  distinct()


lake_coords <- sediment %>%
  filter(
    sitename %in% keep_lakes,
    !grepl("^(bchron_draw|bacon_draw)", id)
  ) %>%
  inner_join(
    tibble::tribble(
      ~sitename,            ~long_ref,    ~lat_ref,
      "FerryLake",         -92.125,      46.013,
      "LilyLake",          -92.273,      45.901,
      "HellHoleLake",      -92.21898,    45.786565,
      "LoneStarLake",      -92.365,      45.932,
      "DarkLakeLake",          -91.475,      45.275,
      "RubyLakeLake",          -91.458333,   45.283333,
      "LittlePineLake",   -91.483333,   45.283333,
      "HostageLake",            -94.1333,     46.55,
      "MudLake",                -94.75,       46.8667,
      "OzawindibLake",          -95.2697,     47.2313,
      "PetersonLake",           -95.3167,     46.9667
    ),
    by = "sitename"
  ) %>%
  filter(
    abs(lat  - lat_ref)  < 2,
    abs(long - long_ref) < 2
  ) %>%
  select(sitename, lat, long, state, altitude, id) %>%
  distinct() %>% head()
