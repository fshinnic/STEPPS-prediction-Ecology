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
  select(id, sitename, lat, long, state, altitude, age_bacon, age_bchron,depth , bchron_draw1, bacon_draw1, ASH, BEECH ,BIRCH  ,
         ELM, HEMLOCK, MAPLE,   OAK, OTHER.CONIFER, OTHER.HARDWOOD,  PINE ,SPRUCE ,TAMARACK) %>%  # optional: drop draw columns
  pull(sitename) %>% 
  unique

keep_lakes <- c(
  "HostageLake",
  "FerryLake",
  "LilyLake",             # or "LilyLake(US:Minnesota)"
  "LoneStarLake",
  "HellHoleLake",
  "OzawindibLake",
  "MudLake"               # or "MudLake(US:Minnesota:Hubbard)"
)

library(dplyr)

lake_coords <- sediment %>%
  filter(sitename %in% keep_lakes & 
           !grepl("^(bchron_draw|bacon_draw)", id)) %>%
  select(sitename, lat, long, state, altitude) %>%
  distinct()

lake_coords

# correct information:
Factor	Ferry Lake	Lily Lake	Hellhole Lake	Lonestar Lake	Dark Lake	Ruby Lake	Little Pine Lake	Hostage	Mud	Ozawindib	Peterson
Longitude	-92.125	-92.273	-92.21898	-92.365	-91.475	-91.458333	-91.483333	-94.1333	-94.75	-95.2697	-95.3167
Lattitude	46.013	45.901	45.786565	45.932	45.275	45.283333	45.283333	46.55	46.8667	47.2313	46.9667
