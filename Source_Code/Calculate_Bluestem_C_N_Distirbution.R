### Source script for calculating the distribution of C and N among litter and SOM fractions and the proportion of C and N.
### Adapted from Bluestem_C_Distribution.R
### Author: Sam Leuthold (sam.leuthold@gmail.com)
### Last update: 12/12/2024

###-----------------------------------------------------------------------------

require(pacman)

p_load(tidyverse,
       patchwork,
       ggpubr,
       gdata,
       conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

###-----------------------------------------------------------------------------

## Create a tibble with the litter C and N stocks ------------------------------

tibble(Harvest_Num = rep(c(0, 7, 12, 18, 24, 36, 120), 4),
       Fraction = "Litter",
       Replicate = rep(seq(1, 4), each = 7),
       Carbon_g_m2 = rep(c(8.2/0.0314, 225.8, 111.2, 89.2, 55.6, 3.4, 0), 4),
       Nitrogen_g_m2 = rep(c(0.27048/0.0314, 4.5, 3.4, 3.0, 2.0, 0.1, 0), 4)) -> Applied_Litter_Data

## Pull the time series data in and combine ------------------------------------

Litter_Time_Series_Data %>% 
  filter(Depth_cm == "Total",
         Fumigation_Treatment == "F") %>% 
  select(Harvest_Num,
         Fraction,
         Replicate,
         C_Mass,
         N_Mass) %>%  
  rename(Carbon_g_m2 = C_Mass,
         Nitrogen_g_m2 = N_Mass) %>% 
  bind_rows(tibble(Harvest_Num = rep(0, 12),
                   Replicate = rep(seq(1,4),3),
                   Fraction = rep(c("MAOM", "lPOM", "CHAOM"), each = 4),
                   Carbon_g_m2 = rep(0, 12),
                   Nitrogen_g_m2 = rep(0, 12))) %>% 
  bind_rows(Applied_Litter_Data) %>% 
  pivot_wider(names_from = Fraction,
              values_from = c(Carbon_g_m2, Nitrogen_g_m2)) %>% 
## Sum across rows for each time step ------------------------------------------
  rowwise() %>% 
  mutate(Carbon_Lost = round(261.146 - sum(Carbon_g_m2_lPOM, Carbon_g_m2_CHAOM, Carbon_g_m2_MAOM, Carbon_g_m2_Litter), 4),
         Nitrogen_Lost = round(8.614013 - sum(Nitrogen_g_m2_lPOM, Nitrogen_g_m2_CHAOM, Nitrogen_g_m2_MAOM, Nitrogen_g_m2_Litter), 4)) -> Intermediate 

## Pull the C out, and the rejoin it with the N to make a long tibble-----------

Intermediate %>% 
  ungroup() %>% 
  select(Harvest_Num,
         Replicate,
         Carbon_g_m2_lPOM,
         Carbon_g_m2_CHAOM,
         Carbon_g_m2_MAOM,
         Carbon_g_m2_Litter,
         Carbon_Lost) %>% 
  pivot_longer(cols = c(Carbon_g_m2_lPOM,
                        Carbon_g_m2_CHAOM,
                        Carbon_g_m2_MAOM,
                        Carbon_g_m2_Litter,
                        Carbon_Lost), 
               names_to = "Fraction", 
               names_prefix = "Carbon_g_m2_",
               values_to = "Carbon_g_m2") %>% 
  mutate(Carbon_Mass_Prop = (Carbon_g_m2/261.146)*100,
         Fraction = case_when(Fraction == "Carbon_Lost" ~ "Lost",
                              TRUE ~ Fraction),
         Fraction = as.factor(Fraction),
         Fraction = fct_relevel(Fraction, "Lost", "Litter", "lPOM", "CHAOM", "MAOM")) %>% 
  left_join(., Intermediate %>% 
              ungroup() %>% 
              select(Harvest_Num,
                     Replicate,
                     Nitrogen_g_m2_lPOM,
                     Nitrogen_g_m2_CHAOM,
                     Nitrogen_g_m2_MAOM,
                     Nitrogen_g_m2_Litter,
                     Nitrogen_Lost) %>% 
              pivot_longer(cols = c(Nitrogen_g_m2_lPOM,
                                    Nitrogen_g_m2_CHAOM,
                                    Nitrogen_g_m2_MAOM,
                                    Nitrogen_g_m2_Litter,
                                    Nitrogen_Lost), 
                           names_to = "Fraction", 
                           names_prefix = "Nitrogen_g_m2_",
                           values_to = "Nitrogen_g_m2") %>% 
              mutate(Nitrogen_Mass_Prop = (Nitrogen_g_m2/8.614013)*100,
                     Fraction = case_when(Fraction == "Nitrogen_Lost" ~ "Lost",
                                          TRUE ~ Fraction),
                     Fraction = as.factor(Fraction),
                     Fraction = fct_relevel(Fraction, "Lost", "Litter", "lPOM", "CHAOM", "MAOM"))) -> C_N_Mass_Dist


