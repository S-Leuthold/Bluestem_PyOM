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

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

###-----------------------------------------------------------------------------

T0_Char_Mass <- 5.92
T0_Char_Pct_C <- 55.31
T0_Char_C_Mass <- (T0_Char_Mass * (T0_Char_Pct_C/100))/0.0314
T0_Char_Pct_N <- 1.685
T0_Char_N_Mass <- (T0_Char_Mass * (T0_Char_Pct_N/100))/0.0314


BlackC_Time_Series_Data %>% 
  filter(Depth_cm == "Total",
         Cover_Treatment == "C") %>% 
  select(Harvest_Num,
         Fraction,
         Replicate,
         C_Mass,
         N_Mass) %>% 
  rename(Carbon_g_m2 = C_Mass,
         Nitrogen_g_m2 = N_Mass) %>% 
  pivot_wider(names_from = Fraction,
              values_from = c(Carbon_g_m2, Nitrogen_g_m2)) %>% 
  rowwise() %>% 
  mutate(Carbon_Lost = round(104.2787 - sum(Carbon_g_m2_lPOM, Carbon_g_m2_CHAOM, Carbon_g_m2_MAOM), 4),
         Nitrogen_Lost = round(3.176815 - sum(Nitrogen_g_m2_lPOM, Nitrogen_g_m2_CHAOM, Nitrogen_g_m2_MAOM), 4)) -> Intermediate

Intermediate %>% 
  ungroup() %>% 
  select(Harvest_Num,
         Replicate,
         Carbon_g_m2_lPOM,
         Carbon_g_m2_CHAOM,
         Carbon_g_m2_MAOM,
         Carbon_Lost) %>% 
  pivot_longer(cols = c(Carbon_g_m2_lPOM,
                        Carbon_g_m2_CHAOM,
                        Carbon_g_m2_MAOM,
                        Carbon_Lost), 
               names_to = "Fraction", 
               names_prefix = "Carbon_g_m2_",
               values_to = "Carbon_g_m2") %>% 
  mutate(Carbon_Mass_Prop = (Carbon_g_m2/104.2787)*100,
         Fraction = case_when(Fraction == "Carbon_Lost" ~ "Lost",
                              TRUE ~ Fraction),
         Fraction = as.factor(Fraction),
         Fraction = fct_relevel(Fraction, "Lost", "lPOM", "CHAOM", "MAOM")) %>% 
  left_join(., Intermediate %>% 
              ungroup() %>% 
              select(Harvest_Num,
                     Replicate,
                     Nitrogen_g_m2_lPOM,
                     Nitrogen_g_m2_CHAOM,
                     Nitrogen_g_m2_MAOM,
                     Nitrogen_Lost) %>% 
              pivot_longer(cols = c(Nitrogen_g_m2_lPOM,
                                    Nitrogen_g_m2_CHAOM,
                                    Nitrogen_g_m2_MAOM,
                                    Nitrogen_Lost), 
                           names_to = "Fraction", 
                           names_prefix = "Nitrogen_g_m2_",
                           values_to = "Nitrogen_g_m2") %>% 
              mutate(Nitrogen_Mass_Prop = (Nitrogen_g_m2/3.176815)*100,
                     Fraction = case_when(Fraction == "Nitrogen_Lost" ~ "Lost",
                                          TRUE ~ Fraction),
                     Fraction = as.factor(Fraction),
                     Fraction = fct_relevel(Fraction, "Lost", "lPOM", "CHAOM", "MAOM"))) -> BlackC_C_N_Mass_Dist



