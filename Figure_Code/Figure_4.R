### Script to calculate the CV of variation for the litter derived C and N for POM, CHAOM, and MAOM
### Written by: Sam Leuthold (sam.leuthold@gmail.com)
### Last Updated: 12/16/2024

### Setup ----------------------------------------------------------------------

rm(list = ls())

require(pacman)

p_load(tidyverse,
       patchwork,
       lmerTest,
       emmeans)

### ----------------------------------------------------------------------------

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_Data.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_N_Distirbution.R"))

### Look at the data -----------------------------------------------------------

Silly_Function <- function(X){
  y <- X * 100
  return(y)
}

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  filter(Harvest_Num != 0) %>% 
  group_by(Fraction, Replicate) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  group_by(Fraction) %>% 
  summarize(Mean_CV_Carbon = mean(CV_Carbon),
            SE_CV_Carbon = sd(CV_Carbon)/2,
            Mean_CV_Nitrogen = mean(CV_Nitrogen),
            SE_CV_Nitrogen = sd(CV_Nitrogen/2)) %>% 
  mutate(across(where(is.numeric), Silly_Function))


## ANOVA and pairwise comparisons ----------------------------------------------

# Carbon -----------------------------------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  filter(Harvest_Num != 0) %>% 
  group_by(Fraction, Replicate) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  ungroup() %>% 
  lmer(CV_Carbon ~ Fraction + (1|Replicate), data = .) %>% 
  emmeans(., pairwise ~ Fraction)

# Nitrogen ---------------------------------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  filter(Harvest_Num != 0) %>% 
  group_by(Fraction, Replicate) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  ungroup() %>% 
  lmer(CV_Nitrogen ~ Fraction + (1|Replicate), data = .) %>% 
  emmeans(., pairwise ~ Fraction)

## Create tibble of CLD values -------------------------------------------------

tibble(Fraction = rep(c("POM", "CHAOM", "MAOM"),2),
       Analyte = c(rep("Carbon", 3), rep("Nitrogen", 3)),
       Label = c("a", "a", "b",
                 "b", "a", "c"),
       Y = c(0.92, 1.05, 0.40,
             0.694, 0.823, 0.419)) %>% 
  mutate(Fraction = fct_relevel(Fraction, "POM", "CHAOM", "MAOM")) -> CV_Letters

## Create plot -----------------------------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  filter(Harvest_Num != 0) %>% 
  group_by(Fraction, Replicate) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  pivot_longer(cols = c(CV_Carbon, 
                        CV_Nitrogen),
               names_to = "Analyte",
               values_to = "CV",
               names_prefix = "CV_") %>% 
  group_by(Fraction, Analyte) %>% 
  summarize(Mean_CV = mean(CV),
            SE_CV = sd(CV)/2) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),
         Fraction = as.factor(Fraction),
         Fraction = fct_relevel(Fraction, "POM", "CHAOM", "MAOM")) %>% 
  ggplot() +
  geom_col(aes(x = Fraction,
               y = Mean_CV,
               fill = Analyte),
           color = "black",
           position = "dodge") +
  geom_errorbar(aes(x = Fraction,
                    ymin = Mean_CV - SE_CV,
                    ymax = Mean_CV + SE_CV,
                    group = Analyte),
                width = 0.5,
                position = position_dodge(width = 0.9)) +
  scale_y_continuous(name = "Temporal CV (%)",
                     breaks = seq(0, 1.25, 0.25),
                     labels = seq(0, 125, 25),
                     limits = c(0, 1.10),
                     expand = c(0,0)) +
  geom_text(data = CV_Letters,
            aes(x = Fraction,
                y = Y,
                label = Label,
                group = Analyte),
            position = position_dodge(width = 0.9),
            size = 8) +
  scale_fill_manual(name = "Analyte",
                    values = c("Carbon" = "#C8C8A9",
                               "Nitrogen" = "#FC9D9A")) +
  ggtitle("Litter Derived SOM") +
  labs(tag = "a.") +
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none",
        title = element_text(size = 15, face = "bold")) -> Bluestem_Litter_CV

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_PyOM_Data.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_PyOM_C_N_Distribution.R"))

### Look at the data -----------------------------------------------------------

Silly_Function <- function(X){
  y <- X * 100
  return(y)
}

BlackC_C_N_Mass_Dist %>% 
  group_by(Fraction, Replicate) %>% 
  filter(Fraction != "Lost",
         Carbon_Mass_Prop < 100) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2, na.rm = TRUE)/mean(Carbon_g_m2, na.rm = TRUE),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  group_by(Fraction) %>% 
  summarize(Mean_CV_Carbon = mean(CV_Carbon),
            SE_CV_Carbon = sd(CV_Carbon)/2,
            Mean_CV_Nitrogen = mean(CV_Nitrogen),
            SE_CV_Nitrogen = sd(CV_Nitrogen/2)) %>% 
  mutate(across(where(is.numeric), Silly_Function))


## ANOVA and pairwise comparisons ----------------------------------------------

# Carbon -----------------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  group_by(Fraction, Replicate) %>% 
  filter(Fraction != "Lost",
         Carbon_Mass_Prop < 100) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  ungroup() %>% 
  lmer(CV_Carbon ~ Fraction + (1|Replicate), data = .) %>% 
  emmeans(., pairwise ~ Fraction)

# Nitrogen ---------------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  filter(Harvest_Num != 0) %>% 
  group_by(Fraction, Replicate) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  ungroup() %>% 
  lmer(CV_Nitrogen ~ Fraction + (1|Replicate), data = .) %>% 
  emmeans(., pairwise ~ Fraction)

## Create tibble of CLD values -------------------------------------------------

tibble(Fraction = rep(c("POM", "CHAOM", "MAOM"),2),
       Analyte = c(rep("Carbon", 3), rep("Nitrogen", 3)),
       Label = c("a", "a", "a",
                 "a", "ab", "b"),
       Y = c(0.348 + 0.03, 0.779 + 0.03, 0.563 + 0.03,
             0.307 + 0.03, 0.551 + 0.03, 0.570 + 0.03)) %>% 
  mutate(Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")) -> CV_Letters_PyOM

## Create plot -----------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  filter(Harvest_Num != 0) %>% 
  group_by(Fraction, Replicate) %>% 
  summarise(CV_Carbon = sd(Carbon_g_m2)/mean(Carbon_g_m2),
            CV_Nitrogen = sd(Nitrogen_g_m2)/mean(Nitrogen_g_m2)) %>% 
  pivot_longer(cols = c(CV_Carbon, 
                        CV_Nitrogen),
               names_to = "Analyte",
               values_to = "CV",
               names_prefix = "CV_") %>% 
  group_by(Fraction, Analyte) %>% 
  summarize(Mean_CV = mean(CV),
            SE_CV = sd(CV)/2) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
         TRUE ~ Fraction),
         Fraction = as.factor(Fraction),
         Fraction = fct_relevel(Fraction, "POM", "CHAOM", "MAOM")) %>% 
  ggplot() +
  geom_col(aes(x = Fraction,
               y = Mean_CV,
               fill = Analyte),
           color = "black",
           position = "dodge") +
  geom_errorbar(aes(x = Fraction,
                    ymin = Mean_CV - SE_CV,
                    ymax = Mean_CV + SE_CV,
                    group = Analyte),
                width = 0.5,
                position = position_dodge(width = 0.9)) +
  scale_y_continuous(name = "Temporal CV (%)",
                     breaks = seq(0, 1.25, 0.25),
                     labels = seq(0, 125, 25),
                     limits = c(0, 1.10),
                     expand = c(0,0)) +
  geom_text(data = CV_Letters_PyOM,
            aes(x = Fraction,
                y = Y,
                label = Label,
                group = Analyte),
            position = position_dodge(width = 0.9),
            size = 8) +
  scale_fill_manual(name = "Analyte",
                    values = c("Carbon" = "#C8C8A9",
                               "Nitrogen" = "#FC9D9A")) +
  ggtitle("PyOM Derived SOM") +
  labs(tag = "b.") +
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.title.y = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1,1),
        legend.text = element_text(size = 15),
        legend.background = element_rect(color = "black"),
        title = element_text(size = 15, face = "bold")) -> PyOM_CV

## Combine with litter CV Plot

Temporal_CV <- Bluestem_Litter_CV + PyOM_CV

Temporal_CV

## Save the plot to file -------------------------------------------------------

ggsave(filename = "./Final_Code/Plots/Figure_4.pdf",
       plot = Temporal_CV,
       device = "pdf",
       width = 180,
       height = 90,
       units = "mm",
       scale = 2,
       dpi = "print")

