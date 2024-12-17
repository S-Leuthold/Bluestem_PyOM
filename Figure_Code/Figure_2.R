### Script to build the figures that show the relative proportion of litter derived C and N in each fraction for the Bluestem litter
### Written by: Sam Leuthold (sam.leuthold@gmail.com)
### Last Updated: 12/14/2024

### Setup ----------------------------------------------------------------------

rm(list = ls())

require(pacman)

p_load(patchwork,
       lmerTest,
       emmeans,
       multcomp,
       tidyverse,
       conflicted)

select <- dplyr::select
filter <- dplyr::filter

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_Data.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_N_Distirbution.R"))

### Look at the data -----------------------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  group_by(Harvest_Num, Replicate) %>% 
  summarize(C_Remaining = sum(Carbon_g_m2), 
            N_Remaining = sum(Nitrogen_g_m2)) -> CN_Remaining

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  group_by(Harvest_Num, Fraction, Replicate) %>% 
  left_join(., 
            CN_Remaining) %>% 
  mutate(Proportion_C = (Carbon_g_m2/C_Remaining)*100,
         Proportion_N = (Nitrogen_g_m2/N_Remaining)*100) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarise(Mean_Proportion_C = mean(Proportion_C),
            SE_Proportion_C = sd(Proportion_C)/2,
            Mean_Proportion_N = mean(Proportion_N),
            SE_Proportion_N = sd(Proportion_N)/2) -> C_N_in_SOM


print(C_N_in_SOM, n = 21)


C_N_Mass_Dist %>% 
  filter(Harvest_Num %in% c(0, 12, 24, 36, 120)) -> C_N_Mass_Dist

## Statistics ------------------------------------------------------------------

### POM

C_N_Mass_Dist %>% 
  filter(Fraction == "lPOM",
         Harvest_Num != 0) %>% 
  lmer(Carbon_g_m2 ~ as.factor(Harvest_Num) + (1|Replicate), data = .) %>% 
  anova(.)

C_N_Mass_Dist %>% 
  filter(Fraction == "lPOM",
         Harvest_Num != 0) %>% 
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest_Num = "Tukey")) %>% 
  cld()

### CHAOM

C_N_Mass_Dist %>% 
  filter(Fraction == "CHAOM") %>% 
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  anova(.)

C_N_Mass_Dist %>% 
  filter(Fraction == "CHAOM",
         Harvest_Num != 0) %>% 
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest_Num = "Tukey")) %>% 
  cld()

### MAOM

C_N_Mass_Dist %>% 
  filter(Fraction == "MAOM",
         Harvest_Num != 0) %>% 
  lmer(Carbon_g_m2 ~ as.factor(Harvest_Num) + (1|Harvest_Num/Replicate), data = .) %>% 
  anova(.)

C_N_Mass_Dist %>% 
  filter(Fraction == "MAOM",
         Harvest_Num != 0) %>%
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest_Num = "Tukey")) %>% 
  cld()

## Create dataframe for errorbars ----------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(C_xmin = mean(Carbon_g_m2) - sd(Carbon_g_m2)/2,
            C_xmax = mean(Carbon_g_m2) + sd(Carbon_g_m2)/2,
            N_xmin = mean(Nitrogen_g_m2) - sd(Nitrogen_g_m2)/2,
            N_xmax = mean(Nitrogen_g_m2) + sd(Nitrogen_g_m2)/2) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),,
         Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM", "Litter")) -> Error_Bars


## Create text labels for proportion in fraction -------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(Mean_C_Proportion = round((mean(Carbon_g_m2)/261.164) * 100, 1),
            Mean_N_Proportion = round((mean(Nitrogen_g_m2)/8.615) * 100, 1)) %>% 
  ungroup() %>% 
  bind_cols(.,
            Error_Bars %>% 
              ungroup() %>% 
              mutate(C_Pos = C_xmax + (0.15 * 30),
                     N_Pos = N_xmax + (0.15 * 2.5)) %>% 
              dplyr::select(C_Pos,
                            N_Pos)) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction))-> C_N_Proportion_Text 

### Create N Plot --------------------------------------------------------------

## Fiddle with the months label ::eyeroll::

Months_Label <- paste0(unique(C_N_Mass_Dist$Harvest_Num), " Months")

Months_Label <- Months_Label[1:6]

names(Months_Label) <- c(12, 24, 36, 120)


## Create C plot ---------------------------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(Mean_C_g_m2 = mean(Carbon_g_m2),
            Mean_N_g_m2 = mean(Nitrogen_g_m2)) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),
         Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")) %>% 
  ggplot() +
  geom_col(aes(x = Mean_C_g_m2,
               y = Fraction,
               fill = Fraction),
           color = "black",
           show.legend = FALSE) +
  scale_fill_manual(values = c("Litter" = "#FC9D9A",
                               "POM" = "#83AF9B", 
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD")) +
  geom_errorbarh(data = Error_Bars,
                 aes(y = Fraction,
                     xmin = C_xmin,
                     xmax = C_xmax),
                 height = 0.4) +
  geom_point(data = C_N_Mass_Dist %>% 
               filter(Fraction != "Lost",
                      Fraction != "Litter",
                      Harvest_Num != 0) %>% 
               mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                                           TRUE ~ Fraction),
                      Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")),
             aes(x = Carbon_g_m2,
                 y = Fraction,
                 fill = Fraction),
             shape = 21,
             alpha = 0.75,
             color = "black") +
  facet_wrap(Harvest_Num~.,
             ncol = 1,
             labeller = labeller(Harvest_Num = Months_Label)) +
  geom_text(data = C_N_Proportion_Text,
            aes(x = C_Pos,
                y = Fraction,
                label = paste0(Mean_C_Proportion, "%"))) +
  scale_color_manual(values = c("Litter" = "#FC9D9A",
                               "POM" = "#83AF9B", 
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD")) +
  scale_x_reverse(name = expression("Litter derived carbon (g C "*m^-2*")"),
                  limits = c(30, 0),
                  expand = c(0,0),
                  breaks = seq(0,30,5)) +
  labs(tag = "a.") +
  ggtitle("Litter Derived") +
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),    
        axis.text.y = element_blank(),
        axis.ticks.y.right = element_blank(),
        plot.tag = element_text(face = "bold", size = 18, color = "black"),
        axis.ticks.y.left = element_blank(),
        legend.position = "none",
        legend.justification = c(1,1),
        strip.text = element_text(hjust = 0, size = 12),
        strip.background = element_blank(),
        legend.background = element_rect(color = "black")) -> C_Stocks_Plot_Litter


## Create plot -----------------------------------------------------------------

C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(Mean_C_g_m2 = mean(Carbon_g_m2),
            Mean_N_g_m2 = mean(Nitrogen_g_m2)) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),
         Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")) %>% 
  ggplot() +
  geom_col(aes(x = Mean_N_g_m2,
               y = Fraction,
               fill = Fraction),
           color = "black",
           show.legend = FALSE) +
  scale_fill_manual(values = c("Litter" = "#FC9D9A",
                               "POM" = "#83AF9B", 
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD")) +
  geom_errorbarh(data = Error_Bars,
                 aes(y = Fraction,
                     xmin = N_xmin,
                     xmax = N_xmax),
                 height = 0.4) +
  geom_point(data = C_N_Mass_Dist %>% 
               filter(Fraction != "Lost",
                      Fraction != "Litter",
                      Harvest_Num != 0) %>% 
               mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                                           TRUE ~ Fraction),
                      Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")),
             aes(x = Nitrogen_g_m2,
                 y = Fraction,
                 fill = Fraction),
             shape = 21,
             alpha = 0.75,
             color = "black") +
  facet_wrap(Harvest_Num~.,
             ncol = 1) +
  geom_text(data = C_N_Proportion_Text,
            aes(x = N_Pos,
                y = Fraction,
                label = paste0(Mean_N_Proportion, "%"))) +
  scale_color_manual(values = c("Litter" = "#FC9D9A",
                                "POM" = "#83AF9B", 
                                "CHAOM" = "#C8C8A9",
                                "MAOM" = "#F9CDAD")) +
  scale_x_continuous(name = expression("Litter derived nitrogen (g N "*m^-2*")"),
                     limits = c(0, 2.5),
                      expand = c(0, 0)) +
  scale_y_discrete(position = "left") +
  ggtitle("") +
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black", 
                                   hjust = 0.5,
                                   margin = margin(r = 10)),
        axis.ticks.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.position = "none",
        legend.justification = c(1,1),
        strip.background = element_blank(),
        legend.background = element_rect(color = "black"),
        strip.text = element_text(color = "white")) -> N_Stocks_Plot_Litter

C_Stocks_Plot_Litter|N_Stocks_Plot_Litter -> Stocks_Plot

Stocks_Plot

### ----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_PyOM_Data.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_PyOM_C_N_Distribution.R"))

### Look at the data -----------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  group_by(Harvest_Num, Replicate) %>% 
  summarize(C_Remaining = sum(Carbon_g_m2), 
            N_Remaining = sum(Nitrogen_g_m2)) -> CN_Remaining

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter") %>% 
  group_by(Harvest_Num, Fraction, Replicate) %>% 
  left_join(., 
            CN_Remaining) %>% 
  mutate(Proportion_C = (Carbon_g_m2/C_Remaining)*100,
         Proportion_N = (Nitrogen_g_m2/N_Remaining)*100) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarise(Mean_Proportion_C = mean(Proportion_C),
            SE_Proportion_C = sd(Proportion_C)/2,
            Mean_Proportion_N = mean(Proportion_N),
            SE_Proportion_N = sd(Proportion_N)/2) -> C_N_in_SOM


print(C_N_in_SOM, n = 21)


C_N_Mass_Dist %>% 
  filter(Harvest_Num %in% c(0, 12, 24, 36, 120)) -> C_N_Mass_Dist

## Statistics ------------------------------------------------------------------

### POM

C_N_Mass_Dist %>% 
  filter(Fraction == "lPOM",
         Harvest_Num != 0) %>% 
  lmer(Carbon_g_m2 ~ as.factor(Harvest_Num) + (1|Replicate), data = .) %>% 
  anova(.)

C_N_Mass_Dist %>% 
  filter(Fraction == "lPOM",
         Harvest_Num != 0) %>% 
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest_Num = "Tukey")) %>% 
  cld()

### CHAOM

C_N_Mass_Dist %>% 
  filter(Fraction == "CHAOM") %>% 
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  anova(.)

C_N_Mass_Dist %>% 
  filter(Fraction == "CHAOM",
         Harvest_Num != 0) %>% 
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest_Num = "Tukey")) %>% 
  cld()

### MAOM

C_N_Mass_Dist %>% 
  filter(Fraction == "MAOM",
         Harvest_Num != 0) %>% 
  lmer(Carbon_g_m2 ~ as.factor(Harvest_Num) + (1|Replicate), data = .) %>% 
  anova(.)

C_N_Mass_Dist %>% 
  filter(Fraction == "MAOM",
         Harvest_Num != 0) %>%
  mutate(Harvest_Num = as.factor(Harvest_Num)) %>% 
  lmer(Carbon_g_m2 ~ Harvest_Num + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest_Num = "Tukey")) %>% 
  cld()

## Create dataframe for errorbars ----------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(C_xmin = mean(Carbon_g_m2) - sd(Carbon_g_m2)/2,
            C_xmax = mean(Carbon_g_m2) + sd(Carbon_g_m2)/2,
            N_xmin = mean(Nitrogen_g_m2) - sd(Nitrogen_g_m2)/2,
            N_xmax = mean(Nitrogen_g_m2) + sd(Nitrogen_g_m2)/2) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),
         Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "lPOM")) %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num)) -> Error_Bars


## Create text labels for proportion in fraction -------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(Mean_C_Proportion = round((mean(Carbon_g_m2)/104.2787) * 100, 1),
            Mean_N_Proportion = round((mean(Nitrogen_g_m2)/3.1768) * 100, 1)) %>% 
  ungroup() %>% 
  bind_cols(.,
            Error_Bars %>% 
              ungroup() %>% 
              mutate(C_Pos = C_xmax + (0.125 * 135),
                     N_Pos = N_xmax + (0.125 * 3.8)) %>% 
              dplyr::select(
                C_Pos,
                N_Pos)) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction)) %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num)) -> C_N_Proportion_Text 

### Create N Plot --------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num))  -> BlackC_C_N_Mass_Dist
## Fiddle with the months label ::eyeroll::

Months_Label <- paste0(unique(BlackC_C_N_Mass_Dist$Harvest_Num), " Months")

Months_Label <- Months_Label[c(2, 1, 3)]

names(Months_Label) <- c(4, 11, 101)

## Create C plot ---------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(Mean_C_g_m2 = mean(Carbon_g_m2),
            Mean_N_g_m2 = mean(Nitrogen_g_m2)) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),
         Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")) %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num)) %>% 
  ggplot() +
  geom_col(aes(x = Mean_C_g_m2,
               y = Fraction,
               fill = Fraction),
           color = "black",
           show.legend = FALSE) +
  scale_fill_manual(values = c("Litter" = "#FC9D9A",
                               "POM" = "#83AF9B", 
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD")) +
  geom_errorbarh(data = Error_Bars,
                 aes(y = Fraction,
                     xmin = C_xmin,
                     xmax = C_xmax),
                 height = 0.4) +
  geom_point(data = BlackC_C_N_Mass_Dist %>% 
               filter(Fraction != "Lost",
                      Fraction != "Litter",
                      Harvest_Num != 0) %>% 
               mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                                           TRUE ~ Fraction),
                      Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")),
             aes(x = Carbon_g_m2,
                 y = Fraction,
                 fill = Fraction),
             shape = 21,
             color = "black",
             alpha = 0.75) +
  facet_wrap(Harvest_Num~.,
             ncol = 1,
             labeller = labeller(Harvest_Num = Months_Label)) +
  geom_text(data = C_N_Proportion_Text,
            aes(x = C_Pos,
                y = Fraction,
                label = paste0(Mean_C_Proportion, "%"))) +
  scale_color_manual(values = c("Litter" = "#FC9D9A",
                                "POM" = "#83AF9B", 
                                "CHAOM" = "#C8C8A9",
                                "MAOM" = "#F9CDAD")) +
  scale_x_reverse(name = expression("PyOM derived carbon (g C "*m^-2*")"),
                  limits = c(145, 0),
                  expand = c(0,0),
                  breaks = seq(0, 120, 30)) +
  labs(tag = "b.") +
  ggtitle("PyOM Derived") +
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),    
        axis.ticks.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),
        plot.tag = element_text(face = "bold", size = 18, color = "black"),
        axis.text.y = element_blank(),
        legend.position = "none",
        legend.justification = c(1,1),
        strip.text = element_text(hjust = 0, size = 12),
        strip.background = element_blank(),
        legend.background = element_rect(color = "black")) -> C_Stocks_Plot_PyOM

C_Stocks_Plot_PyOM


## Create plot -----------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost",
         Fraction != "Litter",
         Harvest_Num != 0) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarize(Mean_C_g_m2 = mean(Carbon_g_m2),
            Mean_N_g_m2 = mean(Nitrogen_g_m2)) %>% 
  mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                              TRUE ~ Fraction),
         Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")) %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num)) %>% 
  ggplot() +
  geom_col(aes(x = Mean_N_g_m2,
               y = Fraction,
               fill = Fraction),
           color = "black",
           show.legend = FALSE) +
  scale_fill_manual(values = c("Litter" = "#FC9D9A",
                               "POM" = "#83AF9B", 
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD")) +
  geom_errorbarh(data = Error_Bars,
                 aes(y = Fraction,
                     xmin = N_xmin,
                     xmax = N_xmax),
                 height = 0.4) +
  geom_point(data = BlackC_C_N_Mass_Dist %>% 
               filter(Fraction != "Lost",
                      Fraction != "Litter",
                      Harvest_Num != 0) %>% 
               mutate(Fraction = case_when(Fraction == "lPOM" ~ "POM",
                                           TRUE ~ Fraction),
                      Fraction = fct_relevel(Fraction, "MAOM", "CHAOM", "POM")),
             aes(x = Nitrogen_g_m2,
                 y = Fraction,
                 fill = Fraction),
             shape = 21,
             color = "black",
             alpha = 0.75) +
  facet_wrap(Harvest_Num~.,
             ncol = 1) +
  geom_text(data = C_N_Proportion_Text,
            aes(x = N_Pos,
                y = Fraction,
                label = paste0(Mean_N_Proportion, "%"))) +
  scale_color_manual(values = c("Litter" = "#FC9D9A",
                                "POM" = "#83AF9B", 
                                "CHAOM" = "#C8C8A9",
                                "MAOM" = "#F9CDAD")) +
  scale_x_continuous(name = expression("PyOM derived nitrogen (g N "*m^-2*")"),
                     limits = c(0, 3.8),
                     expand = c(0, 0)) +
  scale_y_discrete(position = "left") +
  ggtitle("") +
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black", 
                                   hjust = 0.5,
                                   margin = margin(r = 10)),
        axis.ticks.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.position = "none",
        legend.justification = c(1,1),
        strip.background = element_blank(),
        legend.background = element_rect(color = "black"),
        strip.text = element_text(color = "white")) -> N_Stocks_Plot_PyOM

PyOM_Stocks_Plot <- C_Stocks_Plot_PyOM + N_Stocks_Plot_PyOM

PyOM_Stocks_Plot

Stocks_Plot/PyOM_Stocks_Plot + plot_layout(nrow = 2,
                                           heights = c(4,3)) -> Figure_2
###

ggsave(filename = "./Final_Code/Plots/Figure_2.pdf",
       plot = Figure_2,
       device = "pdf",
       width = 90,
       height = 225,
       scale = 1.75,
       units = "mm",
       dpi = "print")


