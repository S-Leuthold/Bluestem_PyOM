### Script to build the figures that show the CN ratios across time of the three different SOM fractions
### Written by: Sam Leuthold (sam.leuthold@gmail.com)
### Last Updated: 12/12/2024

### Setup ----------------------------------------------------------------------

rm(list = ls())

require(pacman)

p_load(patchwork,
       lmerTest,
       emmeans,
       multcomp,
       tidyverse,
       ggbreak,
       conflicted)

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_Data.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_N_Distirbution.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_PyOM_Data.R"))

suppressMessages(source("./Final_Code/Source_Scripts/Calculate_PyOM_C_N_Distribution.R"))

### Litter data ----------------------------------------------------------------

C_N_Mass_Dist %>%
  filter(Harvest_Num %in% c(12, 24, 36, 120),
         Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>% 
  mutate(Time_Group = case_when(Harvest_Num < 100 ~ "Early",
                                TRUE ~ "Late")) %>% 
  group_by(Harvest_Num, Fraction, Time_Group) %>% 
  summarise(CN = mean(Carbon_g_m2/Nitrogen_g_m2)) %>% 
  ggplot(aes(x = Harvest_Num,
             y = CN)) +
  geom_errorbar(data = C_N_Mass_Dist %>% 
                  filter(Harvest_Num %in% c(12, 24, 36, 120),
                         Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>%  
                  mutate(Time_Group = case_when(Harvest_Num < 100 ~ "Early",
                                                TRUE ~ "Late")) %>% 
                  group_by(Harvest_Num, Fraction, Time_Group) %>% 
                  summarise(CN_Mean = mean(Carbon_g_m2/Nitrogen_g_m2),
                            CN_SE = sd(Carbon_g_m2/Nitrogen_g_m2)/2) %>% 
                  mutate(CN_Upper = CN_Mean + CN_SE,
                         CN_Lower = CN_Mean - CN_SE),
                aes(x = Harvest_Num, 
                    ymin = CN_Lower,
                    ymax = CN_Upper),
                width = 0.5,
                inherit.aes = F) +
  geom_line(aes(group = Fraction)) +
  geom_point(aes(fill = Fraction,
                 shape = Fraction),
             color = "black", 
             size = 7) +
  scale_fill_manual(name = "SOM Fraction",
                    values = c("lPOM" = "#83AF9B",
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD"),
                    labels = c("lPOM" = "POM",
                               "CHAOM" = "CHAOM",
                               "MAOM" = "MAOM")) +
  scale_shape_manual(name = "SOM Fraction",
                     values = c("lPOM" = 21,
                                "CHAOM" = 22,
                                "MAOM" = 23),
                     labels = c("lPOM" = "POM",
                                "CHAOM" = "CHAOM",
                                "MAOM" = "MAOM")) +
  coord_cartesian(xlim = c(10,60),
                  ylim = c(5, 20),
                  clip = "on") +
  ylab("Litter Derived C:N Ratio") +
  scale_x_continuous(name = "Harvest Month",
                     breaks = c(12, 24, 36)) +
  theme(legend.position = c(0.01,0.01),
        legend.justification = c(0,0),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 20, color = "black", face  = "bold"),
        axis.text = element_text(size = 17, color = "black"),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.background = element_rect(fill = "white"))  -> Litter_0_36

C_N_Mass_Dist %>%
  filter(Harvest_Num %in% c(12, 24, 36, 120),
         Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>% 
  mutate(Time_Group = case_when(Harvest_Num < 100 ~ "Early",
                                TRUE ~ "Late")) %>% 
  group_by(Harvest_Num, Fraction, Time_Group) %>% 
  summarise(CN = mean(Carbon_g_m2/Nitrogen_g_m2)) %>% 
  ggplot(aes(x = Harvest_Num,
             y = CN)) +
  geom_errorbar(data = C_N_Mass_Dist %>% 
                  filter(Harvest_Num %in% c(12, 24, 36, 120),
                         Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>%  
                  mutate(Time_Group = case_when(Harvest_Num < 100 ~ "Early",
                                                TRUE ~ "Late")) %>% 
                  group_by(Harvest_Num, Fraction, Time_Group) %>% 
                  summarise(CN_Mean = mean(Carbon_g_m2/Nitrogen_g_m2),
                            CN_SE = sd(Carbon_g_m2/Nitrogen_g_m2)/2) %>% 
                  mutate(CN_Upper = CN_Mean + CN_SE,
                         CN_Lower = CN_Mean - CN_SE),
                aes(x = Harvest_Num, 
                    ymin = CN_Lower,
                    ymax = CN_Upper),
                width = 0.5,
                inherit.aes = F) +
  geom_line(aes(group = Fraction)) +
  geom_point(aes(fill = Fraction,
                 shape = Fraction),
             color = "black", 
             size = 7) +
  scale_fill_manual(name = "SOM Fraction",
                    values = c("lPOM" = "#83AF9B",
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD"),
                    labels = c("lPOM" = "POM",
                               "CHAOM" = "CHAOM",
                               "MAOM" = "MAOM")) +
  scale_shape_manual(name = "SOM Fraction",
                     values = c("lPOM" = 21,
                                "CHAOM" = 22,
                                "MAOM" = 23),
                     labels = c("lPOM" = "POM",
                                "CHAOM" = "CHAOM",
                                "MAOM" = "MAOM")) +
  coord_cartesian(xlim = c(100, 125),
                  ylim = c(5, 20),
                  clip = "on") +
  scale_x_continuous(name = "",
                     breaks = c(108, 120),
                     expand = c(0, 0)) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        axis.line.y.left = element_line(color = NULL),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 17, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> Litter_100_120

ggplot() +
  scale_y_continuous(limits = c(5,20),
                     expand = c(0,0)) +
  scale_x_continuous(
                     expand = c(0,0)) +
  coord_cartesian(xlim = c(0,1),
                  clip = "off") +
  annotate(geom = "text",
           x = 1.8,
           y = 5.2,
           label = "/",
           size = 10) +
  annotate(geom = "text",
           x = -0.6,
           y = 5.2,
           label = "/",
           size = 10) +
  annotate(geom = "text",
           x = 1.8,
           y = 19.95,
           label = "/",
           size = 10) +
  annotate(geom = "text",
           x = -0.6,
           y = 19.95,
           label = "/",
           size = 10) +
  annotate(geom = "segment",
           x = -0.6,
           xend = 1.8,
           y = 15,
           yend = 14.2,
           lty = 2) +
  annotate(geom = "segment",
           x = -0.6,
           xend = 1.8,
           y = 14.0,
           yend = 12.5,
           lty = 2) +
  annotate(geom = "segment",
           x = -0.6,
           xend = 1.8,
           y = 10.0,
           yend = 9.35,
           lty = 2) +
  theme(panel.background = element_rect(fill = "white", color = NULL),
        plot.background = element_rect(fill = "white", color = NULL),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> Intermediate_Patch_Litter

Litter_0_36 + inset_element(
  Litter_100_120,
  left = 3/5,
  bottom = -0.075,
  right = 1.01,
  top = 1.01
) +
  inset_element(
    Intermediate_Patch_Litter,
    left = 3/5 - 0.02,
    bottom = -0.025,
    right = 3/5 + 0.02,
    top = 1.01
  ) -> Litter_CN_Plot
  
Litter_CN_Plot

### PyOM -----------------------------------------------------------------------

BlackC_C_N_Mass_Dist %>%
  filter(Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarise(CN = mean(Carbon_g_m2/Nitrogen_g_m2)) %>% 
  ggplot(aes(x = Harvest_Num,
             y = CN)) +
  geom_errorbar(data = BlackC_C_N_Mass_Dist %>% 
                  filter(Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>%  
                  group_by(Harvest_Num, Fraction) %>% 
                  summarise(CN_Mean = mean(Carbon_g_m2/Nitrogen_g_m2),
                            CN_SE = sd(Carbon_g_m2/Nitrogen_g_m2)/2) %>% # Need to check this SE calculation...
                  mutate(CN_Upper = CN_Mean + CN_SE,
                         CN_Lower = CN_Mean - CN_SE),
                aes(x = Harvest_Num, 
                    y = CN_Mean,
                    ymin = CN_Lower,
                    ymax = CN_Upper),
                width = 0.5,
                inherit.aes = F) +
  geom_line(aes(group = Fraction)) +
  geom_point(aes(fill = Fraction,
                 shape = Fraction),
             color = "black", 
             size = 7) +
  scale_fill_manual(name = "SOM Fraction",
                    values = c("lPOM" = "#83AF9B",
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD"),
                    labels = c("lPOM" = "POM",
                               "CHAOM" = "CHAOM",
                               "MAOM" = "MAOM")) +
  scale_shape_manual(name = "SOM Fraction",
                     values = c("lPOM" = 21,
                                "CHAOM" = 22,
                                "MAOM" = 23),
                     labels = c("lPOM" = "POM",
                                "CHAOM" = "CHAOM",
                                "MAOM" = "MAOM")) +
  coord_cartesian(xlim = c(2, 20),
                  ylim = c(10, 40),
                  clip = "on") +
  ylab("PyOM Derived C:N Ratio") +
  scale_x_continuous(name = "Harvest Month",
                     breaks = c(4, 11)) +
  theme(legend.position = "none",
        legend.justification = c(0,0),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 20, color = "black", face = "bold"),
        axis.text = element_text(size = 17, color = "black"),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.background = element_rect(fill = "white"))  -> PyOM_0_12

PyOM_0_12

###-----------------------------------------------------------------------------

BlackC_C_N_Mass_Dist %>%
  filter(Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>% 
  group_by(Harvest_Num, Fraction) %>% 
  summarise(CN = mean(Carbon_g_m2/Nitrogen_g_m2)) %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num)) %>% 
  ggplot(aes(x = Harvest_Num,
             y = CN)) +
  geom_errorbar(data = BlackC_C_N_Mass_Dist %>% 
                  filter(Fraction %in% c("lPOM", "MAOM", "CHAOM")) %>%  
                  group_by(Harvest_Num, Fraction) %>% 
                  summarise(CN_Mean = mean(Carbon_g_m2/Nitrogen_g_m2),
                            CN_SE = sd(Carbon_g_m2/Nitrogen_g_m2)/2) %>% # Need to check this SE calculation...
                  mutate(CN_Upper = CN_Mean + CN_SE,
                         CN_Lower = CN_Mean - CN_SE),
                aes(x = Harvest_Num, 
                    ymin = CN_Lower,
                    ymax = CN_Upper),
                width = 0.5,
                inherit.aes = F) +
  geom_line(aes(group = Fraction)) +
  geom_point(aes(fill = Fraction,
                 shape = Fraction),
             color = "black", 
             size = 7) +
  scale_fill_manual(name = "SOM Fraction",
                    values = c("lPOM" = "#83AF9B",
                               "CHAOM" = "#C8C8A9",
                               "MAOM" = "#F9CDAD"),
                    labels = c("lPOM" = "POM",
                               "CHAOM" = "CHAOM",
                               "MAOM" = "MAOM")) +
  scale_shape_manual(name = "SOM Fraction",
                     values = c("lPOM" = 21,
                                "CHAOM" = 22,
                                "MAOM" = 23),
                     labels = c("lPOM" = "POM",
                                "CHAOM" = "CHAOM",
                                "MAOM" = "MAOM")) +
  coord_cartesian(xlim = c(99, 103),
                  ylim = c(10, 40),
                  clip = "on") +
  scale_x_continuous(name = "",
                     breaks = c(101),
                     expand = c(0, 0)) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        axis.line.y.left = element_line(color = NULL),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 17, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> PyOM_72_84

PyOM_72_84

ggplot() +
  scale_y_continuous(limits = c(10,40),
                     expand = c(0,0)) +
  scale_x_continuous(
    expand = c(0,0)) +
  coord_cartesian(xlim = c(0,1),
                  clip = "off") +
  annotate(geom = "text",
           x = 1.8,
           y = 10.2,
           label = "/",
           size = 10) +
  annotate(geom = "text",
           x = -0.6,
           y = 10.2,
           label = "/",
           size = 10) +
  annotate(geom = "text",
           x = 1.8,
           y = 39.95,
           label = "/",
           size = 10) +
  annotate(geom = "text",
           x = -0.6,
           y = 39.95,
           label = "/",
           size = 10) +
  annotate(geom = "segment",
           x = -0.6,
           xend = 1.78,
           y = 28.70,
           yend = 31.35
           ,
           lty = 2) +
  annotate(geom = "segment",
           x = -0.6,
           xend = 1.8,
           y = 24.10,
           yend = 21.55,
           lty = 2) +
  annotate(geom = "segment",
           x = -0.6,
           xend = 1.8,
           y = 21.50,
           yend = 24.30,
           lty = 2) +
  theme(panel.background = element_rect(fill = "white", color = NULL),
        plot.background = element_rect(fill = "white", color = NULL),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> Intermediate_Patch_PyOM

PyOM_0_12 + inset_element(
  PyOM_72_84,
  left = 2/3,
  bottom = -0.075,
  right = 1.01,
  top = 1.01
) +
  inset_element(
    Intermediate_Patch_PyOM,
    left = 2/3 - 0.02,
    bottom = -0.025,
    right = 2/3 + 0.02,
    top = 1.01
  ) -> PyOM_CN_Plot

PyOM_CN_Plot

Litter_CN_Plot / PyOM_CN_Plot -> CN_Figure

ggsave(filename = "./Final_Code/Plots/Figure_3.pdf",
       plot = CN_Figure,
       device = "pdf",
       units = "px",
       height = 2366,
       width = 1233,
       scale = 2,
       dpi = "print")


