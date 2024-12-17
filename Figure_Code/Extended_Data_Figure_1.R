### Script for evidence of migration of PyOM (Ext. Data Figure 1)
### Written by: Sam Leuthold (sam.leuthold@gmail.com)
### Last Updated: 12/16/2024

### Set-up ---------------------------------------------------------------------

library(pacman)

p_load("tidyverse")

### ----------------------------------------------------------------------------

source("./Final_Code/Source_Scripts/Calculate_PyOM_Data.R")
source("./Final_Code/Source_Scripts/Calculate_PyOM_C_N_Distribution.R")

### ----------------------------------------------------------------------------

T0_Char_Pct_C <- 55.31
T0_Char_Mass <- 5.92
T0_Char_C_Mass <- T0_Char_Mass * (T0_Char_Pct_C/100)

### ----------------------------------------------------------------------------

read_csv("./Final_Code/Data/Leuthold et al., 2025 - PyOM Original Data.csv") %>% 
  filter(TIME == "11mts",
         COVER == "C",
         FRACTION == "bulk") %>% 
  select(REP,
         TIME,
         DEPTH,
         Total_derived_mg_C) %>% 
  rename(Replicate = REP,
         Harvest = TIME,
         Depth_cm = DEPTH) %>% 
# Rename the variables -------------------------------------------------------
  mutate(Replicate = as.numeric(str_sub(Replicate, 2, 2)),
         Derived_C_Mass = Total_derived_mg_C/1000) %>% 
  select(-Total_derived_mg_C) %>% 
  filter(Depth_cm %in% c("0-2", 
                         "2-5",
                         "5-10")) %>% 
  group_by(Depth_cm, Replicate) %>% 
## Get a proportion remaining value ------------------------------------------
  mutate(Substrate_Derived_C_Remaining = Derived_C_Mass /T0_Char_C_Mass * 100) %>%  
  group_by(Depth_cm) %>% 
  summarize(SE = sd(Substrate_Derived_C_Remaining)/sqrt(n()),
            mean = mean(Substrate_Derived_C_Remaining)) %>% 
  mutate(Depth_cm = fct_relevel(Depth_cm, "0-2", "2-5", "5-10"),
         Sampling = "11 months") %>% 
  rbind(.,
        BlackC_Derived_Data %>% 
          group_by(Depth_cm) %>% 
          summarize(SE = sd(Bulk_Initial_Char_C_Remaining)/sqrt(n()),
                    mean = mean(Bulk_Initial_Char_C_Remaining)) %>% 
          mutate(Sampling = "101 months")) %>% 
  mutate(Sampling = fct_relevel(Sampling, "11 months", "101 months")) %>% 
  ungroup() %>% 
  mutate(Test = c(3,2,1,6,5,4)) %>% 
  arrange(Test) %>% 
  group_by(Sampling) %>% 
  mutate(Position = cumsum(mean),
         Depth_cm = fct_relevel(Depth_cm, "0-2", "2-5", "5-10")) -> Bars_Data
  

read_csv("./Final_Code/Data/Leuthold et al., 2025 - PyOM Original Data.csv") %>% 
  filter(TIME == "11mts",
         COVER == "C",
         FRACTION == "bulk") %>% 
  select(REP,
         TIME,
         DEPTH,
         Total_derived_mg_C) %>% 
  rename(Rep = REP,
         Sampling = TIME,
         Depth_cm = DEPTH) %>% 
  mutate(Rep = as.numeric(str_sub(Rep, 2, 2)),
         Derived_C_Mass = Total_derived_mg_C/1000) %>% 
  select(-Total_derived_mg_C) %>% 
  filter(Depth_cm %in% c("0-2", 
                         "2-5",
                         "5-10")) %>% 
  group_by(Depth_cm, 
           Rep,
           Sampling) %>%
  summarise(Derived_C_Mass = mean(Derived_C_Mass, na.rm = TRUE)) %>% 
  ## Get a proportion remaining value ------------------------------------------
  mutate(Substrate_Derived_C_Remaining = Derived_C_Mass /T0_Char_C_Mass * 100) %>%  
  mutate(Depth_cm = as.factor(Depth_cm),
         Depth_cm = fct_relevel(Depth_cm, "0-2", "2-5", "5-10"),
         Sampling = "11 months") %>% 
  select(Rep,
         Depth_cm,
         Sampling,
         Substrate_Derived_C_Remaining) %>% 
  rbind(.,
        BlackC_Derived_Data %>% 
          select(Rep,
                 Depth_cm,
                 Bulk_Initial_Char_C_Remaining) %>%
          rename(Substrate_Derived_C_Remaining = Bulk_Initial_Char_C_Remaining) %>% 
          mutate(Sampling = "101 months") %>% 
          select(Rep, 
                 Depth_cm,
                 Sampling,
                 Substrate_Derived_C_Remaining)) %>% 
  ungroup() %>% 
  mutate(Sampling = as.factor(Sampling),
         Sampling = fct_relevel(Sampling, "11 months", "101 months"),
         Position = case_when(Sampling == "11 months" & Depth_cm == "0-2" ~ Substrate_Derived_C_Remaining + as.numeric(Bars_Data[2,6]),
                              Sampling == "11 months" & Depth_cm == "2-5" ~ Substrate_Derived_C_Remaining + as.numeric(Bars_Data[1,6]),
                              Sampling == "11 months" & Depth_cm == "5-10" ~ Substrate_Derived_C_Remaining + 0,
                              Sampling == "101 months" & Depth_cm == "0-2" ~ Substrate_Derived_C_Remaining + as.numeric(Bars_Data[5,6]),
                              Sampling == "101 months" & Depth_cm == "2-5" ~ Substrate_Derived_C_Remaining + as.numeric(Bars_Data[4,6]),
                              Sampling == "101 months" & Depth_cm == "5-10" ~ Substrate_Derived_C_Remaining + 0)) -> Point_Data



ggplot(data = Bars_Data) +
  geom_col(aes(x = Sampling, 
               y = mean,
               fill = Depth_cm),
           position = "stack",
           color = "black") +
  geom_errorbar(aes(x = Sampling,
                    ymin = Position - SE,
                    ymax = Position + SE),
                width = 0.2) +
  geom_point(data = Point_Data,
             aes(x = Sampling,
                 y = Position,
                 fill = Depth_cm),
             shape = 21,
             size = 4,
             show.legend = FALSE) +
  scale_y_continuous(name = "Proportion of initial PyOM C (%)",
                     limits = c(0, 82),
                     expand = c(0,0)) +
  scale_fill_manual(name = "Depth (cm)",
                    values = c("0-2" = "#83AF9B", 
                               "2-5" = "#C8C8A9",
                               "5-10" = "#F9CDAD")) +
  theme_void() +
  theme(panel.border = element_rect(fill = "transparent",
                                    color = "black",
                                    size = 1),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_text(color = "black",
                                    size = 15,
                                    angle = 90,
                                    face = "bold",
                                    margin = margin(r = 15, l = 15)),
        axis.title.x = element_text(color = "black",
                                    size = 15,
                                    face = "bold",
                                    margin = margin(t = 15, b = 15)),
        axis.text = element_text(color = "black",
                                 size = 12,
                                 margin = margin(r = 5, t = 5)),
        legend.position = c(0.92,0.92),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold"),
        legend.key = element_rect(color = "black"))  -> Ext_Data_Figure_1


### ----------------------------------------------------------------------------

ggsave("./Final_Code/Plots/Extended_Data_Figure_1.tiff",
       plot = Ext_Data_Figure_1,
       units = "in",
       width = 4,
       height = 5,
       device = "tiff",
       dpi = 300,
       scale = 2.2)
