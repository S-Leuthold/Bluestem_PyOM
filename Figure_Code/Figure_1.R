### Script to calculate the CV of variation for the litter derived C and N for POM, CHAOM, and MAOM
### Written by: Sam Leuthold (sam.leuthold@gmail.com)
### Last Updated: 12/14/2024

### Setup ----------------------------------------------------------------------

require(pacman)

p_load("tidyverse",
       "lmerTest")


### Source files ----------------------------------------------------------------

source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_Data.R")

source("./Final_Code/Source_Scripts/Calculate_Bluestem_C_N_Distirbution.R")

source("./Final_Code/Source_Scripts/Calculate_PyOM_Data.R")

source("./Final_Code/Source_Scripts/Calculate_PyOM_C_N_Distribution.R")

### Data manipulation ----------------------------------------------------------

C_N_Mass_Dist %>%
  filter(Fraction != "Lost") %>% 
  group_by(Harvest_Num, Replicate) %>% 
  summarise(Carbon_g_m2 = sum(Carbon_g_m2, na.rm = TRUE),
            Carbon_Mass_Prop = sum(Carbon_Mass_Prop, na.rm = TRUE)/100) -> Bluestem_Intermediate 

### ----------------------------------------------------------------------------
  
ggplot(data = Bluestem_Intermediate,
       aes(x = Harvest_Num,
           y = Carbon_Mass_Prop)) +
  geom_point()

### ----------------------------------------------------------------------------
  
BS_Decay_Model <- nls(Carbon_Mass_Prop ~ exp(-k * as.numeric(Harvest_Num)), data = Bluestem_Intermediate, start = list(k = 0.001))

BS_Decay_Model

### ----------------------------------------------------------------------------

predict(BS_Decay_Model, 
        newdata = data.frame(Harvest_Num = seq(0,120,5))) %>% 
  bind_cols(.,
            data.frame(Harvest_Num = seq(0,120,5))) %>% 
  rename(Carbon_Mass_Prop = `...1`) -> Bluestem_Predicted

## Check fit -------------------------------------------------------------------

ggplot(data = Bluestem_Intermediate,
       aes(x = Harvest_Num,
           y = Carbon_Mass_Prop)) +
  geom_point() +
  geom_line(data = Bluestem_Predicted,
            aes(x = Harvest_Num,
                y = Carbon_Mass_Prop))
### ----------------------------------------------------------------------------

Bluestem_Intermediate %>% 
  ungroup() %>% 
  group_by(Harvest_Num) %>% 
  summarize(Carbon_g_m2_SE = sd(Carbon_g_m2, na.rm = TRUE)/sqrt(n()),
            Carbon_g_m2 = mean(Carbon_g_m2, na.rm = TRUE),
            Carbon_Mass_Prop_SE = sd(Carbon_Mass_Prop, na.rm = TRUE)/sqrt(n()),
            Carbon_Mass_Prop = mean(Carbon_Mass_Prop, na.rm = TRUE)) %>% 
  mutate(Grouping = "Litter") -> Bluestem_Intermediate

###-----------------------------------------------------------------------------
### ----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

BlackC_C_N_Mass_Dist %>% 
  filter(Fraction != "Lost") %>% 
  group_by(Harvest_Num, Replicate) %>% 
  summarise(Carbon_g_m2 = sum(Carbon_g_m2, na.rm = TRUE),
            Carbon_Mass_Prop = sum(Carbon_Mass_Prop)/100) %>% 
  bind_rows(.,
            tibble(Grouping = "PyOM",
                   Harvest_Num = 0,
                   Carbon_g_m2 = 104.2787,
                   Carbon_Mass_Prop = 0.62189)) %>% 
  mutate(Carbon_Mass_Prop = (Carbon_g_m2/167.68),
         Input = "PyOM",
         Grouping = "PyOM") %>% 
  mutate(Harvest_Num = case_when(Harvest_Num == 108 ~ 101,
                                 TRUE ~ Harvest_Num)) -> PyOM_Intermediate

PyOM_Intermediate %>% 
  group_by(Harvest_Num) %>% 
  summarize(Mean = mean(Carbon_Mass_Prop))

### ----------------------------------------------------------------------------

ggplot(data = PyOM_Intermediate,
       aes(x = Harvest_Num,
           y = Carbon_Mass_Prop)) +
  geom_point()

### ----------------------------------------------------------------------------
Py_Decay_Model <- nls(Carbon_Mass_Prop ~ A + (B-A) * exp(-k * as.numeric(Harvest_Num)), data = PyOM_Intermediate, 
                      start = list(k = 0.05,
                                   A = 0.4,
                                   B = 1))

Py_Decay_Model

### ----------------------------------------------------------------------------

predict(Py_Decay_Model, 
        newdata = data.frame(Harvest_Num = seq(0,120,5))) %>% 
  bind_cols(.,
            data.frame(Harvest_Num = seq(0,120,5))) %>% 
  rename(Carbon_Mass_Prop = `...1`) -> PyOM_Predicted

## Check fit -------------------------------------------------------------------

ggplot(data = PyOM_Intermediate %>% 
         filter(Harvest_Num != 4 | Replicate != 2),
       aes(x = Harvest_Num,
           y = Carbon_Mass_Prop)) +
  geom_point() +
  geom_line(data = PyOM_Predicted,
            aes(x = Harvest_Num,
                y = Carbon_Mass_Prop))


### ----------------------------------------------------------------------------  
  PyOM_Intermediate %>% 
  ungroup() %>% 
  group_by(Grouping, Harvest_Num) %>% 
  summarize(Carbon_g_m2_SE = sd(Carbon_g_m2, na.rm = TRUE)/sqrt(n()),
            Carbon_g_m2 = mean(Carbon_g_m2, na.rm = TRUE),
            Carbon_Mass_Prop_SE = sd(Carbon_Mass_Prop, na.rm = TRUE)/sqrt(n()),
            Carbon_Mass_Prop = mean(Carbon_Mass_Prop, na.rm = TRUE)) %>% 
  bind_rows(.,
            tibble(Grouping = "PyOM",
                   Harvest_Num = -1,
                   Carbon_g_m2 = 167.68,
                   Carbon_Mass_Prop = 1)) -> PyOM_Intermediate

### ----------------------------------------------------------------------------

bind_rows(Bluestem_Intermediate,
          PyOM_Intermediate) -> Lehman_Data

### ----------------------------------------------------------------------------

ggplot() +
  geom_line(data = Lehman_Data %>% 
              filter(Harvest_Num <= 0),
            aes(x = Harvest_Num,
                y = Carbon_Mass_Prop,
                group = Grouping),
            lty = 2) +
  geom_line(data = Bluestem_Predicted,
            aes(x = Harvest_Num,
                y = Carbon_Mass_Prop)) +
  geom_line(data = PyOM_Predicted,
            aes(x = Harvest_Num,
                y = Carbon_Mass_Prop),
            lty = 2) +
  geom_errorbar(data = Lehman_Data,
                aes(x = Harvest_Num,
                    ymax = Carbon_Mass_Prop + Carbon_Mass_Prop_SE,
                    ymin = Carbon_Mass_Prop - Carbon_Mass_Prop_SE,
                    group = Grouping)) +
  geom_point(data = Lehman_Data,
                 aes(x = Harvest_Num,
                     y = Carbon_Mass_Prop,
                     shape = Grouping),
             fill = "black",
             size = 4) + 
  scale_shape_manual(name = "C Source",
                    values = c("PyOM" = 21,
                               "Litter" = 24)) +
  scale_y_continuous(name = "Carbon input remaining (% of initial)",
                     breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(name = "Harvest Month",
                     breaks = seq(0, 120, 15)) +
  annotate(geom = "text",
           x = 1.5, 
           y = 0.8,
           label = "Pyrolosis Loss",
           size = 5,
           angle = 270) +
  theme_classic2() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 21),
        axis.text = element_text(color = "black", size = 15),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.title = element_text(face = "bold", color = "black", size = 17),
        legend.text = element_text(color = "black", size = 14),
        legend.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_rect(fill = "transparent", color = "black"))

ggsave("./Final_Code/Plots/Figure_1.pdf",
       plot = last_plot(),
       device = "pdf",
       dpi = "print",
       height = 7,
       width = 11,
       units = "in")

