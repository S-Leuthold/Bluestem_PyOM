read_csv("./Final_Code/Data/Leuthold et al., 2025 - Bluestem 10 Year Sampling.csv") %>% 
  select(-c(Full_name)) %>% 
  filter(Project == "Bluestem") -> Bluestem_Data


read_csv("./Final_Code/Data/Leuthold et al., 2025 - Bluestem Raw.csv") %>% 
  filter(`Litter type (B, E, NL)` == "E") %>% 
  select(Replicate,
         Harvest,
         Fraction, 
         `Soil depth (cm)`, 
         `bulk density (g/cm3)`,
         `volume (cm3)`,
         `% Carbon`,
         `% Nitrogen`) %>% 
  filter(Fraction == "BULK") %>% 
  mutate(`Soil depth (cm)` = case_when(`Soil depth (cm)` == "5-Feb" ~ "2-5",
                                       `Soil depth (cm)` == "10-May" ~ "5-10",
                                       `Soil depth (cm)` == "20-Oct" ~ "10-20",
                                       TRUE ~ `Soil depth (cm)`),
         Replicate = as.numeric(str_sub(Replicate, 2, 2))) %>% 
  filter(`Soil depth (cm)` %in% c("0-2", "2-5")) %>% 
  mutate(Soil_Height = case_when(`Soil depth (cm)` == "2-5" ~ 3,
                                 `Soil depth (cm)` == "0-2" ~ 2),
         Soil_C_Stock = (Soil_Height * 10000 *`bulk density (g/cm3)` * `% Carbon`/100),
         Soil_N_Stock = (Soil_Height * 10000 * `bulk density (g/cm3)` * `% Nitrogen`/100)) %>% 
  group_by(Harvest, Replicate) %>%
  summarize(Soil_C_Stock = sum(Soil_C_Stock),
            Soil_N_Stock = sum(Soil_N_Stock)) %>% 
  group_by(Harvest, Replicate) %>%
  summarize(Mean_SOC_Stock = mean(Soil_C_Stock),
            Mean_STN_Stock = mean(Soil_N_Stock)) %>% 
  bind_rows(.,
            Bluestem_Data %>% 
              mutate(Soil_Height = (Bottom_Depth_cm - Top_Depth_cm),
                     Bulk_Density_g_cm3_Adj = case_when(Depth_cm == "0-2" ~ 0.702,
                                                        Depth_cm == "2-5" ~ 0.923,
                                                        TRUE ~ Bulk_Density_g_cm3),
                     Soil_C_Stock = Soil_Height * (10000) * Bulk_Density_g_cm3_Adj * Bulk_C_Pct/100,
                     Soil_N_Stock = Soil_Height * (10000) * Bulk_Density_g_cm3 * Bulk_N_Pct/100) %>% 
              select(Rep, 
                     Fumigation_Treatment,
                     Depth_cm,
                     Soil_C_Stock,
                     Soil_N_Stock) %>% 
              rename(Replicate = Rep) %>% 
              filter(Depth_cm %in% c("0-2", "2-5")) %>% 
              group_by(Replicate, Fumigation_Treatment) %>% 
              summarize(Soil_C_Stock = sum(Soil_C_Stock),
                        Soil_N_Stock = sum(Soil_N_Stock)) %>% 
              mutate(Harvest = "H7") %>% 
              ungroup() %>% 
              group_by(Harvest, Replicate) %>% 
              summarize(Mean_SOC_Stock = mean(Soil_C_Stock),
                        Mean_STN_Stock = mean(Soil_N_Stock))) %>% 
  mutate(Bulk_CN = Mean_SOC_Stock/Mean_STN_Stock) -> Bulk_SOM_Stocks


Bulk_SOM_Stocks %>% 
  lme4::lmer(Mean_SOC_Stock ~ Harvest + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest = "Tukey")) %>% 
  cld()

Bulk_SOM_Stocks %>% 
  lme4::lmer(Bulk_CN ~ Harvest + (1|Replicate), data = .) %>% 
  glht(., mcp(Harvest = "Tukey")) %>% 
  cld()

Bulk_SOM_Stocks %>% 
  group_by(Harvest) %>% 
  summarize(SE = sd(Bulk_CN)/sqrt(n()),Bulk_CN = mean(Bulk_CN)) %>% 
  ggplot() +
  geom_col(aes(x = Harvest,
               y = Bulk_CN))

Bulk_SOM_Stocks %>% 
  ungroup() %>% 
  summarize(C_SE = sd(Mean_SOC_Stock)/sqrt(n()),
            C = mean(Mean_SOC_Stock),
            N_SE = sd(Mean_STN_Stock)/sqrt(n()),
            N = mean(Mean_STN_Stock))


Bulk_SOM_Stocks %>% 
  group_by(Harvest) %>% 
  summarize(Mean_C = mean(Mean_SOC_Stock ),
            SE_C = sd(Mean_SOC_Stock )/2) %>% 
  mutate(Harvest = case_when(Harvest == "H1" ~ "7 Months",
                             Harvest == "H2" ~ "12 Months",
                             Harvest == "H3" ~ "18 Months",
                             Harvest == "H4" ~ "24 Months",
                             Harvest == "H6" ~ "36 Months",
                             Harvest == "H7" ~ "120 Months"),
         Harvest = fct_relevel(Harvest, "7 Months", "12 Months", "18 Months", "24 Months", "36 Months", "120 Months")) %>% 
  ggplot() +
  geom_col(aes(x = Harvest, 
               y = Mean_C),
           fill = "wheat3",
           color = "black") +
  geom_errorbar(aes(x = Harvest,
                    ymin = Mean_C - SE_C,
                    ymax = Mean_C + SE_C),
                width = 0.5) +
  geom_point(data = Bulk_SOM_Stocks %>% 
                     mutate(Harvest = case_when(Harvest == "H1" ~ "7 Months",
                                                Harvest == "H2" ~ "12 Months",
                                                Harvest == "H3" ~ "18 Months",
                                                Harvest == "H4" ~ "24 Months",
                                                Harvest == "H6" ~ "36 Months",
                                                Harvest == "H7" ~ "120 Months"),
                            Harvest = as.factor(Harvest),
                            Harvest = fct_relevel(Harvest, "7 Months", "12 Months", "18 Months", "24 Months", "36 Months", "120 Months")),
             aes(x = Harvest, 
                 y = Mean_SOC_Stock),
             shape = 21,
             size = 3,
             fill = "wheat2") +
  scale_y_continuous(name =  expression("Bulk Soil C Stock (g C "*m^-2*")"),
                     limits = c(0, 4500),
                     expand = c(0,0)) + 
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 20, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = c(1, 1),
        legend.justification = c(1,1),
        legend.text = element_text(size = 15),
        legend.background = element_rect(color = "black"),
        title = element_text(size = 15, face = "bold")) -> Ext_Data_Figure_3

ggsave(filename = "./Final_Code/Plots/Extended_Data_Figure_3.tiff",
       plot = Ext_Data_Figure_3,
       units = "px",
       height = 1200,
       width = 1200,
       scale = 2,
       dpi = "print")


Bulk_SOM_Stocks %>% 
  mutate(CN = Mean_SOC_Stock/Mean_STN_Stock) %>% 
  group_by(Harvest) %>% 
  summarize(Mean_CN = mean(CN),
            SE_CN = sd(CN)/2) %>% 
  mutate(Harvest = case_when(Harvest == "H1" ~ "7 Months",
                             Harvest == "H2" ~ "12 Months",
                             Harvest == "H3" ~ "18 Months",
                             Harvest == "H4" ~ "24 Months",
                             Harvest == "H6" ~ "36 Months",
                             Harvest == "H7" ~ "120 Months"),
         Harvest = fct_relevel(Harvest, "7 Months", "12 Months", "18 Months", "24 Months", "36 Months", "120 Months")) %>% 
  ggplot() +
  geom_col(aes(x = Harvest, 
               y = Mean_CN),
           fill = "wheat3",
           color = "black") +
  geom_errorbar(aes(x = Harvest,
                    ymin = Mean_CN - SE_CN,
                    ymax = Mean_CN + SE_CN),
                width = 0.5) +
  geom_point(data = Bulk_SOM_Stocks %>% 
               mutate(CN = Mean_SOC_Stock/Mean_STN_Stock,
                      Harvest = case_when(Harvest == "H1" ~ "7 Months",
                                          Harvest == "H2" ~ "12 Months",
                                          Harvest == "H3" ~ "18 Months",
                                          Harvest == "H4" ~ "24 Months",
                                          Harvest == "H6" ~ "36 Months",
                                          Harvest == "H7" ~ "120 Months"),
                      Harvest = as.factor(Harvest), 
                      Harvest = fct_relevel(Harvest, "7 Months", "12 Months", "18 Months", "24 Months", "36 Months", "120 Months")),
             aes(x = Harvest, 
                 y = CN),
             shape = 21,
             size = 3,
             fill = "wheat2") +
  scale_y_continuous(name =  "Bulk Soil C:N",
                     limits = c(0, 15.5),
                     expand = c(0,0)) + 
  theme_classic2() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA),
        axis.title = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = c(1, 1),
        legend.justification = c(1,1),
        legend.text = element_text(size = 15),
        legend.background = element_rect(color = "black"),
        title = element_text(size = 15, face = "bold")) -> Ext_Data_Figure_4


ggsave(filename = "./Final_Code/Plots/Extended_Data_Figure_4.tiff",
       plot = Ext_Data_Figure_4,
       units = "px",
       height = 1200,
       width = 1200,
       scale = 2,
       dpi = "print")
