### Processing script for FTIR Data from bluestem/PyOM work, with figure
### Written by: Sam Leuthold (sam.leuthold@gmail.com)
### Last Updated: 12/16/2024

### Setup ----------------------------------------------------------------------

rm(list = ls())

require(pacman)

p_load("prospectr",
       "tidyverse",
       "progress",
       "opusreader",
       "Bolstad2",
       "lmerTest",
       "patchwork",
       "conflicted")

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("resample", "prospectr")
conflict_prefer("lmer", "lmerTest")

### ----------------------------------------------------------------------------

list.files("./Final_Code/Data/FTIR_Data/") %>% 
  tibble() %>% 
  filter(!grepl("NDN", .),
         !grepl("Scans", .)) %>%
  pull(.) -> PyOM_Files

### ----------------------------------------------------------------------------

Spectral_Output <- NULL

pb <- progress_bar$new(total = length(PyOM_Files))

for(i in 1:length(PyOM_Files)){

  opus_read(file = paste0("./Final_Code/Data/FTIR_Data/", PyOM_Files[i]),
            simplify = FALSE,
            progress = FALSE) %>% 
    pluck(., "spec") %>% 
    as_tibble(.) %>% 
    pivot_longer(cols = everything(),
                 names_to = "Wavenumber",
                 values_to = "Absorbance") %>% 
    mutate(Wavenumber = as.numeric(Wavenumber),
           File_ID = PyOM_Files[i],
           Year = str_split_i(string = File_ID,
                              pattern = "_",
                              i = 2),
           Depth = str_split_i(string = File_ID,
                               pattern = "_",
                               i = 3),
           Replicate = str_split_i(string = File_ID,
                                     pattern = "_",
                                     i = 4),
           Scan_ID = str_split_i(string = File_ID,
                                 pattern = "_",
                                 i = 5),
           Well_ID = str_split_i(string = File_ID,
                                 pattern = "_",
                                 i = 6),
           Well_ID = str_split_i(string = Well_ID,
                                 pattern = "\\.",
                                 i = 1),
           Sample_Name = paste0(Replicate, "_", Depth, "cm_Year-", Year)) %>% 
    filter(round(Wavenumber) >= 600,
           round(Wavenumber) <= 4000) %>% 
    select(Sample_Name,
           Replicate,
           Depth,
           Year,
           Scan_ID,
           Well_ID,
           Wavenumber,
           Absorbance)  -> Spectral_Information
  
  
  ## Baseline correct the spectra ----------------------------------------------
  
    resample(X = Spectral_Information$Absorbance,
            wav = Spectral_Information$Wavenumber,
            new.wav = seq(600, 4000, 2),
            interpol = "spline") %>% 
    unname() %>% 
    tibble(Absorbance = .,
           Wavenumber = seq(600, 4000, 2)) %>% 
    mutate(Sample_Name = unique(Spectral_Information$Sample_Name),
           Replicate = unique(Spectral_Information$Replicate),
           Depth = unique(Spectral_Information$Depth),
           Year = unique(Spectral_Information$Year),
           Scan_ID = unique(Spectral_Information$Scan_ID),
           Well_ID = unique(Spectral_Information$Well_ID)) %>% 
    select(Sample_Name,
           Replicate,
           Depth,
           Year,
           Scan_ID,
           Well_ID,
           Wavenumber,
           Absorbance) -> Spectral_Information
    
  ## Smooth the spectra --------------------------------------------------------
  
  Window_Size  <- 9
  Start_Row <-  1 + ((Window_Size-1)/2)
  End_Row <- nrow(Spectral_Information) - ((Window_Size - 1)/2)
  
  Spectral_Information %>% 
    pull(Absorbance) %>% 
    as.matrix(.) %>% 
    t(.) %>% 
    savitzkyGolay(X = ., 
                  m = 0, 
                  p = 1,
                  w = Window_Size) %>% 
    t(.) %>% 
    as_tibble(.) %>% 
    bind_cols(Spectral_Information %>% 
                slice(Start_Row:End_Row),
              .) %>% 
    rename(SG_Absorbance = V1)  -> Spectral_Information
  
  ## Remove light scatter effects ----------------------------------------------
  
  Spectral_Information %>% 
    pull(SG_Absorbance) %>% 
    as.matrix(.) %>% 
    t(.) %>% 
    standardNormalVariate() %>% 
    t(.) %>% 
    as_tibble(.) %>% 
    bind_cols(Spectral_Information,
              .) %>% 
    rename(SNV_Absorbance = V1) %>% 
    mutate(SNV_Absorbance = as.numeric(SNV_Absorbance)) -> Spectral_Information

  ## Output results ------------------------------------------------------------
  
  bind_rows(Spectral_Output, 
            Spectral_Information) -> Spectral_Output
  
  pb$tick()
  
}

### Clean up -------------------------------------------------------------------

Spectral_Output %>% 
  mutate(Sample_Type = case_when(grepl("R", Sample_Name) ~ "POM_Char",
                                 grepl("Q", Sample_Name) ~ "POM_Bare",
                               TRUE ~ "PyOM")) %>% 
  group_by(Sample_Name,
           Sample_Type,
           Replicate,
           Depth,
           Year,
           Wavenumber) %>% 
  summarise(Absorbance = mean(SNV_Absorbance)) %>% 
  ungroup() -> Spectral_Output

  
### Spectral Subtraction -------------------------------------------------------

## Grab the spectra that didn't have PyOM --------------------------------------

Spectral_Output %>%
  filter(Sample_Type == "POM_Bare") %>%
  group_by(Depth,
           Wavenumber) %>%
  summarise(Absorbance = mean(Absorbance)) %>%
  pivot_wider(names_from = Depth,
              values_from = Absorbance) %>%
  rename(POM_Bare_0_2 = `0-2`,
         POM_Bare_2_5 = `2-5`) -> POM_Bare

POM_Bare %>% 
  pull(POM_Bare_0_2) %>% 
  as.matrix(.) %>% 
  t(.) %>% 
  baseline(.) %>% 
  as_tibble(.) %>% 
  pull(value) -> Baseline_POM_0_2

POM_Bare %>% 
  pull(POM_Bare_2_5) %>% 
  as.matrix(.) %>% 
  t(.) %>% 
  baseline(.) %>% 
  as_tibble(.) %>% 
  pull(value) -> Baseline_POM_2_5
  

## Rearrange remaining data ----------------------------------------------------

Spectral_Output %>%
  filter(Sample_Type == "POM_Char") %>%
  mutate(Year_Depth = paste0(Year, "_", Depth)) %>%
  group_by(Year,
           Year_Depth,
           Replicate,
           Wavenumber) %>%
  summarise(Absorbance = mean(Absorbance)) %>%
  pivot_wider(names_from = Year_Depth,
              values_from = Absorbance) %>%
  rename(POM_Char_1_0_2 = `1_0-2`,
         POM_Char_1_2_5 = `1_2-5`,
         POM_Char_10_0_2 = `10_0-2`,
         POM_Char_10_2_5 = `10_2-5`) %>% 
  mutate(Rep_Year = paste0(Year, "_", Replicate)) -> POM_Char

## Join the data, and subtract -------------------------------------------------

Subtraction_Output <- NULL

unique(POM_Char$Rep_Year)[2]

for(i in unique(POM_Char$Rep_Year)){
  
  POM_Char %>% 
    filter(Rep_Year == i) -> temp_Char

  
  if(unique(temp_Char$Year) == 1){
    
    temp_Char %>% 
      pull(POM_Char_1_0_2) %>% 
      as.matrix(.) %>% 
      t(.) %>% 
      baseline(.) %>% 
      as_tibble(.) %>% 
      rename(Baseline_PyOM_1_0_2 = value) %>% 
      mutate(Subtracted_1_0_2 = Baseline_PyOM_1_0_2 - Baseline_POM_0_2) %>% 
      pull(Subtracted_1_0_2) -> Subtracted_0_2

    temp_Char %>% 
      pull(POM_Char_1_2_5) %>% 
      as.matrix(.) %>% 
      t(.) %>% 
      baseline(.) %>% 
      as_tibble(.) %>% 
      rename(Baseline_PyOM_1_2_5 = value) %>% 
      mutate(Subtracted_1_2_5 = Baseline_PyOM_1_2_5 - Baseline_POM_2_5) %>% 
      pull(Subtracted_1_2_5) -> Subtracted_2_5
    
      } 
  
  if(unique(temp_Char$Year) == 10){
    
    temp_Char %>% 
      pull(POM_Char_10_0_2) %>% 
      as.matrix(.) %>% 
      t(.) %>% 
      baseline(.) %>% 
      as_tibble(.) %>% 
      rename(Baseline_PyOM_10_0_2 = value) %>% 
      mutate(Subtracted_10_0_2 = Baseline_PyOM_10_0_2 - Baseline_POM_0_2) %>% 
      pull(Subtracted_10_0_2) -> Subtracted_0_2
    
    temp_Char %>% 
      pull(POM_Char_10_2_5) %>% 
      as.matrix(.) %>% 
      t(.) %>% 
      baseline(.) %>% 
      as_tibble(.) %>% 
      rename(Baseline_PyOM_10_2_5 = value) %>% 
      mutate(Subtracted_10_2_5 = Baseline_PyOM_10_2_5 - Baseline_POM_2_5) %>% 
      pull(Subtracted_10_2_5) -> Subtracted_2_5
    
  } 
  
  
  tibble(Rep_Year = i,
         Wavenumber = temp_Char$Wavenumber,
         Subtracted_0_2 = Subtracted_0_2,
         Subtracted_2_5 = Subtracted_2_5) -> output
  
  Subtraction_Output <- bind_rows(Subtraction_Output,
                                  output)
    
  }
  
## Clean up --------------------------------------------------------------------

Subtraction_Output %>% 
  pivot_longer(names_to = "Depth",
               cols = Subtracted_0_2:Subtracted_2_5) %>% 
  mutate(Replicate = str_split_i(Rep_Year,
                                 "_",
                                 2),
         Year = str_split_i(Rep_Year,
                            "_",
                            1),
         Depth = str_split_i(Depth,
                             "_",
                             2),
         Depth = case_when(Depth == 0 ~ "0-2",
                           Depth == 2 ~ "2-5"),
         Sample_Name = paste0("POM_", Year, "_", Depth, "_", Replicate)) %>% 
  rename(Subtracted_Absorbance = value) %>% 
  select(Sample_Name,
         Year,
         Replicate,
         Depth,
         Wavenumber,
         Subtracted_Absorbance)  -> POM_Corrected_Spectra

### Create a dataframe with specified peak areas --------------------------------

Spectral_Output %>%
  filter(Sample_Type != "PyOM") %>% 
  group_by(Depth,
           Wavenumber,
           Replicate,
           Year) %>% 
  summarise(Absorbance = mean(Absorbance)) %>% 
  filter(Wavenumber %in% c(seq(1680, 1730, 2),
                           seq(1570, 1610, 2),
                           seq(1500, 1510, 2),
                           seq(1380, 1430, 2),
                           seq(1210, 1260, 2),
                           seq(1020, 1060, 2),
                           seq(878, 882, 2),
                           seq(804, 806, 2),
                           seq(744, 746, 2))) %>% 
  ungroup() %>% 
  mutate(Assignment = case_when(Wavenumber %in% seq(1680, 1730, 2) ~ "Aromatic Carbonyl/Carboxyl C=O Stretch",
                                Wavenumber %in% seq(1570, 1610, 2) ~ "C=C Stretch",
                                Wavenumber %in% seq(1500, 1510, 2) ~ "Lignin, Aromatic C=C Stretch",
                                Wavenumber %in% seq(1380, 1430, 2) ~ "Aromatic C=C stretch",
                                Wavenumber %in% seq(1210, 1260, 2) ~ "Cellulose",
                                Wavenumber %in% seq(1020, 1060, 2) ~ "Aliphatic C-O- and Alcohol C-O stretch",
                                Wavenumber %in% seq(878, 882, 2) ~   "C-H aromatic bending deformation",
                                Wavenumber %in% seq(804, 806, 2) ~ "C-H aromatic bending deformation",
                                Wavenumber %in% seq(744, 746, 2) ~ "C-H aromatic bending deformation")) %>% 
  mutate(ID = paste0(Depth, Replicate, Year, Assignment))-> temp.integration.data


AUC_Data <- NULL

for(i in 1:length(unique(temp.integration.data$ID))){
  
  ## Do a local baseline correction on the wavenumbers of interest.
  
  temp.integration.data %>% 
    filter(ID == unique(temp.integration.data$ID)[i]) %>% 
    pull(Absorbance) %>%
    as.matrix(.) %>%
    t(.) %>%
    baseline(.) %>%
    as_tibble(.) %>%
    mutate(Absorbance = as.numeric(value)) %>% 
    select(Absorbance) %>% 
    bind_cols(tibble(Wavenumber = temp.integration.data %>% 
                                    filter(ID == unique(temp.integration.data$ID)[i]) %>% 
                                    pull(Wavenumber)),
              .) -> temp_A
  
  ## Intergrate under the curve ------------------------------------------------
  
  sintegral(x = temp_A$Wavenumber, 
            fx = temp_A$Absorbance) %>% 
    pluck(., "int") -> AUC
  
  temp.integration.data %>% 
    filter(ID == unique(temp.integration.data$ID)[i])  -> temp.info
  
  ## Pull results together in a dataframe --------------------------------------
  
  tibble(Assignment = temp.info$Assignment,
         Depth = temp.info$Depth,
         Replicate = temp.info$Replicate,
         Year = temp.info$Year,
         AUC = AUC) %>% 
    distinct() -> output
  
  ## Output results ------------------------------------------------------------
  
  bind_rows(AUC_Data,
            output) -> AUC_Data

}

### Compute the relative peak areas --------------------------------------------

AUC_Data %>%
  group_by(Replicate, Year) %>% 
  mutate(RPA = (AUC/sum(AUC)) * 100) %>% 
  ungroup() -> AUC_Data

### Do paired T-Tests on the AUCs from year 1 and year 10 -----------------------

## Create somewhere for the data to end up -------------------------------------

T_Test_Results_02 <- NULL
T_Test_Results_25 <- NULL

## Loop through the combinations -----------------------------------------------

for(i in unique(AUC_Data$Assignment)){
  
  ## Pull the AUC data for an AUC assignment -----------------------------------
  
  AUC_Data %>%
    filter(Assignment == i,
           Depth == "0-2") %>% 
    select(-RPA) %>% 
    pivot_wider(names_from = Year,
                values_from = AUC) %>% 
    rename(one = `1`, 
           ten = `10`)  -> temp_AUC 
    
  ## Compare -------------------------------------------------------------------
  
  t.test(x = temp_AUC$one,
         y = temp_AUC$ten,
         paired = TRUE) -> results
  
  ## Pull results together and output ------------------------------------------
  
  results %>% 
    pluck(., "p.value") -> results
  
  tibble(Assignment = i,
         p_value = results,
         Depth = "0-2") -> results
  
  
  
  T_Test_Results_02 <- bind_rows(T_Test_Results_02,
                                 results)
  
  ## Do the same thing for 2-5 cm ----------------------------------------------
  
  AUC_Data %>%
    filter(Assignment == i,
           Depth == "2-5") %>% 
    select(-RPA) %>% 
    pivot_wider(names_from = Year,
                values_from = AUC) %>% 
    rename(one = `1`, 
           ten = `10`)  -> temp_AUC 
  
  t.test(x = temp_AUC$one,
         y = temp_AUC$ten,
         paired = TRUE) -> results
  
  results %>% 
    pluck(., "p.value") -> results
  
  tibble(Assignment = i,
         p_value = results,
         Depth = "2-5") -> results
  
  T_Test_Results_25 <- bind_rows(T_Test_Results_25,
                                 results)
  
}

## Clean up results from tests -------------------------------------------------

bind_rows(T_Test_Results_02,
          T_Test_Results_25) %>% 
  mutate(Significant = case_when(p_value < 0.05 ~ "*",
                                 TRUE ~ " "),
         Depth = case_when(Depth == "0-2" ~ "0 - 2 cm",
                           Depth == "2-5" ~ "2 - 5 cm")) %>% 
  mutate(Assignment = fct_relevel(Assignment,
                                  "Aromatic Carbonyl/Carboxyl C=O Stretch",
                                  "C=C Stretch",
                                  "Lignin, Aromatic C=C Stretch",
                                  "Aromatic C=C stretch",
                                  "Cellulose",
                                  "Aliphatic C-O- and Alcohol C-O stretch",
                                  "C-H aromatic bending deformation")) -> T_Test_Results

print(T_Test_Results,
      n = nrow(T_Test_Results))

### Plot the RPA Results -------------------------------------------------------

AUC_Data %>% 
  mutate(Assignment = fct_relevel(Assignment,
                                  "Aromatic Carbonyl/Carboxyl C=O Stretch",
                                  "C=C Stretch",
                                  "Lignin, Aromatic C=C Stretch",
                                  "Aromatic C=C stretch",
                                  "Cellulose",
                                  "Aliphatic C-O- and Alcohol C-O stretch",
                                  "C-H aromatic bending deformation")) -> AUC_Data

AUC_Data %>% 
  group_by(Year, Depth, Assignment) %>% 
  summarise(SE_AUC = sd(AUC)/sqrt(n()),
            Mean_AUC = mean(AUC)) %>% 
  mutate(Depth = case_when(Depth == "0-2" ~ "0 - 2 cm",
                           Depth == "2-5" ~ "2 - 5 cm")) %>% 
  filter(Assignment !=  "C-H aromatic stretch (increase)") %>% 
  ggplot(data = .) +
  geom_col(aes(x = Assignment,
               y = Mean_AUC,
               fill = Year),
           position = "dodge",
           color = "black") +
  
  facet_grid(Depth~.) +
  geom_text(data = T_Test_Results %>% 
              filter(Assignment !=  "C-H aromatic stretch (increase)"),
            aes(x = Assignment,
                y = 3.25,
                label = Significant),
            size = 10) +
  scale_fill_manual(name = "Year",
                     values = c("1" = "#83AF9B",
                                "10" = "#F9CDAD")) +
  geom_errorbar(aes(x = Assignment,
                    ymin = Mean_AUC - SE_AUC,
                    ymax = Mean_AUC + SE_AUC,
                    group = Year),
                position = position_dodge(width = 0.85),
                width = 0.25)  +
  scale_y_continuous(name = "Average Area Under Curve (Abs. Units)",
                     limits = c(0, 5.5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(transform = ~ . * 10, name = derive())) +
  
  theme_void() +
  theme(panel.border = element_rect(fill = "transparent",
                                    color = "black",
                                    size = 1),
        panel.background = element_rect(fill = "white", color = "transparent"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = "transparent"),
        axis.title.y = element_text(color = "black",
                                    size = 15,
                                    angle = 90,
                                    face = "bold",
                                    margin = margin(r = 15, l = 15)),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black",
                                 size = 12,
                                 margin = margin(r = 5, t = 5)),
        axis.text.x = element_text(angle = 50,
                                   hjust = 1,
                                   vjust = 1),
        legend.position = "none",
        axis.title.y.right = element_text(angle = 270),
        plot.margin = margin(t = 7, r = 7, b = 7, l = 100)) -> Right_Panel

  Right_Panel

  
  
### Do a baseline correction on the original data for plotting -----------------


Spectral_Output_B <- NULL
  
for(i in unique(Spectral_Output$Sample_Name)){
  
  Spectral_Output %>% 
    filter(Sample_Name == i) %>% 
    pull(Absorbance) %>% 
    as.matrix(.) %>% 
    t(.) %>% 
    baseline(.) %>% 
    as_tibble(.) %>% 
    rename(BL_Absorbance = value) %>% 
    bind_cols(Spectral_Output %>% 
                filter(Sample_Name == i),
              .) -> output
  
  bind_rows(Spectral_Output_B,
            output) -> Spectral_Output_B
}  
  
### ----------------------------------------------------------------------------

Spectral_Output_B %>%
  filter(Sample_Type != "PyOM") %>% 
  group_by(Depth,
           Year,
           Wavenumber) %>% 
  summarise(Absorbance = mean(BL_Absorbance)) %>% 
  ggplot(data = .) +
  geom_line(aes(x = Wavenumber,
                y = Absorbance, 
                color = Year)) +
  scale_color_manual(name = "Year",
                     values = c("1" = "#83AF9B",
                                "10" = "#F9CDAD")) +
  facet_grid(Depth~.) +
  scale_x_reverse(name = expression("Wavenumber (c"*m^-1*")", face = "bold"),) +
  scale_y_continuous("Absorbance") +
  geom_hline(yintercept = 0) +
  theme_void() +
  theme(panel.border = element_rect(fill = "transparent",
                                    color = "black",
                                    size = 1),
        panel.background = element_rect(fill = "white", color = "transparent"),
        strip.text = element_text(color = "black",
                                  angle = 270,
                                  size = 13),
        plot.background = element_rect(fill = "white", color = "transparent"),
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
        legend.position = "bottom",
        plot.margin = margin(t = 7, r = 7, b = 7, l = 10)) -> Left_Panel

### ----------------------------------------------------------------------------

Left_Panel + Right_Panel +
  plot_layout(ncol = 1,
              heights = c(3/8, 5/8)) -> Ext_Data_Figure_2

### ----------------------------------------------------------------------------

ggsave(filename = "./Final_Code/Plots/Extended_Data_Figure_2.tiff",
       plot = Ext_Data_Figure_2,
       device = "tiff",
       units = "in",
       scale = 2,
       height = 8,
       width = 5)




  

