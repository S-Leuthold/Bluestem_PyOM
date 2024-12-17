### Source script for calculating time-series of litter derived C and N.
### Adapted from Bluestem_Litter_Calculations.R
### Author: Sam Leuthold (sam.leuthold@gmail.com)
### Last update: 12/12/2024

### Notes and update log -------------------------------------------------------

###-----------------------------------------------------------------------------

require(pacman)

p_load(tidyverse,
       patchwork,
       ggpubr,
       gdata,
       conflicted)

select <- dplyr::select

filter <- dplyr::filter

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

## Read in Bluestem Data--------------------------------------------------------

read_csv("./Final_Code/Data/Leuthold et al., 2025 - Bluestem 10 Year Sampling.csv") %>% 
  select(-c(Full_name)) %>% 
  filter(Project == "Bluestem") -> Bluestem_Data

## Set initial values ----------------------------------------------------------

# Isotopic values --------------------------------------------------------------

T0_Litter_d13C <- 2112.96
T0_Litter_13C_Atom_Pct <- 3.38037
T0_Litter_d15N <- 10308.7
T0_Litter_15N_Atom_Pct <- 3.9917

# Mass values ------------------------------------------------------------------

T0_litter_Mass <- 18.4
T0_litter_C_Mass <- 8.1569
T0_litter_N_Mass <- 0.269717

###----------------------------------------------------------------------------- 
###-----------------------------------------------------------------------------

## Make the mixing models and calculate ----------------------------------------

Bluestem_Data %>%
  select(-Sample_No) %>% 
  mutate(Bulk_Density_g_cm3_Adj = case_when(Depth_cm == "0-2" ~ 0.702,
                                            Depth_cm == "2-5" ~ 0.923,
                                            TRUE ~ Bulk_Density_g_cm3)) %>% 
  mutate(Soil_Mass_g = Bulk_Density_g_cm3 * ((3.14 * 10 * 10) * (Bottom_Depth_cm - Top_Depth_cm)), # Calculate mass of soil in collar.
         C_Mass_g = Soil_Mass_g * (Bulk_C_Pct/100), # Calculate mass of C in collar.
         N_Mass_g = Soil_Mass_g * (Bulk_N_Pct/100)) %>% # Calculate mass of N in collar.
  pivot_wider(names_from = Cover_Treatment, # Pivot data for mixing model. 
              values_from = c(Soil_Mass_g, C_Mass_g, N_Mass_g,
                              Bulk_N_Pct, Bulk_15N_Atom_Pct,
                              Bulk_C_Pct, Bulk_d13C,
                              lPOM_N_Pct, lPOM_15N_Atom_Pct,
                              lPOM_C_Pct, lPOM_d13C,
                              Proportion_lPOM_Mass,
                              CHOM_N_Pct, CHOM_15N_Atom_Pct,
                              CHOM_C_Pct, CHOM_d13C,
                              Proportion_CHOM_Mass,
                              MAOM_N_Pct, MAOM_15N_Atom_Pct,
                              MAOM_C_Pct, MAOM_d13C,
                              Proportion_MAOM_Mass)) %>%
  mutate(Bulk_C_F_Value = ((Bulk_d13C_E - Bulk_d13C_B) / (T0_Litter_d13C - Bulk_d13C_B)), # Calculate C f-value (Enriched Soil - Bare Soil / Enriched litter - Bare Soil)
         Bulk_N_F_Value = ((Bulk_15N_Atom_Pct_E - Bulk_15N_Atom_Pct_B) / (T0_Litter_15N_Atom_Pct - Bulk_15N_Atom_Pct_B)), # Calculate N f-value
         Bulk_Litter_Derived_C_Mass = ((Soil_Mass_g_E * (Bulk_C_Pct_E/100)) * Bulk_C_F_Value), # Calculate the litter derived C in the bulk soil (g).
         Bulk_Litter_Derived_N_Mass = ((Soil_Mass_g_E * (Bulk_N_Pct_E/100)) * Bulk_N_F_Value), # Calculate the litter derived N in the bulk soil (g). 
         Bulk_Initial_Litter_C_Remaining = ((Bulk_Litter_Derived_C_Mass/T0_litter_C_Mass) * 100), # Calculate the proportion of the original litter remaining (%).
         Bulk_Initial_Litter_N_Remaining = ((Bulk_Litter_Derived_N_Mass/T0_litter_N_Mass) * 100), # Calculate the proportion of the original litter remaining (%).
         ###
         ### Repeat for the lPOM, CHAOM, and MAOM Fractions --------------------
         ###
         lPOM_C_F_Value = ((lPOM_d13C_E - lPOM_d13C_B) / (T0_Litter_d13C - lPOM_d13C_B)),
         lPOM_N_F_Value = ((lPOM_15N_Atom_Pct_E - lPOM_15N_Atom_Pct_B) / (T0_Litter_15N_Atom_Pct - lPOM_15N_Atom_Pct_B)),
         lPOM_Litter_Derived_C_Mass = (((Soil_Mass_g_E * (Proportion_lPOM_Mass_E/100)) * (lPOM_C_Pct_E/100)) * lPOM_C_F_Value),
         lPOM_Litter_Derived_N_Mass = (((Soil_Mass_g_E * (Proportion_lPOM_Mass_E/100)) * (lPOM_N_Pct_E/100)) * lPOM_N_F_Value),
         lPOM_Initial_Litter_C_Remaining = ((lPOM_Litter_Derived_C_Mass/T0_litter_C_Mass) * 100), 
         lPOM_Initial_Litter_N_Remaining = ((lPOM_Litter_Derived_N_Mass/T0_litter_N_Mass) * 100),
         ###
         CHAOM_C_F_Value = ((CHOM_d13C_E - CHOM_d13C_B) / (T0_Litter_d13C - CHOM_d13C_B)),
         CHAOM_N_F_Value = ((CHOM_15N_Atom_Pct_E - CHOM_15N_Atom_Pct_B) / (T0_Litter_15N_Atom_Pct - CHOM_15N_Atom_Pct_B)),
         CHAOM_Litter_Derived_C_Mass = (((Soil_Mass_g_E * (Proportion_CHOM_Mass_E/100)) * (CHOM_C_Pct_E/100)) * CHAOM_C_F_Value),
         CHAOM_Litter_Derived_N_Mass = (((Soil_Mass_g_E * (Proportion_CHOM_Mass_E/100)) * (CHOM_N_Pct_E/100)) * CHAOM_N_F_Value),
         CHAOM_Initial_Litter_C_Remaining = ((CHAOM_Litter_Derived_C_Mass/T0_litter_C_Mass) * 100), 
         CHAOM_Initial_Litter_N_Remaining = ((CHAOM_Litter_Derived_N_Mass/T0_litter_N_Mass) * 100),
         ###
         MAOM_C_F_Value = ((MAOM_d13C_E - MAOM_d13C_B) / (T0_Litter_d13C - MAOM_d13C_B)),
         MAOM_N_F_Value = ((MAOM_15N_Atom_Pct_E - MAOM_15N_Atom_Pct_B) / (T0_Litter_15N_Atom_Pct - MAOM_15N_Atom_Pct_B)),
         MAOM_Litter_Derived_C_Mass = (((Soil_Mass_g_E * (Proportion_MAOM_Mass_E/100)) * (MAOM_C_Pct_E/100)) * MAOM_C_F_Value),
         MAOM_Litter_Derived_N_Mass = (((Soil_Mass_g_E * (Proportion_MAOM_Mass_E/100)) * (MAOM_N_Pct_E/100)) * MAOM_N_F_Value),
         MAOM_Initial_Litter_C_Remaining = ((MAOM_Litter_Derived_C_Mass/T0_litter_C_Mass) * 100), 
         MAOM_Initial_Litter_N_Remaining = ((MAOM_Litter_Derived_N_Mass/T0_litter_N_Mass) * 100)) -> Bluestem_Derived_Data

### Summarize 10 year data -----------------------------------------------------

## Fumigated treatments --------------------------------------------------------

Bluestem_Derived_Data %>%
  filter(Fumigation_Treatment == "F") %>% 
  rowwise() %>% 
  group_by(Rep) %>% 
  summarize(Total_Litter_C = sum(Bulk_Litter_Derived_C_Mass),
            Fraction_Original_C = sum(Bulk_Litter_Derived_C_Mass/T0_litter_C_Mass) * 100,
            Total_Litter_N = sum(Bulk_Litter_Derived_N_Mass),
            Fraction_Original_N = sum(Bulk_Litter_Derived_N_Mass/T0_litter_N_Mass) * 100) %>% 
  summarise(Mean_Litter_C = mean(Total_Litter_C),
            Mean_Fraction_Original_C = mean(Fraction_Original_C),
            Mean_Litter_N = mean(Total_Litter_N),
            Mean_Fraction_Original_N = mean(Fraction_Original_N))

## Non-fumigated treatments ----------------------------------------------------

Bluestem_Derived_Data %>%
  filter(Fumigation_Treatment == "NF") %>% 
  rowwise() %>% 
  group_by(Rep) %>% 
  summarize(Total_Litter_C = sum(Bulk_Litter_Derived_C_Mass),
            Fraction_Original_C = sum(Bulk_Litter_Derived_C_Mass/T0_litter_C_Mass) * 100,
            Total_Litter_N = sum(Bulk_Litter_Derived_N_Mass),
            Fraction_Original_N = sum(Bulk_Litter_Derived_N_Mass/T0_litter_N_Mass) * 100) %>% 
  summarise(Mean_Litter_C = mean(Total_Litter_C),
            Mean_Fraction_Original_C = mean(Fraction_Original_C),
            Mean_Litter_N = mean(Total_Litter_N),
            Mean_Fraction_Original_N = mean(Fraction_Original_N))

### Add the original study data to create a time-series ------------------------

## Set this up to also include N

# Read in original data --------------------------------------------------------
read_csv("./Bluestem_OG_R_Data.csv") %>% 
  # Filter for the litter treatment --------------------------------------------
filter(`Litter type (B, E, NL)` == "E") %>% 
  select(Replicate,
         `Fauna Treatment`,
         Harvest,
         `Soil depth (cm)`,
         Fraction,
         `mg litter-derived C in total section`,
         `mg litter-derived N in total section`) %>% 
  # Rename the variables for ease ----------------------------------------------
rename(Replicate = Replicate,
       Fumigation_Treatment = `Fauna Treatment`,
       Harvest = Harvest,
       Depth_cm = `Soil depth (cm)`,
       Fraction  = Fraction,
       mg_Litter_C_per_Collar = `mg litter-derived C in total section`,
       mg_Litter_N_per_Collar = `mg litter-derived N in total section`) %>% 
  # Fix depths error from Excel ------------------------------------------------
mutate(Depth_cm = case_when(Depth_cm == "5-Feb" ~ "2-5",
                            Depth_cm == "10-May" ~ "5-10",
                            Depth_cm == "20-Oct" ~ "10-20",
                            TRUE ~ Depth_cm),
       Replicate = as.numeric(str_sub(Replicate, 2, 2)),
       # Convert to grams ---------------------------------------------------
       Litter_Derived_Carbon_Mass = mg_Litter_C_per_Collar / 1000,
       Litter_Derived_Nitrogen_Mass = mg_Litter_N_per_Collar / 1000) %>% 
  select(-mg_Litter_C_per_Collar,
         -mg_Litter_N_per_Collar) %>% 
  filter(Depth_cm %in% c("0-2", "2-5"),
         Fraction != "BULK") %>% 
  group_by(Fraction, Depth_cm, Harvest, Fumigation_Treatment, Replicate) %>% 
  # Get an average value per meter squared -------------------------------------
summarize(C_Mass = mean(Litter_Derived_Carbon_Mass)/0.0314,
          N_Mass = mean(Litter_Derived_Nitrogen_Mass)/0.0314) %>% 
  ungroup() %>% 
  # Rename fractions to current set of names -----------------------------------
mutate(Fraction = case_when(Fraction == "BULK" ~ "Bulk",
                            Fraction == "CLAY" ~ "MAOM",
                            Fraction == "SILT" ~ "MAOM",
                            Fraction == "LF" ~ "lPOM",
                            Fraction == "POM" ~ "CHAOM")) %>% 
  group_by(Harvest, Depth_cm, Fraction, Fumigation_Treatment, Replicate) %>% 
  # Add the MAOMs together ---------------------------------------------------
summarize(C_Mass = sum(C_Mass),
          N_Mass = sum(N_Mass)) %>% 
  ungroup() %>% 
  mutate(Fraction = fct_relevel(Fraction, "lPOM", "CHAOM", "MAOM")) %>% 
  group_by(Harvest, Fraction, Fumigation_Treatment, Replicate) %>%
  group_split() %>% 
  # Create a row for total C, adding the depths together -----------------------
map_dfr(~ add_row(.x, 
                  Harvest = unique(.x$Harvest), 
                  Fraction = unique(.x$Fraction),
                  Fumigation_Treatment = unique(.x$Fumigation_Treatment),
                  Replicate = unique(.x$Replicate),
                  Depth_cm = "Total", 
                  C_Mass = sum(.x$C_Mass),
                  N_Mass = sum(.x$N_Mass))) %>% 
  # Add data to the 10-year data -----------------------------------------------
bind_rows(., Bluestem_Derived_Data %>% 
            select(Rep, 
                   Fumigation_Treatment, 
                   Depth_cm, 
                   lPOM_Litter_Derived_C_Mass,
                   CHAOM_Litter_Derived_C_Mass,
                   MAOM_Litter_Derived_C_Mass) %>% 
            rename(Replicate = Rep) %>% 
            drop_na() %>% 
            # Manipulate the data in the way described above -----------------
          pivot_longer(cols = c(lPOM_Litter_Derived_C_Mass,
                                CHAOM_Litter_Derived_C_Mass,
                                MAOM_Litter_Derived_C_Mass),
                       names_to = "Fraction",
                       values_to = "C_Mass") %>% 
            mutate(Fraction = str_extract(Fraction, "^[^_]+")) %>%
            group_by(Fraction, Depth_cm, Fumigation_Treatment, Replicate) %>% 
            summarise(C_Mass = mean(C_Mass)/0.0314) %>%
            left_join(., Bluestem_Derived_Data %>%
                        select(Rep, 
                               Fumigation_Treatment, 
                               Depth_cm, 
                               lPOM_Litter_Derived_N_Mass,
                               CHAOM_Litter_Derived_N_Mass,
                               MAOM_Litter_Derived_N_Mass) %>% 
                        rename(Replicate = Rep) %>% 
                        drop_na() %>% 
                        pivot_longer(cols = c(lPOM_Litter_Derived_N_Mass,
                                              CHAOM_Litter_Derived_N_Mass,
                                              MAOM_Litter_Derived_N_Mass),
                                     names_to = "Fraction",
                                     values_to = "N_Mass") %>% 
                        mutate(Fraction = str_extract(Fraction, "^[^_]+")) %>% 
                        group_by(Fraction, Depth_cm, Fumigation_Treatment, Replicate) %>% 
                        summarise(N_Mass = mean(N_Mass)/0.0314),
                      by = c("Fumigation_Treatment", "Depth_cm", "Fraction", "Replicate")) %>% 
            ungroup() %>% 
            group_by(Fraction, Fumigation_Treatment, Replicate) %>% 
            group_split() %>% 
            map_dfr(~ add_row(.x, 
                              Fraction = unique(.x$Fraction),
                              Fumigation_Treatment = unique(.x$Fumigation_Treatment),
                              Replicate = unique(.x$Replicate),
                              Depth_cm = "Total", 
                              C_Mass = sum(.x$C_Mass),
                              N_Mass = sum(.x$N_Mass))) %>%  
            mutate(Fraction = fct_relevel(Fraction, "lPOM", "CHAOM", "MAOM"),
                   Harvest = "H7") %>% 
            select(Harvest,
                   Depth_cm,
                   Fraction,
                   Fumigation_Treatment,
                   Replicate,
                   C_Mass,
                   N_Mass)) %>% 
  # Change the name of the harvests --------------------------------------------
mutate(Harvest = case_when(Harvest == "H1" ~ "7 Months",
                           Harvest == "H2" ~ "12 Months",
                           Harvest == "H3" ~ "18 Months",
                           Harvest == "H4" ~ "24 Months",
                           Harvest == "H6" ~ "36 Months",
                           Harvest == "H7" ~ "120 Months"),
       Harvest = fct_relevel(Harvest, "7 Months", "12 Months", "18 Months", "24 Months", "36 Months", "120 Months"),
       Harvest_Num = case_when(Harvest == "7 Months" ~ 7,
                               Harvest == "12 Months" ~ 12,
                               Harvest == "18 Months" ~ 18,
                               Harvest == "24 Months" ~ 24,
                               Harvest == "36 Months" ~ 36,
                               Harvest == "120 Months" ~ 120))  -> Litter_Time_Series_Data



