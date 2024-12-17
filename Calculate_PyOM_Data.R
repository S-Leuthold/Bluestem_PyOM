### Script for calculating the litter and char derived C and N from the Black C project.
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

###----------------------------------------------------------------------------
###----------------------------------------------------------------------------

### Read in Black C Data-------------------------------------------------------

read_csv("./Final_Code/Data/Leuthold et al., 2025 - Bluestem 10 Year Sampling.csv") %>% 
  select(-c(Full_name, 
            Fumigation_Treatment)) %>% 
  filter(Project == "BlackC") -> BlackC_Data

### Set initial values ---------------------------------------------------------

## Isotopic values -------------------------------------------------------------

# Char -------------------------------------------------------------------------

T0_Char_d13C <- 3097
T0_Char_Pct_C <- 55.31
T0_Char_d15N <- 17729
T0_Char_15N_Atom_Pct <- ((((T0_Char_d15N / 1000) + 1) * 0.0036765) / ((((T0_Char_d15N / 1000) + 1) * 0.0036765) + 1)) * 100
T0_Char_Pct_N <- 1.685

# Litter -----------------------------------------------------------------------

T0_Litter_d13C <- 3045
T0_Litter_Pct_C  <- 41.92
T0_Litter_d15N <- 19253
T0_Litter_15N_Atom_Pct <- ((((T0_Litter_d15N / 1000) + 1) * 0.0036765) / ((((T0_Litter_d15N / 1000) + 1) * 0.0036765) + 1)) * 100
T0_Litter_Pct_N <- 1.06

## Mass Values -----------------------------------------------------------------

# Char -------------------------------------------------------------------------

T0_Char_Mass <- 5.92
T0_Char_C_Mass <- T0_Char_Mass * (T0_Char_Pct_C/100)
T0_Char_N_Mass <- T0_Char_Mass * (T0_Char_Pct_N/100)

# Litter -----------------------------------------------------------------------

T0_Litter_Mass  <- 17
T0_Litter_C_Mass <- T0_Litter_Mass * (T0_Litter_Pct_C/100)
T0_Litter_N_Mass <- T0_Litter_Mass * (T0_Litter_Pct_N/100)

###----------------------------------------------------------------------------- 
###-----------------------------------------------------------------------------

BlackC_Data %>%
  select(-Sample_No) %>% 
  mutate(Soil_Mass_g = Bulk_Density_g_cm3 * ((3.14 * 10 * 10) * (Bottom_Depth_cm - Top_Depth_cm)),
         C_Mass_g = Soil_Mass_g * (Bulk_C_Pct/100),
         N_Mass_g = Soil_Mass_g * (Bulk_N_Pct/100)) %>%
  pivot_wider(names_from = Cover_Treatment,
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
  ## First, the litter. --------------------------------------------------------

mutate(Bulk_C_F_Value = ((Bulk_d13C_L - Bulk_d13C_BC_B) / (T0_Litter_d13C - Bulk_d13C_BC_B)),
       Bulk_N_F_Value = ((Bulk_15N_Atom_Pct_L - Bulk_15N_Atom_Pct_BC_B) / (T0_Litter_15N_Atom_Pct - Bulk_15N_Atom_Pct_BC_B)),
       Bulk_Litter_Derived_C_Mass = ((Soil_Mass_g_L * (Bulk_C_Pct_L/100)) * Bulk_C_F_Value),
       Bulk_Litter_Derived_N_Mass = ((Soil_Mass_g_L * (Bulk_N_Pct_L/100)) * Bulk_N_F_Value),
       Bulk_Initial_Litter_C_Remaining = ((Bulk_Litter_Derived_C_Mass/T0_Litter_C_Mass) * 100), 
       Bulk_Initial_Litter_N_Remaining = ((Bulk_Litter_Derived_N_Mass/T0_Litter_N_Mass) * 100),
       ###
       lPOM_C_F_Value = ((lPOM_d13C_L - lPOM_d13C_BC_B) / (T0_Litter_d13C - lPOM_d13C_BC_B)),
       lPOM_N_F_Value = ((lPOM_15N_Atom_Pct_L - lPOM_15N_Atom_Pct_BC_B) / (T0_Litter_15N_Atom_Pct - lPOM_15N_Atom_Pct_BC_B)),
       lPOM_Litter_Derived_C_Mass = (((Soil_Mass_g_L * (Proportion_lPOM_Mass_L/100)) * (lPOM_C_Pct_L/100)) * lPOM_C_F_Value),
       lPOM_Litter_Derived_N_Mass = (((Soil_Mass_g_L * (Proportion_lPOM_Mass_L/100)) * (lPOM_N_Pct_L/100)) * lPOM_N_F_Value),
       lPOM_Initial_Litter_C_Remaining = ((lPOM_Litter_Derived_C_Mass/T0_Litter_C_Mass) * 100), 
       lPOM_Initial_Litter_N_Remaining = ((lPOM_Litter_Derived_N_Mass/T0_Litter_N_Mass) * 100),
       ###
       CHAOM_C_F_Value = ((CHOM_d13C_L - CHOM_d13C_BC_B) / (T0_Litter_d13C - CHOM_d13C_BC_B)),
       CHAOM_N_F_Value = ((CHOM_15N_Atom_Pct_L - CHOM_15N_Atom_Pct_BC_B) / (T0_Litter_15N_Atom_Pct - CHOM_15N_Atom_Pct_BC_B)),
       CHAOM_Litter_Derived_C_Mass = (((Soil_Mass_g_L * (Proportion_CHOM_Mass_L/100)) * (CHOM_C_Pct_L/100)) * CHAOM_C_F_Value),
       CHAOM_Litter_Derived_N_Mass = (((Soil_Mass_g_L * (Proportion_CHOM_Mass_L/100)) * (CHOM_N_Pct_L/100)) * CHAOM_N_F_Value),
       CHAOM_Initial_Litter_C_Remaining = ((CHAOM_Litter_Derived_C_Mass/T0_Litter_C_Mass) * 100), 
       CHAOM_Initial_Litter_N_Remaining = ((CHAOM_Litter_Derived_N_Mass/T0_Litter_N_Mass) * 100),
       ###
       MAOM_C_F_Value = ((MAOM_d13C_L - MAOM_d13C_BC_B) / (T0_Litter_d13C - MAOM_d13C_BC_B)),
       MAOM_N_F_Value = ((MAOM_15N_Atom_Pct_L - MAOM_15N_Atom_Pct_BC_B) / (T0_Litter_15N_Atom_Pct - MAOM_15N_Atom_Pct_BC_B)),
       MAOM_Litter_Derived_C_Mass = (((Soil_Mass_g_L * (Proportion_MAOM_Mass_L/100)) * (MAOM_C_Pct_L/100)) * MAOM_C_F_Value),
       MAOM_Litter_Derived_N_Mass = (((Soil_Mass_g_L * (Proportion_MAOM_Mass_L/100)) * (MAOM_N_Pct_L/100)) * MAOM_N_F_Value),
       MAOM_Initial_Litter_C_Remaining = ((MAOM_Litter_Derived_C_Mass/T0_Litter_C_Mass) * 100), 
       MAOM_Initial_Litter_N_Remaining = ((MAOM_Litter_Derived_N_Mass/T0_Litter_N_Mass) * 100)) %>% 
  
  ## Then, the pyrolized material ----------------------------------------------

mutate(Bulk_C_F_Value = ((Bulk_d13C_C - Bulk_d13C_BC_B) / (T0_Char_d13C - Bulk_d13C_BC_B)),
       Bulk_N_F_Value = ((Bulk_15N_Atom_Pct_C - Bulk_15N_Atom_Pct_BC_B) / (T0_Char_15N_Atom_Pct - Bulk_15N_Atom_Pct_BC_B)),
       Bulk_Char_Derived_C_Mass = ((Soil_Mass_g_C * (Bulk_C_Pct_C/100)) * Bulk_C_F_Value),
       Bulk_Char_Derived_N_Mass = ((Soil_Mass_g_C * (Bulk_N_Pct_C/100)) * Bulk_N_F_Value),
       Bulk_Initial_Char_C_Remaining = ((Bulk_Char_Derived_C_Mass/T0_Char_C_Mass) * 100), 
       Bulk_Initial_Char_N_Remaining = ((Bulk_Char_Derived_N_Mass/T0_Char_N_Mass) * 100),
       ###
       lPOM_C_F_Value = ((lPOM_d13C_C - lPOM_d13C_BC_B) / (T0_Char_d13C - lPOM_d13C_BC_B)),
       lPOM_N_F_Value = ((lPOM_15N_Atom_Pct_C - lPOM_15N_Atom_Pct_BC_B) / (T0_Char_15N_Atom_Pct - lPOM_15N_Atom_Pct_BC_B)),
       lPOM_Char_Derived_C_Mass = (((Soil_Mass_g_C * (Proportion_lPOM_Mass_C/100)) * (lPOM_C_Pct_C/100)) * lPOM_C_F_Value),
       lPOM_Char_Derived_N_Mass = (((Soil_Mass_g_C * (Proportion_lPOM_Mass_C/100)) * (lPOM_N_Pct_C/100)) * lPOM_N_F_Value),
       lPOM_Initial_Char_C_Remaining = ((lPOM_Char_Derived_C_Mass/T0_Char_C_Mass) * 100), 
       lPOM_Initial_Char_N_Remaining = ((lPOM_Char_Derived_N_Mass/T0_Char_N_Mass) * 100),
       ###
       CHAOM_C_F_Value = ((CHOM_d13C_C - CHOM_d13C_BC_B) / (T0_Char_d13C - CHOM_d13C_BC_B)),
       CHAOM_N_F_Value = ((CHOM_15N_Atom_Pct_C - CHOM_15N_Atom_Pct_BC_B) / (T0_Char_15N_Atom_Pct - CHOM_15N_Atom_Pct_BC_B)),
       CHAOM_Char_Derived_C_Mass = (((Soil_Mass_g_C * (Proportion_CHOM_Mass_C/100)) * (CHOM_C_Pct_C/100)) * CHAOM_C_F_Value),
       CHAOM_Char_Derived_N_Mass = (((Soil_Mass_g_C * (Proportion_CHOM_Mass_C/100)) * (CHOM_N_Pct_C/100)) * CHAOM_N_F_Value),
       CHAOM_Initial_Char_C_Remaining = ((CHAOM_Char_Derived_C_Mass/T0_Char_C_Mass) * 100), 
       CHAOM_Initial_Char_N_Remaining = ((CHAOM_Char_Derived_N_Mass/T0_Char_N_Mass) * 100),
       ###
       MAOM_C_F_Value = ((MAOM_d13C_C - MAOM_d13C_BC_B) / (T0_Char_d13C - MAOM_d13C_BC_B)),
       MAOM_N_F_Value = ((MAOM_15N_Atom_Pct_C - MAOM_15N_Atom_Pct_BC_B) / (T0_Char_15N_Atom_Pct - MAOM_15N_Atom_Pct_BC_B)),
       MAOM_Char_Derived_C_Mass = (((Soil_Mass_g_C * (Proportion_MAOM_Mass_C/100)) * (MAOM_C_Pct_C/100)) * MAOM_C_F_Value),
       MAOM_Char_Derived_N_Mass = (((Soil_Mass_g_C * (Proportion_MAOM_Mass_C/100)) * (MAOM_N_Pct_C/100)) * MAOM_N_F_Value),
       MAOM_Initial_Char_C_Remaining = ((MAOM_Char_Derived_C_Mass/T0_Char_C_Mass) * 100), 
       MAOM_Initial_Char_N_Remaining = ((MAOM_Char_Derived_N_Mass/T0_Char_N_Mass) * 100))  -> BlackC_Derived_Data

## Look at recoveries ----------------------------------------------------------

BlackC_Derived_Data %>% 
  rowwise() %>% 
  group_by(Rep) %>% 
  summarize(Total_Litter_C = sum(Bulk_Litter_Derived_C_Mass),
            Fraction_Original_Litter_C = sum(Bulk_Litter_Derived_C_Mass/T0_Litter_C_Mass),
            Total_Char_C = sum(Bulk_Char_Derived_C_Mass),
            Fraction_Original_Char_C = sum(Bulk_Char_Derived_C_Mass/T0_Char_C_Mass)) %>% 
  ungroup() %>% 
  summarise(Mean_Char_Remaining = mean(Fraction_Original_Char_C),
            Mean_Litter_Remaining = mean(Fraction_Original_Litter_C))

BlackC_Derived_Data %>% 
  rowwise() %>% 
  group_by(Rep) %>% 
  summarize(Total_Litter_N = sum(Bulk_Litter_Derived_N_Mass),
            Fraction_Original_Litter_N = sum(Bulk_Litter_Derived_N_Mass/T0_Litter_N_Mass),
            Total_Char_N = sum(Bulk_Char_Derived_N_Mass),
            Fraction_Original_Char_N = sum(Bulk_Char_Derived_N_Mass/T0_Char_N_Mass)) %>% 
  ungroup() %>% 
  summarise(Mean_Char_Remaining = mean(Fraction_Original_Char_N),
            Mean_Litter_Remaining = mean(Fraction_Original_Litter_N))




# Read in original data --------------------------------------------------------
read_csv("./Final_Code/Data/Leuthold et al., 2025 - PyOM Original Data.csv") %>% 
  # Filter for the pyrolized treatment --------------------------------------------
  filter(COVER == "C") %>% 
  select(REP,
         COVER,
         TIME,
         DEPTH,
         FRACTION,
         Total_derived_mg_C,
         Total_derived_mg_N) %>% 
  # Rename the variables for ease ----------------------------------------------
  rename(Replicate = REP,
         Cover_Treatment = COVER,
         Harvest = TIME,
         Depth_cm = DEPTH,
         Fraction  = FRACTION,
         mg_Litter_C_per_Collar = Total_derived_mg_C,
         mg_Litter_N_per_Collar = Total_derived_mg_N) %>% 
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
         Fraction != "bulk") %>% 
  group_by(Fraction, Depth_cm, Harvest,Cover_Treatment, Replicate) %>% 
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
  group_by(Harvest, Depth_cm, Fraction, Cover_Treatment, Replicate) %>% 
  # Add the MAOMs together ---------------------------------------------------
  summarize(C_Mass = sum(C_Mass, na.rm = TRUE),
            N_Mass = sum(N_Mass, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Fraction = fct_relevel(Fraction, "lPOM", "CHAOM", "MAOM")) %>% 
  group_by(Harvest, Fraction, Cover_Treatment, Replicate) %>%
  group_split() %>% 
  # Create a row for total C, adding the depths together -----------------------
  map_dfr(~ add_row(.x, 
                    Harvest = unique(.x$Harvest), 
                    Fraction = unique(.x$Fraction),
                    Cover_Treatment = unique(.x$Cover_Treatment),
                    Replicate = unique(.x$Replicate),
                    Depth_cm = "Total", 
                    C_Mass = sum(.x$C_Mass),
                    N_Mass = sum(.x$N_Mass))) %>% 
  # Add data to the 10-year data -----------------------------------------------
  bind_rows(., BlackC_Derived_Data %>% 
            select(Rep, 
                   Depth_cm, 
                   lPOM_Char_Derived_C_Mass,
                   CHAOM_Char_Derived_C_Mass,
                   MAOM_Char_Derived_C_Mass) %>% 
            rename(Replicate = Rep) %>% 
            drop_na() %>% 
            # Manipulate the data in the way described above -----------------
            pivot_longer(cols = c(lPOM_Char_Derived_C_Mass,
                                  CHAOM_Char_Derived_C_Mass,
                                  MAOM_Char_Derived_C_Mass),
                         names_to = "Fraction",
                         values_to = "C_Mass") %>% 
            mutate(Fraction = str_extract(Fraction, "^[^_]+")) %>%
            group_by(Fraction, Depth_cm, Replicate) %>% 
            summarise(C_Mass = mean(C_Mass)/0.0314) %>%
            left_join(., BlackC_Derived_Data %>%
                        select(Rep, 
                               Depth_cm, 
                               lPOM_Char_Derived_N_Mass,
                               CHAOM_Char_Derived_N_Mass,
                               MAOM_Char_Derived_N_Mass) %>% 
                        rename(Replicate = Rep) %>% 
                        drop_na() %>% 
                        pivot_longer(cols = c(lPOM_Char_Derived_N_Mass,
                                              CHAOM_Char_Derived_N_Mass,
                                              MAOM_Char_Derived_N_Mass),
                                     names_to = "Fraction",
                                     values_to = "N_Mass") %>% 
                        mutate(Fraction = str_extract(Fraction, "^[^_]+")) %>% 
                        group_by(Fraction, Depth_cm, Replicate) %>% 
                        summarise(N_Mass = mean(N_Mass)/0.0314),
                      by = c("Depth_cm", "Fraction", "Replicate")) %>% 
            mutate(Cover_Treatment = "C") %>% 
            ungroup() %>% 
            group_by(Fraction, Cover_Treatment, Replicate) %>% 
            group_split() %>% 
            map_dfr(~ add_row(.x, 
                              Fraction = unique(.x$Fraction),
                              Cover_Treatment = (.x$Cover_Treatment),
                              Replicate = unique(.x$Replicate),
                              Depth_cm = "Total", 
                              C_Mass = sum(.x$C_Mass),
                              N_Mass = sum(.x$N_Mass))) %>%  
            mutate(Fraction = fct_relevel(Fraction, "lPOM", "CHAOM", "MAOM"),
                   Harvest = "108mts") %>% 
            select(Harvest,
                   Depth_cm,
                   Fraction,
                   Cover_Treatment,
                   Replicate,
                   C_Mass,
                   N_Mass) %>% 
              unique(.)) %>% 
  # Change the name of the harvests --------------------------------------------
  mutate(Harvest = case_when(Harvest == "4mts" ~ "4 Months",
                             Harvest == "11mts" ~ "11 Months",
                             Harvest == "108mts" ~ "108 Months"),
         Harvest = fct_relevel(Harvest, "4 Months", "11 Months", "108 Months"),
         Harvest_Num = case_when(Harvest == "4 Months" ~ 4,
                                 Harvest == "11 Months" ~ 11,
                                 Harvest == "108 Months" ~ 108)) %>% 
  filter(Replicate != 3 & Harvest != 108) -> BlackC_Time_Series_Data


