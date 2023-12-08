# Function to wrangle data

install.packages("tidyverse")
library(tidyverse)
library(tidyr)
library(stringi)
library(MoMAColors)

path = "/Users/ssackey/Documents/Organoid_Project/6_donor_Tonsil_infection_in_triplicate_timecourse.csv"

try(if (grepl(".xlsx", path)) stop("File needs to be a .csv"))
raw_data = read.csv("/Users/ssackey/Documents/Organoid_Project/6_donor_Tonsil_infection_in_triplicate_timecourse.csv")
View(raw_data)

#Check if the column names have "+" or "-", if they do then print out a warning
if(grepl(colnames() == "+" | colnames() == "-")) {
  print("Need to clean column names")
}

#If the columns have "+", "-", clean the column names by replacing the symbol with "pos", "neg" or "_'
Clean_col_names = ifelse(grepl(colnames() == "+" | colnames() == "-"),
                         str_replace_all(data, pattern=c('\\+', 'Nonneg', '-'), 
                                         replacement=c('pos', '_', 'neg'), stop()))

#Change MOI = 0 to Uninfected
raw_data$MOI[raw_data$MOI == "0"] = "Uninfected"
View(raw_data)

subset_cal = raw_data |> mutate(Uninfected_CD4Total_CM = (CD4Total_p24Total_CM / CD4Total_p24Total_Non_TFH_cells)*100,
                                Uninfected_CD4Total_EM = (CD4Total_p24Total_EM / CD4Total_p24Total_Non_TFH_cells)*100,
                                Uninfected_CD4Total_Naive = (CD4Total_p24Total_Naive / CD4Total_p24Total_Non_TFH_cells)*100,
                                Uninfected_CD4Total_TEMRA = (CD4Total_p24Total_TEMRA / CD4Total_p24Total_Non_TFH_cells)*100,
                                CD4_p24neg_CM = (CD4Total_p24neg_CM / CD4Total_p24neg_Non_TFH_cells)*100,
                                CD4_p24neg_EM = (CD4Total_p24neg_EM / CD4Total_p24neg_Non_TFH_cells)*100,
                                CD4_p24neg_Naive = (CD4Total_p24neg_Naive / CD4Total_p24neg_Non_TFH_cells)*100,
                                CD4_p24neg_TEMRA = ( CD4Total_p24neg_TEMRA/ CD4Total_p24neg_Non_TFH_cells)*100,
                                CD4_p24pos_CM = (CD4Total_p24pos_CM / CD4Total_p24pos_Non_TFH_cells)*100,
                                CD4_p24pos_EM = (CD4Total_p24pos_EM / CD4Total_p24pos_Non_TFH_cells)*100,
                                CD4_p24pos_Naive = (CD4Total_p24pos_Naive / CD4Total_p24pos_Non_TFH_cells)*100,
                                CD4_p24pos_TEMRA = (CD4Total_p24pos_TEMRA / CD4Total_p24pos_Non_TFH_cells)*100,
                                CD8_CM = (CD8pos_CM / CD8pos_Non_Follicular_cells)*100,
                                CD8_EM = (CD8pos_EM / CD8pos_Non_Follicular_cells)*100,
                                CD8_Naive = (CD8pos_Naive / CD8pos_Non_Follicular_cells)*100,
                                CD8_TEMRA = (CD8pos_TEMRA / CD8pos_Non_Follicular_cells)*100, .keep = "unused")
colnames(subset_cal)

#Making a vector with all the columns I want to use for CD4 cells analysis
CD4_col_of_interest = c("CD4Total_p24Total_TFH_cells", "Uninfected_CD4Total_CM", "Uninfected_CD4Total_EM", "Uninfected_CD4Total_Naive", "Uninfected_CD4Total_TEMRA", "CD4Total_p24neg_TFH_cells", "CD4_p24neg_CM", "CD4_p24neg_EM", "CD4_p24neg_Naive", "CD4_p24neg_TEMRA", "CD4Total_p24pos_TFH_cells", "CD4_p24pos_CM", "CD4_p24pos_EM","CD4_p24pos_Naive", "CD4_p24pos_TEMRA")

#Select the column I indicated in my vector from my data frame subset_cal
subset_cal |> select(all_of(CD4_col_of_interest))

#Making a vector with all the metadate columns I want to use for CD4 cells analysis
Meta_data = c("Donor_ID", "Timepoint", "Replicate", "MOI", "Maturity")
subset_cal |> select(all_of(Meta_data))

#Make a new data frame to combine all the columns I want in the new data frame with empty column I will later populate
CD4_df = subset_cal[,c(Meta_data, CD4_col_of_interest)]
CD4_datalong = gather(CD4_df, cell_type, freq, CD4_col_of_interest)
CD4_datalong

CD4_datalong$p24_type = ifelse(grepl("p24neg", CD4_datalong$cell_type), "Infected_p24neg", "Uninfected_CD4Total")
CD4_datalong$p24_type = ifelse(grepl("p24pos", CD4_datalong$cell_type), "Infected_p24pos", CD4_datalong$p24_type)
CD4_datalong$MOI |> as.factor()
CD4_datalong

CD4_datalong$cell_type2 = ifelse(grepl("CM",CD4_datalong$cell_type),"CM", "EM")
CD4_datalong$cell_type2 = ifelse(grepl("Naive",CD4_datalong$cell_type),"Naive", CD4_datalong$cell_type2)
CD4_datalong$cell_type2 = ifelse(grepl("TEMRA",CD4_datalong$cell_type),"TEMRA", CD4_datalong$cell_type2)
CD4_datalong$cell_type2 = ifelse(grepl("TFH",CD4_datalong$cell_type),"TFH", CD4_datalong$cell_type2)
CD4_datalong

#Making a vector with all the columns I want to use for CD8 cells analysis
CD8_col_of_interest = c("CD8pos_Follicular_cells", "CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA")

#Select the column I indicated in my vector from my data frame subset_cal
subset_cal |> select(all_of(CD8_col_of_interest))

#Making a vector with all the metadate columns I want to use for CD8 cells analysis
Meta_data = c("Donor_ID", "Timepoint", "Replicate", "MOI", "Maturity")
subset_cal |> select(all_of(Meta_data))

#Make a new data frame to combine all the columns I want in the new data frame with empty column I will later populate
CD8_df = subset_cal[,c(Meta_data, CD8_col_of_interest)]
CD8_datalong = gather(CD8_df, cell_type, freq, CD8_col_of_interest)
CD8_datalong

CD8_datalong$p24_type = ifelse(grepl("p24neg", CD8_datalong$cell_type), "Infected_p24neg", "Uninfected_CD4Total")
CD8_datalong$p24_type = ifelse(grepl("p24pos", CD8_datalong$cell_type), "Infected_p24pos", CD8_datalong$p24_type)
CD8_datalong$MOI |> as.factor()
CD8_datalong

CD8_datalong$cell_type2 = ifelse(grepl("CM",CD8_datalong$cell_type),"CM", "EM")
CD8_datalong$cell_type2 = ifelse(grepl("Naive",CD8_datalong$cell_type),"Naive", CD8_datalong$cell_type2)
CD8_datalong$cell_type2 = ifelse(grepl("TEMRA",CD8_datalong$cell_type),"TEMRA", CD8_datalong$cell_type2)
CD8_datalong$cell_type2 = ifelse(grepl("Follicular",CD8_datalong$cell_type),"Follicular_cells", CD8_datalong$cell_type2)
CD8_datalong

#Making a vector with all the columns I want to use for NK cells analysis
NK_col_of_interest = c("Follicular_NK_cells","Non_Follicular_NKs" )

#Select the column I indicated in my vector from my data frame subset_cal
subset_cal |> select(all_of(NK_col_of_interest))

#Making a vector with all the metadate columns I want to use for NK cells analysis
Meta_data = c("Donor_ID", "Timepoint", "Replicate", "MOI", "Maturity")
subset_cal |> select(all_of(Meta_data))

#Make a new data frame to combine all the columns I want in the new data frame with empty column I will later populate
NK_df = subset_cal[,c(Meta_data, NK_col_of_interest)]
NK_datalong = gather(NK_df, cell_type, freq, NK_col_of_interest)
NK_datalong

NK_datalong$p24_type = ifelse(grepl("p24neg", NK_datalong$cell_type), "Infected_p24neg", "Uninfected_CD4Total")
NK_datalong$p24_type = ifelse(grepl("p24pos", NK_datalong$cell_type), "Infected_p24pos", NK_datalong$p24_type)
NK_datalong$MOI |> as.factor()
NK_datalong

NK_datalong$cell_type2 = ifelse(grepl("Follicular_NK_cells",NK_datalong$cell_type),"Follicular_NK_cells", "Non_Follicular_NKs") 
NK_datalong

