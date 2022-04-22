#load libraries
library(tidyverse)

#define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
Output_File_Path <- choose.dir(caption = "Select (open on Mac) Pulse_Chase_Analysis folder, where compiled data will be stored")

#define type of cells used in experiment - either "mESCs" or "NPCs"
Cell_Type <- "mESCs"

#define saving path
save_path <- paste(Output_File_Path, Cell_Type, "Molecule_Count", "All_Cell_Lines", sep="/")

#create new folder to save individual data in
dir.create(save_path)

#file to tidy
file <- choose.files() %>% 
  read_csv

#create new columns in table
file_new_cols <- file %>% 
  mutate(Pulse = NA,
         Data_Set = NA,
         Time = NA,
         Cloud_Name = NA) %>% 
  transmute(Cell_Line = Cell_Line,
            Data_Set = Data_Set, 
            Phase = Phase, 
            Pulse = Pulse,
            Time = Time, 
            Cloud_Name = Cloud_Name,
            Total_Centroids = Total_Centroids)

#cell line names
names <- file_new_cols %>% 
   select(Cell_Line) %>% 
   distinct()

print(names)

#change names if required
new_names <- file_new_cols %>% 
   mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT",
                                Cell_Line == "Mettl3_dTAG" ~ "Mettl3_dTAG",
                                Cell_Line == "Transgenic" ~ "Transgenic",
                                Cell_Line == "Ciz1_KO" ~ "Ciz1_KO",
                                Cell_Line == "SPEN_RRM_del" ~ "SPEN_RRM_del",
                                Cell_Line == "SPOC_mut" ~ "SPOC_mut",
                                Cell_Line == "CIZ1_KO" ~ "Ciz1_KO",
                                Cell_Line == "SPEN_KO" ~ "SPEN_RRM_del",
                                Cell_Line == "SPOC_KO" ~ "SPOC_mut",)) 

#list of cell lines
cell_lines <- c("WT", "Ciz1_KO", "SPEN_RRM_del", "SPOC_mut", "Mettl3_dTAG", 
                "Transgenic")

#separate out data
 for (cell in cell_lines){
   individual <- new_names %>% 
     filter(Cell_Line == cell)
   
   file_name <- paste(cell, "Total_Molecule_Count_Compile.csv", sep="_")
   
   write_csv(individual, paste(save_path, file_name, sep="/"))
 }

#check no data has been lost

files <- list.files(path = save_path, pattern = ".csv", full.names  = TRUE)

file_lists <- lapply(files, read_csv, col_types = "ccccncn")

#concatenates all the individual tibbles into 1 tibble
compile_files <- bind_rows(file_lists) 