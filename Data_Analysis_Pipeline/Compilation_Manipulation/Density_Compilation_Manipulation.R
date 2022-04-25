####################################################################################################################
#SCRIPT COMPILES ALL DATA REQUIRED TO CALACULATE DENSITY
#please contact Holly Roach at hmroach@hotmail.co.uk if you have any questions
####################################################################################################################
#OUTPUTS: New_Line_Name_Density_Compile.csv 
#compiles all the Median_nn_dist rows of data for all the C1-C1.csv and C2-C2.csv files into 1 table for the new cell line (New_Line_Name_Density_Compile.csv)
#makes sure tidyverse, stringr and tcltk packaged have been installed - install.packages("package_name")

#load libraries
library(tidyverse)
library(stringr)
library(tcltk)

#create function to replace choose.dir so that it is compatible with mac OS
choose_dir <- function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

################################################################################
#USER INPUT REQUIRED

#set name of new cell line - should be spelled the same as folder name within the directory
#make sure not to use any spaces or symbols other than an underscore (_), dash (-), or full stop (.)
#eg "Mettl3_dTAG" or "SPEN_RRM_del" or "Ciz1_KO"
New_Line_Name <- "Test"           

#define type of cells used in experiment - either "mESCs" or "NPCs"
Cell_Type <- "mESCs"

#define time points used in both datasets
#make sure that the lists include all the possible time-points for your data
#it is ok if the lists contain extra time-points, as these lists are a superset of time points
#it doesn't matter if initiation/maintenance have different time points
Dynamic_Time <- c(10, 20, 30, 40, 50, 60)
Turnover_Time <- c(0, 60, 80, 100, 120, 140, 160, 180, 200, 220)


################################################################################
#OTHER INPUTS
#the below inputs should not need changing

#creates list of datasets being used for compile
Data_Set_List <- c("nascent_Xist_dynamics",
                   "Xist_turnover_on_chromatin")

#creates list of phases being used for compile
Phase_List <- c("Initiation", "Maintenance")

#define file path to where Watershed Algorithm "results" are found - needs to be located in documents due to long file path names
Input_File_Path <- choose_dir(caption = "Select (open on Mac) Watershed_Algorithm_Results folder, where raw data is found")

#define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
Output_File_Path <- choose_dir(caption = "Select (open on Mac) Pulse_Chase_Analysis folder, where compiled data will be stored")


################################################################################
#VALIDATE INPUTS

#checks inputted name of new cell line - compares name to REGEX
if (grepl("\\W", New_Line_Name)) {
  stop(paste("Invalid character in name of new cell line (", New_Line_Name, " - New_Line_Name) - must not contain any spaces or symbols", sep=""))
}

#checks that Watershed_Algorithm_Results folder exists for new cell line
if (!dir.exists(paste(Input_File_Path, Cell_Type, New_Line_Name, sep="/"))) {
  stop(paste("Results folder for new cell line (", New_Line_Name, " - New_Line_Name) does not exist in Watershed_Algorithm_Results folder", sep=""))
}

#checks name of dataset
if (!("nascent_Xist_dynamics" %in% Data_Set_List |"Xist_turnover_on_chromatin" %in% Data_Set_List)) {
  stop("Invalid names entered in Data_Set_List - must contain either nascent_Xist_dynamics or Xist_turnover_on_chromatin or both")
}

#checks name of phase
if (!("Initiation" %in% Phase_List |"Maintenance" %in% Phase_List)) {
  stop("Invalid names entered in Phase_List - must contain either Initiation or Maintenance or both")
}

#checks name of cell type
if (!(Cell_Type == "mESCs" | Cell_Type == "NPCs")) {
  stop("Invalid name entered in Cell_Type - must be mESCs or NPCs")
}

#checks correct folder has been selected for input file path
if (!(str_sub(Input_File_Path, -27, -1) == "Watershed_Algorithm_Results")) {
  stop("Watershed_Algorithm_Results folder not selected for Input_File_Path")
}

#checks correct folder has been selected for output file path
if (!(str_sub(Output_File_Path, -20, -1) == "Pulse_Chase_Analysis")) {
  stop("Pulse_Chase_Analysis folder not selected for Output_File_Path")
}

#checks that user has permissions to read and write files in Watershed_Algorithm_Results folder
if (!(file.access(Input_File_Path, 2) == 0 && file.access(Input_File_Path, 4) == 0 )) {
  stop("Invalid file permissions for Watershed_Algorithm_Results folder - check in folders properties that you have permission to read&write")
}

#checks that user has permissions to read and write files in Pulse_Chase_Analysis folder
if (!(file.access(Output_File_Path, 2) == 0 && file.access(Output_File_Path, 4) == 0 )) {
  stop("Invalid file permissions for Pulse_Chase_Analysis folder - check in folders properties that you have permission to read&write")
}


################################################################################
#creates table which will store table for new cell line data

New_Cell_Line <- tibble(Cell_Line = character(),
                        Data_Set = character(), 
                        Phase = character(), 
                        Pulse = character(),
                        Time = numeric(),
                        Cloud_Name = character(),
                        Distance = numeric())

################################################################################
#loops through all the datasets, phases and time points to collect all the raw C1-C1.csv and C2-C2.csv files

for (Data_Set in Data_Set_List) {
  
  #based on data set being used - defines name of the data
  #based on data set being used - defines vector of Times, each corresponds to a folder in the directory
  if (Data_Set == "nascent_Xist_dynamics") {
    Data_Name <- "Dynamic"
    Time_points <- Dynamic_Time
  } else {
    Data_Name <- "Turnover"
    Time_points <- Turnover_Time
  }
  
  for (Phase in Phase_List) {
    
    #create table to store compiled data for C1-C1_NNA from all Time points - data will be progressively added
    P1_Density_Compiled <- tibble(Time = Time_points)  #Creates table with initial values in Time column of 0, 60, etc to collect Pulse_X data (CX-CX_NNA)
    
    #create table to store compiled data for C2-C2_NNA from all Time points - data will be progressively added
    P2_Density_Compiled <- tibble(Time = Time_points)
    
    #creates a loop to add all the files for each time point
    #Loop through all Time points - adds file path to locate CX_CX_NNA raw data to Density_Compiled table
    for (point in Time_points) {
      #defines file path to the sub directory where raw data is found - use Xist_turnover_on_chromatin or nascent_Xist_dynamics /Initiation or Maintenance
      path <- paste(Input_File_Path, Cell_Type, New_Line_Name, Data_Set, Phase, point, sep="/") 
      
      #tibble_of_files updates/re-writes when each Time_point is looped through to find C1-C1 files
      P1_tibble_of_files <- tibble(File_Path = list.files(path = path, pattern="_C1-C1_NNA", full.names=TRUE), #table with file paths for CX-CX_NNA files with associated Time point
                                   Time = point) %>%                                                           #uses full.names so the working directory doesn't need to change
        mutate(Pulse = "Pulse_1")
      
      #tibble_of_files updates/re-writes when each Time_point is looped through to find C2-C2 files
      P2_tibble_of_files <- tibble(File_Path = list.files(path = path, pattern="_C2-C2_NNA", full.names=TRUE), #table with file paths for CX-CX_NNA files with associated Time point
                                   Time = point) %>%                                                           #uses full.names so the working directory doesn't need to change
        mutate(Pulse = "Pulse_2")
      
      
      #Join new table for pulse 1 (P1_tibble_of_files) to Density_Compiled tibble
      P1_Density_Compiled<- drop_na(full_join(P1_Density_Compiled, P1_tibble_of_files))
      
      #Join new table for pulse 2 (P2_tibble_of_files) to Density_Compiled table
      P2_Density_Compiled<- drop_na(full_join(P2_Density_Compiled, P2_tibble_of_files))
      
      #combine data for pulse 1 and pulse 2
      Density_Compiled <- P1_Density_Compiled %>%
        bind_rows(P2_Density_Compiled)
      
    }
    
    #add a column to Density_Compiled which contains the raw CX-CX_NNA.csv file
    Density_Compiled$Raw_Data <- lapply(Density_Compiled$File_Path, read_csv)
    
    #add a column to Density_Compiled which contains the name of the cloud the Raw_Data file corresponds to
    Density_Compiled<- Density_Compiled %>%
      mutate(Cloud_Name = str_extract(File_Path, r"(\d{0,2}_Cloud\-\d{1,3})")) #extracts the cloud name from the File_path based on XX_Cloud-XX pattern
    
    #opens up raw data stored in Raw_Data column
    unnested_Density <- unnest(Density_Compiled, "Raw_Data")
    
    #removes 1st row of each Time point where volume = 0 as x,y,z = 0
    unnested_Density <- filter(unnested_Density, volume > 0) %>%
      mutate(Data_Set = Data_Name) %>%   #add the reference of the data set use
      select(-File_Path)                 #removes column as unnessassery
    
    #create table containing the useful info only
    Density_Pulse_1_and_2_Compile <- unnested_Density %>%
      transmute(Cell_Line = New_Line_Name,
                Data_Set = Data_Set,
                Phase = Phase,
                Pulse = Pulse,
                Time = Time,
                Cloud_Name = Cloud_Name,
                Distance = Median_nn_dist)
    
    New_Cell_Line <- bind_rows(New_Cell_Line, Density_Pulse_1_and_2_Compile)
  }  
}

################################################################################
#save tables

#if folder for new cell line does not exit, create folder to save density data
save_path <- paste(Output_File_Path, Cell_Type, "Density", New_Line_Name, sep="/")
dir.create(save_path)

#define file path to location of data for all other cell lines
Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Density", "All_Cell_Lines", sep="/")

#save data in both locations
File_Name <- paste(New_Line_Name, "Density_Compile.csv", sep="_")

write_csv(New_Cell_Line, paste(save_path, File_Name, sep="/"))
write_csv(New_Cell_Line, paste(Other_Data_Path, File_Name, sep="/"))

