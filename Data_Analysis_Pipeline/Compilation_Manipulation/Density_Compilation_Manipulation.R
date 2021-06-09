####################################################################################################################
#SCRIPT COMPILES ALL DATA REQUIRED TO CALACULATE DENSITY
####################################################################################################################
#OUTPUTS: New_Line_Name_Density_Compile.csv and New_All_Cell_Lines_Density_Compile.csv 
#compiles all the Median_nn_dist rows of data for all the C1-C1.csv and C2-C2.csv files into 1 table for the new cell line (New_Line_Name_Density_Compile.csv)
#adds this newly compiled data to a table containing the density measurements for all existing cell lines (New_All_Cell_Lines_Density_Compile.csv)
#makes sure tidyverse and stringr packaged have been installed - install.packages("package_name")

#load libraries
library(tidyverse)
library(stringr)

################################################################################
#USER INPUT REQUIRED

#set name of new cell line - should be spelled the same as folder name within the directory
New_Line_Name <- "Mettl3_dTAG"           

#define type of cells used in experiment - either "mESCs" or "NPCs"
Cell_Type <- "mESCs"
  
#creates list of datasets being used for compile - if cell line has stable Xist can include Turnover and Dynamic datasets
#either "Xist_turnover_on_chromatin" or "nascent_Xist_dynamics" or both
Data_Set_List <- c("nascent_Xist_dynamics",
                   "Xist_turnover_on_chromatin")

#define time points used for this specific cell line
#doesn't matter if initiation/maintenance have different time points as long as all time points are included
Dynamic_Time <- c(10, 20, 30,40, 50, 60)
Turnover_Time <- c(0, 60, 80, 100, 120, 140, 160, 180, 200, 220)

#creates list of phases being used for compile
Phase_List <- c("Initiation", "Maintenance")
                
#define file path to where Watershed Algorithm "results" are found - needs to be located in documents due to long file path names
Input_File_Path <- choose.dir(default = "", caption = "Select Watershed_Algorithm_Results folder, where raw data is found")

#define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
Output_File_Path <- choose.dir(default = "", caption = "Select Pulse_Chase_Analysis folder, where compiled data will be stored")


################################################################################
#validate inputs

#checks inputted name of new cell line - compares name to REGEX
if (grepl("\\W", New_Line_Name)) {
  print(paste("New_Line_Name =", New_Line_Name, sep =" "))
  stop("Invalid character in name of new cell line (New_Line_Name)",
       " - must not contain any spaces or symbols")
}

#checks name of dataset
if ("nascent_Xist_dynamics" %in% Data_Set_List |"Xist_turnover_on_chromatin" %in% Data_Set_List) {
} else {
  stop("Invalid names entered in Data_Set_List - must contain either nascent_Xist_dynamics or Xist_turnover_on_chromatin or both")
}

#checks name of phase
if ("Initiation" %in% Phase_List |"Maintenance" %in% Phase_List) {
} else {
  stop("Invalid names entered in Phase_List - must contain either Initiation or Maintenance or both")
}

#checks name of cell type
stopifnot((Cell_Type == "mESCs" | Cell_Type == "NPCs"))
if (Cell_Type == "mESCs" | Cell_Type == "NPCs") {
} else {
  stop("Invalid name entered in Cell_Type - must be mESCs or NPCs")
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
  if (Data_Set == "nascent_Xist_dynamics") {
    Data_Name <- "Dynamic"
  } else {
    Data_Name <- "Turnover"
  }
  
  #based on data set being used - defines vector of Times, each corresponds to a folder in the directory
  if (Data_Set == "nascent_Xist_dynamics" ) {
    Time_points <- Dynamic_Time
  } else {
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
      mutate(Cloud_Name = str_sub(File_Path, -44, -35))  #indexes used to extract name of cloud from the File_Path (constant in all data sets)
    
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
#load table for data that exists for all other cell lines

#create file path to location of data for all other cell lines
Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Density", "All_Cell_Lines_Merged", sep="/")

#opens either template or table containing density data for all other cell lines
if (file.exists(paste(Other_Data_Path, "New_All_Cell_Lines_Density_Compile.csv", sep="/"))) { 
  Other_Cell_Lines <- read_csv(paste(Other_Data_Path, "New_All_Cell_Lines_Density_Compile.csv", sep="/"))          #contains Density data for other cell lines
} else{
  Other_Cell_Lines <- read_csv(paste(Other_Data_Path, "TEMPLATE_Density_Compile.csv", sep="/"), col_types = "ccd") #Read in template, sets col types to character and double
}

################################################################################
#add new cell line data to other existing cloud volume data

if (any(Other_Cell_Lines$Cell_Line == New_Line_Name)) {
  warning("This cell line already exists in Main Dataframe")
  All_Data <- Other_Cell_Lines
} else {
  All_Data <- Other_Cell_Lines %>%
    bind_rows(New_Cell_Line)
} 

################################################################################
#save tables

#create file path to save data just for new cell line
Save_Path_1 <- paste(Output_File_Path, Cell_Type, "Density", New_Line_Name, sep="/")

#create name of file which contains only the new cell line data
File_Name_1 <- paste(New_Line_Name, "Density_Compile.csv", sep="_")

#save data for new cell line
write_csv(New_Cell_Line, paste(Save_Path_1, File_Name_1, sep="/"))


#create file path to save the new cell line data added to other cell line data
Save_Path_2 <- Other_Data_Path

#create name of file which contains all cell line data
File_Name_2 <- "New_All_Cell_Lines_Density_Compile.csv"

#save/update data for all cell lines
write_csv(All_Data, paste(Save_Path_2, File_Name_2, sep="/"))
