####################################################################################################################
#SCRIPT COMPILES ALL DATA REQUIRED TO PERFORM NNA (NNA = nearest neighbor analysis)
#please contact Holly Roach at hmroach@hotmail.co.uk if you have any questions
####################################################################################################################
#OUTPUTS: New_Line_Name_NNA_Compile.csv 
#compiles all the NN_dist_2to1 rows of data for all the C2-C1.csv files into 1 table for the new cell line (New_Line_Name_NNA_Compile.csv)
#makes sure tidyverse, stringr and tcltk packaged have been installed - install.packages("package_name")
#IMPORTANT to make sure you have created a folder specific for your new cell line in the Pulse_Chase_Analysis folder
# --> this folder should contain a nascent_Xist_dynamics and a Xist_turnover_on_chromatin folder
# --> within each dataset folder, there should be an EdU_Control and Random_Control folder which contain files corresponding to these controls

#load libraries
# library(tidyverse)
# library(stringr)
# library(tcltk)

#create function to replace choose.dir so that it is compatible with mac OS
# choose_dir <- function(caption = 'Select data directory') {
#   if (exists('utils::choose.dir')) {
#     choose.dir(caption = caption)
#   } else {
#     tk_choose.dir(caption = caption)
#   }
# }

################################################################################
# #USER INPUT REQUIRED
# 
# #set name of new cell line - should be spelled the same as folder name within the directory
# #make sure not to use any spaces or symbols other than an underscore (_), dash (-), or full stop (.)
# #eg "Mettl3_dTAG" or "SPEN_RRM_del" or "Ciz1_KO"
# New_Line_Name <- "Test"           
# 
# #define type of cells used in experiment - either "mESCs" or "NPCs"
# Cell_Type <- "mESCs"
# 
# #define time points used in both datasets
# #make sure that the lists include all the possible time-points for your data
# #it is ok if the lists contain extra time-points, as these lists are a superset of time points
# #it doesn't matter if initiation/maintenance have different time points
# Dynamic_Time <- c(10, 20, 30, 40, 50, 60)
# Turnover_Time <- c(0, 60, 80, 100, 120, 140, 160, 180, 200, 220)
# 
# 
# ################################################################################
# #OTHER INPUTS
# #the below inputs should not need changing
# 
# #creates list of datasets being used for compile
# Data_Set_List <- c("nascent_Xist_dynamics",
#                    "Xist_turnover_on_chromatin")
# 
# #creates list of phases being used for compile
# Phase_List <- c("Initiation", "Maintenance")
# 
# #define file path to where Watershed Algorithm "results" are found - needs to be located in documents due to long file path names
# Input_File_Path <- choose_dir(caption = "Select (open on Mac) Watershed_Algorithm_Results folder, where raw data is found")
# 
# #define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
# Output_File_Path <- choose_dir(caption = "Select (open on Mac) Pulse_Chase_Analysis folder, where compiled data will be stored")


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
  stop("Invalid names entered in Data_Set_List - must contain either
       nascent_Xist_dynamics or Xist_turnover_on_chromatin or both")
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
  stop("Invalid file permissions for Watershed_Algorithm_Results folder - 
       check in folders properties that you have permission to read&write")
}

#checks that user has permissions to read and write files in Pulse_Chase_Analysis folder
if (!(file.access(Output_File_Path, 2) == 0 && file.access(Output_File_Path, 4) == 0 )) {
  stop("Invalid file permissions for Pulse_Chase_Analysis folder - 
       check in folders properties that you have permission to read&write")
}

#checks that Pulse_Chase_Analysis folder contains folder for new cell line
if (!dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, sep="/"))) {
  stop(paste("Folder for the new cell line (", New_Line_Name, " - New_Line_Name) does not exist in Pulse_Chase_Analysis folder", sep=""))
}

#checks that Pulse_Chase_Analysis folder contains subfolders for either dataset for new cell line
if (!dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "nascent_Xist_dynamics", sep="/")) 
    || !dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "Xist_turnover_on_chromatin", sep="/"))
    ) {
  stop(paste("nascent_Xist_dynamics or Xist_turnover_on_chromatin subfolders do not exits for new cell line (", New_Line_Name, " - New_Line_Name) within the Pulse_Chase_Analysis folder", sep=""))
}

#checks that Pulse_Chase_Analysis folder contains the EdU controls for new cell line
if (!dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "nascent_Xist_dynamics", "EdU_Control", sep="/")) 
    || !dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "Xist_turnover_on_chromatin", "EdU_Control", sep="/"))) {
  stop(paste("EdU_Control subfolder does not exist for new cell line (", New_Line_Name, " - New_Line_Name) within the Pulse_Chase_Analysis folder", sep=""))
}

#checks that Pulse_Chase_Analysis folder contains the random controls for new cell line
if (!dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "nascent_Xist_dynamics", "Random_Control", sep="/")) 
    || !dir.exists(paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "Xist_turnover_on_chromatin", "Random_Control", sep="/"))) {
  stop(paste("Random_Control subfolder does not exist for new cell line (", New_Line_Name, " - New_Line_Name) within the Pulse_Chase_Analysis folder", sep=""))
}


################################################################################
#creates table which will store table for new cell line data

New_Cell_Line <- tibble(Cell_Line = character(),
                        Data_Set = character(), 
                        Phase = character(),
                        Time = numeric(),
                        Cloud_Name = character(),
                        Distance = numeric())


################################################################################
#Compile Raw data for expansion and steady-state phases

#loops through all the datasets, phases and time points to collect all the raw C2-C1.csv files
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
    
    #create tibble to store compiled data for C2-C1_NNA from all Time points - data will be progressively added
    NNA_Compiled <- tibble(Time = Time_points) #Creates tibble with initial values in Time column of 0, 60, etc to collect C2-C1_NNA data 
    
    #creates a loop to add all the files for each time point
    #Loop through all Time points - adds file path to locate C2_C1_NNA raw data to NNA_Compiled table
    for (point in Time_points) {
      
      #defines file path to the sub directory where raw data is found - use Xist_turnover_on_chromatin or nascent_Xist_dynamics /Initiation or Maintenance
      path <- paste(Input_File_Path, Cell_Type, New_Line_Name, Data_Set, Phase, point, sep="/")
      
      #tibble_of_file updates/re-writes when each Time_point is looped through
      tibble_of_files <- tibble(File_Path = list.files(path = path, pattern="_C2-C1_NNA", full.names=TRUE),  #tibble with file paths for C2-C1_NNA files with associated Time point
                                Time = point)                                                                 #uses full.names so the working directory doesn't need to change
      
      #Join new tibble (tibble_of_files) to NNA_Compiled
      NNA_Compiled<- drop_na(full_join(NNA_Compiled, tibble_of_files))
      
    }
    
    # Check if any data existed for this phase, if not then quit
    if (nrow(NNA_Compiled) == 0) {
      print(paste("Warning: No data found for", Data_Set, Phase, "for the new cell line (", New_Line_Name, ")", sep=" "))
      next
    }
    
    #add a column to NNA_Compiled which contains the raw C2-C1_NNA.csv file
    NNA_Compiled$Raw_Data <- lapply(NNA_Compiled$File_Path, read_csv)
    
    #add a column to NNA_Compiled which contains the name of the cloud the Raw_Data file corresponds to
    NNA_Compiled<- NNA_Compiled %>%
      mutate(Cloud_Name = str_extract(File_Path, r"(\d{0,2}_Cloud\-\d{1,3})")) #extracts the cloud name from the File_path based on XX_Cloud-XX pattern
    
    #opens up raw data stored in Raw_Data column
    unnested_NNA <- unnest(NNA_Compiled, "Raw_Data")
    
    #removes 1st row of each Time point where volume = 0 as x,y,z = 0
    unnested_NNA <- filter(unnested_NNA, volume > 0)
    
    #create file with only useful info
    NNA_2_to_1_Compile <- unnested_NNA %>%
      transmute(Cell_Line = New_Line_Name,
                Data_Set = Data_Name,
                Phase = Phase,
                Time = Time,
                Cloud_Name = Cloud_Name,
                Distance =NN_dist_2to1)
    
    New_Cell_Line <- bind_rows(New_Cell_Line, NNA_2_to_1_Compile)
  }

  
  #compile data for EdU control
  #filepath to the location of the file EdU control for new cell line 
  EdU_path <-paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, Data_Set, "EdU_Control",  sep="/")
  
  #store name of all EdU files, for the new cell line, into a vector
  EdU_Files <- list.files(path = EdU_path, pattern = "EdU_CAL.*csv", full.names  = TRUE)   #all EdU control files have this sub string
  
  #creates a list containing the EdU tibble files
  EdU_List <- lapply(EdU_Files, read_csv)
  
  # Check if any data existed for this phase, if not then quit
  if (length(EdU_List) == 0) {
    print(paste("Warning: Either EdU or random control data not found for", Data_Set, Phase, "for the new cell line (", New_Line_Name, ")", sep=" "))
    next
  }
  
  #concatenates all the individual EdU tibbles into 1 tibble
  EdU_Compile <- bind_rows(EdU_List) %>%
    filter(volume > 0)                    #removes 1st row of separate EdU table where volume = 0 as x,y,z = 0
  
  #transform EdU Compile to be in same format as original data
  EdU_Compile <- EdU_Compile %>%
    transmute(Cell_Line = New_Line_Name,
              Data_Set = Data_Name,
              Phase = "EdU",
              Time = NA,
              Cloud_Name = NA,
              Distance =NN_dist_2to1)
  
  #add EdU control data to all data
  New_Cell_Line <- bind_rows(New_Cell_Line, EdU_Compile)

  
  #compile data for random controls
  #filepath to the location of the file random control for new cell line 
  Random_path <-paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, Data_Set, "Random_Control",  sep="/")
  
  #store name of all Random files, for the new cell line, into a vector
  Random_Files <- list.files(path = Random_path, pattern = "RNA_SPLIT.*csv", full.names  = TRUE)   #all random control files have this sub string 0 usually use time point = 80/100mins
  
  #creates a list containing the Random tibble files
  Random_List <- lapply(Random_Files, read_csv)
  
  # Check if any data existed for this phase, if not then quit
  if (length(Random_List) == 0) {
    print(paste("Warning: No random control data found for", Data_Set, Phase, "for the new cell line (", New_Line_Name, ")", sep=" "))
    next
  }
  
  #concatenates all the individual Random tibbles into 1 tibble
  Random_Compile <- bind_rows(Random_List) %>%
    filter(volume > 0)                    #removes 1st row of separate Random table where volume = 0 as x,y,z = 0
  
  #transform Random Compile to be in same format as original data
  Random_Compile <- Random_Compile %>%
    transmute(Cell_Line = New_Line_Name,
              Data_Set = Data_Name,
              Phase = "random_sample",
              Time = NA,
              Cloud_Name = NA,
              Distance =NN_dist_2to1)
  
  #add random control data to all data
  New_Cell_Line <- bind_rows(New_Cell_Line, Random_Compile)
}

#separate out data into respective data_sets
turnover_NNA <- New_Cell_Line %>% 
  filter(Data_Set == "Turnover")

dynamic_NNA <- New_Cell_Line %>% 
  filter(Data_Set == "Dynamic")


################################################################################
#save data 

#define file path to location of folder specific for new cell line
save_path <- paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, sep="/")

#define file path to location of data for all other cell lines
Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", "All_Cell_Lines", sep="/")

#define file name for each dataset compile
File_Name_Dynamic <- paste(New_Line_Name, "Dynamic_Nearest_Neighbour_Compile.csv", sep="_")
File_Name_Turnover <- paste(New_Line_Name, "Turnover_Nearest_Neighbour_Compile.csv", sep="_")

#save dynamic data in both locations
write_csv(dynamic_NNA, paste(save_path, File_Name_Dynamic, sep="/"))
write_csv(dynamic_NNA, paste(Other_Data_Path, File_Name_Dynamic, sep="/"))

#save turnover data in both locations
write_csv(turnover_NNA, paste(save_path, File_Name_Turnover, sep="/"))
write_csv(turnover_NNA, paste(Other_Data_Path, File_Name_Turnover, sep="/"))

