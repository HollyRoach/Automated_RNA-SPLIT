####################################################################################################################
#SCRIPT CREATES NNA_2_to_1_Compile (NNA = nearest neighbor analysis)
####################################################################################################################
#OUTPUTS: New_Line_Name_NNA_Compile.csv and New_All_Cell_Lines_NNA_Compile.csv 
#compiles all the NN_dist_2to1 rows of data for all the C2-C1.csv files into 1 table for the new cell line (New_Line_Name_NNA_Compile.csv)
#adds this newly compiled data to a table containing the NNA measurements for all existing cell lines (New_All_Cell_Lines_NNA_Compile.csv)
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
Data_Set_List <- c(#"nascent_Xist_dynamics",
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
                        Time = numeric(),
                        Cloud_Name = character(),
                        Distance = numeric())


################################################################################
#Compile Raw data for expansion and steady-state phases

#loops through all the datasets, phases and time points to collect all the raw C2-C1.csv files
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
    
    #add a column to NNA_Compiled which contains the raw C2-C1_NNA.csv file
    NNA_Compiled$Raw_Data <- lapply(NNA_Compiled$File_Path, read_csv)
    
    #add a column to NNA_Compiled which contains the name of the cloud the Raw_Data file corresponds to
    NNA_Compiled<- NNA_Compiled %>%
      mutate(Cloud_Name = str_sub(File_Path, -44, -35))  #indexes used to extract name of cloud from the File_Path (constant in all data sets)
    
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
}

################################################################################
#compile data for EdU and Random control

#loops through all the datasets, phases and time points to collect all the EdU control files
for (Data_Set in Data_Set_List) {
  
  #based on data set being used - defines name of the data
  if (Data_Set == "nascent_Xist_dynamics") {
    Data_Name <- "Dynamic"
  } else {
    Data_Name <- "Turnover"
  }
  
  for (Phase in Phase_List) {
    
        #filepath to the location of the file EdU control for new cell line 
    EdU_path <-paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, Data_Set, "EdU_Control",  sep="/")
    
    #store name of all EdU files, for the new cell line, into a vector
    EdU_Files <- list.files(path = EdU_path, pattern = "EdU_CAL.*csv", full.names  = TRUE)   #all EdU control files have this sub string
    
    #creates a list containing the EdU tibble files
    EdU_List <- lapply(EdU_Files, read_csv)
    
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
    
    EdU_Controls <- bind_rows(EdU_Compile)
  }
}

#loops through all the datasets, phases and time points to collect all the random control files
for (Data_Set in Data_Set_List) {
  
  #based on data set being used - defines name of the data
  if (Data_Set == "nascent_Xist_dynamics") {
    Data_Name <- "Dynamic"
  } else {
    Data_Name <- "Turnover"
  }
  
  for (Phase in Phase_List) {
    
    #filepath to the location of the file EdU control for new cell line 
    Random_path <-paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, Data_Set, "Random_Control",  sep="/")
    
    #store name of all Random files, for the new cell line, into a vector
    Random_Files <- list.files(path = Random_path, pattern = "RNA_SPLIT.*csv", full.names  = TRUE)   #all random control files have this sub string 0 usually use time point = 80/100mins

    #creates a list containing the Random tibble files
    Random_List <- lapply(Random_Files, read_csv)
    
    #concatenates all the individual Random tibbles into 1 tibble
    Random_Compile <- bind_rows(Random_List) %>%
      filter(volume > 0)                    #removes 1st row of separate Random table where volume = 0 as x,y,z = 0
    
    #transform Random Compile to be in same format as original data
    Random_Compile <- Random_Compile %>%
      transmute(Cell_Line = New_Line_Name,
                Phase = "random_sample",
                Time = NA,
                Cloud_Name = NA,
                Distance =NN_dist_2to1)
    
    Random_Controls <- bind_rows(Random_Compile, Random_Compile)
  }
} 

#add internal control data for the new cell line
New_Cell_Line <- New_Cell_Line %>% 
  bind_rows(EdU_Controls,
            Random_Controls)

################################################################################
#load table for data that exists for all other cell lines

#create file path to location of data for all other cell lines
Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", "All_Cell_Lines_Merged", sep="/")

#opens either template or table containing density data for all other cell lines
if (file.exists(paste(Other_Data_Path, "New_All_Cell_Lines_NNA_Compile.csv", sep="/"))) { 
  Other_Cell_Lines <- read_csv(paste(Other_Data_Path, "New_All_Cell_Lines_NNA_Compile.csv", sep="/"))          #contains Density data for other cell lines
} else{
  Other_Cell_Lines <- read_csv(paste(Other_Data_Path, "TEMPLATE_NNA_Compile.csv", sep="/"), col_types = "ccd") #Read in template, sets col types to character and double
}

#add new cell line data to other existing cloud volume data
if (any(Other_Cell_Lines$Cell_Line == New_Line_Name)) {
  warning("This cell line already exists in Main Dataframe")
  All_Data <- Other_Cell_Lines
} else {
  All_Data <- New_Cell_Line %>%
    select(Cell_Line, Phase, Distance) %>%  #only want to add these rows of data
    bind_rows(Other_Cell_Lines)
} 


################################################################################
#save tables

#create file path to save data just for new cell line
Save_Path_1 <- paste(Output_File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, sep="/")

#create name of file which contains only the new cell line data
File_Name_1 <- paste(New_Line_Name, "NNA_Compile.csv", sep="_")

#save data for new cell line
write_csv(New_Cell_Line, paste(Save_Path_1, File_Name_1, sep="/"))


#create file path to save the new cell line data added to other cell line data
Save_Path_2 <- Other_Data_Path

#create name of file which contains all cell line data
File_Name_2 <- "New_All_Cell_Lines_NNA_Compile.csv"

#save/update data for all cell lines
write_csv(All_Data, paste(Save_Path_2, File_Name_2, sep="/"))