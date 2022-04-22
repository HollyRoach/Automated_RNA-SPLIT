################################################################################
#SCRIPT TO COMPILE CLOUD VOLUME DATA  
################################################################################

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

#define file path to where "Pulse_Chase_Analysis" is located - results from cloud_volume_calculation.imj are
File_Path <- choose.dir(default = "", caption = "Select Pulse_Chase_Analysis folder, where compiled data will be stored")


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
                        Cloud_Volume = numeric())


################################################################################
#loops through all the datasets to collect all the cloud volume tables for each time point

for (Data_Set in Data_Set_List) {
  
  #based on data set being used - defines name of the data
  if (Data_Set == "nascent_Xist_dynamics") {
    Data_Name <- "Dynamic"
  } else {
    Data_Name <- "Turnover"
  }

  #set file path to the location of the Cloud_Volume files 
  Input_Path <- paste(File_Path, Cell_Type, "Cloud_Volume", New_Line_Name, Data_Set, "Individual_Time_Points", sep="/")

  #stores name of all Cloud_Volume files, in the directory, into a vector
  Volume_Files <- list.files(path = Input_Path, pattern = "0.*csv", full.names  = TRUE)   #all Cloud_Volume files end in a 0
  
    #creates a list containing the Cloud_Volume tibble files
  Volume_List <- lapply(Volume_Files, read_csv)
  
  #concatenates all the individual Cloud_Volume tibbles into 1 tibble
  Total_Compile <- bind_rows(Volume_List) %>%
    select(-X1)
  
  #create new tibble containing the relevant information for both phase
  Phases_Complie <- Total_Compile %>%
    transmute(Cell_Line = New_Line_Name,
              Data_Set = Data_Name, 
              Phase = Phase, 
              Time = Time,
              Cloud_Name = Cloud,
              Cloud_Volume = Volume)
  
  New_Cell_Line <- bind_rows(New_Cell_Line, Phases_Complie)
  
}


################################################################################
#load table for data that exists for all other cell lines

#create file path to location of data for all other cell lines
Other_Data_Path <- paste(File_Path, Cell_Type, "Cloud_Volume", "All_Cell_Lines_Merged", sep="/")

#opens either template or table containing volume data for all other cell lines
if (file.exists(paste(Other_Data_Path, "New_All_Cell_Lines_Cloud_Volume_Compile.csv", sep="/"))) { 
  Other_Cell_Lines <- read_csv(paste(Other_Data_Path, "New_All_Cell_Lines_Cloud_Volume_Compile.csv", sep="/"))          #contains volume data for other cell lines
} else{
  Other_Cell_Lines <- read_csv(paste(Other_Data_Path, "TEMPLATE_Cloud_Volume_Compile.csv", sep="/"), col_types = "ccd") #Read in template, sets col types to character and double
}

################################################################################
#add new cell line data to other existing cloud volume data

if (any(Other_Cell_Lines$Cell_Line == New_Line_Name)) {
  warning("This cell line already exists in Main Dataframe") #prevents accidental duplication of data
  All_Data <- Other_Cell_Lines
} else {
  All_Data <- Other_Cell_Lines %>%
    bind_rows(New_Cell_Line)
} 


################################################################################
#save tables

#create file path to save data just for new cell line
Save_Path_1 <- paste(File_Path, Cell_Type, "Cloud_Volume", New_Line_Name, sep="/")

#create name of file which contains only the new cell line data
File_Name_1 <- paste(New_Line_Name, "Cloud_Volume_Compile.csv", sep="_")

#save data for new cell line
write_csv(New_Cell_Line, paste(Save_Path_1, File_Name_1, sep="/"))


#create file path to save the new cell line data added to other cell line data
Save_Path_2 <- Other_Data_Path

#create name of file which contains all cell line data
File_Name_2 <- "New_All_Cell_Lines_Cloud_Volume_Compile.csv"

#save/update data for all cell lines
write_csv(All_Data, paste(Save_Path_2, File_Name_2, sep="/"))