####################################################################################################################
#SCRIPT COMPILES C1/C2-CENTROIDS AND CALCULATES TOTAL MOLECULE COUNT, TURNOVER AND TRANSCRIPTION DYNAMICS
####################################################################################################################
#OUTPUTS: New_Line_Name_Molecule_Count/Turnover/Transcription_Dynamics_Compile.csv
#compiles all data from the C1-Centroids.csv and C2-Centroids.csv files into 1 table for the new cell line (New_Line_Name_Centroids_Compile.csv)
#uses number of rows in each C1/2-Centroids.csv to calculate number of centroids for each cloud
#adds up number of pulse 1 and 2 centroids to find total molecule count 
#uses number of pulse 2 centroids to calculate transcription dynamics (from dynamic dataset only)
#normalises the number of pulse 1 centroids to calculate turnover (from turnover dataset only)
#makes sure tidyverse, stringr and tcltk packaged have been installed - install.packages("package_name")

#load libraries
library(tidyverse)
library(stringr)
library(tcltk)

#create function to replace choose.dir
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
New_Line_Name <- "Test"           

#define type of cells used in experiment - either "mESCs" or "NPCs"
Cell_Type <- "mESCs"

#creates list of datasets being used for compile - if cell line has stable Xist can include Turnover and Dynamic datasets
Data_Set_List <- c("nascent_Xist_dynamics",
                   "Xist_turnover_on_chromatin")

#define time points used in both datasets
#doesn't matter if initiation/maintenance have different time points as long as all time points are included
Dynamic_Time <- c(10, 20, 30,40, 50, 60)
Turnover_Time <- c(0, 60, 80, 100, 120, 140, 160, 180, 200, 220)

#creates list of phases being used for compile
Phase_List <- c("Initiation", "Maintenance")

#define file path to where Watershed Algorithm "results" are found - needs to be located in documents due to long file path names
Input_File_Path <- choose_dir(caption = "Select (open on Mac) Watershed_Algorithm_Results folder, where raw data is found")

#define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
Output_File_Path <- choose_dir(caption = "Select (open on Mac) Pulse_Chase_Analysis folder, where compiled data will be stored")


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

#checks correct folder has been selected for input file path
if (str_sub(Input_File_Path, -27, -1) == "Watershed_Algorithm_Results") {
} else {
  stop("Watershed_Algorithm_Results folder not selected for Input_File_Path")
}

#checks correct folder has been selected for output file path
if (str_sub(Output_File_Path, -20, -1) == "Pulse_Chase_Analysis") {
} else {
  stop("Pulse_Chase_Analysis folder not selected for Output_File_Path")
}



################################################################################
#creates table which will store table for new cell line data

New_Cell_Line <- tibble(Cell_Line = character(),
                        Data_Set = character(), 
                        Phase = character(), 
                        Pulse = character(),
                        Time = numeric(),
                        Cloud_Name = character(),
                        No_Centroids = numeric())

################################################################################
#calculates number of pulse 1 and 2 centroids for every cloud 
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
    
    #create tibble to store compiled data for C1/C2-centroids from all Time points - data will be progressively added
    compiled_tibble<- tibble(Time = Time_points) #Creates tibble with initial values in Time column of 0, 60, etc
    
    
    #creates a loop to add all the files for each time point
    #Loop through all Time points - adds file path to locate CX_CX_NNA raw data to Density_Compiled table
    for (point in Time_points) {
      
      #defines file path to the sub directory where raw data is found - use Xist_turnover_on_chromatin or nascent_Xist_dynamics /Initiation or Maintenance
      path <- paste(Input_File_Path, Cell_Type, New_Line_Name, Data_Set, Phase, point, sep="/") 
      
      #tibble_of_file updates/re-writes when each Time_point is looped through
      tibble_of_files <- tibble(Pulse_1 = list.files(path = path, pattern="_C1-centroids", full.names=TRUE),  #tibble with file paths for C1/2 files with associated Time point
                                Pulse_2 =  list.files(path = path, pattern="_C2-centroids", full.names=TRUE), #uses full.names so the working directory doesn't need to change
                                Time = point) %>%                                                            
        gather('Pulse_1','Pulse_2',key='Pulse',value='File_Path') %>%                                          #combines Pulse_1/2 file columns to make Pulse/File_Path columns to ensure 1 observation per row
        rename(Pulse = Pulse)
      
      #Join new tibble (tibble_of_files) to compiled_tibble
      compiled_tibble<- drop_na(full_join(compiled_tibble, tibble_of_files))
      
    }
    
    #add a column to compiled_tibble which contains the raw C1/2-centroid.csv file
    compiled_tibble$Raw_Data <- lapply(compiled_tibble$File_Path, read_csv)
    
    #add a column to compiled_tibble which contains the number of rows (centroids) in each Raw_Data file
    compiled_tibble$No_Centroids <- lapply(compiled_tibble$Raw_Data, nrow)
    
    #add a column to compiled_tibble which contains the name of the cloud the Raw_Data file corresponds to
    compiled_tibble<- compiled_tibble %>%
      # mutate(Cloud_Name = str_sub(File_Path, -47, -38))  #indexes used to extract name of cloud from the File_Path (constant in all data sets)
      mutate(Cloud_Name = str_extract(File_Path, r"(\d{0,2}_Cloud\-\d)"))
    
    #creates new tibble (compiled_Pulse_1_and_2_Centroids) - removes columns that are not useful (eg File_Path, Raw_Data etc)
    Pulse_1_and_2_Centroids_compile <- compiled_tibble %>%
      transmute(Cell_Line = New_Line_Name,
                Data_Set = Data_Name,
                Phase = Phase,
                Pulse = Pulse,
                Time = Time,
                Cloud_Name = Cloud_Name,
                No_Centroids = No_Centroids)
    
    #need to unnest the list-column containing the number of centroids
    Pulse_1_and_2_Centroids_compile <- unnest(Pulse_1_and_2_Centroids_compile, "No_Centroids") 
    
    New_Cell_Line <- bind_rows(New_Cell_Line, Pulse_1_and_2_Centroids_compile)
  }  
}

#create file path to location of data for all other cell lines
Cen_Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Turnover", "All_Cell_Lines_Merged", sep="/")

#opens either template or table containing density data for all other cell lines
if (file.exists(paste(Cen_Other_Data_Path, "New_All_Cell_Lines_Centroids_Compile.csv", sep="/"))) { 
  Cen_Other_Cell_Lines <- read_csv(paste(Cen_Other_Data_Path, "New_All_Cell_Lines_Centroids_Compile.csv", sep="/"))          #contains Density data for other cell lines
} else{
  Cen_Other_Cell_Lines <- read_csv(paste(Cen_Other_Data_Path, "TEMPLATE_Centroids_Compile.csv", sep="/"), col_types = "ccd") #Read in template, sets col types to character and double
}

#add new cell line data to other existing cloud volume data
if (any(Cen_Other_Cell_Lines$Cell_Line == New_Line_Name)) {
  warning("This cell line already exists in Main Dataframe")
  All_Centroids <- Cen_Other_Cell_Lines
} else {
  All_Centroids <- New_Cell_Line %>%
    select(-Cloud_Name) %>%
    bind_rows(Cen_Other_Cell_Lines)
} 


#######################################################################################################
#calculate total molecule count

#calculate total number of centroids per cloud
#separate into pulse 1 and pulse 2 data
Pulse_1 <- New_Cell_Line %>%
  filter(Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Pulse_2 <- New_Cell_Line %>%
  filter(Pulse == "Pulse_2") %>%
  rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Total_Count <- Pulse_1 %>%
  mutate(P2_No_Centroids = Pulse_2$P2_No_Centroids,
         Total_Centroids = Pulse_1$P1_No_Centroids + Pulse_2$P2_No_Centroids) %>%
  transmute(Cell_Line = New_Line_Name,                                               
            Phase = Phase,
            Total_Centroids = Total_Centroids) 

#create file path to location of data for all other cell lines
TMC_Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Molecule_Count", "All_Cell_Lines_Merged", sep="/")

#opens either template or table containing density data for all other cell lines
if (file.exists(paste(TMC_Other_Data_Path, "New_All_Cell_Lines_Molecule_Count_Compile.csv", sep="/"))) { 
  TMC_Other_Cell_Lines <- read_csv(paste(TMC_Other_Data_Path, "New_All_Cell_Lines_Molecule_Count_Compile.csv", sep="/"))          #contains Density data for other cell lines
} else{
  TMC_Other_Cell_Lines <- read_csv(paste(TMC_Other_Data_Path, "TEMPLATE_Molecule_Count_Compile.csv", sep="/"), col_types = "ccd") #Read in template, sets col types to character and double
}

#add new cell line data to other existing total molecule count data
if (any(TMC_Other_Cell_Lines$Cell_Line == New_Line_Name)) {
  warning("This cell line already exists in Main Dataframe")
  All_Total_Count <- TMC_Other_Cell_Lines
} else {
  All_Total_Count <- TMC_Other_Cell_Lines %>%
    bind_rows(Total_Count)
} 


#######################################################################################################
#calculate transcription dynamics

Transcription <- New_Cell_Line %>%
  filter(Data_Set == "Dynamic") %>%              #only want to include nascent_Xist_dynamics data
  filter(Pulse == "Pulse_2") %>%                 #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  mutate(Time = case_when(Time == 10 ~ 20,       #new data needs to grouped by times - 20, 40, 60
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                 #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

#create file path to location of data for all other cell lines
Trans_Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines_Merged", sep="/")

#opens either template or table containing density data for all other cell lines
if (file.exists(paste(Trans_Other_Data_Path, "New_All_Cell_Lines_Transcription_Dynamics_Compile.csv", sep="/"))) { 
  Trans_Other_Cell_Lines <- read_csv(paste(Trans_Other_Data_Path, "New_All_Cell_Lines_Transcription_Dynamics_Compile.csv", sep="/"))          #contains Density data for other cell lines
} else{
  Trans_Other_Cell_Lines <- read_csv(paste(Trans_Other_Data_Path, "TEMPLATE_Transcription_Dynamics_Compile.csv", sep="/"), col_types = "ccd") #Read in template, sets col types to character and double
}

#add new cell line data to other existing transcription dynamics data
if (any(Trans_Other_Cell_Lines$Cell_Line == New_Line_Name)) {
  warning("This cell line already exists in Main Dataframe")
  All_Transcription <- Trans_Other_Cell_Lines
} else {
  All_Transcription <- Trans_Other_Cell_Lines %>%
    bind_rows(Transcription)
} 

#######################################################################################################
#calculate turnover
#find the maximum Avg_No_Centroids for C1 data points per time point
#use this maximum to normalise No_Centroids for C1-centroids data

#Find Maximum Avg_No_Centroids for C1-Centroids
Turnover <- New_Cell_Line %>%                        
  filter(Data_Set == "Turnover") %>%                                        #only want to include Xist_turnover_on_chromatin data                                                              #create tibble which shows just Avg_No_Centroids for C1 data
  filter(Pulse == "Pulse_1") %>%                                            #removes C2 Centroids Data
  group_by(Time) %>%                                                        #groups rows in compiled_tibble based on Time
  mutate(Avg_No_Centroids = mean(unlist(No_Centroids)))                     #finds average number of centroids for each Time (eg 0, 60, 100 etc)
  
Max_C1_Centroids <- max(Turnover$Avg_No_Centroids, na.rm = TRUE) #creates variable equal to maximum average number of C1 centroids

#normalise C1 No_Centroids using Max_C1_Centroids
Turnover <- Turnover %>%                                         
  mutate(Normalised_No_Centroids = (No_Centroids/Max_C1_Centroids)*100) %>% #adds a column which normalises No_Centroids based on maximum for C1 data
  select(-Avg_No_Centroids, -Cloud_Name) %>%
  mutate(Data_Set = Data_Name)

#as turnover data is exported to graph pad prism, it does not get added to other cell line data

#######################################################################################################
#save all data

#save centroid data
Save_Path_1 <- paste(Output_File_Path, Cell_Type, "Turnover", New_Line_Name, sep="/")

#create name of file which contains only the new cell line data
File_Name_1 <- paste(New_Line_Name, "Centroids_Compile.csv", sep="_")

#save data for new cell line
write_csv(New_Cell_Line, paste(Save_Path_1, File_Name_1, sep="/"))
write_csv(All_Centroids, paste(Cen_Other_Data_Path, "New_All_Cell_Lines_Centroids_Compile.csv", sep="/"))


#save total molecule count data
Save_Path_2 <- paste(Output_File_Path, Cell_Type, "Molecule_Count", New_Line_Name, sep="/")

File_Name_2 <- paste(New_Line_Name, "Molecule_Count_Compile.csv", sep="_")

write_csv(Total_Count, paste(Save_Path_2, File_Name_2, sep="/"))
write_csv(All_Total_Count, paste(TMC_Other_Data_Path, "New_All_Cell_Lines_Molecule_Count_Compile.csv", sep="/"))

#save transcription dynamics data
Save_Path_3 <- paste(Output_File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="/")

File_Name_3 <- paste(New_Line_Name, "Transcription_Dynamics_Compile.csv", sep="_")

write_csv(Transcription, paste(Save_Path_3, File_Name_3, sep="/"))
write_csv(All_Transcription, paste(Trans_Other_Data_Path, "New_All_Cell_Lines_Transcription_Dynamics_Compile.csv", sep="/"))

#save turnover data
File_Name_4 <- paste(New_Line_Name, "Turnover_Compile.csv", sep="_")

write_csv(Turnover, paste(Save_Path_1, File_Name_4, sep="/"))