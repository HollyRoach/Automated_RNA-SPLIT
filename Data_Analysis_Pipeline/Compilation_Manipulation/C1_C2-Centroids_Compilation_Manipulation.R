####################################################################################################################
#SCRIPT COMPILES C1/C2-CENTROIDS AND CALCULATES TOTAL MOLECULE COUNT, TURNOVER AND TRANSCRIPTION DYNAMICS
#please contact Holly Roach at hmroach@hotmail.co.uk if you have any questions
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
                        No_Centroids = numeric())

################################################################################
#calculates number of pulse 1 and 2 centroids for every cloud 
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
    
    #create tibble to store compiled data for C1/C2-centroids from all Time points - data will be progressively added
    compiled_tibble<- tibble(Time = Time_points) #Creates tibble with initial values in Time column of 0, 60, etc
    
    
    #creates a loop to add all the files for each time point
    #Loop through all Time points - adds file path to locate C1-centroids and C2-centroids raw data to compiled table
    for (point in Time_points) {
      
      #defines file path to the sub directory where raw data is found - use Xist_turnover_on_chromatin or nascent_Xist_dynamics /Initiation or Maintenance
      path <- paste(Input_File_Path, Cell_Type, New_Line_Name, Data_Set, Phase, point, sep="/") 
      
      #tibble_of_file updates/re-writes when each Time_point is looped through
      tibble_of_files <- tibble(Pulse_1 = list.files(path = path, pattern="_C1-centroids", full.names=TRUE),  #tibble with file paths for C1/2 files with associated Time point
                                Pulse_2 =  list.files(path = path, pattern="_C2-centroids", full.names=TRUE), #uses full.names so the working directory doesn't need to change
                                Time = point) %>%                                                            
        gather('Pulse_1','Pulse_2',key='Pulse',value='File_Path') %>%                                         #combines Pulse_1/2 file columns to make Pulse/File_Path columns to ensure 1 observation per row
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
     mutate(Cloud_Name = str_extract(File_Path, r"(\d{0,2}_Cloud\-\d{1,3})")) #extracts the cloud name from the File_path
    
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
    
    #add data for this time point, phase and dataset into a final table (New_Cell_Line)
    New_Cell_Line <- bind_rows(New_Cell_Line, Pulse_1_and_2_Centroids_compile)
  }  
}

#if folder for new cell line does not exit, create folder to save centroid data
Cen_path <- paste(Output_File_Path, Cell_Type, "Centroids", New_Line_Name, sep="/")
dir.create(Cen_path)

#define file path to location of data for all other cell lines
Cen_Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Centroids", "All_Cell_Lines", sep="/")

#save data in both locations
Cen_File_Name <- paste(New_Line_Name, "Centroids_Compile.csv", sep="_")

write_csv(New_Cell_Line, paste(Cen_path, Cen_File_Name, sep="/"))
write_csv(New_Cell_Line, paste(Cen_Other_Data_Path, Cen_File_Name, sep="/"))


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

#puts pulse 1 and pulse 2 columns in the same tibble
#uses these columns to find total count
Total_Count <- Pulse_1 %>%
  full_join(Pulse_2,
            by= c("Cell_Line", "Data_Set", "Phase", "Time", "Cloud_Name"))%>% 
  transmute(Cell_Line = New_Line_Name,
            Data_Set = Data_Set,
            Phase = Phase,
            Time = Time,
            Cloud_Name = Cloud_Name,
            P1_No_Centroids = P1_No_Centroids,
            P2_No_Centroids = P2_No_Centroids,
            Total_Centroids = P1_No_Centroids + P2_No_Centroids) 

#if folder for new cell line does not exit, create folder to save molecule count data
TMC_path <- paste(Output_File_Path, Cell_Type, "Molecule_Count", New_Line_Name, sep="/")
dir.create(TMC_path)

#define file path to location of data for all other cell lines
TMC_Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Molecule_Count", "All_Cell_Lines", sep="/")

#save total molecule count data in both locations
TMC_File_Name <- paste(New_Line_Name, "Total_Molecule_Count_Compile.csv", sep="_")

write_csv(Total_Count, paste(TMC_path, TMC_File_Name, sep="/"))
write_csv(Total_Count, paste(TMC_Other_Data_Path, TMC_File_Name, sep="/"))


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
         Cloud_Name,
         No_Centroids = No_Centroids)

#if folder for new cell line does not exit, create folder to save transcription data
Trans_path <- paste(Output_File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="/")
dir.create(Trans_path)

#define file path to location of data for all other cell lines
Trans_Other_Data_Path <- paste(Output_File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines", sep="/")

#save total molecule count data in both locations
Trans_File_Name <- paste(New_Line_Name, "Transcription_Dynamics_Compile.csv", sep="_")

write_csv(Transcription, paste(Trans_path, Trans_File_Name, sep="/"))
write_csv(Transcription, paste(Trans_Other_Data_Path, Trans_File_Name, sep="/"))

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
  select(-Avg_No_Centroids) %>%
  mutate(Data_Set = Data_Name)

#if folder for new cell line does not exit, create folder to save turnover data
Turn_path <- paste(Output_File_Path, Cell_Type, "Turnover", New_Line_Name, sep="/")
dir.create(Turn_path)

#save data in both locations
Turn_File_Name <- paste(New_Line_Name, "Turnover_Compile.csv", sep="_")

write_csv(Turnover, paste(Turn_path, Turn_File_Name, sep="/"))

