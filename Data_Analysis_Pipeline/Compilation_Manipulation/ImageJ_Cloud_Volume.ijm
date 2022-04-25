///////////////////////////////////////////////////////////////////////////////////
//SCRIPT TO RECORD X, Y, Z DIMENSIONS + VOLUME FOR EACH CROPPED IMAGE IN DIRECTORY
//please email Holly Rpach at hmroach@hotmail.co.uk if you have any questions
//////////////////////////////////////////////////////////////////////////////////
//script records x, y and z dimensions for each images
//based on the x, y and z dimension it can calculate the volume


////USER INPUTS REQUIRED
//set name of new cell line - should be spelled the same as folder name within the directory
//make sure not to use any spaces or symbols other than an underscore (_), dash (-), or full stop (.)
//eg "Mettl3_dTAG" or "SPEN_RRM_del" or "Ciz1_KO"
New_Line_Name = "Test";     

//define type of cells used in experiment - either "mESCs" or "NPCs"
Cell_Type = "mESCs"

//define time points used in both datasets
//make sure that the lists include all the possible time-points for your data
//it is ok if the lists contain extra time-points, as these lists are a superset of time points
//it doesn't matter if initiation/maintenance have different time points
Dynamic_Time = newArray(10, 20, 30, 40, 50, 60);
Turnover_Time = newArray(0, 60, 80, 100, 120, 140, 160, 180, 200, 220);


////OTHER INPUTS
//The below inputs should not need changing

//creates list of datasets being used for compile
Data_Set_List = newArray("nascent_Xist_dynamics", "Xist_turnover_on_chromatin");

//creates list of phases being used for compile
Phase_List = newArray("Initiation", "Maintenance");

//define file path to where foldering containing cropped cloud images for new cell line are stored
Input_File_Path = getDirectory("Select (open on Mac) folder containing all images for new cell line ("+New_Line_Name+")");

//define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
Output_File_Path = getDirectory("Select (open on Mac) Pulse_Chase_Analysis folder, where cloud volume data will be stored");


//////////////////////////////////////////////////////////////////////////////////
//main loop of analysis to caclulate cloud volume for each image per time-point, phase and dataset
setBatchMode(false);

//loop through for each dataset
for (i=0; i<Data_Set_List.length; i++) {
	
	//based on data set being used - defines name of the data
	//based on data set being used - defines vector of Times, each corresponds to a folder in the directory
  	if (Data_Set_List[i] == "nascent_Xist_dynamics") {
    	Data_Name = "Dynamic";
    	Time_points = Dynamic_Time;
  	} else {
    	Data_Name = "Turnover";
    	Time_points = Turnover_Time;
  	}
  	
  	//loop through for each Phase
  	for (j=0; j<Phase_List.length; j++) {
  		
  		//loop through for each time-point
  		for (k=0; k<Time_points.length; k++){
  			
  			//define file path to specific subfolder of images
  			file_path = Input_File_Path + File.separator + Data_Set_List[i] + File.separator + Phase_List[j] + File.separator + Time_points[k] + File.separator;
  			
  			//create list of all image files located in specific subfolder
  			list = getFileList(file_path); 
  			
  			//loop to open only the cropped cloud images within the directory
  			for (l=0; l<list.length; l++) {				
  				//create list of file paths for each image within the specific subfolder
  				Cropped_Image = file_path + list[l];
  				//only open the cropped images by specifying that the file path names ending in ALN.tif
  				if (endsWith(Cropped_Image, "ALN.tif")) {
  					open(Cropped_Image);
  				}
  			}
  			
  			//loop to calculate x, y, z and volume for each cropped image
  			for (l=0; l<nImages; l++) {
				selectImage(l+1);
				
				//get name of open image file 
				imgName = getTitle();	
				
				//define which dataset, phase and time the image comes from 
				Data_Set = Data_Set_List[i];
				Phase = Phase_List[j];
 				Time = Time_points[k];	
 				 				
				//find the x, y, z dimensions of the cropped image
				width = getWidth(); 			//returns x dimension of image in pixels
				x_dimension = width*0.0410000;	//convert width from pixels to microns
				
				height = getHeight(); 			//returns y dimension of image in pixels
				y_dimension = height*0.0410000; //convert height from pixels to microns
				
				depth = nSlices()/2; 			// returns double the number of z stacks why was it double?
				z_dimension = depth*0.125; 		//converts depth to mircons as 0.125microns between each z stack
				
				//record data into results table
				setResult("Phase", l, Phase);									//records the phase the cloud belongs to
				setResult("Time", l, Time);			 							//records the time point the Cloud belongs to
				setResult("Image_Name", l, imgName);				    		//records the name of the image 
				setResult("x_Dimension", l, x_dimension);						//records the x_dimension of the image
				setResult("y_Dimension", l, y_dimension);						//records the y_dimension of the image
				setResult("z_Dimension", l, z_dimension);						//records the z_dimension of the image
				setResult("No_Z_Stacks", l, depth);								//records number of Z stacks in the image
				setResult("Volume", l, x_dimension*y_dimension*z_dimension);	//records volume of cloud
			}
			
			//open and update the result table for each time-point
			setOption("ShowRowNumbers", false); 
			updateResults;
			
			////create files to store cloud volume data
			//create new subfolder for new cell line located in Pulse_Chase_Analysis folder
			new_cell_dir = Output_File_Path + Cell_Type + File.separator +"Cloud_Volume"+ File.separator + New_Line_Name + File.separator;
			File.makeDirectory(new_cell_dir);
			
			//within new cell line subfolder create folder for each dataset
			dataset_dir = Output_File_Path + Cell_Type + File.separator +"Cloud_Volume"+ File.separator + New_Line_Name + File.separator + Data_Set_List[i] + File.separator;
			File.makeDirectory(dataset_dir);
			
			//within new dataset subfolder create folder for each timepoint
			time_dir = Output_File_Path + Cell_Type + File.separator +"Cloud_Volume"+ File.separator + New_Line_Name + File.separator + Data_Set_List[i] + File.separator + "Individual_Time_Points" + File.separator;
			File.makeDirectory(time_dir);
			
			//if time-point is present in this dataset, save and close the results table
			if (nResults > 0) {
				selectWindow("Results");
				saveAs("Results", time_dir+"Cloud_Volume_"+Phase_List[j]+"_"+Time_points[k]+".csv");
				run("Close");
			}
			
			//closes all images
			close("*"); 			  			
  		}  		
  	}	
}

print("Finishing running macro to calculate cloud volumes for each cropped image");
print("Next Run scipt in R to compile all the Cloud Volume tables");