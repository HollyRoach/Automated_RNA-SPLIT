///////////////////////////////////////////////////////////////////////////////////
//SCRIPT TO RECORD X, Y, Z DIMENSIONS + VOLUME FOR EACH CROPPED IMAGE IN DIRECTORY
//////////////////////////////////////////////////////////////////////////////////
//need to check/update lines 31,32,35 + 63 to create correct cloud name and table name


//set directories
dir1 = getDirectory("Choose Source Directory ");		//choose directory where cropped images are loacted
dir2 = getDirectory("Choose Destination Directory ");	//choose directory where to save results table
list = getFileList(dir1);     
setBatchMode(false);

//For Loop which only opens up the cropped images within the directory
for (i=0; i<list.length; i++)  {				

	Cropped_Image = dir1 + list[i]; 			//creates list of files within the directory

	if (endsWith(Cropped_Image, "ALN.tif")) {	//speficies cropped images as file name ends in ALN.tif
   	
    	//opens all the cropped images
    	open(Cropped_Image);
	}
}

//For Loop which calculates data to be stored in the results table
for (i=0; i<nImages; i++) {
	selectImage(i+1);

	//get name of file open - change all indexes by 
	Image = getTitle();
	Cloud = substring(Image, 45, 55);	//43,53 initation-0 \44,54 initation-60 \45,55 initation-100-220 \44,54 maintenace-0 \45,55 maintenance-60 \46,56 maintenace-100-220
										//44,44 iniation-10-60 Xist_dynamics \45,55 maintenance-10-60 Xist_dynamics
	Time = substring(Image, 42, 44);	//41,42 initation-0 \41,43 initation-60 \41,44 initation-100-220 \42,43 maintenance-0 \42,44 maintenance-60 \42,45 maintenance-100-220
										//41,43 initiation-10-60 Xist_dynamics \42,44 maintenance-10-60 Xist_dynamics
	Phase = substring(Image, 30, 41);	//30,40 initiation \30,41 maintenance for Xist_turnover + Xist_Dynamics
								
	//find the x, y, z dimensions of the cropped image
    width = getWidth(); 			//returns x dimension of image in pixels
	x_dimension = width*0.0410000;	//convert width from pixels to microns

	height = getHeight(); 			//returns y dimension of image in pixels
	y_dimension = height*0.0410000; //convert height from pixels to microns

	depth = nSlices()/2; 			// returns double the number of z stacks why was it double?
	z_dimension = depth*0.125; 		//converts depth to mircons as 0.125microns between each z stack

    //record data into results table
	setResult("Phase", i, Phase);									//records the phase the cloud belongs to
	setResult("Time", i, Time);			 							//records the time point the Cloud belongs to
	setResult("Cloud", i, Cloud);				    				//records the Cloud Name 
    setResult("x_Dimension", i, x_dimension);						//records the x_dimension of the image
    setResult("y_Dimension", i, y_dimension);						//records the y_dimension of the image
    setResult("z_Dimension", i, z_dimension);						//records the z_dimension of the image
 	setResult("No_Z_Stacks", i, depth);								//records number of Z stacks in the image
	setResult("Volume", i, x_dimension*y_dimension*z_dimension);	//records volume of cloud
}

//opens and updates the result table
setOption("ShowRowNumbers", false);
updateResults;

//renames and saves the results table
Data_Set = substring(Image, 30, 44);	//42 initation-0 \43 initation-60 \44 initation-100-220 \43 maintenace-0 \44 maintenance-60 \45 maintenace-100-220 
										//43 initation-10-60 Xist_dynamics \44 maintenance-10-60 Xist_dynamics 
Table.rename("Results", "Cloud_Volume_"+Data_Set);
saveAs("Results", dir2+"Cloud_Volume_"+Data_Set+".csv");

close("*"); //closes all images

print("Successfully created table of Cloud Volumes");
print("Next Run scipt in R to compile all the Cloud Volume tables");