setBatchMode(true);

baseFolderPath = "/home/valentin/Desktop/hepatovac/tiles/";
resultsBasePath = "/home/valentin/Desktop/hepatovac/results/masks/";
folders = getFileList(baseFolderPath);

for (f = 0; f < folders.length; f++) {
	folder = folders[f];
	selectedTilesPath = baseFolderPath + folder + "/selected_tiles/";
    maskTilesPath = baseFolderPath + folder + "/mask_tiles/";
    maskedTilesPath = baseFolderPath + folder + "/masked_tiles/";
    resultsPath = resultsBasePath + folder + "/";
    
    list = getFileList(selectedTilesPath);
    for (i = 0; i < list.length; i++) {
    	filePath = selectedTilesPath + list[i];
    	open(filePath);
    	run("Set Scale...", "distance=1 known=0.2643 unit=µm");
    	run("HSB Stack");
    	run("Convert Stack to Images");
    	selectWindow("Hue");
    	rename("0");
    	selectWindow("Saturation");
    	rename("1");
    	selectWindow("Brightness");
    	rename("2");
    	min=newArray(3);
    	max=newArray(3);
    	filter=newArray(3);
    	min[0]=0;   max[0]=255; filter[0]="pass";  
    	min[1]=0;   max[1]=15;  filter[1]="pass";  
    	min[2]=243;   max[2]=255; filter[2]="pass"; 
    	for (j=0; j<3; j++) {
        	selectWindow(""+j);
        	setThreshold(min[j], max[j]);
        	run("Convert to Mask");
        	if (filter[j]=="stop") run("Invert");
    	}
    	imageCalculator("AND create", "0", "1");
    	imageCalculator("AND create", "Result of 0", "2");
    	for (j=0; j<3; j++) {
        	selectWindow(""+j);
        	close();
    	}
    	selectWindow("Result of 0");
    	close();
    	selectWindow("Result of Result of 0");
    	a = getTitle();
    	rename(a);
    	run("Convert to Mask");
    	run("Find Edges");
    	run("Analyze Particles...", "size=5-500 circularity=0.45-1.00 show=Masks display clear include overlay");
    	resultFilePath = resultsPath + list[i] + ".csv";
    	saveAs("Results", resultFilePath);
    	close("Results");
    	maskFilePath = maskTilesPath + list[i];
    	saveAs("PNG", maskFilePath);
    	close(); 
    	selectWindow(a); 
    	close(); 
    	open(filePath);
    	rename("original_img");
    	open(maskFilePath);
    	rename("mask_img");
    	imageCalculator("AND create", "original_img", "mask_img");
    	selectWindow("Result of original_img");
    	maskedFilePath = maskedTilesPath + list[i];
    	saveAs("PNG", maskedFilePath);
    	close("original_img");
    	close("mask_img");
    	close("Result of original_img");
    	call("java.lang.System.gc");
    }
}

setBatchMode(false);
