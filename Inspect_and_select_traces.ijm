
var	operationMode = "Single trace onMouseOver";
var smoothTraces = 0;
var selectedTraceLineWidth = 3;
var useTraceColor = false;
var playMovie = false;
 
macro "Select Traces Tool - C314DbeDcdC314DbfC315DddC315DbdC325DadC326DdeC326C327D02D11C327C328DccC338D13C339D03C339C44aC44bC45bD21C45cD22C45cC45dC46dD14C46dC46eDaeDbcC47eC47fC48fD04D0cC48fDedC48fDfbDfeC48fD0eD8eC49fD9fC49fD1eC49fDe5DebC49fDdcC39fD0dD1cDffC39fD2fC3afD15D2bDcbC3afD0fD8fDf4C3afD7eDdaC3afD7fC3afD1bC3afD1fD2eD3cD3fDe2Df5DfcC2bfDf3C2bfD00D2cD34C2bfD1dD5fDd9C2beD0bD6eDefC2beD20DafDdfDeaC2beD6fDf2C2ceD23D4fC2ceD76De3Df1DfaC2ceC2cdD6dD9eDe6C1cdD2dDecC1cdDe9Df9C1cdD3bDf0C1ddD5eD7dC1ddD01De4C1ddD3eC1dcD3dDe8DeeC1dcD0aD5dD9aC1dcD06D1aD30D40D58Df6C1dcD44D4dDe7C1dcD57Dc3De0C1ebD16D25D9dDd3DdbDf8C1ebD4eD67C1ebD24Dd2C1ebD10C1ebD33D68C1eaC2eaDd1Df7C2eaD4cD66D69D70C2eaD05D31D8dC2eaD2aC2eaD17Da9C2f9D59Db6C2f9De1C3f9D19DfdC3f9D60C3f9D92C3f8D41D82C4f8D77C4f8D51C4f8D6aDc2Dd0C4f7Dd8C5f7D32C5f7D8aC5f7C5f6C6f6D50C6f6D47C6f6D43C6f6Da6C7f5D7aC7f5D35D78C7f5D48Db7C7f5D8bDc1C8f5D29D71D86Da8C8f5D07D79C8f4D61Db3Dc0C8f4D80C8f4D18D95Dd7C9f4D09C9f4Dd4C9f3D53D54D5aCaf3D7bDa2Db8Caf3D81D99Caf3D63D72Db2Caf3D42Da7Caf3Cbf3Da5Cbf3D62Cbf3D73Cbf3D08Cbf3D52D96Ccf3Cce3D85Cce3D5cCce3Cde3Db1Cde3D87DacCde3D83Cde3Cdd3D89Ced3D97Db0Ced3D98DaaCed3D26D88Ced3D91DcaDd5Ced3Da3Cec3Da1Cec3Cfc3D90Cfc3Da0Dd6Cfc3D49Cfb3D28Cfb3D56Cfb3Dc4Cfb3D6bCfb3D93Cfb3Cfa3D6cCfa3D27Cfa3D9bCfa3Cf92DbbCf92D45Cf92D64Cf92D7cCf82D4bDb5Cf82D38D75Cf82D3aCf82Cf72D37Cf71Cf61Db9Cf61D5bCf61D8cCf61Cf51Ce51Ce50Ce40Cd30Cd20Dc9Cc20D4aD9cCc20Cb20Cb10Db4Cb10D46Cb10Ca10Dc7Ca10Dc6Ca10D65Ca10C910C900D74C900D94C900C800D39Dc8C700D84C700D36D55Da4DabDbaDc5"{

//Initialize
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();
run("Clear Results");
close("Results");
roiManager("deselect");
roiManager("Set Color", "gray");
roiManager("Set Line Width", 1);
use_coordinates = false;

//Close previous plots, if present
window = select_image_containing_string("selected cells");
if(window != 0) close();
window = select_image_containing_string("highlighted");
if(window != 0) close();
window = select_image_containing_string("kymograph_smoothed");
if(window != 0) close();
window = select_image_containing_string("kymograph_sorted_smoothed");
if(window != 0) close();
window = select_image_containing_string("lifetime traces plot HITS");
if(window != 0) close();

plot = select_image_containing_string("lifetime traces plot");
setLocation(100, 50);
getLocationAndSize(x, y, width, height);
Overlay.remove;
plotname = getTitle();
basename = substring(plotname, 0, indexOf(plotname, " (lifetime traces plot)"));
plot_id = getImageID();
selectImage(plot_id);
print("\\Clear");
run("NKI get plot styles to log window");		//Run the groovy script to get the plot styles
logWindow = getInfo("log");
styles = split(logWindow, "\n");
print("\\Clear");

//Get nr_cells from kymograph width
kymograph = select_image_containing_string("(kymograph)");
getDimensions(kwidth, kheight, channels, slices, frames);
nr_cells = kwidth;

selectedTraces = newArray(0);
print("Retrieving values from "+nr_cells+" traces...");

//By doing this only once, time is saved when running multiple times.
//BUT: if the smooth radius is changed by the user Plot Values has to be closed manually in order to have effect on the selected traces plot!
PlotValuesTable = "Plot Values";
if(!isOpen(PlotValuesTable)) {
	selectImage(plotname);
waitForUser(getTitle());
	Plot.showValues();	//Put all traces into the Results table
	Table.rename("Results", PlotValuesTable);
//	selectWindow("Results");
//	IJ.renameResults(PlotValuesTable);
}

kymograph_sorted = select_image_containing_string("kymograph sorted");
Overlay.remove;
setLocation(x + width + 2, 50);
getLocationAndSize(x, y, width, height);
kymograph_sorted_smoothed = select_image_containing_string("kymograph_sorted_smoothed");
Overlay.remove;
getDimensions(width, height, channels, slices, frames);
kymograph = select_image_containing_string("(kymograph)");
setLocation(screenWidth-200, screenHeight-200);
rank_vector = select_image_containing_string("rank vector");
rank_array = image1DToArray(rank_vector);
setLocation(screenWidth-200, screenHeight-200);
labelmap = select_image_containing_string("labelmap_cells");
getLut(reds, greens, blues);
setLocation(screenWidth-200, screenHeight-200);
Overlay.remove;
run("Select None");
RGB_overlay = select_image_containing_string("RGB");
setLocation(x, y + height + 70);
Overlay.remove;
run("Select None");
roiManager("Show All without labels");

frameInterval = Stack.getFrameInterval();
selectImage(kymograph);
getDimensions(width, height, channels, slices, frames);
timeArray = Array.getSequence(height);
timeArray = multiplyArraywithScalar(timeArray, frameInterval);
timepoints = Table.size(PlotValuesTable);

//create smooth kymographs
if(smoothTraces > 0) {
	kymograph_smoothed = "kymograph_smoothed";
	Ext.CLIJ2_push(kymograph);
	Ext.CLIJ2_mean2DBox(kymograph, kymograph_smoothed, 0, smoothTraces);
	Ext.CLIJ2_pull(kymograph_smoothed);
	Ext.CLIJ2_release(kymograph_smoothed);
	selectImage(kymograph_smoothed);
	setLocation(screenWidth-200, screenHeight-200);

	kymograph_sorted_smoothed = "kymograph_sorted_smoothed";
	selectImage(kymograph_sorted);
	getLut(reds_kymo, greens_kymo, blues_kymo);
	Ext.CLIJ2_push(kymograph_sorted);
	Ext.CLIJ2_mean2DBox(kymograph_sorted, kymograph_sorted_smoothed, 0, smoothTraces);
	Ext.CLIJ2_pull(kymograph_sorted_smoothed);
	Ext.CLIJ2_release(kymograph_sorted_smoothed);
	selectImage(kymograph_sorted_smoothed);
	setLut(reds_kymo, greens_kymo, blues_kymo);
	setLocation(x + width + 20, y);
	
}
else kymograph_smoothed = kymograph;


call("ij.gui.ImageWindow.setNextLocation", 100, 50);
plot = plot_traces_from_kymograph(plot, kymograph_smoothed, newArray(0));	//Recreate the plot from the kymograph (faster than resetting the styles for many traces)
selectImage(plot_id);
getLocationAndSize(x, y, width, height);
run("Close");

selectImage(plot);
selectTracesToolID = toolID;

if(operationMode == "Select multiple traces") {
	setTool("freehand");
	
	run("Colors...", "foreground=white background=black selection=red");
	run("Roi Defaults...", "color=red stroke=5 group=0");
	waitForUser("Create \n* a ROI in the traces plot\n* a ROI or multi-point selection in the RGB overlay\n* a ROI in the sorted kymograph image\n to highlight cells");
	run("Roi Defaults...", "color=red stroke=1 group=0");
	getLocationAndSize(x_coord, y_coord, width, height);
	if(selectionType == -1) exit("No Selection found");
	if(selectionType == 10) run("Enlarge...", "enlarge="+cellSelectionSize);	//Enlarge (multi)point selections
	selectedTraces = newArray(nr_cells);
	
	setTool(selectTracesToolID);
	setBatchMode(true);
	
	
	//Check where the ROI is and get the selected cells
	if(getTitle() == plot) {
		selection = "plot";
		run("Add Selection...");
		run("Properties... ", "  fill=#220096ff");
		run("Add Selection...");
		run("Properties... ", "  stroke=black");
		n = 0;									//counter for selected traces
		for (i = 0; i <nr_cells; i++) {			//Loop over all traces
			yValues = Table.getColumn("Y"+i, PlotValuesTable);
			for (t = 0; t < timepoints; t++) {	//Check for all timepoints if there are data points inside the selection
				x = timeArray[t];
				y = yValues[t];
				toUnscaled(x, y);
				if (selectionContains(x, y)) {
					selectedTraces[n] = i;
					n++;
					break;
				}
			}
		}
		selectedTraces = Array.trim(selectedTraces, n);
		print(selectedTraces.length+" traces in selection:");
		Array.print(selectedTraces);
	}
	else if(getTitle() == RGB_overlay) {
		selection = "RGB";
		run("Properties... ", "  stroke=cyan width=3");
		run("Add Selection...");
		Overlay.setPosition(0);
		roiManager("add");						//Temporarily add the selection to the ROI manager
		run("Clear Results");
		selectImage(labelmap);
		roiManager("select", nr_cells);			//Remove the selection again
		roiManager("delete");
		Ext.CLIJ2_push(labelmap);
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap, labelmap);
		Ext.CLIJ2_release(labelmap);
		X = Table.getColumn("CENTROID_X", "Results");
		Y = Table.getColumn("CENTROID_Y", "Results");
		n = 0;									//counter for selected traces
			for (i = 0; i <nr_cells; i++) {
				if (selectionContains(X[i], Y[i])) {
					selectedTraces[n] = i;
					n++;
			}
		}
		selectedTraces = Array.trim(selectedTraces, n);
		print(selectedTraces.length+" traces in selection:");
		Array.print(selectedTraces);
	}
	else if(getTitle() == kymograph_sorted || getTitle() == kymograph_sorted_smoothed) {
		selection = "kymograph";
		if(selectionType()>3) exit("In the kymograph an area ROI is required");
		getSelectionBounds(x0, y, width, height);
		x1 = x0 + width;
		selectedTraces = Array.slice(rank_array, x0, x1);
		selectImage(kymograph_sorted);
		getDimensions(width, height, channels, slices, frames);
		makeRectangle(x0, 0, x1-x0, height);
		Overlay.addSelection("", 0, "#55000000");
		run("Select None");
		if(smoothTraces>0) {
			selectImage(kymograph_sorted_smoothed);
			getDimensions(width, height, channels, slices, frames);
			makeRectangle(x0, 0, x1-x0, height);
			Overlay.addSelection("", 0, "#55000000");
			run("Select None");
		}
	}
	else exit("No selection found. Please select a region in the Plot, in the RGB overlay, or in the sorted kymograph.");

	run("Colors...", "foreground=white background=black selection=cyan");
	
	plot_highlighted = plot_traces_from_kymograph(plot, kymograph_smoothed, selectedTraces);
	rename(plot + " highlighted");
	getLocationAndSize(x, y, width, height);
	
	//Create Plot of selected traces
	call("ij.gui.ImageWindow.setNextLocation", 100, y + height + 10);
	plot_selectedOnly = plot_traces(plot, PlotValuesTable, nr_cells, true, selectedTraces);
	rename(plot + " selected cells");
	
	highlight_cells_in_RGB_overlay();
	create_selected_cells_table(selection);
	close(kymograph_smoothed);
}


else if (operationMode == "Single trace onMouseOver") {
	//Runs when clicking left mouse button
	setOption("DisablePopupMenu", true);
	leftButton=16;
	x1 = -1;
	y1 = -1;
	getDimensions(width, height, channels, slices, frames);
	getCursorLoc(x, y, z, flags);
	currentCell = -1;
	oldCell = -1;
	highlightedCell = -1;
	traceNr = -1;
	plot_single_cell = "";
	hidden=false;

	//Add an empty trace that will become the highlighted trace
	selectImage(plot);
	Plot.add("line", newArray(timepoints), newArray(timepoints));

	selectImage(kymograph_sorted);
	setTool(selectTracesToolID);
	while(toolID == selectTracesToolID) {
		if(getTitle() == RGB_overlay) {
			selectedImage = getTitle();
			setBatchMode(true);
			Stack.getPosition(channel, slice, frame);
			selectImage(labelmap);
			currentCell = getPixel(x, y) - 1;	//labels start at 1
			selectImage(RGB_overlay);
			if(currentCell != -1) showStatus("cell "+currentCell + 1);
			else showStatus("background");
			style = getStyle(0, styles);
			if(currentCell != -1) { 
				highlight_ROI_in_RGB_overlay(currentCell);
				if(currentCell != oldCell && oldCell != -1) reset_ROI_to_gray(oldCell);
				if(isOpen(kymograph_sorted_smoothed)) show_line_on_sorted_kymograph(kymograph_sorted_smoothed, rank_array, currentCell, 2);
				show_line_on_sorted_kymograph(kymograph_sorted, rank_array, currentCell, 2);
//				show_horizontal_line_on_sorted_kymograph(kymograph_sorted, frame, 1);
//				show_horizontal_line_on_sorted_kymograph(kymograph_sorted_smoothed, frame, 1);

				selectImage(RGB_overlay);
				if(flags&leftButton != 0) {	//leftclick to show a single trace
					if(currentCell != -1 && highlightedCell != currentCell) {
						selectImage(plot);
						Plot.freeze(true);
						for (i = 0; i < nr_cells; i++) {
							Plot.setStyle(i, "#dddddd" +","+ style[1] +","+ style[2] +","+ style[3]);
						}
						style = getStyle(currentCell, styles);
						style[2] = 3;	//change linewidth
						Plot.replace(nr_cells, style[3], timeArray, Table.getColumn("Y"+currentCell+1, PlotValuesTable));
						Plot.setStyle(nr_cells, style[0] +","+ style[1] +","+ style[2] +","+ style[3]);
	 					Plot.freeze(false);
						selectImage(RGB_overlay);
						highlightedCell = currentCell;
						hidden = true;
					}
				}
				else if(flags&leftButton == 0 && hidden == true) {
					unhide_traces(plot, styles);
					hidden = false;
					highlightedCell = -1;
				}
			}
			else {
				roiManager("Deselect");
				roiManager("Set Color", "gray");
				roiManager("Set Line Width", 1);
				run("Select None");
				selectImage(kymograph_sorted);
				Overlay.remove;
				if(isOpen(kymograph_sorted_smoothed)) {
					selectImage(kymograph_sorted_smoothed);
					Overlay.remove;
				}
				selectImage(RGB_overlay);
				unhide_traces(plot, styles);
				hidden = false;
				highlightedCell = -1;
			}
			if(hidden == true && flags&leftButton == 0) {
				unhide_traces(plot, styles);
				hidden = false;
			}
		}
		else if(getTitle() == plot) {
			selectedImage = getTitle();
			setBatchMode(true);
			Plot.getLimits(xMin, xMax, yMin, yMax);
			Plot.getFrameBounds(plotX, plotY, plotWidth, plotHeight);
			calibrationX = (xMax - xMin) / plotWidth;
			calibrationY = (yMax - yMin) / plotHeight;
			xLoc = (x - plotX)*calibrationX + xMin;
			yLoc = -(y - plotY)*calibrationY + yMax;
			xLocPixel = round(xLoc / frameInterval);
			selectWindow(kymograph);
			makeRectangle(0, xLocPixel, nr_cells, 1);
			data = getProfile();
			closest = get_closest_value_in_array(data, yLoc);
			traceNrPrevious = traceNr;
			traceNr = nr_cells - closest[0] - 1;	//This depends on the orientation of the kymograph image
			showStatus("cell "+traceNr+1+" | t = "+xLoc+" s | tau = "+closest[1]+" ns");
			print("cell "+traceNr+1+" | t = "+xLoc+" s | tau = "+closest[1]+" ns");
			
			if(traceNr != -1) highlight_ROI_in_RGB_overlay(traceNr);
			if(traceNr != traceNrPrevious && traceNrPrevious != -1) reset_ROI_to_gray(traceNrPrevious);
			
			if(isOpen(kymograph_sorted_smoothed)) show_line_on_sorted_kymograph(kymograph_sorted_smoothed, rank_array, traceNr, 2);
			show_line_on_sorted_kymograph(kymograph_sorted, rank_array, traceNr, 2);
//			show_horizontal_line_on_sorted_kymograph(kymograph_sorted, frame, 1);
//			show_horizontal_line_on_sorted_kymograph(kymograph_sorted_smoothed, frame, 1);

			selectImage(plot);
			if(flags&leftButton != 0) {	//leftclick to highlight a single trace
				labelColor = getLabelColorFromStyles(traceNr, styles);
				Plot.setStyle(traceNr, labelColor+", none,5.0,Line");
				if(traceNrPrevious != -1 && traceNrPrevious != traceNr) {
					Plot.setStyle(traceNrPrevious, labelColor+", none,1.0,Line");
					labelColor = getLabelColorFromStyles(traceNrPrevious, styles);
				}
				if(traceNrPrevious == traceNr) {
					Plot.setStyle(traceNr, labelColor+", none,1.0,Line");
					close("Cell "+traceNr+1);
					selectWindow(labelmap);
					run("Select None");
					selectWindow(RGB_overlay);
					run("Select None");
					updateDisplay();
					traceNr = -1;
				}
				setBatchMode(true);
				if(traceNr != -1) {
					selectWindow(labelmap);
					roiManager("select", traceNr);
					close("Cell "+traceNrPrevious+1);
		//			selectWindow(kymograph);
		//			makeRectangle(traceNr, 0, 1, maxTime);
					selectWindow(RGB_overlay);
					roiManager("select", traceNr);
					getLocationAndSize(xWindow, yWindow, widthWindow, heightWindow);
					call("ij.gui.ImageWindow.setNextLocation", xWindow + widthWindow-15, yWindow);
					run("Enlarge...", "enlarge=5");	//expand ROI with 5 µm
					run("Duplicate...", "title=[Cell "+traceNr+1+"] duplicate");
					run("Select None");
					setBatchMode(false);
					wait(50);
					run("Set... ", "zoom=600");
					setBatchMode(true);
					if(playMovie) doCommand("Start Animation [\\]");

					setBatchMode("show");
					selectImage(RGB_overlay);

				}
				setBatchMode(false);
				wait(50);
			}
			
			selectWindow(RGB_overlay);
			if(playMovie) run("Stop Animation");
			if(traceNr != -1) {
				if(isOpen("Cell "+traceNr+1)) selectImage("Cell "+traceNr+1);
				if(playMovie) run("Stop Animation");
			}
		}
		else if(getTitle() == kymograph_sorted || getTitle() == kymograph_sorted_smoothed) {
			selectedImage = getTitle();
			setBatchMode(true);
			selectImage(RGB_overlay);
			Stack.getPosition(channel, slice, frame);
			selectImage(rank_vector);
			currentCell = getPixel(x, 0);	//rank_vector image starts at x=0
			if(currentCell != 0) showStatus("cell "+currentCell + 1);
			else showStatus("");
			show_line_on_sorted_kymograph(kymograph_sorted, rank_array, currentCell, 2);
			if(isOpen(kymograph_sorted_smoothed)) show_line_on_sorted_kymograph(kymograph_sorted_smoothed, rank_array, currentCell, 2);
			show_horizontal_line_on_sorted_kymograph(kymograph_sorted, frame, 1);
			if(isOpen(kymograph_sorted_smoothed)) show_horizontal_line_on_sorted_kymograph(kymograph_sorted_smoothed, frame, 1);
//			show_timeline_on_plot(plot, frame, 2);

			if(flags == 16) {	//leftclick to show a single trace
				selectImage(plot);
				Plot.freeze(true);
				for (i = 0; i < nr_cells; i++) {
					Plot.setStyle(i, "#dddddd" +","+ style[1] +","+ style[2] +","+ style[3]);
				}
				style = getStyle(currentCell, styles);
				style[2] = 3;	//change linewidth
				Plot.replace(nr_cells, style[3], timeArray, Table.getColumn("Y"+currentCell+1, PlotValuesTable));
				Plot.setStyle(nr_cells, style[0] +","+ style[1] +","+ style[2] +","+ style[3]);
				Plot.freeze(false);

				selectImage(RGB_overlay);
				hidden = true;
			}

			style = getStyle(currentCell, styles);
			if(currentCell != oldCell && currentCell != 0) {
				style = getStyle(currentCell, styles);
				selectImage(RGB_overlay);
				roiManager("select", currentCell);
				if(useTraceColor == true) roiManager("Set Color", style[0]);
				else roiManager("Set Color", "white");
				roiManager("Set Line Width", 3);
				if(oldCell != -1) {
					roiManager("select", oldCell);
					roiManager("Set Color", "gray");
					roiManager("Set Line Width", 1);
					run("Select None");
				}		
			}
			if(hidden == true && flags&leftButton == 0) {
				unhide_traces(plot, styles);
				hidden = false;
			}

			//left+ctrl to select multiple traces 
			selectImage(kymograph_sorted);
			xa = x;
			xb = x;
			ctrl = false;
			while(flags == 18) {
				getCursorLoc(xb, y, z, flags);
				ctrl = true;
				Overlay.remove;
				if(xa != xb) {
					makeRectangle(minOf(xa, xb), 0, abs(xb-xa), height);
					Overlay.addSelection("", 0, "#55000000");
					run("Select None");
					wait(10);
				}
			}
			if(ctrl == true && xa != xb) {
				print(xa, xb);
				if(xa < xb) selectedTraces = Array.slice(rank_array, maxOf(0,xa), minOf(styles.length, xb));	//prevent crash if cursor is outside the image
				else selectedTraces = Array.slice(rank_array, maxOf(0,xb), minOf(styles.length, xa));
				Array.print(selectedTraces);
				selectImage(plot);
				Plot.freeze(true);
				for (i = 0; i < nr_cells; i++) {
					Plot.setStyle(i, styles[i]+",hidden");
				}
				for (i = 0; i < selectedTraces.length; i++) {
					style = getStyle(selectedTraces[i], styles);
					style[2] = 3;	//change linewidth
					Plot.setStyle(selectedTraces[i], style[0] +","+ style[1] +","+ style[2] +","+ style[3]);
				}
				Plot.freeze(false);
				
				highlight_cells_in_RGB_overlay();
				selectImage(RGB_overlay);
				getLocationAndSize(x, y, width, height);
				selection = "kymograph";
				create_selected_cells_table(selection);
				setLocation(x + width + 10, y);

				hidden = true;
				while(!isKeyDown("space")) {
					wait(10);
				}
				ctrl = false;
				selectImage(RGB_overlay);
				roiManager("Deselect");
				roiManager("Set Color", "gray");
				roiManager("Set Line Width", 1);
				run("Select None");
				selectImage(kymograph_sorted);
			}	
		}
		oldCell = currentCell;
		selectImage(RGB_overlay);
		Stack.getPosition(channel, slice, frame);
//		roiManager("Deselect");
//		roiManager("Set Color", "gray");
//		roiManager("Set Line Width", 1);
		selectImage(kymograph_sorted);
		show_horizontal_line_on_sorted_kymograph(kymograph_sorted, frame, 1);
		if(isOpen(kymograph_sorted_smoothed)) show_horizontal_line_on_sorted_kymograph(kymograph_sorted_smoothed, frame, 1);
		setBatchMode(false);
//		selectImage(selectedImage);
		wait(25);
		getCursorLoc(x, y, z, flags);
	}
	run("Select None");
	setOption("DisablePopupMenu", false);
}
//END

function generate_rgn_file(name, hitList) {
	preamble = "<StageOverviewRegions>\n<Regions>\n<ShapeList>\n<Items>\n<Item0>\n<Name>hits</Name>\n<Identifier>"+generate_UUID()+"</Identifier>\n<Type>CompoundShape</Type>\n<Font />\n<Verticies>\n<Items />\n</Verticies>\n<DecoratorColors>\n<Items />\n</DecoratorColors>\n<ExtendedProperties>\n<Items />\n</ExtendedProperties>\n<Children>\n<Items>\n";
	postamble = "</Items>\n</Children>\n</Item0>\n</Items>\n<FillMaskMode>None</FillMaskMode>\n<VertexUnitMode>Pixels</VertexUnitMode>\n</ShapeList>\n</Regions>\n<StackList>\n";
	stacklist = "";
	output = File.directory;	//Location of the last opened file - NOT FAIL SAVE!
	if(File.exists(output + "HITS_"+name+".rgn")) File.delete(output + "HITS_"+name+".rgn");
	file = File.open(output + "HITS_"+name+".rgn");
	print(file, preamble);
	selectWindow(hitList);
	nrHits = Table.size;
	for (i = 0; i < nrHits; i++) {
		UUID = generate_UUID();
		start = "<Item"+i+">\n<Number>"+i+1+"</Number>\n<Tag>Cell_"+i+1+"</Tag>\n<Identifier>"+UUID+"</Identifier>\n<Type>Point</Type>\n<Fill>R:1,G:0,B:0,A:0</Fill>\n<Font />\n<Verticies>\n<Items>\n<Item0>\n";
		middle = "\<X\>" + d2s(Table.get("AbsPosX (m)", i),10) + "\</X\>\n\<Y\>"+ d2s(Table.get("AbsPosY (m)", i),10) + "\</Y\>\n";
		end = "</Item0>\n</Items>\n</Verticies>\n<DecoratorColors>\n<Items />\n</DecoratorColors>\n<ExtendedProperties>\n<Items />\n</ExtendedProperties>\n</Item"+i+">\n";
		print(file, start+middle+end);
		stacklist += "<Entry Identifier=\""+UUID+"\" Begin=\"0.0000000000\" End=\"0.0000000000\" SectionCount=\"0\" ReferenceX=\"0.0000000000\" ReferenceY=\"0.0000000000\" FocusStabilizerOffset=\"0.0000000000\" FocusStabilizerOffsetFixed=\"false\" StackValid=\"false\" Marked=\"false\" />\n";
	}
	print(file, postamble);
	print(file, stacklist);
	print(file, "</StackList>\n</StageOverviewRegions>\n");
	File.close(file);
	print(".RGN file saved as " + output + "HITS_"+name+".rgn");
}


function generate_UUID() {
	string = "";
	for (i = 0; i < 36; i++) {
		if(i==8 || i==13 || i==18 || i==23) string += "-";
		else string += toHex(random*16);
	}
	return string;
}


function highlight_cells_in_RGB_overlay() {
	selectImage(RGB_overlay);
	roiManager("Set Color", "gray");
	roiManager("Set Line Width", 1);
	for (n = 0; n < selectedTraces.length; n++) {
		roiManager("select", selectedTraces[n]);
		if(useTraceColor == true) {
			style = getStyle(selectedTraces[n], styles);
			roiManager("Set Color", style[0]);
		}
		else roiManager("Set Color", "white");
		roiManager("Set Line Width", 3);
	}
	roiManager("Show All without labels");
	getLocationAndSize(x, y, width, height);
}

function highlight_ROI_in_RGB_overlay(cell) {
	selectImage(RGB_overlay);
	style = getStyle(cell, styles);
	roiManager("select", cell);
	if(useTraceColor == true) roiManager("Set Color", style[0]);
	else roiManager("Set Color", "white");
	roiManager("Set Line Width", 3);
	run("Select None");
}

function reset_ROI_to_gray(cell) {
		selectImage(RGB_overlay);
		roiManager("Select", cell);
		roiManager("Set Color", "gray");
		roiManager("Set Line Width", 1);
		run("Select None");
}

function create_selected_cells_table(selection) {
	//Create a table with selected cells
	cellCoordinatesTable = "Cell_coordinates.tsv";
	if(isOpen(cellCoordinatesTable)) use_coordinates = true;
	selectedCellsTable = "Selected_cells";
	Table.create("Selected_cells");
	if(selection == "plot" || selection == "kymograph") {
		Ext.CLIJ2_push(labelmap);
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap, labelmap);
		Ext.CLIJ2_release(labelmap);
		X = Table.getColumn("CENTROID_X", "Results");
		Y = Table.getColumn("CENTROID_Y", "Results");
		selectWindow("Results");
		run("Close");
	}
	for (n = 0; n < selectedTraces.length; n++) {
		if(use_coordinates == true) {
			selectedCell = Table.get("Cell total", selectedTraces[n], cellCoordinatesTable);
			absPosX = Table.get("AbsPosX (m)", selectedTraces[n], cellCoordinatesTable);
			absPosY = Table.get("AbsPosY (m)", selectedTraces[n], cellCoordinatesTable);
			Table.set("Cell", n, selectedCell+1, selectedCellsTable);
			Table.set("AbsPosX (m)", n, absPosX, selectedCellsTable);
			Table.set("AbsPosY (m)", n, absPosY, selectedCellsTable);
		}
		else {
			Table.set("Cell", n, selectedTraces[n]+1, selectedCellsTable);
			Table.set("PosX (px)", n, X[n], selectedCellsTable);
			Table.set("PosY (px)", n, Y[n], selectedCellsTable);
		}
	}
	Table.update;
	Table.setLocationAndSize(x + width + 10, y, 250, 400);
	if(use_coordinates == true) generate_rgn_file(basename, selectedCellsTable);
}


function show_line_on_sorted_kymograph(image, rank_array, currentCell, thickness) {
	currentWindow = getTitle();
	selectImage(image);
	Overlay.remove;
	currentSortedCell = firstIndexOfArray(rank_array, currentCell);
	makeLine(currentSortedCell, 0, currentSortedCell, getHeight, thickness);
	Roi.setStrokeColor("black");
//	Roi.setStrokeWidth(2);
	run("Add Selection...");
//	Roi.setStrokeWidth(0);
	run("Select None");
	selectWindow(currentWindow);
}


function show_horizontal_line_on_sorted_kymograph(image, position, thickness) {
	currentWindow = getTitle();
	selectImage(image);
	makeLine(0, position, getWidth(), position, thickness);
	Roi.setStrokeColor("black");
//	Roi.setStrokeWidth(2);
	run("Add Selection...");
//	Roi.setStrokeWidth(0);
	run("Select None");
	selectWindow(currentWindow);
}


function show_timeline_on_plot(plot, frame, thickness) {
	currentWindow = getTitle();
	selectImage(plot);
	Plot.getLimits(xMin, xMax, yMin, yMax);
	Plot.drawLine(frame*frameInterval/xMax, yMin, frame*frameInterval/xMax, yMax);
	Plot.update();
	selectWindow(currentWindow);
}


function unhide_traces(plot, styles) {
	currentWindow = getTitle();
	selectImage(plot);
	Plot.freeze(true);
	for (i = 0; i < nr_cells; i++) {
		style = getStyle(i, styles);
		unhiddenStyle = style[0] +","+ style[1] +","+ style[2] +","+ style[3];
		Plot.setStyle(i, unhiddenStyle);
	}
	Plot.replace(nr_cells, 0, newArray(timepoints), newArray(timepoints));
	Plot.freeze(false);
	selectWindow(currentWindow);
}


//Construct plot of all traces from the kymograph image
function plot_traces_from_kymograph(template_plot, kymograph_all, selectedTraces) {
	selectImage(kymograph_all);
	nr_cells = getWidth();
	y_axis = "Lifetime (ns)";

	getDimensions(width, height, channels, slices, frames);
	timeArray = Array.getSequence(height);
	timeArray = multiplyArraywithScalar(timeArray, frameInterval);
//	plotName = substring(saveName, 0, saveName.length-4) + " (lifetime traces plot)";
	plot_new = template_plot;
	selectImage(template_plot);
	Plot.getLimits(xMin, xMax, yMin, yMax);
	call("ij.gui.ImageWindow.setNextLocation", 100, 50);
	Plot.create(plot_new, "", "");
	Plot.useTemplate(template_plot);
	Plot.setLimits(xMin, xMax, yMin, yMax);

	//setBatchMode(true);
	selectImage(kymograph_all);

	run("Rotate 90 Degrees Left");
//	run("Flip Vertically", "stack");
	getDimensions(width, height, channels, slices, frames);
	n=0;
	for (s = 0; s < slices; s++) {
	    if(slices>1) Stack.setSlice(s+1);
		for (i = 0; i < height; i++) {
			makeRectangle(0, i, width, 1);
			lifetimeData = getProfile();
			if(lifetimeData[0] != 0) {
				color = getLabelColor(i, height);
				Plot.setColor(color);
				if(occursInArray(selectedTraces, n)) Plot.setLineWidth(selectedTraceLineWidth);
				else Plot.setLineWidth(1.0);
				Plot.add("line", timeArray, lifetimeData);
				n++;
			}
		}
	}
	Plot.show();
	setBatchMode("show");
	selectImage(kymograph_all);
	run("Rotate 90 Degrees Right");
//	run("Flip Vertically", "stack");
	getDimensions(width, height, channels, slices, frames);
	run("Select None");
	setBatchMode("show");
	return plot_new;
}


//Construct plot of all traces from a Results table
function plot_traces(plotName, resultsTable, nr_cells, create_selected_plot, selectedTraces) {
	getLocationAndSize(xloc, yloc, width, height);
	//timepoints = Table.size(resultsTable);
	plotName2 = plotName+"_temp";
	selectImage(plotName);
	Plot.getLimits(xMin, xMax, yMin, yMax);
	Plot.create(plotName2, "", "");
	Plot.useTemplate(plotName);
	Plot.setLimits(xMin, xMax, yMin, yMax);
//	Plot.setFrameSize(900, 600);
	if(create_selected_plot == false) {
		for (i = 0; i < nr_cells; i++) {
			lifetimeData = Table.getColumn("Y"+i, resultsTable);
			if(lifetimeData[0] != 0) {
				style = getStyle(i, styles);
				Plot.setColor(style[0]);
				Plot.add("line", timeArray, lifetimeData);
			}
		}
	}
	else {
		for (i = 0; i < nr_cells; i++) {
			if(occursInArray(selectedTraces, i)) {
				lifetimeData = Table.getColumn("Y"+i, resultsTable);
				if(lifetimeData[0] != 0) {
					style = getStyle(i, styles);
					Plot.setColor(style[0]);
					Plot.add("line", timeArray, lifetimeData);
				}
			}
		}		
	}
	if(create_selected_plot == false) close(plotName);
	Plot.show();
	selectImage(plotName2);
	rename(plotName);
	setLocation(xloc, yloc);
	setBatchMode("show");
	return plotName;
}


function plot_single_cell_trace(plotName, resultsTable, cellToPlot) {
	selectImage(plot);	//The plot with all traces
	getLocationAndSize(xloc, yloc, width, height);
	Plot.getLimits(xMin, xMax, yMin, yMax);
	call("ij.gui.ImageWindow.setNextLocation", xloc, yloc);
	Plot.create(plotName, "", "");
	Plot.useTemplate(plot);
	Plot.setLimits(xMin, xMax, yMin, yMax);

	lifetimeData = Table.getColumn("Y"+cellToPlot, resultsTable);
	Plot.setColor(style[0]);
	Plot.setLineWidth(selectedTraceLineWidth);
	Plot.add("line", timeArray, lifetimeData);
	Plot.show();
	setBatchMode("show");
	return plotName;
}


//Returns the index and value of the closest element to the parameter value, using k-d tree
function get_closest_value_in_array(array, value) {
	p = 0;
	sorted = Array.copy(array);
	Array.sort(sorted);
//	Array.print(sorted);
	
	median = sorted[floor(sorted.length)/2];
	while(sorted.length > 2) {
		if(value < median) {
			sorted = Array.slice(sorted, 0, floor((sorted.length)/2)+1);
			median = sorted[floor((sorted.length)/2)];
		}
		else {
			sorted = Array.slice(sorted, floor((sorted.length)/2), sorted.length);
			median = sorted[floor(sorted.length)/2];
		}
	}
	if(value - sorted[0] < sorted[1] - value) closestValue = sorted[0];
	else closestValue = sorted[1];
	index = firstIndexOfArray(array, closestValue);
	//print("Closest value to "+ value +" is "+ sorted[0] +" at index "+ index +", found in "+k+" steps");
	return newArray(index, closestValue);	//return the (first?) index and the value as an array
}


//Gets the 'content style' of the dataset 'label' from the styles array, and return it as an array
function getStyle(label, styles) {
	style_array = split(styles[label], ",");
	return style_array;
}

//Select an image with a name that contains the string argument
function select_image_containing_string(name) {
	imageList = getList("image.titles");
	for(i=0 ; i<imageList.length ; i++) {
		if(matches(imageList[i],".*"+name+".*")) {
			selectImage(imageList[i]);
			return getTitle();
		}
	}
	return;
}

function getLabelColor(label, nrOfLabels) {
	color1 = IJ.pad(toHex(reds[label/nrOfLabels*255]),2);
	color2 = IJ.pad(toHex(greens[label/nrOfLabels*255]),2);
	color3 = IJ.pad(toHex(blues[label/nrOfLabels*255]),2);
	labelColor = "#"+color1+color2+color3;
	return labelColor;
}

function getLabelColorFromStyles(label, styles) {
	style_array = split(styles[label], ",");
//	print(style_array[0]);
	labelcolor = style_array[0];
	return labelColor;
}


function occursInArray(array, value) {
	for(i=0; i<array.length; i++) {
		if(array[i] == value) return true;
	}
	return false;
}

//Multiplies all elements of an array with a scalar
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]= (array[a]) * (scalar);
	}
	return multiplied_array;
}

//Returns the first index at which a value occurs in an array
function firstIndexOfArray(array, value) {
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}

//Returns the values in a 1D image (m,1,1) or (1,m,1) as an array
function image1DToArray(image) {
	getDimensions(width, height, channels, slices, frames);
	array = newArray(maxOf(width, height));
	if(height == 1) {
		for(x=0; x<width; x++) {
			array[x] = getPixel(x, 0);
		}
	}
	else if (width == 1) {
		for(y=0; y<height; y++) {
			array[y] = getPixel(0, y);
		}
	}
	else exit("Error in function 'image1DToArray': 1D image expected");
	return array;
}

}


macro "Select Traces Tool Options" {
	Dialog.create("Select Traces Options");
	Dialog.addRadioButtonGroup("", newArray("Single trace onMouseOver", "Select multiple traces"), 2, 1, "Single trace onMouseOver");
	Dialog.addMessage("\n");
	Dialog.addNumber("Smooth traces (mean filter radius)", smoothTraces);
	Dialog.addNumber("Line thickness of selected traces", selectedTraceLineWidth);
	Dialog.addCheckbox("Use trace color to highlight cells in RGB overlay? (otherwise white)", useTraceColor);
	Dialog.addCheckbox("Play single cell movie when selecting traces in the graph?", playMovie);

	Dialog.show();
	operationMode = Dialog.getRadioButton();
	smoothTraces = Dialog.getNumber();
    selectedTraceLineWidth = Dialog.getNumber();
    useTraceColor = Dialog.getCheckbox();
    playMovie = Dialog.getCheckbox();
}
