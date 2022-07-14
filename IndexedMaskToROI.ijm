run("ROI Manager...");
roiManager("reset");

getStatistics(area, mean, min, max, std, histogram);
for (i=1; i<=max; i++) {
	setThreshold(i, i);
	run("Create Selection");
	roiManager("Add");
}
resetThreshold();