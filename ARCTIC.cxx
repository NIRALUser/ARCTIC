#include "ARCTICCLP.h"
#include "ARCTICComputation.h"
#include "ARCTICParameters.h"




int main(int argc, char** argv)
{
	PARSE_ARGS;
	Parameters CortThickParameters;
	
	//Parameters Initialization
	CortThickParameters.setT1(T1);
	CortThickParameters.setT2(T2);
	CortThickParameters.setpd(pd);
	CortThickParameters.setrawImage(rawImage);
	CortThickParameters.setlabel(label);
	CortThickParameters.setoutputRepCortThick(outputRepCortThick);
	CortThickParameters.setoutputDir();
	CortThickParameters.setBMSFile();
	CortThickParameters.setBMSFile2();
	CortThickParameters.setprefix(prefix);
	CortThickParameters.setInterpolation(InterpOff,Threshold);
	CortThickParameters.setsegmentationMaskingOff(segmentationMaskingOff);
	/*CortThickParameters.setsaveSkullStripping(saveSkullStripping);
	CortThickParameters.setsaveAtlasRegistered(saveAtlasRegistered);
	CortThickParameters.setsaveParcellationRegistered(saveParcellationRegistered);*/
	CortThickParameters.setatlas(atlas);
	CortThickParameters.setparcellation(parcellation);
	CortThickParameters.setsegAtlasDir(segAtlasDir);
	/*CortThickParameters.setorientation(orientation);
	CortThickParameters.setrawImageOrientation(rawImageOrientation);
	CortThickParameters.setatlasOrientation(atlasOrientation);
	CortThickParameters.setsegAtlasOrientation(segAtlasOrientation);*/
	CortThickParameters.setatlasType(atlasType);
	CortThickParameters.setparcellationCase(caseParcellation);
	CortThickParameters.setLabels(WMLabel,GMLabel,CSFLabel);
	CortThickParameters.setdilateOn(dilateOn);
	CortThickParameters.setfilterIteration(filterIteration);
	CortThickParameters.setfilterTimeStep(filterTimeStep);
	CortThickParameters.setfilterMethod(filterMethod);
	CortThickParameters.setmaxBiasDegree(maxBiasDegree);
	CortThickParameters.setPrior(WMPrior,GMPrior,CSFPrior,OtherPrior);
	CortThickParameters.setatlasWarping(atlasWarpingOff);
	CortThickParameters.setgridSize(gridSizeX,gridSizeY,gridSizeZ);
	CortThickParameters.setRegImagesOptions(initialization);
	CortThickParameters.setGMCortThick(GMCortThick);
	CortThickParameters.setWMCortThick(WMCortThick);
	Computation ComputationPipeline(CortThickParameters);
	ComputationPipeline.Compute();
	return 0;
}



