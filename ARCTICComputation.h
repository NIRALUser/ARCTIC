#ifndef ARCTICCOMPUTATION_H
#define ARCTICCOMPUTATION_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <itksys/Process.h>
#include "ARCTICParameters.h"


#include <itkImage.h>
#include <itkSize.h>
#include <itkPoint.h>
#include <itkImageFileReader.h> 
#include <itkImageFileWriter.h>

#include "bmScriptParser.h"
#include "BMString.h"

class Computation
{

public:
	Computation(Parameters _Parameters);
	~Computation();
	void setpipelineMode();
	void WriteBMSFile();
	void WriteBMSFile2();
	void ExecuteBatchmakeScript();
	void ExecuteBatchmakeScript2();
	void Compute();
	void WriteHeaderMRMLScene();
	void WriteImageMRMLScene(string ImageName,string ImagePath,bool LabelFlag,int VolumeNodeNumber,int VolumeDisplayNodeNumber);
	void WriteEndMRMLScene();
	void WriteSnapshotBeginingMRMLScene(string SceneName,int SceneSnapshotNumber,int ActiveVolumeNumber,string ActiveVolumeLabelID);
	void WriteSnapshotEndMRMLScene();
	void WriteMeshSnapshotsMRMLScene(int VolumeNodeNumber1, int VolumeNodeNumber2);
	void CenterImage(string Input, string Output);

private:
	Parameters *m_Parameters;
	ofstream CortThickPipelineBMS;
	ofstream CortThickPipelineBMS2;
	int m_pipelineMode;
	bm::ScriptParser m_Parser;
	bm::ScriptParser m_Parser2;
	string m_labelFileCentered;
	string m_refName;



};
#endif

