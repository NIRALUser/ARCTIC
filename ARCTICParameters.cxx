#include "ARCTICParameters.h"
#include <cstdlib>

#if defined(_WIN32) || defined(WIN32)
#define sep		'\\'
#else
#define sep		'/'
#endif

Parameters::Parameters()
{}

Parameters::~Parameters()
{}

void Parameters::setoutputDir()
{

	if(m_outputRepCortThick.empty())
	{
		size_t found;
		if(!(m_T1.empty()))
		{	
			found=m_T1.rfind(sep);
			if(found!=string::npos)
				m_outputDir=m_T1.substr(0,found+1);
			else
				m_outputDir.append(".");
		}
		else
		{
			if(!(m_label.empty()))
			{
				found=m_label.rfind(sep);
				if(found!=string::npos)
					m_outputDir=m_label.substr(0,found+1);
			}	
		}
	}
	else
		m_outputDir.append(m_outputRepCortThick);		
	//cout<<"Data Directory : "<<getoutputDir()<<endl;
}

void Parameters::setoutputRepCortThick(string _outputRepCortThick)
{
	m_outputRepCortThick.append(_outputRepCortThick);
}

void Parameters::setprefix(string _prefix)
{
	m_prefix.append(_prefix);
	
}

void Parameters::setsaveSkullStripping(string _saveSkullStripping)
{
	m_saveSkullStripping.append(_saveSkullStripping);
}

void Parameters::setsaveAtlasRegistered(string _saveAtlasRegistered)
{
	m_saveAtlasRegistered.append(_saveAtlasRegistered);
}

void Parameters::setsaveParcellationRegistered(string _saveParcellationRegistered)
{
	m_saveParcellationRegistered.append(_saveParcellationRegistered);
}

void Parameters::setT1(string _T1)
{
	m_T1.append(_T1);
}

void Parameters::setT2(string _T2)
{
	m_T2.append(_T2);
}

void Parameters::setpd(string _pd)
{
	m_pd.append(_pd);
}

void Parameters::setrawImage(string _rawImage)
{
	m_rawImage.append(_rawImage);
}

void Parameters::setBMSFile()
{
	
	m_BMSFile.append(getoutputDir());
	if(!m_prefix.empty())
		m_BMSFile=m_BMSFile+sep+m_prefix+"_ARCTIC.bms";
	else
		m_BMSFile=m_BMSFile+sep+"ARCTIC.bms";
	//cout<<"BMSFILE  "<<m_BMSFile<<endl;
}

void Parameters::setBMSFile2()
{
	
	m_BMSFile2.append(getoutputDir());

	if(!m_prefix.empty())
		m_BMSFile2=m_BMSFile2+sep+m_prefix+"_ARCTIC2.bms";
	else
		m_BMSFile2=m_BMSFile2+sep+"ARCTIC2.bms";
	//cout<<"BMSFILE2  "<<m_BMSFile2<<endl;
}

void Parameters::setatlas(string _atlas)
{
	m_atlas.append(_atlas);
}

void Parameters::setparcellation(string _parcellation)
{
	m_parcellation.append(_parcellation);
}

void Parameters::setsegAtlasDir(string _segAtlasDir)
{
	m_segAtlasDir.append(_segAtlasDir);
}

/*void Parameters::setorientation(string _orientation)
{
	m_orientation.append(_orientation);
}

void Parameters::setatlasOrientation(string _atlasOrientation)
{
	m_atlasOrientation.append(_atlasOrientation);
}

void Parameters::setrawImageOrientation(string _rawImageOrientation)
{
	m_rawImageOrientation.append(_rawImageOrientation);
}

void Parameters::setsegAtlasOrientation(string _segAtlasOrientation)
{
	m_segAtlasOrientation.append(_segAtlasOrientation);
}
*/
void Parameters::setatlasType(string _atlasType)
{
	m_atlasType.append(_atlasType);
}

void Parameters::setlabel(string _label)
{
	m_label.append(_label);
}

void Parameters::setparcellationCase(string _parcellationCase)
{
	m_parcellationCase.append(_parcellationCase);
}

void Parameters::setRegImagesOptions(string _initialization)
{

	m_initialization.append(_initialization);
}

void Parameters::setLabels(int _WMLabel,int _GMLabel,int _CSFLabel)
{
	m_WMLabel=_WMLabel;
	m_GMLabel=_GMLabel;
	m_CSFLabel=_CSFLabel;
}

void Parameters::setdilateOn(bool _dilateOn)
{
	m_dilateOn=_dilateOn;
}

void Parameters::setsegmentationMaskingOff(bool _segmentationMaskingOff)
{
	m_segmentationMaskingOff=_segmentationMaskingOff;
}

void Parameters::setfilterIteration(int _filterIteration)
{
	m_filterIteration=_filterIteration;
}

void Parameters::setfilterTimeStep(float _filterTimeStep)
{
	m_filterTimeStep=_filterTimeStep;
}

void Parameters::setfilterMethod(string _filterMethod)
{
	m_filterMethod.append(_filterMethod);
}

void Parameters::setmaxBiasDegree(int _maxBiasDegree)
{
	m_maxBiasDegree=_maxBiasDegree;
}

void Parameters::setPrior(float _WMPrior,float _GMPrior, float _CSFPrior,float _OtherPrior)
{
	m_WMPrior=_WMPrior;
	m_GMPrior=_GMPrior;
	m_CSFPrior=_CSFPrior;
	m_OtherPrior=_OtherPrior;
}

void Parameters::setatlasWarping(bool _atlasWarpingOff)
{
	m_atlasWarping=!(_atlasWarpingOff);
}

void Parameters::setgridSize(int _gridSizeX,int _gridSizeY,int _gridSizeZ)
{
	m_gridSizeX=_gridSizeX;
	m_gridSizeY=_gridSizeY;
	m_gridSizeZ=_gridSizeZ;
}

void Parameters::setGMCortThick(string _GMCortThick)
{
	m_GMCortThick.append(_GMCortThick);
}

void Parameters::setWMCortThick(string _WMCortThick)
{
	m_WMCortThick.append(_WMCortThick);
}

void Parameters::setInterpolation(bool _InterpOff, float _Threshold)
{
	m_InterpOff=_InterpOff;
	m_Threshold=_Threshold;
}

string Parameters::getT1()
{
	return m_T1;
}

string Parameters::getT2()
{
	return m_T2;
}

string Parameters::getpd()
{
	return m_pd;
}

string Parameters::getrawImage()
{
	return m_rawImage;
}

string Parameters::getprefix()
{
	return m_prefix;
}

string Parameters::getoutputDir()
{
	return m_outputDir;
}

string Parameters::getBMSFile()
{
	return m_BMSFile;
}

string Parameters::getBMSFile2()
{
	return m_BMSFile2;
}

string Parameters::getatlas()
{
	return m_atlas;
}

string Parameters::getparcellation()
{
	return m_parcellation;
}

string Parameters::getsegAtlasDir()
{
	return m_segAtlasDir;
}

/*string Parameters::getorientation()
{
	return m_orientation;
}

string Parameters::getatlasOrientation()
{
	return m_atlasOrientation;
}

string Parameters::getsegAtlasOrientation()
{
	return m_segAtlasOrientation;
}

string Parameters::getrawImageOrientation()
{
	return m_rawImageOrientation;
}*/

string Parameters::getatlasType()
{
	return m_atlasType;
}

string Parameters::getlabel()
{
	return m_label;
}

string Parameters::getparcellationCase()
{
	return m_parcellationCase;
}

string Parameters::getinitialization()
{
	return m_initialization;
}

string Parameters::getoutputRepCortThick()
{
	return m_outputRepCortThick;
}

string Parameters::getsaveSkullStripping()
{
	return m_saveSkullStripping;
}

string Parameters::getsaveAtlasRegistered()
{
	return m_saveAtlasRegistered;
}

string Parameters::getsaveParcellationRegistered()
{
	return m_saveParcellationRegistered;
}

int Parameters::getWMLabel()
{
	return m_WMLabel;
}

int Parameters::getGMLabel()
{
	return m_GMLabel;
}

int Parameters::getCSFLabel()
{
	return m_CSFLabel;
}

bool Parameters::getdilateOn()
{
	return m_dilateOn;
}

bool Parameters::getsegmentationMaskingOff()
{
	return m_segmentationMaskingOff;
}

int Parameters::getfilterIteration()
{
	return m_filterIteration;
}

float Parameters::getfilterTimeStep()
{
	return m_filterTimeStep;
}

string Parameters::getfilterMethod()
{
	return m_filterMethod;
}
int Parameters::getmaxBiasDegree()
{
	return m_maxBiasDegree;
}

float Parameters::getWMPrior()
{
	return m_WMPrior;
}

float Parameters::getGMPrior()
{
	return m_GMPrior;
}

float Parameters::getCSFPrior()
{
	return m_CSFPrior;
}

float Parameters::getOtherPrior()
{
	return m_OtherPrior;
}

bool Parameters::getatlasWarping()
{
	return m_atlasWarping;
}

int Parameters::getgridSizeX()
{
	return m_gridSizeX;
}

int Parameters::getgridSizeY()
{
	return m_gridSizeY;
}

int Parameters::getgridSizeZ()
{
	return m_gridSizeZ;
}

string Parameters::getGMCortThick()
{
	return m_GMCortThick;
}

string Parameters::getWMCortThick()
{
	return m_WMCortThick;
}
float Parameters::getThreshold()
{
	return m_Threshold;
}

bool Parameters::getInterpOff()
{
	return m_InterpOff;
}



