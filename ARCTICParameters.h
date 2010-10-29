#ifndef ARCTICPARAMETERS_H
#define ARCTICPARAMETERS_H



#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

class Parameters
{

public:
	Parameters();
	~Parameters();
	void setT1(string _T1);
	void setT2(string _T2);
	void setpd(string _pd);
	void setrawImage(string _rawImage);
	void setBMSFile();
	void setBMSFile2();
	void setoutputDir();
	void setprefix(string _prefix);
	void setoutputRepCortThick(string _outputRepCortThick);
	void setsaveSkullStripping(string _saveSkullStripping);
	void setsaveAtlasRegistered(string _saveAtlasRegistered);
	void setsaveParcellationRegistered(string _saveParcellationRegistered);
	void setatlas(string _atlas);
	void setparcellation(string _parcellation);
	void setsegAtlasDir(string _segAtlasDir);
	void setsegmentationMaskingOff(bool _segmentationMaskingOff);
	/*void setorientation(string _orientation);
	void setatlasOrientation(string _atlasOrientation);
	void setrawImageOrientation(string _rawImageOrientation);
	void setsegAtlasOrientation(string _segAtlasOrientation);*/
	void setatlasType(string _atlasType);
	void setlabel(string _label);
	void setparcellationCase(string _parcellationCase);
	void setLabels(int _WMLabel,int _GMLabel,int _CSFLabel);
	void setdilateOn(bool _dilateOn);
	void setfilterIteration(int _filterIteration);
	void setfilterTimeStep(float _filterTimeStep);
	void setfilterMethod(string _filterMethod);
	void setmaxBiasDegree(int _maxBiasDegree);
	void setPrior(float _WMPrior,float _GMPrior, float _CSFPrior,float _OtherPrior);
	void setatlasWarping(bool _atlasWarpingOff);
	void setgridSize(int _gridSizeX,int _gridSizeY,int _gridSizeZ);
	//void setCortThickOptions(bool _WmOption,bool _GMMapsOption,bool _SdmOption);
	//void setRegImagesOptions(string _metric,string _initialization);
	void setRegImagesOptions(string _initialization);
	void setGMCortThick(string _GMCortThick);
	void setWMCortThick(string _WMCortThick);
	void setInterpolation(bool _InterpOff, float _Threshold);

	string getT1();
	string getT2();
	string getpd();
	string getrawImage();
	string getBMSFile();
	string getBMSFile2();
	string getoutputDir();
	string getprefix();
	string getatlas();
	string getparcellation();
	string getsegAtlasDir();
	/*string getorientation();
	string getatlasOrientation();
	string getsegAtlasOrientation();
	string getrawImageOrientation();*/
	string getatlasType();
	string getlabel();
	string getparcellationCase();
	//string getmetric();
	string getinitialization();
	string getoutputRepCortThick();
	string getsaveSkullStripping();
	string getsaveAtlasRegistered();
	string getsaveParcellationRegistered();
	int getWMLabel();
	int getGMLabel();
	int getCSFLabel();
	bool getdilateOn();
	bool getsegmentationMaskingOff();
	int getfilterIteration();
	float getfilterTimeStep();
	string getfilterMethod();
	int getmaxBiasDegree();
	float getWMPrior();
	float getGMPrior();
	float getCSFPrior();
	float getOtherPrior();
	bool getatlasWarping();
	int getgridSizeX();	
	int getgridSizeY();
	int getgridSizeZ();
	string getGMCortThick();
	string getWMCortThick();
	float getThreshold();
	bool getInterpOff();
	//bool getWmOption();
	//bool getGMMapsOption();
	//bool getSdmOption();


private:
	string m_T1;
	string m_T2;
	string m_pd;
	string m_rawImage;
	string m_outputDir;
	string m_BMSFile;
	string m_BMSFile2;
	string m_prefix;
	string m_atlas;
	string m_parcellation;
	string m_segAtlasDir;
	/*string m_orientation;
	string m_atlasOrientation;
	string m_rawImageOrientation;
	string m_segAtlasOrientation;*/
	string m_atlasType;
	string m_label;
	string m_parcellationCase;
	//string m_metric;
	string m_initialization;
	string m_saveSkullStripping;
	string m_saveAtlasRegistered;
	string m_saveParcellationRegistered;
	string m_outputRepCortThick;
	int m_GMLabel;
	int m_WMLabel;
	int m_CSFLabel;
	int m_filterIteration;
	float m_filterTimeStep;
	string m_filterMethod;
	int m_maxBiasDegree;
	float m_WMPrior;
	float m_GMPrior;
	float m_CSFPrior;
	float m_OtherPrior;
	bool m_atlasWarping;
	bool m_segmentationMaskingOff;
	int m_gridSizeX;	
	int m_gridSizeY;
	int m_gridSizeZ;
	bool m_dilateOn;
	string m_WMCortThick;
	string m_GMCortThick;
	float m_Threshold;
	bool m_InterpOff;
	//bool m_WmOption;
	//bool m_GMMapsOption;
	//bool m_SdmOption;

};

#endif
