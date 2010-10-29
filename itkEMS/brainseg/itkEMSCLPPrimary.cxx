
#include "mu.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "EMSParameters.h"
#include "EMSParametersXMLFile.h"
#include "runEMS.h"

#include <exception>
#include <iostream>
#include <string>
#include <fstream>

#include "itkEMSCLPCLP.h"

//TODO:  Remove the using namespace command because it is dangerous when using with other c++ libraries
using namespace std;

static int EMSMode(const std::string &T1,const std::string &T2,const std::string &PD)
{
	int Mode=0;
	if(!T1.empty() && T2.empty() && PD.empty())
		Mode=1;
	if(T1.empty() && !T2.empty() && PD.empty())
		Mode=2;
	if(T1.empty() && T2.empty() && !PD.empty())
		Mode=3;
	if(!T1.empty() && !T2.empty() && PD.empty())
		Mode=4;
	if(!T1.empty() && T2.empty() && !PD.empty())
		Mode=5;
	if(T1.empty() && !T2.empty() && !PD.empty())
		Mode=6;
	if(!T1.empty() && !T2.empty() && !PD.empty())
		Mode=7;
	return Mode;
}

static void WriteXMLFile(string _XMLFile , string T1 , string T2 , string PD, string  SegAtlasDir , string Orientation , string OutputDir , string OutputFormat, int FilterIteration, float FilterTimeStep , string FilterMethod , bool AtlasWarping , int MaxBiasDegree , float WMPrior, float GMPrior , float CSFPrior , float OtherPrior , int GridSizeX , int GridSizeY , int GridSizeZ , string AtlasType)
{
	const int Mode=EMSMode(T1,T2,PD);

	ofstream XMLFile(_XMLFile.c_str());
	XMLFile<<"<?xml version=\"1.0\"?>\n"<<endl;
	XMLFile<<"<!DOCTYPE SEGMENTATION-PARAMETERS>\n"<<endl;
	XMLFile<<"<SEGMENTATION-PARAMETERS>\n"<<endl;
	XMLFile<<"<SUFFIX>EMS</SUFFIX>\n"<<endl;
	XMLFile<<"<ATLAS-DIRECTORY>"<<SegAtlasDir<<"</ATLAS-DIRECTORY>\n"<<endl;
	XMLFile<<"<ATLAS-ORIENTATION>"<<Orientation<<"</ATLAS-ORIENTATION>\n"<<endl;
	XMLFile<<"<OUTPUT-DIRECTORY>"<<OutputDir<<"/</OUTPUT-DIRECTORY>\n"<<endl;
	XMLFile<<"<OUTPUT-FORMAT>"<<OutputFormat<<"</OUTPUT-FORMAT>\n"<<endl;

	if(Mode==1)
	{	
		XMLFile<<"<IMAGE>\n"<<endl;
		XMLFile<<"	<FILE>"<<T1<<"</FILE>\n"<<endl;
		XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
		XMLFile<<"</IMAGE>\n"<<endl;
	}	
	if(Mode==2)
	{	
		XMLFile<<"<IMAGE>\n"<<endl;
		XMLFile<<"	<FILE>"<<T2<<"</FILE>\n"<<endl;
		XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
		XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==3)
	{	
		XMLFile<<"<IMAGE>\n"<<endl;
		XMLFile<<"	<FILE>"<<PD<<"</FILE>\n"<<endl;
		XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
		XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==4 && AtlasType=="T1")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==4 && AtlasType=="T2")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==5 && AtlasType=="T1")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==5 && AtlasType=="PD")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==6 && AtlasType=="T2")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==6 && AtlasType=="PD")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==7 && AtlasType=="T1")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==7 && AtlasType=="T2")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	if(Mode==7 && AtlasType=="PD")
	{
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<PD<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T1<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	        XMLFile<<"<IMAGE>\n"<<endl;
	        XMLFile<<"	<FILE>"<<T2<<"</FILE>')\n"<<endl;
	        XMLFile<<"	<ORIENTATION>"<<Orientation<<"</ORIENTATION>\n"<<endl;
	        XMLFile<<"</IMAGE>\n"<<endl;
	}
	XMLFile<<"<FILTER-ITERATIONS>"<<FilterIteration<<"</FILTER-ITERATIONS>\n"<<endl;
	XMLFile<<"<FILTER-TIME-STEP>"<<FilterTimeStep<<"</FILTER-TIME-STEP>\n"<<endl;
	XMLFile<<"<FILTER-METHOD>"<<FilterMethod<<"</FILTER-METHOD>\n"<<endl;
	XMLFile<<"<MAX-BIAS-DEGREE>"<<MaxBiasDegree<<"</MAX-BIAS-DEGREE>\n"<<endl;
	XMLFile<<"<PRIOR-1>"<<WMPrior<<"</PRIOR-1>\n"<<endl;
	XMLFile<<"<PRIOR-2>"<<GMPrior<<"</PRIOR-2>\n"<<endl;
	XMLFile<<"<PRIOR-3>"<<CSFPrior<<"</PRIOR-3>\n"<<endl;
	XMLFile<<"<PRIOR-4>"<<OtherPrior<<"</PRIOR-4>\n"<<endl;
	XMLFile<<"<DO-ATLAS-WARP>"<<!AtlasWarping<<"</DO-ATLAS-WARP>\n"<<endl;
	XMLFile<<"<ATLAS-WARP-GRID-X>"<<GridSizeX<<"</ATLAS-WARP-GRID-X>\n"<<endl;
	XMLFile<<"<ATLAS-WARP-GRID-Y>"<<GridSizeY<<"</ATLAS-WARP-GRID-Y>\n"<<endl;
	XMLFile<<"<ATLAS-WARP-GRID-Z>"<<GridSizeZ<<"</ATLAS-WARP-GRID-Z>\n"<<endl;
	XMLFile<<"</SEGMENTATION-PARAMETERS>\n"<<endl;
	
	XMLFile.close();
	
}

int itkEMSCLPPrimary(int argc, char** argv)
{
	PARSE_ARGS;

  	itk::OutputWindow::SetInstance(itk::TextOutput::New());

  	if(EMSFile.empty())
	{	
		if(OutputDir.empty())
		{
			cerr <<"Not enough parameters"<< std::endl;
			return -1;
		}
		EMSFile.append(OutputDir);
		EMSFile.append("/EMS-Param.xml");
		WriteXMLFile(EMSFile,T1,T2,PD,SegAtlasDir,Orientation,OutputDir,OutputFormat,FilterIteration,FilterTimeStep,FilterMethod,AtlasWarping,MaxBiasDegree,WMPrior,GMPrior,CSFPrior,OtherPrior,GridSizeX,GridSizeY,GridSizeZ,AtlasType);	
	}
	
	
	try
	{
	EMSParameters::Pointer emsp = readEMSParametersXML(EMSFile.c_str());
	runEMSFull(emsp, debugflag, !writeflag, templateVolume, priorsList, priorsWeightList);
	}
	catch (itk::ExceptionObject& e)
	{
	cerr << e << std::endl;
	return -1;
	}
	catch (std::exception& e)
	{
	cerr << "Exception: " << e.what() << std::endl;
	return -1;
	}
	catch (std::string& s)
	{
	cerr << "Exception: " << s << std::endl;
	return -1;
	}
	catch (...)
	{
	cerr << "Unknown exception" << std::endl;
	return -1;
	}

 	return 0;
}
