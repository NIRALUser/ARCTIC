 #include "ARCTICComputation.h"

typedef unsigned char InputPixelType;

const unsigned int Dim = 3;

typedef itk::Image< InputPixelType, Dim > ImageType;
typedef itk::ImageFileReader <ImageType> VolumeReaderType;
typedef itk::ImageFileWriter< ImageType > VolumeWriterType;

typedef itk::Point<int,Dim> PointType;
typedef itk::Size<Dim> SizeType;



using namespace std;

#if defined(_WIN32) || defined(WIN32)
#define sep	'\\'
#else
#define sep	'/'
#endif


Computation::Computation(Parameters _Parameters)
{
	m_pipelineMode=0;
	m_Parameters=&(_Parameters);
}

Computation::~Computation()
{

}

void Computation::setpipelineMode()
{
	if((m_Parameters->getlabel()).empty()) 
	{
		if((m_Parameters->getparcellation()).empty() && (m_Parameters->getparcellationCase()).empty())
			m_pipelineMode=1;
		else
		{
			if(!((m_Parameters->getparcellation()).empty()))
				m_pipelineMode=5;
			if(!((m_Parameters->getparcellationCase()).empty()))
				m_pipelineMode=3;
		}
	}
	else
	{
		if((m_Parameters->getparcellation()).empty() && (m_Parameters->getparcellationCase()).empty())
			m_pipelineMode=2;
		else
		{
			if(!((m_Parameters->getparcellation()).empty()))
				m_pipelineMode=6;
			if(!((m_Parameters->getparcellationCase()).empty()))
				m_pipelineMode=4;
		}
	}	
}

void Computation::WriteBMSFile()
{
	CortThickPipelineBMS.open((m_Parameters->getBMSFile()).c_str());

	//script
	CortThickPipelineBMS<<"#APPLICATIONS"<<endl;
	CortThickPipelineBMS<<"SetApp(segmentation @itkEMSCLP)"<<endl;
	CortThickPipelineBMS<<"SetApp(skullstripping @SegPostProcessCLP)"<<endl;
	CortThickPipelineBMS<<"SetApp(skullstrippingT1 @SegPostProcessCLP)"<<endl;
	CortThickPipelineBMS<<"SetApp(skullstrippingT2 @SegPostProcessCLP)"<<endl;
	CortThickPipelineBMS<<"SetApp(skullstrippingpd @SegPostProcessCLP)"<<endl;
	CortThickPipelineBMS<<"SetApp(registration @RegisterImages)"<<endl;
	CortThickPipelineBMS<<"SetApp(parcelRegistration @ResampleVolume2)"<<endl;
	CortThickPipelineBMS<<"SetApp(corticalThickness @CortThickCLP)"<<endl;
	CortThickPipelineBMS<<"SetApp(imageMath @ImageMath)"<<endl;
	CortThickPipelineBMS<<"SetApp(imageStat @ImageStat)"<<endl;
	CortThickPipelineBMS<<"SetApp(orientImage @OrientImage)"<<endl;
	CortThickPipelineBMS<<""<<endl;

	CortThickPipelineBMS<<"#DATA"<<endl;
	CortThickPipelineBMS<<"Set(T1File '"<<m_Parameters->getT1()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(T2File '"<<m_Parameters->getT2()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(pdFile '"<<m_Parameters->getpd()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(rawImage '"<<m_Parameters->getrawImage()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(pdFile '"<<m_Parameters->getpd()<<"')"<<endl;
	CortThickPipelineBMS<<"GetFileName(T1OriginalName ${T1File} NAME_WITHOUT_EXTENSION)"<<endl;
	CortThickPipelineBMS<<"GetFileName(T2OriginalName ${T2File} NAME_WITHOUT_EXTENSION)"<<endl;
	CortThickPipelineBMS<<"GetFileName(pdOriginalName ${pdFile} NAME_WITHOUT_EXTENSION)"<<endl;
	CortThickPipelineBMS<<"GetFileName(rawImageOriginalName ${rawImage} NAME_WITHOUT_EXTENSION)"<<endl;
	CortThickPipelineBMS<<"Set(atlasFile '"<<m_Parameters->getatlas()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(parcellationFile '"<<m_Parameters->getparcellation()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(atlasType '"<<m_Parameters->getatlasType()<<"')"<<endl;
	CortThickPipelineBMS<<"GetFileName(atlasDir ${atlasFile} PATH)"<<endl;
	CortThickPipelineBMS<<"Set(atlasSegDir '"<<m_Parameters->getsegAtlasDir()<<"')"<<endl;
	/*CortThickPipelineBMS<<"Set(orientation '"<<m_Parameters->getorientation()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(atlasOrientation '"<<m_Parameters->getatlasOrientation()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(segAtlasOrientation '"<<m_Parameters->getsegAtlasOrientation()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(rawImageOrientation '"<<m_Parameters->getrawImageOrientation()<<"')"<<endl;*/
	CortThickPipelineBMS<<"Set(outputFormat 'nrrd')"<<endl;
	CortThickPipelineBMS<<"Set(labelWhite '1')"<<endl;
	CortThickPipelineBMS<<"Set(labelGrey '2')"<<endl;
	CortThickPipelineBMS<<""<<endl;

	CortThickPipelineBMS<<"#MODES"<<endl;
	CortThickPipelineBMS<<"Set(Mode '0')"<<endl;
	CortThickPipelineBMS<<"#Mode1 : T1 only"<<endl;
	CortThickPipelineBMS<<"#Mode2 : T1 and T2"<<endl;
	CortThickPipelineBMS<<"#Mode3 : T1 and pd"<<endl;
	CortThickPipelineBMS<<"#Mode4 : T1, T2 and pd"<<endl;
	CortThickPipelineBMS<<"If(${T1File} != '' && ${T2File} == '') "<<endl;
	CortThickPipelineBMS<<"   If(${pdFile} == '') "<<endl;
	CortThickPipelineBMS<<"      Set(Mode '1')"<<endl;
	CortThickPipelineBMS<<"   EndIf(${pdFile})"<<endl;
	CortThickPipelineBMS<<"   If(${pdFile} != '') "<<endl;
	CortThickPipelineBMS<<"      Set(Mode '3')"<<endl;
	CortThickPipelineBMS<<"   EndIf(${pdFile})"<<endl;
	CortThickPipelineBMS<<"EndIf(${T1File})"<<endl;
	CortThickPipelineBMS<<"If(${T1File} != '' && ${T2File} != '') "<<endl;
	CortThickPipelineBMS<<"   If(${pdFile} == '') "<<endl;
	CortThickPipelineBMS<<"      Set(Mode '2')"<<endl;
	CortThickPipelineBMS<<"   EndIf(${pdFile})"<<endl;
	CortThickPipelineBMS<<"   If(${pdFile} != '') "<<endl;
	CortThickPipelineBMS<<"      Set(Mode '4')"<<endl;
	CortThickPipelineBMS<<"   EndIf(${pdFile})"<<endl;
	CortThickPipelineBMS<<"EndIf(${T1File})"<<endl;
	CortThickPipelineBMS<<""<<endl;

	if(!((m_Parameters->getprefix()).empty()))
	{
		CortThickPipelineBMS<<"Set(Prefix '"<<m_Parameters->getprefix()<<"')"<<endl;
		CortThickPipelineBMS<<"Set(RefName '"<<m_Parameters->getprefix()<<"')"<<endl;
	}
	else
	{
		if(!((m_Parameters->getrawImage()).empty()))
		{
			CortThickPipelineBMS<<"Set(RefName ${rawImageOriginalName})"<<endl;

		}
		else
		{
			CortThickPipelineBMS<<"If(${Mode} == '1' || ${Mode} == '3')"<<endl;
			CortThickPipelineBMS<<"	   Set(RefName ${T1OriginalName})"<<endl;
			CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		
			CortThickPipelineBMS<<"If(${Mode} == '2' || ${Mode} == '4')"<<endl;
			CortThickPipelineBMS<<"    If(${atlasType} == 'T1')"<<endl;
			CortThickPipelineBMS<<"		Set(RefName ${T1OriginalName})"<<endl;
			CortThickPipelineBMS<<"    EndIf(${atlasType})"<<endl;
			CortThickPipelineBMS<<"    If(${atlasType} == 'T2')"<<endl;
			CortThickPipelineBMS<<"		Set(RefName ${T2OriginalName})"<<endl;
			CortThickPipelineBMS<<"    EndIf(${atlasType})"<<endl;
			CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		}
	}
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"#DIRECTORIES"<<endl;
	CortThickPipelineBMS<<"Set(dataDir '"<<m_Parameters->getoutputDir()<<"')"<<endl;
	CortThickPipelineBMS<<"Set(pipelineDir ${dataDir}/ARCTIC)"<<endl;
	CortThickPipelineBMS<<"echo('Pipeline results directory : ' ${pipelineDir})"<<endl;
	CortThickPipelineBMS<<"Set(EMSDir ${pipelineDir}/TissueSegmentation)"<<endl;
	CortThickPipelineBMS<<"Set(RegistrationDir ${pipelineDir}/Registration)"<<endl;
	CortThickPipelineBMS<<"Set(CortthickDir ${pipelineDir}/CorticalThickness)"<<endl;
	CortThickPipelineBMS<<"Set(StatDir ${pipelineDir}/Stat)"<<endl;
	CortThickPipelineBMS<<"Set(MeshDir ${pipelineDir}/Mesh)"<<endl;
	CortThickPipelineBMS<<"MakeDirectory(${pipelineDir})"<<endl;
	CortThickPipelineBMS<<"MakeDirectory(${EMSDir})"<<endl;
	if(m_pipelineMode!=2)
	{
		CortThickPipelineBMS<<"MakeDirectory(${RegistrationDir})"<<endl;
	}
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"#OPTIONS"<<endl;
	CortThickPipelineBMS<<"Set(dilateOn '"<<m_Parameters->getdilateOn()<<"')"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"Set(T1correctedFile ${EMSDir}/${T1OriginalName}_corrected_EMS.${outputFormat}) "<<endl;
	CortThickPipelineBMS<<"Set(T2correctedFile ${EMSDir}/${T2OriginalName}_corrected_EMS.${outputFormat}) "<<endl;
	CortThickPipelineBMS<<"Set(pdcorrectedFile ${EMSDir}/${pdOriginalName}_corrected_EMS.${outputFormat}) "<<endl;

	if(m_pipelineMode==1 || m_pipelineMode==3 || m_pipelineMode==5)
	{

		/*CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##############################Set SegAtlas files Orientation######################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		


		if(strcmp( (m_Parameters->getsegAtlasOrientation()).c_str() , (m_Parameters->getorientation()).c_str() ) != 0)
		{
			CortThickPipelineBMS<<"Set(EMSFiles template.gipl white.gipl gray.gipl csf.gipl rest.gipl)"<<endl;
			
			CortThickPipelineBMS<<"ForEach(File ${EMSFiles})"<<endl;
			CortThickPipelineBMS<<"   SetAppOption(orientImage.inputVolume1 ${atlasSegDir}/${File})"<<endl;
			CortThickPipelineBMS<<"   SetAppOption(orientImage.outputVolume ${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"   SetAppOption(orientImage.orientation.orientation ${orientation})"<<endl;
			CortThickPipelineBMS<<"   Run(output ${orientImage})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
	
			CortThickPipelineBMS<<"Set(atlasSegDir ${EMSDir})"<<endl;
			CortThickPipelineBMS<<"Set(segAtlasOrientation ${orientation})"<<endl;
		}*/


		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##############################itkEMS : Tissue segmentation########################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"echo('Tissue segmentation ...')"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testT1corrected ${EMSDir} ${T1OriginalName}_corrected_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testT2corrected ${EMSDir} ${T2OriginalName}_corrected_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testpdcorrected ${EMSDir} ${pdOriginalName}_corrected_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(runEMSOn '0')"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '1')"<<endl;
		CortThickPipelineBMS<<"    If(${testT1corrected} == '')"<<endl;
		CortThickPipelineBMS<<"        Set(runEMSOn '1')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${testT1corrected})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '2')"<<endl;
		CortThickPipelineBMS<<"    If(${testT1corrected} == '' || ${testT2corrected} == '')"<<endl;
		CortThickPipelineBMS<<"        Set(runEMSOn '1')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${testT1corrected})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '3')"<<endl;
		CortThickPipelineBMS<<"    If(${testT1corrected} == '' || ${testpdcorrected} == '')"<<endl;
		CortThickPipelineBMS<<"        Set(runEMSOn '1')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${testT1corrected})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '4')"<<endl;
		CortThickPipelineBMS<<"    If(${testT1corrected} == '' || ${testT2corrected} == '')"<<endl;
		CortThickPipelineBMS<<"    	If(${testpdcorrected} == '')"<<endl;
		CortThickPipelineBMS<<"       	  Set(runEMSOn '1')"<<endl;
		CortThickPipelineBMS<<"  	EndIf(${testpdcorrected})"<<endl;
		CortThickPipelineBMS<<"    EndIf(${testT1corrected})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testLabelsEMS ${EMSDir} *_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"If(${testLabelsEMS} == '')"<<endl;
		CortThickPipelineBMS<<"    Set(runEMSOn '1')"<<endl;
		CortThickPipelineBMS<<"EndIf(${testLabelsEMS})"<<endl;
		CortThickPipelineBMS<<"If(${runEMSOn} == '1')"<<endl;
		CortThickPipelineBMS<<"    Set(EMSfile ${EMSDir}/EMS-param.xml)"<<endl;
		CortThickPipelineBMS<<"    #CREATE EMS FILE"<<endl;
		CortThickPipelineBMS<<"    WriteFile(${EMSfile} '<?xml version=\"1.0\"?>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<!DOCTYPE SEGMENTATION-PARAMETERS>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<SEGMENTATION-PARAMETERS>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<SUFFIX>EMS</SUFFIX>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<ATLAS-DIRECTORY>'${atlasSegDir}'</ATLAS-DIRECTORY>\\n')"<<endl;
		/*if(!((m_Parameters->getsegAtlasOrientation()).empty()))
			CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<ATLAS-ORIENTATION>'${segAtlasOrientation}'</ATLAS-ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<OUTPUT-DIRECTORY>'${EMSDir}/'</OUTPUT-DIRECTORY>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<OUTPUT-FORMAT>Nrrd</OUTPUT-FORMAT>\\n')"<<endl;
		CortThickPipelineBMS<<"    If(${Mode} == '1')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T1File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"    If(${Mode} == '2' && ${atlasType} == 'T1')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T1File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T2File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"    If(${Mode} == '2' && ${atlasType} == T2)"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T2File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T1File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"    If(${Mode} == '3')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T1File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${pdFile}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"    If(${Mode} == '4' && ${atlasType} == T1)"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T1File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T2File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${pdFile}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"    If(${Mode} == '4' && ${atlasType} == 'T2')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T2File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${T1File}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>')\\n"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '<IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <FILE>'${pdFile}'</FILE>\\n')"<<endl;
		/*if(!((m_Parameters->getorientation()).empty()))
			CortThickPipelineBMS<<"        AppendFile(${EMSfile} '  <ORIENTATION>'${orientation}'</ORIENTATION>\\n')"<<endl;*/
		CortThickPipelineBMS<<"        AppendFile(${EMSfile} '</IMAGE>\\n')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<FILTER-ITERATIONS>"<<m_Parameters->getfilterIteration()<<"</FILTER-ITERATIONS>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<FILTER-TIME-STEP>"<<m_Parameters->getfilterTimeStep()<<"</FILTER-TIME-STEP>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<FILTER-METHOD>"<<m_Parameters->getfilterMethod()<<"</FILTER-METHOD>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<MAX-BIAS-DEGREE>"<<m_Parameters->getmaxBiasDegree()<<"</MAX-BIAS-DEGREE>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<PRIOR-1>"<<m_Parameters->getWMPrior()<<"</PRIOR-1>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<PRIOR-2>"<<m_Parameters->getGMPrior()<<"</PRIOR-2>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<PRIOR-3>"<<m_Parameters->getCSFPrior()<<"</PRIOR-3>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<PRIOR-4>"<<m_Parameters->getOtherPrior()<<"</PRIOR-4>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<DO-ATLAS-WARP>"<<m_Parameters->getatlasWarping()<<"</DO-ATLAS-WARP>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<ATLAS-WARP-GRID-X>"<<m_Parameters->getgridSizeX()<<"</ATLAS-WARP-GRID-X>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<ATLAS-WARP-GRID-Y>"<<m_Parameters->getgridSizeY()<<"</ATLAS-WARP-GRID-Y>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '<ATLAS-WARP-GRID-Z>"<<m_Parameters->getgridSizeZ()<<"</ATLAS-WARP-GRID-Z>\\n')"<<endl;
		CortThickPipelineBMS<<"    AppendFile(${EMSfile} '</SEGMENTATION-PARAMETERS>\\n')"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"    SetAppOption(segmentation.EMSFile.EMSFile ${EMSfile})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"    Run(outputseg ${segmentation})"<<endl;
		CortThickPipelineBMS<<"EndIf(${runEMS})"<<endl;
		CortThickPipelineBMS<<""<<endl;

		/*if(strcmp( (m_Parameters->getsegAtlasOrientation()).c_str() , (m_Parameters->getorientation()).c_str() ) != 0)
		{
		
			CortThickPipelineBMS<<"ForEach(File ${EMSFiles})"<<endl;
			CortThickPipelineBMS<<"   DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
	
			CortThickPipelineBMS<<"Set(atlasSegDir ${EMSDir})"<<endl;
			CortThickPipelineBMS<<"Set(segAtlasOrientation ${orientation})"<<endl;
		}*/
	}

	if(m_pipelineMode==1 || m_pipelineMode==3)
	{
		CortThickPipelineBMS<<"If(${atlasType} == 'T2')"<<endl;
		CortThickPipelineBMS<<"	  If(${Mode} == '2' || ${Mode} = '4')"<<endl;
		CortThickPipelineBMS<<"       Set(labelFile ${EMSDir}/${T2OriginalName}_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"   EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T1')"<<endl;
		CortThickPipelineBMS<<"    Set(labelFile ${EMSDir}/${T1OriginalName}_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '1' || ${Mode} = '3')"<<endl;
		CortThickPipelineBMS<<"    Set(labelFile ${EMSDir}/${T1OriginalName}_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
	}

	if(m_pipelineMode==4|| m_pipelineMode==2)
		CortThickPipelineBMS<<"Set(labelFile '"<<m_Parameters->getlabel()<<"')"<<endl;

	if(m_pipelineMode==5)
	{
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"#############################SegPostProcess : Skullstripping######################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"echo('Registration ...')"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"Set(T1strippedFile ${RegistrationDir}/${T1OriginalName}_corrected_EMS_stripped.${outputFormat})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T2')"<<endl;
		CortThickPipelineBMS<<"	  If(${Mode} == '2' || ${Mode} = '4')"<<endl;
		CortThickPipelineBMS<<"     Set(labelFile ${EMSDir}/${T2OriginalName}_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"   EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T1')"<<endl;
		CortThickPipelineBMS<<"    Set(labelFile ${EMSDir}/${T1OriginalName}_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '1' || ${Mode} = '3')"<<endl;
		CortThickPipelineBMS<<"    Set(labelFile ${EMSDir}/${T1OriginalName}_labels_EMS.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"GetFileName(T1StrippedName ${T1strippedFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"#T1"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testT1stripped ${RegistrationDir} ${T1StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(runSkullStripT1On '0')"<<endl;
		CortThickPipelineBMS<<"If(${testT1stripped} == '')"<<endl;
		CortThickPipelineBMS<<"    Set(runSkullStripT1On '1')"<<endl;
		CortThickPipelineBMS<<"EndIf(${testT1stripped})"<<endl;
		CortThickPipelineBMS<<"If(${runSkullStripT1On} == '1')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstrippingT1.input ${labelFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstrippingT1.output ${T1strippedFile})"<<endl;	
		CortThickPipelineBMS<<"    SetAppOption(skullstrippingT1.skullstripping.skullstripping ${T1correctedFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstrippingT1.GM.GM '"<<m_Parameters->getGMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstrippingT1.WM.WM '"<<m_Parameters->getWMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstrippingT1.CSF.CSF '"<<m_Parameters->getCSFLabel()<<"')"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"    If(${dilateOn} == '1')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT1.dilate '1') "<<endl;
		CortThickPipelineBMS<<"    EndIf(${dilateOn})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"    Run(outputskullstripT1 ${skullstrippingT1})"<<endl;
		CortThickPipelineBMS<<"EndIf(${runSkullStripT1On})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"#T2"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '2' || ${Mode} == '4')"<<endl;
		CortThickPipelineBMS<<"    Set(T2strippedFile ${RegistrationDir}/${T2OriginalName}_corrected_EMS_stripped.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    GetFileName(T2StrippedName ${T2strippedFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"    ListFileInDir(testT2stripped ${RegistrationDir} ${T2StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    Set(runSkullStripT2On '0')"<<endl;
		CortThickPipelineBMS<<"    If(${testT2stripped} == '')"<<endl;
		CortThickPipelineBMS<<"        Set(runSkullStripT2On '1')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${testT2stripped})"<<endl;
		CortThickPipelineBMS<<"    If(${runSkullStripT2On} == '1')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT2.input ${labelFile})"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT2.output ${T2strippedFile})"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT2.skullstripping.skullstripping ${T2correctedFile})"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT2.GM.GM '"<<m_Parameters->getGMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT2.WM.WM '"<<m_Parameters->getWMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingT2.CSF.CSF '"<<m_Parameters->getCSFLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"        If(${dilateOn} == '1')"<<endl;
		CortThickPipelineBMS<<"            SetAppOption(skullstrippingT2.dilate '1') "<<endl;
		CortThickPipelineBMS<<"        EndIf(${dilateOn})"<<endl;
		CortThickPipelineBMS<<"        Run(outputskullstripT2 ${skullstrippingT2})"<<endl;
		CortThickPipelineBMS<<"    EndIf(${runSkullStripT2On})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"#pd"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '3' || ${Mode} == '4')"<<endl;
		CortThickPipelineBMS<<"    Set(pdstrippedFile ${RegistrationDir}/${pdOriginalName}_corrected_EMS_stripped.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    GetFileName(pdStrippedName ${pdstrippedFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"    ListFileInDir(testpdstripped ${RegistrationDir} ${pdStrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    Set(runSkullStrippdOn '0')"<<endl;
		CortThickPipelineBMS<<"    If(${testpdstripped} == '')"<<endl;
		CortThickPipelineBMS<<"        Set(runSkullStrippdOn '1')"<<endl;
		CortThickPipelineBMS<<"    EndIf(${testpdstripped})"<<endl;
		CortThickPipelineBMS<<"    If(${runSkullStrippdOn} == '1')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingpd.input ${labelFile})"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingpd.output ${pdstrippedFile})"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingpd.skullstripping.skullstripping ${pdcorrectedFile})"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingpd.GM.GM '"<<m_Parameters->getGMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingpd.WM.WM '"<<m_Parameters->getWMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"        SetAppOption(skullstrippingpd.CSF.CSF '"<<m_Parameters->getGMLabel()<<"')"<<endl;
		CortThickPipelineBMS<<"        If(${dilateOn} == '1')"<<endl;
		CortThickPipelineBMS<<"            SetAppOption(skullstrippingpd.dilate '1') "<<endl;
		CortThickPipelineBMS<<"        EndIf(${dilateOn})"<<endl;
		CortThickPipelineBMS<<"        Run(outputskullstrippd ${skullstrippingpd})"<<endl;
		CortThickPipelineBMS<<"    EndIf(${runSkullStrippdOn})     "<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<""<<endl;

		/*CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##########################Set Atlas/Parcellation Orientation######################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		


		if(strcmp( (m_Parameters->getatlasOrientation()).c_str() , (m_Parameters->getorientation()).c_str() ) != 0)
		{
		//set atlas orientation
		CortThickPipelineBMS<<"SetAppOption(orientImage.inputVolume1 ${atlasFile})"<<endl;
		CortThickPipelineBMS<<"GetFileName(atlasName ${atlasFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"Set(orientedAtlas ${RegistrationDir}/${atlasName}_${orientation}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.outputVolume ${orientedAtlas})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.orientation.orientation ${orientation})"<<endl;
		CortThickPipelineBMS<<"Run(output ${orientImage})"<<endl;
		
		//set parcellation orientation
		CortThickPipelineBMS<<"SetAppOption(orientImage.inputVolume1 ${parcellationFile})"<<endl;
		CortThickPipelineBMS<<"GetFileName(parcellationName ${parcellationFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"Set(orientedParcellation ${RegistrationDir}/${parcellationName}_${orientation}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.outputVolume ${orientedParcellation})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.orientation.orientation ${orientation})"<<endl;
		CortThickPipelineBMS<<"Run(output ${orientImage})"<<endl;

		CortThickPipelineBMS<<"Set(atlasFile ${orientedAtlas})"<<endl;
		CortThickPipelineBMS<<"Set(parcellationFile ${orientedParcellation})"<<endl;
		}*/

		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##########################RegisterImages : Atlas Registration#####################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T2')"<<endl;
		CortThickPipelineBMS<<"	  If(${Mode} == '2' || ${Mode} = '4')"<<endl;
		CortThickPipelineBMS<<"     Set(atlasRegisteredFile ${RegistrationDir}/AtlasRegistered_${T2StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"     Set(transformFile ${RegistrationDir}/AtlasTransform_${T2StrippedName}.txt)"<<endl;
		CortThickPipelineBMS<<"     SetAppOption(registration.fixedImage ${T2strippedFile})"<<endl;
		CortThickPipelineBMS<<"   EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T1')    "<<endl;
		CortThickPipelineBMS<<"    Set(atlasRegisteredFile ${RegistrationDir}/AtlasRegistered_${T1StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    Set(transformFile ${RegistrationDir}/AtlasTransform_${T1StrippedName}.txt)"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.fixedImage ${T1strippedFile})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '1' || ${Mode} = '3')    "<<endl;
		CortThickPipelineBMS<<"    Set(atlasRegisteredFile ${RegistrationDir}/AtlasRegistered_${T1StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    Set(transformFile ${RegistrationDir}/AtlasTransform_${T1StrippedName}.txt)"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.fixedImage ${T1strippedFile})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"GetFileName(atlasName ${atlasRegisteredFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"GetFileName(transformName ${transformFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testatlasreg ${RegistrationDir} ${atlasName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testtransform ${RegistrationDir} ${transformName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(runAtlasRegOn '0')"<<endl;
		CortThickPipelineBMS<<"If(${testatlasreg} == '' && ${testtransform} == '')"<<endl;
		CortThickPipelineBMS<<"    Set(runAtlasRegOn '1')"<<endl;
		CortThickPipelineBMS<<"EndIf(${testatlasreg})"<<endl;
		CortThickPipelineBMS<<"If(${runAtlasRegOn} == '1')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.movingImage ${atlasFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.saveTransform.saveTransform ${transformFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.resampledImage.resampledImage ${atlasRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.registration.registration 'PipelineBSpline')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.initialization.initialization "<<m_Parameters->getinitialization()<<")"<<endl;
		CortThickPipelineBMS<<"    Run(outputreg ${registration})"<<endl;
		CortThickPipelineBMS<<"EndIf(${runAtlasRegOn})"<<endl;
		CortThickPipelineBMS<<""<<endl;

		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"########################ResampleVolume2 : Parcellation registration###############################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T2')"<<endl;
		CortThickPipelineBMS<<"	  If(${Mode} == '2' || ${Mode} = '4')"<<endl;
		CortThickPipelineBMS<<"    Set(parcellationRegisteredFile ${RegistrationDir}/ParcellationRegistered_${T2StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Reference.Reference ${T2strippedFile})"<<endl;
		CortThickPipelineBMS<<"   EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${atlasType} == 'T1')"<<endl;
		CortThickPipelineBMS<<"    Set(parcellationRegisteredFile ${RegistrationDir}/ParcellationRegistered_${T1StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Reference.Reference ${T1strippedFile})"<<endl;
		CortThickPipelineBMS<<"EndIf(${atlasType})"<<endl;
		CortThickPipelineBMS<<"If(${Mode} == '1' || ${Mode} == '3')"<<endl;
		CortThickPipelineBMS<<"    Set(parcellationRegisteredFile ${RegistrationDir}/ParcellationRegistered_${T1StrippedName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Reference.Reference ${T1strippedFile})"<<endl;
		CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"GetFileName(parcelregName ${parcellationRegisteredFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testparcelreg ${RegistrationDir} ${parcelregName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(runParcelRegOn '0')"<<endl;
		CortThickPipelineBMS<<"If(${testparcelreg} == '')"<<endl;
		CortThickPipelineBMS<<"    Set(runParcelRegOn '1')"<<endl;
		CortThickPipelineBMS<<"EndIf(${testparcelreg})"<<endl;
		CortThickPipelineBMS<<"If(${runParcelRegOn} == '1')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Infile ${parcellationFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Outfile ${parcellationRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.transformationFile.transformationFile ${transformFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.interpolation.interpolation 'nn')"<<endl;
		CortThickPipelineBMS<<"    Run(outputparcelreg ${parcelRegistration})"<<endl;
		CortThickPipelineBMS<<"EndIf(${runParcelRegOn})"<<endl;
		CortThickPipelineBMS<<""<<endl;
	}
	if(m_pipelineMode==6)
	{
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"#############################SegPostProcess : Skullstripping######################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"echo('Registration ...')"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"Set(labelFile "<<m_Parameters->getlabel()<<")"<<endl;
		CortThickPipelineBMS<<"GetFileName(rawImageName ${rawImage} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"Set(strippedFile ${RegistrationDir}/${rawImageName}_corrected_EMS_stripped.${outputFormat})"<<endl;
		if(!((m_Parameters->getsaveSkullStripping()).empty()))
		{
			CortThickPipelineBMS<<"SetAppOption(skullstripping.output "<<m_Parameters->getsaveSkullStripping()<<")"<<endl;
			CortThickPipelineBMS<<"Set(strippedFile "<<m_Parameters->getsaveSkullStripping()<<")"<<endl;
			CortThickPipelineBMS<<"GetFileName(Name ${strippedFile} NAME_WITHOUT_EXTENSION)"<<endl;
		}
		else
			CortThickPipelineBMS<<"    SetAppOption(skullstripping.output ${strippedFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstripping.input ${labelFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstripping.skullstripping.skullstripping ${rawImage})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstripping.GM.GM "<<m_Parameters->getGMLabel()<<")"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstripping.WM.WM "<<m_Parameters->getWMLabel()<<")"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstripping.CSF.CSF "<<m_Parameters->getCSFLabel()<<")"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"If(${dilateOn} == '1')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(skullstripping.dilate '1') "<<endl;
		CortThickPipelineBMS<<"EndIf(${dilateOn})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"Run(outputskullstrip ${skullstripping})"<<endl;
		CortThickPipelineBMS<<"GetFileName(rawImageName ${strippedFile} NAME_WITHOUT_EXTENSION)"<<endl;

	/*	CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##########################Set Atlas/Parcellation Orientation######################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		


		if(strcmp( (m_Parameters->getatlasOrientation()).c_str() , (m_Parameters->getrawImageOrientation()).c_str() ) != 0)
		{
		//set atlas orientation
		CortThickPipelineBMS<<"SetAppOption(orientImage.inputVolume1 ${atlasFile})"<<endl;
		CortThickPipelineBMS<<"GetFileName(atlasName ${atlasFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"Set(orientedAtlas ${RegistrationDir}/${atlasName}_${rawImageOrientation}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.outputVolume ${orientedAtlas})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.orientation.orientation ${rawImageOrientation})"<<endl;
		CortThickPipelineBMS<<"Run(output ${orientImage})"<<endl;
		
		//set parcellation orientation
		CortThickPipelineBMS<<"SetAppOption(orientImage.inputVolume1 ${parcellationFile})"<<endl;
		CortThickPipelineBMS<<"GetFileName(parcellationName ${parcellationFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"Set(orientedParcellation ${RegistrationDir}/${parcellationName}_${rawImageOrientation}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.outputVolume ${orientedParcellation})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(orientImage.orientation.orientation ${rawImageOrientation})"<<endl;
		CortThickPipelineBMS<<"Run(output ${orientImage})"<<endl;

		CortThickPipelineBMS<<"Set(atlasFile ${orientedAtlas})"<<endl;
		CortThickPipelineBMS<<"Set(parcellationFile ${orientedParcellation})"<<endl;
		}*/


		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##########################RegisterImages : Atlas Registration#####################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"Set(atlasRegisteredFile ${RegistrationDir}/AtlasRegistered_${rawImageName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(transformFile ${RegistrationDir}/AtlasTransform_${rawImageName}.txt)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(registration.fixedImage ${strippedFile})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"GetFileName(atlasName ${atlasRegisteredFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"GetFileName(transformName ${transformFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testatlasreg ${RegistrationDir} ${atlasName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testtransform ${RegistrationDir} ${transformName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(runAtlasRegOn '0')"<<endl;
		CortThickPipelineBMS<<"If(${testatlasreg} == '' && ${testtransform} == '')"<<endl;
		CortThickPipelineBMS<<"    Set(runAtlasRegOn '1')"<<endl;
		CortThickPipelineBMS<<"EndIf(${testatlasreg})"<<endl;
		CortThickPipelineBMS<<"If(${runAtlasRegOn} == '1')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.movingImage ${atlasFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.saveTransform.saveTransform ${transformFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.resampledImage.resampledImage ${atlasRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.registration.registration 'PipelineBSpline')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(registration.initialization.initialization "<<m_Parameters->getinitialization()<<")"<<endl;
		CortThickPipelineBMS<<"    Run(outputreg ${registration})"<<endl;
		CortThickPipelineBMS<<"EndIf(${runAtlasRegOn})"<<endl;
		CortThickPipelineBMS<<""<<endl;
	
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"########################ResampleVolume2 : Parcellation registration###############################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"    Set(parcellationRegisteredFile ${RegistrationDir}/ParcellationRegistered_${rawImageName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Reference.Reference ${strippedFile})"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"GetFileName(parcelregName ${parcellationRegisteredFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"ListFileInDir(testparcelreg ${RegistrationDir} ${parcelregName}.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(runParcelRegOn '0')"<<endl;
		CortThickPipelineBMS<<"If(${testparcelreg} == '')"<<endl;
		CortThickPipelineBMS<<"    Set(runParcelRegOn '1')"<<endl;
		CortThickPipelineBMS<<"EndIf(${testparcelreg})"<<endl;
		CortThickPipelineBMS<<"If(${runParcelRegOn} == '1')"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Infile ${parcellationFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.Outfile ${parcellationRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.transformationFile.transformationFile ${transformFile})"<<endl;
		CortThickPipelineBMS<<"    SetAppOption(parcelRegistration.interpolation.interpolation 'nn')"<<endl;
		CortThickPipelineBMS<<"    Run(outputparcelreg ${parcelRegistration})"<<endl;
		CortThickPipelineBMS<<"EndIf(${runParcelRegOn})"<<endl;
		CortThickPipelineBMS<<""<<endl;
	}

	if((m_pipelineMode==5 || m_pipelineMode==6) && !m_Parameters->getsegmentationMaskingOff())
	{
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"##############ImageMath : Mask the segmentation with the parcellation registered##################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		CortThickPipelineBMS<<"SetApp(imMath @ImageMath)"<<endl;
		CortThickPipelineBMS<<"GetFileName(labelFileName ${labelFile} NAME_WITHOUT_EXTENSION)"<<endl;
		CortThickPipelineBMS<<"GetFileName(labelFilePath ${labelFile} PATH)"<<endl;
		CortThickPipelineBMS<<"Set(labelFileMasked ${labelFilePath}/${labelFileName}_masked.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imMath.infile ${labelFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imMath.mask.mask ${parcellationRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imMath.outfile.outfile ${labelFileMasked})"<<endl;
		CortThickPipelineBMS<<"Run(outputImageMath ${imMath})"<<endl;
		CortThickPipelineBMS<<"Set(labelFileNotMasked ${labelFile})"<<endl;
		CortThickPipelineBMS<<"Set(labelFile ${labelFile})"<<endl;
	}

	CortThickPipelineBMS<<"##################################################################################################"<<endl;
	CortThickPipelineBMS<<"##########################CortThick : Generate cortical thickness#################################"<<endl;
	CortThickPipelineBMS<<"##################################################################################################"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"echo('Cortical thickness computation...')"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"MakeDirectory(${CortthickDir})"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"SetAppOption(corticalThickness.outputDir ${CortthickDir})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(corticalThickness.WMLabel.WMLabel "<<m_Parameters->getWMLabel()<<")"<<endl;
	CortThickPipelineBMS<<"SetAppOption(corticalThickness.GMLabel.GMLabel "<<m_Parameters->getGMLabel()<<")"<<endl;
	if(m_pipelineMode==1 || m_pipelineMode==3 || m_pipelineMode==5)
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.inputSeg.inputSeg ${labelFile})"<<endl;
	if(m_pipelineMode==2 || m_pipelineMode==4 || m_pipelineMode==6)
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.inputSeg.inputSeg "<<m_Parameters->getlabel()<<")"<<endl;
	if(m_pipelineMode==3 || m_pipelineMode==4)
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.par.par "<<m_Parameters->getparcellationCase()<<")"<<endl;
	if(m_pipelineMode==5 || m_pipelineMode==6)
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.par.par ${parcellationRegisteredFile})"<<endl;
	if(!((m_Parameters->getWMCortThick()).empty()))
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.SaveWM.SaveWM "<<m_Parameters->getWMCortThick()<<")"<<endl;
	if(!((m_Parameters->getGMCortThick()).empty()))
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.SaveGM.SaveGM "<<m_Parameters->getGMCortThick()<<")"<<endl;
	if(!m_Parameters->getInterpOff())
	{
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.Interp '1')"<<endl;
		CortThickPipelineBMS<<"SetAppOption(corticalThickness.Threshold.Threshold "<<m_Parameters->getThreshold()<<")"<<endl;
	}
	CortThickPipelineBMS<<""<<endl;
	if(m_pipelineMode==3 || m_pipelineMode==4 || m_pipelineMode==5 || m_pipelineMode==6)
		CortThickPipelineBMS<<"    ListFileInDir(testcortthick ${CortthickDir} '*WhiteMatDistanceMap_par.csv')"<<endl;
	else
		CortThickPipelineBMS<<"    ListFileInDir(testcortthick ${CortthickDir} '*WhiteMatDistanceMap.csv')"<<endl;
	CortThickPipelineBMS<<"    Set(runCortThickOn '0')"<<endl;
	CortThickPipelineBMS<<"    If(${testcortthick} == '')"<<endl;
	CortThickPipelineBMS<<"        Set(runCortThickOn '1')"<<endl;
	CortThickPipelineBMS<<"    EndIf(${testcortthick})"<<endl;
	CortThickPipelineBMS<<"If(${runCortThickOn} == '1')"<<endl;
	CortThickPipelineBMS<<"    Run(outputCorticalThickness ${corticalThickness})"<<endl;
	CortThickPipelineBMS<<"EndIf(${runCortThickOn})"<<endl;
	CortThickPipelineBMS<<""<<endl;

	CortThickPipelineBMS<<"##################################################################################################"<<endl;
	CortThickPipelineBMS<<"######################ImageMath, ImageStat : Create CSV Spreadsheets##############################"<<endl;
	CortThickPipelineBMS<<"##################################################################################################"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"echo('Stat ...')"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"MakeDirectory(${StatDir})"<<endl;
	CortThickPipelineBMS<<""<<endl;
	CortThickPipelineBMS<<"Set(WMFile ${EMSDir}/${RefName}_WMMap.${outputFormat})"<<endl;
	CortThickPipelineBMS<<"Set(GMFile ${EMSDir}/${RefName}_GMMap.${outputFormat})"<<endl;
	CortThickPipelineBMS<<"Set(CSFFile ${EMSDir}/${RefName}_CSFMap.${outputFormat})"<<endl;

	CortThickPipelineBMS<<"SetAppOption(imageMath.infile ${labelFile})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.extractLabel.extractLabel "<<m_Parameters->getWMLabel()<<")"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.outfile.outfile ${WMFile})"<<endl;
	CortThickPipelineBMS<<"Run(outputImageMath ${imageMath})"<<endl;
	//CortThickPipelineBMS<<"Run(output '${ImageMathCmd} ${labelFile} -extractLabel 1 -outfile ${WMFile}')"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.infile ${labelFile})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.extractLabel.extractLabel "<<m_Parameters->getGMLabel()<<")"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.outfile.outfile ${GMFile})"<<endl;
	CortThickPipelineBMS<<"Run(outputImageMath ${imageMath})"<<endl;
	//CortThickPipelineBMS<<"Run(output '${ImageMathCmd} ${labelFile} -extractLabel 2 -outfile ${GMFile}')"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.infile ${labelFile})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.extractLabel.extractLabel "<<m_Parameters->getCSFLabel()<<")"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageMath.outfile.outfile ${CSFFile})"<<endl;
	CortThickPipelineBMS<<"Run(outputImageMath ${imageMath})"<<endl;
	//CortThickPipelineBMS<<"Run(output '${ImageMathCmd} ${labelFile} -extractLabel 3 -outfile ${CSFFile}')"<<endl;

	CortThickPipelineBMS<<"SetAppOption(imageStat.infile ${WMFile})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageStat.histo '1')"<<endl;
	CortThickPipelineBMS<<"Run(output ${imageStat})"<<endl;
	//CortThickPipelineBMS<<"Run(output '${ImageStatCmd} ${WMFile} -histo')"<<endl;
	CortThickPipelineBMS<<"Set(VolumeWM ${output})"<<endl;		
	CortThickPipelineBMS<<"SetAppOption(imageStat.infile ${GMFile})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageStat.histo '1')"<<endl;
	CortThickPipelineBMS<<"Run(output ${imageStat})"<<endl;
	//CortThickPipelineBMS<<"Run(output '${ImageStatCmd} ${GMFile} -histo')"<<endl;
	CortThickPipelineBMS<<"Set(VolumeGM ${output})"<<endl;	

	CortThickPipelineBMS<<"SetAppOption(imageStat.infile ${CSFFile})"<<endl;
	CortThickPipelineBMS<<"SetAppOption(imageStat.histo '1')"<<endl;
	CortThickPipelineBMS<<"Run(output ${imageStat})"<<endl;
	//CortThickPipelineBMS<<"Run(output '${ImageStatCmd} ${CSFFile} -histo')"<<endl;
	CortThickPipelineBMS<<"Set(VolumeCSF ${output})"<<endl;	
	
	CortThickPipelineBMS<<"Set(VolumeFile ${StatDir}/${RefName}_TissueSegmentationVolumes.csv)"<<endl;

	CortThickPipelineBMS<<"WriteFile(${VolumeFile} 'Volume analysis\\t(in cubic mm)\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '\\n\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '##############\\t##################\\t#\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '#   WHITE MATTER\\tWhite matter volume\\t#\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '#   GREY MATTER\\tGrey Matter volume\\t#\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '#   CSF\\tCerebrospinal fluid volume\\t#\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '##############\\t##################\\t#\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '\\n\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} 'WHITE MATTER\\t'${VolumeWM})"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} 'GREY MATTER\\t'${VolumeGM})"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} 'CSF\\t'${VolumeCSF})"<<endl;
	CortThickPipelineBMS<<"AppendFile(${VolumeFile} '\\n')"<<endl;

	if(m_pipelineMode>=3 && m_pipelineMode<=6)//if there is a parcellation map
	{
		CortThickPipelineBMS<<"Set(ParcVolumeFile ${StatDir}/${RefName}_ParcellationMapVolumes.csv)"<<endl;
		CortThickPipelineBMS<<"WriteFile(${ParcVolumeFile} 'Volume analysis\\t\\t\\t(in cubic mm)\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '\\n\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#########\\t\\t\\t###########################\\t\\t####\\t\\t####\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '# Fields:\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#   TISSUE\\t\\t\\tTissue: WM GM CSF\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#   LABEL\\t\\t\\tLabel number\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#   VOLUME\\t\\t\\tVolume in cubic mm\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#   MEAN\\t\\t\\tRelated mean intensity in skull-stripped image\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#   STD\\t\\t\\tStandard deviation of those voxels\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '#########\\t\\t\\t###########################\\t\\t####\\t\\t####\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} '\\n\\n')"<<endl;

		CortThickPipelineBMS<<"Set(ParcellationWM ${RegistrationDir}/${RefName}_parcellation_WM.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(ParcellationGM ${RegistrationDir}/${RefName}_parcellation_GM.${outputFormat})"<<endl;
		CortThickPipelineBMS<<"Set(ParcellationCSF ${RegistrationDir}/${RefName}_parcellation_CSF.${outputFormat})"<<endl;

		if(m_pipelineMode==3 || m_pipelineMode==4)
			CortThickPipelineBMS<<"Set(parcellationRegisteredFile "<<m_Parameters->getparcellationCase()<<")"<<endl;
		if(m_pipelineMode==6 || m_pipelineMode==4)
			CortThickPipelineBMS<<"Set(File ${rawImage})"<<endl;
		else
			CortThickPipelineBMS<<"Set(File ${T1File})"<<endl;

		CortThickPipelineBMS<<"SetApp(imageMath @ImageMath)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.infile ${parcellationRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.mask.mask ${WMFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.outfile.outfile ${ParcellationWM})"<<endl;
		CortThickPipelineBMS<<"Run(output ${imageMath})"<<endl;
		//CortThickPipelineBMS<<"Run(output '${ImageMathCmd} ${parcellationRegisteredFile} -mask ${WMFile} -outfile ${ParcellationWM}')"<<endl;

		CortThickPipelineBMS<<"Set(OutputFile ${RegistrationDir}/${RefName}_parcellation_WM_stat.txt)"<<endl;
		CortThickPipelineBMS<<"Set(OutputBase ${RegistrationDir}/${RefName}_parcellation_WM)"<<endl;
		
		CortThickPipelineBMS<<"SetApp(imageStat @ImageStat)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.infile ${File})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.label.label ${ParcellationWM})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.outbase.outbase ${OutputBase})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.display '1')"<<endl;
		CortThickPipelineBMS<<"Run(output ${imageStat})"<<endl;
		//CortThickPipelineBMS<<"Run(output '${ImageStatCmd} ${File} -label ${ParcellationWM} -outbase ${OutputBase}')"<<endl;
		//CortThickPipelineBMS<<"Run(output 'tail -n +11 '${OutputFile})"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} 'WM\\n')"<<std::endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} 'LABEL\\t\\t\\tVOLUME\\t\\tMEAN\\t\\tSTD\\n')"<<std::endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} ${output}'\\n\\n')"<<std::endl;
	
		CortThickPipelineBMS<<"SetApp(imageMath @ImageMath)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.infile ${parcellationRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.mask.mask ${GMFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.outfile.outfile ${ParcellationGM})"<<endl;
		CortThickPipelineBMS<<"Run(output ${imageMath})"<<endl;
		//CortThickPipelineBMS<<"Run(output '${ImageMathCmd} ${parcellationRegisteredFile} -mask ${GMFile} -outfile ${ParcellationGM}')"<<endl;

		CortThickPipelineBMS<<"Set(OutputFile ${RegistrationDir}/${RefName}_parcellation_GM_stat.txt)"<<endl;
		CortThickPipelineBMS<<"Set(OutputBase ${RegistrationDir}/${RefName}_parcellation_GM)"<<endl;
	
		CortThickPipelineBMS<<"SetApp(imageStat @ImageStat)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.infile ${File})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.label.label ${ParcellationGM})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.outbase.outbase ${OutputBase})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.display '1')"<<endl;
		CortThickPipelineBMS<<"Run(output ${imageStat})"<<endl;
		//CortThickPipelineBMS<<"Run(output '${ImageStatCmd} ${File} -label ${ParcellationGM} -outbase ${OutputBase}')"<<endl;
		//CortThickPipelineBMS<<"Run(output 'tail -n +11 '${OutputFile})"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} 'GM\\n')"<<std::endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} 'LABEL\\t\\t\\tVOLUME\\t\\tMEAN\\t\\tSTD\\n')"<<std::endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} ${output}'\\n\\n')"<<std::endl;
	
		CortThickPipelineBMS<<"SetApp(imageMath @ImageMath)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.infile ${parcellationRegisteredFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.mask.mask ${CSFFile})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageMath.outfile.outfile ${ParcellationCSF})"<<endl;
		CortThickPipelineBMS<<"Run(output ${imageMath})"<<endl;
		//CortThickPipelineBMS<<"Run(output '${ImageMathCmd} ${parcellationRegisteredFile} -mask ${CSFFile} -outfile ${ParcellationCSF}')"<<endl;

		CortThickPipelineBMS<<"Set(OutputFile ${RegistrationDir}/${RefName}_parcellation_CSF_stat.txt)"<<endl;
		CortThickPipelineBMS<<"Set(OutputBase ${RegistrationDir}/${RefName}_parcellation_CSF)"<<endl;
	
		CortThickPipelineBMS<<"SetApp(imageStat @ImageStat)"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.infile ${File})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.label.label ${ParcellationCSF})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.outbase.outbase ${OutputBase})"<<endl;
		CortThickPipelineBMS<<"SetAppOption(imageStat.display '1')"<<endl;
		CortThickPipelineBMS<<"Run(output ${imageStat})"<<endl;
		//CortThickPipelineBMS<<"Run(output '${ImageStatCmd} ${File} -label ${ParcellationCSF} -outbase ${OutputBase}')"<<endl;
		//CortThickPipelineBMS<<"Run(output 'tail -n +11 '${OutputFile})"<<endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} 'CSF\\n')"<<std::endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} 'LABEL\\t\\t\\tVOLUME\\t\\tMEAN\\t\\tSTD\\n')"<<std::endl;
		CortThickPipelineBMS<<"AppendFile(${ParcVolumeFile} ${output}'\\n\\n')"<<std::endl;
	}	


	CortThickPipelineBMS<<"##################################################################################################"<<endl;
	CortThickPipelineBMS<<"###############################Output Directory Organization######################################"<<endl;
	CortThickPipelineBMS<<"##################################################################################################"<<endl;
	CortThickPipelineBMS<<""<<endl;

	//save in the outputdir
	if(!((m_Parameters->getprefix()).empty()))
	{
	//copy the itkEMS files
		if(m_pipelineMode==1 || m_pipelineMode==3 || m_pipelineMode==5)
		{
			if(m_pipelineMode==5  && !m_Parameters->getsegmentationMaskingOff())
			{
				CortThickPipelineBMS<<"CopyFile(${labelFile} ${EMSDir}/${RefName}_labels_EMS_masked.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${labelFile})"<<endl;
				CortThickPipelineBMS<<"Set(labelFile ${EMSDir}/${RefName}_labels_EMS_masked.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${labelFileNotMasked} ${EMSDir}/${RefName}_labels_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${labelFileNotMasked})"<<endl;
				CortThickPipelineBMS<<"Set(labelFileNotMasked ${EMSDir}/${RefName}_labels_EMS.${outputFormat})"<<endl;
			}
			else
			{
				CortThickPipelineBMS<<"CopyFile(${labelFile} ${EMSDir}/${RefName}_labels_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${labelFile})"<<endl;
				CortThickPipelineBMS<<"Set(labelFile ${EMSDir}/${RefName}_labels_EMS.${outputFormat})"<<endl;
			}
			CortThickPipelineBMS<<"Set(CorrectedEMS ${EMSDir}/${RefName}_T1_corrected_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"If(${Mode} == 2 || ${Mode} == 4)"<<endl;
			CortThickPipelineBMS<<"   If(${atlasType} == T2)"<<endl;
			CortThickPipelineBMS<<"      Set(CorrectedEMS ${EMSDir}/${RefName}_T2_corrected_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"   EndIf(${atlasType})"<<endl;
			CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(Posterior0 ${EMSDir} *posterior0_EMS*)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${Posterior0})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${EMSDir}/${File} ${EMSDir}/${RefName}_posterior0_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(Posterior1 ${EMSDir} *posterior1_EMS*)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${Posterior1})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${EMSDir}/${File} ${EMSDir}/${RefName}_posterior1_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(Posterior2 ${EMSDir} *posterior2_EMS*)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${Posterior2})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${EMSDir}/${File} ${EMSDir}/${RefName}_posterior2_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
			if(!((m_Parameters->getT1()).empty()))
			{
				CortThickPipelineBMS<<"CopyFile(${EMSDir}/${T1OriginalName}_registered_EMS.${outputFormat} ${EMSDir}/${RefName}_T1_registered_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${T1OriginalName}_registered_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${T1correctedFile} ${EMSDir}/${RefName}_T1_corrected_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${T1correctedFile})"<<endl;
				CortThickPipelineBMS<<"Set(T1correctedFile ${EMSDir}/${RefName}_T1_corrected_EMS.${outputFormat})"<<endl;
			}
			if(!((m_Parameters->getT2()).empty()))
			{
				CortThickPipelineBMS<<"CopyFile(${EMSDir}/${T2OriginalName}_registered_EMS.${outputFormat} ${EMSDir}/${RefName}_T2_registered_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${T2OriginalName}_registered_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${T2correctedFile} ${EMSDir}/${RefName}_T2_corrected_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${T2correctedFile})"<<endl;
				CortThickPipelineBMS<<"Set(T2correctedFile ${EMSDir}/${RefName}_T2_corrected_EMS.${outputFormat})"<<endl;
			}
			if(!((m_Parameters->getpd()).empty()))
			{
				CortThickPipelineBMS<<"CopyFile(${EMSDir}/${pdOriginalName}_registered_EMS.${outputFormat} ${EMSDir}/${RefName}_PD_registered_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${pdOriginalName}_registered_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${pdcorrectedFile} ${EMSDir}/${RefName}_PD_corrected_EMS.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${pdcorrectedFile})"<<endl;
				CortThickPipelineBMS<<"Set(pdcorrectedFile ${EMSDir}/${RefName}_PD_corrected_EMS.${outputFormat})"<<endl;
			}
			CortThickPipelineBMS<<"ListFileInDir(TemplateAffine ${EMSDir} *template_affine_EMS*)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${TemplateAffine})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${EMSDir}/${File} ${EMSDir}/${RefName}_template_affine_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(TemplateWarped ${EMSDir} *template_warped_EMS*)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${TemplateWarped})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${EMSDir}/${File} ${EMSDir}/${RefName}_template_warped_EMS.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(AffineFile ${EMSDir} *.affine)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${AffineFile})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(BsplineFile ${EMSDir} *.bspline)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${BsplineFile})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${EMSDir}/${File} ${EMSDir}/${RefName}_template_EMS.bspline)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${EMSDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
		}
	//copy the skullstripped file
		if(m_pipelineMode==5)
		{
			if(!((m_Parameters->getT1()).empty()))
			{
				CortThickPipelineBMS<<"CopyFile(${T1strippedFile} ${RegistrationDir}/${RefName}_T1_stripped.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${T1strippedFile})"<<endl;
				CortThickPipelineBMS<<"Set(SkullStripped ${RegistrationDir}/${RefName}_T1_stripped.${outputFormat})"<<endl;
			}
			if(!((m_Parameters->getT2()).empty()))
			{
				CortThickPipelineBMS<<"CopyFile(${T2strippedFile} ${RegistrationDir}/${RefName}_T2_stripped.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${T2strippedFile})"<<endl;
				CortThickPipelineBMS<<"Set(SkullStripped ${RegistrationDir}/${RefName}_T2_stripped.${outputFormat})"<<endl;
			}
			if(!((m_Parameters->getpd()).empty()))
			{
				CortThickPipelineBMS<<"CopyFile(${pdstrippedFile} ${RegistrationDir}/${RefName}_PD_stripped.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${pdstrippedFile})"<<endl;
				CortThickPipelineBMS<<"Set(SkullStripped ${RegistrationDir}/${RefName}_PD_stripped.${outputFormat})"<<endl;
			}

		}
		if(m_pipelineMode==6)
		{
			CortThickPipelineBMS<<"CopyFile(${strippedFile} ${RegistrationDir}/${RefName}_stripped.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${strippedFile})"<<endl;
			CortThickPipelineBMS<<"Set(SkullStripped ${RegistrationDir}/${RefName}_stripped.${outputFormat})"<<endl;
			if(!m_Parameters->getsegmentationMaskingOff())
			{
				CortThickPipelineBMS<<"CopyFile(${labelFile} ${EMSDir}/${RegistrationDir}_labels_masked.${outputFormat})"<<endl;
				CortThickPipelineBMS<<"Set(labelFile ${RegistrationDir}/${RefName}_labels_masked.${outputFormat})"<<endl;
			}
		}
	//copy the registered files (atlas+parcellation)
		if(m_pipelineMode==6 || m_pipelineMode==5)
		{
			CortThickPipelineBMS<<"CopyFile(${atlasRegisteredFile} ${RegistrationDir}/${RefName}_AtlasRegistered.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${atlasRegisteredFile})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${transformFile} ${RegistrationDir}/${RefName}_TransformFile.txt)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${transformFile})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${parcellationRegisteredFile} ${RegistrationDir}/${RefName}_ParcellationRegistered.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${parcellationRegisteredFile})"<<endl;
			CortThickPipelineBMS<<"Set(AtlasRegistered ${RegistrationDir}/${RefName}_AtlasRegistered.${outputFormat})"<<endl;
			CortThickPipelineBMS<<"Set(ParcellationRegistered ${RegistrationDir}/${RefName}_ParcellationRegistered.${outputFormat})"<<endl;
		}
	//copy the cortthick directory files
		if(!((m_Parameters->getWMCortThick()).empty()))
		{
			CortThickPipelineBMS<<"CopyFile("<<m_Parameters->getWMCortThick()<<" ${CortthickDir}/${RefName}_AvgCortThickOnWMBoundary.${outputFormat})"<<endl;
		}
		if(!((m_Parameters->getGMCortThick()).empty()))
		{
			CortThickPipelineBMS<<"CopyFile("<<m_Parameters->getGMCortThick()<<" ${CortthickDir}/${RefName}_AvgCortThicOnGMBoundary.${outputFormat})"<<endl;
		}
		if(m_pipelineMode==1 || m_pipelineMode==2)
		{
			CortThickPipelineBMS<<"ListFileInDir(WhiteMapDistanceMap ${CortthickDir} *WhiteMatDistanceMap.csv)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${WhiteMapDistanceMap})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${CortthickDir}/${File} ${CortthickDir}/${RefName}_WhiteMatDistanceMap.csv)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(File)"<<endl;
		}
		if(m_pipelineMode>=3 && m_pipelineMode<=6)
		{
			CortThickPipelineBMS<<"CopyFile(${CortthickDir}/parcellationOnWhiteInterp.nrrd ${CortthickDir}/${RefName}_parcellationOnWhiteInterp.nrrd)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/parcellationOnWhiteInterp.nrrd)"<<endl;
			CortThickPipelineBMS<<"CopyFile(${CortthickDir}/parcellationOnWhite.nrrd ${CortthickDir}/${RefName}_parcellationOnWhite.nrrd)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/parcellationOnWhite.nrrd)"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(WhiteMapDistanceMap ${CortthickDir} *WhiteMatDistanceMap_par.csv)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${WhiteMapDistanceMap})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${CortthickDir}/${File} ${CortthickDir}/${RefName}_WhiteMatDistanceMap_par.csv)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(${File})"<<endl;
			CortThickPipelineBMS<<"ListFileInDir(WhiteMapDistanceMap ${CortthickDir} *WhiteMatDistanceMap_par_array.csv)"<<endl;
			CortThickPipelineBMS<<"ForEach(File ${WhiteMapDistanceMap})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${CortthickDir}/${File} ${CortthickDir}/${RefName}_WhiteMatDistanceMap_par_array.csv)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/${File})"<<endl;
			CortThickPipelineBMS<<"EndForEach(${File})"<<endl;
			CortThickPipelineBMS<<"CopyFile(${CortthickDir}/parcellationOnWhite.nrrd ${CortthickDir}/${RefName}_parcellationOnWhite.nrrd)"<<endl;
			CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/parcellationOnWhite.nrrd)"<<endl;
			if(!m_Parameters->getInterpOff())
			{
				CortThickPipelineBMS<<"ListFileInDir(WhiteMapDistanceMapInterp ${CortthickDir} *WhiteMatDistanceMap_par_interp.csv)"<<endl;
				CortThickPipelineBMS<<"ForEach(File ${WhiteMapDistanceMapInterp})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${CortthickDir}/${File} ${CortthickDir}/${RefName}_WhiteMatDistanceMap_par_interp.csv)"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/${File})"<<endl;
				CortThickPipelineBMS<<"EndForEach(${File})"<<endl;
				CortThickPipelineBMS<<"ListFileInDir(WhiteMapDistanceMapInterp ${CortthickDir} *WhiteMatDistanceMap_par_array_interp.csv)"<<endl;
				CortThickPipelineBMS<<"ForEach(File ${WhiteMapDistanceMapInterp})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${CortthickDir}/${File} ${CortthickDir}/${RefName}_WhiteMatDistanceMap_par_array_interp.csv)"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/${File})"<<endl;
				CortThickPipelineBMS<<"EndForEach(${File})"<<endl;
				CortThickPipelineBMS<<"CopyFile(${CortthickDir}/parcellationOnWhiteInterp.nrrd ${CortthickDir}/${RefName}_parcellationOnWhiteInterp.nrrd)"<<endl;
				CortThickPipelineBMS<<"DeleteFile(${CortthickDir}/parcellationOnWhiteInterp.nrrd)"<<endl;
			}
		}

	//remove files with Node (Slicer3 use)
		CortThickPipelineBMS<<"ListDirInDir(ARCTICFolder ${pipelineDir} *)"<<endl;
		CortThickPipelineBMS<<"ForEach(Folder ${ARCTICFolder})"<<endl;
		//CortThickPipelineBMS<<"If(${Folder} != 'CorticalThickness')"<<endl;
		CortThickPipelineBMS<<"	 ListFileInDir(NodeFile ${pipelineDir}/${Folder} *ScalarVolumeNode*)"<<endl;
		CortThickPipelineBMS<<"	 If(${NodeFile} != '')"<<endl;
		CortThickPipelineBMS<<"	 ForEach(File ${pipelineDir}/${Folder}/${NodeFile})"<<endl;
		CortThickPipelineBMS<<"	   DeleteFile(${File})"<<endl;
		CortThickPipelineBMS<<"  EndForEach(File)"<<endl;
		CortThickPipelineBMS<<"	 EndIf(${NodeFile})"<<endl;
		//CortThickPipelineBMS<<"EndIf(${Folder}}
		CortThickPipelineBMS<<"EndForEach(${Folder})"<<endl;
	}
	else
	{
		if(m_pipelineMode==1 || m_pipelineMode==3 || m_pipelineMode==5)
		{
			CortThickPipelineBMS<<"Set(CorrectedEMS ${T1correctedFile})"<<endl;
			CortThickPipelineBMS<<"If(${Mode} == '2' || ${Mode} == '4')"<<endl;
			CortThickPipelineBMS<<"   If(${atlasType} == 'T2')"<<endl;
			CortThickPipelineBMS<<"      Set(CorrectedEMS ${T2correctedFile})"<<endl;
			CortThickPipelineBMS<<"   EndIf(${atlasType})"<<endl;
			CortThickPipelineBMS<<"EndIf(${Mode})"<<endl;
		}
		if(m_pipelineMode==5)
		{			
			if(m_Parameters->getatlasType()=="T1")
				CortThickPipelineBMS<<"Set(SkullStripped ${T1strippedFile})"<<endl;
			if(m_Parameters->getatlasType()=="T2")
				CortThickPipelineBMS<<"Set(SkullStripped ${T2strippedFile})"<<endl;
			CortThickPipelineBMS<<"Set(AtlasRegistered ${atlasRegisteredFile})"<<endl;
			CortThickPipelineBMS<<"Set(ParcellationRegistered ${parcellationRegisteredFile})"<<endl;
		}
		if(m_pipelineMode==6)
		{
			CortThickPipelineBMS<<"Set(SkullStripped ${strippedFile})"<<endl;
			CortThickPipelineBMS<<"Set(AtlasRegistered ${atlasRegisteredFile})"<<endl;
			CortThickPipelineBMS<<"Set(ParcellationRegistered ${parcellationRegisteredFile})"<<endl;
		}
	}

	/*if(strcmp((m_Parameters->getatlasOrientation()).c_str(),(m_Parameters->getrawImageOrientation()).c_str())!= 0 || strcmp((m_Parameters->getatlasOrientation()).c_str(),(m_Parameters->getorientation()).c_str() ) != 0)
	{
		CortThickPipelineBMS<<"DeleteFile(${parcellationFile})"<<endl;
		CortThickPipelineBMS<<"DeleteFile(${atlasFile})"<<endl;
	}*/
	
	//set the modelmaker output files
	CortThickPipelineBMS<<"Set(WMSurface ${MeshDir}/${RefName}_WM_Surface.vtk)"<<endl;
	CortThickPipelineBMS<<"Set(GMSurface ${MeshDir}/${RefName}_GM_Surface.vtk)"<<endl;
	CortThickPipelineBMS<<"GetFileName(LabelPath ${labelFile} PATH)"<<endl;
	CortThickPipelineBMS<<"GetFileName(LabelName ${labelFile} NAME_WITHOUT_EXTENSION)"<<endl;
	if(m_pipelineMode!=2 && m_pipelineMode!=4)
	{
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<"####################################MRML Scene Creation###########################################"<<endl;
		CortThickPipelineBMS<<"##################################################################################################"<<endl;
		CortThickPipelineBMS<<""<<endl;
		WriteHeaderMRMLScene();
		if(m_pipelineMode==1 || m_pipelineMode==3)
		{
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Labels_EMS","labelFile",1,2,1);
			WriteSnapshotBeginingMRMLScene("Segmentation",1,1,"vtkMRMLScalarVolumeNode2");
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Label_EMS","labelFile",1,2,1);
			WriteSnapshotEndMRMLScene();

			WriteMeshSnapshotsMRMLScene(2,3);
		}
		if(m_pipelineMode==5)
		{			
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Labels_EMS","labelFile",1,2,1);
			WriteImageMRMLScene("Skull-stripped","SkullStripped",0,3,2);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,4,3);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,5,2);

			WriteSnapshotBeginingMRMLScene("Segmentation",1,1,"vtkMRMLScalarVolumeNode2");
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Labels_EMS","labelFile",1,2,1);
			WriteImageMRMLScene("Skull-stripped","SkullStripped",0,3,2);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,4,3);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,5,2);
			WriteSnapshotEndMRMLScene();

			WriteSnapshotBeginingMRMLScene("Skull_stripping",2,3,"");
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Labels_EMS","labelFile",1,2,1);
			WriteImageMRMLScene("Skull-stripped","SkullStripped",0,3,2);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,4,3);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,5,3);
			WriteSnapshotEndMRMLScene();

			WriteSnapshotBeginingMRMLScene("Atlas_registration",3,4,"");
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Labels_EMS","labelFile",1,2,1);
			WriteImageMRMLScene("Skull-stripped","SkullStripped",0,3,2);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,4,3);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,5,2);
			WriteSnapshotEndMRMLScene();

			WriteSnapshotBeginingMRMLScene("Parcellation_registration",4,3,"vtkMRMLScalarVolumeNode5");
			WriteImageMRMLScene("Corrected_EMS","CorrectedEMS",0,1,1);
			WriteImageMRMLScene("Labels_EMS","labelFile",1,2,1);
			WriteImageMRMLScene("Skullstripped","SkullStripped",0,3,2);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,4,3);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,5,2);
			WriteSnapshotEndMRMLScene();

			WriteMeshSnapshotsMRMLScene(5,6);

		}
		if(m_pipelineMode==6)
		{
			WriteImageMRMLScene("Labels","labelFile",1,1,1);
			WriteImageMRMLScene("Skullstripped","SkullStripped",0,2,1);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,3,2);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,4,2);

			WriteSnapshotBeginingMRMLScene("Skull_stripping",1,2,"");
			WriteImageMRMLScene("Labels","labelFile",1,1,1);
			WriteImageMRMLScene("Skull-stripped","SkullStripped",0,2,1);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,3,2);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,4,2);
			WriteSnapshotEndMRMLScene();

			WriteSnapshotBeginingMRMLScene("Atlas registration",2,3,"");
			WriteImageMRMLScene("Labels","labelFile",1,1,1);
			WriteImageMRMLScene("Skull stripped","SkullStripped",0,2,1);
			WriteImageMRMLScene("Atlas registered","AtlasRegistered",0,3,2);
			WriteImageMRMLScene("Parcellation registered","ParcellationRegistered",1,4,2);
			WriteSnapshotEndMRMLScene();

			WriteSnapshotBeginingMRMLScene("Parcellation_registration",3,2,"vtkMRMLScalarVolumeNode4");
			WriteImageMRMLScene("Labels","labelFile",1,1,1);
			WriteImageMRMLScene("Skull-stripped","SkullStripped",0,2,1);
			WriteImageMRMLScene("Atlas-registered","AtlasRegistered",0,3,2);
			WriteImageMRMLScene("Parcellation-registered","ParcellationRegistered",1,4,2);
			WriteSnapshotEndMRMLScene();

			WriteMeshSnapshotsMRMLScene(4,5);

		}
		WriteEndMRMLScene();
	}
	CortThickPipelineBMS.close();

}

void Computation::WriteBMSFile2()
{
	CortThickPipelineBMS2.open((m_Parameters->getBMSFile2()).c_str());

	CortThickPipelineBMS2<<"##################################################################################################"<<endl;
	CortThickPipelineBMS2<<"##########################ModelMaker : Generate GM and WM volumes#################################"<<endl;
	CortThickPipelineBMS2<<"##################################################################################################"<<endl;

	CortThickPipelineBMS2<<"SetApp(volumeCreation @ModelMaker)"<<endl;

	CortThickPipelineBMS2<<"Set(dataDir '"<<m_Parameters->getoutputDir()<<"')"<<endl;
	CortThickPipelineBMS2<<"Set(pipelineDir ${dataDir}/ARCTIC)"<<endl;
	CortThickPipelineBMS2<<"Set(MeshDir ${pipelineDir}/Mesh)"<<endl;
	CortThickPipelineBMS2<<"MakeDirectory(${MeshDir})"<<endl;

	CortThickPipelineBMS2<<"Set(RefName '"<<m_refName<<"')"<<endl;
	CortThickPipelineBMS2<<"Set(labelFile '"<<m_labelFileCentered<<"')"<<endl;
	CortThickPipelineBMS2<<"GetFileName(LabelPath ${labelFile} PATH)"<<endl;
	CortThickPipelineBMS2<<"GetFileName(LabelName ${labelFile} NAME_WITHOUT_EXTENSION)"<<endl;	
	CortThickPipelineBMS2<<"SetAppOption(volumeCreation.labelFile ${labelFile})"<<endl;
	CortThickPipelineBMS2<<"SetAppOption(volumeCreation.generateAll '1')"<<endl;
	CortThickPipelineBMS2<<"Set(tmpMRMLFile ${LabelPath}/tmp.mrml)"<<endl;
 	CortThickPipelineBMS2<<"WriteFile(${tmpMRMLFile} '<MRML >\\n')"<<endl;
	CortThickPipelineBMS2<<"AppendFile(${tmpMRMLFile} '</MRML>\\n')"<<endl;
	CortThickPipelineBMS2<<"SetAppOption(volumeCreation.modelSceneFile.modelSceneFile ${tmpMRMLFile})"<<endl;
	CortThickPipelineBMS2<<"SetAppOption(volumeCreation.jointsmooth 1)"<<endl;
	CortThickPipelineBMS2<<"Run(outputModelMaker ${volumeCreation})"<<endl;

	CortThickPipelineBMS2<<"CopyFile(${LabelPath}/Model_1.vtk ${MeshDir}/${RefName}_WM_Surface.vtk)"<<endl;
	CortThickPipelineBMS2<<"CopyFile(${LabelPath}/Model_2.vtk ${MeshDir}/${RefName}_GM_Surface.vtk)"<<endl;

	CortThickPipelineBMS2<<"DeleteFile(${LabelPath}/Model_1.vtk)"<<endl;
	CortThickPipelineBMS2<<"DeleteFile(${LabelPath}/Model_2.vtk)"<<endl;
	CortThickPipelineBMS2<<"DeleteFile(${LabelPath}/Model_3.vtk)"<<endl;
 	CortThickPipelineBMS2<<"DeleteFile(${LabelPath}/tmp.mrml)"<<endl;
	
	CortThickPipelineBMS2.close();
}

void Computation::ExecuteBatchmakeScript()
{
	cout<<"\tExecuting BatchMake..."<<endl;
	
	vector<const char*> args;
	char * envpath = getenv("BatchmakeWrapper_Dir");
	std::string applicationPath;	
	if (! envpath) 
	{
		std::cerr<<"The environment variable 'BatchmakeWrapper_Dir' needs to be set"<<std::endl;
		std::cerr<<"bash usage : export BatchmakeWrapper_Dir=<Batchmake Wrapper Directory>"<<std::endl;
		std::cerr<<"tcsh usage : setenv BatchmakeWrapper_Dir <Batchmake Wrapper Directory>"<<std::endl;
		exit(0);
	}
	else
		applicationPath = string(envpath);
	m_Parser.LoadWrappedApplication( applicationPath.c_str());
	bool ComputationSuccess = m_Parser.Execute((m_Parameters->getBMSFile()).c_str());
	
	if (!ComputationSuccess)
		cerr<<"\tExecuting BatchMake: Error!"<<endl;
}

void Computation::ExecuteBatchmakeScript2()
{
	vector<const char*> args;
	char * envpath = getenv("BatchmakeWrapper_Dir");
	std::string applicationPath;
	applicationPath = string(envpath);
	m_Parser2.LoadWrappedApplication( applicationPath );
	bool ComputationSuccess = m_Parser2.Execute((m_Parameters->getBMSFile2()).c_str());
	cout<<"\tExecuting BatchMake: Done!"<<endl;
	if (ComputationSuccess)
	  cout<<"\tExecuting BatchMake: Done!"<<endl;
	else
	  cerr<<"\tExecuting BatchMake: Error!"<<endl;
}

void Computation::WriteHeaderMRMLScene()
{
	CortThickPipelineBMS<<"Set(MRMLfile ${pipelineDir}/${RefName}_ARCTIC-Scene.mrml)"<<endl;
	CortThickPipelineBMS<<"WriteFile(${MRMLfile} '<MRML userTags=\"\">\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<Selection\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} 'id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  activeVolumeID=\"vtkMRMLScalarVolumeNode1\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"NULL\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<Interaction\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<View\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} 'id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"true\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"true\" backgroundColor=\"0.70196 0.70196 0.90588\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<Camera\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  position=\"0 500 0\" focalPoint=\"0 0 0\" viewUp=\"0 0 1\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	<Slice\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	      id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  fieldOfView=\"385.443 258.983 1.01562\" dimensions=\"445 299 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 -0 0 129.492 -0 -0 1 97.4995 0 1 0 129.492 -0 0 -0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"0\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	<SliceComposite\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '		       id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\" crosshairMode=\"ShowIntersection\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	<Slice\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	      id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  fieldOfView=\"289.673 193.983 1.01562\" dimensions=\"445 298 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 -0 129.492 0 1 0 96.9917 -0 0 1 129.999 0 -0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"0\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	<SliceComposite\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '		       id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" crosshairMode=\"ShowIntersection\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	<Slice\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	      id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  fieldOfView=\"385.443 258.983 1.01562\" dimensions=\"445 299 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-0 0 1 129.999 -1 -0 0 96.9917 -0 1 -0 129.492 0 -0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"0\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '	<SliceComposite\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" crosshairMode=\"ShowIntersection\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n')"<<endl;
}

void Computation::WriteEndMRMLScene()
{
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '</MRML>\\n')"<<endl;
}


void Computation::WriteSnapshotBeginingMRMLScene(string SceneName,int SceneSnapshotNumber,int ActiveVolumeNumber,string ActiveVolumeLabelID)
{
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<SceneSnapshot\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLSceneSnapshotNode"<<SceneSnapshotNumber<<"\" name=\""<<SceneName<<"\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\" > <Selection\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSelectionNode"<<SceneSnapshotNumber<<"\"   name=\"vtkMRMLSelectionNode"<<SceneSnapshotNumber<<"\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    activeVolumeID=\"vtkMRMLScalarVolumeNode"<<ActiveVolumeNumber<<"\"   activeLabelVolumeID=\""<<ActiveVolumeLabelID<<"\"   activeFiducialListID=\"NULL\"   activeROIListID=\"NULL\"   activeCameraID=\"NULL\"   activeViewID=\"NULL\"   activeLayoutID=\"vtkMRMLLayoutNode1\" ></Selection>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Slice\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSliceNode1\"   name=\"Green\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    fieldOfView=\"173.966 166.563 1.01563\"   dimensions=\"282 270 1\"   activeSlice=\"0\"   layoutGridRows=\"1\"   layoutGridColumns=\"1\"   sliceToRAS=\"-1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1\"   layoutName=\"Green\"   orientation=\"Coronal\"   jumpMode=\"0\"   sliceVisibility=\"true\"   widgetVisibility=\"false\"   useLabelOutline=\"false\"   sliceSpacingMode=\"0\"   prescribedSliceSpacing=\"1 1 1\" ></Slice>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Slice\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSliceNode2\"   name=\"Red\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    fieldOfView=\"200.485 191.954 1.01562\"   dimensions=\"282 270 1\"   activeSlice=\"0\"   layoutGridRows=\"1\"   layoutGridColumns=\"1\"   sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\"   layoutName=\"Red\"   orientation=\"Axial\"   jumpMode=\"0\"   sliceVisibility=\"true\"   widgetVisibility=\"false\"   useLabelOutline=\"false\"   sliceSpacingMode=\"0\"   prescribedSliceSpacing=\"1 1 1\" ></Slice>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Slice\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSliceNode3\"   name=\"Yellow\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    fieldOfView=\"191.954 183.786 1.01562\"   dimensions=\"282 270 1\"   activeSlice=\"0\"   layoutGridRows=\"1\"   layoutGridColumns=\"1\"   sliceToRAS=\"0 0 1 0 -1 0 0 0 0 1 0 0 0 0 0 1\"   layoutName=\"Yellow\"   orientation=\"Sagittal\"   jumpMode=\"0\"   sliceVisibility=\"true\"   widgetVisibility=\"false\"   useLabelOutline=\"false\"   sliceSpacingMode=\"0\"   prescribedSliceSpacing=\"1 1 1\" ></Slice>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <SliceComposite\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSliceCompositeNode3\"   name=\"vtkMRMLSliceCompositeNode3\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    backgroundVolumeID=\"vtkMRMLScalarVolumeNode"<<ActiveVolumeNumber<<"\"   foregroundVolumeID=\"\"   labelVolumeID=\""<<ActiveVolumeLabelID<<"\"   compositing=\"0\"   labelOpacity=\"0.35\"   linkedControl=\"1\"   foregroundGrid=\"0\"   backgroundGrid=\"0\"   labelGrid=\"1\"   fiducialVisibility=\"1\"   fiducialLabelVisibility=\"1\"   sliceIntersectionVisibility=\"0\"   layoutName=\"Yellow\"   annotationMode=\"All\"   crosshairMode=\"ShowIntersection\"   crosshairBehavior=\"Normal\"   crosshairThickness=\"Fine\" ></SliceComposite>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <SliceComposite\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSliceCompositeNode2\"   name=\"vtkMRMLSliceCompositeNode2\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    backgroundVolumeID=\"vtkMRMLScalarVolumeNode"<<ActiveVolumeNumber<<"\"   foregroundVolumeID=\"\"   labelVolumeID=\""<<ActiveVolumeLabelID<<"\"   compositing=\"0\"   labelOpacity=\"0.35\"   linkedControl=\"1\"   foregroundGrid=\"0\"   backgroundGrid=\"0\"   labelGrid=\"1\"   fiducialVisibility=\"1\"   fiducialLabelVisibility=\"1\"   sliceIntersectionVisibility=\"0\"   layoutName=\"Red\"   annotationMode=\"All\"   crosshairMode=\"ShowIntersection\"   crosshairBehavior=\"Normal\"   crosshairThickness=\"Fine\" ></SliceComposite>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <SliceComposite\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSliceCompositeNode1\"   name=\"vtkMRMLSliceCompositeNode1\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    backgroundVolumeID=\"vtkMRMLScalarVolumeNode"<<ActiveVolumeNumber<<"\"   foregroundVolumeID=\"\"   labelVolumeID=\""<<ActiveVolumeLabelID<<"\"   compositing=\"0\"   labelOpacity=\"0.35\"   linkedControl=\"1\"   foregroundGrid=\"0\"   backgroundGrid=\"0\"   labelGrid=\"1\"   fiducialVisibility=\"1\"   fiducialLabelVisibility=\"1\"   sliceIntersectionVisibility=\"0\"   layoutName=\"Green\"   annotationMode=\"All\"   crosshairMode=\"ShowIntersection\"   crosshairBehavior=\"Normal\"   crosshairThickness=\"Fine\" ></SliceComposite>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ScriptedModule\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLScriptedModuleNode1\"   name=\"vtkMRMLScriptedModuleNode1\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"  parameter0= \"label 1\" ></ScriptedModule>\\n')"<<endl;
}

void Computation::WriteImageMRMLScene(string ImageName,string ImagePath,bool LabelFlag,int VolumeNodeNumber,int VolumeDisplayNodeNumber)
{
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<VolumeArchetypeStorage\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLVolumeArchetypeStorageNode"<<VolumeNodeNumber<<"\" name=\"vtkMRMLVolumeArchetypeStorageNode"<<VolumeNodeNumber<<"\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  fileName=\"'${"<<ImagePath<<"}'\" useCompression=\"1\" readState=\"0\" writeState=\"0\" centerImage=\"1\" singleFile=\"0\" UseOrientationFromFile=\"1\"></VolumeArchetypeStorage>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<Volume\\n')"<<endl;
	if(!LabelFlag)
	{
		CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLScalarVolumeNode"<<VolumeNodeNumber<<"\" name=\""<<ImageName<<"\" hideFromEditors=\"false\" selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLVolumeArchetypeStorageNode"<<VolumeNodeNumber<<"\" userTags=\"\" displayNodeRef=\"vtkMRMLScalarVolumeDisplayNode"<<VolumeDisplayNodeNumber<<"\"  labelMap=\"0\"></Volume>\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<VolumeDisplay\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLScalarVolumeDisplayNode"<<VolumeDisplayNodeNumber<<"\" name=\"vtkMRMLScalarVolumeDisplayNode"<<VolumeDisplayNodeNumber<<"\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" colorNodeRef=\"vtkMRMLColorTableNodeGrey\" interpolate=\"1\" autoWindowLevel=\"1\" applyThreshold=\"0\" autoThreshold=\"0\"></VolumeDisplay>\\n')"<<endl;
	}
	else
	{
		CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLScalarVolumeNode"<<VolumeNodeNumber<<"\" name=\""<<ImageName<<"\" hideFromEditors=\"false\" selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLVolumeArchetypeStorageNode"<<VolumeNodeNumber<<"\" userTags=\"\" displayNodeRef=\"vtkMRMLLabelMapVolumeDisplayNode"<<VolumeDisplayNodeNumber<<"\"  labelMap=\"1\"></Volume>\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<LabelMapVolumeDisplay\\n')"<<endl;
		CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLLabelMapVolumeDisplayNode"<<VolumeDisplayNodeNumber<<"\" name=\"vtkMRMLLabelMapVolumeDisplayNode"<<VolumeDisplayNodeNumber<<"\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"0.35\" visibility=\"true\" clipping=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeLabels\" ></LabelMapVolumeDisplay>\\n')"<<endl;
	}
}

void Computation::WriteSnapshotEndMRMLScene()
{
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '</SceneSnapshot>\\n')"<<endl;
}

void Computation::WriteMeshSnapshotsMRMLScene(int VolumeNodeNumber1, int VolumeNodeNumber2)
{
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelHierarchyNode1\" name=\"Model Maker Model1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  displayNodeRef=\"vtkMRMLModelDisplayNode1\" expanded=\"true\"></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelDisplayNode1\" name=\"vtkMRMLModelDisplayNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelHierarchyNode2\" name=\"vtkMRMLModelHierarchyNode2\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  parentNodeRef=\"vtkMRMLModelHierarchyNode1\" modelNodeRef=\"vtkMRMLModelNode1\" expanded=\"true\"></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<Model\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelNode1\" name=\"WM-Surface\" hideFromEditors=\"false\" selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLModelStorageNode1\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode2\"></Model>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelDisplayNode2\" name=\"vtkMRMLModelDisplayNode2\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  color=\"0.2 0.501961 0.8\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"false\" clipping=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelHierarchyNode3\" name=\"vtkMRMLModelHierarchyNode3\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  parentNodeRef=\"vtkMRMLModelHierarchyNode1\" modelNodeRef=\"vtkMRMLModelNode2\" expanded=\"true\"></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<Model\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelNode2\" name=\"GM-Surface\" hideFromEditors=\"false\" selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLModelStorageNode2\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode3\"></Model>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelDisplayNode3\" name=\"vtkMRMLModelDisplayNode3\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  color=\"1 0.8 0.701961\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"false\" clipping=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<SceneSnapshot\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLSceneSnapshotNode"<<VolumeNodeNumber1<<"\" name=\"WM-Mesh\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\" > <Selection\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSelectionNode"<<VolumeNodeNumber1<<"\"   name=\"vtkMRMLSelectionNode"<<VolumeNodeNumber1<<"\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    activeVolumeID=\"NULL\"   activeLabelVolumeID=\"NULL\"   activeFiducialListID=\"NULL\"   activeROIListID=\"NULL\"   activeCameraID=\"NULL\"   activeViewID=\"NULL\"   activeLayoutID=\"vtkMRMLLayoutNode1\" ></Selection>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelHierarchyNode1\"   name=\"Model Maker Model1\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    displayNodeRef=\"vtkMRMLModelDisplayNode1\"   expanded=\"true\" ></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelDisplayNode1\"   name=\"vtkMRMLModelDisplayNode1\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    color=\"0.5 0.5 0.5\"   selectedColor=\"1 0 0\"   selectedAmbient=\"0.4\"   ambient=\"0\"   diffuse=\"1\"   selectedSpecular=\"0.5\"   specular=\"0\"   power=\"1\"   opacity=\"1\"   visibility=\"true\"   clipping=\"false\"   backfaceCulling=\"true\"   scalarVisibility=\"false\"   vectorVisibility=\"false\"   tensorVisibility=\"false\"   autoScalarRange=\"true\"   scalarRange=\"0 100\"   colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\"   activeScalarName=\"\"  ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelHierarchyNode2\"   name=\"vtkMRMLModelHierarchyNode2\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    parentNodeRef=\"vtkMRMLModelHierarchyNode1\"   modelNodeRef=\"vtkMRMLModelNode1\"   expanded=\"true\" ></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Model\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelNode1\"   name=\"WM-Surface\"   hideFromEditors=\"false\"   selectable=\"true\"    selected=\"false\"    userTags=\"\"   displayNodeRef=\"vtkMRMLModelDisplayNode2\" ></Model>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelDisplayNode2\"   name=\"vtkMRMLModelDisplayNode2\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    color=\"0.2 0.501961 0.8\"   selectedColor=\"1 0 0\"   selectedAmbient=\"0.4\"   ambient=\"0\"   diffuse=\"1\"   selectedSpecular=\"0.5\"   specular=\"0\"   power=\"1\"   opacity=\"1\"   visibility=\"true\"   clipping=\"false\"   backfaceCulling=\"true\"   scalarVisibility=\"false\"   vectorVisibility=\"false\"   tensorVisibility=\"false\"   autoScalarRange=\"true\"   scalarRange=\"0 100\"  ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelHierarchyNode3\"   name=\"vtkMRMLModelHierarchyNode3\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    parentNodeRef=\"vtkMRMLModelHierarchyNode1\"   modelNodeRef=\"vtkMRMLModelNode2\"   expanded=\"true\" ></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Model\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelNode2\"   name=\"GM-Surface\"   hideFromEditors=\"false\"   selectable=\"true\"    selected=\"false\"    userTags=\"\"   displayNodeRef=\"vtkMRMLModelDisplayNode3\" ></Model>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelDisplayNode3\"   name=\"vtkMRMLModelDisplayNode3\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    color=\"1 0.8 0.701961\"   selectedColor=\"1 0 0\"   selectedAmbient=\"0.4\"   ambient=\"0\"   diffuse=\"1\"   selectedSpecular=\"0.5\"   specular=\"0\"   power=\"1\"   opacity=\"1\"   visibility=\"false\"   clipping=\"false\"   backfaceCulling=\"true\"   scalarVisibility=\"false\"   vectorVisibility=\"false\"   tensorVisibility=\"false\"   autoScalarRange=\"true\"   scalarRange=\"0 100\"  ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '</SceneSnapshot>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<SceneSnapshot\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLSceneSnapshotNode"<<VolumeNodeNumber2<<"\" name=\"GM-Mesh\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\" > <Selection\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLSelectionNode"<<VolumeNodeNumber2<<"\"   name=\"vtkMRMLSelectionNode"<<VolumeNodeNumber2<<"\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    activeVolumeID=\"NULL\"   activeLabelVolumeID=\"NULL\"   activeFiducialListID=\"NULL\"   activeROIListID=\"NULL\"   activeCameraID=\"NULL\"   activeViewID=\"NULL\"   activeLayoutID=\"vtkMRMLLayoutNode1\" ></Selection>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelHierarchyNode1\"   name=\"Model Maker Model1\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    displayNodeRef=\"vtkMRMLModelDisplayNode1\"   expanded=\"true\" ></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelDisplayNode1\"   name=\"vtkMRMLModelDisplayNode1\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    color=\"0.5 0.5 0.5\"   selectedColor=\"1 0 0\"   selectedAmbient=\"0.4\"   ambient=\"0\"   diffuse=\"1\"   selectedSpecular=\"0.5\"   specular=\"0\"   power=\"1\"   opacity=\"1\"   visibility=\"true\"   clipping=\"false\"   backfaceCulling=\"true\"   scalarVisibility=\"false\"   vectorVisibility=\"false\"   tensorVisibility=\"false\"   autoScalarRange=\"true\"   scalarRange=\"0 100\"   colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\"   activeScalarName=\"\"  ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelHierarchyNode2\"   name=\"vtkMRMLModelHierarchyNode2\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    parentNodeRef=\"vtkMRMLModelHierarchyNode1\"   modelNodeRef=\"vtkMRMLModelNode1\"   expanded=\"true\" ></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Model\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelNode1\"   name=\"WM-Surface\"   hideFromEditors=\"false\"   selectable=\"true\"    selected=\"false\"    userTags=\"\"   displayNodeRef=\"vtkMRMLModelDisplayNode2\" ></Model>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelDisplayNode2\"   name=\"vtkMRMLModelDisplayNode2\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    color=\"0.2 0.501961 0.8\"   selectedColor=\"1 0 0\"   selectedAmbient=\"0.4\"   ambient=\"0\"   diffuse=\"1\"   selectedSpecular=\"0.5\"   specular=\"0\"   power=\"1\"   opacity=\"1\"   visibility=\"false\"   clipping=\"false\"   backfaceCulling=\"true\"   scalarVisibility=\"false\"   vectorVisibility=\"false\"   tensorVisibility=\"false\"   autoScalarRange=\"true\"   scalarRange=\"0 100\"  ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelHierarchy\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelHierarchyNode3\"   name=\"vtkMRMLModelHierarchyNode3\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    parentNodeRef=\"vtkMRMLModelHierarchyNode1\"   modelNodeRef=\"vtkMRMLModelNode2\"   expanded=\"true\" ></ModelHierarchy>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <Model\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelNode2\"   name=\"GM-Surface\"   hideFromEditors=\"false\"   selectable=\"true\"    selected=\"false\"    userTags=\"\"   displayNodeRef=\"vtkMRMLModelDisplayNode3\" ></Model>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' <ModelDisplay\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '   id=\"vtkMRMLModelDisplayNode3\"   name=\"vtkMRMLModelDisplayNode3\"   hideFromEditors=\"true\"   selectable=\"true\"    selected=\"false\"    color=\"1 0.8 0.701961\"   selectedColor=\"1 0 0\"   selectedAmbient=\"0.4\"   ambient=\"0\"   diffuse=\"1\"   selectedSpecular=\"0.5\"   specular=\"0\"   power=\"1\"   opacity=\"1\"   visibility=\"true\"   clipping=\"false\"   backfaceCulling=\"true\"   scalarVisibility=\"false\"   vectorVisibility=\"false\"   tensorVisibility=\"false\"   autoScalarRange=\"true\"   scalarRange=\"0 100\"  ></ModelDisplay>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '</SceneSnapshot>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelStorage\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelStorageNode1\" name=\"vtkMRMLModelStorageNode1\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  fileName=\"'${WMSurface}'\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} '<ModelStorage\\n')"<<endl;
	CortThickPipelineBMS<<"AppendFile(${MRMLfile} ' id=\"vtkMRMLModelStorageNode2\" name=\"vtkMRMLModelStorageNode2\" hideFromEditors=\"true\" selectable=\"true\"  selected=\"false\"  fileName=\"'${GMSurface}'\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>\\n')"<<endl;
}

void Computation::CenterImage(string Input,string Output)
{
	//INPUT VOLUME READER
        VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
        imageReader->SetFileName(Input) ;
        imageReader->Update();

	ImageType::Pointer volume = imageReader->GetOutput( );


	//GET SIZE OF THE INPUT VOLUME
	ImageType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

	//SET THE VALUE THE NEW ORIGIN
	ImageType::PointType outputOrigin;
    	for(unsigned int i = 0; i<Dim; i++)
	{
		outputOrigin[i]  = volumeSize[i]/2;
		outputOrigin[i]  = -outputOrigin[i];
	}
	ImageType::Pointer OutputImage = ImageType::New();
	volume->SetOrigin(outputOrigin);

	//OUTPUT VOLUME WRITER
	VolumeWriterType::Pointer imageWriter = VolumeWriterType::New();
	imageWriter->SetFileName(Output);
	imageWriter->SetInput(volume);
	imageWriter->UseCompressionOn();
	try
	{
        	imageWriter->Update();
	}
	catch(itk::ExceptionObject Exp)
	{
		cerr<<"Error writer "<<Exp<<endl;
	}
}



void Computation::Compute()
{
	cout<<"\n\nComputing RegionalCortThickPipeline..."<<endl;
	if(!((m_Parameters->getlabel()).empty()) && ((m_Parameters->getrawImage()).empty()))
	{
		cerr<<"Error : a raw image is needed with the label image."<<endl;
		exit(0);
	}
	if(((m_Parameters->getlabel()).empty()) && ((m_Parameters->getsegAtlasDir()).empty()))
	{
		cerr<<"Error : a segmentation atlas directory is needed."<<endl;
		exit(0);
	}
	if(!((m_Parameters->getparcellation()).empty()) && ((m_Parameters->getatlas()).empty()))
	{
		cerr<<"Error : an space coordinate atlas is needed with the parcellation file."<<endl;
		exit(0);
	}
	setpipelineMode();
	WriteBMSFile();
	ExecuteBatchmakeScript();
	vector<BMString> labelFile = m_Parser.GetScriptActionManager()->GetVariable( "labelFile" );
	vector<BMString> labelPath = m_Parser.GetScriptActionManager()->GetVariable( "LabelPath" );
	vector<BMString> labelName = m_Parser.GetScriptActionManager()->GetVariable( "LabelName" );
	vector<BMString> outputFormat = m_Parser.GetScriptActionManager()->GetVariable( "outputFormat" );
	vector<BMString> refName = m_Parser.GetScriptActionManager()->GetVariable( "RefName" );
	m_labelFileCentered = (labelPath.front()).fromVariable().GetValue() + sep + (labelName.front()).fromVariable().GetValue() + "_centered." + (outputFormat.front()).fromVariable().GetValue();
	m_refName = (refName.front()).fromVariable().GetValue();
	CenterImage((labelFile.front()).fromVariable().GetValue(),m_labelFileCentered);
	WriteBMSFile2();
	ExecuteBatchmakeScript2();
	
	cout<<"Computing ARCTIC : Done!"<<endl;
}

