//#include "itkImageFileWriter.h"

//#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

//#include "CLIModuleTemplateCLP.h"

//add 0814
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkPasteImageFilter.h"

#include "itkImageDuplicator.h"
#include "itkImageToVTKImageFilter.h" 


// vtkITK includes
#include "vtkITKArchetypeImageSeriesScalarReader.h"

#include "vtkPoissonReconstruction.h"
#include "SurfaceReconstructionCLP.h"



#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkMergePoints.h>
#include <vtkPointSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkMath.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangleStrip.h>
#include <vtkMarchingCubes.h>
#include <vtkKdTreePointLocator.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkImageAccumulate.h>
#include <vtkMetaImageReader.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkImageThreshold.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkImageChangeInformation.h>
#include <vtkNew.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkStripper.h>
#include <vtkReverseSense.h>
#include <vtkCleanPolyData.h>
#include <vtkAppendFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkBooleanOperationPolyDataFilter.h>


#include <iostream>






// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

//typedef itk::Image< PixelType,  3> InputImageType;
//typedef itk::Image< unsigned char,  3> OutputImageType;

typedef itk::Image< unsigned char, 3> MidImageType;
typedef itk::Image< unsigned char, 2> SliceImageType;

vtkSmartPointer<vtkPoints>  getVolumeContourPoints (const MidImageType::Pointer midImage) 
{
	typedef itk::ExtractImageFilter< MidImageType,SliceImageType > 
		extractFilterType;
	typedef itk::BinaryContourImageFilter <SliceImageType,SliceImageType >
		binaryContourImageFilterType;
	typedef itk::CastImageFilter< SliceImageType, MidImageType > 
		CastFilterType;
	typedef itk::PasteImageFilter<MidImageType,MidImageType > 
		PasteFilterType;


		
	//Get Resion, size and start of input midImage	
	MidImageType::RegionType inputRegion =
		midImage->GetLargestPossibleRegion();
	MidImageType::SizeType inputSize =
		inputRegion.GetSize();
	MidImageType::IndexType start = inputRegion.GetIndex();
	MidImageType::RegionType desiredRegion;

	//New a output volume used to contain contour volume
	MidImageType::Pointer  outputContour = MidImageType::New();
	outputContour->SetRegions(inputRegion);
	outputContour->Allocate();
	outputContour = midImage;
	

	extractFilterType::Pointer extractFilter = extractFilterType::New();
	extractFilter->SetDirectionCollapseToIdentity();
		

	binaryContourImageFilterType::Pointer binaryContourFilter = 
		binaryContourImageFilterType::New ();
	CastFilterType::Pointer castFilter = CastFilterType::New();
	PasteFilterType::Pointer pasteFilter = PasteFilterType::New();

	unsigned int zSize=inputSize[2];
	inputSize[2]=0;
	for (unsigned int k=0; k<zSize; k++)
	{
		//Extract one slice from midImage
		start[2] = k;
		desiredRegion.SetSize(  inputSize  );
		desiredRegion.SetIndex( start );
		extractFilter->SetExtractionRegion( desiredRegion );
		extractFilter->SetInput(midImage );

		//Extract a contour of one slice
		binaryContourFilter->SetInput(extractFilter->GetOutput() );

		//Case one contour slice of 2D to a slice of 3D which can be a input of pasteFilter
		castFilter->SetInput(binaryContourFilter->GetOutput());
		castFilter->UpdateLargestPossibleRegion();

	
		//Paste a contour slice to outputContour volume
		pasteFilter->SetSourceImage(castFilter->GetOutput() );
		pasteFilter->SetDestinationImage( outputContour );
		pasteFilter->SetDestinationIndex( start );

		//sliceImage3D = castFilter->GetOutput();
		pasteFilter->SetSourceRegion( castFilter->GetOutput()->GetBufferedRegion() );
		pasteFilter->Update();
	
		outputContour=pasteFilter->GetOutput();

	}

	inputSize[2]=zSize;
	
	//debug
	//typedef itk::ImageFileWriter<MidImageType> WriterType3D;
	//WriterType3D::Pointer writeroutconture = WriterType3D::New();
	//writeroutconture->SetFileName("D:\\outputcontour.nrrd");
	//writeroutconture->SetInput(outputContour);
	//writeroutconture->Update();
	
	
	//Get points at boundary  	
	MidImageType::IndexType index;
	vtkSmartPointer<vtkPoints> pointsInit = vtkSmartPointer<vtkPoints>::New();
	MidImageType::PointType contourPoint;
  

	float pointXYZ[3];
	const float matrixVTK2Slicer [3][3]={{-1, 0, 0}, {0,-1, 0}, {0, 0,1}};
	float pointXYZ_S[3]; //slicer point
	int pointNum=0;

	for(unsigned int k=0; k<inputSize[2]; k++)
	{
		for(unsigned int j=0; j<inputSize[1]; j++)
		{
			for(unsigned int i=0; i<inputSize[0]; i++)
			{
				index[0]=i;
				index[1]=j;
				index[2]=k;
		
				if ( outputContour->GetPixel(index)!=0)
				{

					pointNum++;
					outputContour->TransformIndexToPhysicalPoint(index,contourPoint);
					pointXYZ[2]=contourPoint[2];
					pointXYZ[1]=contourPoint[1];
					pointXYZ[0]=contourPoint[0];
					vtkMath::Multiply3x3(matrixVTK2Slicer, pointXYZ,pointXYZ_S); 
					pointsInit->InsertNextPoint(pointXYZ_S); 


				}
				
			}
	
		}
	
	}
	

	return pointsInit;
}

std::vector<int> getLabels(vtkSmartPointer<vtkImageData> midImageVTK)
{
	vtkSmartPointer<vtkImageAccumulate> hist = vtkSmartPointer<vtkImageAccumulate>::New();
	hist->SetInput(midImageVTK);
	int extentMax = 0;
	double dImageScalarMax = midImageVTK->GetScalarTypeMax();
	extentMax = (int)(floor(dImageScalarMax - 1.0));

	int biggestBin = 1000000;     // VTK_INT_MAX - 1;
	if (extentMax < 0 || extentMax > biggestBin)
	{
		std::cout << "\nWARNING: due to lack of color label information and an image with a scalar maximum of "
					<< dImageScalarMax << ", using  " << biggestBin << " as the histogram number of bins" << endl;
		extentMax = biggestBin;
	}
	else
	{
		std::cout
		<<
		"\nWARNING: due to lack of color label information, using the full scalar range of the input image when calculating the histogram over the image: "
		<< extentMax << endl;
	}
	
	hist->SetComponentExtent(0, extentMax, 0, 0, 0, 0);
	hist->SetComponentOrigin(0, 0, 0);
	hist->SetComponentSpacing(1, 1, 1);
	hist->Update();
		
	double *max = hist->GetMax();
	double *min = hist->GetMin();
	if (min[0] == 0 || min[0] <0)
		min[0]=1;

	int StartLabel = (int)floor(min[0]);
	int EndLabel = (int)floor(max[0]);

	int numModelsToGenerate = 0;
	for(int i = StartLabel; i <= EndLabel; i++)
	{
		if ((int)floor((((hist->GetOutput())->GetPointData())->GetScalars())->GetTuple1(i)) > 0)
		{
			numModelsToGenerate++;
		}
	}

	std::cout<<
	"\nNumber of models to be generated "<< numModelsToGenerate << endl;

	// ModelMakerMarch
	double      labelFrequency = 0.0;;
	std::vector<int> loopLabels;
	std::vector<int> skippedModels;
	std::vector<int> madeModels;
	// set up the loop list with all the labels between start and end
	for(int i = StartLabel; i <= EndLabel; i++)
	{
		loopLabels.push_back(i);
	}
    
	for(::size_t l = 0; l < loopLabels.size(); l++)
	{
		// get the label out of the vector
		int i = loopLabels[l];
			
		labelFrequency = (((hist->GetOutput())->GetPointData())->GetScalars())->GetTuple1(i);
		if (labelFrequency == 0.0)
		{
			skippedModels.push_back(i);
			continue;
		}
			else
		{
			madeModels.push_back(i);
		}
	}

	return madeModels;


}

vtkSmartPointer<vtkPolyData>  getMarchingcubePoly (vtkSmartPointer<vtkImageData> midImageVTK, int labelValue, vtkSmartPointer<vtkTransform> transformIJKtoRAS ) 
{
	vtkSmartPointer <vtkImageThreshold> imageThreshold = vtkSmartPointer<vtkImageThreshold>::New();
	imageThreshold->SetInput(midImageVTK);
	imageThreshold->SetReplaceIn(1);
	imageThreshold->SetReplaceOut(1);
	imageThreshold->SetInValue(200);
	imageThreshold->SetOutValue(0);

	imageThreshold->ThresholdBetween(labelValue, labelValue);
	(imageThreshold->GetOutput())->ReleaseDataFlagOn();
	imageThreshold->ReleaseDataFlagOn();
	vtkSmartPointer <vtkImageToStructuredPoints> imageToStructuredPoints = vtkSmartPointer<vtkImageToStructuredPoints>::New();
	imageToStructuredPoints->SetInput(imageThreshold->GetOutput());
	try
	{
		imageToStructuredPoints->Update();
	}
	catch(...)
	{
		std::cerr << "ERROR while updating image to structured points for label " << labelValue << std::endl;
		return NULL;
	}
	imageToStructuredPoints->ReleaseDataFlagOn();
      

	vtkSmartPointer <vtkMarchingCubes> mcubes = vtkSmartPointer<vtkMarchingCubes>::New();

	mcubes->SetInput(imageToStructuredPoints->GetOutput());
	mcubes->SetValue(0, 100.5);
	mcubes->ComputeScalarsOff();
	mcubes->ComputeGradientsOff();
	mcubes->ComputeNormalsOff();
	(mcubes->GetOutput())->ReleaseDataFlagOn();
	try
	{
		mcubes->Update();
	}
	catch(...)
	{
		std::cerr << "ERROR while running marching cubes, for label " << labelValue << std::endl;
		return NULL;
	}


	if ((mcubes->GetOutput())->GetNumberOfPolys()  == 0)
	{
		std::cout << "Cannot create a model from label " << labelValue
		<< "\nNo polygons can be created,\nthere may be no voxels with this label in the volume." << endl;

		if (transformIJKtoRAS)
		{
			transformIJKtoRAS = NULL;
		}
		
		if (imageThreshold)
		{
      
			imageThreshold->SetInput(NULL);
			imageThreshold->RemoveAllInputs();
			imageThreshold = NULL;

		}
		if (imageToStructuredPoints)
		{
			imageToStructuredPoints->SetInput(NULL);
			imageToStructuredPoints = NULL;
		}
		if (mcubes)
		{
			mcubes->SetInput(NULL);
			mcubes = NULL;
		}

		std::cout << "...continuing" << endl;

		//return EXIT_FAILURE;
		return NULL;
	}
	  

	vtkSmartPointer <vtkDecimatePro> decimator = vtkSmartPointer<vtkDecimatePro>::New();
	decimator->SetInput(mcubes->GetOutput());
	decimator->SetFeatureAngle(60);
	decimator->SplittingOff();
	decimator->PreserveTopologyOn();

	decimator->SetMaximumError(1);
	decimator->SetTargetReduction(0.01);
	(decimator->GetOutput())->ReleaseDataFlagOff();

	try
	{
	decimator->Update();
	}
	catch(...)
	{
	std::cerr << "ERROR decimating model " << labelValue << std::endl;
	return NULL;
	}
	 
	 
	vtkSmartPointer<vtkReverseSense> reverser = vtkSmartPointer<vtkReverseSense>::New();
	 
	if ((transformIJKtoRAS->GetMatrix())->Determinant() < 0)
	{
	reverser->SetInput(decimator->GetOutput());
	reverser->ReverseNormalsOn();
	(reverser->GetOutput())->ReleaseDataFlagOn();
	}  
	 
	 
	//Smooth
	vtkSmartPointer <vtkSmoothPolyDataFilter> smootherPoly = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smootherPoly->SetRelaxationFactor(0.33);
	smootherPoly->SetFeatureAngle(60);
	smootherPoly->SetConvergence(0);
	 
	if ((transformIJKtoRAS->GetMatrix())->Determinant() < 0)
	{
		smootherPoly->SetInput(reverser->GetOutput());
	}
	else
	{
		smootherPoly->SetInput(decimator->GetOutput());
	}
	 
	 

	smootherPoly->SetNumberOfIterations(98);
	smootherPoly->FeatureEdgeSmoothingOff();
	smootherPoly->BoundarySmoothingOff();
	(smootherPoly->GetOutput())->ReleaseDataFlagOn();

	try
	{
		smootherPoly->Update();
	}
	catch(...)
	{
		std::cerr << "ERROR updating Poly smoother for model " << labelValue << std::endl;
		return NULL;
	}
      
	vtkSmartPointer<vtkTransformPolyDataFilter> transformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformer->SetInput(smootherPoly->GetOutput());
	transformer->SetTransform(transformIJKtoRAS);
	(transformer->GetOutput())->ReleaseDataFlagOn();
		
	vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->ComputePointNormalsOn();
	normals->SetInput(transformer->GetOutput());
	normals->SetFeatureAngle(60);
	normals->SetSplitting(true);
	(normals->GetOutput())->ReleaseDataFlagOn();

	vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
	stripper->SetInput(normals->GetOutput());
	(stripper->GetOutput())->ReleaseDataFlagOff();
	  
	   
	// the poly data output from the stripper can be set as an input to a
	// model's polydata
	try
	{
		(stripper->GetOutput())->Update();
	}
	catch(...)
	{
		std::cerr << "ERROR updating stripper for model " << labelValue << std::endl;
		return NULL;
	}

	return stripper->GetOutput();

}

MidImageType::Pointer getMask (const MidImageType::Pointer midImage, int labelValue)
{

	//Set un-zero pixel to 255, this is a limit of binaryContourImageFilterType
	MidImageType::SizeType midImageSize =
	midImage->GetLargestPossibleRegion().GetSize();
	MidImageType::IndexType indexMid;
	unsigned char pixelvalue;
			
	typedef itk::ImageDuplicator< MidImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(midImage);
	duplicator->Update();
	MidImageType::Pointer clonedImage = duplicator->GetOutput();
 


	for(unsigned int n=0; n<midImageSize[2]; n++)
	{
		for(unsigned int m=0; m<midImageSize[1]; m++)
		{
			for(unsigned int l=0; l<midImageSize[0]; l++)
			{
				indexMid[0]=l;
				indexMid[1]=m;
				indexMid[2]=n;
				pixelvalue=clonedImage->GetPixel(indexMid);

				if (pixelvalue == unsigned char (labelValue))
					clonedImage->SetPixel(indexMid, 255);
				else
					clonedImage->SetPixel(indexMid, 0);
			}
		}
	}
	
    return clonedImage;

}

double max (double a, double b)
{
	if (a>=b)
		return a;
	else
		return b;
}

double min (double a, double b)
{
	if (a<=b)
		return a;
	else
		return b;
}
void combinePoints (vtkSmartPointer<vtkPolyData> pointsetA, vtkSmartPointer<vtkPolyData> pointsetB)

{

	vtkIdType id;
	vtkSmartPointer<vtkMergePoints> mergePoints = 
	vtkSmartPointer<vtkMergePoints>::New();
	mergePoints->SetDataSet(pointsetA);
	mergePoints->SetDivisions(10,10,10);
	
	double boundsA[6];
	double boundsB[6];
	
	pointsetA->GetBounds(boundsA);
	pointsetB->GetBounds(boundsB);
	double bounds[6];
	bounds[0]=min(boundsA[0], boundsB[0]);
	bounds[2]=min(boundsA[2], boundsB[2]);
	bounds[4]=min(boundsA[4], boundsB[4]);
	bounds[1]=max(boundsA[1], boundsB[1]);
	bounds[3]=max(boundsA[3], boundsB[3]);
	bounds[5]=max(boundsA[5], boundsB[5]);
	
	mergePoints->InitPointInsertion(pointsetA->GetPoints(), bounds);
    
	for (vtkIdType i = 0; i < pointsetA->GetNumberOfPoints(); i++)
	{
	mergePoints->InsertUniquePoint(pointsetA->GetPoint(i), id);
	}

	for (vtkIdType i = 0; i < pointsetB->GetNumberOfPoints(); i++)
	{
	mergePoints->InsertUniquePoint(pointsetB->GetPoint(i), id);

	}

	std::cout << "There are now "
		<< pointsetA->GetNumberOfPoints() << " points. inside function" << std::endl;



	
	mergePoints = NULL;

}
vtkSmartPointer<vtkPolyData> combinePolys (vtkSmartPointer<vtkPolyData> polyA, vtkSmartPointer<vtkPolyData> polyB)
{
	
	vtkSmartPointer<vtkPolyData> polyMarchingCombine = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkAppendFilter> appendFilterAS =   vtkSmartPointer<vtkAppendFilter>::New();
	vtkSmartPointer<vtkAppendFilter> appendFilterASC =   vtkSmartPointer<vtkAppendFilter>::New();


		appendFilterAS->AddInputConnection(polyA->GetProducerPort());
		appendFilterAS->AddInputConnection(polyB->GetProducerPort());
		appendFilterAS->Update();


			polyMarchingCombine->ShallowCopy(appendFilterAS->GetOutput()); 
	
	return polyMarchingCombine;
}


template <class T>
int DoIt( int argc, char * argv[], T )
{
	PARSE_ARGS;

	typedef    T InputPixelType;
	typedef    T OutputPixelType;

	typedef itk::Image<InputPixelType,  3> InputImageType;
	typedef itk::Image<OutputPixelType, 3> OutputImageType;


	typedef itk::ImageFileReader<InputImageType>  ReaderType;
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	typedef itk::CastImageFilter< InputImageType, MidImageType > CastFilterType;
	typedef itk::ImageToVTKImageFilter<MidImageType> ConnectorType;
	
	//
	vtkSmartPointer<vtkPoints> pointsInit = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New(); 
	vtkSmartPointer<vtkPolyData> polyDataCombine = vtkSmartPointer<vtkPolyData>::New();  //contain points on contours
	vtkSmartPointer<vtkPolyData> marchingPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> marchingPolyPre = vtkSmartPointer<vtkPolyData>::New();

	//

	vtkSmartPointer<vtkKdTreePointLocator> kDTree;
	vtkSmartPointer<vtkPolyData> polyMarchingCombine; //Marching cube result
	vtkSmartPointer<vtkPolyData> surfacePre = vtkSmartPointer<vtkPolyData>::New();

	//
	std::string Name="Surface";
	std::vector<int> madeModels; //The labels for surface reconstruction
	std::vector< MidImageType::Pointer > midImage_V;
	std::vector<std::string> inputFileName;
	std::vector< vtkSmartPointer<vtkImageData> > midImageVTK_V;
	std::vector< vtkSmartPointer<vtkTransform> > transformIJKtoRAS_V;

	
	
	//
	if (InputVolumeAxial.size())
	{
     inputFileName.push_back(InputVolumeAxial.c_str());
    }
	if (InputVolumeSag.size())
	{
    inputFileName.push_back(InputVolumeSag.c_str());
    }
	if (InputVolumeCor.size())
	{
    inputFileName.push_back(InputVolumeCor.c_str());
    }
	

	//Read input volumes and get labels
	for(::size_t l=0; l <inputFileName.size(); l++)
	{
 
		std::string filename=inputFileName[l];
		std::cout << "Now File name is  "
				<< inputFileName[l] << std::endl;
	
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(inputFileName[l].c_str());
		reader->Update();

		//typename WriterType::Pointer writer = WriterType::New();
		//writer->SetFileName(OutputVolumeAxial.c_str());
		//writer->SetInput(reader->GetOutput());
		//writer->Update();


		//Cast input volume to a volume with unsigned char pixel
		typename CastFilterType::Pointer castInputFilterS = CastFilterType::New();
		castInputFilterS->SetInput(reader->GetOutput());
		castInputFilterS->Update();
		
		MidImageType::Pointer midImage= MidImageType ::New();
		midImage= castInputFilterS->GetOutput();
		midImage_V.push_back(midImage);
	
		//itk to vtk
	
		typename ConnectorType::Pointer connector =ConnectorType::New();
		connector->SetInput(midImage);
	
		vtkNew<vtkImageChangeInformation> ici;
		ici->SetInput(connector->GetOutput());
		ici->SetOutputSpacing(1, 1, 1);
		ici->SetOutputOrigin(0, 0, 0);
		ici->Update();

		vtkSmartPointer<vtkImageData> midImageVTK = vtkSmartPointer<vtkImageData>::New();
		
		midImageVTK = ici->GetOutput();
		midImageVTK->Update();

		midImageVTK_V.push_back(midImageVTK);

		//for acquiring transform
		vtkSmartPointer<vtkITKArchetypeImageSeriesScalarReader> scalarReader = 
			vtkSmartPointer<vtkITKArchetypeImageSeriesScalarReader>::New();
		vtkSmartPointer<vtkTransform> transformIJKtoRAS = vtkSmartPointer<vtkTransform>::New();
		
		scalarReader->SetArchetype(filename.c_str());
		scalarReader->SetOutputScalarTypeToNative();
		scalarReader->SetDesiredCoordinateOrientationToNative();
		scalarReader->SetUseNativeOriginOn();
		scalarReader->Update();
		transformIJKtoRAS->SetMatrix(scalarReader->GetRasToIjkMatrix());
		transformIJKtoRAS->Inverse();
		
		transformIJKtoRAS_V.push_back(transformIJKtoRAS);


		//Generate a histogram of the labels
		if (l==0)
			madeModels = getLabels(midImageVTK);
		else
			continue;
	
	}
	
	

	
	for (::size_t m =0; m < madeModels.size(); m++)
	{
		int labelValue= madeModels[m];
		std::stringstream lable;
		lable<<labelValue;

		for(::size_t l=0; l <inputFileName.size(); l++)
		{

			pointsInit=getVolumeContourPoints(getMask(midImage_V[l], labelValue));
			polydata->SetPoints(pointsInit); 
			std::cout << "There are  "
				<< polydata->GetNumberOfPoints() << " points in " << inputFileName[l] << std::endl;
	
			//Get marching cube result
			marchingPoly=getMarchingcubePoly(midImageVTK_V[l],labelValue, transformIJKtoRAS_V[l]);

			std::cout << "There are  "
				<< marchingPoly->GetNumberOfPoints() << " marchingpoy in " << inputFileName[l] << std::endl;

			if(l==0)
			{
				polyDataCombine->ShallowCopy(polydata);
				
				if(polyMarchingCombine)
				{
					polyMarchingCombine=NULL;
				}
				polyMarchingCombine = vtkSmartPointer<vtkPolyData>::New();
				polyMarchingCombine->ShallowCopy(marchingPoly); 
			}
			else
			{
				//Combine  points
				std::cout << "There are  "
				<< polyDataCombine->GetNumberOfPoints() << "input points in polydata before" << " l= " <<l<< std::endl;
				
				combinePoints (polyDataCombine, polydata);

				std::cout << "There are  "
				<< polyDataCombine->GetNumberOfPoints() << "input points in polydata after" << " l= " <<l<< std::endl;


				//combine marching cube results
				std::cout << "There are  "
				<< marchingPoly->GetNumberOfPoints() << "input marching result before combine."  << " l= " <<l<< std::endl;
				
				if(polyMarchingCombine)
				{
					polyMarchingCombine=NULL;
				}
				polyMarchingCombine = vtkSmartPointer<vtkPolyData>::New();
				polyMarchingCombine->ShallowCopy(combinePolys(marchingPolyPre, marchingPoly)); 
				
				std::cout << "There are  "
				<< polyMarchingCombine->GetNumberOfPoints() << " input marching result after combine." << std::endl;

			
			}

			marchingPolyPre->ShallowCopy(polyMarchingCombine);
		}
		
		if(kDTree)
		{
			kDTree->SetDataSet(NULL);		
			kDTree=NULL;
		}
		kDTree=	vtkSmartPointer<vtkKdTreePointLocator>::New();
		kDTree->SetDataSet(polyMarchingCombine);
		kDTree->BuildLocator();

		vtkSmartPointer<vtkFloatArray> pointNormalsRetrieved = 
		vtkFloatArray::SafeDownCast(kDTree->GetDataSet()->GetPointData()->GetNormals());
 
		double pTTrans[3];
		vtkIdType iD;
		double closestPoint[3];
		double pNdouble[3];
		float pN[3];
		double dotRslt;
		double normpNPre;
		double normpN;
		ofstream normalOut("normalsExtract.txt");

		vtkSmartPointer<vtkFloatArray> pointNormalsArray = vtkSmartPointer<vtkFloatArray>::New();
		pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
		pointNormalsArray->SetNumberOfTuples(polyMarchingCombine->GetNumberOfPoints());

  
		for (unsigned int i=0; i< polyDataCombine->GetNumberOfPoints(); i++)
		{
		polyDataCombine->GetPoint(i, pTTrans);
  
		// Find the closest points to TestPoint
		iD = kDTree->FindClosestPoint(pTTrans);
		normalOut << "The closest point is point " << iD << std::endl;
 
		//Get the coordinates of the closest point
 
		kDTree->GetDataSet()->GetPoint(iD, closestPoint);
		normalOut << "Coordinates:    " << closestPoint[0] << " " << closestPoint[1] << " " << closestPoint[2] << std::endl;
		normalOut <<"points on contour" <<pTTrans[0] << " " <<pTTrans[1]<<" " << pTTrans[2] << std::endl;


		pointNormalsRetrieved->GetTuple(iD, pNdouble);
		normalOut << "Point normal " << iD << ": "  << pNdouble[0] << " " << pNdouble[1] << " " << pNdouble[2] << std::endl;

		pN[0]=(float)pNdouble[0];
		pN[1]=(float)pNdouble[1];
		pN[2]=(float)pNdouble[2];

		pointNormalsArray->SetTuple(i, pN);  

		}

		// Add the normals to the points in the polydata
		polyDataCombine->GetPointData()->SetNormals(pointNormalsArray); 


		//Make a vtkPolyData with a vertex on each point.
		vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
		vertexFilter->SetInputConnection(polyDataCombine->GetProducerPort());
		vertexFilter->Update();


		vtkSmartPointer<vtkPolyData> polydataNew = vtkSmartPointer<vtkPolyData>::New();
		polydataNew->ShallowCopy(vertexFilter->GetOutput()); 
  


		vtkSmartPointer<vtkXMLPolyDataWriter> writerPoly = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writerPoly->SetInput(polydataNew);
		writerPoly->SetFileName(outputPolyFile.c_str());
		writerPoly->Update();


		vtkSmartPointer<vtkXMLPolyDataReader> readerPoly = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		readerPoly->SetFileName(outputPolyFile.c_str());
		readerPoly->Update();

		//PoissonReconstruction
		vtkSmartPointer<vtkPoissonReconstruction> poissonFilter = 
		vtkSmartPointer<vtkPoissonReconstruction>::New();
		poissonFilter->SetDepth(12);
		poissonFilter->SetInputConnection(readerPoly->GetOutputPort());
		poissonFilter->Update();
  
  
  
		////Write the file
		//std::string fileName;
		//fileName = "f:" + std::string("/") + Name + lable.str()+ std::string(".vtp");
		//vtkSmartPointer<vtkXMLPolyDataWriter> writerSurface =
		//vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		//writerSurface->SetInputConnection(poissonFilter->GetOutputPort());
		//writerSurface->SetFileName(fileName.c_str());
		//writerSurface->Update();

		
		if (m==0)
		{
		surfacePre->ShallowCopy(poissonFilter->GetOutput());
		}
		else
		{

		vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
		vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();	
		booleanOperation->SetInputConnection( 0, surfacePre->GetProducerPort() );
		booleanOperation->SetInputConnection( 1, poissonFilter->GetOutputPort() );
		booleanOperation->SetOperationToUnion();
		booleanOperation->Update();
		surfacePre->ShallowCopy(booleanOperation->GetOutput());

		}
		
	}
	
	
	vtkSmartPointer<vtkXMLPolyDataWriter> writerSurface =
	vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writerSurface->SetInput(surfacePre);
	writerSurface->SetFileName(outputSurfaceFile.c_str());
	writerSurface->Update();
	

	return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(InputVolumeAxial, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
