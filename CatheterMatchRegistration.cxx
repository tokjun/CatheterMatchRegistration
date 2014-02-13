/*=========================================================================

  Program:   Registration stand-alone
  Module:    $HeadURL$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
=========================================================================*/

#include "CatheterMatchRegistrationCLP.h"

#include "itkOrientImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkMaskImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageMaskSpatialObject.h"

#include "itkQuaternionRigidTransformGradientDescentOptimizer.h"
#include "itkCatheterMatchRigidTransformGradientDescentOptimizer.h"
#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkQuaternionRigidTransform.h"
#include "itkCatheterMatchRigidTransform.h"
#include "itkResampleImageFilter.h"
#include "itkBinomialBlurImageFilter.h"

#include "itkPluginUtilities.h"

#include "itkTimeProbesCollectorBase.h"

#include "vnl/vnl_cross.h"

#define TESTMODE_ERROR_TOLERANCE 0.1

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

class ScheduleCommand : public itk::Command
{
public:
  typedef ScheduleCommand         Self;
  typedef itk::SmartPointer<Self> Pointer;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScheduleCommand, itk::Command);

  void SetLearningRates( std::vector<double> & LearningRates )
  {
    m_LearningRates = LearningRates;
  }

  void SetNumberOfIterations( std::vector<int> & NumberOfIterations )
  {
    m_NumberOfIterations = NumberOfIterations;
    this->m_NextChange = NumberOfIterations[0];
  }

  void SetSchedule( std::vector<int> & NumberOfIterations,
                    std::vector<double> & LearningRates )
  {
    this->SetNumberOfIterations(NumberOfIterations);
    this->SetLearningRates(LearningRates);
  }

  void SetRegistration(itk::ProcessObject* reg)
  {
    m_Registration = reg;
  }

  void DoExecute( itk::GradientDescentOptimizer* optimizer )
  {
    if( m_Schedule < m_NumberOfIterations.size() - 1 )
      {
      if( static_cast<int>(optimizer->GetCurrentIteration() )
          >= this->m_NumberOfIterations[m_Schedule] )
        {
        m_Schedule++;
        optimizer->SetLearningRate( this->m_LearningRates[m_Schedule] );
        this->m_NextChange = optimizer->GetCurrentIteration()
          + this->m_NumberOfIterations[m_Schedule];
        // std::cout << "Iteration: " << optimizer->GetCurrentIteration()
        //          << " Parameters: " << optimizer->GetCurrentPosition()
        //          << " LearningRate: " << optimizer->GetLearningRate()
        //          << std::endl;
//         std::cout << "<filter-comment>"
//                   << "Iteration: "
//                   << optimizer->GetCurrentIteration() << ", "
//                   << " Switching LearningRate: "
//                   << optimizer->GetLearningRate() << " "
//                   << "</filter-comment>"
//                   << std::endl;
        }
      }
  }

  void Execute( itk::Object *caller, const itk::EventObject & event )
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute( const itk::Object *caller, const itk::EventObject & itkNotUsed(event) )
  {
    itk::GradientDescentOptimizer* optimizer = (itk::GradientDescentOptimizer *)(const_cast<itk::Object *>(caller) );

    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetCurrentPosition() << "   ";
    std::cout << optimizer->GetValue() << std::endl;
    if( m_Registration )
      {
      // for our purposes, an iteration even is a progress event
      // m_Registration->UpdateProgress( static_cast<double>(optimizer->GetCurrentIteration()) /
      //                                static_cast<double>(optimizer->GetNumberOfIterations()));
      }

    if( optimizer->GetCurrentIteration() >= this->m_NextChange )
      {
      this->DoExecute( optimizer );
      }
  }

protected:
  std::vector<int>            m_NumberOfIterations;
  std::vector<double>         m_LearningRates;
  unsigned int                m_Schedule;
  unsigned int                m_NextChange;
  itk::ProcessObject::Pointer m_Registration;
  ScheduleCommand()
  {
    m_Schedule = 0;
    m_Registration = 0;
  }

  ~ScheduleCommand()
  {
  }
};

// typedef itk::Image<signed short, 3> Volume;

template <class T1, class T2>
int DoIt2( int argc, char * argv[], const T1 &, const T2 & )
{
  //
  // Command line processing
  //
  PARSE_ARGS;

  const    unsigned int ImageDimension = 3;
  typedef  T1                                        FixedPixelType; // ##
  typedef itk::Image<FixedPixelType, ImageDimension> FixedImageType; // ##
  typedef itk::Image< unsigned char, ImageDimension >  MaskImageType;

  typedef itk::ImageFileReader<FixedImageType>                   FixedFileReaderType;   // ##
  typedef itk::OrientImageFilter<FixedImageType, FixedImageType> FixedOrientFilterType; // ##

  typedef  T2                                         MovingPixelType; // ##
  typedef itk::Image<MovingPixelType, ImageDimension> MovingImageType; // ##

  typedef itk::ImageFileReader<MovingImageType>                    MovingFileReaderType;   // ##
  typedef itk::OrientImageFilter<MovingImageType, MovingImageType> MovingOrientFilterType; // ##

  typedef itk::MattesMutualInformationImageToImageMetric<FixedImageType, MovingImageType> MetricType;    // ##
  typedef itk::QuaternionRigidTransformGradientDescentOptimizer                           OptimizerType;
  typedef itk::CatheterMatchRigidTransformGradientDescentOptimizer                        CathMatchOptimizerType;
  
  typedef itk::LinearInterpolateImageFunction<MovingImageType, double>                    InterpolatorType; // ##
  typedef itk::ImageRegistrationMethod<FixedImageType, MovingImageType>                   RegistrationType; // ##
  typedef itk::QuaternionRigidTransform<double>                                           TransformType;
  typedef itk::CatheterMatchRigidTransform<double>                                        CathMatchTransformType;

  typedef OptimizerType::ScalesType                                                       OptimizerScalesType;
  typedef CathMatchOptimizerType::ScalesType                                              CathMatchOptimizerScalesType;

  typedef itk::ResampleImageFilter<MovingImageType, MovingImageType>                      ResampleType;             // ##
  typedef itk::LinearInterpolateImageFunction<MovingImageType, double>                    ResampleInterpolatorType; // ##
  typedef itk::ImageFileWriter<MovingImageType>                                           WriterType;               // ##
  typedef itk::ImageFileWriter<FixedImageType>                                            FixedWriterType;          // ##
  typedef itk::ContinuousIndex<double, 3>                                                 ContinuousIndexType;

  typedef itk::ThresholdImageFilter<FixedImageType>  FixedThresholdFilterType;
  typedef itk::ThresholdImageFilter<MovingImageType>  MovingThresholdFilterType;
  typedef itk::ResampleImageFilter<FixedImageType, MaskImageType>     FixedResampleType;
  typedef itk::ResampleImageFilter<MovingImageType, MaskImageType>    MovingResampleType;
  typedef itk::NearestNeighborInterpolateImageFunction<FixedImageType> FixedNNInterpolator;
  typedef itk::NearestNeighborInterpolateImageFunction<MovingImageType> MovingNNInterpolator;

  typedef itk::ImageMaskSpatialObject< ImageDimension >   MaskType;

  // bool DoInitializeTransform = false;
  // int RandomSeed = 1234567;

  // Add a time probe
  itk::TimeProbesCollectorBase collector;

  // Read the fixed and moving volumes
  //
  //
  typename FixedFileReaderType::Pointer fixedReader = FixedFileReaderType::New();
  fixedReader->SetFileName( FixedImageFileName.c_str() );

  try
    {
    collector.Start( "Read fixed volume" );
    fixedReader->Update();
    collector.Stop( "Read fixed volume" );
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Error Reading Fixed image: " << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  typename MovingFileReaderType::Pointer movingReader = MovingFileReaderType::New();
  movingReader->SetFileName( MovingImageFileName.c_str() );

  try
    {
    collector.Start( "Read moving volume" );
    movingReader->Update();
    collector.Stop( "Read moving volume" );
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Error Reading Moving image: " << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // CathMatch option (enabled only if both fixed and moving
  // catheter transforms are specified
  bool CathMatch = false;

  // Check if any catheter transforms are specified
  typedef itk::TransformFileReader TransformReaderType;
  TransformReaderType::TransformType::Pointer fixedCatheterTransform;
  TransformReaderType::TransformType::Pointer movingCatheterTransform;
  TransformReaderType::Pointer fixedCatheterTransformReader;
  TransformReaderType::Pointer movingCatheterTransformReader;
  
  if( FixedCatheterTransform != "" && MovingCatheterTransform != "")
    {
    fixedCatheterTransformReader = TransformReaderType::New();
    movingCatheterTransformReader = TransformReaderType::New();
    fixedCatheterTransformReader->SetFileName( FixedCatheterTransform );
    movingCatheterTransformReader->SetFileName( MovingCatheterTransform );

    try
      {
      fixedCatheterTransformReader->Update();
      movingCatheterTransformReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }

    if (fixedCatheterTransformReader->GetTransformList()->size() != 0 &&
        movingCatheterTransformReader->GetTransformList()->size() != 0)
      {
      fixedCatheterTransform = *(fixedCatheterTransformReader->GetTransformList()->begin() );
      movingCatheterTransform = *(movingCatheterTransformReader->GetTransformList()->begin() );
      CathMatch = true;
      }
    }


  // Read the fixed and moving mask images
  typename FixedFileReaderType::Pointer fixedMaskReader;
  typename MovingFileReaderType::Pointer movingMaskReader;

  // Flag to use ROI. Enabled when fixed and moving mask images are specified.
  bool fMask = false;

  if ( FixedMaskImageFileName != "" && MovingMaskImageFileName != "")
    {
    fMask = true;

    fixedMaskReader = FixedFileReaderType::New();
    fixedMaskReader->SetFileName( FixedMaskImageFileName.c_str() );
    
    try
      {
      collector.Start( "Read fixed mask volume" );
      fixedMaskReader->Update();
      collector.Stop( "Read fixed mask volume" );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Error Reading Fixed mask image: " << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    
    movingMaskReader = MovingFileReaderType::New();
    movingMaskReader->SetFileName( MovingMaskImageFileName.c_str() );
    
    try
      {
      collector.Start( "Read moving mask volume" );
      movingMaskReader->Update();
      collector.Stop( "Read moving mask volume" );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Error Reading Moving image: " << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    }


  /*
  typename FixedWriterType::Pointer fixedWriter = FixedWriterType::New();
  fixedWriter->SetFileName ( "/tmp/fixed.nrrd" );
  fixedWriter->SetInput ( fixedReader->GetOutput() );
  fixedWriter->Write();

  typename WriterType::Pointer movingWriter = WriterType::New();
  movingWriter->SetFileName ( "/tmp/moving.nrrd" );
  movingWriter->SetInput ( movingReader->GetOutput() );
  movingWriter->Write();
  */

  // This was added by Fedorov:
  // In testing mode, the initial transform is assumed to be derived by
  // registering the provided input with the same parameters on *some*
  // platform. The testing mode is essentially testing *reproducibility*, not
  // correctness of the registration.
  //
  // The measure of reproducibility used is the error vector for a (0,0,0)
  // image index subject to the input transform vs. the recovered transform.
  // This measure is caclulated and reported after registration.
  //
  // NOTE: Testing mode can only be invoked from Command line interface: the
  // "TestingMode" parameter is "hidden" in GUI.
  typename TransformType::Pointer groundTruthTransform = NULL;
  typename MovingImageType::Pointer movingImageTestingMode = NULL;
  if( TestingMode )
    {
    typedef itk::TransformFileReader TransformReaderType;
    TransformReaderType::Pointer groundTruthTransformReader;
    groundTruthTransform = TransformType::New();

    if( InitialTransform != "" )
      {
      groundTruthTransformReader = TransformReaderType::New();
      groundTruthTransformReader->SetFileName( InitialTransform );
      try
        {
        groundTruthTransformReader->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
        }
      if( groundTruthTransformReader->GetTransformList()->size() == 0 )
        {
        std::cerr << "Non-empty transform should be specified in Testing Mode" << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
      {
      std::cerr << "Transform must be specified in TestingMode!" << std::endl;
      return EXIT_FAILURE;
      }

    TransformReaderType::TransformType::Pointer initial
      = *(groundTruthTransformReader->GetTransformList()->begin() );

    // most likely, the transform coming in is a subclass of
    // MatrixOffsetTransformBase
    typedef itk::MatrixOffsetTransformBase<double, 3, 3> DoubleMatrixOffsetType;
    typedef itk::MatrixOffsetTransformBase<float, 3, 3>  FloatMatrixOffsetType;

    DoubleMatrixOffsetType::Pointer da
      = dynamic_cast<DoubleMatrixOffsetType *>(initial.GetPointer() );
    FloatMatrixOffsetType::Pointer fa
      = dynamic_cast<FloatMatrixOffsetType *>(initial.GetPointer() );

    if( da )
      {
      vnl_svd<double> svd(da->GetMatrix().GetVnlMatrix() );

      groundTruthTransform->SetMatrix( svd.U() * vnl_transpose(svd.V() ) );
      groundTruthTransform->SetOffset( da->GetOffset() );
      }
    else if( fa )
      {
      vnl_matrix<double> t(3, 3);
      for( int i = 0; i < 3; ++i )
        {
        for( int j = 0; j < 3; ++j )
          {
          t.put(i, j, fa->GetMatrix().GetVnlMatrix().get(i, j) );
          }
        }

      vnl_svd<double> svd( t );

      groundTruthTransform->SetMatrix( svd.U() * vnl_transpose(svd.V() ) );
      groundTruthTransform->SetOffset( fa->GetOffset() );
      }
    else
      {
      std::cout << "Initial transform is an unsupported type.\n";
      return EXIT_FAILURE;
      }

    std::cout << "Testing mode ground truth transform: ";
    groundTruthTransform->Print( std::cout );
    }


  //// Masking
  // Threshold -- have to threshold the mask volume first
  //typename FixedImageType::Pointer fixedMaskedImage;
  //typename MovingImageType::Pointer movingMaskedImage;
  MaskType::Pointer  fixedSpatialObjectMask = MaskType::New();
  MaskType::Pointer  movingSpatialObjectMask = MaskType::New();

  if ( fMask )
    {
    typename FixedThresholdFilterType::Pointer fixedThresholdFilter = FixedThresholdFilterType::New();
    typename MovingThresholdFilterType::Pointer movingThresholdFilter = MovingThresholdFilterType::New();

    // Threshold to remove labels other than 1
    fixedThresholdFilter->SetInput(0, fixedMaskReader->GetOutput() );
    fixedThresholdFilter->SetOutsideValue(0);
    fixedThresholdFilter->ThresholdOutside(1, 255);
    fixedThresholdFilter->ReleaseDataFlagOn();

    movingThresholdFilter->SetInput(0, movingMaskReader->GetOutput() );
    movingThresholdFilter->SetOutsideValue(0);
    movingThresholdFilter->ThresholdOutside(1, 255);
    movingThresholdFilter->ReleaseDataFlagOn();

    // Resample
    typename FixedNNInterpolator::Pointer fixednninterp = FixedNNInterpolator::New();
    typename MovingNNInterpolator::Pointer movingnninterp = MovingNNInterpolator::New();

    typename FixedResampleType::Pointer  fixedResample = FixedResampleType::New();
    fixedResample->SetInput( fixedThresholdFilter->GetOutput() );
    fixedResample->SetOutputParametersFromImage( fixedReader->GetOutput() );
    fixedResample->SetInterpolator( fixednninterp );
    fixedResample->SetDefaultPixelValue( 0 );
    fixedResample->ReleaseDataFlagOn();
    fixedResample->Update();

    typename MovingResampleType::Pointer  movingResample = MovingResampleType::New();
    movingResample->SetInput( movingThresholdFilter->GetOutput() );
    movingResample->SetOutputParametersFromImage( movingReader->GetOutput() );
    movingResample->SetInterpolator( movingnninterp );
    movingResample->SetDefaultPixelValue( 0 );
    movingResample->ReleaseDataFlagOn();
    movingResample->Update();

    fixedSpatialObjectMask->SetImage( fixedResample->GetOutput());
    movingSpatialObjectMask->SetImage( movingResample->GetOutput());
    }    

  // user decide if the input images need to be smoothed

  // Reorient to axials to avoid issues with registration metrics not
  // transforming image gradients with the image orientation in
  // calculating the derivative of metric wrt transformation
  // parameters.
  //
  // Forcing image to be axials avoids this problem. Note, that
  // reorientation only affects the internal mapping from index to
  // physical coordinates.  The reoriented data spans the same
  // physical space as the original data.  Thus, the registration
  // transform calculated on the reoriented data is also the
  // transform forthe original un-reoriented data.

  typename FixedOrientFilterType::Pointer orientFixed = FixedOrientFilterType::New(); // ##
  itk::PluginFilterWatcher watchOrientFixed(orientFixed,   "Orient Fixed Image",  CLPProcessInformation,  1.0 / 5.0,
                                            0.0);
  orientFixed->UseImageDirectionOn();
  orientFixed->SetDesiredCoordinateOrientationToAxial();

  if( FixedImageSmoothingFactor != 0 )
    {
    typedef itk::BinomialBlurImageFilter<FixedImageType, FixedImageType> BinomialFixedType;
    typename BinomialFixedType::Pointer BinomialFixed = BinomialFixedType::New();
    BinomialFixed->SetInput( fixedReader->GetOutput() );
    BinomialFixed->SetRepetitions( FixedImageSmoothingFactor * 2);
    itk::PluginFilterWatcher watchfilter(BinomialFixed, "Binomial Filter Fixed",  CLPProcessInformation, 1.0 / 5.0,
                                         1.0 / 5.0);
    BinomialFixed->Update();
    orientFixed->SetInput(BinomialFixed->GetOutput() );
    }
  else
    {
    orientFixed->SetInput( fixedReader->GetOutput() );
    }
  collector.Start( "Orient fixed volume" );
  orientFixed->Update();
  collector.Stop( "Orient fixed volume" );

  typename MovingOrientFilterType::Pointer orientMoving = MovingOrientFilterType::New(); // ##
  itk::PluginFilterWatcher watchOrientMoving(orientMoving,  "Orient Moving Image", CLPProcessInformation,  1.0 / 5.0,
                                             2.0 / 5.0);
  orientMoving->UseImageDirectionOn();
  orientMoving->SetDesiredCoordinateOrientationToAxial();

  if( MovingImageSmoothingFactor != 0 )
    {
    typedef itk::BinomialBlurImageFilter<MovingImageType,  MovingImageType> BinomialMovingType;
    typename BinomialMovingType::Pointer BinomialMoving = BinomialMovingType::New();
    BinomialMoving->SetInput(  movingReader->GetOutput() );
    BinomialMoving->SetRepetitions( MovingImageSmoothingFactor * 2);
    itk::PluginFilterWatcher watchfilter(BinomialMoving, "Binomial Filter Moving",  CLPProcessInformation, 1.0 / 5.0,
                                         3.0 / 5.0);
    BinomialMoving->Update();
    orientMoving->SetInput(BinomialMoving->GetOutput() );
    }
  else
    {
    orientMoving->SetInput( movingReader->GetOutput() );
    }

  collector.Start( "Orient moving volume" );
  orientMoving->Update();
  collector.Stop( "Orient moving volume" );

  // if an initial transform was specified, read it
  typedef itk::TransformFileReader TransformReaderType;
  TransformReaderType::Pointer initialTransform;

  if( InitialTransform != "" && !TestingMode )
    {
    initialTransform = TransformReaderType::New();
    initialTransform->SetFileName( InitialTransform );
    try
      {
      initialTransform->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Set up the optimizer
  //
  //
  ScheduleCommand::Pointer Schedule = ScheduleCommand::New();
  Schedule->SetSchedule( Iterations, LearningRate );

  int sum = 0;
  for( size_t i = 0; i < Iterations.size(); i++ )
    {
    sum += Iterations[i];
    }

  typename OptimizerType::Pointer optimizer;
  typename CathMatchOptimizerType::Pointer cathMatchOptimizer;
  typename TransformType::Pointer transform = TransformType::New();
  typename CathMatchTransformType::Pointer cathMatchTransform = CathMatchTransformType::New();

  if (CathMatch)
    {
    cathMatchOptimizer     = CathMatchOptimizerType::New();
    cathMatchOptimizer->SetNumberOfIterations( sum );
    cathMatchOptimizer->SetLearningRate( LearningRate[0] );
    cathMatchOptimizer->AddObserver( itk::IterationEvent(), Schedule );
    
    cathMatchTransform = CathMatchTransformType::New();
    //typedef CathMatchOptimizerType::ScalesType CathMatchOptimizerScalesType;
    CathMatchOptimizerScalesType scales( cathMatchTransform->GetNumberOfParameters() );
    scales.Fill( 1.0 );
    scales[0] = 1.0 / vnl_math_sqr(TranslationScale);
    cathMatchOptimizer->SetScales( scales );
    cathMatchOptimizer->SetMinimize( true );
    }
  else
    {
    optimizer     = OptimizerType::New();
    optimizer->SetNumberOfIterations( sum );
    optimizer->SetLearningRate( LearningRate[0] );
    optimizer->AddObserver( itk::IterationEvent(), Schedule );
    
    transform = TransformType::New();
    //typedef OptimizerType::ScalesType OptimizerScalesType;
    OptimizerScalesType scales( transform->GetNumberOfParameters() );
    scales.Fill( 1.0 );
    for( unsigned j = 4; j < 7; j++ )
      {
      scales[j] = 1.0 / vnl_math_sqr(TranslationScale);
      }
    optimizer->SetScales( scales );
    optimizer->SetMinimize( true );
    }


  // Initialize the transform
  //
  //
  if (CathMatch)
    {
    CathMatchTransformType::InputVnlVectorType axisOrigin;
    CathMatchTransformType::InputVnlVectorType axisNormal;
    CathMatchTransformType::InputVnlVectorType baseTranslation;
    CathMatchTransformType::VnlQuaternionType  baseRotation;
    CathMatchTransformType::InputVnlVectorType movingCatheterOrigin;
    CathMatchTransformType::InputVnlVectorType movingCatheterNormal;

    typedef itk::MatrixOffsetTransformBase<double, 3, 3> DoubleMatrixOffsetType;
    typedef itk::MatrixOffsetTransformBase<float, 3, 3>  FloatMatrixOffsetType;
    DoubleMatrixOffsetType::Pointer fda
      = dynamic_cast<DoubleMatrixOffsetType *>(fixedCatheterTransform.GetPointer() );
    FloatMatrixOffsetType::Pointer ffa
      = dynamic_cast<FloatMatrixOffsetType *>(fixedCatheterTransform.GetPointer() );
    DoubleMatrixOffsetType::Pointer mda
      = dynamic_cast<DoubleMatrixOffsetType *>(movingCatheterTransform.GetPointer() );
    FloatMatrixOffsetType::Pointer mfa
      = dynamic_cast<FloatMatrixOffsetType *>(movingCatheterTransform.GetPointer() );
    

    ////
    // NOTE: Since the CLI imports inverse transforms, we need to calculate inverse of
    //  the imported transforms for catheters. 
    // 
    // TODO: Check file type before setting the pointer to AssineTransformType
    //
    typedef itk::AffineTransform< double, 3 > AffineTransformType;
    AffineTransformType::Pointer atrans = AffineTransformType::New();
    AffineTransformType::InverseTransformBaseType::Pointer itrans;
    AffineTransformType::Pointer iatrans;
    AffineTransformType::MatrixType matrix;
    AffineTransformType::InverseMatrixType imatrix;
    AffineTransformType::OffsetType offset;
    
    atrans->SetIdentity();
    atrans->SetParameters(fixedCatheterTransform->GetParameters());
    std::cerr << "transform: " << atrans->GetMatrix() << std::endl;
    itrans = atrans->GetInverseTransform();
    iatrans = dynamic_cast<AffineTransformType *>(itrans.GetPointer() );
    matrix = iatrans->GetMatrix();
    offset = iatrans->GetOffset();
    for (int i = 0; i < 3; i ++)
      {
      axisNormal[i] = matrix[i][2];
      axisOrigin[i] = offset[i];
      }
    
    atrans->SetIdentity();
    atrans->SetParameters(movingCatheterTransform->GetParameters());
    itrans = atrans->GetInverseTransform();
    iatrans = dynamic_cast<AffineTransformType *>(itrans.GetPointer() );
    matrix = iatrans->GetMatrix();
    offset = iatrans->GetOffset();
    for (int i = 0; i < 3; i ++)
      {
      movingCatheterNormal[i] = matrix[i][2];
      movingCatheterOrigin[i] = offset[i];
      }

    axisNormal.normalize();
    movingCatheterNormal.normalize();

    // Compute base rotation in quaternion

    //  - 1. Find a vector orthogonal to the fixed and moving catheters by
    //       cross product
    CathMatchTransformType::InputVnlVectorType baseRotationAxis;
    baseRotationAxis = vnl_cross_3d(movingCatheterNormal, axisNormal);
    double baseAngleSin = baseRotationAxis.magnitude();
    double baseAngleCos = dot_product(movingCatheterNormal, axisNormal);
    double baseAngle = vcl_atan2(baseAngleSin, baseAngleCos);

    baseRotationAxis.normalize();
    baseRotation = CathMatchTransformType::VnlQuaternionType(baseRotationAxis, baseAngle);
    baseTranslation = axisOrigin - movingCatheterOrigin;

    std::cout << "baseAngle:        " << baseAngle << std::endl;
    std::cout << "Transform Axis:   " << axisNormal << std::endl;
    std::cout << "Transform Origin: " << axisOrigin << std::endl;
    std::cout << "moving Catheter Axis:   " << movingCatheterNormal << std::endl;
    std::cout << "moving Catheter Origin: " << movingCatheterOrigin << std::endl;
    std::cout << "Base Translation: " << baseTranslation << std::endl;
    std::cout << "Base Rotation:    " << baseRotation << std::endl;
    cathMatchTransform->SetTransformAxis(axisOrigin, axisNormal);
    cathMatchTransform->SetBaseTransform(baseTranslation, baseRotation);
    TransformType::ParametersType param = cathMatchTransform->GetParameters();
    param[0] = 2.0;
    param[1] = 1.0;
    cathMatchTransform->SetParameters(param);
    std::cerr << "returned parameter = " << cathMatchTransform->GetParameters() << std::endl;

    param[0] = 0.0;
    param[1] = 0.0;
    cathMatchTransform->SetParameters(param);
    std::cerr << "Matrix = "
              << cathMatchTransform->GetMatrix() << std::endl;
    std::cerr << "Offset = "
              << cathMatchTransform->GetOffset() << std::endl;
    std::cerr << "Center = "
              << cathMatchTransform->GetCenter() << std::endl;
    }
  else
    {
    typename TransformType::InputPointType centerFixed;
    typename FixedImageType::RegionType::SizeType sizeFixed =
      orientFixed->GetOutput()->GetLargestPossibleRegion().GetSize();
    // Find the center
    ContinuousIndexType indexFixed;
    for( unsigned j = 0; j < 3; j++ )
      {
      indexFixed[j] = (sizeFixed[j] - 1) / 2.0;
      }
    orientFixed->GetOutput()->TransformContinuousIndexToPhysicalPoint( indexFixed, centerFixed );
    
    typename TransformType::InputPointType centerMoving;
    typename MovingImageType::RegionType::SizeType sizeMoving =
      orientMoving->GetOutput()->GetLargestPossibleRegion().GetSize();
    // Find the center
    ContinuousIndexType indexMoving;
    for( unsigned j = 0; j < 3; j++ )
      {
      indexMoving[j] = (sizeMoving[j] - 1) / 2.0;
      }
    orientMoving->GetOutput()->TransformContinuousIndexToPhysicalPoint( indexMoving, centerMoving );
    
    transform->SetCenter( centerFixed );
    transform->Translate(centerMoving - centerFixed);
    std::cout << "Centering transform: "; transform->Print( std::cout );
    
    // If an initial transformation was provided, then grab the rigid
    // part and use it instead of the the centering transform.
    // (Should this be instead of the centering transform or composed
    // with the centering transform?
    //
    if( !TestingMode )
      {
      if( InitialTransform != ""
          && initialTransform->GetTransformList()->size() != 0 )
        {
        TransformReaderType::TransformType::Pointer initial
          = *(initialTransform->GetTransformList()->begin() );
    
        // most likely, the transform coming in is a subclass of
        // MatrixOffsetTransformBase
        typedef itk::MatrixOffsetTransformBase<double, 3, 3> DoubleMatrixOffsetType;
        typedef itk::MatrixOffsetTransformBase<float, 3, 3>  FloatMatrixOffsetType;
    
        DoubleMatrixOffsetType::Pointer da
          = dynamic_cast<DoubleMatrixOffsetType *>(initial.GetPointer() );
        FloatMatrixOffsetType::Pointer fa
          = dynamic_cast<FloatMatrixOffsetType *>(initial.GetPointer() );
    
        if( da )
          {
          vnl_svd<double> svd(da->GetMatrix().GetVnlMatrix() );
    
          transform->SetMatrix( svd.U() * vnl_transpose(svd.V() ) );
          transform->SetOffset( da->GetOffset() );
          }
        else if( fa )
          {
          vnl_matrix<double> t(3, 3);
          for( int i = 0; i < 3; ++i )
            {
            for( int j = 0; j < 3; ++j )
              {
              t.put(i, j, fa->GetMatrix().GetVnlMatrix().get(i, j) );
              }
            }
    
          vnl_svd<double> svd( t );
    
          transform->SetMatrix( svd.U() * vnl_transpose(svd.V() ) );
          transform->SetOffset( fa->GetOffset() );
          }
        else
          {
          std::cout << "Initial transform is an unsupported type.\n";
          }
    
        std::cout << "Initial transform: "; transform->Print( std::cout );
        }
      }
    }
  

  // Set up the metric
  //
  typename MetricType::Pointer  metric        = MetricType::New();
  metric->SetNumberOfHistogramBins( HistogramBins );
  metric->SetNumberOfSpatialSamples( SpatialSamples );
  metric->ReinitializeSeed(314159265);

  if (fMask)
    {
    metric->SetFixedImageMask( fixedSpatialObjectMask );
    metric->SetMovingImageMask( movingSpatialObjectMask );
    }

  // Create the interpolator
  //
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  // Set up the registration
  //
  typename RegistrationType::Pointer registration = RegistrationType::New();
  if (CathMatch)
    {
    registration->SetTransform( cathMatchTransform );
    registration->SetInitialTransformParameters( cathMatchTransform->GetParameters() );
    registration->SetMetric( metric );
    registration->SetOptimizer( cathMatchOptimizer );
    }
  else
    {
    registration->SetTransform( transform );
    registration->SetInitialTransformParameters( transform->GetParameters() );
    registration->SetMetric( metric );
    registration->SetOptimizer( optimizer );
    }

  registration->SetInterpolator( interpolator );
  registration->SetFixedImage( orientFixed->GetOutput() );
  registration->SetMovingImage( orientMoving->GetOutput() );

  Schedule->SetRegistration( registration );

  try
    {
    itk::PluginFilterWatcher watchRegistration(registration,
                                               "Registering",
                                               CLPProcessInformation,
                                               1.0 / 5.0, 4.0 / 5.0);
    collector.Start( "Register" );
    registration->Update();
    collector.Stop( "Register" );
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << err << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    return EXIT_FAILURE;
    }


  if (CathMatch)  
    {
    std::cerr << "Last parameter = "
              << registration->GetLastTransformParameters() << std::endl;
    TransformType::ParametersType parameters = registration->GetLastTransformParameters();
    parameters[0] += TestOffset;
    parameters[1] += TestAngleOffset/3.14159265358979;
    cathMatchTransform->SetParameters( parameters);
    std::cerr << "Matrix = "
              << cathMatchTransform->GetMatrix() << std::endl;
    std::cerr << "Offset = "
              << cathMatchTransform->GetOffset() << std::endl;
    std::cerr << "Center = "
              << cathMatchTransform->GetCenter() << std::endl;
    transform->SetMatrix( cathMatchTransform->GetMatrix() );
    transform->SetTranslation( cathMatchTransform->GetTranslation() );
    transform->SetOffset( cathMatchTransform->GetOffset() );
    transform->SetCenter( cathMatchTransform->GetCenter() );
    }
  else
    {
    transform->SetParameters( registration->GetLastTransformParameters() );
    }
    
  // compute eular angle and displacement
  double Ephi;          // roll
  double Etheta;        // pitch
  double Epsi;          // yaw
  double Edisplacement; // displacement
  
  //TransformType::ParametersType Eparams = registration->GetLastTransformParameters();
  TransformType::ParametersType Eparams = transform->GetParameters();
  Edisplacement = sqrt( Eparams[4] * Eparams[4] + Eparams[5] * Eparams[5] + Eparams[6] * Eparams[6] );
  Ephi =
    atan(2
         * (Eparams[0] * Eparams[1] + Eparams[2]
            * Eparams[3]) / (1 - 2 * (Eparams[1] * Eparams[1] + Eparams[2] * Eparams[2]) ) );
  Etheta = asin(2 * (Eparams[0] * Eparams[2] - Eparams[1] * Eparams[3]) );
  Epsi =
    atan(2
         * (Eparams[0] * Eparams[3] + Eparams[1]
            * Eparams[2]) / (1 - 2 * (Eparams[2] * Eparams[2] + Eparams[3] * Eparams[3]) ) );
  
  Ephi *= 180 / 3.14159265;
  Etheta *= 180 / 3.14159265;
  Epsi *= 180 / 3.14159265;
  
  std::cout << "\nRoll: " << Ephi;
  std::cout << " Pitch: " << Etheta;
  std::cout << " Yaw: " << Epsi;
  std::cout << " Displacement: " << Edisplacement << std::endl;


  if( OutputTransform != "" )
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer outputTransformWriter;
    
    outputTransformWriter = TransformWriterType::New();
    outputTransformWriter->SetFileName( OutputTransform );
    outputTransformWriter->SetInput( transform );

    try
      {
      outputTransformWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    }

  transform->Print(std::cout);

  // Resample to the original coordinate frame (not the reoriented
  // axial coordinate frame) of the fixed image
  //
  if( ResampledImageFileName != "" )
    {
    typename ResampleType::Pointer resample = ResampleType::New();
    typename ResampleInterpolatorType::Pointer Interpolator = ResampleInterpolatorType::New();
    itk::PluginFilterWatcher watchResample(resample,
                                           "Resample",
                                           CLPProcessInformation,
                                           1.0 / 5.0, 4.0 / 5.0);

    resample->SetInput( movingReader->GetOutput() );
    if (CathMatch)
      {
      resample->SetTransform( cathMatchTransform );
      }
    else
      {
      resample->SetTransform( transform );
      }
    resample->SetInterpolator( Interpolator );

    // Set the output sampling based on the fixed image.
    // ResampleImageFilter needs an image of the same type as the
    // moving image.
    typename MovingImageType::Pointer fixedInformation = MovingImageType::New();
    fixedInformation->CopyInformation( fixedReader->GetOutput() );
    resample->SetOutputParametersFromImage( fixedInformation );

    collector.Start( "Resample" );
    resample->Update();
    collector.Stop( "Resample" );

    typename WriterType::Pointer resampledWriter = WriterType::New();
    resampledWriter->SetFileName( ResampledImageFileName.c_str() );
    resampledWriter->SetInput( resample->GetOutput() );
    try
      {
      collector.Start( "Write volume" );
      resampledWriter->Write();
      collector.Stop( "Write volume" );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Report the time taken by the registration
  collector.Report();

  // In testing mode, take the image corner, and compute the error as the
  // difference between the locations of points after transforming with the
  // recovered transform and the ground truth transform.
  if( TestingMode )
    {
    // estimate the error of registration at index (0,0,0)
    ContinuousIndexType testIndex;
    typename TransformType::InputPointType testPoint;
    typename TransformType::OutputVectorType errorVector;

    testIndex[0] = 0.;
    testIndex[1] = 0.;
    testIndex[2] = 0.;

    orientFixed->GetOutput()->TransformContinuousIndexToPhysicalPoint(testIndex, testPoint);
    typename TransformType::OutputPointType groundTruthPoint, recoveredPoint;
    groundTruthPoint = groundTruthTransform->TransformPoint(testPoint);
    recoveredPoint = transform->TransformPoint(testPoint);
    errorVector[0] = groundTruthPoint[0] - recoveredPoint[0];
    errorVector[1] = groundTruthPoint[1] - recoveredPoint[1];
    errorVector[2] = groundTruthPoint[2] - recoveredPoint[2];
    std::cout << "Magnitude of error vector: " << errorVector.GetNorm() << std::endl;

    if( errorVector.GetNorm() > TESTMODE_ERROR_TOLERANCE )
      {
      std::cerr << "Registration error in testing mode exceeds threshold "
                << TESTMODE_ERROR_TOLERANCE << std::endl;

      std::cerr << "Recovered transform: " << std::endl;
      transform->Print(std::cerr);

      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}

template <class T>
int DoIt( int argc, char * argv[], const T& targ)
{

  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(MovingImageFileName, pixelType, componentType);

    // This filter handles all types

    switch( componentType )
      {
      case itk::ImageIOBase::CHAR:
      case itk::ImageIOBase::UCHAR:
      case itk::ImageIOBase::SHORT:
        return DoIt2( argc, argv, targ, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::INT:
        return DoIt2( argc, argv, targ, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::ULONG:
        return DoIt2( argc, argv, targ, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt2( argc, argv, targ, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt2( argc, argv, targ, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt2( argc, argv, targ, static_cast<float>(0) );
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
  return EXIT_FAILURE;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{

  PARSE_ARGS;

  // this line is here to be able to see the full output on the dashboard even
  // when the test succeeds (to see the reproducibility error measure)
  std::cout << std::endl << "ctest needs: CTEST_FULL_OUTPUT" << std::endl;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(FixedImageFileName, pixelType, componentType);

    // This filter handles all types

    switch( componentType )
      {
      case itk::ImageIOBase::CHAR:
      case itk::ImageIOBase::UCHAR:
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::UINT:
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
        return DoIt( argc, argv, static_cast<float>(0) );
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
