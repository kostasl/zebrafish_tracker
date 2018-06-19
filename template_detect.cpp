///*
/// Uses Template Image To detect Matching and orientation of fish Body
/// \todo Need to Look into feature descriptors for matching as a faster method.
/// OpenCV includes some ready made, ORB being free. I Found an intuitive and external to OPENCV called FREAK
/// which could be implemented (uses a retina inspired sampling pattern to encode a binary string of the template)
///*



#include <template_detect.h>
#include <larvatrack.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
#include <opencv2/highgui/highgui.hpp>


/// CUDA //

#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    #include <opencv2/opencv_modules.hpp> //THe Cuda Defines are in here
    #include "opencv2/cudaimgproc.hpp"
    #include "opencv2/cudaarithm.hpp"
    #include <opencv2/core/cuda.hpp>
    #include <opencv2/photo/cuda.hpp>
    #include <opencv2/core/cuda_types.hpp>
#endif

#include <QString>
#include <random>
#include <QDirIterator>
#include <QDir>
#include <QDebug>





//#include <opencv2/cudaimgproc.hpp> //Template Matching \todo Compile OpenCv With CUDA Support

extern double gTemplateMatchThreshold;
extern int gFishTemplateAngleSteps;
extern int gnumberOfTemplatesInCache;
extern cv::Mat gFishTemplateCache;
extern MainWindow* pwindow_main;
extern bool bTemplateSearchThroughRows;

#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    extern cv::Ptr<cv::cuda::TemplateMatching> gpu_MatchAlg;
#endif

static cv::Mat loadImage(const std::string& name)
{
    cv::Mat image = cv::imread(name, cv::IMREAD_GRAYSCALE);
    if (image.empty())
    {
        std::cerr << "Can't load image - " << name << std::endl;
        exit(-1);
    }
    return image;
}


///
/// \brief makeTemplateCache Creates a rotational replicate of the template across angles of size
/// stores them in consecutive blocks. Saves time from having to rotate the template EveryTime
///  Handles Symmetry Issues, If Template Is not Square, It copies to a square canvas before Rotation
/// \param templateIn
/// \param imgTemplateOut
/// \param iAngleStepDeg (minimum 1 degree)
///
void makeTemplateVar(cv::Mat& templateIn,cv::Mat& imgTemplateOut, int iAngleStepDeg)
{
    int iAngleIncrements = 360.0/iAngleStepDeg;
    int ifishtemplateAngle = 0;
    //Allocate Dark/Black Mat of appropriate Size to Fit All Replications of the templ.
    //Space Them Out Along Long Dim So no cutting Occurs
    int mxDim = std::max(templateIn.cols,templateIn.rows);
    imgTemplateOut = cv::Mat::zeros(mxDim,mxDim*iAngleIncrements,CV_8UC1);
    cv::Point tempCentre = cv::Point(templateIn.cols/2,templateIn.rows/2);
    cv::Point ptbottomRight = cv::Point(mxDim,mxDim); //Square fitting Largest Dimension
    cv::Point cntrdiff = cv::Point(mxDim/2,mxDim/2)-tempCentre ;
    cv::Rect templRegion(cv::Point(0,0),ptbottomRight); //Define A Rect region to snatch SubMatrices
    cv::Point rotCentre = cv::Point(templRegion.width/2,templRegion.height/2);

    //Find Alignment OffSet - Correct Template In Case Its Tilded Already //
    int iAngleOffset = 0;
    cv::Point ptCentreCorrection(0,0);
    /// Detect edges using Threshold
    cv::Mat templ_thres,template_blur;
    std::vector<std::vector<cv::Point> > contours;
    std::vector<cv::Vec4i> hierarchy;

    //Segment Fish Body Blob
    //cv::threshold( templateIn, templ_thres,thres, 255, cv::THRESH_BINARY );
    cv::GaussianBlur(templateIn,template_blur,cv::Size(7,7),5,5);
    cv::adaptiveThreshold(template_blur, templ_thres,255,cv::ADAPTIVE_THRESH_MEAN_C,cv::THRESH_BINARY,21,0); //Last Param Is const substracted from mean
    //ADAPTIVE_THRESH_MEAN_C

    /// Find contours
    cv::findContours( templ_thres, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

    //Should be one contour
    if (contours.size() == 1)
    {
        cv::RotatedRect rectBody = fitEllipse(contours[0]);
        // Make Slight Angle Corrections 0.5 Corrected
        iAngleOffset = rectBody.angle; //Make 0 Angle  the vertical image axis

        if  ((iAngleOffset) > 90) iAngleOffset-=180;
        if  ((iAngleOffset) < -90) iAngleOffset+=180;


        //Debug Show Shape Centre //
        cv::circle(templ_thres,rectBody.center,3,100,1);
        //Correct For Centre Offset

        ptCentreCorrection.x = tempCentre.x - rectBody.center.x;
       // ptCentreCorrection.y //= tempCentre.x - rectBody.center.x;
        std::clog << "[info] makeTemplateVar: correct Template DAngle : " << iAngleOffset << " DX:" << ptCentreCorrection.x <<  std::endl;

    }else
    {
        std::clog << "[warning] makeTemplateVar: multiple contours detected on template " << std::endl;
    }


    // DEBUG //
    //std::string winname= QString(QString("TemplShape")+QString::number(gnumberOfTemplatesInCache)).toStdString();
    //cv::imshow(winname ,templ_thres);



    //Got through Each Angle
    for (int i=0;i<iAngleIncrements;i++)
    {
        //Make Rotation MAtrix
        cv::Mat Mrot = cv::getRotationMatrix2D(rotCentre,360.0-ifishtemplateAngle+iAngleOffset*0.5,1.0);

        //Copy to Larger Square

        //Get icon SubMatrix from Large Template Mat
        cv::Mat templ_rot = imgTemplateOut(templRegion);
        //ReLocate Template To Centre Of Large Canvas, Before Rotation
        templateIn.copyTo(templ_rot(cv::Rect(0,0,templateIn.cols,templateIn.rows)+cntrdiff+ptCentreCorrection ));
        //Make Rotation Transformation
        cv::warpAffine(templ_rot,templ_rot,Mrot,templRegion.size());

        //Slide To Next Empty region box for the Next Angle
        templRegion.x +=templRegion.width;
        ifishtemplateAngle += iAngleStepDeg;
    }



}


///
/// \brief templatefindFishInImage Scans input image for templates and returns the best location
///  which that exceed a threshold value for a match
/// \param templRegion Rect of template img Size to look for within the larger template Cache
/// \param startRow - Optimization So search begins from the most likely Template as of the last one
/// \param startCol - Optimization So search begins from the most likely Template Angle
/// \param findFirstMatch if true It Looks for 1st template row that exceeds threshold - otherwise it looks for best match through all cache
/// \note The calling Function needts reposition maxLoc To the global Frame, if imgGreyIn is a subspace of the image
/// if Row scanning is disabled when bTemplateSearchThroughRows is not set
/// Use of UMat for matchTemplate is superfluous , as the GPU is not Utilized and it causes a dealloc mem bug to trigger -  Removed UMat
/// A gpuAssisted  function for this is included in the bottom of the file.
int templatefindFishInImage(cv::Mat& imgRegionIn,cv::Mat& imgtemplCache,cv::Size templSz, double& matchScore, cv::Point& locations_tl,int& startRow,int& startCol,bool findFirstMatch)
{
  const int iIdxAngleMargin = 3; //Offset Of Angle To begin Before LastKnownGood Angle
  int matchColIdx = 0;
  int Colidx = 0; //Current Angle Index Being tested in the loop

  assert(!imgRegionIn.empty());
  assert(templSz.height*startRow <= imgtemplCache.rows);
  assert(templSz.width*startCol <= imgtemplCache.cols);

  //startRow = 0;
  cv::Mat templ_rot; //The Matched template
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
        cv::cuda::GpuMat dimgRegionIn(imgRegionIn);
#endif


  int idRow = startRow;
  int ibestMatchRow  = startRow;
  double minVal;
  double maxVal = 0.0;
  double maxGVal = 0.0;
  cv::Point ptmaxLoc,ptminLoc;
  cv::Point ptGmaxLoc,ptGminLoc;
  cv::Mat outMatchConv; //Convolution Result
  //iAngles = imgtempl.cols/templRegion.width;
  //Slide over each little template icon


  cv::Point ptbottomRight = cv::Point(templSz.width,templSz.height);
  cv::Rect templRegion(cv::Point(0,0),ptbottomRight);

  // Start from Angle Region as set from last search
  if (startCol > iIdxAngleMargin)
      startCol -=iIdxAngleMargin; //Move to Template 3Angle Steps anticlockwise

 //Best Match Flag Forces A Full Search Through the Cache
  if (findFirstMatch) //Reset Search To Start From Top
  {
      startRow = 0;
      idRow = 0;
      startCol = 0;
  }

  //Initialiaze At last known Good Location
  templRegion.x = std::max(0,std::min(templSz.width*startCol,imgtemplCache.cols));
  templRegion.y = std::max(0,std::min(templSz.height*startRow,imgtemplCache.rows));
  Colidx = startCol;

  ///Run Through All rotated Templates - optional starting row for optimization
  /// \note For Speed Up - We Remove Running Through All Rows - Stick to LastKnown Good Row- And Let Random switch when this fails take Care Of Changing Row
  int iScanRowLimit = imgtemplCache.rows;

  if (!bTemplateSearchThroughRows) //Do not Search Subsequent Template Rows
      iScanRowLimit = std::min(templSz.height*startRow + templRegion.height,imgtemplCache.rows) ;

  for (int j=templSz.height*startRow; j<iScanRowLimit;j+=templRegion.height) //Remove for Speed Optimization.
  {
      templRegion.y    = j;
       /// Run Throught each  *Columns/Angle* (Ie Different Angles of this template
      //for (int i=templSz.width*startCol; i<imgtempl.cols;i+=templRegion.width)

      ///Limit Search To Within FIXED (15) Degrees/Templates from Starting Col.Point If Not Searching From The top
      while(templRegion.x < imgtemplCache.cols && ((Colidx-startCol) < TEMPLATE_COL_SEARCH_REGION || startCol==0)){
        //Obtain next Template At Angle
         imgtemplCache(templRegion).copyTo(templ_rot) ;

        //Run Lib MatchTemplate COnvolution Either On CPU OR GPU //
         // WIth Current Setup Using GPU for multiple template Matching/scanning is too slow
#if defined(USE_CUDA) && defined(USE_CUDA_FOR_TEMPLATE_MATCHING) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
        maxVal = gpu_matchTemplate(templ_rot,dimgRegionIn,ptmaxLoc);
#else
         cv::matchTemplate(imgRegionIn,templ_rot,outMatchConv, CV_TM_CCORR_NORMED  ); // CV_TM_CCOEFF_NORMED ,TM_SQDIFF_NORMED
         //Find Min Max Location
         cv::minMaxLoc(outMatchConv,&minVal,&maxVal,&ptminLoc,&ptmaxLoc);
#endif         //Convolution  // CV_TM_SQDIFF_NORMED Poor Matching

        //Assume Value < 0.7 is non Fish,

        if (maxGVal < maxVal)
        {
            maxGVal         = maxVal;
            ptGmaxLoc       = ptmaxLoc; //The calling Function needts reposition maxLoc To the global Frame
            matchColIdx     = Colidx;
            ibestMatchRow   = idRow;
            matchScore      = maxVal; //Save Score Of Best Match
            locations_tl    = ptGmaxLoc;

        }

        //Shift Region To Next Block
        if (findFirstMatch) //Skip Cols If Just Scanning
        {
            templRegion.x +=3*templSz.width;
            Colidx+=3;
        }
        else
        {
            templRegion.x +=templSz.width;
            Colidx++;
        }
      } //Loop Through Columns/Angle


   idRow++;             //We Start Again From The Next Row
   Colidx               = 0;
   templRegion.x        = 0; //ReStart from 1st col


    //Dont scan all Rows Just Check If Found on this One before Proceeding

   ///Check If Matching Exceeeds threshold And Stop Loops - Return Found Column Idx//
   if (maxGVal >= gTemplateMatchThreshold && !findFirstMatch && maxGVal > 0.1)
   {
       //Save Results To Output
       //matchScore    = maxGVal;
       //locations_tl  = ptGmaxLoc;
       startCol = Colidx; //Save AS Best - For next Iteration - Only If Exceeds Threshold
       break; ///Stop The loop Rows search Here

    }else //Reset StartCol Hint, So Next Time The whole Cache Row is scanned
        startCol = 0;

//   else{ //Nothing Found YEt-- Proceed To Next Template variation
//       matchColIdx  = 0;
//       matchScore   = maxGVal;
//       locations_tl = cv::Point(0,0);
//       //Didnt Find Template /Try Next Row
//       //startRow = 0;//Start From Top Of All Templates On Next Search


//   }



 } //Loop Through Rows - Starting From Last Call Best Match Row


  if (startRow != ibestMatchRow) //Store Best Row Match
  {
      std::stringstream ss;
      // Log As Message //
      ss << "Ch. Templ. Row:" << startRow << " -> "  << ibestMatchRow;
      pwindow_main->LogEvent(QString::fromStdString(ss.str()));

      startRow = ibestMatchRow;
      //matchColIdx = ibestMatch;
      //cv::imshow("Templ",templ_rot);
  }


 if (maxGVal < gTemplateMatchThreshold)
 {

     startRow = (rand() % static_cast<int>(gnumberOfTemplatesInCache - 0 + 1));//Start From RANDOM rOW On Next Search
     startCol = 0;

     std::stringstream ss;
     // Log As Message //
     if (ibestMatchRow < gnumberOfTemplatesInCache)
        ss << "Found row:" << ibestMatchRow << " but gives Low Match-pick next Randomly ->"  << startRow;
     else
     {
         ss << "Reached End of Templates row:" << ibestMatchRow << " Start Over on next";
         startRow = 0;
     }
     pwindow_main->LogEvent(QString::fromStdString(ss.str()));

 }

 if (templ_rot.cols > 0 && templ_rot.rows > 0)
         //cv::imshow("MTemplate",templ_rot);
        pwindow_main->showInsetTemplateimg(templ_rot); //Show The template Image Used
 else
     std::cerr << "template_detect: templ_rot is empty - matchcol:" << matchColIdx << " Row: " << ibestMatchRow << std::endl;

  return matchColIdx;
}


///
///\brief addTemplateToCache
///\note assumes all Templates are the same size
///
int addTemplateToCache(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl)
{

    //Make Variations And store in template Cache
    cv::Mat fishTemplateVar;
    cv::Mat mtCacheRow,mtEnlargedCache;
    makeTemplateVar(imgTempl,fishTemplateVar, gFishTemplateAngleSteps);

    ///Initialize The Cache if this the 1st Template added
    if (idxTempl == 0)
        FishTemplateCache = cv::Mat::zeros(fishTemplateVar.rows,fishTemplateVar.cols,CV_8UC1);
    else{
        //Copy COntents To Enlarged Cache and replace pointer
        mtEnlargedCache = cv::Mat::zeros(FishTemplateCache.rows+fishTemplateVar.rows,fishTemplateVar.cols,CV_8UC1);
        //Get Ref To Old Sized Cache Only
        mtCacheRow = mtEnlargedCache(cv::Rect(0,0,FishTemplateCache.cols,FishTemplateCache.rows));
        FishTemplateCache.copyTo(mtCacheRow); //Copy Old Cache into New replacing that part of empty cache
        mtEnlargedCache.copyTo(FishTemplateCache); //Copy Back So gFishTemplateCache = mtEnlargedCache;

        mtEnlargedCache.release();
        mtEnlargedCache.deallocate();
    }
     //Fill The Last (New Row) In The Cache
    mtCacheRow = FishTemplateCache(cv::Rect(0,fishTemplateVar.rows*(idxTempl),fishTemplateVar.cols,fishTemplateVar.rows));
    fishTemplateVar.copyTo(mtCacheRow); //Copy To Row In CAche
    gnumberOfTemplatesInCache++; //Increment Count

    // DEBUG //
    //cv::imshow("Fish Template",FishTemplateCache(cv::Rect(0,0,std::max(imgTempl.cols,imgTempl.rows),FishTemplateCache.rows) ));
    //cv::imshow("Templ",imgTempl);
    //  //

    std::clog << "New Template added, Templ. Count now:" << gnumberOfTemplatesInCache << std::endl;


   return ++idxTempl;
}

///
/// \brief deleteLastTemplateRow
/// \param FishTemplateCache
/// \param idxTempl - Row To Remove
/// \return
///
int deleteTemplateRow(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl)
{
    //Draw Black
    int mxDim = std::max(imgTempl.cols,imgTempl.rows);
    //\note RECT constructor takes starting point x,y, size_w, size_h)
    cv::Rect rectblankcv(0,mxDim*(idxTempl),FishTemplateCache.cols,mxDim);
    cv::Mat mFishTemplate_local;// = FishTemplateCache.getMat(cv::ACCESS_WRITE);
    cv::rectangle(mFishTemplate_local,rectblankcv,CV_RGB(0,0,0),CV_FILLED); //Blank It OUt

    //If Removing Last Row, Then Its Simple
    //Shrink Template
    if (idxTempl == (gnumberOfTemplatesInCache-1))
    {
        FishTemplateCache  = FishTemplateCache(cv::Rect(0,0,FishTemplateCache.cols,mxDim*(idxTempl)));
        gnumberOfTemplatesInCache--;
    }else //Other Wise, We need to Cut and stich
    {
        //Cut In 2- Halves and rejoin
        cv::Mat mTop;
        cv::Mat mBottom;
        FishTemplateCache(cv::Rect(0,0,FishTemplateCache.cols,mxDim*(idxTempl))).copyTo(mTop) ;
        //cv::imshow("Template Top",mTop(cv::Rect(0,0,mxDim,mTop.rows)));

        FishTemplateCache(cv::Rect(0,mTop.rows+mxDim,mTop.cols,FishTemplateCache.rows-mTop.rows-mxDim)).copyTo(mBottom);
        //cv::imshow("Template Bottom",mBottom(cv::Rect(0,0,mxDim,mBottom.rows)));

        cv::Mat mtShrankCache   = cv::Mat::zeros(FishTemplateCache.rows-mxDim,FishTemplateCache.cols,CV_8UC1);

        gnumberOfTemplatesInCache--;
        mTop.copyTo(mtShrankCache(cv::Rect(0,0,FishTemplateCache.cols,mxDim*(idxTempl))));
        mBottom.copyTo( mtShrankCache( cv::Rect(0,mTop.rows,FishTemplateCache.cols,mtShrankCache.rows-mTop.rows) ));
        mtShrankCache.copyTo(FishTemplateCache);
        mtShrankCache.deallocate();

    }

    // DEBUG //
    cv::imshow("Updated Fish Template Cache",FishTemplateCache(cv::Rect(0,0,mxDim,FishTemplateCache.rows)));

    return 1;
}



int loadTemplatesFromDirectory(QString strDir)
{
    QDir dirTempl(strDir);
    cv::Mat templFrame;
    int fileCount = 0;
    if (!dirTempl.exists())
    {
        qWarning("Cannot find the a template directory");
        return 0;
    }

        QStringList fileFilters; fileFilters << "*.png" << "*.tiff" << "*.pgm" << "*.png";
        QStringList imgFiles = QDir(strDir).entryList(fileFilters,QDir::Files,QDir::Name);
        strDir.append('/');
        QListIterator<QString> itfile (imgFiles);
        while (itfile.hasNext())
        {
          QString filename = itfile.next();
          std::string filepath = filename.prepend(strDir ).toStdString();

          qDebug() << "*Load Template: " << filename;
          templFrame  = loadImage(filepath);
          addTemplateToCache(templFrame,gFishTemplateCache,gnumberOfTemplatesInCache);
          fileCount++;
        }


         qDebug() << "Loaded # " << fileCount << "Templates";
        return fileCount;
}



//////////////////////// MATCH TEMPLATE EXAMPLE FOR GPU ////////////////////////
/// \brief gpu_matchTemplate
/// \param templ_h
/// \param image_h IMage region TO search In
/// \param ptBestMatch
/// \return MatchTemplateSCore
///
#if defined(USE_CUDA) && defined(HAVE_OPENCV_CUDAARITHM) && defined(HAVE_OPENCV_CUDAIMGPROC)
    double gpu_matchTemplate(cv::Mat templ_h,cv::cuda::GpuMat& dimage,cv::Point& ptBestMatch)
    {
    //cv::cuda::setDevice(0);	//initialize CUDA

    ////cv::Mat image_h = cv::imread(	"/home/buddy/Documents/workspace/OpenCVTemplateMatch1/src/input_image.jpg");
    ////cv::Mat templ_h = cv::imread(				"/home/buddy/Documents/workspace/OpenCVTemplateMatch1/src/template_image.jpg");

    cv::cuda::GpuMat dtempl(templ_h); //upload image on gpu

    cv::cuda::GpuMat dresult;

    //dtempl.upload(image_h);
    //cv::Ptr<cv::cuda::TemplateMatching> alg = cv::cuda::createTemplateMatching(templ_h.type(), CV_TM_CCORR_NORMED);


    //cv::cuda::GpuMat dst;
    gpu_MatchAlg->match(dimage, dtempl, dresult);

    //cv::cuda::normalize(dresult, dresult, 0, 1, cv::NORM_MINMAX, -1);
    double max_value;

    //cv::Point location;
    //Find Best Match
    cv::cuda::minMaxLoc(dresult, 0, &max_value, 0, &ptBestMatch);


    ////copying back to host memory for display
    //cv::Mat result_h;
    //result.download(result_h);

    //Return The Match Value
    return max_value;

    }
#endif
