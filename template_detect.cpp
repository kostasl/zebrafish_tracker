
#include <template_detect.h>

extern double gMatchShapeThreshold;
extern int gFishTemplateAngleSteps;
extern int gnumberOfTemplatesInCache;


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

    //Got through Each Angle
    for (int i=0;i<iAngleIncrements;i++)
    {
        //Make Rotation MAtrix
        cv::Mat Mrot = cv::getRotationMatrix2D(rotCentre,360.0-ifishtemplateAngle,1.0);

        //Copy to Larger Square

        //Get icon SubMatrix from Large Template Mat
        cv::Mat templ_rot = imgTemplateOut(templRegion);
        //ReLocate Template To Centre Of Large Canvas, Before Rotation
        templateIn.copyTo(templ_rot(cv::Rect(0,0,templateIn.cols,templateIn.rows)+cntrdiff ));
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
/// \note The calling Function needts reposition maxLoc To the global Frame, if imgGreyIn is a subspace of the image
///
int templatefindFishInImage(cv::Mat& imgGreyIn,cv::Mat& imgtempl,cv::Size templSz, double& matchScore, cv::Point& locations_tl,int& startRow,int& startCol)
{
  int matchIdx;
  int idx = 0; //Current Angle Index Being tested in the loop
  int idRow = 0;
  double minVal, maxVal;
  double maxGVal = 0.0;
  cv::Point ptmaxLoc,ptminLoc;
  cv::Point ptGmaxLoc,ptGminLoc;
  cv::Mat outMatchConv; //Convolution Result
  //iAngles = imgtempl.cols/templRegion.width;
  //Slide over each little template icon
  assert(!imgGreyIn.empty());

  cv::Point ptbottomRight = cv::Point(templSz.width,templSz.height);
  cv::Rect templRegion(cv::Point(0,0),ptbottomRight);

  // Start from Angle Region as set from last search
  if (startCol > 3)
      startCol -=3; //Move to Template 3Angle Steps anticlockwise

  templRegion.x = templSz.height*startCol;
  //Run Through All rotated Templates - optional starting row for optimization
  for (int j=templSz.height*startRow; j<imgtempl.rows;j+=templRegion.height)
  {
      templRegion.y    = j;

      for (int i=templSz.height*startCol; i<imgtempl.cols;i+=templRegion.width)
      {
        //Obtain next Template At Angle
        cv::Mat templ_rot(imgtempl,templRegion);
        //Convolution
        cv::matchTemplate(imgGreyIn,templ_rot,outMatchConv,CV_TM_CCOEFF_NORMED);
        //Find Min Max Location
        cv::minMaxLoc(outMatchConv,&minVal,&maxVal,&ptminLoc,&ptmaxLoc);
        //Assume Value < 0.7 is non Fish,
        if (maxGVal < maxVal)
        {
            maxGVal     = maxVal;
            ptGmaxLoc   = ptmaxLoc; //The calling Function needts reposition maxLoc To the global Frame
            matchIdx   = idx;
        }

        //Shift Region To Next Block
        templRegion.x +=templSz.width;
        idx++;
      } //Loop Through Columns
   ///Check If Matching Exceeeds threshold
   if (maxGVal > gMatchShapeThreshold)
   {
       //Save Results To Output
       matchScore    = maxGVal;
       locations_tl  = ptGmaxLoc;
       startRow = idRow;
       //Stop Search
       return matchIdx; //Done Searching Stop Going Through loop

       //break;
    }else { //Nothing Found YEt-- Proceed To Next Template variation
       matchIdx = 0;
       matchScore = 0.0;
       locations_tl = cv::Point(0,0);
   }

   idRow++;
   idx              = 0;
   templRegion.x    = 0; //Start from 1st col again


 } //Loop Through Rows

  startRow = 0;//Start From Top Of All Templates On Next Search
  return matchIdx;
}


///
///\brief addTemplateToCache
///\note assumes all Templates are the same size
///
int addTemplateToCache(cv::Mat& imgTempl,cv::Mat& FishTemplateCache,int idxTempl)
{
    //Make Variations And store in template Cache
    cv::Mat fishTemplateVar,mtCacheRow,mtEnlargedCache;
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
    }
     //Fill The Last (New Row) In The Cache
    mtCacheRow = FishTemplateCache(cv::Rect(0,fishTemplateVar.rows*(idxTempl),fishTemplateVar.cols,fishTemplateVar.rows));
    fishTemplateVar.copyTo(mtCacheRow); //Copy To Row In CAche
    gnumberOfTemplatesInCache++; //Increment Count
    // DEBUG //
    cv::imshow("Fish Template",FishTemplateCache);

   return ++idxTempl;
}
