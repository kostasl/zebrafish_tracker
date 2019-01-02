///*
/// Implements Algorithm based on A NEW EFFICIENT ELLIPSE DETECTION METHOD, Yonghong Xie   ,IEEE, 2002
///
/// (1) Store all edge pixels in a one dimensional array.
/// (2) Clear the accumulator array .
/// (3) For each pixel (x1, y1 ), carry out the following steps from (4) to (14).
/// (4) For each other pixel (x2, y2), if the distance between (x1, y1) and (x 2, y2)
/// is greater than the required least distance  for  a  pair  of  pixels  to  be  considered  then
/// carry out the following steps from (5) to (14).
///
/// (5) From  the  pair  of  pixels  (x1,  y1) and  (x2,  y2),  using
/// equations   (1)   to   (4)   to   calculate   the   center,
/// orientation and length of major axis for the assumed ellipse.
///
/// (6) For  each  third  pixel  (x,  y),  if  the  distance  between
/// (x,  y)  and  (x0,  y0)   is  ?greater?  than  the  required  least
/// distance  for  a  pair  of  pixels  to  be  considered  :
///
/// "The distance between (x, y) and (x 0 , y 0 ) should be less than the distance between (x 1 , y 1 ) and (x 0 ,y 0 ) or between (x 2 , y 2 ) and (x 0 , y 0 ) ."
/// *found in MATlab implementation : ie 3rd point distance <= a; % (otherwise the formulae in paper do not work)
///  then carry out the following steps from (7) to (9).
/// (7)  Using  equations  (5)  and  (6)  to  calculate  the  length  of minor axis.
/// (8)  Increment  the  accumulator  for  this  length  of  minor  axis by 1.
/// (9)  Loop  until  all  pixels  are  computed  for  this  pair  of  pixels.
/// (10) Find the maxium element in accumulator array.
/// The related  length  is  the  possible  length  of  minor  axis
/// for  assumed  ellipse.  If  the  vote  is  greater  than  the
/// required   least   number   for   assumed   ellipse,   one  ellipse is detected.
/// (11)   Output ellipse parameters.
/// (12)   Remove the pixels on the detected ellipse from edge pixel array.
/// (13)   Clear accumulator array.
/// (14)   Loop until all pairs of pixels are computed.
/// (15)   Superimpose   detected   ellipses   on   the   original  image.
/// (16)   End.
///

/// Eqns:
/// x 0 = (x 1 + x 2 )/2  --(1)
/// y 0 = (y 1 + y 2 )/2  --(2)
/// a = [(x 2 – x 1 ) + (y 2 – y 1 ) ] /2 ---(3)
/// α = atan [(y 2 – y 1 )/(x 2 – x 1 )], (4)
/// b2 = (a 2 d 2 sin 2 τ)/( a 2 -d 2 cos 2 τ ) (5)
/// cos τ = ( a 2 + d 2 – f 2 )/(2ad) (6)

///Summary : Algorithm Checks a candidate ellipse with major axis between to pair of test points,
///  then estimates minor axis by testing all 3rd points and uses a voting procedure to check for possible minor axis and ellipse
#include <ellipse_detect.h>
#include <template_detect.h>
#include <larvatrack.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip> //for setprecision
#include <limits>
#include <string>
#include <random>
#include <config.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>

//extern MainWindow window_main;
extern MainWindow* pwindow_main;

extern bool bUseHistEqualization;
extern int gi_CannyThresSmall;
extern int gi_CannyThres;
extern int gi_VotesEllipseThres;
extern int gi_minEllipseMajor;
extern int gi_maxEllipseMajor;
extern int g_BGthresh;
extern int gEyeTemplateAngleSteps;
extern int giHeadIsolationMaskVOffset; //V Distance When Drawing Arc In getEyeSegThreshold
//cv::Mat imgDebug;

extern cv::Mat kernelOpenfish;
extern cv::Mat frameDebugC;
extern cv::Mat gEyeTemplateCache;

extern int gthresEyeSeg;

// Static Memory Buffers //
static cv::Mat imgIn_thres; // Crash Here  Frame:55200 RSS: 1100.57MB
static cv::Mat imgEdge_local; //Crash Here
static cv::Mat imgUpsampled_gray;
static cv::Mat img_colour;

inline int getMax(int* darray,int length,double& votes)
{
    double max=darray[0];
    int maxIdx = 0;
    //find max and mins
    for(int j=0; j<length; j++)
    {
        if(max<=darray[j])
        {
            max=darray[j];
            maxIdx = j;
        }
    }
    votes = max;
    return maxIdx;
}

/// Fills A list with  point coords where pixels (edges image) are above a threshold (non-zero)
void getEdgePoints(cv::Mat& imgEdgeIn,tEllipsoidEdges& vedgepoint)
{
   const float pxThres = 100.0; //threshold is non-zero
   //vedgepoint.clear();


//Split Image In Two
  for(int i=0; i<imgEdgeIn.rows; i++)
      for(int j=0; j<imgEdgeIn.cols; j++)
      { cv::Point pt(j,i); //x,y
           if ( imgEdgeIn.at<uchar>(pt) >  pxThres)
           {
               vedgepoint.push_back(tEllipsoidEdge(pt));
               //imgDebug.at<uchar>(pt) = 125;
           }
      }

}


/// Fills A list with  point coords where pixels (edges image) are above a threshold (non-zero)
void getEdgePoints(std::vector<cv::Point>& contour,tEllipsoidEdges& vedgepoint)
{
   //vedgepoint.clear();

//Split Image In Two
  for(int i=0; i<contour.size(); i++)
      {
               vedgepoint.push_back(tEllipsoidEdge(contour[i]));
               //imgDebug.at<uchar>(contour[i]) = 155;
       }

}

/// \todo Check Image Bounds
void drawEllipse(cv::Mat imgOut,tDetectedEllipsoid ellipse)
{

    cv::ellipse(imgOut,ellipse.rectEllipse,CV_RGB(250,50,50),1,cv::LINE_8);
    cv::circle(imgOut,ellipse.rectEllipse.center,1,CV_RGB(0,0,255),1);

    //    assert(ellipse.ptAxisMj2.y <= imgOut.rows && ellipse.ptAxisMj2.y >= 0);
    //    assert(ellipse.ptAxisMj2.x <= imgOut.cols && ellipse.ptAxisMj2.x >= 0);
    //    assert(ellipse.ptAxisMj1.y <= imgOut.rows && ellipse.ptAxisMj1.y >= 0);
    //    assert(ellipse.ptAxisMj1.x <= imgOut.cols && ellipse.ptAxisMj1.x >= 0);
    // Assertion Was Failing So I imposed the limits to avoid Seg Faults //
    if (ellipse.ptAxisMj2.y > imgOut.rows || ellipse.ptAxisMj2.y < 0)
        ellipse.ptAxisMj2.y = 0;

    if (ellipse.ptAxisMj2.x > imgOut.cols || ellipse.ptAxisMj2.x < 0)
        ellipse.ptAxisMj2.x = 0;

    if (ellipse.ptAxisMj1.y > imgOut.rows || ellipse.ptAxisMj1.y < 0)
        ellipse.ptAxisMj1.y = 0;

    if (ellipse.ptAxisMj1.x > imgOut.cols || ellipse.ptAxisMj1.x < 0)
        ellipse.ptAxisMj1.x = 0;


    imgOut.at<cv::Vec3b>(ellipse.ptAxisMj1)[1] = 255; imgOut.at<cv::Vec3b>(ellipse.ptAxisMj1)[2] = 255;
    imgOut.at<cv::Vec3b>(ellipse.ptAxisMj2)[1] = 255; imgOut.at<cv::Vec3b>(ellipse.ptAxisMj2)[2] = 255;
    //cv::circle(img_colour,ptxy1,1,CV_RGB(0,255,255),1);
    //cv::circle(,ptxy2,1,CV_RGB(0,255,255),1);
    //Debug Mark As Good Pair
    //imgDebug.at<uchar>(ellipse.ptAxisMj1) = 255;
    //imgDebug.at<uchar>(ellipse.ptAxisMj2) = 255;


}

int deleteUsedEdges( )
{

//    ///Step 12 - Remove the points from the image Before Restarting
//    for (std::vector<tEllipsoidEdges::iterator>::iterator itd = vedgePoints_trial.begin(); itd !=vedgePoints_trial.end(); )
//    {
//        tEllipsoidEdge* pEdge = &(*(*itd)); //Pickout Stored Iterator Pointers to Main list
//        //If this edge Is on The winning Ellipse's Minor Axis - Then Its been Used /Remove
//        if (pEdge->minorAxisLength == idx)
//        {
//            imgDebug.at<uchar>(pEdge->ptEdge) = 5; //Debug

//            pEdge->ptEdge.x = 0;
//            pEdge->ptEdge.y = 0;
//            itd = vedgePoints_trial.erase(itd);
//        }else {
//            ++itd;
//        }
//    } //Loop Through Used Points
//    //Invalidate Pair of Points
//    ptxy1.x = 0; ptxy1.y = 0;
//    ptxy2.x = 0; ptxy2.y = 0;

}


//Operator for Priority Ordering
bool operator<(const tDetectedEllipsoid& a,const tDetectedEllipsoid& b) {
  return a.fitscore < b.fitscore; //Max Heap
}

///
/// \brief detectEllipse Implements the Efficient ellipse Detection Algorithm -
/// \param vedgePoints_all
/// \param qEllipsoids
/// \notes The min Votes Threszhold is not fixed but continuously adapted to be just below the highest, see  gi_VotesEllipseThres
/// \return
///
/// \todo check image bounds
int detectEllipse(tEllipsoidEdges& vedgePoints_all, std::priority_queue<tDetectedEllipsoid>& qEllipsoids)
{
    const int minEllipseMajor   = gi_minEllipseMajor;
    const int maxEllipseMajor   = gi_maxEllipseMajor;
    const int minMinorEllipse   = gi_minEllipseMajor/1.8;
    int thresMinVotes     = gi_VotesEllipseThres;

    const int accLength = vedgePoints_all.size();
    int accumulator[accLength]; //The Score Holding (Histogram ) Array - Each index is a Minor Axis Length
    double HighestVotes = 0;
    double Highest2dVotes = 0;

    if (accLength < 3)
        return 0;

    std::vector<tEllipsoidEdges::iterator> vedgePoints_trial; //Containts edge points of specific ellipsoid trial
    vedgePoints_trial.reserve(10);

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator


    /// Begin Ellipsoid Detection ///
    memset(accumulator,0,sizeof(int)*(accLength)); //Reset Accumulator MAtrix
    //std::clog << "== Start === "  << std::endl;
    ///Loop through All Edge points (3)
    for (tEllipsoidEdges::iterator it1 = vedgePoints_all.begin();it1 != vedgePoints_all.end();++it1)
    {
        cv::Point2f ptxy1 = (*it1).ptEdge;
        if (ptxy1.x == 0 && ptxy1.y == 0)
            continue ; //point has been deleted
        if (ptxy1.x == 1 && ptxy1.y == 1)
            continue ; //point has been deleted

        cv::Point2f ptxy2;
        ///(4)

    /// Random Pair Formation //
        //Copy List Of Edges over and Randomize
        tEllipsoidEdges vedgePoints_pair = vedgePoints_all;
        std::uniform_int_distribution<> distr(1, vedgePoints_pair.size()-1); // define the range

        while (vedgePoints_pair.size() > 0)
        {
            tEllipsoidEdges::iterator it2 = vedgePoints_pair.begin();
            it2 += distr(eng);
            ptxy2 = (*it2).ptEdge;
            it2 = vedgePoints_pair.erase(it2);
    ////End of Random Pair //
//        for (tEllipsoidEdges::iterator it2 = vedgePoints_all.begin();it2 != vedgePoints_all.end(); ++it2 ) {
//            ptxy2 = (*it2).ptEdge;

            if (ptxy2.x == 0 && ptxy2.y == 0)
                continue ; //point has been deleted
            //if (ptxy2.x == 1 && ptxy2.y == 1)
//                continue ; //point has been deleted


            double d = cv::norm(ptxy2-ptxy1);
            if (d < minEllipseMajor || d > maxEllipseMajor)
                continue;

            //Use Eqns 1-4 and calculate Ellipse params
            cv::Point2f ptxy0;
            ptxy0.x = (ptxy2.x + ptxy1.x )/2.0;  //--(1)
            ptxy0.y = (ptxy2.y + ptxy1.y)/2.0;  //--(2)
            double a = d/2.0; //[(x 2 – x 1 )^2 + (y 2 – y 1 )^2 ] /2 //--(3) a the half-length of the major axis

            double alpha = atan2(ptxy2.y - ptxy1.y,ptxy2.x - ptxy1.x);//atan [(y 2 – y 1 )/(x 2 – x 1 )] //--(4) α the orientation of the ellipse

            //double dCntrLScore = round(cv::norm(ptxy0-ptLEyeMid));
            //double dCntrRScore = round(cv::norm(ptxy0-ptREyeMid));
            //double dCntrScore = std::min(dCntrLScore,dCntrRScore);


            ///Step (6) - 3rd Pixel;
            vedgePoints_trial.clear();
            for (tEllipsoidEdges::iterator it3 = vedgePoints_all.begin();it3 != vedgePoints_all.end(); ++it3 )
            {

                cv::Point2f ptxy3 = (*it3).ptEdge;

                if (ptxy3.x == 0 && ptxy3.y == 0)
                    continue ; //point has been deleted

                double d = (cv::norm(ptxy0-ptxy3));
                //Measure Distance From Centre of Eyes (Located at centre of img frame)

                double dd = d*d;

                if (d >= a || d < minMinorEllipse) //Candidate 3rd point of minor axis distance needs to be less than alpha away
                    continue;

                //Calculate Minor Axis
                double aa = a*a;
                double f = cv::norm(ptxy2-ptxy3);
                double ff = f*f;
                //double c = cv::norm(ptxy1-ptxy3);

                ///Step 7 - Calc the length of minor axis
                double costau = ( aa + dd - ff)/(2.0*a*d);
                double coscostau = costau*costau;  //eqn (6)
// b = sqrt( (aSq * thirdPtDistsSq(K) .* sinTauSq) ./ (aSq - thirdPtDistsSq(K) .* cosTau.^2 + eps) );
                double bb = aa*dd*(1.0-coscostau)/(aa - dd * coscostau + 0.00001); //(5)
                int b = std::round((sqrt(bb)));
                ///Step 8
                if (b > 1)
                {
                    //Make A "weighted" Band Of width 3
                    accumulator[b-1]+=1;
                    accumulator[b]  +=10; //increment x10 accumulator for this minor Axis = imgIn.at<uchar>(ptxy3)
                    accumulator[b+1]+=1; //increment x10 accumulator for this minor Axis = imgIn.at<uchar>(ptxy3)
///                 Add Intensity Density In the scoring - Eyes Are brighter Than Other features of the head
//                    double ellArea = M_PI*b*a;
//                    int iellArea = 0;
                    //Foci
//                    cv::Point2f focA,focB;
//                    //Take Direction Pointed By Major Axis
//                    focA = ptxy0+((ptxy1-ptxy0)/cv::norm(ptxy1-ptxy0))*sqrt(aa-bb);
//                    focB = ptxy0+((ptxy2-ptxy0)/cv::norm(ptxy2-ptxy0))*sqrt(aa-bb);
//                    int iIntensityTot = 0;
//                    //Divide Total Intensity By Area
//                    //Calc Pixel Density Find Pixels Inside Ellipse
//                    for(int i=0; i<imgIn.rows; i++)
//                        for(int j=0; j<imgIn.cols; j++)
//                        {    cv::Point2f pt(j,i); //x,y
//                            //Point Is inside Ellipse (Sum Of distances from Foci is less than Major axis)
//                            if ((cv::norm(focA-pt) + cv::norm(focB-pt)) < 2*a )
//                            {
//                                iIntensityTot += imgIn.at<uchar>(pt);
//                                iellArea ++;
//                            }
//                        }
//                    //Calc Normalized intensity Intensity Density and add to score
//                    double idensityScore = 1000.0*((double)iIntensityTot/iellArea)/(255.0*(double)iellArea);
//                    accumulator[b]+= idensityScore;

                    //Add Point to tracked List
                    it3->minorAxisLength = b;
                    vedgePoints_trial.push_back(it3); //Store Pointer To Point
                }

            ///Step 9 Loop Until All Pixels 3rd are computed for this pair of pixes
            }

            ///Step 10 //Find Max In accumulator array. The related length is the possible length of minor axis for assumed ellipse.
            double dvotesMax;
            int idx  = getMax(accumulator,accLength,dvotesMax);

            //Detect If Ellipse Is found /
            //idx Is the size of the minor axis
            if (dvotesMax > thresMinVotes) //Found ellipse
            {
                ///Step 11 Output Ellipse Parameters
                //cv::RotatedRect ellipse(ptxy0,ptxy1,ptxy2);
                //alpha += M_PI/2.0;
                cv::RotatedRect r(ptxy0,cv::Size2f(2*a,2*idx), alpha*(180.0/M_PI));
                tDetectedEllipsoid ellipse(ptxy0,ptxy1,ptxy2,dvotesMax,r);
                //vellipses.push_back(ellipse);
                qEllipsoids.push(ellipse); //Automatically Sorted

                ///Step 12 - Remove the points from the image Before Restarting
                for (std::vector<tEllipsoidEdges::iterator>::iterator itd = vedgePoints_trial.begin(); itd !=vedgePoints_trial.end(); )
                {
                    tEllipsoidEdge* pEdge = &(*(*itd)); //Pickout Stored Iterator Pointers to Main list
                    //If this edge Is on The winning Ellipse's Minor Axis - Then Its been Used /Remove
                    if (abs(pEdge->minorAxisLength - idx) == 0  ) //Delete The bin || pEdge->minorAxisLength == idx-1 || pEdge->minorAxisLength == idx-1
                    {
                        //imgDebug.at<uchar>(pEdge->ptEdge) = 200; //Debug - Show Used

                        pEdge->ptEdge.x = 0;
                        pEdge->ptEdge.y = 0;
                        itd = vedgePoints_trial.erase(itd);
                    }else {
                        ++itd;
                    }
                } //Loop Through Used Points

                //Invalidate  2nd Point of pair before moving to the next
                ptxy2.x = 0; ptxy2.y = 0;

//                if (r.boundingRect().contains(ptLEyeMid) )
//                   ptLEyeMid.x = 0; ptLEyeMid.y = 0;

//                if (r.boundingRect().contains(ptREyeMid) )
//                   ptLEyeMid.x = 0; ptLEyeMid.y = 0;


            }else {//Mark As Dull Pair
                //imgDebug.at<uchar>(ptxy1) = 55;
                //imgDebug.at<uchar>(ptxy2) = 55;
            }

            //Find Max Votes - Used to Re-adjust Threshold
            if (HighestVotes < dvotesMax)
            {
                Highest2dVotes  = HighestVotes;
                HighestVotes = dvotesMax;
                 //std::cout << "mxVot:" << HighestVotes << std::endl;
            }



        //it2 = vedgePoints.erase(it2);
        //
          //  it2->x = 0; it2->y = 0; //Delete Point


        ///Step 13 - Clear Accumulator
        memset(accumulator,0,sizeof(int)*(accLength)); //Reset Accumulator MAtrix

        } //Loop through each 2nd point in pair

//        cv::waitKey(1);
        ptxy1.x = 0; ptxy1.y = 0; //Invalidate pt1
        //it1 = vedgePoints.erase(it1);
        //it1->x = 0; it1->y = 0; //Delete Point
    } //Loop through all  point as 1st point pair (Prob: pairs can be repeated)


    gi_VotesEllipseThres = thresMinVotes = 0.95*Highest2dVotes;//Adapt Threshold To Best Score
//    std::clog << "ThresVot:" << gi_VotesEllipseThres << std::endl;


}

///
/// \brief getEyeSegThreshold Samples the N most intense Pixels in an arc below the estimated position of the eyes given the
/// upsampled head image
/// \param pimgIn //Upsampled Grey Scale HEad Image
/// \param ptcenter //Center Of Head Image around which to estimate Eye Position
/// \param ellipseSample_pts //Holds the Drawn Arc Points around the last spine Point
/// \param minVal - The min Intensity Value Sampled
/// \param maxVal - The min Intensity Value Sampled
/// \return Grey threshold for Eye Segmentation
///
int getEyeSegThreshold(cv::Mat& pimgIn,cv::Point2f ptcenter,std::vector<cv::Point>& ellipseSample_pts,int& minVal,int& maxVal)
{
        const int isampleN = 10;
        const int voffset = giHeadIsolationMaskVOffset+1;

        int iThresEyeSeg = 0;
        minVal = 255;
        maxVal = 0;

        //std::vector<cv::Point> ellipse_pts;
        std::priority_queue<int,std::vector<int>> eyeSegMaxHeap;


        //Construct Elliptical Circle around last Spine Point - of Radius step_size
        //Crash Here Stack
        cv::ellipse2Poly(ptcenter, cv::Size(voffset/2,voffset*0.9), 0, 185,345 , 1, ellipseSample_pts);
        for (int i=0;i<ellipseSample_pts.size();i++)
        {
            //iThresEyeSeg += imgUpsampled_gray.at<uchar>(ellipse_pts[i]);
            ellipseSample_pts[i].x = std::max(1,std::min(pimgIn.cols,ellipseSample_pts[i].x));
            ellipseSample_pts[i].y = std::max(1,std::min(pimgIn.rows,ellipseSample_pts[i].y));

            assert(ellipseSample_pts[i].x >= 0 && ellipseSample_pts[i].x <= pimgIn.cols);
            assert(ellipseSample_pts[i].y >= 0 && ellipseSample_pts[i].y <= pimgIn.rows);
            uchar val = pimgIn.at<uchar>(ellipseSample_pts[i]);
            eyeSegMaxHeap.push(val);

            if (val < minVal && val > 0)
                minVal = val;

            if (val > maxVal)
                maxVal = val;
        }

        for (int i=0;i<isampleN;i++)
        {//Withdraw To N values
            iThresEyeSeg  += eyeSegMaxHeap.top();
            eyeSegMaxHeap.pop();
        }

        //Add the Manual Entry And Divide to Get Mean Value
        iThresEyeSeg = (iThresEyeSeg+gthresEyeSeg)/(isampleN+1);




    return std::min(std::max(3,iThresEyeSeg),255);
}

///
/// \brief detectEllipses - Upsamples image Detects Eyes - Used oN Head Isolated Image
/// \param pimgIn
/// \param vellipses
/// \param outHeadFrameMonitor Image TO report Back Underlying Process of segmentation / The Edge Detection
/// \param outHeadFrameProc The Head Image with the Ellipses drawn
/// \return
///
int detectEllipses(cv::Mat& pimgIn,tEllipsoids& vellipses,cv::Mat& outHeadFrameMonitor,cv::Mat& outHeadFrameProc)
{

     cv::Mat imgFishHead_Lapl;

    int ret = 0;//Return Value Is the Count Of Ellipses Detected (Eyes)
    //assert(pimgIn.cols == imgEdge.cols && pimgIn.rows == imgEdge.rows);
    ///Keep Image processing Arrays Static to avoid memory Alloc On Each Run
    //cv::Mat img_contour;
    assert(pimgIn.rows > 0 && pimgIn.cols > 0);
    //cv::Mat imgEdge_dbg;


    std::priority_queue<tDetectedEllipsoid> qEllipsoids;
    std::vector<std::vector<cv::Point> > contours_canny;
    std::vector<cv::Vec4i> hierarchy_canny; //Contour Relationships  [Next, Previous, First_Child, Parent]
    std::vector<cv::Point> vEyeSegSamplePoints;


    tDetectedEllipsoid lEll,rEll;
    std::vector<cv::Point> vt;
    std::vector<cv::RotatedRect> ve;
    std::vector<cv::Point> vLEyeHull; //Left Eye
    std::vector<cv::Point> vREyeHull; //Left Eye

    cv::Point2f ptLEyeMid,ptREyeMid;


    //Upsamples an image which causes blur/interpolation it.
    cv::pyrUp(pimgIn, imgUpsampled_gray, cv::Size(pimgIn.cols*2,pimgIn.rows*2));

    int lengthLine = 13;
    cv::Point2f ptcentre(imgUpsampled_gray.cols/2,imgUpsampled_gray.rows/3+7);

    /*ptLEyeMid.x = ptcentre.x-lengthLine;
    ptLEyeMid.y = ptcentre.y/2; //y=0 is the top left corner
    ptREyeMid.x = ptcentre.x + lengthLine; //ptcentre.x+lengthLine;
    ptREyeMid.y = ptcentre.y/2; //y=0 is the top left corner *cos((angleDeg-90)*(M_PI/180.0))

    ptLEyeMid.x = std::max(1,std::min(imgUpsampled_gray.cols,(int)ptLEyeMid.x));
    ptLEyeMid.y = std::max(1,std::min(imgUpsampled_gray.rows,(int)ptLEyeMid.y));

    ptREyeMid.x = std::max(1,std::min(imgUpsampled_gray.cols,(int)ptREyeMid.x));
    ptREyeMid.y = std::max(1,std::min(imgUpsampled_gray.rows,(int)ptREyeMid.y));
    */
    cv::GaussianBlur(imgUpsampled_gray,imgUpsampled_gray,cv::Size(3,3),3,3);


    // Locate Eye Points //
    double minVal,maxVal;
    cv::Point ptMax,ptMin;
    ///COVER Right Eye - Find Left EYE //
    cv::Rect rRightMask(imgUpsampled_gray.cols/2,0,imgUpsampled_gray.cols,imgUpsampled_gray.rows);
    cv::Mat imgEyeDiscover = imgUpsampled_gray.clone();
    // Make Body Mask For bOth //
    cv::circle(imgEyeDiscover,cv::Point(imgUpsampled_gray.cols/2,imgUpsampled_gray.rows),giHeadIsolationMaskVOffset,CV_RGB(0,250,50),CV_FILLED); //Mask Body

    ///COVER Right Eye - Find Left EYE //
    cv::Mat imgEyeCover = imgEyeDiscover.clone();
    cv::rectangle(imgEyeCover,rRightMask,cv::Scalar(0),-1);
    //Find Eye On Left Side
    cv::minMaxLoc(imgEyeCover,&minVal,&maxVal,&ptMin,&ptMax);
    ptLEyeMid = ptMax;

    ///COVER Left Eye - Find RIGHT EYE //
    imgEyeCover = imgEyeDiscover.clone();
    cv::Rect rLeftMask(0,0,imgUpsampled_gray.cols/2,imgUpsampled_gray.rows);
    cv::rectangle(imgEyeCover,rLeftMask,cv::Scalar(0),-1);
    cv::minMaxLoc(imgEyeCover,&minVal,&maxVal,&ptMin,&ptMax); //Find Centre of RIght Eye
    ptREyeMid = ptMax;


    /// Make Arc from Which to get Sample Points For Eye Segmentation
    int ilFloodRange,iuFloodRange;
    //Fill BG //
//    cv::floodFill(imgUpsampled_gray, cv::Point(0,0), cv::Scalar(0),0,cv::Scalar(3),cv::Scalar(5));

    //equalize the histogram
     //cv::Mat hist_equalized_image;

    /// Equalize Histogram to Enhance Contrast
    if (bUseHistEqualization)
     cv::equalizeHist(imgUpsampled_gray, imgUpsampled_gray);

    /// Estimate Eye Segmentation threshold from sample points in Image
    int iThresEyeSeg = getEyeSegThreshold(imgUpsampled_gray,ptcentre,vEyeSegSamplePoints,ilFloodRange,iuFloodRange);

    //int ilFloodSeed=imgUpsampled_gray.at<uchar>(ptLEyeMid)+1;
    //int irFloodSeed=imgUpsampled_gray.at<uchar>(ptREyeMid)+1;
    //int stepL = (ilFloodSeed - ilFloodRange)/gi_minEllipseMajor;
    //int stepR = (irFloodSeed - ilFloodRange)/gi_minEllipseMajor;
    //Assist by Filling Holes IN Eye Shape - Use Fill
    //cv::floodFill(imgUpsampled_gray, ptLEyeMid, cv::Scalar(iThresEyeSeg+1),0,cv::Scalar(abs(2*(ilFloodRange+1)-ilFloodSeed)),cv::Scalar(abs(iuFloodRange-ilFloodSeed)),CV_FLOODFILL_FIXED_RANGE);
    //cv::floodFill(imgUpsampled_gray, ptREyeMid, cv::Scalar(iThresEyeSeg+1),0,cv::Scalar(abs(2*(ilFloodRange+1)-irFloodSeed)),cv::Scalar(abs(iuFloodRange-irFloodSeed)),CV_FLOODFILL_FIXED_RANGE);
//    cv::floodFill(imgUpsampled_gray, ptLEyeMid, cv::Scalar(iThresEyeSeg+1),0,cv::Scalar(abs(ilFloodSeed)/10),cv::Scalar(abs(ilFloodSeed)/10),CV_FLOODFILL_FIXED_RANGE);
//    cv::floodFill(imgUpsampled_gray, ptREyeMid, cv::Scalar(iThresEyeSeg+1),0,cv::Scalar(abs(irFloodSeed)/10),cv::Scalar(abs(irFloodSeed)/10),CV_FLOODFILL_FIXED_RANGE);



    //cv::floodFill(imgUpsampled_gray, ptLEyeMid, cv::Scalar(iThresEyeSeg+1),0,cv::Scalar(stepL),cv::Scalar(stepL));
    //cv::floodFill(imgUpsampled_gray, ptREyeMid, cv::Scalar(iThresEyeSeg+1),0,cv::Scalar(stepR),cv::Scalar(stepR));


    //Show Eye Points to User //
    cv::circle(imgUpsampled_gray,ptREyeMid,2,cv::Scalar(255),1);
    cv::circle(imgUpsampled_gray,ptLEyeMid,2,cv::Scalar(255),1);

    /// Show Masks with nominal width//
    cv::line(imgUpsampled_gray,ptcentre,cv::Point(imgUpsampled_gray.cols/2,0),CV_RGB(0,250,50),1);//Split Eyes iEyeMaskSepWidth
    cv::circle(imgUpsampled_gray,cv::Point(imgUpsampled_gray.cols/2,imgUpsampled_gray.rows),giHeadIsolationMaskVOffset,CV_RGB(0,250,50),1); //Mask Body


    // Do Thresholding Of Masked Image to Obtain Segmented Eyes //
    cv::threshold(imgUpsampled_gray, imgIn_thres,iThresEyeSeg,255,cv::THRESH_BINARY); // Log Threshold Image + cv::THRESH_OTSU

    //Try Laplacian CV_8U
    cv::Laplacian(imgIn_thres,imgFishHead_Lapl,imgIn_thres.type(),1);

    imgFishHead_Lapl.copyTo(imgEdge_local);


    //Separate Eyes Mask
    //Add Thick Mid line to erase inner Eye Edges and artefacts
    cv::line(imgEdge_local,ptcentre,cv::Point(imgUpsampled_gray.cols/2,0),CV_RGB(0,0,0),iEyeMaskSepWidth);//Split Eyes
    cv::circle(imgEdge_local,cv::Point(imgUpsampled_gray.cols/2,imgUpsampled_gray.rows),giHeadIsolationMaskVOffset,CV_RGB(0,0,0),-1); //Mask Body

    //cv::adaptiveThreshold(imgIn, imgIn_thres, 255,cv::ADAPTIVE_THRESH_GAUSSIAN_C,cv::THRESH_BINARY,2*(imgIn.cols/2)-1,10 ); // Log Threshold Image + cv::THRESH_OTSU

    outHeadFrameMonitor = imgEdge_local.clone();
    //imgIn_thres.copyTo(outHeadFrameMonitor);

    //cv::erode(imgIn_thres,imgIn_thres,kernelOpen,cv::Point(-1,-1),1);
    cv::morphologyEx(imgIn_thres,imgIn_thres, cv::MORPH_OPEN, kernelOpenfish,cv::Point(-1,-1),1); //Break Connections
    //cv::morphologyEx(imgEdge_local,imgEdge_local, cv::MORPH_CLOSE, kernelOpenfish,cv::Point(-1,-1),1);
    //cv::erode(imgIn_thres,imgIn_thres,kernelOpen,cv::Point(-1,-1),3);

    cv::findContours(imgEdge_local, contours_canny,hierarchy_canny, cv::RETR_CCOMP,cv::CHAIN_APPROX_SIMPLE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE


    //Empty List
    vellipses.clear();
    vellipses.shrink_to_fit();

    //std::vector<std::vector<cv::Point>> vEyes;

    //Find Parent Contour
    int iLEye = findMatchingContour(contours_canny,hierarchy_canny,ptLEyeMid,2,vt,ve);
    int iREye = findMatchingContour(contours_canny,hierarchy_canny,ptREyeMid,2,vt,ve);

    cv::RotatedRect rcLEye,rcREye;
    //Make Debug Img

    cv::cvtColor( imgUpsampled_gray,img_colour, cv::COLOR_GRAY2RGB);
    //cv::cvtColor( imgUpsampled_gray,img_contour, cv::COLOR_GRAY2RGB);

    //for( size_t i = 0; i< contours_canny.size(); i++ )

    if (iLEye != -1) //If Contour Is found
    {

        cv::convexHull( cv::Mat(contours_canny[iLEye]), vLEyeHull, false );
        if (vLEyeHull.size() > 4)
        {
            //vEyes.push_back(vLEyeHull);
            rcLEye =  cv::fitEllipse(vLEyeHull);
            tDetectedEllipsoid dEll(rcLEye,100);
            lEll.fitscore       = dEll.fitscore;
            lEll.rectEllipse    = dEll.rectEllipse;
            qEllipsoids.push(dEll);
            //cv::drawContours( img_contour, vEyes, 0, CV_RGB(10,205,10),1);
            //cv::drawContours( imgEdge_local, vEyes, 0, CV_RGB(255,255,255),1);
        }
        //getEdgePoints(contours_canny.at(iLEye),vedgePoints_all);
    }

    //Check if Std Ellipse Finding Worked / Otherwise Try Method 2
    int mjAxis1 = std::max(lEll.rectEllipse.size.width, lEll.rectEllipse.size.height);
    if (mjAxis1 < gi_minEllipseMajor || mjAxis1 > gi_maxEllipseMajor || lEll.fitscore < 10)
    {
        //Empty
        while (qEllipsoids.size() > 0)
            qEllipsoids.pop(); //Empty All Other Candidates


        tEllipsoidEdges vedgePoints_all; //All edge points from Image Of EDge detection
        vedgePoints_all.clear();

        //If Contour Finding Fails Then Take Raw Edge points and mask L/R half of image
        try
        {
            //imgEdge_local = cv::Mat::zeros(imgUpsampled_gray.rows,imgUpsampled_gray.cols,CV_8UC1);
           // cv::Canny( imgIn_thres, imgEdge_local, gi_CannyThresSmall,gi_CannyThres  );
        }
        catch (char* e)
        {
            pwindow_main->LogEvent("Error detectEllipses  L Eye Canny processing ");
            std::cerr << e << std::endl;

        }
        outHeadFrameMonitor = imgEdge_local.clone();
        //COVER Right Eye
        cv::Rect r(imgEdge_local.cols/2,0,imgIn_thres.cols,imgIn_thres.rows);
        //imgEdge.copyTo(imgEdge_local);
        cv::rectangle(imgEdge_local,r,cv::Scalar(0),-1);
        //cv::imshow("REyeCover",imgEdge_local);

        getEdgePoints(imgEdge_local,vedgePoints_all);
        detectEllipse(vedgePoints_all,qEllipsoids); //Run Ellipsoid fitting Algorithm
        //imgEdge_local.copyTo(imgEdge_dbg);
        if (qEllipsoids.size() == 0 )
        //    qDebug() << " L Eye Backup Ellipse Detection found score: " << qEllipsoids.top().fitscore;
        //else
            qDebug() << " L Eye Backup Ellipse Failed";

    }


    ///Store Left Eye And Draw Detected Ellipsoid
    if (qEllipsoids.size() > 0)
    {
        //Pick Best Match For this Eye from to of Priority List
        lEll = qEllipsoids.top();
        //Draw it
        drawEllipse(img_colour,lEll);

        //Store it To Output Vector
        vellipses.push_back(lEll);
        ret++;
        cv::Point2f featurePnts[4];
        lEll.rectEllipse.points(featurePnts);

        ///Draw Left Eye Rectangle
        for (int j=0; j<4;j++) //Rectangle Eye
               cv::line(img_colour,featurePnts[j],featurePnts[(j+1)%4] ,CV_RGB(10,10,130),1);
        //Draw Line
        cv::line(img_colour,lEll.ptAxisMj1,lEll.ptAxisMj2,CV_RGB(10,10,130),1);
        //Empty
        while (qEllipsoids.size() > 0)
            qEllipsoids.pop(); //Empty All Other Candidates
    }

   ///// End oF LEft Eye Trace ///



    /// - RIGHT EYE - Reset And Redraw - ////
    if (iREye != -1)
    {
        //imgEdge_local = cv::Mat::zeros(imgUpsampled_gray.rows,imgUpsampled_gray.cols,CV_8UC1);
        //Crash Can Occur Here When the Fish Is Rushing too fast -
        // malloc(): memory corruption (fast):
        cv::convexHull( cv::Mat(contours_canny[iREye]), vREyeHull, false );

        if (vREyeHull.size() > 4)
        {
            //vEyes.push_back(vREyeHull);
            rcREye =  cv::fitEllipse(vREyeHull);
            tDetectedEllipsoid dEll(rcREye,100);
            rEll.fitscore       = dEll.fitscore;
            rEll.rectEllipse    = dEll.rectEllipse;

            qEllipsoids.push(dEll); //Index 1 / Right Eye
            //cv::drawContours( img_contour, vEyes, vEyes.size()-1, CV_RGB(10,05,210),1);
            //cv::drawContours( imgEdge_local, vEyes,vEyes.size()-1, CV_RGB(255,255,255),1);
        }
        //getEdgePoints(contours_canny.at(iREye),vedgePoints_all);
    }

    //Check
    int mjAxis2 = std::max(rEll.rectEllipse.size.width, rEll.rectEllipse.size.height);
    ///Check If 1st Method Failed And Run Backup Method If 1st Failed
    if (mjAxis2 < gi_minEllipseMajor || mjAxis2 > gi_maxEllipseMajor || rEll.fitscore < 10)
    {
        //Empty
        while (qEllipsoids.size() > 0)
            qEllipsoids.pop(); //Empty All Other Candidates

        tEllipsoidEdges vedgePoints_all; //All edge points from Image Of EDge detection
        vedgePoints_all.clear();

        //If Contour Finding Fails Then Take Raw Edge points and *MASK* L/R half of image
        try
        {
            //imgEdge_local = cv::Mat::zeros(imgUpsampled_gray.rows,imgUpsampled_gray.cols,CV_8UC1);
            //cv::Canny( imgIn_thres, imgEdge_local, gi_CannyThresSmall,gi_CannyThres  );

        }
        catch (char* e)
        {
            pwindow_main->LogEvent("Error in R Eye Canny processing ");
            std::cerr << e << std::endl;
        }

        outHeadFrameMonitor = imgEdge_local.clone();
        //Cover LEFT Eye Edges
        cv::Rect r(0,0,imgEdge_local.cols/2,imgEdge_local.rows);
        //imgEdge.copyTo(imgEdge_local);
        cv::rectangle(imgEdge_local,r,cv::Scalar(0),-1);
        //cv::imshow("LEyeCover",imgEdge_local);

        getEdgePoints(imgEdge_local,vedgePoints_all);
        detectEllipse(vedgePoints_all,qEllipsoids);
        if (qEllipsoids.size() == 0 )
//            qDebug() << " R Eye Backup Ellipse Failed";
            pwindow_main->LogEvent("R Eye Backup Ellipse Failed ");

            //qDebug() << " R Eye Backup Ellipse Detection found score: " << qEllipsoids.top().fitscore;
        //else
    }

    // Check If Found and Draw R Eye //
    if (qEllipsoids.size() > 0)
    {
        rEll = qEllipsoids.top();
        drawEllipse(img_colour,rEll);
        //Store it To Output Vector
        vellipses.push_back(rEll);
        ret++;

        cv::Point2f featurePnts[4];
        rEll.rectEllipse.points(featurePnts);


        ///Draw Left Eye Rectangle
        for (int j=0; j<4;j++) //Rectangle Eye
               cv::line(img_colour,featurePnts[j],featurePnts[(j+1)%4] ,CV_RGB(130,10,10),1);

        cv::line(img_colour,rEll.ptAxisMj1,rEll.ptAxisMj2 ,CV_RGB(130,10,10),1);

        while (qEllipsoids.size() > 0)
            qEllipsoids.pop(); //Empty All Other Candidates
    }

    /// L And R Eyes Detection is Done- Check Results //

    // Evaluate Detection - Use  Limit Checks on Eye Characteristics ////
    int area1 = (int)lEll.rectEllipse.size.width*lEll.rectEllipse.size.height;
    int area2 = (int)rEll.rectEllipse.size.width*rEll.rectEllipse.size.height;

    ///Check L Eye Again
    mjAxis1 = std::max(lEll.rectEllipse.size.width, lEll.rectEllipse.size.height);
    if (mjAxis1 < gi_minEllipseMajor || mjAxis1 > gi_maxEllipseMajor || lEll.fitscore < 10)
    {
        ret = 0; //SOme Detection Error - Ask To Change Threshold
        //qDebug() << "L eye bound error  MjAxis:" << mjAxis1;
        pwindow_main->LogEvent("L eye MjAxis value" + QString::number(mjAxis1) + " is out of bounds   ");
    }

    /// Check R Eye Again //
    mjAxis2 = std::max(rEll.rectEllipse.size.width, rEll.rectEllipse.size.height);
    if (mjAxis2 < gi_minEllipseMajor || mjAxis2 > gi_maxEllipseMajor || rEll.fitscore < 10)
    {
        ret = 0;
        //qDebug() << "R eye bound error MjAxis :" << mjAxis1;

        pwindow_main->LogEvent("R eye MjAxis value: " + QString::number(mjAxis2) + " is out of bounds   ");
    }

    if (std::abs(area1 - area2) > (std::max(area2,area1)))
    {
        ret = 0; //SOme Detection Error - Ask To Change Threshold
        //qDebug() << "R-L eye Areas Diff too large " << std::abs(area1 - area2);
        pwindow_main->LogEvent("R-L eye Areas Diff too large " + QString::number(std::abs(area1 - area2)));
    }


///////// End of Checks //////////

   /// Debug //
   //cv::bitwise_or(imgEdge,imgEdge_dbg,imgEdge_dbg);
   // cv::imshow("Fish Edges ",imgEdge_dbg);
   // cv::imshow("Fish Edges h",imgEdge);
   //cv::imshow("Debug EllipseFit",imgDebug);
   //cv::imshow("Fish Threshold ",imgIn_thres);
   //cv::imshow("Fish CONTOUR ",img_contour);




   // Show Eye Anchor Points
    img_colour.at<cv::Vec3b>(ptLEyeMid)[0] = 0; img_colour.at<cv::Vec3b>(ptLEyeMid)[1] = 0;img_colour.at<cv::Vec3b>(ptLEyeMid)[2] = 205; //Blue
    img_colour.at<cv::Vec3b>(ptREyeMid)[0] = 0; img_colour.at<cv::Vec3b>(ptREyeMid)[1] = 0;img_colour.at<cv::Vec3b>(ptREyeMid)[2] = 205; //Blue

    // Show Eye Segmentation Arc Sample points
    for (int i=0;i<vEyeSegSamplePoints.size();i++)
    {
        img_colour.at<cv::Vec3b>(vEyeSegSamplePoints[i])[0] = 220;
        img_colour.at<cv::Vec3b>(vEyeSegSamplePoints[i])[1] = 220;
        img_colour.at<cv::Vec3b>(vEyeSegSamplePoints[i])[2] = 50;
    }


    //outHeadFrameProc = img_colour.clone(); //Make A deep Copy - Avoids Seg Faults with C
    img_colour.copyTo(outHeadFrameProc);
    // Memory Crash Here //
    //contours_canny.clear();
    //contours_canny.shrink_to_fit();

return ret;

} //End of DetectEllipses


///
/// \brief calculates and show_histogram and its 1st derivative calculated from Right to left
/// \param name of window to show hist.
/// \param image - source image on which to obtain histogram
///
void show_histogram(std::string const& name, cv::Mat1b const& image)
{
    // Set histogram bins count
    int bins = 256;
    int picoffset   =20;
    int histSize[] = {bins};
    // Set ranges for histogram bins
    float lranges[] = {0, 256};
    const float* ranges[] = {lranges};
    // create matrix for histogram
    cv::Mat hist,hist_smooth,hist_grad;
    int channels[] = {0};

    // create matrix for histogram visualization
    const int  hist_height = 256;
    cv::Mat3b hist_image = cv::Mat3b::zeros(hist_height*2+picoffset, bins);
    //Obtain Histogram
    cv::calcHist(&image, 1, channels, cv::Mat(), hist, 1, histSize, ranges, true, false);

    double max_val      = 0.0;
    double max_val_grad = 0.0;
    double min_val_grad = 0.0;
    cv::Point   max_idx_grad;
    cv::Point   min_idx_grad;

    max_val = 150;

    cv::blur(hist,hist_smooth,cv::Size(1,41));

    hist_grad = cv::Mat::zeros(hist.rows,hist.cols,hist.type());
    //max_idx_grad = new int(hist.dims);
    //min_idx_grad = new int(hist.dims);


    //cv::Sobel( hist_smooth, hist_grad, CV_8UC1, 1, 0, 3, 1, 0, cv::BORDER_DEFAULT );
    //cv::convertScaleAbs( hist_grad, hist_grad );

    //Compute Changes / Grad
    for(int b = bins-1; b >1; b--) {
        float const binVal = hist_smooth.at<float>(b);
        hist_grad.at<float>(b) = hist_smooth.at<float>(b) - hist_smooth.at<float>(b+1);
    }

    //Smooth The Histogram Gradient
    cv::blur(hist_grad,hist_grad,cv::Size(1,91));

    //Locate Max - For normalizing and Min for finding the seg. threshold
    cv::minMaxLoc(hist_grad, &min_val_grad, &max_val_grad, &min_idx_grad, &max_idx_grad);

    /// Find Eye segmentation Threshold - The distribution has a bimodality-
    /// Take a point close/above to the 1 -ve point in hist. starting from highest to lowest grey values

    //Plot Histogram
    for(int b = 0; b < bins; b++) {
        float const binVal = hist_smooth.at<float>(b);
        int   const height = cvRound(binVal*hist_height/max_val);
        cv::line
            ( hist_image
            , cv::Point(b, 2*hist_height-height + picoffset ), cv::Point(b, 2*hist_height)
            , cv::Scalar::all(165)
            );


        //Plot Hist Gradient
        float const sbinVal = hist_grad.at<float>(b);
        int  const sheight = cvRound(sbinVal*hist_height/max_val_grad);
        cv::line
            ( hist_image
            , cv::Point(b, hist_height-sheight + picoffset), cv::Point(b, hist_height)
            , cv::Scalar::all(255)
            );

        //Draw Set Points
        if (min_idx_grad.y == b)
        {
           cv::circle(hist_image,cv::Point(b, hist_height-sheight+ picoffset),6,CV_RGB(80,80,80),CV_FILLED);
           cv::putText(hist_image,"L",cv::Point(b, hist_height-sheight+ picoffset),CV_FONT_HERSHEY_PLAIN,0.7,CV_RGB(200,200,200),1);
        }

        if (max_idx_grad.y == b)
        {
           cv::circle(hist_image, cv::Point(b, hist_height-sheight+ picoffset),6,CV_RGB(100,100,100),CV_FILLED);
           cv::putText(hist_image,"H",cv::Point(b, hist_height-sheight+ picoffset),CV_FONT_HERSHEY_PLAIN,0.7,CV_RGB(200,200,200),1);
        }

        if (gthresEyeSeg == b)
        {
           cv::circle(hist_image, cv::Point(b, hist_height-sheight+ picoffset),8,CV_RGB(150,150,150),CV_FILLED);
           cv::putText(hist_image,"T",cv::Point(b, hist_height-sheight+ picoffset),CV_FONT_HERSHEY_PLAIN,0.7,CV_RGB(200,200,200),1);
        }

        //Adaptive Threshold
        //Find the first low above after the first High mode -
        //Get the 1st peak above this 1st low - Set it to be the thresh. The adaptive

    }

    cv::imshow(name, hist_image);

    hist_image.deallocate();
    hist_grad.deallocate();

}


