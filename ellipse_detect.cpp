///*
/// Implements Algorithm based on A NEW EFFICIENT ELLIPSE DETECTION METHOD,  2002
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
#include <larvatrack.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip> //for setprecision
#include <limits>
#include <string>
#include <random>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/core/utility.hpp"
//#include <opencv2/bgsegm.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>

extern int gi_CannyThresSmall;
extern int gi_CannyThres;
extern int gi_VotesEllipseThres;
extern int gi_minEllipseMajor;
extern int gi_maxEllipseMajor;
extern int g_BGthresh;
cv::Mat imgDebug;

extern cv::Mat kernelOpen;
extern cv::Mat frameDebugC;

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
               imgDebug.at<uchar>(pt) = 125;
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
               imgDebug.at<uchar>(contour[i]) = 155;
       }

}


void drawEllipse(cv::Mat imgOut,tDetectedEllipsoid ellipse)
{

    cv::ellipse(imgOut,ellipse.rectEllipse,CV_RGB(250,50,50),1,cv::LINE_8);
    cv::circle(imgOut,ellipse.rectEllipse.center,1,CV_RGB(0,0,255),1);
    imgOut.at<cv::Vec3b>(ellipse.ptAxisMj1)[1] = 255; imgOut.at<cv::Vec3b>(ellipse.ptAxisMj1)[2] = 255;
    imgOut.at<cv::Vec3b>(ellipse.ptAxisMj2)[1] = 255; imgOut.at<cv::Vec3b>(ellipse.ptAxisMj2)[2] = 255;
    //cv::circle(img_colour,ptxy1,1,CV_RGB(0,255,255),1);
    //cv::circle(img_colour,ptxy2,1,CV_RGB(0,255,255),1);
    //Debug Mark As Good Pair
    imgDebug.at<uchar>(ellipse.ptAxisMj1) = 255;
    imgDebug.at<uchar>(ellipse.ptAxisMj2) = 255;


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

int detectEllipse(tEllipsoidEdges& vedgePoints_all, std::priority_queue<tDetectedEllipsoid>& qEllipsoids)
{
    const int minEllipseMajor   = gi_minEllipseMajor;
    const int maxEllipseMajor   = gi_maxEllipseMajor;
    const int minMinorEllipse   = 1;
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
    std::cout << "== Start === "  << std::endl;
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
                {   //Make A Band Of width 2
                    //accumulator[b+1]++; //imgIn<uchar>.at(ptxy3) ; //increment accumulator for this minor Axis
                    //accumulator[b+1]++; //increment accumulator for this minor Axis
                    //accumulator[b]-= dCntrScore/4; //Add Points for being close to Eye centre

                    accumulator[b-1]+=1;
                    accumulator[b]  +=10; //increment x10 accumulator for this minor Axis = imgIn.at<uchar>(ptxy3)
                    accumulator[b+1]+=1; //increment x10 accumulator for this minor Axis = imgIn.at<uchar>(ptxy3)
                    //accumulator[b-1]+=10; //Add points to smaller ellipsoids so as to smooth out close edge point issues
                    //accumulator[b+1]+=10;
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
                cv::RotatedRect r(ptxy0,cv::Size2f(2*a,2*idx), alpha*(180/M_PI));
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
                        imgDebug.at<uchar>(pEdge->ptEdge) = 200; //Debug - Show Used

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
                imgDebug.at<uchar>(ptxy1) = 55;
                imgDebug.at<uchar>(ptxy2) = 55;
            }

            //Find Max Votes - Used to Re-adjust Threshold
            if (HighestVotes < dvotesMax)
            {
                Highest2dVotes  = HighestVotes;
                HighestVotes = dvotesMax;
                 std::cout << "mxVot:" << HighestVotes << std::endl;
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
    std::cout << "ThresVot:" << gi_VotesEllipseThres << std::endl;


}

int detectEllipses(cv::Mat& imgIn,cv::Mat imgEdge,cv::Mat& imgOut,int angleDeg,tEllipsoids& vellipses)
{
    assert(imgIn.cols == imgEdge.cols && imgIn.rows == imgEdge.rows);

    std::priority_queue<tDetectedEllipsoid> qEllipsoids;

    std::vector<std::vector<cv::Point> > contours_canny;
    std::vector<cv::Vec4i> hierarchy_canny; //Contour Relationships  [Next, Previous, First_Child, Parent]

    cv::Point2f ptLEyeMid,ptREyeMid;
    cv::Point2f ptcentre(imgIn.cols/2,imgIn.rows/2);
    int lengthLine = 4;
    ptLEyeMid.x = ptcentre.x-lengthLine;
    ptLEyeMid.y = ptcentre.y; //y=0 is the top left corner
    ptREyeMid.x = ptcentre.x + lengthLine; //ptcentre.x+lengthLine;
    ptREyeMid.y = ptcentre.y; //y=0 is the top left corner *cos((angleDeg-90)*(M_PI/180.0))


    cv::Mat img_colour,img_contour,imgIn_thres,imgEdge_local,imgEdge_dbg;
    //Debug
    imgDebug = cv::Mat::zeros(imgIn.rows,imgIn.cols,CV_8UC1);

    //cv::GaussianBlur(imgIn,img_blur,cv::Size(3,3),3,3);
    //cv::Laplacian(img_blur,img_edge,CV_8UC1,g_BGthresh);
    //Get Pixel Value Between Eyes
    int thresEyeSeg = imgIn.at<uchar>(ptcentre) + imgIn.at<uchar>(ptcentre.y-1,ptcentre.x) + imgIn.at<uchar>(ptcentre.y+1,ptcentre.x) + imgIn.at<uchar>(ptcentre.y+2,ptcentre.x);
    thresEyeSeg     += (int)imgIn.at<uchar>(imgIn.rows-3,ptcentre.x) + (int)imgIn.at<uchar>(imgIn.rows-1,ptcentre.x);
    thresEyeSeg     = thresEyeSeg/6+10;
    //Report
    std::cout << (int)imgIn.at<uchar>(ptcentre) << "+" << (int)imgIn.at<uchar>(ptcentre.y-1,ptcentre.x) << "+" <<(int) imgIn.at<uchar>(ptcentre.y+1,ptcentre.x) << "+" << (int)imgIn.at<uchar>(ptcentre.y+2,ptcentre.x);
    std::cout << "+" << (int)imgIn.at<uchar>(imgIn.rows-2,ptcentre.x) << "+" <<  (int)imgIn.at<uchar>(imgIn.rows-1,ptcentre.x) << " avg:" <<  thresEyeSeg << std::endl;

    cv::threshold(imgIn, imgIn_thres,thresEyeSeg,255,cv::THRESH_BINARY); // Log Threshold Image + cv::THRESH_OTSU
    //cv::adaptiveThreshold(imgIn, imgIn_thres, 255,cv::ADAPTIVE_THRESH_GAUSSIAN_C,cv::THRESH_BINARY,2*(imgIn.cols/2)-1,10 ); // Log Threshold Image + cv::THRESH_OTSU

    cv::erode(imgIn_thres,imgIn_thres,kernelOpen,cv::Point(-1,-1),1);
    cv::morphologyEx(imgIn_thres,imgIn_thres, cv::MORPH_OPEN, kernelOpen,cv::Point(-1,-1),2);
    //cv::morphologyEx(imgIn_thres,imgIn_thres, cv::MORPH_CLOSE, kernelOpen,cv::Point(-1,-1),1);
    //cv::erode(imgIn_thres,imgIn_thres,kernelOpen,cv::Point(-1,-1),3);

    cv::findContours(imgIn_thres, contours_canny,hierarchy_canny, cv::RETR_CCOMP,cv::CHAIN_APPROX_SIMPLE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE


    //Empty List
    vellipses.clear();
    tEllipsoidEdges vedgePoints_all; //All edge points from Image Of EDge detection

    vedgePoints_all.reserve(imgEdge.cols*imgEdge.rows/2);


    std::vector<cv::Point> vt;
    std::vector<cv::RotatedRect> ve;
    std::vector<std::vector<cv::Point>> vEyes;
    std::vector<cv::Point> vLEyeHull; //Left Eye
    std::vector<cv::Point> vREyeHull; //Left Eye
    //Find Parent Contour
    int iLEye = findMatchingContour(contours_canny,hierarchy_canny,ptLEyeMid,2,vt,ve);
    int iREye = findMatchingContour(contours_canny,hierarchy_canny,ptREyeMid,2,vt,ve);

    cv::RotatedRect rcLEye,rcREye;
    //Make Debug Img
    cv::cvtColor( imgIn,img_colour, cv::COLOR_GRAY2RGB);
    cv::cvtColor( imgIn,img_contour, cv::COLOR_GRAY2RGB);

    //for( size_t i = 0; i< contours_canny.size(); i++ )
    vedgePoints_all.clear();

    if (iLEye != -1)
    {
        imgEdge_local = cv::Mat::zeros(imgIn.rows,imgIn.cols,CV_8UC1);
        cv::convexHull( cv::Mat(contours_canny[iLEye]), vLEyeHull, false );

        vEyes.push_back(vLEyeHull);
        rcLEye =  cv::fitEllipse(vLEyeHull);
        cv::drawContours( img_contour, vEyes, 0, CV_RGB(10,205,10),1);
        cv::drawContours( imgEdge_local, vEyes, 0, CV_RGB(255,255,255),1);
        //getEdgePoints(contours_canny.at(iLEye),vedgePoints_all);
    }
    else {
        //If Contour Finding Fails Then Take Raw Edge points and mask L/R half of image
        cv::Canny( imgIn_thres, imgEdge_local, gi_CannyThresSmall,gi_CannyThres  );
        //Cover Right Eye
        cv::Rect r(ptcentre.x,ptcentre.y,imgEdge.cols/2-1,imgEdge.rows);
        //imgEdge.copyTo(imgEdge_local);
        cv::rectangle(imgEdge_local,r,cv::Scalar(0),-1);

    }

    getEdgePoints(imgEdge_local,vedgePoints_all);
    detectEllipse(vedgePoints_all,qEllipsoids);


    if (qEllipsoids.size() > 0)
    {
        tDetectedEllipsoid dEll = qEllipsoids.top();
        drawEllipse(img_colour,dEll);


        cv::Point2f featurePnts[4];
        rcLEye.points(featurePnts);
        ///Draw Left Eye Rectangle
        for (int j=0; j<4;j++) //Rectangle Eye
               cv::line(img_colour,featurePnts[j],featurePnts[(j+1)%4] ,CV_RGB(10,10,130),1);


        if (qEllipsoids.size() > 0)  qEllipsoids.pop();
        if (qEllipsoids.size() > 0) qEllipsoids.pop();
        if (qEllipsoids.size() > 0) qEllipsoids.pop();
        if (qEllipsoids.size() > 0) qEllipsoids.pop();
    }
   imgEdge_local.copyTo(imgEdge_dbg);
   ///End oF LEft Eye Trace //

    //Reset And Redraw - Now Right Eye
    if (iREye != -1)
    {
        imgEdge_local = cv::Mat::zeros(imgIn.rows,imgIn.cols,CV_8UC1);
        cv::convexHull( cv::Mat(contours_canny[iREye]), vREyeHull, false );
        vEyes.push_back(vREyeHull);
        rcREye =  cv::fitEllipse(vREyeHull);
        cv::drawContours( img_contour, vEyes, vEyes.size()-1, CV_RGB(10,05,210),1);
        cv::drawContours( imgEdge_local, vEyes,vEyes.size()-1, CV_RGB(255,255,255),1);
        //getEdgePoints(contours_canny.at(iREye),vedgePoints_all);
    }
    else {
        //If Contour Finding Fails Then Take Raw Edge points and mask L/R half of image
        cv::Canny( imgIn_thres, imgEdge_local, gi_CannyThresSmall,gi_CannyThres  );
        //Cover LEFT Eye Edges
        cv::Rect r(0,0,imgEdge.cols/2,imgEdge.rows);
        //imgEdge.copyTo(imgEdge_local);
        cv::rectangle(imgEdge_local,r,cv::Scalar(0),-1);


    }
    vedgePoints_all.clear();
    getEdgePoints(imgEdge_local,vedgePoints_all);
    detectEllipse(vedgePoints_all,qEllipsoids);


    if (qEllipsoids.size() > 0)
    {
        tDetectedEllipsoid dEll = qEllipsoids.top();
        drawEllipse(img_colour,dEll);

        cv::Point2f featurePnts[4];
        rcREye.points(featurePnts);
        ///Draw Left Eye Rectangle
        for (int j=0; j<4;j++) //Rectangle Eye
               cv::line(img_colour,featurePnts[j],featurePnts[(j+1)%4] ,CV_RGB(130,10,10),1);

        if (qEllipsoids.size() > 0)  qEllipsoids.pop();
        if (qEllipsoids.size() > 0) qEllipsoids.pop();
        if (qEllipsoids.size() > 0) qEllipsoids.pop();
        if (qEllipsoids.size() > 0) qEllipsoids.pop();

    }



    //DEBUG //
    //cv::bitwise_or(imgEdge,imgEdge_dbg,imgEdge_dbg);
    cv::imshow("Fish Edges ",imgEdge_dbg);
    cv::imshow("Fish Edges h",imgEdge);
    cv::imshow("Fish Threshold ",imgIn_thres);


    //Get All Edge Points Manually
    //Iterate Edge Image and extract edge points

    //cv::circle(img_colour,ptLEyeMid,1,CV_RGB(255,0,0),1);
    //cv::circle(img_colour,ptREyeMid,1,CV_RGB(0,250,0),1);
    img_colour.at<cv::Vec3b>(ptLEyeMid)[0] = 255; img_colour.at<cv::Vec3b>(ptLEyeMid)[1] = 0;
    img_colour.at<cv::Vec3b>(ptREyeMid)[0] = 0; img_colour.at<cv::Vec3b>(ptREyeMid)[1] = 250;


    cv::imshow("Fish CONTOUR ",img_contour);
//detect ellipse

    ///Draw Best 2 Ellipses



    ///Debug//
    cv::imshow("Debug EllipseFit",imgDebug);

    cv::imshow("Ellipse fit",img_colour);
    std::cout << "Done"  << std::endl;

}
