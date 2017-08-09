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
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip> //for setprecision
#include <limits>
#include <string>

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

cv::Mat imgDebug;


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

void getEdgePoints(cv::Mat& imgEdgeIn,tEllipsoidEdges& vedgepoint)
{
   const float pxThres = 1.0;

   vedgepoint.clear();

   imgDebug = cv::Mat::zeros(imgEdgeIn.rows,imgEdgeIn.cols,CV_8UC1);

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

int detectEllipses(cv::Mat& imgIn,cv::Mat& imgOut,tEllipsoids& vellipses)
{

    const int minEllipseMajor   = gi_minEllipseMajor;
    const int maxEllipseMajor   = gi_maxEllipseMajor;
    const int thresMinVotes     = gi_VotesEllipseThres;
    const int minMinorEllipse   = 1;
    int accLength = imgIn.cols+imgIn.rows;
    double HighestVotes = 0;
    cv::Mat img_blur,img_edge,img_colour;
    cv::GaussianBlur(imgIn,img_blur,cv::Size(1,1),1,1);
    cv::Canny( img_blur, img_edge, gi_CannyThresSmall,gi_CannyThres  );
    //cv::findContours(frameCanny, contours_canny,hierarchy_canny, cv::RETR_CCOMP,cv::CHAIN_APPROX_NONE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE
    //Debug
    cv::imshow("Fish Edges ",img_edge);

    //Make Debug Img
    cv::cvtColor( imgIn,img_colour, cv::COLOR_GRAY2RGB);
//    imgIn.convertTo(img_colour,CV_8UC3);

    //Empty List
    vellipses.clear();
    tEllipsoidEdges vedgePoints_all; //All edge points from Image Of EDge detection
    tEllipsoidEdges vedgePoints_trial; //Containts edge points of specific ellipsoid trial


    int accumulator[accLength];


    //Iterate Edge Image and extract edge points
    getEdgePoints(img_edge,vedgePoints_all);
    memset(accumulator,0,sizeof(int)*(accLength)); //Reset Accumulator MAtrix

    std::cout << "== Start === "  << std::endl;
    ///Loop through All Edge points (3)
    for (tEllipsoidEdges::iterator it1 = vedgePoints_all.begin();it1 != vedgePoints_all.end();++it1)
    {
        cv::Point2f ptxy1 = (*it1).ptEdge;
        if (ptxy1.x == 0 && ptxy1.y == 0)
            continue ; //point has been deleted
        cv::Point2f ptxy2;



        ///(4)
        for (tEllipsoidEdges::iterator it2 = vedgePoints_all.begin();it2 != vedgePoints_all.end(); ++it2 )
        {

            ptxy2 = (*it2).ptEdge;

            if (ptxy2.x == 0 && ptxy2.y == 0)
                continue ; //point has been deleted


            double d = cv::norm(ptxy2-ptxy1);
            if (d < minEllipseMajor || d > maxEllipseMajor)
                continue;

            //Use Eqns 1-4 and calculate Ellipse params
            cv::Point2f ptxy0;
            ptxy0.x = (ptxy2.x + ptxy1.x )/2.0;  //--(1)
            ptxy0.y = (ptxy2.y + ptxy1.y)/2.0;  //--(2)
            double a = d/2.0; //[(x 2 – x 1 )^2 + (y 2 – y 1 )^2 ] /2 //--(3) a the half-length of the major axis

            double alpha = atan2(ptxy2.y - ptxy1.y,ptxy2.x - ptxy1.x);//atan [(y 2 – y 1 )/(x 2 – x 1 )] //--(4) α the orientation of the ellipse

            ///Step (6) - 3rd Pixel;
            for (tEllipsoidEdges::iterator it3 = vedgePoints_all.begin();it3 != vedgePoints_all.end(); ++it3 )
            {

                cv::Point2f ptxy3 = (*it3).ptEdge;
                if (ptxy3.x == 0 && ptxy3.y == 0)
                    continue ; //point has been deleted

                double d = cv::norm(ptxy0-ptxy3);
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
                int b = round(sqrt(bb));
                ///Step 8
                if (b > 0)
                {
                    accumulator[b]++; //increment accumulator for this minor Axis
                    //Add Point to tracked List
                    //vedgePoints_trial
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
                tDetectedEllipsoid ellipse(cv::RotatedRect(ptxy0,cv::Size2f(2*a,2*idx), alpha*(180/M_PI)));
                vellipses.push_back(ellipse);
                cv::ellipse(img_colour,ellipse.rectEllipse,CV_RGB(250,50,50),1,cv::LINE_8);

                cv::circle(img_colour,ptxy0,1,CV_RGB(0,0,255),1);
                img_colour.at<cv::Vec3b>(ptxy1)[1] = 255; img_colour.at<cv::Vec3b>(ptxy1)[2] = 255;
                img_colour.at<cv::Vec3b>(ptxy0)[1] = 255; img_colour.at<cv::Vec3b>(ptxy0)[2] = 255;
                //cv::circle(img_colour,ptxy1,1,CV_RGB(0,255,255),1);
                //cv::circle(img_colour,ptxy2,1,CV_RGB(0,255,255),1);
                //Debug Mark As Good Pair
                imgDebug.at<uchar>(ptxy1) = 255;
                imgDebug.at<uchar>(ptxy2) = 255;


            }else
            {//Mark As Dull Pair
                imgDebug.at<uchar>(ptxy1) = 25;
                imgDebug.at<uchar>(ptxy2) = 25;
            }


            if (HighestVotes < dvotesMax)
            {
                HighestVotes = dvotesMax;
                std::cout << "mxVot:" << HighestVotes << std::endl;
            }


        ///Step 12 - Remove the points from the image Before Restarting

        //it2 = vedgePoints.erase(it2);
        //
            //ptxy1.x = 0; ptxy1.y = 0;
            //ptxy2.x = 0; ptxy2.y = 0;
          //  it2->x = 0; it2->y = 0; //Delete Point


        ///Step 13 - Clear Accumulator
        memset(accumulator,0,sizeof(int)*(accLength)); //Reset Accumulator MAtrix

        } //Loop through each 2nd point in pair

//        cv::waitKey(1);
        //it1 = vedgePoints.erase(it1);
        //it1->x = 0; it1->y = 0; //Delete Point
    } //Loop through all  point as 1st point pair (Prob: pairs can be repeated)


    ///Debug//
    cv::imshow("Debug D",imgDebug);

    cv::imshow("Ellipse fit",img_colour);
    std::cout << "Done"  << std::endl;

}
