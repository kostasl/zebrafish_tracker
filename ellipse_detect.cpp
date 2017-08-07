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


inline int getMax(int* darray,int length,double& votes)
{
    double max=darray[0];
    int maxIdx = 0;
    //find max and mins
    for(int j=0; j<length; j++)
    {
        if(max<darray[j])
        {
            max=darray[j];
            maxIdx = j;
        }
    }
    votes = max;
    return maxIdx;
}

void getEdgePoints(cv::Mat& imgEdgeIn,std::vector<cv::Point2f>& vedgepoint)
{
   const float pxThres = 1.0;

   vedgepoint.clear();
   /*
  for (int y = 0; y < imgEdgeIn.rows; ++y)
  {
      const float* row_ptr = imgEdgeIn.ptr<float>(y);
      for (int x = 0; x < imgEdgeIn.cols; ++x)
      {
          float value = row_ptr[x];
          if ( value >  pxThres)

      }
  }
  */
  for(int i=0; i<imgEdgeIn.rows; i++)
      for(int j=0; j<imgEdgeIn.cols; j++)
           if ( imgEdgeIn.at<uchar>(i,j) >  pxThres)
               vedgepoint.push_back(cv::Point(i,j));

}

int detectEllipses(cv::Mat& imgIn,cv::Mat& imgOut,std::vector<cv::RotatedRect>& vellipses)
{
    const int minEllipseMajor = 4;
    const int thresMinVotes = 160;
    int accLength = imgIn.cols+imgIn.rows;
    cv::Mat img_blur,img_edge;
    cv::GaussianBlur(imgIn,img_blur,cv::Size(1,1),1,1);
    cv::Canny( img_blur, img_edge, gi_CannyThresSmall,gi_CannyThres  );
    //cv::findContours(frameCanny, contours_canny,hierarchy_canny, cv::RETR_CCOMP,cv::CHAIN_APPROX_NONE , cv::Point(0, 0) ); //cv::CHAIN_APPROX_SIMPLE
    //Debug
    cv::imshow("Fish Edges ",img_edge);

    std::vector<cv::Point2f> vedgePoints;

    int accumulator[accLength];


    //Iterate Edge Image and extract edge points
    getEdgePoints(img_edge,vedgePoints);


    ///Loop through All Edge points (3)
    for (std::vector<cv::Point2f>::iterator it1 = vedgePoints.begin();it1 != vedgePoints.end(); ++it1 )
    {
        cv::Point2f ptxy1 = *it1;
        cv::Point2f ptxy2;
        memset(accumulator,0,sizeof(int)*(accLength)); //Reset Accumulator MAtrix

        ///(4)
        for (std::vector<cv::Point2f>::iterator it2 = vedgePoints.begin();it2 != vedgePoints.end(); ++it2 )
        {
            ptxy2 = *it2;
            double d = cv::norm(ptxy1-ptxy2);
            if (d < minEllipseMajor)
                continue;

            //Use Eqns 1-4 and calculate Ellipse params
            cv::Point2f ptxy0;
            ptxy0.x = (ptxy1.x + ptxy2.x )/2.0;  //--(1)
            ptxy0.y = (ptxy1.y + ptxy2.y)/2.0;  //--(2)
            double a = d/2.0; //[(x 2 – x 1 )^2 + (y 2 – y 1 )^2 ] /2 //--(3) a the half-length of the major axis

            double alpha = atan(ptxy2.y - ptxy1.y)/(ptxy2.x - ptxy1.x);//atan [(y 2 – y 1 )/(x 2 – x 1 )] //--(4) α the orientation of the ellipse

            ///Step (6) - 3rd Pixel;
            for (std::vector<cv::Point2f>::iterator it3 = vedgePoints.begin();it3 != vedgePoints.end(); ++it3 )
            {
                cv::Point2f ptxy3 = *it3;
                double d = cv::norm(ptxy0-ptxy3);
                double dd = d*d;

                if (d > a || d ==0) //Candidate 3rd point of minor axis distance needs to be less than alpha away
                    continue;

                //Calculate Minor Axis
                double aa = a*a;
                double f = cv::norm(ptxy2-ptxy3);
                double ff = f*f;
                double c = cv::norm(ptxy1-ptxy3);

                ///Step 7 - Calc the length of minor axis
                double costau = ( aa + dd - ff)/(2.0*a*d);
                double coscostau = costau*costau;  //eqn (6)
// b = sqrt( (aSq * thirdPtDistsSq(K) .* sinTauSq) ./ (aSq - thirdPtDistsSq(K) .* cosTau.^2 + eps) );
                double bb = aa*dd*(1.0-coscostau)/(aa - dd * coscostau + 0.00001); //(5)
                int b = round(sqrt(bb));
                ///Step 8
                if (b > 0)
                    accumulator[b]++; //increment accumulator for this minor Axis

            ///Step 9 Loop Until All Pixels 3rd are computed for this pair of pixes
            }

            ///Step 10 //Find Max In accumulator array. The related length is the possible length of minor axis for assumed ellipse.
            double dvotesMax;
            int idx  = getMax(accumulator,accLength,dvotesMax);

            //idx Is the size of the minor axis
            if (dvotesMax > thresMinVotes) //Found ellipse
            {
                ///Step 11 Output Ellipse Parameters
                //cv::RotatedRect ellipse(ptxy0,ptxy1,ptxy2);
                cv::RotatedRect ellipse(ptxy0,cv::Size2f(idx,2*a),alpha);
                vellipses.push_back(ellipse);
                cv::ellipse(imgIn,ellipse,CV_RGB(250,50,50),1,cv::LINE_8);

            }

        ///Step 12 - Remove the points from the image Before Restarting


        } //Loop through each 2nd point in pair


    } //Loop through all  point as 1st point pair (Prob: pairs can be repeated)

cv::imshow("Ellipse",imgIn);
}
