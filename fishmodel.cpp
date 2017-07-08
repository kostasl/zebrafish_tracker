#include "fishmodel.h"

fishModel::fishModel()
{

        c_spinePoints = 8;
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());

        this->mouthPoint.x = 0;
        this->mouthPoint.y = 0;
        this->leftEyeHull.clear();
        this->rightEyeHull.clear();

}

fishModel::fishModel(cvb::CvTrack* track):fishModel()
{

    this->ID    = track->id;
    this->blobLabel = track->label;
    this->track = track; //Copy Pointer
    this->coreTriangle[2].x = track->centroid.x;
    this->coreTriangle[2].y = track->centroid.y;
    this->resetSpine();
}

float fishModel::leftEyeAngle()
{

    if (this->leftEyeRect.size.width < this->leftEyeRect.size.height)
        return (this->leftEyeRect.angle-90.0)*CV_PI/180.0;
    else
        {
         return (this->leftEyeRect.angle)*CV_PI/180.0;
        }


}

/// \brief return (corrected for leading edge from horizontal line -Pi ... +Pi) Rectangle angle in Radians
float fishModel::rightEyeAngle()
{

if (this->rightEyeRect.size.width < this->rightEyeRect.size.height)
    return (this->rightEyeRect.angle-90.0)*CV_PI/180.0;
else
    {
     return (this->rightEyeRect.angle)*CV_PI/180.0;
    }

}

///
/// \brief fishModel::resetSpine make a straight Spline pointing towards Blobs Bearings Angle
///
void fishModel::resetSpine()
{
    //spline.reserve(c_spinePoints);
    this->spline.clear();
    this->splineTheta.clear();
    this->spline.push_back(this->coreTriangle[2]);


    for (int i=1;i<c_spinePoints;i++)
    {

        cv::Point2f sp;
        sp.x = spline[i-1].x - c_spineSegL*cos(this->bearingRads);
        sp.y = spline[i-1].y - c_spineSegL*sin(this->bearingRads);

        splineTheta.push_back(this->bearingRads); //initial Angles are all the same
        spline.push_back(sp); //Place Initial Point of spline
    }

}
///
/// \brief fishModel::getSpine Make Spine Points based on Polar Coords
/// /Calculates Spine Positions
void fishModel::getSpline(std::vector<cv::Point2f>& outspline)
{

    this->spline.clear();
    this->spline.push_back(this->coreTriangle[2]);

    for (int i=1;i<c_spinePoints;i++)
    {
        cv::Point2f sp;
        sp.x = spline[i-1].x - c_spineSegL*cos(this->splineTheta[i-1]);
        sp.y = spline[i-1].y - c_spineSegL*sin(this->splineTheta[i-1]);
        spline.push_back(sp);
    }

    outspline = spline; //Copy To output

}

float fishModel::vergenceAngle()
{

}


///
/// \brief fishModel::fitSpineToContour
/// \param contours_body
/// \param idxInnerContour
/// \param idxOuterContour
/// \return fitness error score
///
double fishModel::fitSpineToContour(std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour)
{

    cv::Scalar TRACKER_COLOURMAP[] ={CV_RGB(150,150,150),
                                     CV_RGB(200,100,100),
                                     CV_RGB(150,200,50),
                                     CV_RGB(50,250,00),
                                     CV_RGB(150,150,00),
                                     CV_RGB(250,250,00),
                                     CV_RGB(200,200,80),
                                     CV_RGB(20,200,180)};



    ///Param sfish model should contain initial spline curve (Hold Last Frame Position)

    //Run Until Convergence Error is below threshold - Or Change is too small

    ///Compute Error terms for all data points/obtain local quadratic approx of fsd
    //For each contour Point
    std::vector<cv::Point> contour = contours_body[idxOuterContour];
    std::vector<cv::Point2f> spline = sfish.spline;

    //Measure squared Distance error to closest Curve(spline) Point
        //Add to total error
    double dfitPtError[contour.size()][spline.size()]; //Per Spine Point Fit Error
    double dfitPtError_last[contour.size()][spline.size()]; //Last Iteration Fit Error per spine point
    double dfitPtError_change[contour.size()][spline.size()]; //Aproximate Error Change per spine point between iteration
    std::vector< std::vector<cv::Point2f> > dfitPtErrorVector(contour.size(),std::vector<cv::Point2f>(spline.size()) ); //Aproximate Error Gradient vector (Used for Corrections)
    double dfitPtError_total = 0.0;
    double dfitPtError_total_last = 0.0;
    double dDifffitPtError_total = 1000.0;


    //Init Vectors
    //std::fill(dDifffitPtError.begin(), dDifffitPtError.end(), 0); //Difference in Spine Errors
    //std::fill(dfitPtError_last.begin(), dfitPtError_last.end(), 0); //Last Spine Error Vector
    memset(dfitPtError,0.0,contour.size()*spline.size()*sizeof(double));


 //while (dDifffitPtError_total > 100.0) //While Error Change Above Threshold
 {
     dfitPtError_total_last = dfitPtError_total;
     dfitPtError_total = 0.0; //Reset
     //dfitPtError_last = dfitPtError;

     //std::fill(dfitPtError.begin(), dfitPtError.end(), 0);
     //std::fill(dfitPtErrorVector.begin(), dfitPtErrorVector.end(), cv::Point(0,0));


    for (int i=0;i<contour.size();i++) //For Each Data point
    {
        int idxNear = 0;
        cv::Point2f ptsrc = contour[i];
        cv::Point2f ptNear = spline[idxNear];
        double distX = pow(ptsrc.x-ptNear.x,2);
        double distY = pow(ptsrc.y-ptNear.y,2);
        double mindist = distX + distY;
        double dist;
        ///Find Closest Curve POint to contour point
        for (int j=1; j<spline.size();j++)
        {
            dfitPtErrorVector[i][j].x = 0.0;
            dfitPtErrorVector[i][j].y = 0.0;
            //cv::Point vecSrcToSpline = ptsrc-spline[j];
            //dist = pow(vecSrcToSpline.x,2)+pow(vecSrcToSpline.y,2); //Squared euclidian distance error
            distX = ptsrc.x-spline[j].x;
            distY = ptsrc.y-spline[j].y;
            dist  = pow(distX,2) + pow(distY,2);
            if (dist < mindist)
            {
                idxNear = j;
                ptNear = spline[idxNear];
                mindist = dist;
            }

        }
        ///Measure Error
        //Add to error for that spine point
        dfitPtError_last[i][idxNear]   = dfitPtError[i][idxNear]; //Save in matrix of contour / Spine association
        dfitPtError[i][idxNear]        = mindist; //Save in matrix of contour / Spine association
        dfitPtError_change[i][idxNear] = dfitPtError_last[i][idxNear]-dfitPtError[i][idxNear];
        dfitPtErrorVector[i][idxNear].x  = (ptsrc.x-ptNear.x)*0.1; //Move Closer Scaled By Square Error
        dfitPtErrorVector[i][idxNear].y  = 0; //ptsrc.y-spline[idxNear].y;
        dfitPtError_total              +=mindist;

    }    //Total Spine Fitness Has been measured
    dDifffitPtError_total = abs(dfitPtError_total -  dfitPtError_total_last);

    //Update Spline Positions by  Approx Error Gradient
    for (int i=0; i<contour.size();i++)
        for (int j=0; j<spline.size();j++)
    {

        if (abs(dfitPtError_change[i][j]) > 0.0)
        {
            //dfitPtErrorVector[j]= dfitPtErrorVector[j]/dfitPtError_total;
            spline[j].x += ((double)dfitPtErrorVector[i][j].x);//*dfitPtError_change[i][j]/dfitPtError_total; //Weighted Update
            spline[j].y += ((double)dfitPtErrorVector[i][j].y);//*dfitPtError_change[i][j]/dfitPtError_total; //Weighted Update
            qDebug() << "sj"<< j << " ex:" << ((double)dfitPtErrorVector[i][j].x)*dfitPtError_change[i][j]/dfitPtError_total << " ey:" << ((double)dfitPtErrorVector[i][j].y)*dfitPtError_change[i][j]/dfitPtError_total;
            assert(!(isnan(spline[j].y) || isnan(spline[j].x)));
        }
    }

    for (int j=0; j<8;j++) //Rectangle Eye
    {
        cv::circle(frameDebugC,spline[j],2,TRACKER_COLOURMAP[j],1);
    }
    cv::imshow("Debug C",frameDebugC);


    qDebug() << "D err:" << dDifffitPtError_total;
    //cv::waitKey(100);
} //While Error Change Is larger Than

sfish.spline = spline;
}



