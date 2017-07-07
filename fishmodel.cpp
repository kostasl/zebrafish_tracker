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

void fishModel::resetSpine()
{
//spline.reserve(c_spinePoints);
this->spline.clear();
this->spline.push_back(this->coreTriangle[2]);
spline[0] = this->coreTriangle[2];
int len = 6;
for (int i=1;i<c_spinePoints;i++)
{
    cv::Point sp;
    sp.x = spline[i-1].x + len*sin(this->bearingRads);
    sp.y = spline[i-1].y + len*cos(this->bearingRads);
    spline.push_back(sp);
}

}

float fishModel::vergenceAngle()
{

}


