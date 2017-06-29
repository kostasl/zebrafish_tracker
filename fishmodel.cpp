#include "fishmodel.h"

fishModel::fishModel()
{

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

float fishModel::vergenceAngle()
{

}


