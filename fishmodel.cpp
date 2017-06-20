#include "fishmodel.h"

fishModel::fishModel()
{

        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());

        this->leftEyeHull.clear();
        this->rightEyeHull.clear();

}

fishModel::fishModel(cvb::CvTrack* track):fishModel()
{

    this->ID    = track->id;
    this->track = track; //Copy Pointer
}


