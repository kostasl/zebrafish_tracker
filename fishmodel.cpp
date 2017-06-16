#include "fishmodel.h"

fishModel::fishModel()
{

        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());

}

fishModel::fishModel(cvb::CvTrack* track):fishModel()
{

    this->ID track->id;
    this->track = track; //Copy Pointer
}


