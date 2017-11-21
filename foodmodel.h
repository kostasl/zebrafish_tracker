#ifndef FOODMODEL_H
#define FOODMODEL_H

#include <string>
#include <QDebug>
//#include <QApplication>

//#include "larvatrack.h" //If included here it causes circular search if fishModel Defs.
#include "config.h"
#include "ellipse_detect.h"
#include "zfttracks.h"

extern float gDisplacementThreshold;


typedef cv::KeyPoint zfdblob;
typedef std::vector<zfdblob> zfdblobs;


class foodModel
{
public:
    foodModel();
    foodModel(zfdblob blob,zfdID ID);
    void updateState(zfdblob* fblob,int Angle, cv::Point2f bcentre,unsigned int nFrame);
    zfdID ID;
    int inactiveFrames; //Count Of Number Of Frames That this model Has not Been Matched To Any Fish
    zfdblob  zfoodblob;
    zftTrack zTrack;
    unsigned int nLastUpdateFrame; ///<-Holds the frame Number of the last State Update
};

typedef std::map<zfdID,foodModel* > foodModels;
typedef std::pair<zfdID, foodModel* > IDFoodModel;
#endif // FOODMODEL_H
