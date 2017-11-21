#include "foodmodel.h"
#include "ellipse_detect.h"
#include "config.h"
#include "zfttracks.h"

#include <string>
#include <QDebug>
//#include <QApplication>

//#include "larvatrack.h" //If included here it causes circular search if fishModel Defs.


static int lastFoodID = 0;

foodModel::foodModel()
{
lastFoodID++;
inactiveFrames = 0;
this->ID = lastFoodID;

}

foodModel::foodModel(zfdblob blob,zfdID ID)
{
    inactiveFrames = 0;
    this->ID = ID;

    zTrack.id   = this->ID;
    //zTrack.colour = CV_RGB(0,10,200);
    zTrack.centroid = blob.pt;
    this->zfoodblob = blob;
}


void foodModel::updateState(zfdblob* fblob,int Angle, cv::Point2f bcentre,unsigned int nFrame)
{

    nLastUpdateFrame = nFrame; //Set Last Update To Current Frame
    this->zfoodblob      = *fblob;
    this->zTrack.pointStack.push_back(bcentre);
    this->zTrack.effectiveDisplacement = cv::norm(fblob->pt-this->zTrack.centroid);
    this->zTrack.centroid = bcentre;//fblob->pt; //Or Maybe bcentre
    inactiveFrames = 0;
    ///Optimization only Render Point If Displaced Enough from Last One
    if (this->zTrack.effectiveDisplacement > gDisplacementThreshold)
    {
        this->zTrack.pointStackRender.push_back(bcentre);
        this->zTrack.active++;

        this->zTrack.inactive = 0;
    }else {
        this->zTrack.inactive++;
    }


}
