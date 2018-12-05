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
activeFrames = 0;
blobMatchScore = 0;
this->ID = lastFoodID;
ROIID = 0;
isTargeted = false; //Saves Location To Data File When True

}

foodModel::foodModel(zfdblob blob,zfdID ID):foodModel()
{
    inactiveFrames = 0;
    activeFrames = 0;
    blobMatchScore = 0;
    this->ID = ID;

    zTrack.id   = this->ID;
    zTrack.colour = CV_RGB(20,170,20);
    zTrack.centroid = blob.pt;




    zTrack.boundingBox.x = blob.pt.x - 6;
    zTrack.boundingBox.y = blob.pt.y - 6;
    zTrack.boundingBox.width = 12;
    zTrack.boundingBox.height = 12;

//    zTrack.boundingBox = cv::Rect(blob.pt.x - 5,blob.pt.y - 5,5,5);
    this->zfoodblob = blob;
}


foodModel::~foodModel()
{
    zTrack.pointStack.clear();
    zTrack.pointStackRender.clear();

}

void foodModel::updateState(zfdblob* fblob,int Angle, cv::Point2f bcentre,unsigned int nFrame,int matchScore,float szradius)
{

    blobMatchScore = matchScore;
    nLastUpdateFrame = nFrame; //Set Last Update To Current Frame
    this->zfoodblob      = *fblob;
    this->zTrack.pointStack.push_back(bcentre);
    this->zTrack.effectiveDisplacement = cv::norm(fblob->pt-this->zTrack.centroid);
    this->zTrack.centroid = bcentre;//fblob->pt; //Or Maybe bcentre
    this->blobRadius = szradius;
    zTrack.boundingBox.x = bcentre.x - 6;
    zTrack.boundingBox.y = bcentre.y - 6;
    zTrack.boundingBox.width = 12;
    zTrack.boundingBox.height = 12;


    inactiveFrames = 0;
    ///Optimization only Render Point If Displaced Enough from Last One
    if (this->zTrack.effectiveDisplacement > gDisplacementThreshold)
    {
        this->zTrack.pointStackRender.push_back(bcentre);
        //this->zTrack.active++;

        //this->zTrack.inactive = 0;
    }else {
        //this->zTrack.inactive++;
    }

    if (isTargeted )
        zTrack.colour = CV_RGB(240,210,10); //YelloW For Tracking
    else
        zTrack.colour = CV_RGB(0,200,10);

 //cv::Scalar(0,120,200)

}


std::ostream& operator<<(std::ostream& out, const foodModels& v)
{

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)

    //Check Through Models And Find The Closest Food To This FoodBlob
    //foodModels::iterator ft = v.begin();

    //while (ft != v.end())
    for (int i =0;i<v.size();i++)
    {
        foodModel* pfood = v.at(i);

         if (pfood->isTargeted) //Only Log The Marked Food
         {
            out << pfood->ID << "\t" << pfood->zTrack << "\n";

            pfood->zTrack.pointStack.clear();
            pfood->zTrack.pointStack.shrink_to_fit();
         }
    //++ft;
    }

    return out;
}


///
/// \brief operator << //Overloaded Stream Operator
/// Output Marked tracked foodModels State to Log File
/// \param out
/// \param h
/// \return
///
QTextStream& operator<<(QTextStream& out, const foodModels& v)
{

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)
    out.setRealNumberNotation(QTextStream::RealNumberNotation::FixedNotation );
    out.setRealNumberPrecision(2);


    //Set Global 1st Spine Direction (Helps to detect Errors)
    //Output Spine Point Angular Deviations from the previous spine/tail Segment in Degrees

    //foodModels::iterator ft = vfoodmodels.begin();
    //while (ft != vfoodmodels.end())
    for (int i =0;i<v.size();i++)
    {
        foodModel* pfood = (foodModel*)v.at(i);

         if (pfood->isTargeted) //Only Log The Marked Food
         {
            out << pfood->ID << "\t" << pfood->zTrack << "\n";

            pfood->zTrack.pointStack.clear();
            pfood->zTrack.pointStack.shrink_to_fit();
         }
    //++ft;
    }

    return out;
}
///
/// \brief getActiveFoodCount Aux. function returning the usuable food count - instead of just the number of instances given by .size()
/// \param vfoodmodels
/// \return
///
int foodModel::getActiveFoodCount(foodModels& vfoodmodels)
{
    int retNfood = 0;
    foodModels::iterator ft = vfoodmodels.begin();

    while (ft != vfoodmodels.end())
    {
        foodModel* pfood = ft->second;
        assert(pfood);

        // Render Food that has been on for A Min of Active frames / Skip unstable Detected Food Blob - Except If Food is being Tracked
        if (pfood->activeFrames < gcMinFoodModelActiveFrames && (!pfood->isTargeted))
        {
            ++ft; //Item Is not Counted
            continue;
        }

        ++ft;
        retNfood++; //only count the rendered Food Items ie. Active Ones
    }

return retNfood;
}




