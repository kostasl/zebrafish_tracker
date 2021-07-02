#ifndef FOODMODEL_H
#define FOODMODEL_H

#include <string>
#include <QDebug>
//#include <QApplication>

//#include "larvatrack.h" //If included here it causes circular search if fishModel Defs.
//#include "config.h"
#include "ellipse_detect.h"
#include "zfttracks.h"

extern float gDisplacementThreshold;
//extern uint gi_MaxFoodID;
extern trackerState gTrackerState;

typedef cv::KeyPoint zfdblob;
typedef std::vector<zfdblob> zfdblobs;

class preyModel; //fwd definition

typedef std::map<zfdID,preyModel* > foodModels;
typedef std::pair<zfdID, preyModel* > IDFoodModel;


class preyModel
{
public:
    preyModel();
    preyModel(zfdblob blob,zfdID ID);
    ~preyModel();

    cv::Point2f  predictMove();// Draws Prediction Of next position
    void updateState(zfdblob fblob,int Angle, cv::Point2f bcentre,unsigned int nFrame,int matchScore,float szradius);
    static int getActiveFoodCount(foodModels& vfoodmodels);
    bool isUnused(); //Contains the logic of when to delete food item
    zfdID ID;
    int ROIID;
    int inactiveFrames; //Count Of Number Of Frames That this model Has not Been Matched To Any Fish
    int activeFrames; //Count Of Last Consecutive Active Frames
    float blobRadius; //To Use For Size Estimation
    zfdblob  zfoodblob;
    int blobMatchScore;
    zftTrack zTrack;


    zfdblob previous_position;
    cv::Point2f velocity;
    float omegaDeg;

    float headingTheta;
    float muPropulsion;
    float sigmaPropulsion;
    float df_propulsion;
    float df_friction;
    float mass;
    float dTheta;
    float muTurn;
    float sigmaTurn;
    uint16_t color;
    uint16_t BGcolor;

    unsigned int nLastUpdateFrame; ///<-Holds the frame Number of the last State Update
    bool isTargeted;
    bool isActive;
    bool isNew;

    private:
    cv::Point2f alpha_beta_TrackingFilter_step(cv::Point2f blobPt); //Filter Position Updates based on a simple tuned g-h filter
    cv::Point2f ptEstimated; //Filter Internal Vars
    cv::Point2f ptPredicted;

    double g_FPS_scaling = 400.0/gTrackerState.gfVidfps;
    double dx = 0.01*g_FPS_scaling;
    double dy = 0.01*g_FPS_scaling;
    // These values have been tested in R using pre-recorded Prey Data
    double dt = g_FPS_scaling*1.0; // Filter TimeStep (a dt 1.0 works when fps 410)

    //constexpr static
    double g = g_FPS_scaling*1.0/100.0; //Measurement Scaling - 2orders Larger than h works best
    double h = g_FPS_scaling*1.0/1000.0; //Prediction Scaling - Gain Smaller when highly noisy environment


};


//To Write Selected Food Vector Data To File
std::ostream& operator<<(std::ostream& out, const foodModels& h);
QTextStream& operator<<(QTextStream& out, const foodModels& h);


//Top Item Is one With Least Penalty (Score)
class CompareFoodScore {
    public:
    bool operator()(preyModel*& t1, preyModel*& t2) // Returns true if t1 is less than t2 /Ordering Highest 1st
    {
       return t1->blobMatchScore < t2->blobMatchScore;
    }
};

typedef std::priority_queue<preyModel*,std::vector<preyModel*>,CompareFoodScore> qfoodModels;

/// \fn inline void cvReleaseFoodModels(fishModels &fishes)
/// \brief Clear Fish LIst
/// \param fishmodles List
/// \see
inline void ReleaseFoodModels(foodModels &vfood)
{
  for (foodModels::iterator it=vfood.begin(); it!=vfood.end(); ++it)
  {
      preyModel* pfood = (*it).second;
        //Let ReleaseTracks Handle This
//      if (fish->track)
//      {
//         fish->track->pointStack.clear();
//         delete fish->track;
//      }

      delete pfood;
  }
  vfood.clear();
  gTrackerState.gi_MaxFoodID = 1; //Reset COunter from Where IDs are drawn
}


#endif // FOODMODEL_H
