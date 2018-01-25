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
extern uint gi_MaxFoodID       = 0;

typedef cv::KeyPoint zfdblob;
typedef std::vector<zfdblob> zfdblobs;


class foodModel
{
public:
    foodModel();
    foodModel(zfdblob blob,zfdID ID);
    ~foodModel();
    void updateState(zfdblob* fblob,int Angle, cv::Point2f bcentre,unsigned int nFrame,int matchScore);
    zfdID ID;
    int inactiveFrames; //Count Of Number Of Frames That this model Has not Been Matched To Any Fish
    int activeFrames; //Count Of Last Consecutive Active Frames
    zfdblob  zfoodblob;
    int blobMatchScore;
    zftTrack zTrack;
    unsigned int nLastUpdateFrame; ///<-Holds the frame Number of the last State Update
};

typedef std::map<zfdID,foodModel* > foodModels;
typedef std::pair<zfdID, foodModel* > IDFoodModel;


//Top Item Is one With Least Penalty (Score)
class CompareFoodScore {
    public:
    bool operator()(foodModel*& t1, foodModel*& t2) // Returns true if t1 is less than t2 /Ordering Highest 1st
    {
       return t1->blobMatchScore < t2->blobMatchScore;
    }
};

typedef std::priority_queue<foodModel*,std::vector<foodModel*>,CompareFoodScore> qfoodModels;

/// \fn inline void cvReleaseFoodModels(fishModels &fishes)
/// \brief Clear Fish LIst
/// \param fishmodles List
/// \see
inline void ReleaseFoodModels(foodModels &vfood)
{
  for (foodModels::iterator it=vfood.begin(); it!=vfood.end(); ++it)
  {
      foodModel* pfood = (*it).second;
        //Let ReleaseTracks Handle This
//      if (fish->track)
//      {
//         fish->track->pointStack.clear();
//         delete fish->track;
//      }

      delete pfood;
  }
  vfood.clear();
  gi_MaxFoodID = 1; //Reset COunter from Where IDs are drawn
}


#endif // FOODMODEL_H
