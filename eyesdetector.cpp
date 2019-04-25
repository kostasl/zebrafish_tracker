#include "eyesdetector.h"
#include <random>

/// implements a reinforcement learning technique to discover parameters for
/// best extraction of eye information from video frames
tEyeDetectorState EyesDetector::DrawNextState(tEyeDetectorState currentState)
{

   int idxVal_i = currentState.iSegThres1-baseIdxRow;
   int idxVal_j = currentState.VergenceAngle-baseIdxCol;
   tEyeDetectorState nextState = currentState;
    //Draw if we are exploring or Greedy
   if (drand48() < pExplore)
   {
       bExploreMove = true;
   }

    //Choose Action that determines Next State
    // Increase/Decrease Threshold
    if (drand48()<0.5 )
    {
        nextState.iSegThres1++;
        nextState.iSegThres1 = std::max(nextState.iSegThres1,baseIdxRow);
    }else {
        nextState.iSegThres1--;
    }

    // Explore Next State Value


    //Get Next states Value and update current - Greedy
   int idxNextVal_i = nextState.iSegThres1-baseIdxRow;
   int idxNextVal_j = nextState.VergenceAngle-baseIdxCol;
    //Update The Current State value by propagating value of next state backwards
    mStateValue[idxVal_i][idxVal_j] = mStateValue[idxVal_i][idxVal_j] + 0.1*(mStateValue[idxNextVal_i][idxNextVal_j]-mStateValue[idxVal_i][idxVal_j] );

}


EyesDetector::EyesDetector(int RangeValThres_min,int RangeValThres_max,int AngleVal_min,int AngleVal_max): EyesDetector(RangeValThres_max-RangeValThres_min, (AngleVal_max-AngleVal_min))
{
///Keep The min values so We can translate, locatre State  Values
baseIdxRow = RangeValThres_min;
baseIdxCol = AngleVal_min;

    for(int i = 0;i < mStateValue.size();++i)
    {
          for(int j = 0;j < mStateValue[i].size();++j)
          {
            mStateValue[i][j] = 100; // start with high value on all - so we get to explore more
          }
    }


}

EyesDetector::EyesDetector()
{


}

EyesDetector::~EyesDetector()
{

}


