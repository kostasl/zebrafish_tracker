#include "eyesdetector.h"
#include <random>

/// implements a reinforcement learning technique to discover parameters for
/// best extraction of eye information from video frames
tEyeDetectorState EyesDetector::DrawNextState(tEyeDetectorState currentState,double RewardScore)
{

   ulong idxVal_i = currentState.iSegThres1-baseIdxRow;
   ulong idxVal_j = currentState.VergenceAngle-baseIdxCol;
   tEyeDetectorState nextState = currentState;
    //Draw if we are exploring or Greedy
   if (drand48() < pExplore)
   {
       bExploreMove = true;
   }

    //Explore: Choose Action that determines Next State
    // Todo : Maybe policy could be modified to context, ie  EyeVergence
    // Increase/Decrease Threshold
    if (drand48()<0.5 & bExploreMove)
    {
        nextState.iSegThres1++;
    }else {
        nextState.iSegThres1--;
    }

    //Policy: Take greedy action with some prob
    if (!bExploreMove)
    {

        ulong idxState_ActUp = std::min(baseIdxRow,idxVal_i++);
        ulong idxState_ActDown = std::max((ulong)0,idxVal_i--);

        //Calc Marginal Value For each Action taken

        if (mStateValue[idxState_ActUp][idxVal_j] > mStateValue[idxState_ActDown][idxVal_j] )
                nextState.iSegThres1 = idxState_ActUp;
        else {

        }
    }


    nextState.iSegThres1 = std::max((ulong)nextState.iSegThres1,baseIdxRow);
    nextState.iSegThres1 = std::min((ulong)nextState.iSegThres1,baseIdxRow);


    // Explore Next State Value


}



void EyesDetector::UpdateStateValue(tEyeDetectorState nextState)
{

    ulong idxVal_i = currentState.iSegThres1-baseIdxRow;
    ulong idxVal_j = currentState.VergenceAngle-baseIdxCol;

    //Get Next states Value and update current - Greedy
    ulong idxNextVal_i = nextState.iSegThres1-baseIdxRow;
    ulong idxNextVal_j = nextState.VergenceAngle-baseIdxCol;
    //Update The Current State value by propagating value of next state backwards

    //TD learning - Use immediate Reward and add discounted future state rewards
    mStateValue[idxVal_i][idxVal_j] = mStateValue[idxVal_i][idxVal_j] +
                                      alpha*(RewardScore + gamma*mStateValue[idxNextVal_i][idxNextVal_j]-mStateValue[idxVal_i][idxVal_j] );

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


