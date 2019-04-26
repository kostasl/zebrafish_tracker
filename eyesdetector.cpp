#include "eyesdetector.h"
#include <random>
#include <QDebug>
#include <assert.h>
/// \brief implements a reinforcement learning technique to discover parameters for
/// best extraction of eye information from video frames
///

/// \brief DrawNextAction- Implements policy and returns action as a move in state space (ie increases or decreases the threshold)
/// this action must then be passed on to the "envinment" here the fitting of ellipsoids -
/// the fit score and eye v angle determine the next state - which needs to be passed to UpdateStateValue before drawing the next Action
/// \returns a modified currentState reflecting the action taken in the state space based on (greedy policy)
tEyeDetectorState EyesDetector::DrawNextAction(tEyeDetectorState currentState)
{

   ulong idxVal_i = currentState.iSegThres1-baseIdxRow;
   ulong idxVal_j = currentState.VergenceState-baseIdxCol;
   tEyeDetectorState nextState = currentState;
    //Draw if we are exploring or Greedy
   if (drand48() < pExplore)
   {
       bExploreMove = true;
   }
   else
       bExploreMove = false;

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

        ulong idxState_ActUp = std::min( mStateValue.size()-1,idxVal_i+1);
        ulong idxState_ActDown = std::max((ulong)0,idxVal_i-1);
        double qUp   = 0; //Marginal Value of going action up thresh.
        double qDown = 0;//Marginal Value of going action down thresh.

        //Calc Marginal Value For each Action taken (the state can shift the eye vergence so marginalize across a range )
        for (ulong i=0;i < mStateValue[idxState_ActUp].size();i++ )
                qUp += mStateValue[idxState_ActUp][i];

        for (ulong i=0;i < mStateValue[idxState_ActDown].size();i++ )
                qDown += mStateValue[idxState_ActDown][i];
        // if equal value
        if (qUp == qDown ) //Add noise to randomly choose direction
            qUp += (drand48()-0.5);

        if (qUp > qDown )
                nextState.iSegThres1++;
        else {
                nextState.iSegThres1--;
        }
    } // Taking greedy action / not exploring

    //Do not exceed state space limits - bounce off walls;
    nextState.iSegThres1 = std::max((ulong)nextState.iSegThres1,baseIdxRow);
    nextState.iSegThres1 = std::min((ulong)nextState.iSegThres1,baseIdxRow+mStateValue.size()-1);


    //Return next state to the environment - which will evaluate the fit score to update the value funct.
    return nextState;


}


/// \brief call it after drawing a new state and evaluating its fitness
/// \param RewardScore the immediate value of the current state and
double EyesDetector::UpdateStateValue(tEyeDetectorState nextState,double RewardScore)
{

    ulong idxVal_i = (ulong)currentState.iSegThres1-baseIdxRow;
    ulong idxVal_j = (ulong)currentState.VergenceState;

    //Get Next states Value and update current - Greedy
    ulong idxNextVal_i = (ulong)nextState.iSegThres1-baseIdxRow;
    ulong idxNextVal_j = (ulong)nextState.VergenceState;
    //Update The Current State value by propagating value of next state backwards

    //TD learning - Use immediate Reward and add discounted future state rewards
    mStateValue[idxVal_i][idxVal_j] = mStateValue[idxVal_i][idxVal_j] +
                                      alpha*(RewardScore + gamma*mStateValue[idxNextVal_i][idxNextVal_j]-mStateValue[idxVal_i][idxVal_j] );

    return mStateValue[idxVal_i][idxVal_j];
}


EyesDetector::EyesDetector(int RangeValThres_min,int RangeValThres_max,int AngleVal_min,int AngleVal_max): EyesDetector(RangeValThres_max-RangeValThres_min, 4)
{
///Keep The min values so We can translate, locatre State  Values
baseIdxRow = RangeValThres_min;
//baseIdxCol = AngleVal_min;

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

tEyeDetectorState EyesDetector::getCurrentState()
{
    return currentState;
}


