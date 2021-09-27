#include "config.h"
#include "eyesdetector.h"
#include <random>
#include <QDebug>
#include <assert.h>


#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include "cereal/types/vector.hpp"
#include <fstream>

/// \brief implements a reinforcement learning technique to discover parameters for
/// best extraction of eye information from video frames
///

/// \brief DrawNextAction- Implements policy and returns action as a move in state space (ie increases or decreases the threshold)
/// this action must then be passed on to the "environment" here the fitting of ellipsoids -
/// the fit score and eye v angle determine the next state - which needs to be passed to UpdateStateValue before drawing the next Action
/// \returns a modified currentState reflecting the action taken in the state space based on (greedy policy)
tEyeDetectorState EyesDetector::DrawNextAction(tEyeDetectorState currentState)
{
    std::uniform_int_distribution<> distrA(0, mStateValue.size()-1);
    std::uniform_int_distribution<> distrB(0, mStateValue[0].size()-1);

   int idxVal_i = std::min((int)mStateValue.size()-1, std::max(0,(int)currentState.iSegThres1-baseIdxRow));
   int idxVal_j = std::min((int)mStateValue[0].size()-1, std::max(0,(int)currentState.iDSegThres2));
   int idxVal_k = (int)currentState.VergenceState;

   tEyeDetectorState nextState = currentState;
    //Draw if we are exploring or Greedy
   if (drand48() < pExplore)
   {
       bExploreMove = true;
   }
   else
       bExploreMove = false;

    //bExploreMove = false; //Never Explore
    //Explore: Choose Action that determines Next State
    // Todo : Maybe policy could be modified to context, ie  EyeVergence
    //Jump to random threshold
    if (bExploreMove)
    {
        nextState.iSegThres1  =distrA(generator)+baseIdxRow;
        nextState.iDSegThres2 =distrB(generator);
    }

    //Policy: Take greedy action with some prob
    if (!bExploreMove)
    {

        ulong idxState_ActUp = std::min((int) mStateValue.size()-1,idxVal_i+1);
        ulong idxState_ActDown = std::max(0,idxVal_i-1);
        double qUp   = 0; //Marginal Value of going action up thresh.
        double qDown = 0;//Marginal Value of going action down thresh.
        double qStay = 0;

        //Calc Marginal Value For each Action taken (the state can shift the eye vergence so marginalize across a range )
        for (ulong v=0;v < mStateValue[idxVal_i][idxVal_j].size();v++ )
                qStay += mStateValue[idxVal_i][idxVal_j][v];

        for (ulong v=0;v < mStateValue[idxState_ActUp][idxVal_j].size();v++ )
                qUp += mStateValue[idxState_ActUp][idxVal_j][v];

        for (ulong v=0;v < mStateValue[idxState_ActDown][idxVal_j].size();v++ )
                qDown += mStateValue[idxState_ActDown][idxVal_j][v];

        if (qUp > qDown && qUp > qStay)
                nextState.iSegThres1++;

        if (qUp < qDown && qDown > qStay)
                nextState.iSegThres1--;

        // Check Limits: Do not exceed state space limits - bounce off walls;
        nextState.iSegThres1 = std::max((int)nextState.iSegThres1,baseIdxRow);
        nextState.iSegThres1 = std::min((int)nextState.iSegThres1,baseIdxRow+(int)mStateValue.size()-1);

        //Otherwise stay on the same Main threshold
        /// Calc action Value of Changing 2nd threshold //
        idxVal_i = (int)nextState.iSegThres1-baseIdxRow;
        qUp =qDown = 0.0;
        int idx2ndState_ActUp = std::min((int) mStateValue[0].size()-1,idxVal_j+1);
        int idx2ndState_ActDown = std::max(0,idxVal_j-1);
        // Action taken to move main thesh, then value of stay changed, recalc
        if (nextState.iSegThres1 != currentState.iSegThres1)
        {   qStay = 0.0;
            for (ulong v=0;v < mStateValue[idxVal_i][idxVal_j].size();v++ )
                    qStay += mStateValue[idxVal_i][idxVal_j][v];
        }

        for (ulong v=0;v < mStateValue[idxVal_i][idx2ndState_ActUp].size();v++ )
                qUp += mStateValue[idxVal_i][idx2ndState_ActUp][v];

        for (ulong v=0;v < mStateValue[idxVal_i][idx2ndState_ActDown].size();v++ )
                qDown += mStateValue[idxVal_i][idx2ndState_ActDown][v];

        if (qUp > qDown && qUp > qStay)
            nextState.iDSegThres2++;

        if (qUp < qDown && qDown > qStay)
            nextState.iDSegThres2--;

    } // Taking greedy action / not exploring

    //Enforce bounds on 2nd threshold
    nextState.iDSegThres2 = std::max((int)nextState.iDSegThres2,0);
    nextState.iDSegThres2= std::min((int)nextState.iDSegThres2,(int)mStateValue[0].size()-1);


    //Return next state to the environment - which will evaluate the fit score to update the value funct.
    return nextState;
}


/// \brief call it after drawing a new state and evaluating its fitness
/// Make transition to to NewState once value has been updated
/// \param RewardScore the immediate value of the current state and
double EyesDetector::UpdateStateValue(tEyeDetectorState toState,double RewardScore)
{

    int idxVal_i = std::max(0, std::min( (int)currentState.iSegThres1-(int)baseIdxRow, (int)mStateValue.size()-1 ));
    int idxVal_j = std::max(0,std::min( currentState.iDSegThres2,(int)mStateValue[0].size()-1) ); //The difference to 2nd threshold/ enfornce lower bound
    int idxVal_k = currentState.VergenceState;

    //Get Next states Value and update current - Greedy
    int idxNextVal_i = std::max(0, std::min( (int)toState.iSegThres1-(int)baseIdxRow, (int)mStateValue.size()-1));
    int idxNextVal_j = std::max(0,std::min( (int)toState.iDSegThres2, (int)mStateValue[0].size()-1) );
    int idxNextVal_k = (int)toState.VergenceState;

    //Update The Current State value by propagating value of next state backwards

    //TD learning - Use immediate Reward and add discounted future state rewards
    if (!bExploreMove) //propagate value from nearby threshold transitions
        mStateValue[idxVal_i][idxVal_j][idxVal_k] = mStateValue[idxVal_i][idxVal_j][idxVal_k] +
                                      alpha*(RewardScore + gamma*mStateValue[idxNextVal_i][idxNextVal_j][idxNextVal_k]-mStateValue[idxVal_i][idxVal_j][idxVal_k] );
    else {//When exploring we do not update based on the transition / as the action taken could not have been  chosen -
         // only account for user/reward in state
        mStateValue[idxVal_i][idxVal_j][idxVal_k] = mStateValue[idxVal_i][idxVal_j][idxVal_k] + alpha*(RewardScore - mStateValue[idxVal_i][idxVal_j][idxVal_k]);
    }

    currentState = toState; //Make transition to new state as set by the environment

    return mStateValue[idxVal_i][idxVal_j][idxVal_k];
}


EyesDetector::EyesDetector(int RangeValThres_min,int RangeValThres_max,int AngleVal_min,int AngleVal_max): EyesDetector(RangeValThres_max-RangeValThres_min,(RangeValThres_max-RangeValThres_min)/2, 4)
{
///Keep The min values so We can translate, locatre State  Values
baseIdxRow = RangeValThres_min;

//baseIdxCol = AngleVal_min;

    for(int i = 0;i < mStateValue.size();++i)
    {
          for(int j = 0;j < mStateValue[i].size();++j)
          {
              for(int k = 0;k < mStateValue[i][j].size();++k)
                 mStateValue[i][j][k] = 100; // start with high value on all - so we get to explore more
          }
    }

 /// Load Archived values if they Exists
 /// Load Saved Learned Behaviour
  assert(gsEyeDetectorFilename.size() > 0);
  qDebug() << "Load EyeDetector State:" << gsEyeDetectorFilename;
  std::ifstream is(gsEyeDetectorFilename.toStdString());
  if (is.is_open())
  {

    try
      {
        cereal::XMLInputArchive archive(is);
        archive(mStateValue); //Load State Value
      }catch (QString e)
      {
              qDebug() << "Failed to open RLEyeDetector Data file:" << e;
      }

  }

}

void EyesDetector::SaveState()
{
    //Save Learned Values to Disk
    std::ofstream os(gsEyeDetectorFilename.toStdString());
    cereal::XMLOutputArchive archive(os);
    this->serialize(archive); //Load State Value

    os.flush();
    //archive.~XMLOutputArchive(); //XMLArchive only flushes if destroyed
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

void EyesDetector::setCurrentState(tEyeDetectorState State)
{
    currentState = State;
}


