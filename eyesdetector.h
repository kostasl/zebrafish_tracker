#ifndef EYESDETECTOR_H
#define EYESDETECTOR_H

#include "config.h"
#include <vector>
#include <math.h>
#include <random>



#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include "cereal/types/vector.hpp"

#include <fstream>

#include "ellipse_detect.h"

/// \class EyesDetector
/// \brief Uses reinforcement learning techniques to identify best way to discriminate eye angle.
///  Uses ellipsedetection score to estimate value of each action.
/// state is
///
/// can be converted to baysian sampler, to fit a model parameters

extern int gthresEyeSeg;

typedef unsigned long ulong;
typedef std::vector<std::vector< std::vector<double> > > tStateValueMatrix;
typedef struct EyeDetectorState
{
    EyeDetectorState(){
        iSegThres1 = gthresEyeSeg;
        iDSegThres2 = 1;
        yEyePosition = 7;
        VergenceState = 1;
    }
    EyeDetectorState(int initSegT1,int initSegT2,int initVAngle):iSegThres1(initSegT1),iDSegThres2(initSegT2){
        setVergenceState(initVAngle);
        yEyePosition = 7;
    }

    //Translate vergence angles to discretized eye vergence states
    void setVergenceState(double vangle){
                //Set A discrete number of eye vergence states - here 4 - ranging from -20 to 100
        if (vangle < 20.0)
            VergenceState = 0;

        if (vangle >= 20.0 && vangle < 40.0)
            VergenceState = 1;

        if (vangle >= 40.0 && vangle < 60.0)
            VergenceState = 2;

        if (vangle >= 60.0)
            VergenceState = 3;

       if (VergenceState > 3)
           throw "Error:setVergenceState";


      //  VergenceState = (int)std::max(0.0, std::min(12.0, round( (vangle)/10.0 )+2.0 )  );
    }

    void setEyePos(tDetectedEllipsoid lEye,tDetectedEllipsoid rEye){

    }


    int iSegThres1;
    int iDSegThres2; //Distance Of T2 from T1
    int yEyePosition; //Centre on Y axis for Eye fitted ellipsoids - We can allow the learning algorithm to fix/move the head cutout to improve fit
    int VergenceState;



} tEyeDetectorState;



class EyesDetector
{
public:
    EyesDetector(ulong  ValD1,ulong  ValD2,ulong ValD3):mStateValue( ValD1,std::vector<std::vector<double>>(ValD2,std::vector<double>(ValD3)) ){} //Rows - number of segthres to explore - Cols - Vergence state descriptors
    EyesDetector(int RangeValThres_min,int RangeValThres_max,int AngleVal_min,int AngleVal_max);
    EyesDetector();
    tEyeDetectorState getCurrentState();
    void setCurrentState(tEyeDetectorState State);
    ~EyesDetector();
        // OVerload the member access functs
    double & operator[](int i)
    {

      return mStateValue[1][1][i];
    }

    /// \todo make 1 iterator design flattening the 3D array
    const double  operator[] (int i) const
    {
      return mStateValue[1][1][i];
    }

    /// overload size operator / return full state object size
    size_t
    size() const _GLIBCXX_NOEXCEPT
    { return mStateValue.size()*mStateValue[1].size()*mStateValue[1][1].size();}

    void resize(int rows, int cols,int index)//resize the 3 Dimensional array .
    {
        mStateValue.resize(rows);
        for(int i = 0;i < rows;++i)
        {
            mStateValue[i].resize(cols);
            for (int j=0;j<cols;j++)
                mStateValue[i][j].resize(index);
        }

    }
    // This method lets cereal know which data members to serialize
      template<class Archive>
      void serialize(Archive & archive)
      {
        archive(CEREAL_NVP(mStateValue)); // serialize things by passing them to the archive
      }


    /// The TD learning transitions function
    /// \param RewardScore the fitness score as evaluated after fitting ellipse with given settings
    tEyeDetectorState DrawNextAction(tEyeDetectorState currentState);
    double UpdateStateValue(tEyeDetectorState nextState,double RewardScore);
    void SaveState();/// Saves learned behaviour to disk in RLEyeDetector.xml;

private:
    const double pExplore = 0.1;
    const double alpha = 0.1; //lqearning rate for step update
    const double gamma = 0.2; //discounting of fut. rewards
    tStateValueMatrix mStateValue;
    int baseIdxRow; //Used to Translate State Value range to value table idxs
    int baseIdxCol;
    bool bExploreMove;
    tEyeDetectorState currentState;
    std::default_random_engine generator;

};

#endif // EYESDETECTOR_H
