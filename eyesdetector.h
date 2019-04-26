#ifndef EYESDETECTOR_H
#define EYESDETECTOR_H

#include <vector>
#include <math.h>
/// \class EyesDetector
/// \brief Uses reinforcement learning techniques to identify best way to discriminate eye angle.
///  Uses ellipsedetection score to estimate value of each action.
/// state is
///
/// can be converted to baysian sampler, to fit a model parameters

extern int gthresEyeSeg;

typedef unsigned long ulong;
typedef struct EyeDetectorState
{
    EyeDetectorState(){
        iSegThres1 = gthresEyeSeg;
        iSegThres2 = gthresEyeSeg+1;
        VergenceState = 1;
    }
    EyeDetectorState(int initSegT1,int initSegT2,int initVAngle):iSegThres1(initSegT1),iSegThres2(initSegT2){
        setVergenceState(initVAngle);
    }

    //Translate vergence angles to discretized eye vergence states
    void setVergenceState(double vangle){
                //Set A discrete number of eye vergence states - here 4 - ranging from -20 to 100
        if (vangle < 20.0)
            VergenceState = 1;
        else if (vangle >= 20.0 & vangle < 40.0)
            VergenceState = 2;
        else if (vangle >= 40.0 & vangle < 60.0)
            VergenceState = 3;
        else if (vangle >= 60.0)
            VergenceState = 4;

       if (VergenceState > 4)
           throw "Error:setVergenceState";


      //  VergenceState = (int)std::max(0.0, std::min(12.0, round( (vangle)/10.0 )+2.0 )  );
    }
    int iSegThres1;
    int iSegThres2;
    int VergenceState;



} tEyeDetectorState;



class EyesDetector
{
public:
    EyesDetector(int VarRows,int ValCol):mStateValue( VarRows,std::vector<double>(ValCol) ){} //Rows - number of segthres to explore - Cols - Vergence state descriptors
    EyesDetector(int RangeValThres_min,int RangeValThres_max,int AngleVal_min,int AngleVal_max);
    EyesDetector();
    tEyeDetectorState getCurrentState();
    ~EyesDetector();

    std::vector<double> & operator[](int i)
    {
      return mStateValue[i];
    }
    const std::vector<double> & operator[] (int i) const
    {
      return mStateValue[i];
    }
    void resize(int rows, int cols)//resize the two dimentional array .
    {
        mStateValue.resize(rows);
        for(int i = 0;i < rows;++i) mStateValue[i].resize(cols);
    }

    /// The TD learning transitions function
    /// \param RewardScore the fitness score as evaluated after fitting ellipse with given settings
    tEyeDetectorState DrawNextAction(tEyeDetectorState currentState);
    double UpdateStateValue(tEyeDetectorState nextState,double RewardScore);


private:
    const double pExplore = 0.1;
    const double alpha = 0.1; //learning rate for step update
    const double gamma = 0.1; //discounting of fut. rewards
    std::vector< std::vector<double> > mStateValue;
    ulong baseIdxRow; //Used to Translate State Value range to value table idxs
    ulong baseIdxCol;
    bool bExploreMove;
    tEyeDetectorState currentState;
};

#endif // EYESDETECTOR_H
