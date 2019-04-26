#ifndef EYESDETECTOR_H
#define EYESDETECTOR_H

#include <vector>
/// \class EyesDetector
/// \brief Uses reinforcement learning techniques to identify best way to discriminate eye angle.
///  Uses ellipsedetection score to estimate value of each action.
/// state is
///
/// can be converted to baysian sampler, to fit a model parameters

typedef unsigned long ulong;

typedef struct
{
    int iSegThres1;
    int iSegThres2;
    int VergenceAngle;
} tEyeDetectorState;



class EyesDetector
{
public:
    EyesDetector(int VarRows,int ValCol):mStateValue( VarRows,std::vector<int>(ValCol) ){};
    EyesDetector(int RangeValThres_min,int RangeValThres_max,int AngleVal_min,int AngleVal_max);
    EyesDetector();
    ~EyesDetector();

    std::vector<int> & operator[](int i)
    {
      return mStateValue[i];
    }
    const std::vector<int> & operator[] (int i) const
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
    tEyeDetectorState DrawNextState(tEyeDetectorState currentState,double RewardScore);
    void UpdateStateValue(tEyeDetectorState nextState);


private:
    const double pExplore = 0.1;
    const double alpha = 0.1; //learning rate for step update
    const double gamma = 0.1; //discounting of fut. rewards
    std::vector< std::vector<int> > mStateValue;
    ulong baseIdxRow; //Used to Translate State Value range to value table idxs
    ulong baseIdxCol;
    bool bExploreMove;
    tEyeDetectorState currentState;
};

#endif // EYESDETECTOR_H
