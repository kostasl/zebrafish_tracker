#include "foodmodel.h"
#include "ellipse_detect.h"
#include "config.h"
#include "zfttracks.h"

#include <string>
#include <QDebug>
//#include <QApplication>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include "larvatrack.h" //If included here it causes circular search if fishModel Defs.


static int lastFoodID = 0;


preyModel::preyModel()
{
lastFoodID++;
inactiveFrames = 0;
activeFrames = 0;
blobMatchScore = 0;

muTurn = 0.0; //Mean Turn Angle De
sigmaTurn = 15.0;
muPropulsion = 1.0;
sigmaPropulsion = 0.01;

this->ID = lastFoodID;
ROIID = 0;
isTargeted = false; //Saves Location To Data File When True
isNew = true;

}

preyModel::preyModel(zfdblob blob,zfdID ID):preyModel()
{
    inactiveFrames = 0;
    activeFrames = 0;
    blobMatchScore = 0;
    this->ID = ID;

    zTrack.id   = this->ID;
    zTrack.colour = CV_RGB(20,170,20);
    zTrack.centroid = blob.pt;

    //Initialize Filter Estimate to x0,y0
    ptEstimated = ptPredicted =  blob.pt;


    headingTheta = M_PI;
    velocity.x = 0.0;
    velocity.y = 0.0;
    df_propulsion =  muPropulsion+gsl_ran_gaussian(gTrackerState.r,sigmaPropulsion ); // g2prop.random(); //(float)random(-50,50)/1000.0f;
    df_friction = 0.75;  /// Water Drag Force Coefficient / (prop to velocity)
    omegaDeg = 0; //Angular Velocity In Degrees
    mass     = 1.0; //Assume unity mass initially

    zTrack.boundingBox.x = blob.pt.x - 6;
    zTrack.boundingBox.y = blob.pt.y - 6;
    zTrack.boundingBox.width = 12;
    zTrack.boundingBox.height = 12;

//    zTrack.boundingBox = cv::Rect(blob.pt.x - 5,blob.pt.y - 5,5,5);
    this->zfoodblob = blob;
}


preyModel::~preyModel()
{
    zTrack.pointStack.clear();
    zTrack.pointStackRender.clear();

}

/// \brief Models prey as a particle with propulsion force F and heading Theta moving in 2D space
///
void preyModel::predictMove()
{

//    //Draw Random Heading Change (Degrees)
//       dTheta = gTurn.random();//(float)random(-2,-2)/100.0;
//       //Update Angular Speed //With turn Drag
//       gTurn.setMean(5.0); //Reset Mean

//       //bounce off the border
//       if (current_position.X >= (matrix.width()-1)|| current_position.X <= 0)
//       {
//         velocity.X = -(df_propulsion/mass)*cos(headingTheta);
//         dTheta += 90;
//         gTurn.move(40);//Change the mean turning angle randomly
//         //gTurn.move(sgn(velocity.X)*0.001 ); //Change the mean turning angle randomly
//       }
//       if (current_position.Y >= (matrix.height()-1) || current_position.Y <= 0 )
//       {

//         velocity.Y = (df_propulsion/mass)*sin(headingTheta);
//         gTurn.move(40);//Change the mean turning angle randomly
//       }
//       omegaDeg +=dTheta-df_friction*omegaDeg*0.5; //Update Angular Speed //Assume A Turn Drag too
//       omegaDeg = max(-25,min(omegaDeg,25)); //Maximum Angular Velocity

//       //Update Heading in Rads
//       headingTheta += omegaDeg*(M_PI/180.0);//% (float)TWO_PI;
//       headingTheta = (float)((int)(headingTheta*100.0) % (int)(TWO_PI*100.0))/100.0;

//       // Draw random in propulsion (Random 2D walk)


//       //Calc velocity update with Drag Forces
//       velocity.X += max(-MAX_VELOCITY,min( (df_propulsion/mass)*cos(headingTheta)-df_friction*pow(velocity.X,1),MAX_VELOCITY) ) ;
//       velocity.Y += max(-MAX_VELOCITY,min( (df_propulsion/mass)*sin(headingTheta) -df_friction*pow(velocity.Y,1),MAX_VELOCITY) ) ;

//       current_position.X += velocity.X ;
//       current_position.Y += velocity.Y;


//       //Bound Position
//       current_position.X = max(0.0,min(current_position.X,matrix.width()-1));
//       current_position.Y = max(0.0,min(current_position.Y,matrix.height()-1));

}

/// \brief Implements a simple position filtering technique using a G-H Filter, such as tracking between subsequent frames is improved,
/// removing any sudden jumps and noise from the estimated prey position
///
cv::Point2f preyModel::alpha_beta_TrackingFilter_step(cv::Point2f blobPt)
{


    cv::Point2f res;
    double dResidual = 0.0;

    //Prediction Step
    ptPredicted.x <- ptEstimated.x +(dx*dt);
    ptPredicted.y <- ptEstimated.y +(dy*dt);

    //#Update step / X
    dResidual = blobPt.x -ptPredicted.x;
    dx = dx +h*(dResidual/dt);
    ptEstimated.x = ptPredicted.x + g*dResidual;

    //#Update step / Y
    dResidual = blobPt.y -ptPredicted.y;
    dy = dy +h*(dResidual/dt);
    ptEstimated.y = ptPredicted.y + g*dResidual;


    // Make New prey item position Estimate
    res.x =  ptEstimated.x;
    res.y =  ptEstimated.y;

    //##print(res)

    return(res);

}

void preyModel::updateState(zfdblob fblob,int Angle, cv::Point2f bcentre,unsigned int nFrame,int matchScore,float szradius)
{

    cv::Point2f ptFiltered = alpha_beta_TrackingFilter_step(bcentre);

    float fDistToNewPosition = cv::norm(ptFiltered,this->zTrack.centroid );
    if (fDistToNewPosition > gTrackerState.gMaxClusterRadiusFoodToBlob)
    {
         qDebug() << "Prey " << this->ID << " match too far d: " << fDistToNewPosition << " Mscore :" << matchScore;
         //return;
    }


    blobMatchScore = matchScore;
    nLastUpdateFrame = nFrame; //Set Last Update To Current Frame
    this->zfoodblob      = fblob;
    this->zTrack.pointStack.push_back(ptFiltered);
    this->zTrack.effectiveDisplacement = fDistToNewPosition;


    this->zTrack.centroid = ptFiltered;//fblob->pt; //Or Maybe bcentre
    this->blobRadius = szradius;
    zTrack.boundingBox.x = ptFiltered.x - 6;
    zTrack.boundingBox.y = ptFiltered.y - 6;
    zTrack.boundingBox.width = 12;
    zTrack.boundingBox.height = 12;
    //Establish stable initial phase before removing new flag
    if (activeFrames > gTrackerState.gcMinFoodModelActiveFrames && isNew)
    {
        isNew = false; //Having succeded to achieve n consec. active frames this food item is established
        inactiveFrames = 0; //Reset Counter Of inactive Frames
    }
    ///Trick - Update is called when fooditem has been matched, yet we use the
    /// the update to check if it has been inactive for too long- if found on next frame it will become active again
    //Although it may have been found here, it is still marked inactive until the next round
    isActive = (inactiveFrames < gTrackerState.gcMaxFoodModelInactiveFrames) && !isNew;



    ///Trick 2: Only mark as active if blob size is > 1 , otherwise we may be just tracking pixel flow
    /// Filter out activity based on optic flow only/where a blob cannot be seen/ but tracking a video pixel nontheless
    if (fblob.size > 1)
        inactiveFrames = 0; //Reset Counter Of inactive Frames
    else {
       inactiveFrames++;
    }
    this->zTrack.inactive = inactiveFrames;
    ///Optimization only Render Point If Displaced Enough from Last One
    if (this->zTrack.effectiveDisplacement > gTrackerState.gDisplacementThreshold)
    {
        this->zTrack.pointStackRender.push_back(ptFiltered);
        //this->zTrack.active++;
        //this->zTrack.inactive = 0;
    }else {
        //this->zTrack.inactive++;
    }


   static cv::Scalar colHighlight =  cv::Scalar(0,250,20,50); //Colour
   static cv::Scalar colPlain = cv::Scalar(255,255,0,30);

    if (isTargeted )
        zTrack.colour = CV_RGB(240,210,10); //YelloW For Tracking
    else
        zTrack.colour = colPlain; //CV_RGB(0,200,10);

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
        preyModel* pfood = v.at(i);

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
        preyModel* pfood = (preyModel*)v.at(i);

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
int preyModel::getActiveFoodCount(foodModels& vfoodmodels)
{
    int retNfood = 0;
    foodModels::iterator ft = vfoodmodels.begin();

    while (ft != vfoodmodels.end())
    {
        preyModel* pfood = ft->second;
        assert(pfood);

        // Render Food that has been on for A Min of Active frames / Skip unstable Detected Food Blob - Except If Food is being Tracked
        if (pfood->activeFrames < gTrackerState.gcMinFoodModelActiveFrames && (!pfood->isTargeted))
        {
            ++ft; //Item Is not Counted
            continue;
        }

        ++ft;
        retNfood++; //only count the rendered Food Items ie. Active Ones
    }

return retNfood;
}


/// Logic for when food item should be deleted
bool preyModel::isUnused()
{
    bool bLost =  (!this->isActive
                 && !this->isNew
                 || this->inactiveFrames > gTrackerState.gcMaxFoodModelInactiveFrames
                 || (this->activeFrames < gTrackerState.gcMinFoodModelActiveFrames && this->inactiveFrames > gTrackerState.gcMaxFoodModelInactiveFrames/2))
                 && (this->isTargeted == false);

            return bLost;
}


