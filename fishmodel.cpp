#include "fishmodel.h"
#include "ellipse_detect.h"
#include "config.h"

#include <opencv2/video/tracking.hpp>

extern cv::Mat frameDebugC;
extern cv::Size gszTemplateImg;
extern cv::Point gptTail,gptHead;

//extern double eyeStepIncrement;

//extern int gFishTailSpineSegmentLength;
//extern int gFitTailIntensityScanAngleDeg;
//extern const int gcFishContourSize; //Fixed number of fish Contour Points

//extern double gTemplateMatchThreshold;

fishModel::fishModel()
{
        bNewModel = true;
        stepUpdate = 1.0; // with fast rate and slow down with updates
        bearingAngle                = 0.0f;
        Delta_bearingAngle          = 0.0f; //last Change In bearing
        lastTailFitError            = 0.0;
        matchScore               = 0.0;
        nFailedEyeDetectionCount    = 0;

        inactiveFrames              = 0;
        matchScore               = 0;
        //coreTriangle.push_back(cv::Point());
        //coreTriangle.push_back(cv::Point());
        //coreTriangle.push_back(cv::Point());

        this->mouthPoint.x = 0;
        this->mouthPoint.y = 0;
        //this->leftEyeHull.clear();
        //this->rightEyeHull.clear();
        this->ID    = 0;
        zTrack.id   = this->ID;
        zTrack.colour = CV_RGB(255,0,0);
        leftEyeTheta          = 0; //In Degrees - A Value that looks wrong to show its not initialized
        rightEyeTheta         = 0; //In Degrees
        c_spineSegL           = gTrackerState.gFishTailSpineSegmentLength;


        //mState = cv::Mat::zeros(stateSize,1, type);

}


/////deprecated tracks And Blobs Here - To be removed
//fishModel::fishModel(cvb::CvTrack* track,cvb::CvBlob* blob):fishModel()
//{


//    this->ID        = track->id;
//    this->blobLabel = track->label;
//    //this->track     = track; //Copy Pointer
//    this->bearingRads = cvb::cvAngle(blob);
//    this->coreTriangle[2].x = track->centroid.x;
//    this->coreTriangle[2].y = track->centroid.y;

//    templateScore           = 0;
//    this->resetSpine();
//}

fishModel::fishModel(zftblob blob,int bestTemplateOrientation,cv::Point ptTemplateCenter):fishModel()
{
    stepUpdate = 1.0; // with fast rate and slow down with updates
    inactiveFrames  = 0;
    this->ID        = blob.hash() ;
    //this->blobLabel = blob.hash();

    zTrack.id       = this->ID;

    this->zfishBlob = blob; //Copy Localy
    //this->track     = NULL;
    this->bearingRads  = (float)bestTemplateOrientation*CV_PI/180.0;
    this->bearingAngle = (float)bestTemplateOrientation;
    this->ptRotCentre  = ptTemplateCenter;
    zTrack.centroid    = ptTemplateCenter;

    matchScore           = 0;
    this->resetSpine();

    qDebug() << "<<KF Init>>";
    // >>>> KF State Initialization
    // intialization of KF...
    KF.init(stateSize, measSize, contrSize, type);

    // KF State Vectors //
    mMeasurement = cv::Mat::zeros(measSize,1,type);
    mState =  cv::Mat::zeros(stateSize,1,type);

    cv::setIdentity(KF.errorCovPre, Scalar::all(1e-2f));
    cv::setIdentity(KF.errorCovPost, Scalar::all(1e-3f)); // default is 0, for smoothing try 0.1


    // [x,y,v_x,v_y,angle,angle_v]
    //qDebug() << "<<KF Set State>>";
    mState.at<float>(0) = (float)ptTemplateCenter.x; //X
    mState.at<float>(1) = (float)ptTemplateCenter.y; //Y
    mState.at<float>(2) = 0.0f; //speed X
    mState.at<float>(3) = 0.0f; //speed Y
    mState.at<float>(4) = 0.0f; //Angle Diff // (float)bestTemplateOrientation;// (Deg)
    mState.at<float>(5) = 0.0f; //Accel V Angle (Deg)
    mState.at<float>(6) = 0.0f; //Not Used
    mState.at<float>(7) = 10.0f;// Left Eye
    mState.at<float>(8) = -10.0f; //Right Eye


    // declare an array of floats to feed into Kalman Filter Transition Matrix, also known as State Transition Model
    // Transition State Matrix A  [x,y,v_x,v_y,angle,angle_v]
        // Note: set dT at each processing step!
        // [ 1 0 dT 0  0 0 ]
        // [ 0 1 0  dT 0 0 ]
        // [ 0 0 1  0  0 0 ]
        // [ 0 0 0  1  0 0 ]
        // [ 0 0 0  0  1 dT ]
        // [ 0 0 0  0  0 1 ]
    //qDebug() << "<<KF Init Transition M>>";
    KF.transitionMatrix = cv::Mat::zeros(measSize, stateSize, type);
    cv::setIdentity(KF.transitionMatrix,cv::Scalar::all(1.0f));

    // Measure Matrix H  [z_x, z_y, angle]
    // [ 1 0 0 0 0 0 ]
    // [ 0 1 0 0 0 0 ]
    // [ 0 0 0 0 1 0 ]
    KF.measurementMatrix = cv::Mat::zeros(measSize, stateSize, CV_32FC1);
    cv::setIdentity(KF.measurementMatrix,cv::Scalar::all(1.0f));

    mMeasurement.at<float>(0) = ptTemplateCenter.x;
    mMeasurement.at<float>(1) = ptTemplateCenter.y;
    mMeasurement.at<float>(2) = 0.0f; //speed X
    mMeasurement.at<float>(3) = 0.0f; //speed Y
    mMeasurement.at<float>(4) = 0.0f;// No change in angle initiallythis->bearingAngle;
    mMeasurement.at<float>(5) = 0.0f; //V Angle Accell
    mMeasurement.at<float>(6) = 0.0f; //Not Used
    mMeasurement.at<float>(7) = 10.0f;//this->leftEye.getEyeAngle();
    mMeasurement.at<float>(8) = -10.0f;//this->rightEye.getEyeAngle();
    //KF.measurementMatrix.at<float>(0) = 1.0f;
    //KF.measurementMatrix.at<float>(7) = 1.0f;
    //KF.measurementMatrix.at<float>(16) = 1.0f;


//    // Process Noise Covariance Matrix Q  [E_x,E_y, E_v_x,E_v_y ,E_angle,Eangle_v]
//        // [ Ex   0   0     0     0    0   0  0]
//        // [ 0    Ey  0     0     0    0   0  0]
//        // [ 0    0   Ev_x  0     0    0   0  0]
//        // [ 0    0   0     Ev_y  0    0   0  0]
//        // [ 0    0   0     0     Ea   0   0  0]
//        // [ 0    0   0     0     0    Ea_v0  0]
          // [ 0    0   0     0     0    0   lE 0]
          // [ 0    0   0     0     0    0   0  rE]
    cv::setIdentity(KF.processNoiseCov, cv::Scalar(6e-5)); // default is 1, for smoothing try 0.0001
    //Maybe Noise Suppression too high at 1e-3 , introduces lag in position
//    KF.processNoiseCov.at<float>(0,0) = 1e-2;
//    KF.processNoiseCov.at<float>(1,1) = 1e-2;
//    KF.processNoiseCov.at<float>(2,2) = 1e-3;
//    KF.processNoiseCov.at<float>(3,3) = 1e-3;
//    KF.processNoiseCov.at<float>(4,4) = 1e-4f; //Angle Change between frames
//    KF.processNoiseCov.at<float>(5,5) = 1e-4f; //Angular Accell (V of Diff)
//    KF.processNoiseCov.at<float>(4,5) = 1e-4f; //Angular Diff to Angle Accell
//    KF.processNoiseCov.at<float>(6,6) = 0.0f; //Not Used
    KF.processNoiseCov.at<float>(7,7) = 1e-4f; //Left Eye
    KF.processNoiseCov.at<float>(8,8) = 1e-4f; //Right Eye
//    //KF.processNoiseCov.at<float>(35) = 1e-1f;

//    // Measures Noise Covariance Matrix R - Set high/low so Filter Follows Measurement more closely
    cv::setIdentity(KF.measurementNoiseCov, cv::Scalar(1e-5f)); // default is 1, increasing should smooth but I get erratic behaviour
    //KF.measurementNoiseCov.at<float>(2,2) = 1e-3f; //Angular Diff (Speed)- per frame
    //KF.measurementNoiseCov.at<float>(3,3) = 1e-3f; //Angle Accell
    //KF.measurementNoiseCov.at<float>(1,2) = 1e-2f; // Y speed X pos
    //KF.measurementNoiseCov.at<float>(0,3) = 1e-2f; //X Speed - X
    //KF.measurementNoiseCov.at<float>(2,4) = 1e-5f; //Y Speed - Y
    //KF.measurementNoiseCov.at<float>(4,4) = 1e-2f; //Angular Diff (Speed)- per frame
    //KF.measurementNoiseCov.at<float>(5,5) = 1e-2f; //Angle Accell
    //KF.measurementNoiseCov.at<float>(4,5) = 1e-3f; //1e-1f; //Angular V_speed - Accell Covar
    //KF.measurementNoiseCov.at<float>(5,6) = 0;//1e-4f;//1e-1f; //Angular Accelleration-

    KF.measurementNoiseCov.at<float>(7,7) = 1e-3f; //Left Eye
    KF.measurementNoiseCov.at<float>(8,8) = 1e-3f; //Right Eye

    KF.statePre = mState;
    KF.statePost = mState.clone();

    //qDebug() << "Fish Model Construct.";
}

fishModel::~fishModel()
{
    //Clear The Vectors Contained
    this->zTrack.pointStack.clear();
    this->zTrack.pointStack.shrink_to_fit();
    this->zTrack.pointStackRender.clear();
    this->zTrack.pointStackRender.shrink_to_fit();

}

float fishModel::leftEyeAngle()
{

//    if (this->leftEyeRect.size.width < this->leftEyeRect.size.height)
//        return (this->leftEyeRect.angle-90.0)*CV_PI/180.0;
//    else
//        {
//         return (this->leftEyeRect.angle)*CV_PI/180.0;
//        }
//These Values Are Kalman Filtered
return leftEyeTheta; //leftEye.rectEllipse.angle;
}

/// \brief return (corrected for leading edge from horizontal line -Pi ... +Pi) Rectangle angle in Radians
float fishModel::rightEyeAngle()
{

//if (this->rightEyeRect.size.width < this->rightEyeRect.size.height)
//    return (this->rightEyeRect.angle-90.0)*CV_PI/180.0;
//else
//    {
//     return (this->rightEyeRect.angle)*CV_PI/180.0;
//    }
//These Values Are Kalman Filtered
return rightEyeTheta; //leftEye.rectEllipse.angle;

}

///
/// \brief fishModel::resetSpine make a straight Spline pointing towards Blobs Bearings Angle
///
void fishModel::resetSpine()
{
    //Reset Legth
    c_spineSegL =  gTrackerState.gc_FishTailSpineSegmentLength_init;

    this->spline.clear();
    spline.reserve(c_spinePoints);

    for (int i=0;i<c_spinePoints;i++)
    {

        splineKnotf sp;
        //1st Spine Is in Opposite Direction of Movement and We align 0 degrees to be upwards (vertical axis)
        //if (this->bearingRads > CV_PI)
        if (this->bearingRads < 0.0f)
            this->bearingRads += 2.0*CV_PI;

            sp.angleRad    = (this->bearingRads)-CV_PI ; //  //Spine Looks In Opposite Direction
            sp.spineSegLength = c_spineSegL;    //Default Size
            if (sp.angleRad < 0.0f)
                sp.angleRad += 2.0*CV_PI;
        //else
//            sp.angleRad    = (this->bearingRads)+CV_PI/2.0; //CV_PI/2 //Spine Looks In Opposite Direcyion

        assert(!std::isnan(sp.angleRad && std::abs(sp.angleRad) <= 2.0f*CV_PI && (sp.angleRad) >= 0 ));

        if (i==0)
        {
            sp.x =  this->ptRotCentre.x;
            sp.y =  this->ptRotCentre.y;
        }
        else
        {
            //0 Degrees Is vertical Axis Looking Up
            sp.x        = spline[i-1].x + ((double)c_spineSegL)*sin(sp.angleRad);
            sp.y        = spline[i-1].y - ((double)c_spineSegL)*cos(sp.angleRad);
        }

        spline.push_back(sp); //Add Knot to spline
    }



//    //    ///DEBUG
//        for (int j=0; j<c_spinePoints;j++) //Rectangle Eye
//        {
//            cv::circle(frameDebugC,cv::Point(spline[j].x,spline[j].y),2,TRACKER_COLOURMAP[j],1);
//        }

       // drawSpine(frameDebugC);
       // cv::waitKey(300);

}

///
/// \brief fishModel::getSpine Recalculates Point positions using Stored Knot Params (Angles)
/// /Calculates Spine Positions assumes initial point x0 y0 stored at 0 index of vector
void fishModel::calcSpline(t_fishspline& outspline)
{

    //this->spline.clear();
    //this->spline.push_back(this->coreTriangle[2]);
    double dspineSegL =  outspline[0].spineSegLength;

    for (int i=1;i<c_spinePoints;i++)
    {

       outspline[i].x = outspline[i-1].x + (dspineSegL)*sin(outspline[i-1].angleRad);
       outspline[i].y = outspline[i-1].y - (dspineSegL)*cos(outspline[i-1].angleRad);
       outspline[i].spineSegLength = dspineSegL;
       assert(!std::isnan(outspline[i].y) && !std::isnan(outspline[i].x));
       assert(!std::isnan(outspline[i].angleRad));
    }

}

///
/// \brief fishModel::getSpine VARIATION Make a Spine variation modifying a specific param value indicated by index
/// \param inspline vector Passed by value so the original is unchanged
/// \param outspline - The variational Spline
/// \return distance of variation in Config space
double fishModel::getdeltaSpline(t_fishspline inspline, t_fishspline& outspline,int idxparam,double sgn)
{
    const double dAngleStep = -sgn*CV_PI/80.0;
    double dvarSpineSeg = (double)inspline[0].spineSegLength;
    double ret = 0.0;
    outspline = inspline;

    //If idxparam = 1,2 then we are varying initial Spline Point x0, y0 params
    if (idxparam == 0)
    {
        ret = sgn*0.5;
        outspline[0].x -= ret;
    }else if (idxparam == 1)
    {
        ret = sgn*0.5;
        outspline[0].y -= ret;
    }// segment size,
    else if (idxparam == 2)
     {

      //  if (dvarSpineSeg < this->c_MaxSpineLengthLimit)
       ret = sgn*0.3;
      //  else
      //      ret = sgn*0.0000001; //Stop

        dvarSpineSeg += ret;
        outspline[0].spineSegLength = dvarSpineSeg;
    }
    else //Index > 2 is spine Angles
    {
        outspline[idxparam-3].angleRad += dAngleStep;// Angle variation for this theta
        ret = dAngleStep*dvarSpineSeg+cos(dAngleStep)*dvarSpineSeg; //rTheta
    }

    //Readjust xi,yi (In variational terms calc x_i = f(q1,q2,q3..), y_i = f(q1,q2,q3..)
    calcSpline(outspline);



return ret; //Return Distance In Q space
}

///
/// \brief fishModel::getSplineParams - Returns vector of Cspace spine coordinates
/// \param inspline
/// \param outparams
/// \return
///
void fishModel::getSplineParams(t_fishspline& inspline,std::vector<double>& outparams)
{
    outparams.clear();
    outparams.reserve(c_spineParamCnt);

    //Add x0 - yo
    //outparams[0] = inspline[0].x;
    //outparams[1] = inspline[0].y;
    outparams.push_back(inspline[0].x);
    outparams.push_back(inspline[0].y);
    outparams.push_back(inspline[0].spineSegLength); //Assume same Length Across Spine, given by 1st knot
    for (int i=0;i<(c_spinePoints);i++)
    {
        outparams.push_back(inspline[i].angleRad);
        //outparams[i+2] = inspline[i].angleRad;

    }
    //outparams[0] = inspline[0].x;
    //outparams[1] = inspline[0].y;
    //for (int i=0;i<(c_spinePoints);i++)
//        outparams[i+2] = inspline[i].angleRad;



}


/// \brief Modifies a Spline according to Cspace params
void fishModel::setSplineParams(t_fishspline& inspline,std::vector<double>& inparams)
{
    double dvarSpineSegLength;
    for (int i=0;i<(c_spineParamCnt);i++)
    {
        if (i==0)
            inspline[0].x = inparams[0];
        if (i==1)
            inspline[0].y = inparams[1];
        if (i==2) //Segment Size
        {
            dvarSpineSegLength = inparams[2];
            inspline[0].spineSegLength = dvarSpineSegLength;
        }
        if (i>2)
            inspline[i-3].angleRad = inparams[i]; //Param 3 is actually 1st spine knot's angle
    }

    //Readjust xi,yi (In variational terms calc x_i = f(q1,q2,q3..), y_i = f(q1,q2,q3..)
    calcSpline(inspline);

}

float fishModel::vergenceAngle()
{

}

/// \brief Implements spine Curve function by combining piecewise elements
cv::Point2f fishModel::getPointAlongSpline(float z,t_fishspline& pspline)
{
    //const float spineLength = this->c_spineSegL*pspline.size(); //The fitted Spine's lentgh is fixed
    const float spineSegLength = pspline[0].spineSegLength;
    const float spineLength = spineSegLength*pspline.size(); //The fitted Spine's lentgh is fixed

    int idx = z/spineSegLength; //Find knot index which is contains point
    double segLen = z - idx*spineSegLength;  //Modulo Find length input var along a linear segment

    if (idx > (pspline.size()-1)) //If Input Exceeds Spine Length
        return  cv::Point2f(pspline[pspline.size()-1].x,pspline[pspline.size()-1].y);

    ///Now construct point using Length along curve and return
    cv::Point2f ptC;
    ptC.x = pspline[idx].x + segLen*sin(pspline[idx].angleRad);
    ptC.y = pspline[idx].y - segLen*cos(pspline[idx].angleRad);

    return ptC;
}


///
/// \brief distancePointToSpline Finds Foot point - Important for measuring fit error between points and spline
/// Different Schemes exists - such as PDM, TDM - SDM - Start with Point Distance Radial Measure
/// \param pt
/// \param spline
/// \return
///
double fishModel::distancePointToSpline(cv::Point2f ptsrc,t_fishspline& pspline)
{

    //const int spineLength = this->c_spineSegL*pspline.size(); //The fitted Spine's lentgh is fixed
    const int spineLength = pspline[0].spineSegLength*pspline.size(); //The fitted Spine's varies
    const double dCStep = 0.3; //Step size on when searching along Spine Curve for closest Point (Foot Point )tk

    int idxNear = 0;
    //Take Spine  Point From Body / Set As foot point
    cv::Point2f ptFoot = cv::Point2f(pspline[idxNear].x,pspline[idxNear].y);
    //double distX = pow(ptsrc.x-ptFoot.x,2);
    //double distY = pow(ptsrc.y-ptFoot.y,2);
    double mindist;// = sqrt(distX + distY);
    mindist= cv::norm(ptsrc-ptFoot); //Start from Foot / Body Point look for point closer than this
    double dist = mindist;

    float fScanC = 0.0;
    ///Find Closest point on Curve to POint
    ///  \todo this is a crude/naive search -ok for small spine lengths
    /// Best to improve by calculating/estimating FootPoint projection
    while (fScanC < spineLength)
    {
        cv::Point2f  ptTest = getPointAlongSpline(fScanC,pspline);
        //distX = pow(ptsrc.x-ptTest.x,2);
        //distY = pow(ptsrc.y-ptTest.y,2);
        //dist = sqrt(distX + distY);
        dist= cv::norm(ptsrc-ptTest);
        //Check if distance minimized -
        if (dist < mindist)
        {
            mindist = dist;
            ptFoot  = ptTest;
        }

        fScanC += dCStep; //Move Along Curve
    }
#ifdef _ZTFDEBUG_
    //Show Foot Points
    cv::circle(frameDebugC,ptFoot,1,CV_RGB(10,10,255),1);
#endif
    return mindist;
}

/// \brief Uses detected ellipsoids to set fish's eye model state / using an incremental update
///\return total Score for fit
int fishModel::updateEyeMeasurement(tEllipsoids& vLeftEll,tEllipsoids& vRightEll)
{
    int retPerfScore = 0;
    double fleftEyeTheta = 0.0f;
    //int ileftEyeSamples = 0;
    double frightEyeTheta = 0.0f;
    //int irightEyeSamples = 0;


    // If we are stuck on same frame then estimate the unbiased empirical mean angle for each eye
    // use an incremental mean calculation
    if (uiFrameIterations > 1)
        stepUpdate = 1.0/std::min(200.0, (double)uiFrameIterations);


    tDetectedEllipsoid mleftEye(vLeftEll);
    tDetectedEllipsoid mrightEye(vRightEll);

    fleftEyeTheta = mleftEye.getEyeAngle();
    frightEyeTheta = mrightEye.getEyeAngle();
    //Incremental Update
    if (std::isnan(fleftEyeTheta) )
        this->nFailedEyeDetectionCount++;
    else
    {
        this->leftEyeTheta  = this->leftEyeTheta + stepUpdate*(fleftEyeTheta - this->leftEyeTheta );
        this->leftEye = mleftEye; //Measurement
    }

    if (std::isnan(frightEyeTheta) )
        this->nFailedEyeDetectionCount++;
    else{
        this->rightEyeTheta = this->rightEyeTheta + stepUpdate*(frightEyeTheta - this->rightEyeTheta );
        this->rightEye = mrightEye;
    }



    if (mleftEye.fitscore > 0 && mrightEye.fitscore > 0)
    {
       this->nFailedEyeDetectionCount = 0; // Reset Error Count
       retPerfScore = this->leftEye.fitscore + this->rightEye.fitscore;
    }else //penalize
    {
       retPerfScore =  (mleftEye.fitscore + mrightEye.fitscore)- 400;
    }

    //Reset Step size to default
     stepUpdate = gTrackerState.eyeStepIncrement;

 return (retPerfScore);


 //    //tDetectedEllipsoid
 //    // Go through All detected ellipsoids,
 //    for (int i=0; i< vLeftEll.size(); i++)
 //    {
 //        // select ones are for left/right eye
 //        tDetectedEllipsoid Eye = vLeftEll[i];
 //        if (Eye.cLabel == 'L') //Left eye
 //        {
 //            fleftEyeTheta += Eye.getEyeAngle();
 //            ileftEyeSamples +=1;
 //        }
 //    // and obtain mean angle for left/right eye from set of detected ellipsoids.
 //    }

 //    for (int i=0; i< vRightEll.size(); i++)
 //    {
 //        tDetectedEllipsoid REye = vRightEll[i];
 //        frightEyeTheta += REye.getEyeAngle();
 //        irightEyeSamples +=1;
 //        this->leftEye.fitscore += REye.fitscore;
 //    }

     //Get Mean sample angles
 //    if (ileftEyeSamples >0)
 //    {
 //        fleftEyeTheta   = fleftEyeTheta/(float)ileftEyeSamples;
 //        this->leftEye.fitscore = this->leftEye.fitscore/(float)ileftEyeSamples;
 //    }else
 //    {
 //        this->nFailedEyeDetectionCount++;
 //    }

 //    if (irightEyeSamples >0)
 //    {
 //       frightEyeTheta  = frightEyeTheta/(float)irightEyeSamples;
 //        this->rightEye.fitscore = this->rightEye.fitscore/(float)irightEyeSamples;
 //    }else
 //    {
 //        this->nFailedEyeDetectionCount++;
 //    }



//    if (vell.size() > 0)
//    {//Left Eye Detected First
//        tDetectedEllipsoid lEye = vell.at(0); //L Eye Is pushed 1st

//        // Update Internal Variable for Eye Angle //
//        // Use an incremental/ recent average rule
//        lEye.rectEllipse.angle = fleftEyeTheta;
//        this->leftEye           = lEye;

//        if (lEye.fitscore > 50)
//            this->leftEyeTheta = this->leftEyeTheta + stepUpdate*(fleftEyeTheta - this->leftEyeTheta );

//    }else
//    { //Set To Not detected - Do not update estimates - set score to 0
//        //this->leftEye       = tDetectedEllipsoid(cv::RotatedRect(),0);
//        //this->leftEyeTheta  = 180;
//        this->leftEye.fitscore = 0;
//        this->nFailedEyeDetectionCount++;
//    }


//   // ss.str(""); //Empty String
//    if (vell.size() > 1)
//    {
//      tDetectedEllipsoid rEye = vell.at(1); //R Eye Is pushed 2nd

//      frightEyeTheta     = rEye.rectEllipse.angle - 90;
//      //Fix Equivalent Angles To Range -50 +30
//      if (frightEyeTheta < -90)
//           frightEyeTheta      = rEye.rectEllipse.angle+90;
//      if (frightEyeTheta > 30)
//          frightEyeTheta       = rEye.rectEllipse.angle-90;

//      rEye.rectEllipse.angle = frightEyeTheta;
//      this->rightEye     = rEye; //Save last fitted ellipsoid struct

//      // Update Internal Variable for Eye Angle //
//      // Use an incremental/ recent average rule
//      if (rEye.fitscore > 50)
//      this->rightEyeTheta = this->rightEyeTheta + stepUpdate*(frightEyeTheta - this->rightEyeTheta );

//    }else
//    { //Set To Not detected
//     //   ss << "R Eye Detection Error - Check Threshold";
//     //   window_main.LogEvent(QString::fromStdString(ss.str()));

//        //this->rightEye       = tDetectedEllipsoid(cv::RotatedRect(),0);
//        //this->rightEyeTheta  = 180;
//        this->rightEye.fitscore = 0;
//        this->nFailedEyeDetectionCount++;
//    }

}


void fishModel::drawAnteriorBox(cv::Mat& frameScene, cv::Scalar colour=CV_RGB(00,00,255))
{
    ///Draw a Red Rotated Frame around Detected Body
    cv::Point2f boundBoxPnts[4];
    bodyRotBound.points(boundBoxPnts);
     for (int j=0; j<4;j++) //Rectangle Body
       cv::line(frameScene,boundBoxPnts[j],boundBoxPnts[(j+1)%4],colour,1,cv::LINE_8);
}

///
/// \brief fishModel::Update - Called On Every Frame Processed To Update Model State
/// The track point is set to the blob position and not the template centre
/// \param fblob
/// \param templatematchScore
/// \param Angle
/// \param bcentre
///

bool fishModel::stepPredict(unsigned int nFrame)
{

    double dT = (double)(nFrame-nLastUpdateFrame);///((double)gTrackerState.gfVidfps+1.0)

    //cout << "T:" << KF.transitionMatrix << endl;
    if (!bNewModel) // First detection!
    { //Add  Speed Contributions
        KF.transitionMatrix.at<float>(0,2) = 1.0;
        KF.transitionMatrix.at<float>(1,3) = 1.0;
        KF.transitionMatrix.at<float>(4,5) = 1.0f;//1e-9f;// 0.01; //Angular Accell Feeds into Angle Diff (Speed)
        KF.transitionMatrix.at<float>(5,6) = 0.0f; //1e-9f;//0.01; //Angular Speed Diff Feeds into Angular Speed

        mState = KF.predict();
        if (std::isnan(mState.at<float>(4)))
        {
            cerr << "[KalmanERROR] Sm:" << mState << endl;
            KF.init(stateSize, measSize, contrSize, type); //Re Init
            return(false);
        }

//        if (abs(mState.at<float>(4) - this->bearingAngle) > 20 )
//            qDebug() << "KF Angle prediction error:"  << this->bearingAngle << "->" << mState.at<float>(4);

        //Integrate
        //this->Delta_bearingAngle = mState.at<float>(4);
        this->bearingAngle   = (int)(zfishBlob.angle + mState.at<float>(4))%360; // Angle;

        this->bearingRads   =  this->bearingAngle*CV_PI/180.0;
        assert(!std::isnan(this->bearingAngle));
        this->ptRotCentre    = cv::Point2f(mState.at<float>(0), mState.at<float>(1)); //bcentre;
        this->leftEyeTheta   = mState.at<float>(7); // Eye Angle Left;
        this->rightEyeTheta   = mState.at<float>(8); // Eye Angle Right;


        bPredictedPosition = true;

        return (true);

    }

    return (false);
}

///
/// \brief fishModel::Update - Called On Every Frame Where Measurements from Blob are available - To Update Model State
/// The track point is set to the Kalman filtered blob state
///
/// \param fblob
/// \param matchScore - The fish blob's classifier score
/// \param Angle
/// \param bcentre
///
bool fishModel::updateState(zftblob* fblob, cv::Point2f bcentre,unsigned int nFrame,int SpineSegLength,int TemplRow, int TemplCol)
{

    //assert(!std::isnan(Angle));

    //Compare displacements to Last Measurements Not To Predicted Position In BearingAngle
    // Note Blob Angles can flip 180 between frames so we need to take the closest to the current orientation
    // \note Issue with Mod numbers Compass - Correcting towards from 355 to 0 the wrong way
    float angleDisplacement;
    float angleDisplacementA = getAngleDiff(zfishBlob.angle,this->bearingAngle);
    float angleDisplacementB = getAngleDiff((int)(zfishBlob.angle+180)%360,this->bearingAngle);
    //Choose the Displacement Closer to Current Angle (Fix Blob Noisy angle inversions)
    angleDisplacement = (abs(angleDisplacementA) < abs(angleDisplacementB))?angleDisplacementA:angleDisplacementB;
    //Angle = (abs(angleDisplacementA) < abs(angleDisplacementB))?Angle:(int)(Angle+180.0f)%360;


    double stepDisplacement = cv::norm(bcentre - this->zTrack.centroid);
    float dT = (float)(nFrame-nLastUpdateFrame);///((double)gTrackerState.gfVidfps+1.0)

    if (bNewModel)
        dT = 0;
    else
        if (angleDisplacement/dT > 150 ) //Correct Large Change - Blob Angle Flipped
        {
            qDebug() << "[W] Large angle change detected";
            //return (false);
        }

    if (!gTrackerState.bDraggingTemplateCentre)
        bUserDrag = false;

    //qDebug() << "Fish-Update M.";
    ///\note mod angles cannot be KF tracked as transitions 0->360 are non linear - instead I track the change in angle and integrate
    //Set to 1 frame minimum time step
    KF.transitionMatrix.at<float>(0,2) = dT;
    KF.transitionMatrix.at<float>(1,3) = dT;
    KF.transitionMatrix.at<float>(4,5) = dT;//0.0f; //  Angle Accell Feeds into Angle V
    KF.transitionMatrix.at<float>(5,6) = 0.0f;//dT;//dT; //Angular V Diff Feeds into Angle V
    //KF.transitionMatrix.at<float>(5,4) = 0;//dT;

    mMeasurement.at<float>(0) = bcentre.x;
    mMeasurement.at<float>(1) = bcentre.y;
    mMeasurement.at<float>(2) = mMeasurement.at<float>(3) = mMeasurement.at<float>(4) = mMeasurement.at<float>(5) = 0.0f;
    mMeasurement.at<float>(7) = this->leftEye.getEyeAngle();
    mMeasurement.at<float>(8) = this->rightEye.getEyeAngle();

    if (dT > 0){
        stepDisplacement = stepDisplacement/dT;
        angleDisplacement = angleDisplacement/dT;
        //Add Speed as measured Blob speed (Do not involve Filter Predictions in Measured Speed)
        mMeasurement.at<float>(2) = (bcentre.x-zfishBlob.pt.x)/dT;
        mMeasurement.at<float>(3) = (bcentre.y-zfishBlob.pt.y)/dT; //Y speed;
        mMeasurement.at<float>(4) = (float)angleDisplacement/dT;
        mMeasurement.at<float>(5) = Delta_bearingAngle - angleDisplacement;//Angle - this->bearingAngle; //angleDisplacement; //min(1.0f,max(-1.0f,(float)(angleDisplacement)/1.0f)); //(geAngleDiff(zfishBlob.angle,Angle)); //Ang Speed
        //mMeasurement.at<float>(6) = 0;//(floqat)(mState.at<float>(5) - (angleDisplacement))/2.0f; // min(0.1f,max(-0.1f,(float)(mState.at<float>(5) - (angleDisplacement))/2.0f)); //Ang Accelleration
    }
    //else
    //    mMeasurement.at<float>(2) = mMeasurement.at<float>(3) = mMeasurement.at<float>(5) = mMeasurement.at<float>(6) =0;
    // >>>> Matrix A -  Note: set dT at each processing step :

    /// Kalman Update - Measurements From Blob //
    /// generate measurement
    //mMeasurement += KF.measurementMatrix*mState;

    /// Kalman FILTER //
    ///  Re-Order - First adjust to measurement - then Predict
    //Reject Updates That Are Beyond Bounds
//    if (stepDisplacement > gTrackerState.gDisplacementLimitPerFrame ||
//        angleDisplacement > gTrackerState.gAngleChangeLimitPerFrame){
//        inactiveFrames++;
//    }else{ //Measurement valid - C0nsume
     mCorrected = KF.correct(mMeasurement); // Kalman Correction
     //mCorrected.copyTo(mState);
//     if (abs(mCorrected.at<float>(4) - Angle) > 20)
//         qDebug() << "KF Angle Meas:" << Angle <<  " Error Pred:" << mState.at<float>(4) << " Corrected " << mCorrected.at<float>(4);
     //Catch KALMAN error and skip update

     if (std::isnan(mCorrected.at<float>(0)))
         return(false);

     inactiveFrames  = 0;
 //   }


    this->zTrack.id     = ID;
    this->matchScore  = fblob->response;// templatematchScore;
    this->Delta_bearingAngle = mCorrected.at<float>(4);
    this->bearingAngle   =   zfishBlob.angle  + mCorrected.at<float>(4); // Integrate Angle Change onto Bearing;
    assert(!std::isnan(this->bearingAngle));
    this->bearingRads   =  this->bearingAngle*CV_PI/180.0;

    this->ptRotCentre    = cv::Point2f(mCorrected.at<float>(0), mCorrected.at<float>(1)); //bcentre;


    this->leftEyeTheta   = mCorrected.at<float>(6); // Eye Angle Left;
    this->rightEyeTheta   = mCorrected.at<float>(7); // Eye Angle Right;

    //Update The Eye Ellipsoids
    this->leftEye.rectEllipse.angle = this->leftEyeTheta;
    this->rightEye.rectEllipse.angle = this->rightEyeTheta;

    this->leftEye = tDetectedEllipsoid(this->leftEye.rectEllipse,this->leftEye.fitscore);
    this->rightEye = tDetectedEllipsoid(this->rightEye.rectEllipse,this->rightEye.fitscore);

    //Blob Position Is not FIltered
    this->zfishBlob      = *fblob;
    //this->c_spineSegL   = SpineSegLength;
    this->zTrack.pointStack.push_back(this->ptRotCentre);

    this->zTrack.effectiveDisplacement = cv::norm(this->ptRotCentre - this->zTrack.centroid);
    this->zTrack.centroid = this->ptRotCentre;//fblob->pt; //Or Maybe bcentre
    ///Optimization only Render Point If Displaced Enough from Last One
    if (this->zTrack.effectiveDisplacement > gTrackerState.gDisplacementThreshold)
    {
        this->zTrack.pointStackRender.push_back(this->ptRotCentre);
        this->zTrack.active++;
    }else {
        this->zTrack.inactive++;
    }

    /// Update Template Box Bound
    //int bestAngleinDeg = fish->bearingAngle;
    cv::RotatedRect fishRotAnteriorBox(fblob->pt,gTrackerState.gszTemplateImg ,this->bearingAngle); //this->ptRotCentre
    /// Save Anterior Bound
    this->bodyRotBound = fishRotAnteriorBox;

    this->idxTemplateCol = TemplCol;
    this->idxTemplateRow = TemplRow;

    //Set Spine Source to Rotation Centre
    this->spline[0].x       = this->ptRotCentre.x;
    this->spline[0].y       = this->ptRotCentre.y;
    //this->spline[0].angleRad   = this->bearingRads+CV_PI; //+180 Degrees so it looks in Opposite Direction


    //Check if frame advanced
    if (nLastUpdateFrame == nFrame)
        uiFrameIterations++; //Increment count of calculation cycles that we are stuck on same frame
    else
    {
        nLastUpdateFrame = nFrame; //Set Last Update To Current Frame
        uiFrameIterations = 0;
    }


    bNewModel = false; //Flag THat this model Has been now positioned
    bPredictedPosition = false; // position is based on corrected measurement

    return(true);

}//End of UpdateState


/// \brief Revised fitSpineContour , V2- Faster as it focuses/iterates around the spine points not the contour points
///  and uses the OpenCV pointPolygonTest to measure point distance to contour.
/// \param contours_body Fish Body Contour
/// \param frameImg_grey
double fishModel::fitSpineToContour2(cv::Mat& frameImg_grey, std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour)
{
    static const int cntParam = this->c_spineParamCnt;
    static const int cFitSpinePointsCount = ZTF_TAILSPINECOUNT;//

    //Parameter Found By Experience for current size fish
    ///Param sfish model should contain initial spline curve (Hold Last Frame Position)

    //Run Until Convergence Error is below threshold - Or Change is too small

    ///Compute Error terms for all data points/obtain local quadratic approx of fsd
    //For each contour Point
    std::vector<cv::Point> contour = contours_body[idxOuterContour];
    t_fishspline tmpspline = this->spline;
    t_fishspline dsSpline; //Variational Spline

    //Measure squared Distance error to closest Curve(spline) Point
        //Add to total error
    double dfitPtError_total = 10000.0;
    double dfitPtError_total_last = 0.0;
    double dDifffitPtError_total = 1000.0;

    double dTemp = 1.0; //Anealling Temperature
    /// \todo Optimize - Make Fish Contour Size Fixed - Then Allocate this as a buffer on the heap and reuse
    static double dJacobian[cFitSpinePointsCount][cntParam];//Vector of \nabla d for error functions
    memset(dJacobian,0.0,cFitSpinePointsCount*(cntParam)*sizeof(double));

    static double dGradf[cntParam];//Vector of Grad F per param
    memset(dGradf,0.0,cntParam*sizeof(double));

    static double dGradi[cntParam];//Vector of Grad Intensity per SPine POint param
    memset(dGradi,0.0,cntParam*sizeof(double));

    static double dResiduals[cFitSpinePointsCount];//Vector of \nabla d for error functions
    memset(dResiduals,0.0,cFitSpinePointsCount*sizeof(double));

    int cntpass     = 0;
    int cntStuck    = 0; //Number of Cycles solution has converge to an unnacceptable solution
    int cntSolved   = 0; //Number of Cycles Solution Is acceptable
    double dVarScale    = 1.0;
    //Do A number of Passes Before  Convergence //&& (dfitPtError_total/contour.size() > 8)
    while (cntpass < gTrackerState.gMaxFitIterations && (cntStuck < 5) && (cntSolved < 3) )
    {
        //Converged But Error Is still Large per Countour Point Then Jolt
        if (std::abs(dDifffitPtError_total) < 0.01 && dfitPtError_total/contour.size() > c_fitErrorPerContourPoint) //Time Out Convergece Count
        {
           cntStuck++;
           dVarScale = -dVarScale*1.0; //*1.2
        }
        else
        {
            cntStuck    = 0;
            dVarScale   = 1.0;
        }

        //Check For Ealy Convergence And Stop Early
        if (std::abs(dDifffitPtError_total) <= c_fitErrorPerContourPoint) //Time Out Convergece Count
        {
            cntSolved++;
            //dVarScale = dVarScale*0.93;
        }
        else
        {
            cntSolved   = 0;
            dVarScale   = 1.0;
        }


        //Reset Grad INfo - Start Pass From Last Point
        memset(dGradi,0.0,cntParam*sizeof(double));
        memset(dGradf,0.0,cntParam*sizeof(double));

        cntpass++;
        ///For Annealing
        dTemp = (double)cntpass/gTrackerState.gMaxFitIterations;
        //Prob Of Acceptance P = exp(-(sn-s)/T) >= drand
        ///
        dfitPtError_total_last  = dfitPtError_total;
        dfitPtError_total       = 0.0; //Reset

        double dq,ds; //Variation In Space And Score Variation
        //For Each Contour Point
       ///\todo invert the problem, go through each spine point and check against contour
        for (uint i=1;i<cFitSpinePointsCount;i+=1) //For Each Data point make a row in Jacobian
        {
            //dResiduals[i] = distancePointToSpline((cv::Point2f)contour[i],tmpspline);
            //Distance to closest Contour Edge, +ve inside, 0 on edge -ve outside
            // Invert so minimum is at centre of contour
            dResiduals[i] = -pointPolygonTest(contour, cv::Point2f(tmpspline[i].x,tmpspline[i].y), true );

            //Add Extra Grad Info/Cost for small segLength / thus pulling to longer spine length
            //dResiduals[i] -= (1.0)*c_MaxSpineLengthLimit/tmpspline[i].spineSegLength;

            // Push to get Tail Pt Coincide with Last spine point
            double distToTailTip = cv::norm( gptTail - cv::Point(tmpspline[tmpspline.size()-1].x,tmpspline[tmpspline.size()-1].y) );
            dResiduals[i] -= (0.01)*distToTailTip ;
           // double penalty = dResiduals[i]*0.10; //Calc Scaled Penalty
           // for (int s=0;s<tmpspline.size();s++)
           // {
           //     int pptTest = pointPolygonTest(contour, cv::Point2f(tmpspline[s].x,tmpspline[s].y), false );
           //     if (pptTest < 0 ) //if spine point is outside contour then increase residuals
           //         dResiduals[i] += penalty;
            //}
            dfitPtError_total       +=dResiduals[i];

            //Add Variation dx to each param and calc derivative
            //Start from param idx 2 thus skipping the 1st point(root ) position and only do angle variations
            for (int k=2;k < cntParam; k++)
            {   /// \note using only +ve dx variations and not -dx - In this C space Ds magnitude should be symmetrical to dq anyway

                dq = getdeltaSpline(tmpspline,dsSpline,k,dVarScale); //Return param variation
                // dsSpline residual of variation spline
                ds = -pointPolygonTest(contour, cv::Point2f(dsSpline[i].x,dsSpline[i].y), true );
                {
                    dJacobian[i][k] = (ds-dResiduals[i])/(dq);
                    //Got towards smaller Distance
                    dGradf[k]           += dResiduals[i]*dJacobian[i][k]; //Error Grad - gives Gradient in Cspace vars to Total error
                }
            }

        }//Loop Through All Contour (Data points)

//        //Now Using Intensity Specific Algorithm
//        ///Add Gradient Of Intensity - GradNow - GradVs - Spine Point Struct indexes Range from 0 To SpinePointCount
//        calcSpline(tmpspline);
//        for (int k=1;k < cFitSpinePointsCount; k++)
//        {
//            //Get Variation Of Spine Point
//            dq = getdeltaSpline(tmpspline,dsSpline,k+2,0.8); //Return Angle param variation
//            calcSpline(dsSpline); //Calc Position Of Next Spine Point Based On Angle Of Preceding one

//            float pxi0 = frameImg_grey.at<uchar>(cv::Point(tmpspline[k].x,tmpspline[k].y));
//            float pxi1 = frameImg_grey.at<uchar>(cv::Point(dsSpline[k].x,dsSpline[k].y));
//            dGradi[k]           += (pxi1 - pxi0)/dq;
/////       double dsi = max(1.0,cv::norm(cv::Point(tmpspline[k-2].x,tmpspline[k-2].y)-cv::Point(dsSpline[k-2].x,dsSpline[k-2].y)));
//        }

        std::vector<double> cparams(c_spineParamCnt);
        //Recover Params from tmp Spline
        getSplineParams(tmpspline,cparams);
        assert(cparams.size() == c_spineParamCnt);

        //Pass and modify CSpace Params with gradient descent
        cparams[2] -= 0.001*dGradf[2]; //- 0.0001*dGradi[2];
        for (int i=3;i<cntParam;i++)
        {   //Go Down Distance to Contour And Up Intensity Gradient
            cparams[i] += 0.01*dGradf[i];// -  0.01*dGradi[i];

#ifdef _ZTFDEBUG_
            qDebug() << "lamda GradF_"<< i << "-:" << 0.01*dGradf[i] << " GradI:" << 0.005*dGradi[i];
#endif
        }
        // Adapt spline Using the modified Params that follow the contour distance gradient
        setSplineParams(tmpspline,cparams);

        dDifffitPtError_total = dfitPtError_total - dfitPtError_total_last; //Total Residual /Error Measure Change

        //if (dfitPtError_total > 1000)
            //this->resetSpine(); //Start over


    }//While Error Change Is larger Than

#ifdef _ZTFDEBUG_
    qDebug() << "ID:" <<  this->ID << cntpass << " EChange:" << dDifffitPtError_total;
#endif

        // Copy Fitted Spline On Fish Spline
        this->spline            = tmpspline;

        // Copy Fitted TailSeg Length - But Impose Limits on Tail Length (in case contour was too poor and tail is too small)
        this->c_spineSegL       = std::max(c_MinSpineLengthLimit, std::min(c_MaxSpineLengthLimit, (double)tmpspline[0].spineSegLength));
        this->lastTailFitError = dfitPtError_total/c_spinePoints;


#ifdef _ZTFDEBUG_
        qDebug() << "Converged in n: " << cntpass;
#endif


///  DEBUG ///
#ifdef _ZTFDEBUG_
    for (int j=0; j<c_spinePoints;j++) //Rectangle Eye
    {
        cv::circle(frameDebugC,cv::Point(spline[j].x,spline[j].y),2,TRACKER_COLOURMAP[j],1);
    }
    cv::drawContours(frameDebugC,contours_body,idxOuterContour,CV_RGB(200,20,20),1);

#endif
///    End Debug ///

    return dfitPtError_total; //Return Total Fit Error
   // qDebug() << "D err:" << dDifffitPtError_total;
}

bool fishModel::isFrameUpdated(uint nFrame)
{
    return (nLastUpdateFrame == nFrame);
}

/// \brief Check if Model Position can be still considered to actively represent a fish's location
bool fishModel::isValid()
{
    //templateScore >= gTrackerState.gTemplateMatchThreshold // inactiveFrames < 2

    return(this->zfishBlob.response >= gTrackerState.fishnet_classifier_thres &&
           inactiveFrames == 0 &&  pointIsInROI(ptRotCentre, bodyRotBound.size.width) ||
           binFocus );


}
void fishModel::drawBodyTemplateBounds(cv::Mat& outframe)
{

    int bestAngleinDeg = this->bearingAngle;
    //cv::RotatedRect fishRotAnteriorBox(centre, cv::Size(gLastfishimg_template.cols,gLastfishimg_template.rows),bestAngleinDeg);

//    stringstream strLbl;
//    strLbl << "A: " << bestAngleinDeg;

    QString strlbl(QString::number(ID));
    cv::Scalar colour;
    colour = CV_RGB(60,60,60); //Grey - Means Not Validated/Fish Detected Region

    if (this->isValid())
    {
        colour = CV_RGB(250,250,0); //Yellow Measured-Filtered Position
        if (bPredictedPosition)
            colour = CV_RGB(50,50,210); //Blue Predicted position
    }




    //if ()
    //    colour = CV_RGB(50,250,0); //Green - FishNet Classified


#ifdef _ZTFDEBUG_
    QString strlbl("A: " + QString::number(bestAngleinDeg));
    cv::putText(outframe,strlbl.toStdString(),this->bodyRotBound.boundingRect().br()+cv::Point(-10,15),CV_FONT_NORMAL,0.4,colour,1);
#endif

    ///Draw a center Point
    cv::circle(outframe,this->bodyRotBound.center,5,colour,1);

    ///Draw a Red Rotated Frame around Detected Body
    drawAnteriorBox(outframe,colour);

    cv::putText(outframe,strlbl.toStdString(),this->bodyRotBound.boundingRect().br()+cv::Point(-10,15),CV_FONT_NORMAL,0.4,colour,1);

}

///
/// \brief fishModel::fitSpineToIntensity inspired by  Giovanni's tail fitting method : starting from initial point on the fish body and an initial direction for the tail it searches for the
///   highest intensity within a angle range -c_tailscanAngle +c_tailscanAngle degrees on a largely blurred scene image , through each spine point incrementally.
///  \note (15/1/18) Modified such that it updates to centre of Mass of sampled points and not to Highest intensity - as this HInt is noisy and ends up curling around the body frequently
///  It pretty fast compared to the fitContour Variational method
/// \param src
/// \param start
/// \param tgt_start Is a vector from point Start towards the initial guess of tail direction
/// \param step_size Segment Lenght between  Spine anchor points
/// \param anchor_pts // The Spine
///
void fishModel::fitSpineToIntensity(cv::Mat &frameimg_Blur,int c_tailscanAngle){
    const size_t AP_N= this->c_spinePoints;
    const int step_size = this->c_spineSegL;

    //const int c_tailscanAngle = gFitTailIntensityScanAngleDeg;
    if (!isValid())
        return; //Do not Calculate For Inactive Fish

    uint pxValMax;

    int angle; //In Deg of Where The spline point is looking towards - Used by Ellipse Arc Drawing

    //cv::Mat frameimg_Blur;
    //cv::GaussianBlur(imgframeIn,frameimg_Blur,cv::Size(5,5),5,5);
    //cv::imshow("IntensitTailFit",frameimg_Blur);

    std::vector<cv::Point> ellipse_pts; //Holds the Drawn Arc Points around the last spine Point

    for(unsigned int k=1;k<AP_N;k++)
    { //Loop Through Spine Points
        ellipse_pts.clear();


        angle = spline[k-1].angleRad/CV_PI*180.0-90.0; //Get Angle In Degrees for Arc Drawing Tranlated Back to 0 horizontal
        //Construct Elliptical Circle around last Spine Point - of Radius step_size
        cv::ellipse2Poly(cv::Point(spline[k-1].x,spline[k-1].y), cv::Size(step_size,step_size), 0, angle-c_tailscanAngle, angle+c_tailscanAngle, 2, ellipse_pts);

        if (ellipse_pts.size() ==0)
        {
            qDebug() << "fitSpineToIntensity: Failed empty ellipse2Poly";
            continue;
        }
        ///Calculate Moment of inertia Sum m theta along arc
        pxValMax                = 0;
        uint iTailArcMoment     = 0;
        uint iPxIntensity       = 0;
        uint iSumPxIntensity    = 1;
        // Loop Integrate Mass //
        for(int idx=0;idx<(int)ellipse_pts.size();idx++){
            //Obtain Value From Image at Point on Arc - Boundit in case it goes outside image
            int x = std::min(frameimg_Blur.cols,std::max(1,ellipse_pts[idx].x));
            int y = std::min(frameimg_Blur.rows,std::max(1,ellipse_pts[idx].y));
            iPxIntensity = frameimg_Blur.at<uchar>(cv::Point(x,y));

            //Use idx As Angle /Position
            iTailArcMoment  += idx*iPxIntensity;
            iSumPxIntensity += iPxIntensity;
        } //Loop Through Arc Sample Points

        //Update Spline to COM (Centre Of Mass) And Set As New Spline Point
        uint comIdx = iTailArcMoment/iSumPxIntensity;
        spline[k].x     = ellipse_pts[comIdx].x;
        spline[k].y     = ellipse_pts[comIdx].y;
        /// Get Arc tan and Translate back to 0 being the Vertical Axis
        if (k==1) //1st point Always points in the opposite direction of the body
            spline[k-1].angleRad    = (this->bearingRads)-CV_PI ; //  //Spine Looks In Opposite Direction
        else
            spline[k-1].angleRad = std::atan2(spline[k].y-spline[k-1].y,spline[k].x-spline[k-1].x)+CV_PI/2.0; // ReCalc Angle in 0 - 2PI range Of previous Spline POint to this New One

        //Set Next point Angle To follow this one - Otherwise Large deviation Spline
        if (k < AP_N)
            spline[k].angleRad = spline[k-1].angleRad;

        //Constrain Large Deviations
        if (std::abs(spline[k-1].angleRad - spline[k].angleRad) > CV_PI/2.0)
           spline[k].angleRad = spline[k-1].angleRad; //Spine Curvature by Initializing next spine point Constraint Next

    }

}



void fishModel::drawSpine(cv::Mat& outFrame)
{
    for (int j=0; j<c_spinePoints-1;j++) //Rectangle Eye
    {
        if (isValid())
           cv::circle(outFrame,cv::Point(spline[j].x,spline[j].y),2,TRACKER_COLOURMAP[j],1);
        else
           cv::circle(outFrame,cv::Point(spline[j].x,spline[j].y),2,CV_RGB(100,100,100),1);

        //Connect Spine Points
        //cv::line(outFrame,cv::Point(spline[j].x,spline[j].y),cv::Point(spline[j+1].x,spline[j+1].y),TRACKER_COLOURMAP[0],1);

            //cv::line(outFrame,cv::Point(spline[j].x,spline[j].y),cv::Point(spline[j+1].x,spline[j+1].y),TRACKER_COLOURMAP[0],1);
//        else
//        { //Draw Terminal (hidden) point - which is not a spine knot
//            cv::Point ptTerm;
//            ptTerm.x = spline[j].x + ((double)c_spineSegL)*sin(spline[j].angleRad);
//            ptTerm.y = spline[j].y - ((double)c_spineSegL)*cos(spline[j].angleRad);

//            cv::line(outFrame,cv::Point(spline[j].x,spline[j].y),ptTerm,TRACKER_COLOURMAP[0],1);
//        }
    }
    //Draw Final Tail (imaginary Spine) Point
    if (isValid())
        cv::circle(outFrame,cv::Point(spline[c_spinePoints-1].x,spline[c_spinePoints-1].y),2,TRACKER_COLOURMAP[c_spinePoints-1],1);
    else
        cv::circle(outFrame,cv::Point(spline[c_spinePoints-1].x,spline[c_spinePoints-1].y),2,CV_RGB(100,100,100),1);

}





///
/// \brief operator << //Overloaded Stream Operator
/// Output fishModel State to Log File
/// \param out
/// \param h
/// \return
///
std::ostream& operator<<(std::ostream& out, const fishModel& h)
{

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)


    out << h.nLastUpdateFrame
        << "\t" << h.ID
        << "\t" << h.bearingAngle
        << "\t" << h.zTrack
        << "\t" << h.leftEyeTheta
        << "\t" << h.rightEyeTheta
        << "\t" << h.matchScore
        << "\t" << h.lastTailFitError
        << "\t" << h.leftEye.fitscore
        << "\t" << h.rightEye.fitscore
        << "\t" << h.nFailedEyeDetectionCount;

    return out;
}


///
/// \brief operator << //Overloaded Stream Operator
/// Output fishModel State to Log File
/// \param out
/// \param h
/// \return
///
QTextStream& operator<<(QTextStream& out, const fishModel& h)
{
    const float Rad2Deg = (180.0/CV_PI);

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)
    out.setRealNumberNotation(QTextStream::RealNumberNotation::FixedNotation );
    out.setRealNumberPrecision(2);
   // assert(!std::isnan( h.leftEyeTheta ) && !std::isnan( h.rightEyeTheta ) );
    out << h.ID <<"\t"<< h.bearingAngle <<"\t" << h.zTrack << "\t" << h.leftEyeTheta << "\t" <<  h.rightEyeTheta;

    //Set Global 1st Spine Direction (Helps to detect Errors)
    assert(h.spline.size() > 0);

    out << "\t" << Rad2Deg* h.spline[0].angleRad;
    //Output Spine Point Angular Deviations from the previous spine/tail Segment in Degrees
    for (int i=1;i<h.c_spinePoints;i++)
    {

       out << "\t" << Rad2Deg*( h.spline[i-1].angleRad - h.spline[i].angleRad);

    }
     out << "\t" << h.c_spineSegL;
     out << "\t" << h.matchScore;
     out << "\t" << h.lastTailFitError;
     //if (h.leftEye)
     out << "\t" << h.leftEye.fitscore;
     //if (h.rightEye)
        out << "\t" << h.rightEye.fitscore;
     out << "\t" << h.nFailedEyeDetectionCount;

    return out;
}










///
/// \brief fishModel::fitSpineToContour / Least squares fit of constrained spline model with params x0,yo, and theta_i-n
///
/// \param contours_body
/// \param idxInnerContour
/// \param idxOuterContour
/// \warning uses arbitrary constants for detecting fit error or convergence
/// \return fitness error score
///
double fishModel::fitSpineToContour(cv::Mat& frameImg_grey, std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour)
{
    static const int cntParam = this->c_spineParamCnt;
    static const int gcFishContourSize = ZTF_FISHCONTOURSIZE+1;//Fixed PLus 1 tail Point

    //Parameter Found By Experience for current size fish
    ///Param sfish model should contain initial spline curve (Hold Last Frame Position)

    //Run Until Convergence Error is below threshold - Or Change is too small

    assert(contours_body.size() >= idxOuterContour && contours_body[idxOuterContour].size() == gcFishContourSize);
    ///Compute Error terms for all data points/obtain local quadratic approx of fsd
    //For each contour Point
    std::vector<cv::Point> contour = contours_body[idxOuterContour];
    t_fishspline tmpspline = this->spline;
    t_fishspline dsSpline; //Variational Spline

    //Measure squared Distance error to closest Curve(spline) Point
        //Add to total error
    double dfitPtError_total = 10000.0;
    double dfitPtError_total_last = 0.0;
    double dDifffitPtError_total = 1000.0;

    double dTemp = 1.0; //Anealling Temperature
    /// \todo Optimize - Make Fish Contour Size Fixed - Then Allocate this as a buffer on the heap and reuse
    static double dJacobian[gcFishContourSize][cntParam];//Vector of \nabla d for error functions
    memset(dJacobian,0.0,gcFishContourSize*(cntParam)*sizeof(double));

    static double dGradf[cntParam];//Vector of Grad F per param
    memset(dGradf,0.0,cntParam*sizeof(double));

    static double dGradi[cntParam];//Vector of Grad Intensity per SPine POint param
    memset(dGradi,0.0,cntParam*sizeof(double));

    static double dResiduals[gcFishContourSize];//Vector of \nabla d for error functions
    memset(dResiduals,0.0,gcFishContourSize*sizeof(double));

    int cntpass     = 0;
    int cntStuck    = 0; //Number of Cycles solution has converge to an unnacceptable solution
    int cntSolved   = 0; //Number of Cycles Solution Is acceptable
    double dVarScale    = 1.0;
    //Do A number of Passes Before  Convergence //&& (dfitPtError_total/contour.size() > 8)
    while (cntpass < gTrackerState.gMaxFitIterations && (cntStuck < 5) && (cntSolved < 3) )
    {
        //Converged But Error Is still Large per Countour Point Then Jolt
        if (std::abs(dDifffitPtError_total) < 0.01 && dfitPtError_total/contour.size() > 10) //Time Out Convergece Count
        {
            cntStuck++;
           dVarScale = -dVarScale*1.0; //*1.2
        }
        else
        {
            cntStuck    = 0;
            dVarScale   = 1.0;
        }

        //Check For Ealy Convergence And Stop Early
        if (std::abs(dDifffitPtError_total) < 0.06 && dfitPtError_total/contour.size() <= 8) //Time Out Convergece Count
        {
            cntSolved++;
            //dVarScale = dVarScale*0.93;
        }
        else
        {
            cntSolved   = 0;
            dVarScale   = 1.0;
        }


        //Reset Grad INfo - Start Pass From Last Point
        memset(dGradi,0.0,cntParam*sizeof(double));
        memset(dGradf,0.0,cntParam*sizeof(double));

        cntpass++;
        ///For Annealing
        dTemp = (double)cntpass/gTrackerState.gMaxFitIterations;
        //Prob Of Acceptance P = exp(-(sn-s)/T) >= drand
        ///
        dfitPtError_total_last  = dfitPtError_total;
        dfitPtError_total       = 0.0; //Reset

        double dq,ds; //Variation In Space And Score Variation
        //For Each Contour Point
       ///\todo invert the problem, go through each spine point and check against contour
        for (uint i=0;i<gcFishContourSize;i+=1) //For Each Data point make a row in Jacobian
        {
            dResiduals[i] = distancePointToSpline((cv::Point2f)contour[i],tmpspline);
           // double penalty = dResiduals[i]*0.10; //Calc Scaled Penalty
           // for (int s=0;s<tmpspline.size();s++)
           // {
           //     int pptTest = pointPolygonTest(contour, cv::Point2f(tmpspline[s].x,tmpspline[s].y), false );
           //     if (pptTest < 0 ) //if spine point is outside contour then increase residuals
           //         dResiduals[i] += penalty;
            //}
            dfitPtError_total       +=dResiduals[i];

            //Add Variation dx to each param and calc derivative
            //Start from param idx 2 thus skipping the 1st point(root ) position and only do angle variations
            for (int k=2;k < cntParam; k++)
            {   /// \note using only +ve dx variations and not -dx - In this C space Ds magnitude should be symmetrical to dq anyway

                dq = getdeltaSpline(tmpspline,dsSpline,k,dVarScale); //Return param variation
                // dsSpline residual of variation spline
                ds = distancePointToSpline((cv::Point2f)contour[i],dsSpline);
                //dsSpline.clear();
                //getdeltaSpline(tmpspline,dsSpline,k,-0.25) ; //add dx
                //ds += distancePointToSpline((cv::Point2f)contour[i],dsSpline); // Add df dsSpline residual of variation spline

                //if (dq > 0.0)
                {
                    dJacobian[i][k] = (ds-dResiduals[i])/(dq);
                    //Got towards smaller Distance
                    dGradf[k]           += dResiduals[i]*dJacobian[i][k]; //Error Grad - gives Gradient in Cspace vars to Total error
                }


            }

        }//Loop Through All Contour (Data points)

        //Now Using Intensity Specific Algorithm
//        ///Add Gradient Of Intensity - GradNow - GradVs - Spine Point Struct indexes Range from 0 To SpinePointCount
//        for (int k=2;k < cntParam; k++)
//        {

//            float pxi0 = frameImg_grey.at<uchar>(cv::Point(tmpspline[k-2].x,tmpspline[k-2].y));
//            float pxi1 = frameImg_grey.at<uchar>(cv::Point(dsSpline[k-2].x,dsSpline[k-2].y));
//            double dsi =std::max(1.0,cv::norm(cv::Point(tmpspline[k-2].x,tmpspline[k-2].y)-cv::Point(dsSpline[k-2].x,dsSpline[k-2].y)));
//            //Go Towards Higher INtensity Pixels
//            dGradi[k]           += (pxi1 - pxi0)/dq;
//        }

        std::vector<double> cparams(c_spineParamCnt);
        getSplineParams(tmpspline,cparams);
        assert(cparams.size() == c_spineParamCnt);

        ///modify CSpace Params with gradient descent
        //cparams[2] -= 0.001*dGradf[2] - 0.01*dGradi[2];
        for (int i=0;i<cntParam;i++)
        {   //Go Down Distance to Contour And Up Intensity Gradient
            cparams[i] -= 0.01*dGradf[i] - 0.01*dGradi[i];
#ifdef _ZTFDEBUG_
            qDebug() << "lamda GradF_"<< i << "-:" << 0.01*dGradf[i] << " GradI:" << 0.005*dGradi[i];
#endif
        }
        ///Modify Spline - ie move closer
        setSplineParams(tmpspline,cparams);


        dDifffitPtError_total = dfitPtError_total - dfitPtError_total_last; //Total Residual /Error Measure Change

        //if (dfitPtError_total > 1000)
            //this->resetSpine(); //Start over


    }//While Error Change Is larger Than

#ifdef _ZTFDEBUG_
    qDebug() << "ID:" <<  this->ID << cntpass << " EChange:" << dDifffitPtError_total;
#endif



        this->spline            = tmpspline;
        this->c_spineSegL       = tmpspline[0].spineSegLength;
        this->lastTailFitError = dfitPtError_total/c_spinePoints;


#ifdef _ZTFDEBUG_
        qDebug() << "Converged in n: " << cntpass;
#endif


///  DEBUG ///
#ifdef _ZTFDEBUG_
    for (int j=0; j<c_spinePoints;j++) //Rectangle Eye
    {
        cv::circle(frameDebugC,cv::Point(spline[j].x,spline[j].y),2,TRACKER_COLOURMAP[j],1);
    }
    cv::drawContours(frameDebugC,contours_body,idxOuterContour,CV_RGB(200,20,20),1);

#endif
///    End Debug ///

    return dfitPtError_total; //Return Total Fit Error
   // qDebug() << "D err:" << dDifffitPtError_total;
}

