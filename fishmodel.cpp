#include "fishmodel.h"
#include "ellipse_detect.h"
#include "config.h"

extern cv::Mat frameDebugC;
extern cv::Size gszTemplateImg;


//extern int gFishTailSpineSegmentLength;
//extern int gFitTailIntensityScanAngleDeg;
//extern const int gcFishContourSize; //Fixed number of fish Contour Points
const int gMaxFitIterations = ZTF_TAILFITMAXITERATIONS;
extern double gTemplateMatchThreshold;

fishModel::fishModel()
{
        inactiveFrames = 0;
        templateScore = 0;
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());
        coreTriangle.push_back(cv::Point());

        this->mouthPoint.x = 0;
        this->mouthPoint.y = 0;
        //this->leftEyeHull.clear();
        //this->rightEyeHull.clear();
        this->ID    = 0;
        zTrack.id   = this->ID;
        zTrack.colour = CV_RGB(255,0,0);

        c_spineSegL           = gFishTailSpineSegmentLength;
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

    inactiveFrames  = 0;
    this->ID        = blob.hash() ;
    this->blobLabel = blob.hash();

    zTrack.id       = this->ID;


    this->zfishBlob = blob; //Copy Localy
    //this->track     = NULL;
    this->bearingRads = bestTemplateOrientation*CV_PI/180.0;
    this->bearingAngle = bestTemplateOrientation;
    this->ptRotCentre    = ptTemplateCenter;
    zTrack.centroid = ptTemplateCenter;
   // this->coreTriangle[2].x = this->zfishBlob.pt.x;
    //this->coreTriangle[2].y = this->zfishBlob.pt.y;


    templateScore           = 0;
    this->resetSpine();

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

return leftEye.rectEllipse.angle;
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

return leftEye.rectEllipse.angle;

}

///
/// \brief fishModel::resetSpine make a straight Spline pointing towards Blobs Bearings Angle
///
void fishModel::resetSpine()
{
    //
    this->spline.clear();
    spline.reserve(c_spinePoints);

    for (int i=0;i<c_spinePoints;i++)
    {

        splineKnotf sp;
        //1st Spine Is in Opposite Direction of Movement and We align 0 degrees to be upwards (vertical axis)
        //if (this->bearingRads > CV_PI)
        if (this->bearingRads < 0)
            this->bearingRads += 2.0*CV_PI;

            sp.angleRad    = (this->bearingRads)-CV_PI ; //  //Spine Looks In Opposite Direction

            if (sp.angleRad < 0)
                sp.angleRad += 2.0*CV_PI;
        //else
//            sp.angleRad    = (this->bearingRads)+CV_PI/2.0; //CV_PI/2 //Spine Looks In Opposite Direcyion

        assert(!std::isnan(sp.angleRad && std::abs(sp.angleRad) <= 2*CV_PI && (sp.angleRad) >= 0 ));

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

    for (int i=1;i<c_spinePoints;i++)
    {

       outspline[i].x = outspline[i-1].x + ((double)c_spineSegL)*sin(outspline[i-1].angleRad);
       outspline[i].y = outspline[i-1].y - ((double)c_spineSegL)*cos(outspline[i-1].angleRad);
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
    const double dAngleStep = sgn*CV_PI/36.0;
    double ret = 0.0;
    outspline = inspline;

    //If idxparam = 0 or 1 then we are varying initial Spline Point
    if (idxparam == 0) //x0, y0 params
    {
        ret = sgn*0.05;
        outspline[0].x -= ret;

    }else if (idxparam == 1)
    {
        ret = sgn*0.05;
        outspline[0].y -= ret;

    }else //Index > 1 is spine Angles
    {
        outspline[idxparam-2].angleRad += dAngleStep;// Angle variation for this theta
        ret = dAngleStep*this->c_spineSegL+cos(dAngleStep)*this->c_spineSegL; //rTheta
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
    for (int i=0;i<(c_spineParamCnt);i++)
    {
        if (i==0)
            inspline[0].x = inparams[0];
        if (i==1)
            inspline[0].y = inparams[1];
        if (i>1)
            inspline[i-2].angleRad = inparams[i]; //Param 3 is actually 1st spine knot's angle
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
    const float spineLength = this->c_spineSegL*pspline.size(); //The fitted Spine's lentgh is fixed
    int idx = z/this->c_spineSegL; //Find knot index which is contains point
    double segLen = z - idx*this->c_spineSegL;  //Modulo Find length input var along a linear segment

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

    const int spineLength = this->c_spineSegL*pspline.size(); //The fitted Spine's lentgh is fixed
    const double dCStep = 0.3; //Step size on when searching along Spine Curve for closest Point (Foot Point )tk

    int idxNear = 0;
    cv::Point2f ptFoot = cv::Point2f(pspline[idxNear].x,pspline[idxNear].y);
    //double distX = pow(ptsrc.x-ptFoot.x,2);
    //double distY = pow(ptsrc.y-ptFoot.y,2);
    double mindist;// = sqrt(distX + distY);
    mindist= cv::norm(ptsrc-ptFoot);
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


///
/// \brief fishModel::Update - Called On Every FrameProcessed To Update Model State
/// The track point is set to the blob position and not the template centre
/// \param fblob
/// \param templatematchScore
/// \param Angle
/// \param bcentre
///

void fishModel::updateState(zftblob* fblob,double templatematchScore,int Angle, cv::Point2f bcentre,unsigned int nFrame)
{
    nCurrentFrame = nFrame; //Set Last Update To Current Frame
    this->templateScore  = templatematchScore;
    this->bearingAngle   = Angle;
    this->bearingRads   =  Angle*CV_PI/180.0;
    this->ptRotCentre    = bcentre;
    this->zfishBlob      = *fblob;
    this->zTrack.pointStack.push_back(fblob->pt);
    this->zTrack.effectiveDisplacement = cv::norm(fblob->pt-this->zTrack.centroid);
    this->zTrack.centroid = bcentre;//fblob->pt; //Or Maybe bcentre
    ///Optimization only Render Point If Displaced Enough from Last One
    if (this->zTrack.effectiveDisplacement > gDisplacementThreshold)
    {
        this->zTrack.pointStackRender.push_back((cv::Point)fblob->pt);
        this->zTrack.active++;
    }else {
        this->zTrack.inactive++;
    }

    ///Update Template Box Bound
    //int bestAngleinDeg = fish->bearingAngle;
    cv::RotatedRect fishRotAnteriorBox(bcentre,gszTemplateImg ,Angle);
    /// Save Anterior Bound
    this->bodyRotBound = fishRotAnteriorBox;



    //Set Spine Source to Rotation Centre
    this->spline[0].x       = bcentre.x;
    this->spline[0].y       = bcentre.y;
    //this->spline[0].angleRad   = this->bearingRads+CV_PI; //+180 Degrees so it looks in Opposite Direction


    assert(!std::isnan(this->bearingRads));

}


///
/// \brief fishModel::fitSpineToContour / Least squares fit of constrained spline model with params x0,yo, and theta_i-n
///
/// \param contours_body
/// \param idxInnerContour
/// \param idxOuterContour
/// \return fitness error score
///
double fishModel::fitSpineToContour(cv::Mat& frameImg_grey, std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour)
{
    static const int cntParam = this->c_spineParamCnt;
    static const int gcFishContourSize = ZTF_FISHCONTOURSIZE+1;//Fixed PLus 1 tail Point

    const int c_fitErrorPerContourPoint = 10; //Parameter Found By Experience for current size fish
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
    int cntStuck    = 0;
    int cntSolved   = 0;
    double dVarScale    = 1.0;
    //Do A number of Passes Before  Convergence //&& (dfitPtError_total/contour.size() > 8)
    while (cntpass < gMaxFitIterations && (cntStuck < 5) && (cntSolved < 3) )
    {
        //Converged But Error Is still Large per Countour Point Then Jolt
        if (std::abs(dDifffitPtError_total) < 0.01 && dfitPtError_total/contour.size() > 10) //Time Out Convergece Count
        {
            cntStuck++;
            dVarScale = dVarScale*1.12;
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
        dTemp = (double)cntpass/gMaxFitIterations;
        //Prob Of Acceptance P = exp(-(sn-s)/T) >= drand
        ///
        dfitPtError_total_last  = dfitPtError_total;
        dfitPtError_total       = 0.0; //Reset

        double dq,ds; //Variation In Space And Score Variation
        //For Each Contour Point
        for (uint i=0;i<gcFishContourSize;i+=1) //For Each Data point make a row in Jacobian
        {
            dResiduals[i] = distancePointToSpline((cv::Point2f)contour[i],tmpspline);

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

    //If Convergece TimedOut Then likely the fit is stuck with High Residual and no gradient
    //Best To reset Spine and Start Over Next Time
    /// \todo Make this parameter threshold formal
    if (dfitPtError_total/contour.size() > c_fitErrorPerContourPoint)
    {
        this->resetSpine(); //No Solution Found So Reset
         qDebug() << "Reset Spine after n:" << cntpass;
    }
    else //Update Spine Model
    {
        this->spline = tmpspline;
#ifdef _ZTFDEBUG_
        qDebug() << "Converged in n: " << cntpass;
#endif
    }



///  DEBUG ///
#ifdef _ZTFDEBUG_
    for (int j=0; j<c_spinePoints;j++) //Rectangle Eye
    {
        cv::circle(frameDebugC,cv::Point(spline[j].x,spline[j].y),2,TRACKER_COLOURMAP[j],1);
    }
    cv::drawContours(frameDebugC,contours_body,idxOuterContour,CV_RGB(200,20,20),1);

#endif
///    End Debug ///

    return dDifffitPtError_total;
   // qDebug() << "D err:" << dDifffitPtError_total;
}


void fishModel::drawBodyTemplateBounds(cv::Mat& outframe)
{

    int bestAngleinDeg = this->bearingAngle;
    //cv::RotatedRect fishRotAnteriorBox(centre, cv::Size(gLastfishimg_template.cols,gLastfishimg_template.rows),bestAngleinDeg);


//    stringstream strLbl;
//    strLbl << "A: " << bestAngleinDeg;
    QString strlbl("A: " + QString::number(bestAngleinDeg));

    cv::Scalar colour;
    if (this->templateScore >= gTemplateMatchThreshold)
        colour = CV_RGB(250,250,0);
    else
        colour = CV_RGB(30,30,250);

    cv::putText(outframe,strlbl.toStdString(),this->bodyRotBound.boundingRect().br(),CV_FONT_NORMAL,0.4,colour,1);


    ///Draw a Red Rotated Frame around Detected Body
    cv::Point2f boundBoxPnts[4];
    this->bodyRotBound.points(boundBoxPnts);
    for (int j=0; j<4;j++) //Rectangle Body
        cv::line(outframe,boundBoxPnts[j],boundBoxPnts[(j+1)%4] ,colour,1);

}

///
/// \brief fishModel::fitSpineToIntensity implements  Giovanni's tail fitting method : starting from initial point on the fish body and an initial direction for the tail it searches for the
///   highest intensity within a angle range -c_tailscanAngle +c_tailscanAngle degrees on a largely blurred scene image , through each spine point incrementally.
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


    uint loc,pxValMax;

    int angle; //In Deg of Where The spline point is looking towards - Used by Ellipse Arc Drawing

    //cv::Mat frameimg_Blur;


    //imgframeIn.convertTo(frameimg_Blur,CV_32F,1./255);
    //draw_inv=ones-draw_inv;

    //cv::GaussianBlur(imgframeIn,frameimg_Blur,cv::Size(5,5),5,5);
    //cv::imshow("IntensitTailFit",frameimg_Blur);

    std::vector<cv::Point> ellipse_pts; //Holds the Drawn Arc Points around the last spine Point
    //cv::imshow("BlurSpine",frameimg_Blur);

    for(unsigned int k=1;k<AP_N;k++){ //Loop Through Spine Points
        ellipse_pts.clear();


        //Get Angl
        //angle = (atan2(tgt.y,tgt.x) +_CV_PI)/CV_PI*180.0; //Find towards Next Point Angle In Degrees

        angle = spline[k-1].angleRad/CV_PI*180.0-90.0; //Get Angle In Degrees for Arc Drawing Tranlated Back to 0 horizontal
        //angle = spline[k].angleRad/CV_PI*180.0; //Get Angle In Degrees for Arc Drawing Around THe point this Spine Was looking At initially
        //Construct Elliptical Circle around last Spine Point - of Radius step_size
        cv::ellipse2Poly(cv::Point(spline[k-1].x,spline[k-1].y), cv::Size(step_size,step_size), 0, angle-c_tailscanAngle, angle+c_tailscanAngle, 1, ellipse_pts);

        ///sSearch for Maximum Intensity on Source Image Along the Arc
        pxValMax=0;
        for(int idx=0;idx<(int)ellipse_pts.size();idx++){
            //Obtain Value From Image at Point on Arc
            int x = std::min(frameimg_Blur.cols,std::max(1,ellipse_pts[idx].x));
            int y = std::min(frameimg_Blur.rows,std::max(1,ellipse_pts[idx].y));
            loc=frameimg_Blur.at<uchar>(cv::Point(x,y));

            //If New Maximum Found THen Update Spline Point to point to this and Update Previous Spline Point's angle
            if(loc>=pxValMax){ //Gt Or Equal - Othewise Points Get Stuck between frames int the black Bg when loc = 0
                spline[k].x     = ellipse_pts[idx].x;
                spline[k].y     = ellipse_pts[idx].y;
                ///Get Arc tan and Translate back to 0 being the Vertical Axis
                spline[k-1].angleRad = std::atan2(spline[k].y-spline[k-1].y,spline[k].x-spline[k-1].x)+CV_PI/2; // ReCalc Angle in 0 - 2PI range Of previous Spline POint to this New One
                //spline[k].angleRad = spline[k-1].angleRad;
                //Constrain Large Deviations
                if (std::abs(spline[k-1].angleRad - spline[k].angleRad) > CV_PI/3.0)
                    spline[k].angleRad = spline[k-1].angleRad; //Spine Curvature by Initializing next spine point Constraint Next


                pxValMax=loc; //Save as New Maximum Point
            }
        }

        //out<<tgt<<' '<<angle<<' '<<tmp_pts[k]<<' '<<loc<<' '<<index<<' '<<ellipse_pts.size()<<endl;

        //anchor_pts[k]=tmp_pts[k];
    }

}




void fishModel::drawSpine(cv::Mat& outFrame)
{
    for (int j=0; j<c_spinePoints-1;j++) //Rectangle Eye
    {
        cv::circle(outFrame,cv::Point(spline[j].x,spline[j].y),2,TRACKER_COLOURMAP[j],2);
        //if (j<(c_spinePoints-1))
            cv::line(outFrame,cv::Point(spline[j].x,spline[j].y),cv::Point(spline[j+1].x,spline[j+1].y),TRACKER_COLOURMAP[0],1);
//        else
//        { //Draw Terminal (hidden) point - which is not a spine knot
//            cv::Point ptTerm;
//            ptTerm.x = spline[j].x + ((double)c_spineSegL)*sin(spline[j].angleRad);
//            ptTerm.y = spline[j].y - ((double)c_spineSegL)*cos(spline[j].angleRad);

//            cv::line(outFrame,cv::Point(spline[j].x,spline[j].y),ptTerm,TRACKER_COLOURMAP[0],1);
//        }
    }
    cv::circle(outFrame,cv::Point(spline[c_spinePoints-1].x,spline[c_spinePoints-1].y),2,TRACKER_COLOURMAP[c_spinePoints-1],1);

}



///
/// \brief operator << //Overloaded Stream Operator // Output Current State Of The Track
/// \param out
/// \param h
/// \return
///
std::ostream& operator<<(std::ostream& out, const zftTrack& h)
{

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)
    if (h.pointStack.size() > 0)
    {
        zftTrackPoint ptt = h.pointStack.back();
        out << ptt.x << "\t" << ptt.y;
    }
    else
        out << 0 << "\t" << 0;

    return out;
}

///
/// \brief operator << //Overloaded Stream Operator // Output Current State Of The Track
/// \param out
/// \param h
/// \return
///
QTextStream& operator<<(QTextStream& out, const zftTrack& h)
{

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)
if (h.pointStack.size() > 0)
{
    zftTrackPoint ptt = h.pointStack.back();
    out << ptt.x << "\t" << ptt.y;
}else
    out << 0 << "\t" << 0;

    return out;
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


    out << h.nCurrentFrame <<"\t"<< h.ID <<"\t"<< h.bearingAngle <<"\t" << h.zTrack;

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
    const double Rad2Deg = (180.0/CV_PI);

    //for (auto it = h.pointStack.begin(); it != h.pointStack.end(); ++it)
    out.setRealNumberNotation(QTextStream::RealNumberNotation::FixedNotation );
    out.setRealNumberPrecision(2);
    out << h.ID <<"\t"<< h.bearingAngle <<"\t" << h.zTrack << "\t" << h.leftEyeTheta << "\t" <<  h.rightEyeTheta;

    //Set Global 1st Spine Direction (Helps to detect Errors)
    assert(h.spline.size() > 0);

    out << "\t" << Rad2Deg* h.spline[0].angleRad;
    //Output Spine Point Angular Deviations from the previous spine/tail Segment in Degrees
    for (int i=1;i<h.c_spinePoints;i++)
    {

       out << "\t" << Rad2Deg*( h.spline[i-1].angleRad - h.spline[i].angleRad);

    }


    return out;
}








