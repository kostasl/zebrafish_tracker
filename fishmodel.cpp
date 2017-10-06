#include "fishmodel.h"
#include "ellipse_detect.h"
fishModel::fishModel()
{
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
}

fishModel::fishModel(cvb::CvTrack* track,cvb::CvBlob* blob):fishModel()
{


    this->ID        = track->id;
    this->blobLabel = track->label;
    this->track     = track; //Copy Pointer
    this->bearingRads = cvb::cvAngle(blob);
    this->coreTriangle[2].x = track->centroid.x;
    this->coreTriangle[2].y = track->centroid.y;

    templateScore           = 0;
    this->resetSpine();
}

fishModel::fishModel(zftblob blob):fishModel()
{

    this->ID        = blob.hash() ;
    this->blobLabel = blob.hash();
    zTrack.id       = this->ID;

    this->zfishBlob = blob; //Copy Localy
    this->track     = NULL;
    this->bearingRads = blob.angle*M_PI/180.0;
    this->coreTriangle[2].x = this->zfishBlob.pt.x;
    this->coreTriangle[2].y = this->zfishBlob.pt.y;


    templateScore           = 0;
    this->resetSpine();

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
        sp.angle    = this->bearingRads;
        assert(!std::isnan(sp.angle));
        if (i==0)
        {
            sp.x =  this->coreTriangle[2].x;
            sp.y =  this->coreTriangle[2].y;
        }
        else
        {
            sp.x        = spline[i-1].x - ((double)c_spineSegL)*cos(sp.angle);
            sp.y        = spline[i-1].y - ((double)c_spineSegL)*sin(sp.angle);
        }

        spline.push_back(sp); //Add Knot to spline
    }

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

       outspline[i].x = outspline[i-1].x - ((double)c_spineSegL)*cos(outspline[i-1].angle);
       outspline[i].y = outspline[i-1].y - ((double)c_spineSegL)*sin(outspline[i-1].angle);
       assert(!std::isnan(outspline[i].y) && !std::isnan(outspline[i].x));
    }

}

///
/// \brief fishModel::getSpine VARIATION Make a Spine variation modifying a specific param value indicated by index
/// \param inspline vector Passed by value so the original is unchanged
/// \param outspline - The variational Spline
/// \return distance of variation in Config space
double fishModel::getdeltaSpline(t_fishspline inspline, t_fishspline& outspline,int idxparam,double sgn)
{
    const static int cntParam = outspline.size()+1;
    const double dAngleStep = sgn*CV_PI/16.0;
    double ret = 0.0;
    outspline = inspline;

    //If idxparam = 0 or 1 then we are varying initial Spline Point
    if (idxparam == 0) //x0, y0 params
    {
        ret = sgn*0.05;
        outspline[0].x += ret;

    }else if (idxparam == 1)
    {
        ret = sgn*0.05;
        outspline[0].y += ret;

    }else
    { // Param INdex is > 1 so it refers to angles starting from 0 idx knot
         //ret = 0.0;//*this->c_spineSegL+cos(dAngleStep)*this->c_spineSegL; //rTheta
         ret = dAngleStep*this->c_spineSegL+cos(dAngleStep)*this->c_spineSegL; //rTheta
         outspline[idxparam-2].angle += dAngleStep;// Angle variation for this theta
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
    outparams[0] = inspline[0].x;
    outparams[1] = inspline[0].y;

    for (int i=0;i<(c_spinePoints-1);i++)
    {
        outparams[i+2] = inspline[i].angle;
    }

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
            inspline[i-2].angle = inparams[i]; //Param 3 is actually 1st spine knot's angle
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
    ptC.x = pspline[idx].x - segLen*cos(pspline[idx].angle);
    ptC.y = pspline[idx].y - segLen*sin(pspline[idx].angle);

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

    return mindist;
}

///
/// \brief fishModel::fitSpineToContour / Least squares fit of constrained spline model with params x0,yo, and theta_i-n
///
/// \param contours_body
/// \param idxInnerContour
/// \param idxOuterContour
/// \return fitness error score
///
double fishModel::fitSpineToContour(std::vector<std::vector<cv::Point> >& contours_body,int idxInnerContour,int idxOuterContour)
{
    const int cntParam = this->c_spineParamCnt;
    const int gMaxFitIterations = 5;

    ///Param sfish model should contain initial spline curve (Hold Last Frame Position)

    //Run Until Convergence Error is below threshold - Or Change is too small

    ///Compute Error terms for all data points/obtain local quadratic approx of fsd
    //For each contour Point
    std::vector<cv::Point> contour = contours_body[idxOuterContour];
    t_fishspline tmpspline = this->spline;
    t_fishspline dsSpline; //Variational Spline

    //Measure squared Distance error to closest Curve(spline) Point
        //Add to total error
    double dfitPtError_total = 0.0;
    double dfitPtError_total_last = 0.0;
    double dDifffitPtError_total = 1000.0;

    double dTemp = 1.0; //Anealling Temperature

    double dJacobian[contour.size()][cntParam];//Vector of \nabla d for error functions
    memset(dJacobian,0.0,contour.size()*(cntParam)*sizeof(double));
    double dGradf[cntParam];//Vector of Grad F per param
    memset(dGradf,0.0,cntParam*sizeof(double));
    double dResiduals[contour.size()];//Vector of \nabla d for error functions
    memset(dResiduals,0.0,contour.size()*sizeof(double));

    int cntpass = 0;
    while (cntpass < gMaxFitIterations && abs(dDifffitPtError_total) > 0.2)
    {
        cntpass++;
        ///For Annealing
        dTemp = (double)cntpass/gMaxFitIterations;
        //Prob Of Acceptance P = exp(-(sn-s)/T) >= drand
        ///
        dfitPtError_total_last  = dfitPtError_total;
        dfitPtError_total = 0.0; //Reset

        for (int i=0;i<contour.size();i+=1) //For Each Data point make a row in Jacobian
        {
            dResiduals[i] = distancePointToSpline((cv::Point2f)contour[i],tmpspline);

            dfitPtError_total       +=dResiduals[i];

            for (int k=0;k<cntParam; k++) //Add Variation dx to each param and calc derivative
            {   /// \note using only +ve dx variations and not -dx - In this C space Ds magnitude should be symmetrical to dq anyway
                double dq = getdeltaSpline(tmpspline,dsSpline,k,+0.75); //Return param variation
                double ds = distancePointToSpline((cv::Point2f)contour[i],dsSpline); // dsSpline residual of variation spline
                //dsSpline.clear();
                //getdeltaSpline(tmpspline,dsSpline,k,-0.25) ; //add dx
                //ds += distancePointToSpline((cv::Point2f)contour[i],dsSpline); // Add df dsSpline residual of variation spline

                //if (dq > 0.0)
                dJacobian[i][k] = (ds-dResiduals[i])/(dq);
                dGradf[k]           += dResiduals[i]*dJacobian[i][k]; //Error Grad - gives Gradient in Cspace vars to Total error

            }

        } //Loop Through Data Points

        std::vector<double> cparams;
        getSplineParams(tmpspline,cparams);

        ///modify CSpace Params with gradient descent
        for (int i=0;i<cntParam;i++)
        {
            cparams[i] -= 0.005*dGradf[i];
            //qDebug() << "lamda GradF_"<< i << "-:" << 0.003*dGradf[i];
        }
        ///Modify Spline - ie move closer
        setSplineParams(tmpspline,cparams);


        dDifffitPtError_total = dfitPtError_total - dfitPtError_total_last; //Total Residual /Error Measure Change

        //if (dfitPtError_total > 1000)
            //this->resetSpine(); //Start over


    }

    qDebug() << "ID:" <<  this->ID << cntpass << " EChange:" << dDifffitPtError_total;

    this->spline = tmpspline;

    return dDifffitPtError_total;
//    ///DEBUG
//    for (int j=0; j<8;j++) //Rectangle Eye
//    {
//        cv::circle(frameDebugC,spline[j],2,TRACKER_COLOURMAP[j],1);
//    }
//    cv::imshow("Debug C",frameDebugC);


   // qDebug() << "D err:" << dDifffitPtError_total;
    //cv::waitKey(100);
} //While Error Change Is larger Than




