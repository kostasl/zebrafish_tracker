
///\file Contains the code to process Blobs and produce Tracks
/// - for fish and food separatelly
/// \author Kostas Lagogiannis
///

#include "zfttracks.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stack>
#include <list>
#include <QDebug>




//// New track.
//maxTrackID++;
//CvTrack *track = new CvTrack;

////Copies Blob data to track
//track->id = maxTrackID;
//track->label = blob->label;
//track->minx = blob->minx;
//track->miny = blob->miny;
//track->maxx = blob->maxx;
//track->maxy = blob->maxy;
//track->centroid = blob->centroid;
//track->effectiveDisplacement = sqrt((double)blob->area); //Set To largest value initially
//track->lifetime = 0;
//track->active = 0;
//track->inactive = 0;
////Set Track Colour
////Random colour
//int c1 =  rand() % 200 + 30;
//int c2 =  rand() % 200 + 30;
//int c3 =  rand() % 200 + 30;
//track->colour = CV_RGB(c1,c2,c3);

//track->pROI = proi; //Set Pointer to ROI containing the 1st blob
//track->pointStack.push_back(pntCentroid); //Add 1st Point to list of Track
//tracks.insert(CvIDTrack(maxTrackID, track));




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
        //If the Point Stack Is empty //Report The stable Point Saved on Centroid/
        out << h.centroid.x << "\t" << h.centroid.y;

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
    //If the Point Stack Is empty //Report The stable Point Saved on Centroid/
    out << h.centroid.x << "\t" << h.centroid.y;


    return out;
}

void zftRenderTrack(zftTrack& track, const cv::Mat& frameIn, cv::Mat& frameOut, unsigned short mode, int fontface,float fontScale )
{
      cv::Scalar colTxt =  cv::Scalar(0,250,20,50); //Colour
      cv::Scalar colCircle = cv::Scalar(255,255,0,30);

     if (mode & CV_TRACK_RENDER_HIGHLIGHT)
     {
       colTxt = CV_RGB(240,210,10);//cv::Scalar(255,255,0,90);
     }

     if (mode&CV_TRACK_RENDER_ID)
     {
//        if (track.inactive < 200)
//        {
          std::stringstream buffer;
          buffer << track.id;
          cv::putText(frameOut, buffer.str().c_str(), (cv::Point)(track.boundingBox.tl())+cv::Point(-5,-2),fontface,fontScale,colTxt);
          //cv::putText(frameOut, buffer.str().c_str(), track.boundingBox.tl(),fontface,fontScale, CV_RGB(10.,220.,0.));

     }
    if (mode & CV_TRACK_RENDER_BOUNDING_CIRCLE)
        cv::circle(frameOut,track.centroid, track.boundingBox.width/2,track.colour,1); //Mark Where Search Is Done


      if (mode & CV_TRACK_RENDER_BOUNDING_BOX)
            cv::rectangle(frameOut,track.boundingBox,track.colour);
        //if (track.inactive > 0)
       //   cv::rectangle(frameOut,track.boundingBox,CV_RGB(0., 40., 150.));
       // else
          //cv::rectangle(frameOut,track.boundingBox,CV_RGB(0., 0., 255.));


      if (mode&CV_TRACK_RENDER_TO_LOG)
      {
        std::clog << "Track " << track.id << std::endl;
        if (track.inactive)
          std::clog << " - Inactive for " << track.inactive << " frames" << std::endl;
        else
        std::clog << " - Lifetime " <<track.lifetime << std::endl;
        std::clog << " - Active " << track.active << std::endl;
        std::clog << " - Bounding box: (" << track.boundingBox.tl() << ", " << track.boundingBox.br()  << ")" << std::endl;
        std::clog << " - Centroid: (" << track.centroid.x << ", " << track.centroid.y << ")" << std::endl;
        std::clog << std::endl;
      }

      if (mode&CV_TRACK_RENDER_TO_STD)
      {
        std::cout << "Track " << track.id << std::endl;
        if (track.inactive)
          std::cout << " - Inactive for " <<track.inactive << " frames" << std::endl;
        else
        std::cout << " - Lifetime " << track.lifetime << std::endl;
        std::cout << " - Active " << track.active << std::endl;
        std::clog << " - Bounding box: (" << track.boundingBox.tl() << ", " << track.boundingBox.br()  << ")" << std::endl;
        std::cout << " - Centroid: (" << track.centroid.x << ", " << track.centroid.y << ")" << std::endl;
        std::cout << std::endl;
      }


      //Render Path
      //cv::Mat img = cv::Mat::zeros(400, 400, CV_8UC3);
      if (mode&CV_TRACK_RENDER_PATH)
      {
          //std::vector<cv::Point> plotPts(track.pointStack.begin(), track.pointStack.end());
          cv::Mat mTrack(track.pointStackRender);
          cv::Point *pts = (cv::Point*) mTrack.data;
          int npts = mTrack.rows;
          cv::polylines(frameOut, &pts,&npts, 1,
                          false, 			// draw open contour (i.e. joint end to start)
                          track.colour ,// colour RGB ordering (here = green)
                          1, 		        // line thickness
                          cv::LINE_AA, 0);
          //delete pts;
          //if (mTrack.u)
          //  qDebug() << "mTrack.u->refcount" << mTrack.u->refcount;

          mTrack.release();
      }
}






