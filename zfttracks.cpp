
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




void zftRenderTrack(zftTrack& track, const cv::Mat& frameIn, cv::Mat& frameOut, unsigned short mode, int fontface,float fontScale )
{

     if (mode&CV_TRACK_RENDER_ID)
     {
//        if (track.inactive < 200)
//        {
          std::stringstream buffer;
          buffer << track.id;
          cv::putText(frameOut, buffer.str().c_str(), (cv::Point)(track.centroid),fontface,fontScale, CV_RGB(0.,255.,0.));
     }

      if (mode&CV_TRACK_RENDER_BOUNDING_BOX)
        if (track.inactive)
          cv::rectangle(frameOut,track.boundingBox,CV_RGB(0., 0., 50.));
        else
          cv::rectangle(frameOut,track.boundingBox,CV_RGB(0., 0., 255.));

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
                          CV_AA, 0);
          //delete pts;
          if (mTrack.u)
            qDebug() << "mTrack.u->refcount" << mTrack.u->refcount;

          mTrack.release();
      }
}






