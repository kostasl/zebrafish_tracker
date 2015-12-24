// Copyright (C) 2007 by Cristóbal Carnero Liñán
// grendel.ccl@gmail.com
//
// This file is part of cvBlob.
//
// cvBlob is free software: you can redistribute it and/or modify
// it under the terms of the Lesser GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// cvBlob is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Lesser GNU General Public License for more details.
//
// You should have received a copy of the Lesser GNU General Public License
// along with cvBlob.  If not, see <http://www.gnu.org/licenses/>.
//


#include <cmath>
#include <iostream>
#include <sstream>
#include <stack>
#include <list>
using namespace std;

#if (defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__) || defined(__WINDOWS__) || (defined(__APPLE__) & defined(__MACH__)))
    #include <cv.h>
#else
    #include <opencv/cv.h>
#endif

//Trying to Find Mat includes
#include <cv.h>   		// open cv general include file
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp" //Draw Polyline
#include "cvblob.h"

namespace cvb
{

  double distantBlobTrack(CvBlob const *b, CvTrack const *t)
  {
    double d1;
    if (b->centroid.x<t->minx)
    {
      if (b->centroid.y<t->miny)
        d1 = MAX(t->minx - b->centroid.x, t->miny - b->centroid.y);
      else if (b->centroid.y>t->maxy)
        d1 = MAX(t->minx - b->centroid.x, b->centroid.y - t->maxy);
      else // if (t->miny < b->centroid.y)&&(b->centroid.y < t->maxy)
        d1 = t->minx - b->centroid.x;
    }
    else if (b->centroid.x>t->maxx)
    {
      if (b->centroid.y<t->miny)
        d1 = MAX(b->centroid.x - t->maxx, t->miny - b->centroid.y);
      else if (b->centroid.y>t->maxy)
        d1 = MAX(b->centroid.x - t->maxx, b->centroid.y - t->maxy);
      else
        d1 = b->centroid.x - t->maxx;
    }
    else // if (t->minx =< b->centroid.x) && (b->centroid.x =< t->maxx)
    {
      if (b->centroid.y<t->miny)
           d1 = t->miny - b->centroid.y;
      else if (b->centroid.y>t->maxy)
            d1 = b->centroid.y - t->maxy;
      else 
	return 0.;
    }

    double d2;
    if (t->centroid.x<b->minx)
    {
      if (t->centroid.y<b->miny)
            d2 = MAX(b->minx - t->centroid.x, b->miny - t->centroid.y);
      else if (t->centroid.y>b->maxy)
            d2 = MAX(b->minx - t->centroid.x, t->centroid.y - b->maxy);
      else // if (b->miny < t->centroid.y)&&(t->centroid.y < b->maxy)
            d2 = b->minx - t->centroid.x;
    }
    else if (t->centroid.x>b->maxx)
    {
      if (t->centroid.y<b->miny)
           d2 = MAX(t->centroid.x - b->maxx, b->miny - t->centroid.y);
      else if (t->centroid.y>b->maxy)
           d2 = MAX(t->centroid.x - b->maxx, t->centroid.y - b->maxy);
      else
           d2 = t->centroid.x - b->maxx;
    }
    else // if (b->minx =< t->centroid.x) && (t->centroid.x =< b->maxx)
    {
      if (t->centroid.y<b->miny)
            d2 = b->miny - t->centroid.y;
      else if (t->centroid.y>b->maxy)
            d2 = t->centroid.y - b->maxy;
      else 
	return 0.;
    }

    return MIN(d1, d2);
  }

  // Access to matrix
#define C(blob, track) close[((blob) + (track)*(nBlobs+2))]
  // Access to accumulators
#define AB(label) C((label), (nTracks))
#define AT(id) C((nBlobs), (id))
  // Access to identifications
#define IB(label) C((label), (nTracks)+1)
#define IT(id) C((nBlobs)+1, (id))
  // Access to registers
#define B(label) blobs.find(IB(label))->second
#define T(id) tracks.find(IT(id))->second

  void getClusterForTrack(unsigned int trackPos, CvID *close, unsigned int nBlobs, unsigned int nTracks, CvBlobs const &blobs, CvTracks const &tracks, list<CvBlob*> &bb, list<CvTrack*> &tt);

  void getClusterForBlob(unsigned int blobPos, CvID *close, unsigned int nBlobs, unsigned int nTracks, CvBlobs const &blobs, CvTracks const &tracks, list<CvBlob*> &bb, list<CvTrack*> &tt)
  {
    for (unsigned int j=0; j<nTracks; j++)
    {
      if (C(blobPos, j))
      {
	tt.push_back(T(j));

	unsigned int c = AT(j);

	C(blobPos, j) = 0;
	AB(blobPos)--;
	AT(j)--;

	if (c>1)
	{
	  getClusterForTrack(j, close, nBlobs, nTracks, blobs, tracks, bb, tt);
	}
      }
    }
  }

  void getClusterForTrack(unsigned int trackPos, CvID *close, unsigned int nBlobs, unsigned int nTracks, CvBlobs const &blobs, CvTracks const &tracks, list<CvBlob*> &bb, list<CvTrack*> &tt)
  {
    for (unsigned int i=0; i<nBlobs; i++)
    {
      if (C(i, trackPos))
      {
        bb.push_back(B(i));

        unsigned int c = AB(i);

        C(i, trackPos) = 0;
        AB(i)--;
        AT(trackPos)--;

        if (c>1)
        {
          getClusterForBlob(i, close, nBlobs, nTracks, blobs, tracks, bb, tt);
        }
      }
    }
  }

  void cvUpdateTracks(CvBlobs const &blobs, CvTracks &tracks, ltROIlist& vRoi, const double thDistance, const unsigned int thInactive, const unsigned int thActive)
  {
    CV_FUNCNAME("cvUpdateTracks");
    __CV_BEGIN__;

    unsigned int nBlobs = blobs.size();
    unsigned int nTracks = tracks.size();

    // Proximity matrix:
    // Last row/column is for ID/label.
    // Last-1 "/" is for accumulation.
    CvID *close = new unsigned int[(nBlobs+2)*(nTracks+2)]; // XXX Must be same type than CvLabel.

    //KL: Note Huge Try Block
    try
    {
      // Inicialization:
      unsigned int i=0;
      for (CvBlobs::const_iterator it = blobs.begin(); it!=blobs.end(); ++it, i++)
      {
        AB(i) = 0;
        IB(i) = it->second->label;
      }

      //KL:Reassign Max Track ID - Search through all trackss
      CvID maxTrackID = 0;
      unsigned int j=0;
      for (CvTracks::const_iterator jt = tracks.begin(); jt!=tracks.end(); ++jt, j++)
      {
        AT(j) = 0;
        IT(j) = jt->second->id;
        if (jt->second->id > maxTrackID)
          maxTrackID = jt->second->id;
      }

      // Proximity matrix calculation and "used blob" list inicialization:
      for (i=0; i<nBlobs; i++)
        for (j=0; j<nTracks; j++)
          if (C(i, j) = (distantBlobTrack(B(i), T(j)) < thDistance))
          {
            AB(i)++;
            AT(j)++;
          }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Detect inactive tracks
      for (j=0; j<nTracks; j++)
      {
        unsigned int c = AT(j);

        if (c==0)
        {
          //cout << "Inactive track: " << j << endl;

          // Inactive track.
          CvTrack *track = T(j);
          track->inactive++;
          track->label = 0;
        }
       }

      // Detect new tracks
      for (i=0; i<nBlobs; i++)
      {
        unsigned int c = AB(i);

        if (c==0)
        {
          CvBlob *blob = B(i);
          //cout << "Blob (new track): " << maxTrackID+1 << endl;
          //cout << *B(i) << endl;
          //Check If Blob is within some ROI
          cv::Point pntCentroid = cv::Point(blob->centroid.x,blob->centroid.y);
          //KL: Detect Which ROI
          ltROI* proi = ltGetFirstROIContainingPoint(vRoi ,pntCentroid);
          if (proi == 0)
              continue; //Ignore this blob its out of ROI

          // New track.
          maxTrackID++;
          CvTrack *track = new CvTrack;

          track->id = maxTrackID;
          track->label = blob->label;
          track->minx = blob->minx;
          track->miny = blob->miny;
          track->maxx = blob->maxx;
          track->maxy = blob->maxy;
          track->centroid = blob->centroid;
          track->lifetime = 0;
          track->active = 0;
          track->inactive = 0;


          track->pROI = proi; //Set Pointer to ROI containing the 1st blob
          track->pointStack.push_back(pntCentroid); //Add 1st Point to list of Track
          tracks.insert(CvIDTrack(maxTrackID, track));
        }
      }

      // Clustering
      for (j=0; j<nTracks; j++)
      {
        unsigned int c = AT(j);

        if (c)
        {
          list<CvTrack*> tt; tt.push_back(T(j));
          list<CvBlob*> bb;

          getClusterForTrack(j, close, nBlobs, nTracks, blobs, tracks, bb, tt);

          // Select track
          //KL :SEG FAULT is caused by these searches failing -low rate occurance)
          CvTrack *track = NULL;
          unsigned int area = 0;
          for (list<CvTrack*>::const_iterator it=tt.begin(); it!=tt.end(); ++it)
          {
            CvTrack *t = *it;

            unsigned int a = (t->maxx-t->minx)*(t->maxy-t->miny);
            if (a>area)
            {
              area = a;
              track = t;
            }
          }

          // Select blob //KL SET TO NULL Detect not found
          CvBlob *blob = NULL;
          area = 0;
          //cout << "Matching blobs: ";
          for (list<CvBlob*>::const_iterator it=bb.begin(); it!=bb.end(); ++it)
          {
            CvBlob *b = *it;

            //cout << b->label << " ";
            //Remove blobs that are not in the same ROI as the track - and those that fail the filter
            if (track != NULL)
            {
               ltROI* blbroi = ltGetFirstROIContainingPoint(vRoi ,cv::Point(b->centroid.x,b->centroid.y) );
                if (blbroi == 0 )
                    continue;

                if (b->area>area && *blbroi == *track->pROI )
                {
                  area = b->area;
                  blob = b;
                }
            }
          }
          //cout << endl;
          //KL: SKip If not found - Think this Matching is what updates the state of the tracks
          if  (track != NULL && blob != NULL )
          {

                  // Update track
                  //cout << "Matching: track=" << track->id << ", blob=" << blob->label << endl;
                  track->label = blob->label;
                  track->centroid = blob->centroid;
                  //KL: Make A point list
                  track->pointStack.push_back(cv::Point(blob->centroid.x,blob->centroid.y)); //KL:Add The new point to the List
                  track->minx = blob->minx;
                  track->miny = blob->miny;
                  track->maxx = blob->maxx;
                  track->maxy = blob->maxy;

                  if (track->inactive)
                    track->active = 0;
                  track->inactive = 0;
           }


          // Others to inactive
          for (list<CvTrack*>::const_iterator it=tt.begin(); it!=tt.end(); ++it)
          {
            CvTrack *t = *it;

            if (t!=track)
            {
              //cout << "Inactive: track=" << t->id << endl;
              t->inactive++;
              t->label = 0;
            }
          }
        }
      }//CLUSTERINg
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      for (CvTracks::iterator jt=tracks.begin(); jt!=tracks.end();)
        if ((jt->second->inactive>=thInactive)||((jt->second->inactive)&&(thActive)&&(jt->second->active<thActive)))
        {
          delete jt->second;
          tracks.erase(jt++);
        }
        else
        {
          jt->second->lifetime++;
          if (!jt->second->inactive)
            jt->second->active++;
          ++jt;
        }
     } //Closes Huge Try block
        catch (...)
        {
          delete[] close;
          throw; // TODO: OpenCV style.
        }

        delete[] close;

    __CV_END__;
  }

  CvFont *defaultFont = NULL;

  void cvRenderTracks(CvTracks const tracks, IplImage *imgSource, IplImage *imgDest, unsigned short mode, CvFont *font )
  {
    CV_FUNCNAME("cvRenderTracks");
    __CV_BEGIN__;

    CV_ASSERT(imgDest&&(imgDest->depth==IPL_DEPTH_8U)&&(imgDest->nChannels==3));

    if ((mode&CV_TRACK_RENDER_ID)&&(!font))
    {
      if (!defaultFont)
      {
	font = defaultFont = new CvFont;
	cvInitFont(font, CV_FONT_HERSHEY_DUPLEX, 0.5, 0.5, 0, 1);
	// Other fonts:
	//   CV_FONT_HERSHEY_SIMPLEX, CV_FONT_HERSHEY_PLAIN,
	//   CV_FONT_HERSHEY_DUPLEX, CV_FONT_HERSHEY_COMPLEX,
	//   CV_FONT_HERSHEY_TRIPLEX, CV_FONT_HERSHEY_COMPLEX_SMALL,
	//   CV_FONT_HERSHEY_SCRIPT_SIMPLEX, CV_FONT_HERSHEY_SCRIPT_COMPLEX
      }
      else
	font = defaultFont;
    }

    if (mode)
    {
        for (CvTracks::const_iterator it=tracks.begin(); it!=tracks.end(); ++it)
        {
            cvRenderTrack(*((*it).second) ,it->first , imgSource, imgDest, mode, font );
        }
    }

    __CV_END__;
  }


  void cvRenderTrack(CvTrack& track,const unsigned int trackID, IplImage *imgSource, IplImage *imgDest, unsigned short mode, CvFont *font )
  {
      CV_FUNCNAME("cvRenderTrack");
      __CV_BEGIN__;

        if (mode&CV_TRACK_RENDER_ID)
          if (!track.inactive)
          {
            stringstream buffer;
            buffer << trackID;
            cvPutText(imgDest, buffer.str().c_str(), cvPoint((int)track.centroid.x, (int)track.centroid.y), font, CV_RGB(0.,255.,0.));
          }

        if (mode&CV_TRACK_RENDER_BOUNDING_BOX)
          if (track.inactive)
            cvRectangle(imgDest, cvPoint(track.minx, track.miny), cvPoint(track.maxx-1,track.maxy-1), CV_RGB(0., 0., 50.));
          else
            cvRectangle(imgDest, cvPoint(track.minx, track.miny), cvPoint(track.maxx-1, track.maxy-1), CV_RGB(0., 0., 255.));

        if (mode&CV_TRACK_RENDER_TO_LOG)
        {
          clog << "Track " << track.id << endl;
          if (track.inactive)
            clog << " - Inactive for " << track.inactive << " frames" << endl;
          else
            clog << " - Associated with blob " << track.label << endl;
          clog << " - Lifetime " <<track.lifetime << endl;
          clog << " - Active " << track.active << endl;
          clog << " - Bounding box: (" << track.minx << ", " <<track.miny << ") - (" << track.maxx << ", " << track.maxy << ")" << endl;
          clog << " - Centroid: (" << track.centroid.x << ", " << track.centroid.y << ")" << endl;
          clog << endl;
        }

        if (mode&CV_TRACK_RENDER_TO_STD)
        {
          cout << "Track " << track.id << endl;
          if (track.inactive)
            cout << " - Inactive for " <<track.inactive << " frames" << endl;
          else
            cout << " - Associated with blobs " << track.label << endl;
          cout << " - Lifetime " << track.lifetime << endl;
          cout << " - Active " << track.active << endl;
          cout << " - Bounding box: (" <<track.minx << ", " << track.miny << ") - (" << track.maxx << ", " <<track.maxy << ")" << endl;
          cout << " - Centroid: (" << track.centroid.x << ", " << track.centroid.y << ")" << endl;
          cout << endl;
        }


        //Render Path
        //cv::Mat img = cv::Mat::zeros(400, 400, CV_8UC3);
        //if (mode&CV_TRACK_RENDER_PATH)
        //{
            std::vector<CvPoint>* pvec = &track.pointStack;
            CvPoint *pts = (CvPoint*) cv::Mat(track.pointStack).data;
            int npts = cv::Mat(track.pointStack).rows;
            //Random colour
            int c1 =  rand() % 200 + 30;
            int c2 =  rand() % 200 + 30;
            int c3 =  rand() % 200 + 30;
            cvPolyLine(imgDest, &pts,&npts, 1,
                            false, 			// draw open contour (i.e. joint end to start)
                            cv::Scalar(c1,c2,c3),// colour RGB ordering (here = green)
                            1, 		        // line thickness
                            CV_AA, 0);
        //}
__CV_END__;
  }


} //END OF NAMESPACE

