#ifndef ZFTTRACKS
#define ZFTTRACKS


/// \brief Type of identification numbers.
typedef unsigned int CvID;


///// \brief Struct that contain information about one track.
///// \see CvID
///// \see CvLabel
///// \see ltROI
//struct CvTrack
//{
//  CvID id; ///< Track identification number.
//  ltROI* pROI; ///< Pointer To Region of Interest structure to which this track belongs
//  CvLabel label; ///< Label assigned to the blob related to this track.
//  CvScalar colour = CV_RGB(255., 0., 0.); ///> Colourwhen drawing Countour

//  unsigned int minx; ///< X min.same as  corresponding blob bounding box
//  unsigned int maxx; ///< X max.
//  unsigned int miny; ///< Y min.
//  unsigned int maxy; ///< y max.

//  CvPoint2D64f centroid; ///< Centroid.
//  CvTrackPoints pointStack; /// <Holds list of past centroid positions along the track
//  double effectiveDisplacement; ///< Used to indicate a px speed measure so as to estimate possible blob distance from track on next frame.
//  unsigned int lifetime; ///< Indicates how much frames the object has been in scene.
//  unsigned int active; ///< Indicates number of frames that has been active from last inactive period.
//  unsigned int inactive; ///< Indicates number of frames that has been missing.
//};

///// \var typedef std::map<CvID, CvTrack *> CvTracks
///// \brief List of tracks.
///// \see CvID
///// \see CvTrack
//typedef std::map<CvID, CvTrack *> CvTracks;

///// \var typedef std::pair<CvID, CvTrack *> CvIDTrack
///// \brief Pair (identification number, track).
///// \see CvID
///// \see CvTrack
//typedef std::pair<CvID, CvTrack *> CvIDTrack;



#endif // ZFTTRACKS

