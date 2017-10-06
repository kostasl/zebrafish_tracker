
///\file Contains the code to process Blobs and produce Tracks
/// - for fish and food separatelly
/// \author Kostas Lagogiannis
///




//bool operator<(const fishmodel& a,const tDetectedEllipsoid& b) {
//  return a.fitscore < b.fitscore; //Max Heap
//}





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

