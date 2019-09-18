# README Zebrafish Behaviour Analysis Software #

@author : Konstantinos Lagogiannis
Developed while postdoc in Martin Meyer's lab, KCL

### What is this repository for? ###

* Quick summary
A software that tracks moving larvae in timelapse videos, exporting 3 data csv files for each region of interest - vial. _N file contains the larvae count per frame, _pos contains position information of identified blobs (most possibly larvae) in the ROI, _pos_track.csv file exports the data from a basic tracking algorithm that attempts to link blobs across frames via track ID. 
The package further contains a subfolder /matlabscripts, that allows for importing and analysing basic statistics of the exported data. 
The data is not assumed clean, and some problems from the tracker are dealt with within the those scripts to clear the data.
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)


### Using the app ###
User:
 * Chooses input video file, then on the second dialogue choose the text file to export track info in CSV format.
The directory where the data csv files are exported must have  format EXP_N_YYYYMMDD_Ssec giving the experiment number N (can be a range or anything really),the date when the embryos were collected, and the timelapse frame period in seconds).

 * The video begins paused -  use to left mouse clicks to define a new region in the image over which you want to count the larvae.
 * Define the regions of interest using two left mouse click for each region. The red box defines the region over which the larvae are counted-tracked and recorded.

### Control keys for tracker during Run time :
 * Press p to pause Image. once paused:
 *  s to save snapshots in CSV outdir pics subfolder.
 *  2 Left Clicks to define the 2 points of region-of interest for tracking.
 *  m to show the masked image of the larva against BG.
 *  t Start Tracking
 *  q Exit Quit application

## Output
The tracker produces N files Vn_XXX.csv - one for each ROI n defined by the order by which the ROI was created. First ROI 1, 2nd ROI gets V2_... For each ROI n it generates the following three types of file:

* Vn_pos.csv -> contains summary counts of the number of blobs, and the number of tracks identified on each video frame.
* Vn_pos.csv -> identified Blob  Positions and area, per frame - file header is : frameN,SerialN,BlobLabel,Centroid_X,Centroid_Y,Area 
* Vn_pos_tracks.csv -> is the main file for using to plot, analyse tracks. It columns info are    frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive. This file allows analysis of tracks by giving positions Centroid_X, Centroid_Y along with TrackID - ie X,Y positions with the same TrackID belong to the same track (as identified by the tracker).
* if the Save-Mode is enabled (via pressing the s key) then a subfolder named "pics" will contain an image per frame of the tracked video. This can be then combined into a new "tracked" video with avcodec or ffmpeg.

*Note:* The package of source files contains example MATLAB scripts that can process these output files, plot and extract statistics.



 
### How do I get set up? ###

* Summary of set up
Install opencv 3.0 / best compiled with with flag -D WITH_QT=ON, improves experience and range of available functions.
http://opencv.org/downloads.html

* Uses a version of cvBlob which is embedded in the package, for independent development.
* Configuration
* Dependencies
 openCV3.0 (compilied with WITH_QT support)
 Qt4 and above

        sudo apt-get install git libopencv-dev qt5-default g++

* Building
        git clone https://github.com/i-git/FIMTrack.git FIMTrack
        cd larvatrack
        qmake larvatrack.pro
        make


* How to run tests
get videos with vials,
save videos in correct folder naming.make output folders, run application from build folder. Choose video(s) (more than one can be selected if they are in sequence), choose output file location, set regions of interest, click r (run) and then t (track).

Wait until the video(s) have been tracked in full, the tracks will be displayed on screen realtime as the tracker is processing the video.

* Deployment instructions


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin (kostasl)
* Other community or team contact
