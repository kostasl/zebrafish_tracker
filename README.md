# README Zebrafish Behaviour Analysis Software #

@author : Konstantinos Lagogiannis
Developed while postdoc in Martin Meyer's lab, KCL
to analyse data for article "Learning steers the ontogeny of hunting behaviour in larval zebrafish"
BioArxiv PrePrint: https://doi.org/10.1101/2019.12.19.883157
### What is this repository for? ###

* Quick summary
This is a larval zebrafish behavioural tracking software that processes monochrome IR DarkField highspeed (100-450 fps) images/videos to extract larva body centroid and orientation, Eye angle, tail spine posture and prey count. It also allows for tracking user indicated prey items.
The tracker exports csv table data files containing pixel coordinates and angles in degrees of tracked items that are located with a ROI.
The package further contains a subfolder /RPlots, that allows for importing and analysing the tracker data.
The tracker deals only with extracting information from image pixels, and so filtering and analysis of tracker data is handled externally by those R scripts.

* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)


### Using the app ###
The application offers a GUI from which to select input video and where to export the CSV data files.
These can also be set by command line options allowing to start the tracker through scripts and process multiple video files automatically.


User:
 * Chooses input video file, then on the second dialogue choose the text file to export track info in CSV format.
The directory where the data csv files are exported must have  format EXP_N_YYYYMMDD_Ssec giving the experiment number N (can be a range or anything really),the date when the embryos were collected, and the timelapse frame period in seconds).

*Command line options:
    "{help h usage ? |    | print this help  message}"
    "{outputdir   o |    | Dir where To save sequence of images }"
    "{invideofile v |    | Behavioural Video file to analyse }"
    "{invideolist f |    | A text file listing full path to video files to process}"
    "{startframe s | 1  | Video Will start by Skipping to this frame}"
    "{stopframe p | 0  | Video Will stop at this frame}"
    "{startpaused P | 0  | Start tracking Paused On 1st Frame/Need to Run Manually}"
    "{duration d | 0  | Number of frames to Track for starting from start frame}"
    "{logtofile l |    | Filename to save clog stream to }"
    "{ModelBG b | 0  | Learn and Substract Stationary Objects from Foreground mask}"
    "{BGThreshold bgthres | 30  | Absolute grey value used to segment BG (g_Segthresh)}"
    "{SkipTracked t | 0  | Skip Previously Tracked Videos}"
    "{PolygonROI r | 0  | Use pointArray for Custom ROI Region}"
    "{ModelBGOnAllVids a | 1  | Only Update BGModel At start of vid when needed}"
    "{FilterPixelNoise pn | 0  | Filter Pixel Noise During Tracking (Note:This has major perf impact so use only when necessary due to pixel noise. BGProcessing does it by default)}"
    "{DisableOpenCL ocl | 0  | Disabling the use of OPENCL can avoid some SEG faults hit when running multiple trackers in parallel}"
    "{EnableCUDA cuda | 0  | Use CUDA for MOG, and mask processing - if available  }"
    "{HideDataSource srcShow | 0  | Do not reveal datafile source, so user can label data blindly  }"
    "{EyeHistEqualization histEq | 0  | Use hist. equalization to enhance eye detection contrast  }"
    "{TrackFish ft | 1  | Track Fish not just the moving prey }"
    "{MeasureMode M | 0 | Click 2 points to measure distance to prey}"


### Control keys for tracker during Run time :
 * Press p to pause Image. once paused:
 *  s to save snapshots in CSV outdir pics subfolder.
 *  2 Left Clicks to define the 2 points of region-of interest for tracking.
 *  m to show the masked image of the larva against BG.
 *  t Start Tracking
 *  q Exit Quit application

 
### How do I get set up? ###
* Dependencies
    - OpenCV
        -Install opencv and compile with flag -D WITH_QT=ON, improves experience and range of available functions.
        http://opencv.org/downloads.html. Latest known opencv version that worked with tracker is OpenCV 3.4.4
    - GSL dev 2.4
        -Random number generator library
    - QT4 or above
    - GCC

```
sudo apt-get install git libopencv-dev qt5-default g++ libgsl-dev
```

* Building the tracker
```
git clone https://github.com/dafishcode/zebrafishtrack
cd zebrafishtrack
qmake zebraprey_track.pro
make
```
## Output
The tracker produces N files Vn_XXX.csv - one for each ROI n defined by the order by which the ROI was created. First ROI 1, 2nd ROI gets V2_... For each ROI n it generates the following three types of file:

* Vn_pos.csv -> contains summary counts of the number of blobs, and the number of tracks identified on each video frame.
* Vn_pos.csv -> identified Blob  Positions and area, per frame - file header is : frameN,SerialN,BlobLabel,Centroid_X,Centroid_Y,Area 
* Vn_pos_tracks.csv -> is the main file for using to plot, analyse tracks. It columns info are    frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive. This file allows analysis of tracks by giving positions Centroid_X, Centroid_Y along with TrackID - ie X,Y positions with the same TrackID belong to the same track (as identified by the tracker).
* if the Save-Mode is enabled (via pressing the s key) then a subfolder named "pics" will contain an image per frame of the tracked video. This can be then combined into a new "tracked" video with avcodec or ffmpeg.

*Note:* The package of source files contains example MATLAB scripts that can process these output files, plot and extract statistics.

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin (kostasl)
* Other community or team contact
