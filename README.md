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


### How do I get set up? 

* Dependencies
    - OpenCV
        -Install opencv and compile with flag -D WITH_QT=ON, improves experience and range of available functions.
        https://opencv.org/releases/. Latest known opencv version that worked with tracker is OpenCV 3.4.4
    - GSL dev 2.4
        -Random number generator library
    - QT4 or above
    - tensorflow 
    - GCC (buildtools)
    

```
sudo apt-get install git libopencv-dev qt5-default g++ libgsl-dev
```

### Tensorflow
Then follow instructions for [ install the tensorflow C library ](https://www.tensorflow.org/install/lang_c)
 
Download & extract :

    - Linux CPU only 	https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-2.7.0.tar.gz
    - Linux GPU support 	https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-gpu-linux-x86_64-2.7.0.tar.gz

On Linux and macOS, you may want to extract to /usr/local/lib:
```
FILENAME=libtensorflow-cpu-linux-x86_64-2.7.0.tar.gz
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/${FILENAME}
sudo tar -C /usr/local -xzf ${FILENAME}
```

then 

On Linux/macOS, if you extract the TensorFlow C library to a system directory, such as /usr/local, configure the linker with ldconfig:
```
sudo ldconfig /usr/local/lib
```

* Building the tracker
```
git clone https://github.com/dafishcode/zebrafishtrack
cd zebrafishtrack
qmake zebraprey_track.pro
make
```

### Using the app ###

The application offers a GUI from which to select input video and where to export the CSV data files.
These can also be set by command line options allowing to start the tracker through scripts and process multiple video files automatically.

User:
 * Chooses input video file, then on the second dialogue choose the text file to export track info in CSV format.
The directory where the data csv files are exported must have  format EXP_N_YYYYMMDD_Ssec giving the experiment number N (can be a range or anything really), the date when the embryos were collected, and the timelapse frame period in seconds.

Command line options:

    - help h usage ? print a help  message
    - outputdir   o Dir where To save sequence of images 
    - invideofile v  Behavioural Video file to analyse 
    - invideolist f A text file listing full path to video files to process
    - startframe s 1  Video Will start by Skipping to this frame
    - stopframe p 0  Video Will stop at this frame
    - startpaused P 0 Start tracking Paused On 1st Frame/Need to Run Manually
    - duration d | 0  | Number of frames to Track for starting from start frame
    - logtofile l |  | Filename to save clog stream to 
    - ModelBG b | 0  | Learn and Substract Stationary Objects from Foreground mask
    - BGThreshold bgthres | 30  | Absolute grey value used to segment BG (g_Segthresh)
    - SkipTracked t | 0  | Skip Previously Tracked Videos
    - PolygonROI r | 0  | Use pointArray for Custom ROI Region
    - ModelBGOnAllVids a | 1  | Only Update BGModel At start of vid when needed
    - FilterPixelNoise pn | 0  | Filter Pixel Noise During Tracking (Note:This has major perf impact so use only when necessary due to pixel noise. BGProcessing does it by default)
    - DisableOpenCL ocl | 0  | Disabling the use of OPENCL can avoid some SEG faults hit when running multiple trackers in parallel
    - EnableCUDA cuda | 0  | Use CUDA for MOG, and mask processing - if available  
    - HideDataSource srcShow | 0  | Do not reveal datafile source, so user can label data blindly  
    - EyeHistEqualization histEq | 0  | Use hist. equalization to enhance eye detection contrast  
    - TrackFish ft | 1  | Track Fish not just the moving prey 
    - MeasureMode M | 0 | Click 2 points to measure distance to prey


### Control keys for tracker during Run time :

Some important shortcut keys: 

 *  Press p to pause Image. once paused:
 *  s to save snapshots in CSV outdir pics subfolder.
 *  2 Left Clicks to define the 2 points of region-of interest for tracking.
 *  m to show the masked image of the larva against BG.
 *  t Start Tracking
 *  q Exit Quit application.
 *  T save detected region as template image.
 
## Output

The tracker produces N files Vn_XXX.csv - one for each ROI n defined by the order by which the ROI was created. First ROI 1, 2nd ROI gets V2_... For each ROI n it generates the following three types of file:

* Vn_pos.csv -> contains summary counts of the number of blobs, and the number of tracks identified on each video frame.
* Vn_pos.csv -> identified Blob  Positions and area, per frame - file header is : frameN, SerialN, BlobLabel, Centroid_X, Centroid_Y, Area 
* Vn_pos_tracks.csv -> is the main file for using to plot, analyse tracks. It columns info are    frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive. This file allows analysis of tracks by giving positions Centroid_X, Centroid_Y along with TrackID - ie X,Y positions with the same TrackID belong to the same track (as identified by the tracker).
* if the Save-Mode is enabled (via pressing the s key) then a subfolder named "pics" will contain an image per frame of the tracked video. This can be then combined into a new "tracked" video with avcodec or ffmpeg.

*Note:* The package of source files contains example MATLAB scripts that can process these output files, plot and extract statistics.

## Fish detection and tracking 

Two methods are used, the first step uses a Deep Neural Network and the second (optional) step uses template matching to accurate detect orientation of the larva's head.
The DNN is trained using a python script called _trainTF_DNN.py_ over a set of prelabelled images found in the tensorDNN/img subirectory.
The labels are taken from the directory names where each set of images is stored. For now I have included 2 sets for fish labels, one is small scale as in early behavioural videos, and the other is large fish - as in my latest recordings (lower fps).
Images for training are extracted from videos by double click on fish to create new Box, and pressing CTRL+T to save as pgm file in the templates subdirectory. Make sure that fish training samples are accuratelly positioned fish anterior box (manually oriented using "[" and "]" keys)
You may also want to save false positives (error templates) to train the classifier on non-fish objects too.
Convert pgm to jpg prior to copying new samples in the training folder by running :

<code>
mogrify -format jpg *.pgm && rm *.pgm
</code>

For the template matching, the tracker uses small set  (<10) images included in the executables resources as default, but also loads the templates found in the templates subdirectory of the tracking output folder. 
This allows to the user to save a few templates that are specific to the video being tracked.
All templates need to be of the same size - (rows,columns) and the the size of the 1st loaded template defines the 

Tracking of fish position is assisted by a Kalman Filter - these preditions and measurements appear as blue and yello box marking position of fish's anterior.

## Prey detection and tracking

I used an G-H filter to track the motion of the prey items.
This simple position filtering technique (G-H Filter) improves tracking between subsequent frames 
 by removing any sudden non-prey like jumps and noise from the estimated prey position.

The filter params where adjusted against prey motion data, and these need to be scaled based on 
video scaling spatial (field of view/resolution)  and time (fps). I have added crude scaling to the video fps and assume the video frame contains a 35 mm circular arena.

The tracker also uses Kalman filtering for larva motion, rotation and eye angles.

### Video Processing 

The program can accept a list of videos to process as a parameter to invideolist pointing to a txt file with a file name on each line.
If you run into problems with broken videos, such as the number of frames not being correctly recognized, try fixing it with :
```
ffmpeg -v error -i source_broken_video.mp4 -c copy fixed_video.mp4
```


### Contribution guidelines ###

* Code review : This is a prototype of a tracker and a new version of the tracker is required with a improved code organization.
* Writing tests
* Speed optimizations for classification and detection of angle
* Improved GUI and user guide.

### Who do I talk to? ###

* Repo owner or admin (kostasl)
* Other community or team contact
