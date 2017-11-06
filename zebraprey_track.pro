#Good to have compiled opencv with qt support so the window can autoresize.
##Suggestion Is though that for Qt5 do not Compile With With_QT=ON, as this will link to Qt4

TEMPLATE = app

QT += widgets gui #qml quick

CONFIG += c++11

SOURCES += main.cpp \
    cvBlob/cvcontour.cpp \
    cvBlob/cvcolor.cpp \
    cvBlob/cvblob.cpp \
    cvBlob/cvaux.cpp \
    GUI/QtOpencvCore.cpp \
    GUI/mainwindow.cpp \
    fishmodel.cpp \
    CSS/CurveCSS.cpp \
    GUI/TrackerScene.cpp \
    #GUI/TrackerScene.cpp
    ellipse_detect.cpp \
    template_detect.cpp \
    zfttracks.cpp

RESOURCES += qml.qrc
#INCLUDEPATH += /usr/include/opencv
#INCLUDEPATH += /home/kostasl/OpenCV/opencv-3.2.0/include
#INCLUDEPATH += /home/kostasl/OpenCV/opencv-3.3.0/include
INCLUDEPATH += /home/kostasl/OpenCV/opencv-3.3.1/include
#INCLUDEPATH += /usr/include/cvblob
#INCLUDEPATH += ~/opencv/

#INCLUDEPATH += /home/kostasl/workspace/cvblobLib
##Note: you can pass multiple items to pkg-config as input, so running
##pkg-config --libs $(pkg-config --print-requires --print-requires-private glfw3)

#`pkg-config opencv cvblob --cflags`
##Figure out VERSION : pkg-config --modversion opencv
##Or Check CV_MAJOR_VERSION, CV_MINOR_VERSION

#LIBS += -L /home/kostasl/OpenCV/opencv-3.2.0/build/lib
#LIBS += -L /home/kostasl/OpenCV/opencv-3.3.0/build-Dbg/lib
#LIBS += -L /home/kostasl/OpenCV/opencv-3.3.0/build-Dbg/lib -lopencv_dnn -lopencv_ml -lopencv_objdetect -lopencv_shape -lopencv_stitching -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_imgproc -lopencv_flann -lopencv_viz -lopencv_core
LIBS += -L /home/kostasl/OpenCV/opencv-3.3.1/build/lib
LIBS += `pkg-config opencv --libs`

# Additional import path used to resolve QML modules in Qt Creator's code model
##QML_IMPORT_PATH =
##
##Assume Libs are copied with the package into
QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN/libs
QMAKE_LFLAGS += -Wl,--rpath=/opt/Qt/5.9/5.9/gcc_64/lib/
QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/Qt/5.9.2/gcc_64/lib/

#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/opencv/build/lib/
#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/OpenCV/opencv-3.3.0/build-Dbg
QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/OpenCV/opencv-3.3.1/build

QMAKE_LFLAGS_RPATH=
###Using command : cp `ldd larvatrack | sed -re s/^.+\=\>// | sed -re 's/^(.+) \(.+\)/\1/'` /libs

# Default rules for deployment.
include(deployment.pri)


##LIBS+=-L/home/kostasl/Qt/5.8/gcc_64/lib/ #Compilation At office DEsktop

#unix: LIBS += -L$$PWD/../cvblob/build-cvBlobLib-Desktop-Release/ -lcvBlobLib
#INCLUDEPATH += $$PWD/../cvblob/cvBlob
#DEPENDPATH += $$PWD/../cvblob/cvBlob

OTHER_FILES += \
    ../cvblob/cvBlob/cvBlobLib.pro.user \
    cvBlob/cvBlobLib.pro.user \
    cvBlob/cvBlobLib.pro~ \
    cvBlob/CMakeLists.txt

SUBDIRS += \
    cvBlob/cvBlobLib.pro

HEADERS += \
    cvBlob/cvblob.h \
    larvatrack.h \
    ltROI.h \
    GUI/QtOpencvCore.hpp \
    #GUI/TrackerScene.hpp \
    GUI/mainwindow.h \
    fishmodel.h \
    CSS/CurveCSS.h \
    CSS/std.h \
    GUI/TrackerScene.hpp \
    ellipse_detect.h \
    template_detect.h \
    zfttracks.h


FORMS += \
    GUI/mainwindow.ui

DISTFILES += \
    img/fishbodyb_tmp.pgm \
    zebraprey_track.supp \
    img/fishbody_tmp9.pgm
