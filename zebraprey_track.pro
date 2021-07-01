#Good to have compiled opencv with qt support so the window can autoresize.
##Suggestion Is though that for Qt5 do not Compile With With_QT=ON, as this will link to Qt4
##Also: May get link error to OpenGL and Qt5Test.so.5 when OpenCV is linked a different version of QT than then one used
# to compile project. SOLUTION: change soft links to Qt5 libs in OS dirs to point to the correct QT lib version ones:
# e.g.:
# sudo rm /usr/lib/x86_64-linux-gnu/libQt5Test.so.5
# sudo ln -s /opt/Qt/5.15.0/gcc_64/lib/libQt5Test.so.5 /usr/lib/x86_64-linux-gnu/libQt5Test.so.5


TEMPLATE = app

QT += widgets gui qml quick testlib
#QTDIR = "/opt/Qt/5.15.0"


CONFIG += c++11
CONFIG += warn_off

SOURCES += main.cpp \
    #cvBlob/cvcontour.cpp \
    #cvBlob/cvcolor.cpp \
    #cvBlob/cvblob.cpp \
    #cvBlob/cvaux.cpp \
    GUI/QtOpencvCore.cpp \
    GUI/mainwindow.cpp \
    fishmodel.cpp \
    CSS/CurveCSS.cpp \
    GUI/TrackerScene.cpp \
    #GUI/TrackerScene.cpp
    ellipse_detect.cpp \
    preymodel.cpp \
    template_detect.cpp \
    zfttracks.cpp \
    fgmaskprocessing.cpp \
    errorhandlers.cpp \
    config.cpp \
    eyesdetector.cpp \
    fishdetector.cpp

RESOURCES += qml.qrc

RC_ICONS = myappico.ico

QT_CONFIG -= no-pkg-config
CONFIG += link_pkgconfig
PKGCONFIG += opencv gsl   #or whatever package here

##pkg-config --libs $(pkg-config --print-requires --print-requires-private glfw3)
#pkg-config --list-all

#INCLUDEPATH += `pkg-config opencv --cflags`

#INCLUDEPATH += /home/kostasl/workspace/cvblobLib
##Note: you can pass multiple items to pkg-config as input, so running

##Figure out VERSION : pkg-config --modversion opencv
##Or Check CV_MAJOR_VERSION, CV_MINOR_VERSION

##LIBS+=-L/home/kostasl/Qt/5.8/gcc_64/lib/ #Compilation At office DEsktop
#LIBS += `pkg-config opencv --libs`
#LIBS +=-lgsl -lgslcblas -lm
#LIBS +=-lm /home/kostasl/Qt/5.15.1/gcc_64/lib/libQt5OpenGL.so.5 /home/kostasl/Qt/5.15.1/gcc_64/lib/libQt5OpenGL.so.5 /home/kostasl/Qt/5.15.1/gcc_64/lib/libQt5Test.so.5
#LIBS +=-lm /opt/Qt/5.15.0/gcc_64/lib/libQt5OpenGL.so.5
#LIBS += -lm /opt/Qt/5.15.0/gcc_64/lib/libQt5Test.so.5

#QMAKE_LIBDIR = /opt/Qt/5.15.0/gcc_64/lib/
QMAKE_LIBDIR_OPENGL = $LD_LIBRARY_PATH #///opt/Qt/5.15.0/gcc_64/lib/
#unix|win32: LIBS += -lQt5OpenGL

QMAKE_CFLAGS_DEBUG += -v -da -Q
QMAKE_CFLAGS += -rdynamic
# Additional import path used to resolve QML modules in Qt Creator's code model
##QML_IMPORT_PATH =
##
##Assume Libs are copied with the package into
QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN/libs
#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/Qt/5.11.1/gcc_64/lib/ ##Office
#QMAKE_LFLAGS += -Wl,--rpath=/media/extStore/opt/Qt3.0.1/5.9.2/gcc_64/lib/ #Home
#QMAKE_LFLAGS += -Wl,--rpath=/opt/Qt/5.9/5.9/gcc_64/lib/ #Home

#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/opencv/build/lib/
#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/OpenCV/opencv-3.3.0/build-Dbg
#QMAKE_LFLAGS += -Wl,--rpath=/home/klagogia1/OpenCV/opencv-3.3.1/build #Home
#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/OpenCV/opencv-3.3.1/build #Office

QMAKE_LFLAGS_RPATH=
###Using command : cp `ldd larvatrack | sed -re s/^.+\=\>// | sed -re 's/^(.+) \(.+\)/\1/'` /libs

# Default rules for deployment.
include(deployment.pri)


#OTHER_FILES += \
#    ../cvblob/cvBlob/cvBlobLib.pro.user \
#    cvBlob/cvBlobLib.pro.user \
#    cvBlob/cvBlobLib.pro~ \
#    cvBlob/CMakeLists.txt

SUBDIRS += \
    cvBlob/cvBlobLib.pro

HEADERS += \
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
    zfttracks.h \
    config.h \
    foodmodel.h \
    fgmaskprocessing.h \
    errorhandlers.h \
    eyesdetector.h \
    fishdetector.h


FORMS += \
    GUI/mainwindow.ui

DISTFILES += \
    README.md \
    img/fishbodyb_tmp.pgm \
    zebraprey_track.supp \
    img/fishbody_tmp9.pgm


