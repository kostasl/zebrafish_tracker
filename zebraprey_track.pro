#Good to have compiled opencv with qt support so the window can autoresize.

TEMPLATE = app

QT += widgets gui #qml quick


SOURCES += main.cpp \
    cvBlob/cvtrack.cpp \
    cvBlob/cvlabel.cpp \
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

RESOURCES += qml.qrc
#INCLUDEPATH += /usr/include/opencv
INCLUDEPATH += /home/kostasl/OpenCV/opencv-3.1.0/include
#INCLUDEPATH += /usr/include/cvblob
#INCLUDEPATH += ~/opencv/

#INCLUDEPATH += /home/kostasl/workspace/cvblobLib
##Note: you can pass multiple items to pkg-config as input, so running
##pkg-config --libs $(pkg-config --print-requires --print-requires-private glfw3)

#`pkg-config opencv cvblob --cflags`
##Figure out VERSION : pkg-config --modversion opencv
##Or Check CV_MAJOR_VERSION, CV_MINOR_VERSION


LIBS += `pkg-config opencv --libs`
#LIBS += -L/usr/local/lib -lcvblob
#LIBS += -L/home/kostasl/workspace -lcvblob
# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =
##
##Assume Libs are copied with the package into
QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN/libs
QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/Qt/5.8/gcc_64/lib/
#QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/opencv/build/lib/
QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/OpenCV/opencv-3.1.0/build
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
    GUI/TrackerScene.hpp


FORMS += \
    GUI/mainwindow.ui
