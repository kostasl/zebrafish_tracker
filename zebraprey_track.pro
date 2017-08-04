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
    fastms/src/examples/example_batchprocessing.cpp \
    fastms/src/examples/example_gui.cpp \
    fastms/src/examples/main.cpp \
    fastms/src/libfastms/solver/solver.cpp \
    fastms/src/libfastms/solver/solver_base.cpp \
    fastms/src/libfastms/solver/solver_host.cpp \
    fastms/src/libfastms/util/has_cuda.cpp \
    fastms/src/libfastms/util/image_mat.cpp \
    fastms/src/mex/fastms_mex.cpp \
    GUI/TrackerScene.cpp \
    SinhaSIFT/MySIFT.cpp \
    SinhaSIFT/SIFT.cpp \
    SinhaSIFT/stdafx.cpp
    #GUI/TrackerScene.cpp

RESOURCES += qml.qrc
INCLUDEPATH += /usr/include/opencv
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
QMAKE_LFLAGS += -Wl,--rpath=/home/kostasl/opencv/build/lib/
QMAKE_LFLAGS_RPATH=
###Using command : cp `ldd larvatrack | sed -re s/^.+\=\>// | sed -re 's/^(.+) \(.+\)/\1/'` /libs

# Default rules for deployment.
include(deployment.pri)


LIBS+=-L/home/kostasl/Qt/5.8/gcc_64/lib/ #Compilation At office DEsktop

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
    fastms/src/examples/example_batchprocessing.h \
    fastms/src/examples/example_gui.h \
    fastms/src/examples/param.h \
    fastms/src/examples/util.h \
    fastms/src/libfastms/solver/solver.h \
    fastms/src/libfastms/solver/solver_base.h \
    fastms/src/libfastms/solver/solver_common_operators.h \
    fastms/src/libfastms/solver/solver_device.h \
    fastms/src/libfastms/solver/solver_host.h \
    fastms/src/libfastms/util/has_cuda.h \
    fastms/src/libfastms/util/image.h \
    fastms/src/libfastms/util/image_access.h \
    fastms/src/libfastms/util/image_access_convert.h \
    fastms/src/libfastms/util/image_mat.h \
    fastms/src/libfastms/util/mem.h \
    fastms/src/libfastms/util/real.h \
    fastms/src/libfastms/util/sum.h \
    fastms/src/libfastms/util/timer.h \
    fastms/src/libfastms/util/types_equal.h \
    fastms/src/mex/mex_util.h \
    GUI/TrackerScene.hpp \
    SinhaSIFT/Descriptor.h \
    SinhaSIFT/KeyPoint.h \
    SinhaSIFT/SIFT.h \
    SinhaSIFT/stdafx.h \
    SinhaSIFT/targetver.h


FORMS += \
    GUI/mainwindow.ui
