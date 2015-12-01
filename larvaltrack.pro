TEMPLATE = app

QT += qml quick widgets gui


SOURCES += main.cpp \
    cvBlob/cvtrack.cpp \
    cvBlob/cvlabel.cpp \
    cvBlob/cvcontour.cpp \
    cvBlob/cvcolor.cpp \
    cvBlob/cvblob.cpp \
    cvBlob/cvaux.cpp

RESOURCES += qml.qrc
INCLUDEPATH += /usr/include/opencv
#INCLUDEPATH += /usr/include/cvblob

#INCLUDEPATH += /home/kostasl/workspace/cvblobLib
#`pkg-config opencv cvblob --cflags`

LIBS += `pkg-config opencv --libs`
#LIBS += -L/usr/local/lib -lcvblob
#LIBS += -L/home/kostasl/workspace -lcvblob
# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Default rules for deployment.
include(deployment.pri)

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
    cvBlob/cvblob.h
