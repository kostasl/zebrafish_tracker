TEMPLATE = app

QT += qml quick widgets gui


SOURCES += main.cpp

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

unix: LIBS += -L$$PWD/../cvblob/build-cvBlobLib-Desktop-Release/ -lcvBlobLib

INCLUDEPATH += $$PWD/../cvblob/cvBlob
DEPENDPATH += $$PWD/../cvblob/cvBlob
