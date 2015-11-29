TEMPLATE = app

QT += qml quick widgets gui


SOURCES += main.cpp

RESOURCES += qml.qrc
INCLUDEPATH += /usr/include/opencv
INCLUDEPATH += /usr/include/cvblob

#INCLUDEPATH += /home/kostasl/workspace/cvblobLib
#`pkg-config opencv cvblob --cflags`

LIBS += `pkg-config opencv --libs`
#LIBS += -L/usr/local/lib -lcvblob
#LIBS += -L/home/kostasl/workspace -lcvblob
# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Default rules for deployment.
include(deployment.pri)


unix:!macx: LIBS += -L$$PWD/../build-cvblobLib-Desktop-Debug/ -lcvblobLib

INCLUDEPATH += $$PWD/../cvblob
DEPENDPATH += $$PWD/../cvblob
