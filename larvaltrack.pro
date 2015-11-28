TEMPLATE = app

QT += qml quick widgets gui


SOURCES += main.cpp

RESOURCES += qml.qrc
INCLUDEPATH += /usr/include/opencv
INCLUDEPATH += /usr/include/cvblog

LIBS += `pkg-config opencv cvblob --libs`
LIBS += -L/usr/local/lib -lcvblob
# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Default rules for deployment.
include(deployment.pri)
