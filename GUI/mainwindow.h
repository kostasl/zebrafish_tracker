#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <larvatrack.h>
#include <QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QGraphicsSceneWheelEvent>
#include <QPixmap>
#include <QEvent>
#include <QKeyEvent>
#include <QStringListModel>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/video/video.hpp>
#include "opencv2/video/background_segm.hpp"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void showCVimg(cv::Mat& img);
    void showInsetimg(cv::Mat& img); //Used to Show Small IMage Next to main scene
    void showInsetTemplateimg(cv::Mat& img); //For Showing The template Image
    void showVideoFrame(cv::Mat& img,unsigned int nFrame);

    void saveScreenShot();
    void saveTemplateImage(cv::Mat& templ);
    void tickProgress();
    void setTotalFrames(uint FrameCount);
    void UpdateTailSegSizeSpinBox(float fTailSize);
    void LogEvent(QString strMessage);
    unsigned int nFrame = 0;
    unsigned int nTotalFrameCount = 0;


    QString stroutDirCSV;
    QString vidFilename;


    ~MainWindow();

public slots:
    void echoChanged(int);
    void changeEvent(QEvent *e);
    void textEdited(const QString);
    void fishvalueChanged(int i);
    void eyevalueChanged(int i);
    void tailSizevalueChanged(int i);
    void maxEllipseSizevalueChanged(int i);
    void minEllipseSizevalueChanged(int i);


protected:
    bool eventFilter(QObject *obj, QEvent *event);
    void handleWheelOnGraphicsScene(QGraphicsSceneWheelEvent* scrollevent);
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseDblClickEvent(QGraphicsSceneMouseEvent * mouseEvent );
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void dragMoveEvent(QGraphicsSceneDragDropEvent* mouseEvent );
    void handleSliderChange(QEvent* event);
    void createSpinBoxes();

    QStringListModel*                     mModelMessageList;
    QStringList                           mMessageList;
    QImage qimg; //SCene Image Updated in ShowCV Image
    QImage qimgHead; //SCene Image Updated in ShowCV Image
    cv::Mat frameScene; //CvMat Last Frame Drawn
private slots:
    void on_spinBoxTemplateThres_valueChanged(int arg1);

    void on_spinBoxMOGBGRatio_valueChanged(int arg1);

    void on_actionTrack_Fish_triggered(bool checked);

    void on_actionTrack_Food_triggered(bool checked);

    void on_actionRecord_Tracks_to_File_w_triggered(bool checked);

    void on_actionQuit_triggered();

    void on_actionPaus_tracking_p_triggered(bool checked);

    void on_actionStart_tracking_triggered();

    void on_actionPaus_tracking_p_triggered();


    void on_checkBoxGPU_toggled(bool checked);

    void on_checkBoxMOG_toggled(bool checked);

    void on_checkBoxNoiseFilter_toggled(bool checked);

    void on_spinBoxFoodThresMin_editingFinished();

    void on_spinBoxSpineSegSize_valueChanged(int arg1);

    void on_spinBoxFoodThresMin_valueChanged(int arg1);

    void on_spinBoxFoodThresMax_valueChanged(int arg1);

private:
    Ui::MainWindow      *ui;

    /**
     * @brief Scene Object for displaying image on Form
     */
    QGraphicsScene      *mScene;
    QGraphicsScene*     mInsetScene;
    QGraphicsScene*     mInsetTemplateScene;

    /**
    * Image to be displayed in  this scene.
    */
    QGraphicsPixmapItem*                            mImage;
    QGraphicsPixmapItem*                            mImageInset;
    QGraphicsPixmapItem*                            mImageTemplateInset;

    cv::Mat*                                        mpLastCVImg;

    ///



};

#endif // MAINWINDOW_H


