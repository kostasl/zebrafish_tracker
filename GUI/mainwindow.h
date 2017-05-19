#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <larvatrack.h>
#include <QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QPixmap>


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

    ~MainWindow();

private:
    Ui::MainWindow      *ui;

    /**
     * @brief Scene Object for displaying image on Form
     */
    QGraphicsScene      *mScene;

    /**
    * Image to be displayed in  this scene.
    */
    QGraphicsPixmapItem*                            mImage;


};

#endif // MAINWINDOW_H
