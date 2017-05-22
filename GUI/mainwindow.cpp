#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "QtOpencvCore.hpp"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    this->mScene = new QGraphicsScene(this->ui->graphicsView);

    this->ui->graphicsView->setScene(this->mScene);
    this->ui->graphicsView->setFixedSize(1280,1024);
    mScene->setSceneRect(this->ui->graphicsView->rect());
    //Add Empty/New PixMap on Which we will set the images onto
    this->mImage = mScene->addPixmap(QPixmap());


    this->mScene->installEventFilter(this);

}

void MainWindow::showCVimg(cv::Mat& img)
{
    QImage qimg = QtOpencvCore::img2qimg(img);

    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimg);
    // scale pixMap image to fit into the QLabel
    //        pixMap = pixMap.scaled(this->ui->graphicsView->size(), Qt::KeepAspectRatio);


    //this->mImage->setPixmap(pixMap);
    this->mScene->setSceneRect(0, 0, qimg.width(), qimg.height());
    this->ui->graphicsView->setFixedSize(qimg.width(),qimg.height());

    mImage->setPixmap(pixMap);
    this->mImage->setPos(-0, -0);

    //this->ui->graphicsView->fitInView(mImage);
    this->ui->graphicsView->show();
    //mImage

}



bool MainWindow::eventFilter(QObject *obj, QEvent *event) {
    if (event->type() == QEvent::GraphicsSceneWheel) {
        handleWheelOnGraphicsScene(dynamic_cast<QGraphicsSceneWheelEvent*> (event));

        // Don't propagate
        event->accept();
        return true;
    }

    if (event->type() == QEvent::KeyPress) {
         QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
         qDebug() << "Ate key press" << keyEvent->key();
         return true;
     } else {
         return false;
     }

    return false;
}



void MainWindow::handleWheelOnGraphicsScene(QGraphicsSceneWheelEvent* scrollevent)
{
  const int degrees = scrollevent->delta()  / 8;
  qDebug() << degrees;

  int steps = degrees / 15;
  qDebug() << steps;

  double scaleFactor = 1.0; //How fast we zoom
  const qreal minFactor = -100.0;
  const qreal maxFactor = 100.0;
  qreal h11 = 1.0, h22 = 0;

  if(steps > 0)
  {
     h11 = (h11 >= maxFactor) ? h11 : (h11 + scaleFactor);
     h22 = (h22 >= maxFactor) ? h22 : (h22 + scaleFactor);
  }
  else
 {
     h11 = (h11 <= minFactor) ? minFactor : (h11 - scaleFactor);
     h22 = (h22 <= minFactor) ? minFactor : (h22 - scaleFactor);
 }
    this->ui->graphicsView->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    this->ui->graphicsView->setTransform(QTransform(h11, 0, 0,0, h22, 0, 0,0,1));

}


MainWindow::~MainWindow()
{
    delete ui;
}
