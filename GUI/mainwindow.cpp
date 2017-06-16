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
    //this->ui->graphicsView->setFixedSize(1280,1024);
    //mScene->setSceneRect(this->ui->graphicsView->rect());

    mScene->addText("Hello, fish!")->setPos(100,100);
    //Add Empty/New PixMap on Which we will set the images onto
    this->mImage = mScene->addPixmap(QPixmap());

    this->ui->graphicsView->setSceneRect(this->frameGeometry()); // set the scene's bounding rect to rect of mainwindow

    this->ui->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    this->ui->graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);

    this->ui->horizontalSlider->setRange(0,10000);
    this->mScene->installEventFilter(this);
    this->installEventFilter(this); //To Capture Resize


}

void MainWindow::showVideoFrame(cv::Mat& img,unsigned int nFrame)
{
    this->ui->horizontalSlider->setValue(nFrame);

    showCVimg(img);
}

void MainWindow::showCVimg(cv::Mat& img)
{
    QImage qimg = QtOpencvCore::img2qimg(img);

    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimg);
    // scale pixMap image to fit into the QLabel
    pixMap = pixMap.scaled(this->ui->graphicsView->size(), Qt::KeepAspectRatio);


    //this->mImage->setPixmap(pixMap);

    //this->ui->graphicsView->setSceneRect(this->frameGeometry()); // set the scene's bounding rect to rect of mainwindow
    /// A Scene contains graphic objects, but rendering them requires a View. The Scenes size can be larger than the View's
    this->ui->graphicsView->setSceneRect(this->ui->graphicsView->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->mScene->setSceneRect(this->ui->graphicsView->geometry());
    QRect bound = this->ui->graphicsView->geometry();
    //this->ui->graphicsView->setFixedSize(qimg.width(),qimg.height());

    mImage->setPixmap(pixMap);
    //this->ui->graphicsView->fitInView(this->mImage->boundingRect(), Qt::KeepAspectRatio);

    this->mImage->setPos(bound.topLeft().x() ,bound.topLeft().y());


    this->ui->graphicsView->fitInView(mImage, Qt::KeepAspectRatio);
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
         char key =  keyEvent->text().toStdString().at(0);
         qDebug() << "Ate key press " << keyEvent->text().toStdString().c_str() << " k: " << key;

         ::keyCommandFlag(this,key,nFrame);
         return true;
     }

    if (event->type() == QEvent::Resize) {
         QResizeEvent *resizeEvent = static_cast<QResizeEvent*>(event);
         qDebug(" Resized (New Size) - Width: %d Height: %d",
                resizeEvent->size().width(),
                resizeEvent->size().height());



          ui->graphicsView->setFixedSize(resizeEvent->size()*0.9);

         //QRectF bounds = mScene->itemsBoundingRect();
         //bounds.setWidth(bounds.width()*0.9);         // to tighten-up margins
         //bounds.setHeight(bounds.height()*0.9);       // same as above

          //ui->graphicsView->fitInView(bounds, Qt::KeepAspectRatio);
        // ui->graphicsView->centerOn(0, 0);


        // this->ui->gridLayout->setGeometry(this->frameGeometry());
         //this->ui->graphicsView->setSceneRect(this->frameGeometry()); // set the scene's bounding rect to rect of mainwindow
        return true;
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

  //ui->graphicsView->setFixedSize(ui->graphicsView->size()+steps);

//  if(steps > 0)
//  {
//     h11 = (h11 >= maxFactor) ? h11 : (h11 + scaleFactor);
//     h22 = (h22 >= maxFactor) ? h22 : (h22 + scaleFactor);
//  }
//  else
// {
//     h11 = (h11 <= minFactor) ? minFactor : (h11 - scaleFactor);
//     h22 = (h22 <= minFactor) ? minFactor : (h22 - scaleFactor);
// }
//    this->ui->graphicsView->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
//    this->ui->graphicsView->setTransform(QTransform(h11, 0, 0,0, h22, 0, 0,0,1));
    //this->mImage->transform().scale(50,50);
    this->mImage->transform().scale(this->mImage->scale()+ steps,this->mImage->scale()+ steps);
}


MainWindow::~MainWindow()
{
    delete ui;
}
