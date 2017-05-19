#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "QtOpencvCore.hpp"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->mScene = new QGraphicsScene(this->ui->graphicsView);
    //mScene->setSceneRect(this->ui->graphicsView->rect());
    this->ui->graphicsView->setScene(this->mScene);
    this->ui->graphicsView->setFixedSize(1024,1200);

    //Add Empty/New PixMap on Which we will set the images onto
    this->mImage = mScene->addPixmap(QPixmap());
    this->mImage->setPos(10, 10);


}

void MainWindow::showCVimg(cv::Mat& img)
{
    QImage qimg = QtOpencvCore::img2qimg(img);

    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimg);
    // scale pixMap image to fit into the QLabel
    //        pixMap = pixMap.scaled(this->ui->graphicsView->size(), Qt::KeepAspectRatio);


    //this->mImage->setPixmap(pixMap);
    //this->mScene->setSceneRect(0, 0, img.width(), img.height());
    //this->mScene->clear();
    mImage->setPixmap(pixMap);
    //this->mScene->addPixmap(pixMap);
    this->ui->graphicsView->show();
    //mImage

}

MainWindow::~MainWindow()
{
    delete ui;
}
