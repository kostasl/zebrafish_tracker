#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "QtOpencvCore.hpp"
#include <QStringListModel>
#include <qlineedit.h>

extern fishModels vfishmodels; //Vector containing live fish models
extern bool bPaused;
extern bool bStoreThisTemplate;
extern bool bDraggingTemplateCentre;
extern bool bStartFrameChanged;
bool bSceneMouseLButtonDown;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    this->mScene = new QGraphicsScene(this->ui->graphicsView);
    this->mInsetScene = new QGraphicsScene(this->ui->graphicsViewHead);

    this->ui->graphicsView->setScene(this->mScene);
    this->ui->graphicsViewHead->setScene(this->mInsetScene);
    //this->ui->graphicsView->setFixedSize(1280,1024);
    //mScene->setSceneRect(this->ui->graphicsView->rect());

    mScene->addText("Zebrafish Tracker Scene")->setPos(100,100);
    //Add Empty/New PixMap on Which we will set the images onto
    this->mImage         = mScene->addPixmap(QPixmap());
    this->mImageInset    = mInsetScene->addPixmap(QPixmap());

    this->ui->graphicsView->setSceneRect(this->ui->graphicsView->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->ui->graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    this->ui->graphicsViewHead->setSceneRect(this->ui->graphicsViewHead->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsViewHead->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->ui->graphicsViewHead->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);


    this->ui->horizontalSlider->setRange(0,50000);
    this->mScene->installEventFilter(this);

    this->installEventFilter(this); //To Capture Resize

    // Log Events on List View //
    mModelMessageList = new QStringListModel(this);
    // Populate our model
    mModelMessageList->setStringList(mMessageList);
    // Glue model and view together
    this->ui->listView->setModel(mModelMessageList);


    //this->ui->txtboxFrameNumber->setInputMask("9999999"); //Just Numbers 0-9
    this->ui->txtboxFrameNumber->setReadOnly(false);
    this->ui->txtboxFrameNumber->setMaxLength(5);
    //this->ui->txtboxFrameNumber->setEchoMode(EchoMode::);
    //this->ui->txtboxFrameNumber->installEventFilter(this); //To Capture Text Change
      this->ui->spinBoxFrame->installEventFilter(this);
}

void MainWindow::showVideoFrame(cv::Mat& img,unsigned int nFrame)
{
    this->ui->horizontalSlider->setValue(nFrame);

    showCVimg(img);
}

void MainWindow::saveScreenShot(QString stroutDirCSV,QString vidFilename)
{
    //std::stringstream frameNumberString; frameNumberString << nFrame;
    QString frameNumberString = QString::number(nFrame);

    ::saveImage(frameNumberString,stroutDirCSV,vidFilename,*this->mpLastCVImg);

}

void MainWindow::setTotalFrames(uint FrameCount)
{
    this->ui->horizontalSlider->setMaximum(FrameCount);
    this->nTotalFrameCount = FrameCount;
}

void MainWindow::tickProgress()
{
    this->ui->horizontalSlider->setValue(this->ui->horizontalSlider->value()+1);
    this->ui->txtboxFrameNumber->setText(QString::number(nFrame));
     this->ui->txtboxFrameNumber->setReadOnly(false);

    ui->spinBoxFrame->setValue(nFrame);
}

void MainWindow::showInsetimg(cv::Mat& img)
{

    qimgHead = QtOpencvCore::img2qimg(img);
    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimgHead);
    QRect bound = this->ui->graphicsViewHead->geometry();
    this->mInsetScene->setSceneRect(bound);
    this->mInsetScene->addPixmap(pixMap);

    this->mImageInset->setPixmap(pixMap);
    this->mImageInset->setPos(bound.topLeft().x() ,bound.topLeft().y());

   // this->ui->graphicsViewHead->show();

}


//Write To Message List Below Tracker
void MainWindow::LogEvent(QString strMessage)
{

    mMessageList.append(QString::number(nFrame) + "." + strMessage);
    this->ui->listView->show();
    mModelMessageList->setStringList(mMessageList);


}

void MainWindow::showCVimg(cv::Mat& img)
{
    frameScene = img;
    qimg = QtOpencvCore::img2qimg(img);

    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimg);

    //Removed Scaling
    // scale pixMap image to fit into the QLabel
    //pixMap = pixMap.scaled(this->ui->graphicsView->size(), Qt::KeepAspectRatio);
    //this->mImage->setPixmap(pixMap);
    //this->ui->graphicsView->setSceneRect(this->frameGeometry()); // set the scene's bounding rect to rect of mainwindow

    /// A Scene contains graphic objects, but rendering them requires a View.
    ///  The Scenes size can be larger than the View's
    this->ui->graphicsView->setSceneRect(this->ui->graphicsView->geometry()); // set the scene's bounding rect to rect of mainwindow
    QRect bound = this->ui->graphicsView->geometry();
    this->mScene->setSceneRect(bound);

    //this->ui->graphicsView->setFixedSize(qimg.width(),qimg.height());

    mImage->setPixmap(pixMap);
    //this->ui->graphicsView->fitInView(this->mImage->boundingRect(), Qt::KeepAspectRatio);

    this->mImage->setPos(bound.topLeft().x() ,bound.topLeft().y());


    //this->ui->graphicsView->fitInView(mImage, Qt::KeepAspectRatio);
    //this->ui->graphicsView->show();

    this->mpLastCVImg = &img; //Save Pointer to frame
    //mImage

}

void MainWindow::textEdited(QString strFrame)
{
    qDebug() << "Txt Edited to:" << strFrame;
}

void MainWindow::echoChanged(int i)
{
    qDebug() << "Echo bx" ;
}

void MainWindow::changeEvent(QEvent *e)
{
    qDebug() << "Change E" << e->type() ;
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event) {
    char key = 0;
    if (event->type() == QEvent::GraphicsSceneWheel) {
        handleWheelOnGraphicsScene(dynamic_cast<QGraphicsSceneWheelEvent*> (event));

        // Don't propagate
        event->accept();
        return true;
    }

    if (event->type() == QEvent::KeyPress) {
         QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
         std::string strkey = keyEvent->text().toStdString();
         if (strkey.length() > 0)
            key =  strkey.at(0);

         qDebug() << "Ate key press " << keyEvent->text().toStdString().c_str() << " k: " << key << " from " << obj->objectName();
        ///Catch Frame Number Edit Enter Press
         if (obj == ui->spinBoxFrame && (event->type() == QEvent::KeyPress) )
         {
             QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);

             //Check If Enter Pressed - Then Change Start Frame Number
             if(Qt::Key_Enter == keyEvent->key() || keyEvent->key() == Qt::Key_Return )
             {
                 qDebug() << "Enter pressed - Seek New Frame";
                 nFrame = ui->spinBoxFrame->value();
                 this->ui->horizontalSlider->setValue( ui->spinBoxFrame->value());

                 bStartFrameChanged = true;
                 event->accept();

                 return true;// We VE handled the event
             }else
             {
                 return false; //Pass THe event To its intented Receipient

             }



             //If Keypress Not from TxtFrameBox
         }else //Propagate KeyPress to Prog Flow Control
         {
             ::keyCommandFlag(this,key,nFrame);
         }
         return true;
     }

    if (event->type() == QEvent::Resize) {
         QResizeEvent *resizeEvent = static_cast<QResizeEvent*>(event);
         //qDebug(" Resized (New Size) - Width: %d Height: %d",                resizeEvent->size().width(),                resizeEvent->size().height());



          //ui->graphicsView->setFixedSize(resizeEvent->size()*0.9);

         //QRectF bounds = mScene->itemsBoundingRect();
         //bounds.setWidth(bounds.width()*0.9);         // to tighten-up margins
         //bounds.setHeight(bounds.height()*0.9);       // same as above

          //ui->graphicsView->fitInView(bounds, Qt::KeepAspectRatio);
        // ui->graphicsView->centerOn(0, 0);


        // this->ui->gridLayout->setGeometry(this->frameGeometry());
         //this->ui->graphicsView->setSceneRect(this->frameGeometry()); // set the scene's bounding rect to rect of mainwindow
        return true;
    }

    if (event->type() == QEvent::GraphicsSceneDragMove) {
        dragMoveEvent(dynamic_cast<QGraphicsSceneDragDropEvent*> (event));
        return true;
    }

    //Detect Drag During Pause
    if (event->type() == QEvent::GraphicsSceneMouseMove) {
        mouseMoveEvent(dynamic_cast<QGraphicsSceneMouseEvent*> (event));
        return true;
    }

    if (event->type() == QEvent::GraphicsSceneMouseDoubleClick)
    {
        mouseDblClickEvent(dynamic_cast<QGraphicsSceneMouseEvent*> (event));

        return true;

    }

    if (event->type() == QEvent::GraphicsSceneMousePress)
    {
        mousePressEvent(dynamic_cast<QGraphicsSceneMouseEvent*> (event));
        return true;
    }

    if (event->type() == QEvent::GraphicsSceneMouseRelease)
    {
        mouseReleaseEvent(dynamic_cast<QGraphicsSceneMouseEvent*> (event));

        qDebug() << "Mouse Up";
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



  //this->mImage->transform().scale(this->mImage->scale()+ steps,this->mImage->scale()+ steps);
}



void MainWindow::dragMoveEvent(QGraphicsSceneDragDropEvent* mouseEvent )
{

    qDebug() << "Drag Event : " <<  mouseEvent->pos().y();
}

void MainWindow::mouseMoveEvent ( QGraphicsSceneMouseEvent* mouseEvent )
{

    //qDebug() << "Mouse Mv";
    if (bDraggingTemplateCentre ) //bDraggingTemplateCentre
    {
         //qDebug() << "Dragging";


        QPointF ptSceneclick = mouseEvent->scenePos();// this->ui->graphicsView->mapToScene( mouseEvent->pos().x(),mouseEvent->pos().y() );
        // get the item that was clicked on
        QGraphicsItem* item = mScene->itemAt( ptSceneclick, this->ui->graphicsView->transform() );
        if (!item)
            return;
        // get the scene pos in the item's local coordinate space
        QPointF ptImg = item->mapFromScene(ptSceneclick);

        cv::Point ptMouse(ptImg.x(),ptImg.y());
        for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
        {

            fishModel* fish = (*it).second;
            if (fish->bodyRotBound.boundingRect().contains(ptMouse)) //Clicked On Fish Box
            {
                //bDraggingTemplateCentre = true; //Still True /Terminates Upon Click

                qDebug() << "Drag to  pos x: " << ptMouse.x << " y:" << ptMouse.y;
                fish->bodyRotBound.center = ptMouse;
                fish->ptRotCentre         = ptMouse;
                ///Draw a Red Rotated Frame around Detected Body
                cv::Point2f boundBoxPnts[4];
                fish->bodyRotBound.points(boundBoxPnts);
                 for (int j=0; j<4;j++) //Rectangle Body
                   cv::line(frameScene,boundBoxPnts[j],boundBoxPnts[(j+1)%4] ,CV_RGB(00,00,255),1,cv::LINE_8);

                //showCVimg(frameScene);

            }
        }//For eAch Fishs

} //Check For Mouse Down And Mouse Moving - Dragging




}//Mouse Move Event


///
/// \brief MainWindow::mousePressEvent
/// Two Cases - Either A Move the Fish Bounding Window Initiates
/// Or It Terminates a previous one - Saving The new template
/// \param mouseEvent
///

void MainWindow::mousePressEvent ( QGraphicsSceneMouseEvent * mouseEvent )
{
    //Already Dragging - Terminate And Save Template
    if (bDraggingTemplateCentre)
    {
        bStoreThisTemplate = true;
        bDraggingTemplateCentre = false;
        //bSceneMouseLButtonDown = false; //Act as Button Released With This Press
        qDebug() << "Store New Template position";
    }
    else
    {
        //bPaused = true;
        bSceneMouseLButtonDown = true;
        qDebug() << "Mouse Down";

    }


}

void MainWindow::mouseReleaseEvent( QGraphicsSceneMouseEvent * mouseEvent )
{
    bSceneMouseLButtonDown = false;
    //bDraggingTemplateCentre = false;
    qDebug() << "Mouse Up";
}

///
/// \brief MainWindow::mouseDblClickEvent Start Dragging Bounding Box Of Fish
/// \param mouseEvent
///
void MainWindow::mouseDblClickEvent( QGraphicsSceneMouseEvent * mouseEvent )
{

    QPointF ptSceneclick = mouseEvent->scenePos();// this->ui->graphicsView->mapToScene( mouseEvent->pos().x(),mouseEvent->pos().y() );
    // get the item that was clicked on
    QGraphicsItem* item = mScene->itemAt( ptSceneclick, this->ui->graphicsView->transform() );
    if (!item)
        return;
    // get the scene pos in the item's local coordinate space
    QPointF ptImg = item->mapFromScene(ptSceneclick);

    cv::Point ptMouse(ptImg.x(),ptImg.y());
    for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
    {

        fishModel* fish = (*it).second;
        if (fish->bodyRotBound.boundingRect().contains(ptMouse)) //Clicked On Fish Box
        {
            bDraggingTemplateCentre = true;
            qDebug() << "Start Dragging Fish Bound from position x: " << ptMouse.x << " y:" << ptMouse.y;

        }
    }


}



MainWindow::~MainWindow()
{
    delete ui;
}
