#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "larvatrack.h" //For resetDataRecording()
#include "QtOpencvCore.hpp"
#include <QStringListModel>
#include <qlineedit.h>

extern QFile outfishdatafile;
extern QFile outfooddatafile;

extern fishModels vfishmodels; //Vector containing live fish models
extern foodModels vfoodmodels; //Vector containing live fish models
extern bool bPaused;
extern bool bStoreThisTemplate;
extern bool bDraggingTemplateCentre;
extern bool bStartFrameChanged;
extern int g_Segthresh;
extern int gthresEyeSeg;
extern int gi_maxEllipseMajor;
extern int gi_minEllipseMajor;
extern QElapsedTimer gTimer;
extern ltROIlist vRoi;
extern int gFishBoundBoxSize;
extern double gTemplateMatchThreshold;
extern double gdMOGBGRatio;
extern bool bTrackFood;
extern bool bTracking;
extern bool bExiting;
extern bool bRecordToFile;

bool bSceneMouseLButtonDown;
bool bDraggingRoiPoint;


extern cv::Ptr<cv::BackgroundSubtractorMOG2> pMOG2; //MOG2 Background subtractor

cv::Point* ptDrag;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ptDrag = 0;
    this->nFrame = 0;
    this->mScene = new QGraphicsScene(this->ui->graphicsView);
    this->mInsetScene = new QGraphicsScene(this->ui->graphicsViewHead);
    this->mInsetTemplateScene = new QGraphicsScene(this->ui->graphicsViewTemplate);

    this->ui->graphicsView->setScene(this->mScene);
    this->ui->graphicsViewHead->setScene(this->mInsetScene);
    this->ui->graphicsViewTemplate->setScene(this->mInsetTemplateScene);

    //this->ui->graphicsView->setFixedSize(1280,1024);
    //mScene->setSceneRect(this->ui->graphicsView->rect());

    mScene->addText("Zebrafish Tracker Scene")->setPos(100,100);
    //Add Empty/New PixMap on Which we will set the images onto
    QPixmap pxMapEmpty1(this->ui->graphicsView->geometry().width(),this->ui->graphicsView->geometry().height());
    QPixmap pxMapEmpty2(this->ui->graphicsViewHead->geometry().width(),this->ui->graphicsViewHead->geometry().height());
    QPixmap pxMapEmpty3(this->ui->graphicsViewTemplate->geometry().width(),this->ui->graphicsViewTemplate->geometry().height());
    this->mImage                 = mScene->addPixmap(pxMapEmpty1);
    this->mImageInset            = mInsetScene->addPixmap(pxMapEmpty2);
    this->mImageTemplateInset    = mInsetTemplateScene->addPixmap(pxMapEmpty3);

    this->ui->graphicsView->setSceneRect(this->ui->graphicsView->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->ui->graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    this->ui->graphicsViewHead->setSceneRect(this->ui->graphicsViewHead->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsViewHead->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->ui->graphicsViewHead->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);


    this->ui->graphicsViewTemplate->setSceneRect(this->ui->graphicsViewTemplate->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsViewTemplate->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->ui->graphicsViewTemplate->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);


    this->ui->horizontalSlider->setRange(0,50000);
    this->mScene->installEventFilter(this);

    this->installEventFilter(this); //To Capture Resize

    // Log Events on List View //
    mModelMessageList = new QStringListModel(this);
    // Populate our model
    mModelMessageList->setStringList(mMessageList);
    // Glue model and view together
    this->ui->listView->setModel(mModelMessageList);



    this->ui->horizontalSlider->installEventFilter(this);

    createSpinBoxes();
    nFrame = 0;
}

void MainWindow::createSpinBoxes()
{

    this->ui->spinBoxFrame->installEventFilter(this);

    this->ui->spinBoxEyeThres->installEventFilter(this); //-Ve Values Allow for lowering Avg Threshold
    this->ui->spinBoxEyeThres->setRange(-100,400); //-Ve Values Allow for lowering Avg Threshold
    this->ui->spinBoxEyeThres->setValue(gthresEyeSeg);



    this->ui->spinBoxFishThres->installEventFilter(this);
    this->ui->spinBoxFishThres->setRange(19,100); //Too low Below 19 App Stalls -- Too many Large Objects Appear
    this->ui->spinBoxFishThres->setValue(g_Segthresh);



    this->ui->spinBoxMinEllipse->installEventFilter(this);
    this->ui->spinBoxMinEllipse->setRange(5,22);
    this->ui->spinBoxMinEllipse->setValue(gi_minEllipseMajor);


    this->ui->spinBoxMaxEllipse->installEventFilter(this);
    this->ui->spinBoxMaxEllipse->setRange(15,35);
    this->ui->spinBoxMaxEllipse->setValue(gi_maxEllipseMajor);


    this->ui->spinBoxSpineSegSize->installEventFilter(this);
    this->ui->spinBoxSpineSegSize->setRange(2,20);
    this->ui->spinBoxSpineSegSize->setValue(gFishTailSpineSegmentLength);

    //These spinBoxes Use Slots For Events
    this->ui->spinBoxTemplateThres->setValue(gTemplateMatchThreshold*100.0);

    this->ui->spinBoxMOGBGRatio->setValue(gdMOGBGRatio*100.0);


    //this->connect(this->ui->spinBoxEyeThres, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),this->ui->spinBoxEyeThres, &QSlider::setValue);

//    QObject::connect(this->ui->spinBoxEyeThres,
//                     static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
//                     this->ui->spinBoxEyeThres,
//                     static_cast<void (MainWindow::*)(int)>(&MainWindow::valueChanged));


    ///\todo Remove These and Make Use of Auto Slots via form editor just like the spinBoxMOGBGRatio or spinBoxTemplateThres
    QObject::connect(this->ui->spinBoxEyeThres,
                     SIGNAL(valueChanged(int)),
                     this,
                     SLOT(eyevalueChanged(int)));

    QObject::connect(this->ui->spinBoxFishThres,
                     SIGNAL(valueChanged(int)),
                     this,
                     SLOT(fishvalueChanged(int)));

    QObject::connect(this->ui->spinBoxSpineSegSize,
                     SIGNAL(valueChanged(int)),
                     this,
                     SLOT(tailSizevalueChanged(int)));


    QObject::connect(this->ui->spinBoxMinEllipse,
                     SIGNAL(valueChanged(int)),
                     this,
                     SLOT(minEllipseSizevalueChanged(int)));


    QObject::connect(this->ui->spinBoxMaxEllipse,
                     SIGNAL(valueChanged(int)),
                     this,
                     SLOT(maxEllipseSizevalueChanged(int)));


}

void MainWindow::showVideoFrame(cv::Mat& img,unsigned int nFrame)
{
    //this->ui->horizontalSlider->setValue(nFrame);

    showCVimg(img);
}

void MainWindow::saveScreenShot()
{
    //std::stringstream frameNumberString; frameNumberString << nFrame;
    QString frameNumberString = QString::number(this->nFrame);

    ::saveImage(frameNumberString,stroutDirCSV,vidFilename,this->frameScene);

}

//Saves Selected Template Images From Running Video To Special templates Subfolder for future Re-Use
void MainWindow::saveTemplateImage(cv::Mat& imgTempl)
{
    //std::stringstream frameNumberString; frameNumberString << nFrame;
    QString dirToSave = stroutDirCSV;
    QString frameNumberString = "Templ_" + QString::number(nFrame);

    //Make ROI dependent File Name
    QFileInfo fiVid(vidFilename);
    QString fileVidCoreName = "templ_"+ fiVid.completeBaseName();

    dirToSave.append("/templates/");
    QString imageToSave =  fileVidCoreName + "_" + frameNumberString + ".pgm";
    imageToSave.prepend(dirToSave);

    if (!QDir(dirToSave).exists())
    {
        std::clog << "Make directory " << dirToSave.toStdString() << std::endl;
        QDir().mkpath(dirToSave);
    }

    bool saved = cv::imwrite(imageToSave.toStdString(), imgTempl);


}

void MainWindow::setTotalFrames(uint FrameCount)
{
    this->ui->horizontalSlider->setMaximum(FrameCount);
    this->nTotalFrameCount = FrameCount;
}

void MainWindow::tickProgress()
{
    //this->ui->horizontalSlider->setValue(this->ui->horizontalSlider->value()+1);
    this->ui->horizontalSlider->setValue(nFrame);
    this->ui->spinBoxFrame->setValue(nFrame);
}


//Write To Message List Below Tracker
void MainWindow::LogEvent(QString strMessage)
{
    //\note THe frame Number seemed to be advanced by 1, in the log file
    mMessageList.append(QString::number(gTimer.elapsed()/60000,'g',4) + " " + QString::number(nFrame-1) + "# " + strMessage);
    this->ui->listView->show();
    mModelMessageList->setStringList(mMessageList);

    std::clog << gTimer.elapsed()/60000 << " #" << nFrame-1 << " " << strMessage.toStdString() << std::endl;

    try{
        //if ((gTimer.elapsed()/60000) % 3 == 0 ) //Scroll To Bottom Every 3 sec
//            this->ui->listView->scrollToBottom();
    }
    catch (char* e)
    {
        std::cerr << gTimer.elapsed()/60000 << " [Error] listView->scrollToBottom(); error: "<< e << std::endl;
    }
}

void MainWindow::showInsetimg(cv::Mat& img)
{

    if (img.cols == 0 || img.rows == 0)
        return;

    qimgHead = QtOpencvCore::img2qimg(img);
    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimgHead);
    QRect bound = this->ui->graphicsViewHead->geometry();
    this->mInsetScene->setSceneRect(bound);

    this->mImageInset->setPixmap(pixMap);
    this->mImageInset->setPos(bound.topLeft().x() ,bound.topLeft().y());

   // this->ui->graphicsViewHead->show();

}

void MainWindow::showInsetTemplateimg(cv::Mat& img)
{

    if (img.cols == 0 || img.rows == 0)
        return;

    qimgHead = QtOpencvCore::img2qimg(img);
    // convert the opencv image to a QPixmap (to show in a QLabel)
    QPixmap pixMap = QPixmap::fromImage(qimgHead);
    QRect bound = this->ui->graphicsViewTemplate->geometry();
    this->mInsetTemplateScene->setSceneRect(bound);

    this->mImageTemplateInset->setPixmap(pixMap);
    this->mImageTemplateInset->setPos(bound.topLeft().x() ,bound.topLeft().y());

}

void MainWindow::showCVimg(cv::Mat& img)
{
    frameScene.release();
    //frameScene = img.clone();
    img.copyTo(frameScene);

    /// Draw / Overlay Info From This Window //
    if (ptDrag)
        cv::circle(frameScene,*ptDrag,5,cv::Scalar(200,200,0),2);

    /////

    qimg = QtOpencvCore::img2qimg(frameScene);


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

    //this->mpLastCVImg = &img; //Save Pointer to frame
    //mImage


}


void MainWindow::tailSizevalueChanged(int i)
{
    qDebug() << "Tails SpinBox gave " << i;
    LogEvent(QString("Tail Segment Size changed:") + QString::number(i));
    gFishTailSpineSegmentLength = i;
}

void MainWindow::eyevalueChanged(int i)
{
    //qDebug() << "Eye SpinBox gave " << i;

    if (!bSceneMouseLButtonDown)
        LogEvent(QString("changed Eye Seg Threshold:") + QString::number(i));

    gthresEyeSeg = i;
}

void MainWindow::fishvalueChanged(int i)
{
    qDebug() << "fish SpinBox gave " << i;
    LogEvent(QString("Changed Fish BG Threshold:") + QString::number(i));
    g_Segthresh = i;
}

void MainWindow::maxEllipseSizevalueChanged(int i)
{
    gi_maxEllipseMajor = i;
    LogEvent(QString("changed Max Ellipse to:") + QString::number(i));
}

void MainWindow::minEllipseSizevalueChanged(int i)
{
    gi_minEllipseMajor = i;
    LogEvent(QString("changed Min Ellipse changed to:") + QString::number(i));
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
    //qDebug() << "Change E" << e->type() ;
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event) {
    char key = 0;

    if (obj->objectName() == "horizontalSlider")
    {
        handleSliderChange(event);

    }
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

         //Cancel Any Drag Event Going On
         if (bDraggingTemplateCentre)
         {
            bDraggingTemplateCentre = false;
            LogEvent("[info] Cancelled Template Adjustment");
         }


         //qDebug() << "Ate key press " << keyEvent->text().toStdString().c_str() << " k: " << key << " from " << obj->objectName();
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
                 bPaused = false; //For 1 Frame and it will pause again at new frame
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


//    if (obj == ui->spinBoxEyeThres)
//        gthresEyeSeg =  ui->spinBoxEyeThres->value();

//    if (obj == ui->spinBoxFishThres)
//        g_Segthresh = ui->spinBoxFishThres->value();

//    if (obj == ui->spinBoxMaxEllipse)
//         gi_maxEllipseMajor = ui->spinBoxMaxEllipse->value();

//    if (obj == ui->spinBoxMinEllipse)
//         gi_minEllipseMajor = ui->spinBoxMinEllipse->value();

//    if (obj == ui->spinBoxSpineSegSize)
//        gFishTailSpineSegmentLength = ui->spinBoxSpineSegSize->value();

    return false;
}

void MainWindow::handleSliderChange(QEvent* event)
{

    if (event->type() == QEvent::MouseButtonPress)
    {
        bPaused = true;
    }

    if (event->type() == QEvent::MouseButtonRelease)
    {
        bStartFrameChanged = true;
        bPaused = false;
        nFrame = this->ui->horizontalSlider->value();
        this->ui->spinBoxFrame->setValue(nFrame);
    }
    if (event->type() == QEvent::HoverMove)
    {
        this->ui->horizontalSlider->setToolTip(QString::number(this->ui->horizontalSlider->value()));
    }
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

    setCursor(Qt::ArrowCursor);

    QPointF ptSceneclick = mouseEvent->scenePos();
    QGraphicsItem* item = mScene->itemAt( ptSceneclick, this->ui->graphicsView->transform() );
    if (!item)
        return;

    // get the scene pos in the item's local coordinate space
    QPointF ptImg = item->mapFromScene(ptSceneclick);
    cv::Point ptMouse(ptImg.x(),ptImg.y());

    //qDebug() << "Mouse Mv";
    if (bDraggingTemplateCentre ) //bDraggingTemplateCentre
    {
         //qDebug() << "Dragging";
        // this->ui->graphicsView->mapToScene( mouseEvent->pos().x(),mouseEvent->pos().y() );
        // get the item that was clicked on

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


    if (bDraggingRoiPoint)
    {   //Update Point - Bound it from Periphery - OtherWise TemplateMatch Fails due to Small Image Crop at boundary
        ptDrag->x = std::max(gFishBoundBoxSize, std::min(frameScene.cols - gFishBoundBoxSize,ptMouse.x));
        ptDrag->y = std::max(gFishBoundBoxSize, std::min(frameScene.rows - gFishBoundBoxSize,ptMouse.y));
    }
    else
        ptDrag = 0;


    //Check If Around ROI Points
    for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {
        ltROI* iroi = &(*it);
        for (std::vector<cv::Point>::iterator it = iroi->vPoints.begin() ; it != iroi->vPoints.end(); ++it)
        {
            //4 pixels Around the ROI Point / Get Moving Cursor
            if (cv::norm(*it-ptMouse) < 4)
            {
                    // HighLight Point //
//                    cv::circle(frameScene,*it,5,cv::Scalar(100,200,0),2);
                    ptDrag = &(*it);
                    setCursor(Qt::CrossCursor);
            }
        }
    }





}//Mouse Move Event


///
/// \brief MainWindow::mousePressEvent
/// Two Cases - Either A Move the Fish Bounding Window Initiates
/// Or It Terminates a previous one - Saving The new template
/// \param mouseEvent
///

void MainWindow::mousePressEvent ( QGraphicsSceneMouseEvent * mouseEvent )
{

    QPointF ptSceneclick = mouseEvent->scenePos();
    QGraphicsItem* item = mScene->itemAt( ptSceneclick, this->ui->graphicsView->transform() );
    if (!item)
        return;

    // get the scene pos in the item's local coordinate space
    QPointF ptImg = item->mapFromScene(ptSceneclick);
    cv::Point ptMouse(ptImg.x(),ptImg.y());

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


    //Check If Around ROI Points
    for (std::vector<ltROI>::iterator it = vRoi.begin(); it != vRoi.end(); ++it)
    {

        ltROI* iroi = &(*it);

        for (std::vector<cv::Point>::iterator it = iroi->vPoints.begin() ; it != iroi->vPoints.end(); ++it)
        {
            //cv::Point  ptR = *it;

            //4 pixels Around the ROI Point / Get Moving Cursor
            if (cv::norm(*it-ptMouse) < 4)
            {
                    //setCursor(Qt::PointingHandCursor);
                    setCursor(Qt::CrossCursor);
                    //If Button Down , THen Move Point
                    if (bSceneMouseLButtonDown)
                    {
                        ptDrag = &(*it); //*Returns Ref, and then AddressOf Operator gives Pointer

                        ptDrag->x = ptMouse.x;
                        ptDrag->y = ptMouse.y;
                        bDraggingRoiPoint = true;
                    }
            }

        }
    }



}

void MainWindow::mouseReleaseEvent( QGraphicsSceneMouseEvent * mouseEvent )
{
    bSceneMouseLButtonDown = false;
    bDraggingRoiPoint = false;

    setCursor(Qt::ArrowCursor);
    ptDrag = 0; //Empty the Dragged Point Pointer

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
            LogEvent("[info] Adjust Template from position ON- Start Dragging");
            //this->statusBar()->set
            qDebug() << "Start Dragging Fish Bound from position x: " << ptMouse.x << " y:" << ptMouse.y;

        }
    }

    ///Check Clicking On Food Item
    for (foodModels::iterator it=vfoodmodels.begin(); it!=vfoodmodels.end(); ++it)
    {

        foodModel* food = (*it).second;
        if (food->zTrack.boundingBox.contains(ptMouse) ) //Clicked On Fish Box
        //if (cv::norm((cv::Point) food->zTrack.centroid - ptMouse) < 5 ) //Clicked On Fish Box
        {
            // Make Targeted
            if (!food->isTargeted)
            {
                food->isTargeted = true;
                qDebug() << "Food Targetting On  x: " << ptMouse.x << " y:" << ptMouse.y;
                LogEvent("[info] Begin Tracking Food Item");
            }else
            {
               LogEvent("[info] END Tracking Food Item");
                food->isTargeted = false;
            }

        }
    }



}



MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_spinBoxTemplateThres_valueChanged(int arg1)
{
 double newTMatchThresh = (double)arg1/100.0;
 gTemplateMatchThreshold = newTMatchThresh;
 LogEvent(QString("[info] Changed Template Match Thres:" ) + QString::number(newTMatchThresh,'g',4) ) ;
}

void MainWindow::on_spinBoxMOGBGRatio_valueChanged(int arg1)
{
    double newBGRatio = (double)arg1/100.0;
    gdMOGBGRatio = newBGRatio; //Updated Value Takes effect in processFrame and in updateBGFrame

    if (pMOG2)
    {
       pMOG2->setBackgroundRatio(gdMOGBGRatio);
       LogEvent(QString("[info] Changed MOG BG Ratio: " ) + QString::number(gdMOGBGRatio,'g',4) ) ;

    }


}

void MainWindow::on_actionTrack_Fish_triggered(bool checked)
{

    if (!bTracking)
    {
        //iLastKnownGoodTemplateRow = 0; //Reset Row
        //iLastKnownGoodTemplateCol = 0;
        LogEvent(QString("Tracking ON"));
    }else
        LogEvent(QString("Tracking OFF"));

    bTracking = checked;


}

void MainWindow::on_actionTrack_Food_triggered(bool checked)
{

    bTrackFood = checked;
}

void MainWindow::on_actionRecord_Tracks_to_File_w_triggered(bool checked)
{
    bRecordToFile = checked;
    if (bRecordToFile)
    {
      LogEvent(QString(">> Recording Tracks ON - New File <<"));

      resetDataRecording(outfishdatafile);
      writeFishDataCSVHeader(outfishdatafile);
      resetDataRecording(outfooddatafile);
      writeFoodDataCSVHeader(outfooddatafile);

    }
    else
      LogEvent(QString("<< Recording Tracks OFF >>"));

}

void MainWindow::on_actionQuit_triggered()
{
    bExiting = true;
    LogEvent("[info] User Terminated - Bye!");
}



void MainWindow::on_actionStart_tracking_triggered()
{
    if (bPaused)
        LogEvent("[info] Running");

    bPaused = false;

}

void MainWindow::on_actionPaus_tracking_p_triggered()
{
    bPaused = true;
    if (bPaused)
    LogEvent("[info] Paused");

}

void MainWindow::on_actionPaus_tracking_p_triggered(bool checked)
{
}
