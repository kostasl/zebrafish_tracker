#include "mainwindow.h"
#include "ui_mainwindow.h"


#include "config.h"
#include <QtMath>
#include <QScrollBar>
#include <QEvent>
#include "larvatrack.h" //For resetDataRecording()
#include "QtOpencvCore.hpp"
#include <QStringListModel>
#include <qlineedit.h>

extern QFile outfishdatafile;
extern QFile outfooddatafile;

extern fishModels vfishmodels; //Vector containing live fish models
extern foodModels vfoodmodels; //Vector containing live fish models
extern pointPairs vMeasureLines; //vector of point pairs/ user defined lines

extern trackerState gTrackerState;

std::pair<cv::Point,cv::Point> userPointPair; //The currently defined point pair prior to adding to list

extern QElapsedTimer gTimer;

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
    mScene->addText("Zebrafish Tracker Scene")->setPos(100,100);

    this->mInsetScene = new QGraphicsScene(this->ui->graphicsViewHead);
    this->mInsetTemplateScene = new QGraphicsScene(this->ui->graphicsViewTemplate);

    this->ui->graphicsView->setScene(this->mScene);
    this->ui->graphicsViewHead->setScene(this->mInsetScene);
    this->ui->graphicsViewTemplate->setScene(this->mInsetTemplateScene);

    this->ui->graphicsViewHead->setSceneRect(this->ui->graphicsViewHead->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsViewTemplate->setSceneRect(this->ui->graphicsViewTemplate->geometry()); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsView->setSceneRect(QRect(0,0,3000,3000)); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsView->horizontalScrollBar()->setValue(1500); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsView->verticalScrollBar()->setValue(1500); // set the scene's bounding rect to rect of mainwindow
    this->ui->graphicsView->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    this->ui->graphicsView->setResizeAnchor(QGraphicsView::AnchorUnderMouse);

    //this->ui->graphicsView->setFixedSize(1280,1024);


    //Add Empty/New PixMap on Which we will set the images onto
    QPixmap pxMapEmpty1(250,250);
    QPixmap pxMapEmpty2(this->ui->graphicsViewHead->geometry().width(),this->ui->graphicsViewHead->geometry().height());
    QPixmap pxMapEmpty3(this->ui->graphicsViewTemplate->geometry().width(),this->ui->graphicsViewTemplate->geometry().height());
    this->mImage                 = mScene->addPixmap(pxMapEmpty1);
    this->mImage->setPos(100,100);
    this->mImageInset            = mInsetScene->addPixmap(pxMapEmpty2);
    this->mImageTemplateInset    = mInsetTemplateScene->addPixmap(pxMapEmpty3);


    //this->ui->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded); //leave option down to  Form Editor
        //this->ui->graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);


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

    this->ui->actionTrack_Fish->setChecked(gTrackerState.bTracking);
    this->ui->checkBoxGPU->setChecked(gTrackerState.bUseGPU);
    this->ui->checkBoxMOG->setChecked(gTrackerState.bUseBGModelling);
    this->ui->checkBoxNoiseFilter->setChecked(gTrackerState.bRemovePixelNoise);

    this->ui->checkBoxHistEqualizer->setChecked(gTrackerState.bUseHistEqualization);

    createSpinBoxes();
    nFrame = 0;
    //Reset Point Pair
    userPointPair.first.x = userPointPair.first.y = -10;
    userPointPair.second.x = userPointPair.second.y = -10;
}

void MainWindow::createSpinBoxes()
{

    this->ui->spinBoxFrame->installEventFilter(this);

    this->ui->spinBoxEyeThres->installEventFilter(this); //-Ve Values Allow for lowering Avg Threshold
    this->ui->spinBoxEyeThres->setRange(-500,500); //-Ve Values Allow for lowering Avg Threshold
    this->ui->spinBoxEyeThres->setValue(gTrackerState.gthresEyeSeg);

    this->ui->spinBoxFoodThresMax->setValue(gTrackerState.g_SegFoodThesMax);
    this->ui->spinBoxFoodThresMin->setValue(gTrackerState.g_SegFoodThesMin);



    this->ui->spinBoxFishThres->installEventFilter(this);
    this->ui->spinBoxFishThres->setRange(1,100); //Too low Below 19 App Stalls -- Too many Large Objects Appear
    this->ui->spinBoxFishThres->setValue(gTrackerState.g_Segthresh);



    this->ui->spinBoxMinEllipse->installEventFilter(this);
    this->ui->spinBoxMinEllipse->setRange(5,22);
    this->ui->spinBoxMinEllipse->setValue(gTrackerState.gi_minEllipseMajor);


    this->ui->spinBoxMaxEllipse->installEventFilter(this);
    this->ui->spinBoxMaxEllipse->setRange(15,35);
    this->ui->spinBoxMaxEllipse->setValue(gTrackerState.gi_maxEllipseMajor);


//    this->ui->spinBoxSpineSegSize->installEventFilter(this);
//    this->ui->spinBoxSpineSegSize->setRange(2,20);
//    this->ui->spinBoxSpineSegSize->setValue(gFishTailSpineSegmentLength);

    //These spinBoxes Use Slots For Events
    this->ui->spinBoxTemplateThres->setValue(gTrackerState.gTemplateMatchThreshold*100.0);

    this->ui->spinBoxMOGBGRatio->setValue(gTrackerState.gdMOGBGRatio*100.0);


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

//    QObject::connect(this->ui->spinBoxSpineSegSize,
//                     SIGNAL(valueChanged(int)),
//                     this,
//                     SLOT(tailSizevalueChanged(int)));


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
        if ((gTimer.elapsed()/60000) % 3 == 0 ) //Scroll To Bottom Every 3 sec
            this->ui->listView->scrollToBottom();
    }
    catch (char* e)
    {
        std::cerr << gTimer.elapsed()/60000 << " [Error] listView->scrollToBottom(); error: "<< e << std::endl;
    }
}

void MainWindow::SetTrackerState(int stateID)
{

    switch (stateID)
    {
        case 0: //paused
        {
            this->statusBar()->showMessage(tr("Paused, press r to resume"));
            gTrackerState.bPaused = true;
        }
        break;
        case 1:
            this->statusBar()->showMessage(tr("Tracking - press p to pause"));
        break;
        case 5:
            this->statusBar()->showMessage(tr("Click to manually set new prey item to track "));
        break;
        case 6:
            {
                this->statusBar()->showMessage(tr("Measure mode, click on 1st source point of measurement"));
             }
        break;
        case 7:
            this->statusBar()->showMessage(tr("Click on 2nd point of measurement"));
        break;
        case 8:
              this->statusBar()->showMessage(tr("New point pair added. "));
        break;
        default:
            this->statusBar()->showMessage(tr("Tracking - press p to pause"));
    }
}

///
/// \brief getFoodItemAtLocation Return Pointer to 1st food item found at clicked (mouse) location
/// \param ptLocation
/// \return
///
preyModel* MainWindow::getFoodItemAtLocation(cv::Point ptLocation)
{
    preyModel* rfood = NULL;

    ///Check First if Clicking On Food Item
    for (foodModels::iterator it=vfoodmodels.begin(); it!=vfoodmodels.end(); ++it)
    {

        preyModel* pfood = (*it).second;
        //if (pfood->zTrack.boundingBox.contains(ptLocation) ) //Clicked On Fish Box
        if (cv::norm(pfood->zfoodblob.pt - (cv::Point2f)ptLocation) < gTrackerState.gMaxClusterRadiusFoodToBlob ) //Clicked On Fish Box
        {
            rfood = pfood; //Found and return clicked Food Item
        }
    }
    return rfood;
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
    QRectF bound = this->ui->graphicsView->sceneRect(); // set the scene's bounding rect to rect of mainwindow

    //this->ui->graphicsView->setFixedSize(qimg.width(),qimg.height());

    mImage->setPixmap(pixMap);
    QRectF imgbound = this->mImage->boundingRect();
    //bound.setHeight(bound.height()*2);
    //bound.setWidth(bound.width()*2);
    //this->mScene->setSceneRect(bound);



    this->mImage->setPos(bound.center().x()-imgbound.width()/2,
                         bound.center().y()-imgbound.height()/2);

    //this->ui->graphicsView->fitInView(this->mImage->boundingRect(), Qt::KeepAspectRatio);

    //this->mpLastCVImg = &img; //Save Pointer to frame
    //mImage


}



void MainWindow::UpdateSpinBoxToValue()
{

    this->ui->spinBoxEyeThres->setValue(gTrackerState.gthresEyeSeg);
    this->ui->spinBoxFishThres->setValue(gTrackerState.g_Segthresh);

    this->ui->spinBoxFoodThresMax->setValue(gTrackerState.g_SegFoodThesMax);
    this->ui->spinBoxFoodThresMin->setValue(gTrackerState.g_SegFoodThesMin);

}

void MainWindow::eyevalueChanged(int i)
{
    //qDebug() << "Eye SpinBox gave " << i;

    if (bSceneMouseLButtonDown)
        LogEvent(QString("Changed Eye Seg Threshold:") + QString::number(i));

    gTrackerState.gthresEyeSeg = i;
}

void MainWindow::fishvalueChanged(int i)
{
    qDebug() << "fish SpinBox gave " << i;
    LogEvent(QString("Changed Fish BG Threshold:") + QString::number(i));
    gTrackerState.g_Segthresh = i;

 }

void MainWindow::maxEllipseSizevalueChanged(int i)
{
    gTrackerState.gi_maxEllipseMajor = i;
    LogEvent(QString("changed Max Ellipse to:") + QString::number(i));
}

void MainWindow::minEllipseSizevalueChanged(int i)
{
    gTrackerState.gi_minEllipseMajor = i;
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

    if (event->type() == QEvent::GraphicsSceneWheel)// QEvent::Wheel
    { //

        QGraphicsSceneWheelEvent *wheelEvent = dynamic_cast<QGraphicsSceneWheelEvent *>(event);
        handleWheelOnGraphicsScene(wheelEvent);

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
         if (gTrackerState.bDraggingTemplateCentre)
         {
            gTrackerState.bDraggingTemplateCentre = false;
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

                 gTrackerState.bStartFrameChanged = true;
                 gTrackerState.bPaused = false; //For 1 Frame and it will pause again at new frame
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
        gTrackerState.bPaused = true;
    }

    if (event->type() == QEvent::MouseButtonRelease)
    {
        gTrackerState.bStartFrameChanged = true;
        gTrackerState.bPaused = false;
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
    //Using Native Graphics Scene Scale capabilities :
    // thx to :https://stackoverflow.com/questions/19113532/qgraphicsview-zooming-in-and-out-under-mouse-position-using-mouse-wheel
  if (scrollevent->modifiers() & Qt::ControlModifier)
      qDebug() << "Pan ";


        // Do a wheel-based zoom about the cursor position

        int angle =  scrollevent->delta();
        double factor = qPow(1.0015, angle);

        qDebug() << " Angle:" << angle << " Zm:"<< factor;


        QPointF targetScenePos = scrollevent->scenePos(); // this->ui->graphicsView->mapToScene(scrollevent->pos().x(),scrollevent->pos().y());
        QPoint targetScreenPos  = scrollevent->screenPos();


        this->ui->graphicsView->scale(factor, factor);
        //Check Where Mouse Is On Scene After Transform
        QPointF nScenePos = this->ui->graphicsView->mapToScene(targetScreenPos);
        QPointF deltaViewportPos = targetScenePos-nScenePos;


        //this->ui->graphicsView->translate(100,100);
        //this->ui->graphicsView->centerOn(targetScenePos);
        //QPointF deltaViewportPos = targetViewportPos - QPointF(this->ui->graphicsView->viewport()->width() / 2.0, this->ui->graphicsView->viewport()->height() / 2.0);
        //QPointF viewportCenter = this->ui->graphicsView->mapFromScene(targetScenePos) - deltaViewportPos;

        //this->ui->graphicsView->centerOn(this->ui->graphicsView->mapToScene(viewportCenter.toPoint()));

//       this->mScene->setSceneRect(bound);
//        this->ui->graphicsView->fitInView(this->mImage->boundingRect(), Qt::KeepAspectRatio);
        //this->mImage->setPos(bound.topLeft().x() ,bound.topLeft().y());


    return;
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
    cv::Point ptMouse((int)ptImg.x(),(int)ptImg.y());

    //qDebug() << "Mouse Mv";
    if (gTrackerState.bDraggingTemplateCentre ) //bDraggingTemplateCentre
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
        ptDrag->x = std::max(gTrackerState.gFishBoundBoxSize/4, std::min(frameScene.cols - gTrackerState.gFishBoundBoxSize/4,ptMouse.x));
        ptDrag->y = std::max(gTrackerState.gFishBoundBoxSize/4, std::min(frameScene.rows - gTrackerState.gFishBoundBoxSize/4,ptMouse.y));
        gTrackerState.bROIChanged = true;
    }
    else
        ptDrag = 0;

    //Check If Around ROI Points
    for (std::vector<ltROI>::iterator it = gTrackerState.vRoi.begin(); it != gTrackerState.vRoi.end(); ++it)
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


/// \brief Outputs the 2 point click locations, for prey and mouth point positions - as set by the user
/// Additionally reports the eye vergence -
/// This function is used to validate the capture strike data from the retracked events- it is used with the script
/// validateCaptureStrikeData.R
///
void MainWindow::reportUserMeasurement(cv::Point ptMouse)
{
 ///
    float fg_EyeVergence = 0.0f;
    fishModels::reverse_iterator rt = vfishmodels.rbegin();
    if (rt != vfishmodels.rend() ) //Pick the the last fish
    {
       fishModel* pfish = rt->second;
        if (pfish)
            fg_EyeVergence = pfish->leftEyeTheta - pfish->rightEyeTheta;
    }

    this->ui->statusbar->showMessage(("Measurement point set"));
    if (userPointPair.first.x < 0)
    {
        userPointPair.first.x = ptMouse.x;
        userPointPair.first.y = ptMouse.y;
        this->SetTrackerState(7);
    }else
    {
        userPointPair.second.x = ptMouse.x;
        userPointPair.second.y = ptMouse.y;
        vMeasureLines.push_back(userPointPair);
        QString strMetro = QString("[INFO] Prey pos X:") + QString::number(userPointPair.first.x) +
                QString(" Y:") + QString::number(userPointPair.first.y) +
                QString(" Distance: ") + QString::number( cv::norm(userPointPair.second-userPointPair.first)) +
                QString(" EyeV: ") + QString::number( fg_EyeVergence);
        //Compose a comma delimeted string contaning the validation bout data - user can copy paste them to the validation script
        QString strDat = QString("[DATA] [") + QString::number(userPointPair.first.x) +
                QString(",") + QString::number(userPointPair.first.y) +
                QString(",") + QString::number(userPointPair.second.x) +
                QString(",") + QString::number(userPointPair.second.y) +
                QString(",") + QString::number(nFrame) +
                QString(",") + QString::number(gTrackerState.uiStopFrame) +
                QString(",") + QString::number(fg_EyeVergence)+
                QString("]");

        this->LogEvent(strMetro );
        this->LogEvent(strDat);
        qDebug() << strMetro;
        this->SetTrackerState(0);
        //Reset Point
        userPointPair.first.x = -10;
        /// Save an image for the records //
        this->saveScreenShot();
    }
}

///
/// \brief MainWindow::mousePressEvent
/// Two Cases - Either A Move the Fish Bounding Window Initiates
/// Or It Terminates a previous one - Saving The new template
/// \param mouseEvent
///

void MainWindow::mousePressEvent ( QGraphicsSceneMouseEvent* mouseEvent )
{

    QPointF ptSceneclick = mouseEvent->scenePos();
    QGraphicsItem* item = mScene->itemAt( ptSceneclick, this->ui->graphicsView->transform() );
    if (!item)
        return;

    // get the scene pos in the item's local coordinate space
    QPointF ptImg = item->mapFromScene(ptSceneclick);
    cv::Point ptMouse(ptImg.x(),ptImg.y());

    //Already Dragging - Terminate And Save Template
    if (gTrackerState.bDraggingTemplateCentre)
    {
        gTrackerState.bStoreThisTemplate = true;
        gTrackerState.bDraggingTemplateCentre = false;
        //bSceneMouseLButtonDown = false; //Act as Button Released With This Press
        qDebug() << "Store New Template position";
    }
    else
    {
        if (mouseEvent->buttons() == Qt::LeftButton)
        {
            bSceneMouseLButtonDown = true;
            if (gTrackerState.bAddPreyManually)
            {
                preyModel* pfood = new preyModel(cv::KeyPoint(ptMouse,1),++gTrackerState.gi_MaxFoodID  );
                pfood->blobMatchScore = 0;
                vfoodmodels.insert(IDFoodModel(pfood->ID,pfood));
            }

            if (gTrackerState.bMeasure2pDistance)
            {

                reportUserMeasurement(ptMouse);
            }
        }

        if (mouseEvent->buttons() == Qt::RightButton){
            preyModel* food = getFoodItemAtLocation(ptMouse);
            if (food) //Only delete non targeted item
                if (!food->isTargeted)
                {
                    food->inactiveFrames = gTrackerState.gcMaxFoodModelInactiveFrames+1;
                    food->isActive = false;
                    LogEvent("[info] Clicked to deactivate Food Item");
                }
        }

            //Delete Food item

        //qDebug() << "Mouse Down";

    }


    //Check If Around ROI Points
    for (std::vector<ltROI>::iterator it = gTrackerState.vRoi.begin(); it != gTrackerState.vRoi.end(); ++it)
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
    //qDebug() << "Mouse Up";
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

    cv::Point ptMouse((int)ptImg.x(),(int)ptImg.y());
    bool bFoodItemClicked = false;
    preyModel* food = getFoodItemAtLocation(ptMouse);

    if (food && mouseEvent->buttons() == Qt::LeftButton)
    {
        bFoodItemClicked = true;
        // Make Targeted
        if (!food->isTargeted)
        {
            food->isTargeted = true;
            qDebug() << "Food Targetting On  x: " << ptMouse.x << " y:" << ptMouse.y;
            LogEvent("[info] Begin Tracking Food Item");
            return;
        }else
        {
           LogEvent("[info] END Tracking Food Item");
            food->isTargeted = false;
        }
    }


    // If No Food Item Found At Location
    // Start Dragging Of Fish Template
    if (!bFoodItemClicked)
    {

        for (fishModels::iterator it=vfishmodels.begin(); it!=vfishmodels.end(); ++it)
        {

            fishModel* fish = (*it).second;
            if (fish->bodyRotBound.boundingRect().contains(ptMouse)) //Clicked On Fish Box
            {
                gTrackerState.bDraggingTemplateCentre = true;
                LogEvent("[info] Adjust Template from position ON- Start Dragging");
                this->statusBar()->showMessage(tr("Adjust fish detection template position"));
                qDebug() << "Start Dragging Fish Bound from position x: " << ptMouse.x << " y:" << ptMouse.y;

            }
        }
    }





}



MainWindow::~MainWindow()
{
    delete ui;
}

// Public method for Base code to update gui with new value of template Threshold
void MainWindow::updateTemplateThres()
{

    //Calls on_spinBoxTemplateThres_valueChanged();
    this->ui->spinBoxTemplateThres->setValue(gTrackerState.gTemplateMatchThreshold*100.0);
}

void MainWindow::on_spinBoxTemplateThres_valueChanged(int arg1)
{
 double newTMatchThresh = (double)arg1/100.0;
 gTrackerState.gTemplateMatchThreshold = newTMatchThresh;
 LogEvent(QString("[info] Changed Template Match Thres:" ) + QString::number(newTMatchThresh,'g',4) ) ;
}

void MainWindow::on_spinBoxMOGBGRatio_valueChanged(int arg1)
{
    double newBGRatio = (double)arg1/100.0;
    gTrackerState.gdMOGBGRatio = newBGRatio; //Updated Value Takes effect in processFrame and in updateBGFrame

    if (pMOG2)
    {
       pMOG2->setBackgroundRatio(gTrackerState.gdMOGBGRatio);
       LogEvent(QString("[info] Changed MOG BG Ratio: " ) + QString::number(gTrackerState.gdMOGBGRatio,'g',4) ) ;

    }


}

void MainWindow::on_actionTrack_Fish_triggered(bool checked)
{

    if (!gTrackerState.bTracking)
    {
        //iLastKnownGoodTemplateRow = 0; //Reset Row
        //iLastKnownGoodTemplateCol = 0;
        LogEvent(QString("Tracking ON"));
        this->ui->actionTrack_Fish->setChecked(true);
    }else
    {
        LogEvent(QString("Tracking OFF"));
        this->ui->actionTrack_Fish->setChecked(false);
    }

    gTrackerState.bTracking=!gTrackerState.bTracking;


}


void MainWindow::on_actionRecord_Tracks_to_File_w_triggered(bool checked)
{
    gTrackerState.bRecordToFile = checked;
    if (gTrackerState.bRecordToFile)
    {
      LogEvent(QString(">> Recording Tracks ON - New File <<"));

      resetDataRecording(outfishdatafile,"tracks");
      writeFishDataCSVHeader(outfishdatafile);
      resetDataRecording(outfooddatafile,"food");
      writeFoodDataCSVHeader(outfooddatafile);
//      closeDataFile(outfishdatafile); //
//      closeDataFile(outfooddatafile); //

//      QFileInfo fileInfFish(outfishdatafile);

//      QFileInfo fileInfFood(outfooddatafile);

//      if ( !openDataFile(fileInfFish.absoluteDir().absolutePath(),fileInfFish.completeBaseName(),outfishdatafile,"_tracks") )
//         LogEvent(QString("[Error] Opening Data Fish Tracks File"));

//      if ( !openDataFile(fileInfFood.absoluteDir().absolutePath(),fileInfFood.completeBaseName(),outfooddatafile,"_food") )
//          LogEvent(QString("[Error] Opening Data Food Tracks File"));
    }
    else
      LogEvent(QString("<< Recording Tracks OFF >>"));

}

void MainWindow::on_actionTrack_Food_triggered(bool checked)
{

    gTrackerState.bTrackFood = checked;
}

void MainWindow::on_actionQuit_triggered()
{
    gTrackerState.bExiting = true;
    LogEvent("[info] User Terminated - Bye!");
}



void MainWindow::on_actionStart_tracking_triggered()
{
    if (gTrackerState.bPaused)
        LogEvent("[info] Running");

    gTrackerState.bPaused = false;

}

void MainWindow::on_actionPaus_tracking_p_triggered()
{
    gTrackerState.bPaused = true;
    if (gTrackerState.bPaused)
    LogEvent("[info] Paused");

}
void MainWindow::on_actionPaus_tracking_p_triggered(bool checked)
{

}


void MainWindow::on_checkBoxGPU_toggled(bool checked)
{
    gTrackerState.bUseGPU = checked;
    if (gTrackerState.bUseGPU)
        LogEvent("[info] GPU use is ON");
    else
        LogEvent("[info] GPU use is OFF");
}

void MainWindow::on_checkBoxMOG_toggled(bool checked)
{
    gTrackerState.bUseBGModelling = checked;
    if (gTrackerState.bUseBGModelling)
            LogEvent("[info] Use Of MOG for BG Model is ON");
        else
            LogEvent("[info] Use Of MOG for BG Model is OFF");
}

void MainWindow::on_checkBoxNoiseFilter_toggled(bool checked)
{
    gTrackerState.bRemovePixelNoise = checked;
    if (gTrackerState.bRemovePixelNoise)
            LogEvent("[info] Pixel Noise Filtering is ON");
        else
            LogEvent("[info] Pixel Noise Filtering is OFF");
}



//void MainWindow::on_checkBoxGPU_toggled(bool checked)
//{
//    bUseGPU = checked;
//    if (bUseGPU)
//        LogEvent("[info] GPU use is OFF");
//    else
//        LogEvent("[info] GPU use is ON");
//}

//void MainWindow::on_checkBoxMOG_toggled(bool checked)
//{
//    gbUseBGModelling = checked;
//    if (gbUseBGModelling)
//            LogEvent("[info] Use Of MOG for BG Model is OFF");
//        else
//            LogEvent("[info] Use Of MOG for BG Model is ON");
//}

//void MainWindow::on_checkBoxNoiseFilter_toggled(bool checked)
//{
//    bRemovePixelNoise = checked;
//    if (bRemovePixelNoise)
//            LogEvent("[info] Pixel Noise Filtering is OFF");
//        else
//            LogEvent("[info] Pixel Noise Filtering is ON");
//}

//Set New Minimum Thrshold Scan range for Food Segmentation

//Main Loop Calls this to update The GUI SpinBox on current fitted Spine Size
void MainWindow::UpdateTailSegSizeSpinBox(float fTailSize)
{
    this->ui->doubleSpinBoxSpineSegSize->blockSignals(true);  //Don't fire change Event (avoid implicit conv to int)
    this->ui->doubleSpinBoxSpineSegSize->setValue(fTailSize);
    this->ui->doubleSpinBoxSpineSegSize->blockSignals(false);  //Don't fire change Event (avoid implicit conv to int)

}

void MainWindow::on_spinBoxFoodThresMin_valueChanged(int arg1)
{
    gTrackerState.g_SegFoodThesMin = arg1;
}

void MainWindow::on_spinBoxFoodThresMax_valueChanged(int arg1)
{
    gTrackerState.g_SegFoodThesMax = arg1;
}

void MainWindow::on_spinBoxSpineSegSize_valueChanged(int arg1)
{
    tailSizevalueChanged(arg1);
}
void MainWindow::tailSizevalueChanged(float i)
{
    qDebug() << "Tails SpinBox gave " << i;
    LogEvent(QString("Tail Segment Size changed:") + QString::number(i));
    gTrackerState.gFishTailSpineSegmentLength = i;
    //Update All fish Model's spine Length
    fishModels::iterator ft = vfishmodels.begin();
    while (ft != vfishmodels.end() ) //Render All Fish
    {
        fishModel* pfish = ft->second;
        pfish->c_spineSegL = i;

        ++ft;
    }

}


void MainWindow::on_spinBoxFoodThresMin_editingFinished()
{

}

void MainWindow::on_doubleSpinBoxSpineSegSize_valueChanged(double arg1)
{
    tailSizevalueChanged(arg1);
}

void MainWindow::on_checkBoxHistEqualizer_clicked(bool checked)
{
    gTrackerState.bUseHistEqualization = checked;
}

void MainWindow::on_checkBoxMOG_stateChanged(int arg1)
{

}

void MainWindow::on_spinBoxEyeThres_valueChanged(const QString &arg1)
{

}

void MainWindow::on_spinBoxEyeThres_valueChanged(int arg1)
{

}

void MainWindow::on_checkBoxHistEqualizer_stateChanged(int arg1)
{

}

/// \brief Change the eye Mask Width betwee eyes
void MainWindow::on_spinBoxEyeMaskW_valueChanged(int arg1)
{
    gTrackerState.iEyeMaskSepWidth = arg1;
}



void MainWindow::on_graphicsView_rubberBandChanged(const QRect &viewportRect, const QPointF &fromScenePoint, const QPointF &toScenePoint)
{

}

