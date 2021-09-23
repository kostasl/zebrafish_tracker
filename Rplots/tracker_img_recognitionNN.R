## Development Code for Willsaw NN image recognition of zfish larva
## KL 2021
## I will implement a MB random kernel NN as in my Ant Navigation paper 
## This Files Trains the weights of  2 Layer Net : Random Expansion Feature Layer (Random Kernel Trick) -> Perceptron Output neuron 
## The connections in L1 Random Sparse Net are data independent, while L2 wights are trained via +ve larva head template samples to filter patterns that look like larva
## The input images are labeled as Fish and non fish according to the folder they are placed.
# Output is saved as pgm images which are then loaded by tracker software so as to implemend the simple classifier net implemended here 

library("pixmap")
library("yaml")



## Initialize the random Weight vector of an L1 (KC) neuron  
init_random_W3 <- function(m, p){
  m <- rnorm(NROW(m),0.0,sd=1) ## larger SD sizes make learning divergent ///#runif(NROW(m))/(10*NROW(m)) ##Weak Synapses
  return (m)
}


## Initialize the random Weight vector of an L1 (KC) neuron  
init_random_W2 <- function(m, p){
  # Draw random Number N of Sampled Inputs From Binomial
  #n <- rbinom(1,length(m),p)
  ## Set N random synapses as inputs to KC
  #idx <- sample(1:length(m),n)
  m <- rnorm(NROW(m),0,sd=1) ## larger SD sizes make learning divergent ///#runif(NROW(m))/(10*NROW(m)) ##Weak Synapses
  ## Likely Stronger Subset
  #m[idx] <- 1#/length(m) #runif(n)/(10*length(m)) #1/NROW(m)
  #print(length(m))
  return (m)
}


## Initialiaze the random Weight vector of an L1 (KC) neuron  
init_sparse_W <- function(m,p){
  # Draw random Number N of Sampled Inputs From Binomial
  n <- rbinom(1,length(m),p)
  ## Set N random synapses as inputs to KC
  idx <- sample(1:length(m),n)
  #m <- runif(NROW(m))/(1000000*NROW(m)) ##Weak Synapses
  ## Likely Stronger Subset
  m[idx] <- 1#/length(m) #runif(n)/(10*length(m)) #1/NROW(m)
  #print(length(m))
  return (m)
}


## Sigmoid//logistic Transfer Function
N_transfer <- function(activation)
{
  return (1/(1 + exp(-activation) ))
}


## Derivative Sigmoid/logistic Transfer Function
N_transfer_D <- function(activation)
{
  return (N_transfer(activation)*(1-N_transfer(activation)) )
}

## Neural Actvation for a layer with Weitgh Matrix W and Input matrix of col vectors -
#  vector of Biases B - This is to a matrix to operate overall all input activations over Inputs in order 
N_activation <- function(X,W,B)
{
  ## If more Than Half input units (based on Avg inputs) is active - then activate KC
  
  #  if ( sum(X) > KC_THRES )
  #    return (1)
  # else
  mat_B <- matrix(B,ncol=ncol(X),nrow=length(B) )
  ## Bias is just another W attached to a fixed Input 1
  return (W%*%X + mat_B)
}

## Copies A Smaller Matrix To the Middle of a larger one - (Pasting Image on larger canvas)
matrix_paste <- function(src,target)
{
  x_offset <- round((ncol(target)- ncol(src))/2 )
  y_offset <- round((nrow(target)- nrow(src))/2 )
  target[ (1+y_offset):(y_offset+nrow(src)),(1+x_offset):(x_offset+ncol(src))] <-src
  return(target)
}


sparse_binarize <- function(X,INPUT_SPARSENESS)
{
  ## Binarize Input Image / Set Sparsensess ##
  pxsparse = 1.0
  thres_bin = mean(X)
  max_iter = 100
  
  while (pxsparse > INPUT_SPARSENESS & max_iter > 0)
  {
    X_bin <-  as.numeric(X > thres_bin)
    pxsparse <- sum(as.numeric(X_bin[X_bin > 0]))/length(X_bin)
    thres_bin = thres_bin + 0.01;
    max_iter = max_iter + 1
  }
  
  
  ##message("Input Sparseness:",pxsparse)
  
  return(X_bin)
}




makeInputMatrix <- function(img_list,inmat_W)
{
  
  ## TRAIN /TEST ##
  img_list <- img_list[,1]
 
  
  
  ## Matrix Of  image Input Vectors  
  mat_X = matrix(0,nrow=ncol(inmat_W[[1]]),ncol=length(img_list) )
  fileidx <- 0
  for (in_img in img_list)
  {
    fileidx= fileidx + 1
    
    imgT <- read.pnm(as.character(in_img) )
    mat_img <- getChannels(imgT)  ## Converted 0..1 real
    X <- as.vector(mat_img)
    ### Add Fixed input 1 - To operate As Adjustable Bias for each input   
    ##X[length(X)] <- 1 # c(as.vector(mat_img),1)
    
    ##message(in_img," Input Dim:", dim(mat_img)[1],"x",dim(mat_img)[2],"=",dim(mat_img)[2]*dim(mat_img)[1])
    
    if (length(X) > n_top_px)
    {
      warning("input sample too big - skipping")
      next
    }
    
    ## If Loaded image is Smaller than INput - Then Resize (paste into) larger Matrix to fit
    if (length(X) < n_top_px)
    {
      ##message("input sample Smaller than Canvas - pasting")
      mat_t <- matrix(0,nrow=img_dim[1],ncol=img_dim[2])
      mat_img <- matrix_paste(mat_img,mat_t)
      ##message("Converted  Dim:", dim(mat_img)[1],"x",dim(mat_img)[2],"=",dim(mat_img)[2]*dim(mat_img)[1])
    }
    
    ##mypic = new("pixmapGrey", size=dim(mat_img),grey = mat_img);plot(mypic)
    ## Convert to Col Vector
    ##X <- sparse_binarize(as.vector(mat_img),INPUT_SPARSENESS)
    X <- as.vector(mat_img)
    
    dim(X) <- c(length(X),1) ## Make into Col Vector
    
    if (length(X) > n_top_px)
    {
      warning(in_img,"Image Too large. Skipping")
      next
    }
    
    stopifnot(length(mat_X[,fileidx]) ==length(X) ) 
    mat_X[,fileidx] <- X # Save In Vector to MAtrix
    
  }## Load All Input Vectors Into Matrix X
  
  return(mat_X)
} ##Make INput List



### EXPORT MAtrix to YAML for !!opencv-matrix##
matrixToYamlForOpenCV <- function(mat)
{
  ##!!opencv-matrix
  strHeader <- "!!opencv-matrix\n rows: %d\n cols: %d\n dt: f\n data:["
  ## Matrix Will Be saved Transposed Because Data serial Reading Order is by row in OpenCV while Here Serialization is By Column
  strHeader <- sprintf(strHeader,nrow(mat),ncol(mat))
  strData = ""
  lineWidth = 100#min(100,length(mat))
  
  ## COMMENTED OUT AS SUDDENLY as.yaml started adding new Lines in data automatically - The following routine does this manually
  # ##Long Strings Need to be broken by New Lines otherwise OPENCV FIle Storage FAils
  # if (length(mat) > lineWidth){
  #   fcon<- tempfile("matrixYaml")
  #   
  #   ## Need to export To TEMP file as string Concat Is too SLOW
  #   for (i in seq(1,length(mat),lineWidth) ){
  #     write(mat[i:min(i+lineWidth-1,length(mat)-1)], file = fcon,ncolumns = 100, append = TRUE, sep = ", ")
  #     ##cat(",",file = fcon, append = TRUE,fill=FALSE)
  #     #strData = paste0(strData, toString(mat[i:min(i+lineWidth,length(mat))],"\n",width=0))
  #   }
  #   ## Join  All Lines Back INto One String - With Newlines and Commas and Intendation
  #   #strData = noquote()
  #   strData = noquote(paste0(
  #                     readLines(con = fcon, n = -1L, ok = TRUE, warn = TRUE, skipNul = FALSE, encoding = "UTF-8"),
  #                            collapse=",\n         ")) ##
  #   #class(strData) <- "verbatim"
  # }else
    
  ## Need to TransposeBecause Data serial Reading Order is by row in OpenCV while Here Serialization is By Column
  strData = toString(t(mat))
  
  
  strRet=  noquote(paste0(strHeader,strData,"]")) 
  
  
  ##Character vectors that have a class of ‘verbatim’ will not be quoted in the output YAML 
  class(strRet) <- "verbatim"
  ##class(strRet) <- "verbatim"
  #write(strRet,file = fcon)  
  return(noquote(strRet) )
}




## Process 2 Layer Network - Return Last Node Output produced for each input image
## Note : Input Layer Is Simplified - No activation function needed or Bias - Input image intentities are taken as activations
## Target_output is vector of desired output for each output Neuron these I chose to be L2_1=1 (Fish) L2_2=1 (Non Fish)
net_proc_images_batch <- function(mat_X,inmat_W,inLayer_Bias,mat_Y,learningRate = 0.0)  
{
  
  L_X     <- list() ## Output Of Layer k
  L_A     <- list() ## Activation of  Layer k
  L_delta <- list()  ## Delta is the "cost attributable to (the value of) that node". 
  L2_out <- list()
  
  scaledEta = learningRate/ ncol(mat_X)
  outError = 0 #'Mean Sq Error Of File Batch'

    ##Input Layer Is Simplified - No activation function needed or Bias - Input image intentities are taken as activations
    ## Due to R hell with numbers ecoming Factors I need to do this tricl
    
    ##Target_output <- mat_Y[fileidx,] 
     
    ## Forward Propagation ##
    ##mat_X[1,] = 1 ## Add COnst Input so W_i operates as a Bias IxW 
   ##  matrix of output Col vectors correspond to output
    for (l in 1:(N_Layers+1) )
    {
      if (l == 1) ## 1st Layer is INput 
      {
        ##Activation 
        L_A[[1]] <- mat_X
        L_X[[1]] <- mat_X##N_transfer(X-0.5)
        
      }else
      {
        #L_X[[l-1]][1,] = 1 ##No Need - Removed for TEST/ Add COnst Input so W_i operates as a Bias IxW
        ## Idx k for Bias and W are advanced +1 - 
        L_A[[l]] <- N_activation(L_X[[l-1]], inmat_W[[l-1]], inLayer_Bias[[l-1]] )
        
        L_X[[l]] <- N_transfer( L_A[[l]])
        
      }
    }
    MSQError =  sum(((mat_Y) - L_X[[N_Layers+1]] )^2)/ ncol(mat_X)

    ## Note Indexes L_X and mat/biases are off by one because LX_1 is considerened an input layer, while neural layer is l=2
    
    ### Back Propagation ###
    for (l in (N_Layers):1 )
    {
      ##On Output Layer
      if (l == (N_Layers))
      {
        ## Element Wise Product (hadamart Product)   // SIGMOID * COST_Derivative ( squared error cost function wrt activations)
        L_delta[[l]] <- ( N_transfer_D( L_X[[l+1]] )   * (L_X[[l+1]] - mat_Y)) ##*N_transfer_D(L_X[[l]])
      }else{
        #L_delta[[l]] <- L_delta[[l+1]] %*% t(mat_W[[l]])*N_transfer_D(L_X[[l]])
        L_delta[[l]] <-   t(t(L_delta[[l+1]]) %*%(inmat_W[[l+1]])*t(N_transfer_D( L_X[[l+1]] )))    ## %*% t(N_transfer_D(L_X[[l-1]]) )
      }
      
      dE <- (L_delta[[l]])%*%t(L_X[[l]])  ##Instead of L_A[
      dW <- (scaledEta)*dE 
      ##Add average change over batch samples
      inmat_W[[l]] <- inmat_W[[l]] -  dW ##length(img_list)

      ## Error Non-Conform
      inLayer_Bias[[l]] <- inLayer_Bias[[l]] - rowSums( scaledEta*(L_delta[[l]]) )  
      ##if (l==1)
      ##  hist(dW)
    }
    
    #hist(Layer_Bias[[1]])
    #if (fileidx %% 10 == 0)  
    #  hist(inmat_W[[1]], main="After")
    
    
    ## Forward Propagation ## 
    # L1 <- N_transfer(N_activation(mat_X, inmat_W[[1]],inmat_B[[1]]  ) ) ## v_Layer_Bias[i] matrix(rep(,nrow(mat_X)),nrow=nrow(mat_X) )  
    # ## L2 matrix of output Col vectors correspond 
    # L2 <-  N_transfer(N_activation(L1,inmat_W[[2]] ,  inmat_B[[2]] )  )  
    # 
    # MSQError =  sum(((mat_Y) - L2 )^2)/ ncol(mat_X)
    # 
    ## Forward Propagation ##
    ##  matrix of output Col vectors correspond to output
    for (l in 1:(N_Layers+1) )
    {
      if (l == 1) ## 1st Layer is INput 
      {
        ##Activation 
        L_A[[1]] <- mat_X
        L_X[[1]] <- mat_X ##N_transfer(X-0.5)
        
      }else
      {
        L_A[[l]] <- N_activation(L_X[[l-1]], inmat_W[[l-1]], inLayer_Bias[[l-1]] )
        stopifnot(!any(is.nan( L_A[[l]]))) ##Check LAST W update did not cause NAN
        L_X[[l]] <- N_transfer( L_A[[l]])
        stopifnot(!any(is.nan(L_X[[l]]))) ##Check For NAN
      }
    }
    MSQError =  sum(((mat_Y) - L_X[[N_Layers+1]] )^2)/ ncol(mat_X)
    
    
    
    
    ##message("Recognition Output for Img ",in_img," is F:",L_X[[3]][1]," non-F:",L_X[[3]][2]," Active KC:",L2_out[[fileidx]]$KC_active/N_KC )
    #dim(X) = img_dim##dim(mat_img)
    
    
    #image(X_bin)
    #title(main = paste(in_img,strRes," R:",L2_Neurons_out[1]-L2_Neurons_out[2]), font.main = 4)
  
  
  lout <- list(X=mat_X,
               W=inmat_W,
               B=inLayer_Bias,
               output= L_X[[N_Layers+1]],
               outputActivation= L_A[[N_Layers+1]],
               Target=mat_Y,
               #out=data.frame( do.call(rbind,L2_out ) ),
               MSQError=MSQError )
  
  #MSQError =  sum(((mat_Y) - L2)^2)/nrow(mat_X)
  
  return(lout  )
}



img_dim <- c(38,28)

dInitialLearningRate    <- dLearningRate <- 0.0001
n_top_px <- img_dim[2]*img_dim[1]
N_KC = round(n_top_px*5) ## Number of Kenyon Cells (Input layer High Dim Coding)
#N_SYN_per_KC <- 50 NOT USED  #ONLY FOR L1-2  n_top_px/500 ## Number of pic Features each KC neuron Codes for
#KC_THRES <- N_SYN_per_KC*0.25 ## Number of INput that need to be active for KC to fire/Activate
#v_Layer_N2 <- c(n_top_px, N_KC, 2)
###THESE SETTINGS L=5 (500,300,100,20, 2) WORKED WITH eta=10e-4
v_Layer_N <- c(n_top_px, 500,300,100,20, 2) ##number of Units per layer (assume fully connected with Normal Dist of Strength)
N_Layers <- length(v_Layer_N)-1
#v_Layer_CON <- c(N_SYN_per_KC/n_top_px,)
Layer_Bias <- list() ## Number of INput that need to be active for Neuron to fire/Activate
INPUT_SPARSENESS = 0.25


batchSize = 132 # Number of Training IMages for Each Learning Episode (which will define error graident )
Nbatches = 1500 ## Number of random batchs (of size batchSize) to repeat training over
trainingN = 100 ## Training Cycles For Each Batch



mat_W <<- list() # List Of Weight Matrices

## Make Sparse Random Synaptic Weight matrix Selecting Inputs for each KC
for (k in 1:N_Layers)
{
  mat_W[[k]] <- matrix(0,ncol=v_Layer_N[k],nrow=v_Layer_N[k+1])
  ## Init Random
  if (k==1)
    mat_W[[k]] <-t(apply(mat_W[[k]],1,init_random_W2,0)) ## ##t(apply(mat_W[[k]],1,init_sparse_W,N_SYN_per_KC/n_top_px)) ##
  
  if (k==2)
    mat_W[[k]] <- t(apply(mat_W[[k]],1,init_random_W2,0)) ##
  
  if (k>=3)
    mat_W[[k]] <- t(apply(mat_W[[k]],1,init_random_W3,0)) ##
  
    
  Layer_Bias[[k]] <- matrix(1,ncol=1,nrow=v_Layer_N[k+1] )  #rep(1,) ## Initialiaze Neural Biases
}

#hist(mat_W[[3]])
hist(colSums(mat_W[[1]]),main="Number of inputs per KC")

##Layer 2 (Output Perceptron)
#L2_Neurons <<- 2

### RUN BATCH TRAINING ###
## Apply Input Image ##
## list training files 
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
sPathTrainingSamples="../img/trainset/fish"
sPathTrainingSamplesB="../img/trainset/fish_B"
sPathTrainingNonSamples="../img/trainset/nonfish/"
sPathTrainingNonSamplesB="../img/trainset/nonfish_B/"
sPathTestingSamplesFish="../img/fish/"
sPathTestingSamplesNonFish="../img/nonfish/"


#load(file="fishNetL3.RData")
#load(file=paste0("fishNetL",5,"-B.RData"))

img_list_train_fish =  cbind(files=list.files(path=sPathTrainingSamples,pattern="*pgm",full.names = T),F=1,NF=0)
img_list_train_fishB =  cbind(files=list.files(path=sPathTrainingSamplesB,pattern="*pgm",full.names = T),F=1,NF=0) 
img_list_test_fish = cbind(files=list.files(path=sPathTestingSamplesFish,pattern="*pgm",full.names = T),F=1,NF=0)
img_list_train_nonfish =   cbind(files=list.files(path=sPathTrainingNonSamples,pattern="*pgm",full.names = T),F=0,NF=1)
img_list_train_nonfishB =   cbind(files=list.files(path=sPathTrainingNonSamplesB,pattern="*pgm",full.names = T),F=0,NF=1)
img_list_test_nonfish =  cbind(files=list.files(path=sPathTestingSamplesNonFish,pattern="*pgm",full.names = T),F=0,NF=1)

## A Small Sample set for Testing the Algorithm Discrimiation between X and 1 
img_list_train_ones <-  cbind(files=list.files(path="../img/trainset/ones/",pattern="*pgm",full.names = T),F=1,NF=0)
img_list_train_X <-  cbind(files=list.files(path="../img/trainset/X/",pattern="*pgm",full.names = T),F=0,NF=1)

img_list_debug_f <-  cbind(files=list.files(path="../img/debug/fish/",pattern="*pgm",full.names = T),F=1,NF=0)
img_list_debug_nf <-  cbind(files=list.files(path="../img/debug/nonfish/",pattern="*pgm",full.names = T),F=0,NF=1)


img_list_all <- rbind.data.frame(img_list_train_fish,
                                 img_list_train_fishB,
                                 img_list_test_fish,
                                 img_list_train_nonfish,
                                 img_list_train_nonfishB,
                                 img_list_test_nonfish,stringsAsFactors = FALSE) #
#img_list_train <- rbind.data.frame(img_list_train_fish,img_list_test_fish,img_list_train_nonfish,img_list_train_nonfish,stringsAsFactors = FALSE)
#img_list_test=  list.files(path=sPathTestingSamples,pattern="*pgm",full.names = T) Samples ##


img_list_train <- rbind.data.frame(img_list_train_fish,
                                   img_list_train_fishB,
                                   img_list_train_nonfish,
                                   img_list_train_nonfishB, stringsAsFactors = FALSE) #
#img_list_train <- rbind.data.frame(img_list_train_ones,img_list_train_X)
#img_list_all <- rbind.data.frame(img_list_train_fish,stringsAsFactors = FALSE)
#img_list_train <-  rbind.data.frame(img_list_debug_f,img_list_debug_nf)


lFitError <- list()
dfitRecord <- data.frame()

vTrainingError <- vector()
## Subset INput LIst Into Batches
rIdx = 1

### Start Main FishNet TRAINING LOOP ###
for (b in 1:Nbatches)
{
  ##Make Matrix -For All net inputs
  ## Balance Random Batch Sample Equally between Classes  ##
  img_list_suffled <- img_list_train[sample(1:nrow(img_list_train)),] #rbind.data.frame((img_list_train[img_list_train$F == 1,])[sample(1:(batchSize/2)),] , (img_list_train[img_list_train$NF == 1,])[sample(1:(batchSize/2)),] )
  ##Select Subset Batch
  img_list_suffled <- head(img_list_suffled,batchSize)
  
  mat_X <- makeInputMatrix(img_list_suffled,mat_W)
  
  label_list <-cbind.data.frame(F=(img_list_suffled[,2]),NF=(img_list_suffled[,3]) )##Target output/labels
  mat_Y <- t(apply(as.matrix(label_list),2,strtoi))  

  message("~~~~ Learning Rate: ",dLearningRate)

  
  for (i in 1:trainingN)
  {  
    
    #dLearningRate    <- dInitialLearningRate/i
    
    
    # TRAIN On Fish 
    dnetout <- net_proc_images_batch(mat_X, mat_W, Layer_Bias, mat_Y, dLearningRate )
    mat_W   <<- dnetout$W
    mat_B   <<- dnetout$B

    #message(fileidx,". MSQERR:",MSQError,"  ", L2_out[[fileidx]]$L2_F, "-", L2_out[[fileidx]]$L2_NF, " ERR: ", L2_out[[fileidx]]$Err ," ",in_img)
    ## Evaluate Fraction Correct - Check If Fish Out Neuron Is Maximum - Compare to Target Vector 
    vclassOut <- apply(dnetout$output,2,function(x){ if(which(x== max(x) )==2) return(0) else return(1) })
    fractCorrect <- sum(dnetout$Target[1,] == vclassOut)/length(vclassOut)
    
    message(rIdx,". MSQERR:",dnetout$MSQError,", % Correct: ",fractCorrect)
  

    vTrainingError[rIdx] = dnetout$MSQError  #plot(unlist(dnetout$out$MSERR),main=paste(i,"Mean SQ Err"))
    rIdx = rIdx +1
    
    plot(vTrainingError,ylim=c(0.00001,1),type="l",log="y") #ylim=c(0,1)xlim=c(0,trainingN*Nbatches)
    
  }## Repeated Training On Batch 
  
  ##Mark End Of Batch
  points(length(vTrainingError),tail(vTrainingError,1),ylim=c(0.00001,1),type="p",log="y") #ylim=c(0,1)xlim=c(0,trainingN*Nbatches)
} ## Different Batch Suffles  

save.image(file=paste0("fishNetL",N_Layers,"-B.RData"))

plot(c) #ylim=c(0,1)





### VALIDATE Calc Final Performance ###

##Select Subset Batch
img_list_suffled <-rbind.data.frame(img_list_debug_f,img_list_debug_nf)## img_list_train_nonfish#img_list_all[sample(1:nrow(img_list_all),50 ),] #img_list_test_nonfish #img_list_test_nonfish # ##img_list_all[sample(1:nrow(img_list_all),100 ),]# 
mat_X <- makeInputMatrix(img_list_suffled,mat_W)
label_list <-cbind.data.frame(F=(img_list_suffled[,2]),NF=(img_list_suffled[,3]) )##Target output/labels
mat_Y <- t(apply(as.matrix(label_list),2,strtoi))
# 
# ## TODO : Move this in Funct - Use One MAtrix For Net -Matrix Of Biases
# mat_B <- list()
# for (k in 1:N_Layers)
#   mat_B[[k]] <- matrix(Layer_Bias[[k]],ncol=ncol(mat_X),nrow=length(Layer_Bias[[k]]) )


dLearningRate    <- 0.0
# TRAIN On Fish 
#for (i in 1:1000)
#{
  dnetoutV <- net_proc_images_batch(mat_X, mat_W, Layer_Bias, mat_Y, dLearningRate )
  mat_W <-      dnetoutV$W
  Layer_Bias <- dnetoutV$B
  vclassOut <- apply(dnetoutV$output,2,function(x){ if(which(x== max(x) )==2) return(0) else return(1) })
  fractCorrect <- sum(dnetoutV$Target[1,] == vclassOut)/length(vclassOut)
  
  message("~~~~~~~~~ Validation. MSQERR:",dnetoutV$MSQError,", % Correct: ",fractCorrect," ~~~~~~~~~~~~~")
#}

  
hist(dnetoutV$Target - N_transfer( dnetoutV$output) )




## FOR EXPORT MATRICES TO YAML SO THEY CAN BE LOADED INTO TRACKER  ##
fishNet <- list(NLayer=length(dnetout$W))
for (l in 1:fishNet$NLayer)
{
  fishNet[[paste0("LW",as.character(l))]] = dnetout$W[[l]]
  fishNet[[paste0("LB",as.character(l))]] = dnetout$B[[l]]
  attr(fishNet[[paste0("LW",as.character(l))]], "tag") <- "!!opencv-matrix" ##Adding tags Also Change The Header to Verbatim, which does not work in OPENCV
  attr(fishNet[[paste0("LB",as.character(l))]], "tag") <- "!!opencv-matrix" ##Adding tags Also Change The Header to Verbatim, which does not work in OPENCV
}



#fishNet <- list(LW1=Layer_Bias[[2]]
#)
## EXPORT TO YAML FOR OPENCV - Custom/hacked exporter routine specific to OPENCV
  filename <- "fishNet.yml"
  con <- file(filename, "w")
  message("Exporting to YAML-Wait for it , this may take a while...")
  write("%YAML:1.0",con) ##Header Is necessary For OPENCV 
  ## Indentation Needs to be 3
  str_yaml<- noquote(as.yaml(fishNet, handlers=list(matrix=matrixToYamlForOpenCV), line.sep="\n",indent=3,unicode=F))
  #  CLEAN UP String FRom Escape Chars
  ## the Var Type tag !!opencv-matrix needs to be on same line as vairable Name - For OpenCV file store - FIX
  str_yaml <- gsub("\n  !!opencv-matrix","!!opencv-matrix",str_yaml,fixed = T)
  str_yaml <- gsub("\\n","\n",str_yaml,fixed = T)## Remove Escaped NewLines
  str_yaml <- gsub("\"","",str_yaml,fixed = T)#Remove Quates

  write(noquote(str_yaml),con,append=TRUE) ##Header Is necessary For OPENCV 
  close(con)





  
  
  
  
  
  
  
  ############## 


  
  