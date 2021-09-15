## Development Code for Willsaw NN image recognition of zfish larva
## KL 2021
## I will implement a MB random kernel NN as in my Ant Navigation paper 
## This Files Trains the weights of  2 Layer Net : Random Expansion Feature Layer (Random Kernel Trick) -> Perceptron Output neuron 
## The connections in L1 Random Sparse Net are data independent, while L2 wights are trained via +ve larva head template samples to filter patterns that look like larva
## The input images are labeled as Fish and non fish according to the folder they are placed.
# Output is saved as pgm images which are then loaded by tracker software so as to implemend the simple classifier net implemended here 

library("pixmap")
library("yaml")

## TODO GIGURE OUT INIT That Works

## Initialiaze the random Weight vector of an L1 (KC) neuron  
init_random_W <- function(m, p){
  # Draw random Number N of Sampled Inputs From Binomial
  #n <- rbinom(1,length(m),p)
  ## Set N random synapses as inputs to KC
  #idx <- sample(1:length(m),n)
  m <- rnorm(NROW(m),0,sd=0.5) #runif(NROW(m))/(10*NROW(m)) ##Weak Synapses
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
  return (1/(1+exp(-activation)))
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
    mat_img <- getChannels(imgT)*255  ##x255 Convert back to uint 8-bit
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
        L_A[[l]] <- N_activation(L_X[[l-1]], inmat_W[[l-1]], inLayer_Bias[[l-1]] )
        L_X[[l]] <- N_transfer( L_A[[l]])
      }
    }
    MSQError =  sum(((mat_Y) - L_X[[N_Layers+1]] )^2)/ ncol(mat_X)

    ## Note Indexes L_X and mat/biases are off by one because LX_1 is considerened an input layer, while neural layer is l=2
    ## Back Propagation ##
    for (l in (N_Layers):1 )
    {
      ##On Output Layer
      if (l == (N_Layers))
      {
        ## Element Wise Product (hadamart Product)
        L_delta[[l]] <- ( N_transfer_D( L_X[[l+1]] )   * (L_X[[l+1]] - mat_Y)) ##*N_transfer_D(L_X[[l]])
      }else{
        #L_delta[[l]] <- L_delta[[l+1]] %*% t(mat_W[[l]])*N_transfer_D(L_X[[l]])
        L_delta[[l]] <-   t(t(L_delta[[l+1]]) %*%(inmat_W[[l+1]])*t(N_transfer_D( L_X[[l+1]] )))    ## %*% t(N_transfer_D(L_X[[l-1]]) )
      }
      
      dE <- (L_delta[[l]])%*%t(L_A[[l]])  
      dW <- learningRate*dE 
      ##Add average change over batch samples
      inmat_W[[l]] <- inmat_W[[l]] -  dW ##length(img_list)

      ## Error Non-Conform
      inLayer_Bias[[l]] <- inLayer_Bias[[l]] - rowSums( learningRate*(L_delta[[l]]) )  
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
        L_X[[1]] <- mat_X##N_transfer(X-0.5)
        
      }else
      {
        L_A[[l]] <- N_activation(L_X[[l-1]], inmat_W[[l-1]], inLayer_Bias[[l-1]] )
        L_X[[l]] <- N_transfer( L_A[[l]])
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
N_Layers <- 2

n_top_px <- img_dim[2]*img_dim[1]
N_KC = n_top_px*6 ## Number of Kenyon Cells (Input layer High Dim Coding)
N_SYN_per_KC <- 5# n_top_px/500 ## Number of pic Features each KC neuron Codes for
#KC_THRES <- N_SYN_per_KC*0.25 ## Number of INput that need to be active for KC to fire/Activate
v_Layer_N <- c(n_top_px, N_KC, 2)
v_Layer_N3 <- c(n_top_px, N_KC,100, 2)

#v_Layer_CON <- c(N_SYN_per_KC/n_top_px,)
Layer_Bias <- list() ## Number of INput that need to be active for Neuron to fire/Activate
INPUT_SPARSENESS = 0.25

mat_W <<- list() # List Of Weight Matrices

## Make Sparse Random Synaptic Weight matrix Selecting Inputs for each KC
for (k in 1:N_Layers)
{
  mat_W[[k]] <- matrix(0,ncol=v_Layer_N[k],nrow=v_Layer_N[k+1])
  ## Init Random
  if (k==1)
    mat_W[[k]] <- t(apply(mat_W[[k]],1,init_sparse_W,N_SYN_per_KC/n_top_px)) ##
  else
    mat_W[[k]] <- t(apply(mat_W[[k]],1,init_random_W,0)) ##
  
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
sPathTrainingNonSamples="../img/trainset/nonfish/"
sPathTestingSamplesFish="../img/fish/"
sPathTestingSamplesNonFish="../img/nonfish/"


img_list_train_fish =  cbind(files=list.files(path=sPathTrainingSamples,pattern="*pgm",full.names = T),F=1,NF=0) 
img_list_test_fish = cbind(files=list.files(path=sPathTestingSamplesFish,pattern="*pgm",full.names = T),F=1,NF=0)
img_list_train_nonfish =   cbind(files=list.files(path=sPathTrainingNonSamples,pattern="*pgm",full.names = T),F=0,NF=1)
img_list_test_nonfish =  cbind(files=list.files(path=sPathTestingSamplesNonFish,pattern="*pgm",full.names = T),F=0,NF=1)

## A Small Sample set for Testing the Algorithm Discrimiation between X and 1 
img_list_train_ones <-  cbind(files=list.files(path="../img/trainset/ones/",pattern="*pgm",full.names = T),F=1,NF=0)
img_list_train_X <-  cbind(files=list.files(path="../img/trainset/X/",pattern="*pgm",full.names = T),F=0,NF=1)

img_list_debug_f <-  cbind(files=list.files(path="../img/debug/fish/",pattern="*pgm",full.names = T),F=1,NF=0)
img_list_debug_nf <-  cbind(files=list.files(path="../img/debug/nonfish/",pattern="*pgm",full.names = T),F=0,NF=1)


img_list_all <- rbind.data.frame(img_list_train_fish,img_list_test_fish,img_list_train_nonfish,img_list_test_nonfish,stringsAsFactors = FALSE) #
#img_list_train <- rbind.data.frame(img_list_train_fish,img_list_test_fish,img_list_train_nonfish,img_list_train_nonfish,stringsAsFactors = FALSE)
#img_list_test=  list.files(path=sPathTestingSamples,pattern="*pgm",full.names = T) Samples ##


img_list_train <- rbind.data.frame(img_list_train_fish,img_list_test_fish,
                                   img_list_test_nonfish,img_list_train_nonfish,
                                   stringsAsFactors = FALSE)
#img_list_train <- rbind.data.frame(img_list_train_ones,img_list_train_X)
#img_list_all <- rbind.data.frame(img_list_train_fish,stringsAsFactors = FALSE)

##img_list_train <-  rbind.data.frame(img_list_debug_f,img_list_debug_nf)

batchSize = 15 # Number of Training IMages for Each Leanring Episode (which will define error graident )
Nbatches = 1500
trainingN = 5 ##Training Cycles For Each Batch



lFitError <- list()
dfitRecord <- data.frame()

vTrainingError <- vector()
## Subset INput LIst Into Batches
rIdx = 1

for (b in 1:Nbatches)
{
  ##Make Matrix -For All net inputs
  img_list_suffled <- img_list_train[sample(1:nrow(img_list_train)),]
  ##Select Subset Batch
  img_list_suffled <- head(img_list_suffled,batchSize)
  
  mat_X <- makeInputMatrix(img_list_suffled,mat_W)
  
  label_list <-cbind.data.frame(F=(img_list_suffled[,2]),NF=(img_list_suffled[,3]) )##Target output/labels
  mat_Y <- t(apply(as.matrix(label_list),2,strtoi))  

  
  for (i in 1:trainingN)
  {  
    
    dLearningRate    <- 0.0001
    
    # TRAIN On Fish 
    dnetout <- net_proc_images_batch(mat_X, mat_W, Layer_Bias, mat_Y, dLearningRate )
    mat_W   <<- dnetout$W
    mat_B   <<- dnetout$B

    #message(fileidx,". MSQERR:",MSQError,"  ", L2_out[[fileidx]]$L2_F, "-", L2_out[[fileidx]]$L2_NF, " ERR: ", L2_out[[fileidx]]$Err ," ",in_img)
    message(rIdx,". MSQERR:",dnetout$MSQError,"  ")
    
        
    vTrainingError[rIdx] = dnetout$MSQError  #plot(unlist(dnetout$out$MSERR),main=paste(i,"Mean SQ Err"))
    rIdx = rIdx +1
    
    plot(vTrainingError) #ylim=c(0,1)xlim=c(0,trainingN*Nbatches)
    
  }## Repeated Training On Batch 

} ## Different Batch Suffles  


plot(c) #ylim=c(0,1)





### Calc Final Performance 

##Select Subset Batch
img_list_suffled <- img_list_all[sample(1:nrow(img_list_all),50 ),]#rbind.data.frame(img_list_debug_f,img_list_debug_nf) img_list_train_nonfish#img_list_all[sample(1:nrow(img_list_all),50 ),] #img_list_test_nonfish #img_list_test_nonfish # 
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
dnetoutV <- net_proc_images_batch(mat_X, mat_W, Layer_Bias, mat_Y, dLearningRate )
            
message("Final MSQERR:",dnetoutV$MSQError)

  
hist(dnetout$Target - N_transfer( dnetout$output) )

hist(dnetout$W[[1]],main="After")




## FOR EXPORT TO YAML

fishNet <- list(LW1=dnetout$W[[1]],
                LW2=dnetout$W[[2]],
                LB1=dnetout$B[[1]],
                LB2=dnetout$B[[2]] 
)

#attr(fishNet$LW1, "tag") <- "!!opencv-matrix" ##Adding tags Also Change The Header to Verbatim, which does not work in OPENCV
attr(fishNet$LW2, "tag") <- "!!opencv-matrix"
attr(fishNet$LB1, "tag") <- "!!opencv-matrix"
attr(fishNet$LB2, "tag") <- "!!opencv-matrix"



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


## Process 2 Layer Network - Return Last Node Output produced for each input image
## Note : Input Layer Is Simplified - No activation function needed or Bias - Input image intentities are taken as activations
## Target_output is vector of desired output for each output Neuron these I chose to be L2_1=1 (Fish) L2_2=1 (Non Fish)
net_proc_images <- function(input_list,inmat_W,inLayer_Bias,learningRate = 0.0)  
{
  
  L_X     <- list() ## Output Of Layer k
  L_A     <- list() ## Activation of  Layer k
  L_delta <- list()  ## Delta is the "cost attributable to (the value of) that node". 
  L2_out <- list()
  
  Target_output = c(1,0)  
  fileidx <- 0
  outError = 0 #'Mean Sq Error Of File Batch'
  ## TRAIN /TEST ##
  img_list <- input_list[,1]
  label_list <-cbind.data.frame(F=(input_list[,2]),NF=(input_list[,3]) )##Target output/labels
  
  ## Matrix Of  image Input Vectors  
  mat_X = matrix(0,nrow=ncol(inmat_W[[1]]),ncol=length(img_list) )
  
  for (in_img in img_list)
  {
    fileidx= fileidx + 1
    
    imgT <- read.pnm(as.character(in_img) )
    mat_img <- getChannels(imgT)
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
    X <- sparse_binarize(as.vector(mat_img),INPUT_SPARSENESS)
    
    
    dim(X) <- c(length(X),1) ## Make into Col Vector
    
    if (length(X) > n_top_px)
    {
      warning(in_img,"Image Too large. Skipping")
      next
    }
    
    stopifnot(length(mat_X[,fileidx]) ==length(X) ) 
    mat_X[,fileidx] <- X # Save In Vector to MAtrix
    
    
    ##Input Layer Is Simplified - No activation function needed or Bias - Input image intentities are taken as activations
    ##Activation 
    L_A[[1]] <- X
    L_X[[1]] <- X##N_transfer(X-0.5)
    ## Due to R hell with numbers ecoming Factors I need to do this tricl
    mat_Y <- apply(as.matrix(label_list),2,strtoi)
    
    Target_output <-mat_Y[fileidx,] 
    ### 
    
    ## Forward Propagation ## 
    for (i in 2:(N_Layers+1) )
    {
      L_A[[i]] <- N_activation(L_X[[i-1]], inmat_W[[i-1]], inLayer_Bias[[i-1]]) ## v_Layer_Bias[i]  
      L_X[[i]] <- N_transfer(L_A[[i]])
    }
    
    outError = outError + sum( ( Target_output - L_X[[i]] ) )
    
    #hist(mat_W[[l-1]],main="before")
    ## Note Indexes L_X and mat/biases are off by one because LX_1 is considerened an input layer, while neural layer is l=2
    ## Back Propagation ##
    for (l in (N_Layers):1 )
    {
      ##On Output Layer
      if (l == (N_Layers))
      {
        ## Element Wise Product (hadamart Product)
        L_delta[[l]] <- ((N_transfer_D(L_X[[l+1]])) * (L_X[[l+1]] - Target_output)) ##*N_transfer_D(L_X[[l]])
      }else{
        #L_delta[[l]] <- L_delta[[l+1]] %*% t(mat_W[[l]])*N_transfer_D(L_X[[l]])
        L_delta[[l]] <-   t(t(L_delta[[l+1]]) %*%(inmat_W[[l+1]])*t(N_transfer_D( L_X[[l+1]] )))    ## %*% t(N_transfer_D(L_X[[l-1]]) )
        
      }
      
      dE <- (L_delta[[l]])%*%t(L_A[[l]])  
      dW <- learningRate*dE 
      ##Add average change over batch samples
      inmat_W[[l]] <- inmat_W[[l]] -  dW ##length(img_list)
      inLayer_Bias[[l]] <- inLayer_Bias[[l]] - learningRate*(L_delta[[l]])  
    }
    
    #hist(Layer_Bias[[1]])
    if (fileidx %% 10 == 0)
      hist(inmat_W[[1]], main="After")
    
    #hist( KC_out_act[1,])
    ## Calc Layer 1 output based on Activation Threshold Funciton for Units
    #KC_output <- apply(KC_out_act, 2, N_transfer)
    #dim(KC_output) <-c(1,length(KC_output)) ## Make into Row Vector
    
    ## Learn Positive Samples - Simple Perceptron Learning Rule - Over an augmented input X which is projected To High dim via sparse random matrix in L1##
    ## Take Active L1 outputs And Set L2 Input synapses of Output Neuron to High
    # Continuous output /
    #L2_Neurons_out <<- (as.numeric(( KC_output %*%  W_L2)/N_KC))
    # Weight Gradient -
    #DW <- t(KC_output)%*% (Target_output - L2_Neurons_out)/2
    #W_L2 <<-  W_L2 + learningRate*DW ## Update Weight Vector 
    
    
    ##W_L2[W_L2 < 1] <<- 1 ##Cap To Limit Values - Saturation Of Synapses
    ##W_L2[,outNeuronIdx][W_L2< -1] <<- -1
    
    #stopifnot(ncol(KC_output) == nrow(W_L2))
    ## Calc Output Neuron - Perceptron 
    ## FWD PROP CAlc MSQE - Fwd Propagate Entire Input Matrix 
    
    
    ## Forward Propagation ## 
    L1 <- N_transfer(N_activation(mat_X, inmat_W[[1]], matrix(inLayer_Bias[[1]],ncol=ncol(mat_X),nrow=length(inLayer_Bias[[1]])) ) ) ## v_Layer_Bias[i] matrix(rep(,nrow(mat_X)),nrow=nrow(mat_X) )  
    ## L2 matrix of output Col vectors correspond 
    L2 <-  N_transfer(N_activation(L1,inmat_W[[2]] , matrix(inLayer_Bias[[2]],ncol=ncol(mat_X),nrow=length(inLayer_Bias[[2]])) )  )  
    
    MSQError =  sum((t(mat_Y) - L2 )^2)/ ncol(mat_X)
    
    
    L2_out[[fileidx]] <- list(Err=0.5*sum((Target_output - L2[,fileidx])^2),
                              MSERR=MSQError,
                              L2_F=L2[1,fileidx],
                              L2_NF=L2[2,fileidx],
                              KC_active = sum(L1[fileidx,][L1[fileidx,]>0.5]),
                              KC_total = length(L1[fileidx,]),
                              # input_sparse= pxsparse,
                              file=in_img
    )
    
    
    
    message(fileidx,". MSQERR:",MSQError,"  ", L2_out[[fileidx]]$L2_F, "-", L2_out[[fileidx]]$L2_NF, " ERR: ", L2_out[[fileidx]]$Err ," ",in_img)
    
    ##message("Recognition Output for Img ",in_img," is F:",L_X[[3]][1]," non-F:",L_X[[3]][2]," Active KC:",L2_out[[fileidx]]$KC_active/N_KC )
    dim(X) = dim(mat_img)
    
    
    #image(X_bin)
    #title(main = paste(in_img,strRes," R:",L2_Neurons_out[1]-L2_Neurons_out[2]), font.main = 4)
  }
  
  lout <- list(X=mat_X,
               W=inmat_W,
               B=inLayer_Bias,
               out=data.frame( do.call(rbind,L2_out ) ))
  
  MSQError =  sum((t(mat_Y) - L2)^2)/nrow(mat_X)
  
  return(lout  )
}



#}

#       ##Test 
#   dnetout_tfish <- net_proc_images(img_list_test_fish,0.3, c(1,0))
#   #hist(unlist(dnetout_fish$L2_F),breaks=10,main="Class Fish Responses to F+ samples" )
#   ##L1 threshold set too high - patterns do not propagate
#   stopifnot(sum(unlist(dnetout_tfish$KC_active)) > 1)
#     
#   SqError_FisF <- sum((1 -unlist(dnetout_tfish$L2_F))^2)/nrow(dnetout_tfish)
#   SqError_FisNotNF <- sum((-1 -unlist(dnetout_tfish$L2_NF))^2)/nrow(dnetout_tfish)
#   
#   message("Sq. Error: F+:",SqError_FisF," F-:",SqError_FisNotNF)
#   
#   dnetout_tnotfish <- net_proc_images(img_list_test_nonfish,dLearningRate,c(0,1))
#   ##hist(unlist(dnetout_notfish$L2_NF),breaks=10,main="Class Not Fish Responses to F- samples" )
#   
#   SqError_NFisNotF <- sum((-1 -unlist(dnetout_tnotfish$L2_F))^2)/nrow(dnetout_tnotfish)
#   SqError_NFisNF <- sum((1 -unlist(dnetout_tnotfish$L2_NF))^2)/nrow(dnetout_tnotfish)
#   
#   #message("Sq. Error: F+:",SqError_NFisNotF," F-:",SqError_NFisNF)
# 
#   
#   SqClassError = sum(SqError_FisF,SqError_FisNotNF,SqError_NFisNotF,SqError_NFisNF)   
#   lFitError[[i]] = list(S_error=SqClassError,FishisFish=SqError_FisF,FishNotFish=SqError_FisNotNF,
#                         NotFishIsFish=SqError_NFisNotF,
#                         NotFishIsNotFish=SqError_NFisNF)
#   
#   dfitRecord <<- data.frame( do.call(rbind,lFitError ) )
#   plot(unlist(dfitRecord$S_error),xlim=c(1,100) ,xlab="Fit Error",ylim=c(0,2))
# 
# 
# ## Test Fish
# dnetout_testfish <- net_proc_images(img_list_test_fish,0,c(-1,1))
# ##hist(unlist(dnetout_notfish$L2_NF),breaks=10,main="Class Not Fish Responses to F- samples" )
# dnetout_testnonfish <- net_proc_images(img_list_test_nonfish,0,c(-1,1))
# 
# 
# SqError_NFisNotF <- sum((-1 -unlist(dnetout_testfish$L2_F))^2)/nrow(dnetout_testfish)
# SqError_NFisNF <- sum((1 -unlist(dnetout_testfish$L2_NF))^2)/nrow(dnetout_testfish)
# 
# 
# 
# 
# plot(unlist(dfitRecord$S_error),xlim=c(1,100) ,xlab="Fit Error",ylim=c(00,800))
# 
# plot(unlist(dfitRecord$FishisFish) )
# plot(unlist(dfitRecord$FishNotFish) )
# plot(unlist(dfitRecord$NotFishIsFish) )
# plot(unlist(dfitRecord$NotFishIsNotFish) )
# 
# 
# 
# ## Test Responses ## 
# dnetout_test <- net_proc_images(img_list_test,0.0,c(0,0))
# dnetout_test_nonfish <- net_proc_images(img_list_train_nonfish,0.0,c(-1,1))
# 
# 
# 
# colSums(W_L2)
# 
# l_out_n <- net_proc_images(img_list_train_nonfish,-dLearningRate/5,1)
# colSums(W_L2)
# 
# colSums(W_L2)
# 
# 
#  ## Save KC Matrix As Image (To Be loaded by Tracker)
#  ## Save Trained Weights As Image (To load by Tracker)
#  KC_sparse_pic <- pixmapGrey(mat_W[[1]], nrow=dim(mat_W[[1]])[1],ncol=dim(mat_W[[1]])[2])
#  plot(KC_sparse_pic)
#  write.pnm(KC_sparse_pic,file="L1_W_SparseNet.pgm")
# # 
# # ## Shift Values
#  L_W2 <- (mat_W[[2]])+1.0
#  L_W2 <- apply(L_W2,2,rev)
#  Lout_pic <- pixmapGrey(L_W2, nrow=dim(mat_W[[2]])[1],ncol=dim(mat_W[[2]])[2],cellres=1)
#  plot(Lout_pic)
#  write.pnm(Lout_pic,file="L2_W_SparseNet.pgm")
# 
#  Lout_pic <- pixmapGrey(Layer_Bias[[1]], nrow=dim(Layer_Bias[[1]])[1],ncol=dim(Layer_Bias[[1]])[2],cellres=1)
#  plot(Lout_pic)
#  write.pnm(Lout_pic,file="L1_B_SparseNet.pgm")
#  
#  
#  Lout_pic <- pixmapGrey(Layer_Bias[[2]],  nrow=dim(Layer_Bias[[2]])[1],ncol=dim(Layer_Bias[[2]])[2],cellres=1)
#  plot(Lout_pic)
#  write.pnm(Lout_pic,file="L2_B_SparseNet.pgm")
#  #
#  # 
# # 
# ## Load Matrix And Test
# imgT <- read.pnm(file="L1_W_SparseNet.pgm")
# mat_WL1 <- getChannels(imgT)
# 
# imgT <- read.pnm(file="L2_W_SparseNet.pgm")
# mat_WL2 <- getChannels(imgT)-1.0 ##Reverse Order

# 
# hist(W_L2)
# 
# cc <- getChannels(Lout_pic)
# plot(pic_outL)