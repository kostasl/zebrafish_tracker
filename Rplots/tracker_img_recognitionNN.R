## Development Code for Willsaw NN image recognition of zfish larva
## KL 2021
## I will implement a MB random kernel NN as in my Ant Navigation paper 
## This Files Trains the weights of  2 Layer Net : Random Expansion Feature Layer (Random Kernel Trick) -> Perceptron Output neuron 
## The connections in L1 Random Sparse Net are data independent, while L2 wights are trained via +ve larva head template samples to filter patterns that look like larva
## The input images are labeled as Fish and non fish according to the folder they are placed.
# Output is saved as pgm images which are then loaded by tracker software so as to implemend the simple classifier net implemended here 

library("pixmap")

## Initialiaze the random Weight vector of an L1 (KC) neuron  
init_random_W <- function(m,p){
  # Draw random Number N of Sampled Inputs From Binomial
  n <- rbinom(1,NROW(m),p)
  ## Set N random synapses as inputs to KC
  idx <- sample(1:NROW(m),n)
  
  m[idx] <- 1
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

N_activation <- function(X,W,B)
{
  ## If more Than Half input units (based on Avg inputs) is active - then activate KC
 
 #  if ( sum(X) > KC_THRES )
#    return (1)
 # else
    return (X%*%W -B)
}

## Copies A Smaller Matrix To the Middle of a larger one - (Pasting Image on larger canvas)
matrix_paste <- function(src,target)
{
  x_offset <- round((ncol(target)- ncol(src))/2 )
  y_offset <- round((nrow(target)- nrow(src))/2 )
  target[ (1+y_offset):(y_offset+nrow(src)),(1+x_offset):(x_offset+ncol(src))] <-src
  return(target)
}



## Process 2 Layer Network - Return Last Node Output produced for each input image
## Target_output is vector of desired output for each output Neuron these I chose to be L2_1=1 (Fish) L2_2=1 (Non Fish)
net_proc_images <- function(img_list,learningRate = 0.0,Target_output = c(1,-1))  
{

  L_X     <- list() ## Output Of Layer k
  L_A     <- list() ## Activation of  Layer k
  L_delta <- list()  ## Delta is the "cost attributable to (the value of) that node". 

  
  fileidx <- 0
  ## TRAIN /TEST ##
  for (in_img in img_list)
  {
    fileidx= fileidx + 1
    
    imgT <- read.pnm(in_img)
    mat_img <- getChannels(imgT)
    X <- as.vector(mat_img)
    message(in_img," Input Dim:", dim(mat_img)[1],"x",dim(mat_img)[2],"=",dim(mat_img)[2]*dim(mat_img)[1])
    
    if (length(X) > n_top_px)
    {
      warning("input sample too big - skipping")
      next
    }
    ## If Loaded image is Smaller than INput - Then Resize (paste into) larger Matrix to fit
    if (length(X) < n_top_px)
    {
      message("input sample Smaller than Canvas - pasting")
      mat_t <- matrix(0,nrow=img_dim[1],ncol=img_dim[2])
      mat_img <- matrix_paste(mat_img,mat_t)
      message("Converted  Dim:", dim(mat_img)[1],"x",dim(mat_img)[2],"=",dim(mat_img)[2]*dim(mat_img)[1])
      X <- as.vector(mat_img)
    }
    
    ##mypic = new("pixmapGrey", size=dim(mat_img),grey = mat_img);plot(mypic)
    dim(X) <- c(1,length(X)) ## Make into Row Vector
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
    ##Layer 1
      ##Output
      L_X[[1]] <- X_bin
      ##Activation 
      L_A[[1]] <- X
      
    message("Input Sparseness:",pxsparse)
    
    ## TRAIN Network ### 
    outError = 0
    ## Forward Propagation ## 
    for (i in 2:(N_Layers+1) )
    {
      L_A[[i]] <- N_activation(L_X[[i-1]], mat_W[[i-1]],v_Layer_Bias[i])  
      L_X[[i]] <- N_transfer(L_A[[i]])
    }
    outError = sum((Target_output - L_X[[i]])^2)
    
    ## Back Prop 
    for (l in (N_Layers+1):2 )
    {
      ##On Output Layer
      if (l == (N_Layers+1))
      {
        ## Element Wise Product (hadamart Product)
        L_delta[[l]] <- ((N_transfer_D(L_X[[l]])) * (L_X[[l]] - Target_output)) ##*N_transfer_D(L_X[[l]])
      }else{
        #L_delta[[l]] <- L_delta[[l+1]] %*% t(mat_W[[l]])*N_transfer_D(L_X[[l]])
        L_delta[[l]] <-   (L_delta[[l+1]]) %*% t(mat_W[[l]]) *N_transfer_D(L_X[[l]]) ## %*% t(N_transfer_D(L_X[[l-1]]) )
      }
      
        dE <- t(L_A[[l-1]]) %*% (L_delta[[l]])  
        dW <- learningRate*dE
        mat_W[[l-1]] = mat_W[[l-1]] -  dW
    }
      
    
    L_delta[[l]]*
    
    #hist( KC_out_act[1,])
    ## Calc Layer 1 output based on Activation Threshold Funciton for Units
    KC_output <- apply(KC_out_act, 2, N_transfer)
    dim(KC_output) <-c(1,length(KC_output)) ## Make into Row Vector
    
    ## Learn Positive Samples - Simple Perceptron Learning Rule - Over an augmented input X which is projected To High dim via sparse random matrix in L1##
    ## Take Active L1 outputs And Set L2 Input synapses of Output Neuron to High
    # Continuous output /
    L2_Neurons_out <<- (as.numeric(( KC_output %*%  W_L2)/N_KC))
    # Weight Gradient -
    DW <- t(KC_output)%*% (Target_output - L2_Neurons_out)/2
    W_L2 <<-  W_L2 + learningRate*DW ## Update Weight Vector 
    
    
    ##W_L2[W_L2 < 1] <<- 1 ##Cap To Limit Values - Saturation Of Synapses
    ##W_L2[,outNeuronIdx][W_L2< -1] <<- -1
    
    stopifnot(ncol(KC_output) == nrow(W_L2))
    ## Calc Output Neuron - Perceptron 
    
    
    L2_out[[fileidx]] <- list(L2_F=L2_Neurons_out[1],
                              L2_NF=L2_Neurons_out[2],
                              KC_active = sum(KC_output),
                              KC_total = length(KC_output),
                              input_sparse= pxsparse,
                              file=in_img
    )
    message("Recognition Output for Img ",in_img," is F:",L2_Neurons_out[1]," non-F:",L2_Neurons_out[2]," Active KC:",L2_out[[fileidx]]$KC_active/N_KC )
    dim(X_bin) = dim(mat_img)
    
    strRes <- "Not a Fish"
    if (L2_Neurons_out[1] > L2_Neurons_out[2])
      strRes <- "A Fish"
    #image(X_bin)
    #title(main = paste(in_img,strRes," R:",L2_Neurons_out[1]-L2_Neurons_out[2]), font.main = 4)
  }
  
  
  return( data.frame( do.call(rbind,L2_out ) ))
}


img_dim <- c(38,28)
N_Layers <- 2

n_top_px <- img_dim[2]*img_dim[1]
N_KC = n_top_px*5 ## Number of Kenyon Cells (Input layer High Dim Coding)
N_SYN_per_KC <- n_top_px/5 ## Number of pic Features each KC neuron Codes for
v_Layer_N <- c(n_top_px, N_KC, 2)
v_Layer_Bias <- c(N_SYN_per_KC*0.25, 0,0) ## Number of INput that need to be active for Neuron to fire/Activate
INPUT_SPARSENESS = 0.20

mat_W <- list() # List Of Weight Matrices
## Make Sparse Random Synaptic Weight matrix Selecting Inputs for each KC
for (k in 1:N_Layers)
{
  mat_W[[k]] <- matrix(0,nrow=v_Layer_N[k],ncol=v_Layer_N[k+1])
  ## Init Random
  mat_W[[k]] <- apply(mat_W[[k]],2,init_random_W, N_SYN_per_KC/n_top_px)
}

  hist(colSums(mat_W[[1]]),main="Number of inputs per KC")

  
  
##Layer 2 (Output Perceptron)
L2_Neurons <<- 2
W_L2 <<- mat_W[[2]] #matrix(0,ncol=L2_Neurons,nrow=N_KC)

## Apply Input Image ##
## list training files 
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
sPathTrainingSamples="../img/trainset/fish"
sPathTrainingNonSamples="../img/trainset/nonfish/"
sPathTestingSamplesFish="../img/fish/"
sPathTestingSamplesNonFish="../img/nonfish/"

img_list_train_fish =  list.files(path=sPathTrainingSamples,pattern="*pgm",full.names = T)
img_list_test_fish =  list.files(path=sPathTestingSamplesFish,pattern="*pgm",full.names = T)
img_list_train_nonfish =  list.files(path=sPathTrainingNonSamples,pattern="*pgm",full.names = T)
img_list_test_nonfish =  list.files(path=sPathTestingSamplesNonFish,pattern="*pgm",full.names = T)

#img_list_test=  list.files(path=sPathTestingSamples,pattern="*pgm",full.names = T)

# TRAIN On Fish Samples ##
KC_THRES <- N_SYN_per_KC*0.25 ## Number of INput that need to be active for KC to fire/Activate

dLearningRate <-1.0
lFitError <- list()
dfitRecord <- data.frame()


for (i in 1:5)
{
  
  dnetout_fish <- net_proc_images(img_list_train_fish,dLearningRate, c(1,-1))
  ##Test 
  dnetout_tfish <- dnetout_fish# net_proc_images(img_list_test_fish,dLearningRate, c(1,-1))
  #hist(unlist(dnetout_fish$L2_F),breaks=10,main="Class Fish Responses to F+ samples" )
  ##L1 threshold set too high - patterns do not propagate
  stopifnot(sum(unlist(dnetout_tfish$KC_active)) > 1)
    
  SqError_FisF <- sum((1 -unlist(dnetout_tfish$L2_F))^2)/nrow(dnetout_tfish)
  SqError_FisNotNF <- sum((-1 -unlist(dnetout_tfish$L2_NF))^2)/nrow(dnetout_tfish)
  
  message("Sq. Error: F+:",SqError_FisF," F-:",SqError_FisNotNF)
  
  dnetout_notfish <- net_proc_images(img_list_train_nonfish,dLearningRate,c(-1,1))
  dnetout_tnotfish <- dnetout_notfish #net_proc_images(img_list_test_nonfish,0,c(-1,1))
  ##hist(unlist(dnetout_notfish$L2_NF),breaks=10,main="Class Not Fish Responses to F- samples" )
  
  SqError_NFisNotF <- sum((-1 -unlist(dnetout_tnotfish$L2_F))^2)/nrow(dnetout_tnotfish)
  SqError_NFisNF <- sum((1 -unlist(dnetout_tnotfish$L2_NF))^2)/nrow(dnetout_tnotfish)
  
  #message("Sq. Error: F+:",SqError_NFisNotF," F-:",SqError_NFisNF)

  
  SqClassError = sum(SqError_FisF,SqError_FisNotNF,SqError_NFisNotF,SqError_NFisNF)   
  lFitError[[i]] = list(S_error=SqClassError,FishisFish=SqError_FisF,FishNotFish=SqError_FisNotNF,
                        NotFishIsFish=SqError_NFisNotF,
                        NotFishIsNotFish=SqError_NFisNF)
  
  dfitRecord <<- data.frame( do.call(rbind,lFitError ) )
  plot(unlist(dfitRecord$S_error),xlim=c(1,100) ,xlab="Fit Error",ylim=c(0,2))
}

## Test Fish
dnetout_testfish <- net_proc_images(img_list_test_fish,0,c(-1,1))
##hist(unlist(dnetout_notfish$L2_NF),breaks=10,main="Class Not Fish Responses to F- samples" )
dnetout_testnonfish <- net_proc_images(img_list_test_nonfish,0,c(-1,1))


SqError_NFisNotF <- sum((-1 -unlist(dnetout_testfish$L2_F))^2)/nrow(dnetout_testfish)
SqError_NFisNF <- sum((1 -unlist(dnetout_testfish$L2_NF))^2)/nrow(dnetout_testfish)




plot(unlist(dfitRecord$S_error),xlim=c(1,100) ,xlab="Fit Error",ylim=c(00,800))

plot(unlist(dfitRecord$FishisFish) )
plot(unlist(dfitRecord$FishNotFish) )
plot(unlist(dfitRecord$NotFishIsFish) )
plot(unlist(dfitRecord$NotFishIsNotFish) )



## Test Responses ## 
dnetout_test <- net_proc_images(img_list_test,0.0,c(0,0))
dnetout_test_nonfish <- net_proc_images(img_list_train_nonfish,0.0,c(-1,1))



colSums(W_L2)

l_out_n <- net_proc_images(img_list_train_nonfish,-dLearningRate/5,1)
colSums(W_L2)

colSums(W_L2)


## Save KC Matrix As Image (To Be loaded by Tracker)
## Save Trained Weights As Image (To load by Tracker)
KC_sparse_pic <- pixmapGrey(sMat, nrow=dim(sMat)[1],ncol=dim(sMat)[2],cellres=1)
plot(KC_sparse_pic)
write.pnm(KC_sparse_pic,file="KC_SparseNet.pgm")

## Shift Values
LW_shifted <- (W_L2/2.0)+0.5
Lout_pic <- pixmapGrey(W_L2, nrow=nrow(W_L2),ncol=ncol(W_L2),cellres=1)
plot(Lout_pic)
write.pnm(Lout_pic,file="outputLayer_trained.pgm")



hist(W_L2)

cc <- getChannels(Lout_pic)
plot(pic_outL)