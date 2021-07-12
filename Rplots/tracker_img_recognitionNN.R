## Development Code for Willsaw NN image recognition of zfish larva
## KL 2021
## I will implement a MB random kernel NN as in my Ant Navigation paper 
## This Files Trains the weights of  2 Layer Net : Random Expansion Feature Layer (Random Kernel Trick) -> Perceptron Output neuron 
## The connections in L1 Random Sparse Net are data independent, while L2 wights are trained via +ve larva head template samples to filter patterns that look like larva
## The input images are labeled as Fish and non fish according to the folder they are placed.
# Output is saved as pgm images which are then loaded by tracker software so as to implemend the simple classifier net implemended here 

library("pixmap")


initKC <- function(m){
  # Draw random Number N of Sampled Inputs From Binomial
  n <- rbinom(1,NROW(m),N_SYN_per_KC/n_top_px)
  ## Set N random synapses as inputs to KC
  idx <- sample(1:NROW(m),n)
  
  m[idx] <- 1
  return (m)
}

KC_activation <- function(X)
{
  ## If more Than Half input units (based on Avg inputs) is active - then activate KC
  if ( sum(X) > KC_THRES )
    return (1)
  else
    return (0)
}
## Copies A Smaller Matrix To the Middle of a larger one - (Pasting Image on larger canvas)
matrix_paste <- function(src,target)
{
  x_offset <- round((ncol(target)- ncol(src))/2 )
  y_offset <- round((nrow(target)- nrow(src))/2 )
  target[ (1+y_offset):(y_offset+nrow(src)),(1+x_offset):(x_offset+ncol(src))] <-src
  return(target)
}

img_dim <- c(38,28)
n_top_px <- img_dim[2]*img_dim[1]
N_KC = n_top_px*5 ## Number of Kenyon Cells (Input layer High Dim Coding)
N_SYN_per_KC <- n_top_px/10 ## Number of pic Features each KC neuron Codes for
KC_THRES <- N_SYN_per_KC*0.55 ## Number of INput that need to be active for KC to fire/Activate

## Make Sparse Random Synaptic Weight matrix Selecting Inputs for each KC
mat_KC <<- matrix(0,nrow=n_top_px,ncol=N_KC)
sMat <<- apply(mat_KC,2,initKC)
hist(colSums(sMat),main="Number of inputs per KC")

##Layer 2 (Output Perceptron)
L2_Neurons <<- 2
W_L2 <<- matrix(0,ncol=L2_Neurons,nrow=N_KC)

## Process 2 Layer Network - Return Last Node Output produced for each input image
net_proc_images <- function(img_list,learningRate = 0.0,outNeuronIdx=1)  
{
  L2_out <- list()
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
    ## Binarize Input Image 
    X_bin <-  as.numeric(X > mean(X))
    KC_out_act = X_bin %*% sMat 
    #hist( KC_out_act[1,])
    ## Calc Layer 1 output based on Activation Threshold Funciton for Units
    KC_output <- apply(KC_out_act, 2, KC_activation)
    dim(KC_output) <-c(1,length(KC_output)) ## Make into Row Vector
    
    ## Learn Positive Samples ##
    ## Take Active L1 outputs And Set L2 Input synapses of Output Neuron to High

    ##W_L2[which(KC_output > 0.5),1] <-W_L2[which(KC_output > 0.5),1] + W_L2[which(KC_output > 0.5),1]*learningRate
    W_L2[,outNeuronIdx] <<-  W_L2[,outNeuronIdx] +  t(learningRate*KC_output)
    ##W_L2[W_L2 < 1] <<- 1 ##Cap To Limit Values - Saturation Of Synapses
    ##W_L2[,outNeuronIdx][W_L2< -1] <<- -1
    
    stopifnot(ncol(KC_output) == nrow(W_L2))
    ## Calc Output Neuron - Perceptron 
    L2_out[[fileidx]] <- list(L2_out=as.numeric(( KC_output %*%  W_L2)/N_KC),
                              KC_active = sum(KC_output),
                              KC_total = nrow(KC_output)
                              )
    message("Recognition Output for Img ",in_img," is F:",L2_out[[fileidx]]$L2_out[1]," non-F:",L2_out[[fileidx]]$L2_out[2]," Active KC:",L2_out[[fileidx]]$KC_active/N_KC )
    dim(X_bin) = dim(mat_img)
    image(X_bin)
    strRes <- "Not a Fish"
    if (L2_out[[fileidx]]$L2_out[1] > L2_out[[fileidx]]$L2_out[2])
      strRes <- "A Fish"
    title(main = paste(in_img,strRes," R:",L2_out[[fileidx]]$L2_out[1]-L2_out[[fileidx]]$L2_out[2]), font.main = 4)
  }
  
  return(L2_out)
}

## Apply Input Image ##
## list training files 
setwd("/home/kostasl/workspace/zebrafishtrack/Rplots")
sPathTrainingSamples="../img/trainset/"
sPathTrainingNonSamples="../img/nonfish/"
sPathTestingSamples="../img/test/"

img_list_train_fish =  list.files(path=sPathTrainingSamples,pattern="*pgm",full.names = T)
img_list_train_nonfish =  list.files(path=sPathTrainingNonSamples,pattern="*pgm",full.names = T)

img_list_test=  list.files(path=sPathTestingSamples,pattern="*pgm",full.names = T)

dLearningRate = 1.0 ## FALSE##TRUE
l_out <- net_proc_images(img_list_train_fish,dLearningRate,1)
l_out <- net_proc_images(img_list_train_fish,-dLearningRate/5,2)
colSums(W_L2)
l_out_n <- net_proc_images(img_list_train_nonfish,dLearningRate,2)
l_out_n <- net_proc_images(img_list_train_nonfish,-dLearningRate/5,1)
colSums(W_L2)
l_out_test <- net_proc_images(img_list_test,0.0)
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