## Development Code for Willsaw NN image recognition of zfish larva
## KL 2021
## I will implement a MB random kernel NN as in my Ant Navigation paper 


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
N_KC = n_top_px*2 ## Number of Kenyon Cells (Input layer High Dim Coding)
N_SYN_per_KC <- n_top_px/100 ## Number of pic Features each KC neuron Codes for
KC_THRES <- 3

## Make Sparse Random Synaptic Weight matrix Selecting Inputs for each KC
mat_KC = matrix(0,nrow=n_top_px,ncol=N_KC)
sMat <- apply(mat_KC,2,initKC)
hist(colSums(sMat),main="Number of inputs per KC")

##Layer 2 (Output Perceptron)
L2_Neurons <- 1
W_L2 <- matrix(0,ncol=L2_Neurons,nrow=N_KC)

## Apply Input Image ##
## list training files 
sPathTrainingSamples="../img"
img_list =  list.files(path=sPathTrainingSamples,pattern="*pgm",full.names = T)

for (in_img in img_list)
{
  imgT <- read.pnm(in_img)
  mat_img <- getChannels(imgT)
  X <- as.vector(mat_img)
  message("Input Dim:", dim(mat_img)[1],"x",dim(mat_img)[2],"=",dim(mat_img)[2]*dim(mat_img)[1])
  
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
  image(mat_img)
  dim(X) <- c(1,length(X)) ## Make into Row Vector
  
  KC_out_act = X %*% sMat 
  hist( KC_out_act[1,])
  ## Calc Layer 1 output based on Activation Threshold Funciton for Units
  KC_output <- apply(KC_out_act, 2, KC_activation)
  dim(KC_output) <-c(1,length(KC_output)) ## Make into Row Vector

  sum(KC_output)

  ## Learn Positive Samples ##
  ## Take Active L1 outputs And Set L2 Input synapses of Output Neuron to High
  W_L2[which(KC_output > 0.5),1] <- 1
  ## Calc Output Neuron - Perceptron 
  L2_out <-  (KC_output %*%  W_L2)/N_KC
  message("Recognition Output for Img ",in_img," is ",L2_out )
}