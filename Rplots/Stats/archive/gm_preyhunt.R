###Generative Model OF Pre Hunt ##
##Author KL Feb 2018
##statistical Repetitions of Exp
nTrials <- 10000

##Number of initial prey Added to Dish / Sample From Binomial With   mean of 15
nPrey <- round(rnorm(mean=15,n=10,sd=5))

pDetection -> 0.15 ##Probability of detecting Prey
pCapture   -> 0.1 ##Probability Of Successfull Hunt Given detection
pObserve    -> 0.60 ##Prob Of Seeing An Event or A rotifer given not all events are recorded


xp <-  round(rnorm(mean=15,n=10000,sd=7))
X11()
hist(xp)