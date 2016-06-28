ROC  <- function(train,test,FittingFunction,time,formula){
    require(dplyr)
    #Check if patient has had an event at a user specified time
    statusAtTime <- with(test,fu < time & dead)
    #Fit a model on train dataset and return prediction on test set
    prediction <- FittingFunction(train=train,test=test,surv=FALSE,time=time,formula=formula) %>% unlist()
    threshold <- seq(0,1,length.out=100)
    ROCData <- lapply(threshold,function(i){
        eventPredict = prediction >= i
        TPR=sum(eventPredict & statusAtTime)/sum(statusAtTime)
        FPR=sum(eventPredict & !statusAtTime)/sum(!statusAtTime)
        data.frame(TPR=TPR,FPR=FPR)
    }) %>% rbind_all() %>% arrange(-row_number())
    ROCData$threshold <- threshold
    #trapisunalgun
    AreaOfTrapezoids=with(ROCData,c(0,(FPR[2:(length(TPR))]-FPR[1:(length(TPR)-1)])*0.5*(TPR[1:(length(TPR)-1)]+TPR[2:(length(TPR))])))
    ROCData$AUC=sum(AreaOfTrapezoids)
    ROCData$cumArea=cumsum(AreaOfTrapezoids)
    return(ROCData)
}