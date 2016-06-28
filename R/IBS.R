IBS <- function(train,test,FittingFunction,formula){
    require(survival)
    #Kaplan Meier estimate of censoring survival function on the train set
    censor_object=with(train,Surv(as.numeric(fu),event=!dead))
    censorFit <- survfit(censor_object ~ 1)
    censorSurvPred <- censorFit$surv
    timeVector <- censorFit$time
    #Fit a model on train dataset and return prediction on test set
    survPred <- FittingFunction(train=train,test=test,surv=TRUE,time=timeVector,formula=formula)
    
    BSVector <- sapply(1:length(timeVector),function(i){
        EventPreTime <- with(test,fu < timeVector[i] & dead)
        AliveAtTime <- with(test,fu > timeVector[i])
        BS=mean((survPred[[i]]^2*EventPreTime/censorSurvPred[i]) + ((1-survPred[[i]]^2)*AliveAtTime/censorSurvPred[i]))
        BS
    })
    #Estimating integral of the Brier score function from 0 to t_max
    h=diff(timeVector)
    areaOfTrapezoids=h*0.5*(BSVector[1:(length(BSVector)-1)]+BSVector[2:(length(BSVector))])
    IntegratedBS <- sum(areaOfTrapezoids,na.rm=T)/max(timeVector)
    return(data.frame(IntegratedBS=IntegratedBS))
}