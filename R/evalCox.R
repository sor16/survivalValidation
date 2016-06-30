#' Predicting Cox Proportional Hazard model
#' @export
#Calculates Cox model and returns list of prediciton for a user specified time
evalCox <- function(train,test=train,formula,time,surv=TRUE){
    require(survival)
    if(all(!c("fu","event") %in% names(train))){
        stop("follow up time in data has to be called fu and the event has to be called event.")
    }
    surv_object <- with(train,Surv(as.numeric(fu),event))
    coxModel <- coxph(as.formula(formula),data=train)
    #Use model to predict outcome in testset
    baselineCumHazard <- basehaz(coxModel)
    LinearComb <- predict(coxModel,newdata=test)
    if(surv){
        predictionFUN <- function(i) exp(-i*exp(LinearComb))
    }else{
        predictionFUN <- function(i) 1-exp(-i*exp(LinearComb))
    }
    prediction <- lapply(baselineCumHazard$hazard[time==baselineCumHazard$time],FUN=predictionFUN)
    return(prediction)
}
