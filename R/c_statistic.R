#' Concordance Statistic
#' @export
#Calculates concordance statistics as defined by Harrel et al.
c_statistic <- function(train,test,FittingFunction,covariates,time,surv=TRUE){
    prediction <- FittingFunction(train=train,test=test,covariates=covariates,time=time,surv=TRUE) %>%
        unlist() %>% as.numeric()
    epsilon=10^-6
    fuDiff <- with(test,outer(fu,fu,'-'))
    predDiff <- outer(prediction,prediction,'-')
    eventSum <- with(test,outer(event,event,'+'))
    eventDiff <- with(test,outer(event,event,'-'))
    permissableMat <- matrix(TRUE,nrow = nrow(test),ncol=nrow(test))
    permissableMat[which(eventDiff*fuDiff > epsilon | (abs(fuDiff) < epsilon & (eventSum < epsilon)) | eventSum < epsilon)] <- FALSE

    concordanceMat <- matrix(0,nrow = nrow(test),ncol=nrow(test))

    concordanceMat[which(permissableMat & abs(fuDiff) > epsilon & predDiff*fuDiff > epsilon)] <- 1

    concordanceMat[which(permissableMat & abs(fuDiff) > epsilon & abs(predDiff) < epsilon)] <- 0.5

    concordanceMat[which(permissableMat & abs(fuDiff) < epsilon & abs(eventSum-2) < epsilon)] <- 0.5

    concordanceMat[which(permissableMat & abs(fuDiff) < epsilon & abs(eventSum-2) < epsilon & abs(predDiff) < epsilon)] <- 1

    concordanceMat[which(permissableMat & abs(fuDiff) < epsilon & abs(eventSum-1) < epsilon)] <- 0.5

    concordanceMat[which(permissableMat & abs(fuDiff) < epsilon & abs(eventSum-1) < epsilon & eventDiff*predDiff < 0)] <- 1

    sum(concordanceMat[upper.tri(concordanceMat)])/sum(permissableMat[upper.tri(permissableMat)])
}
