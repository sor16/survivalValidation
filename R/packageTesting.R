library(survival) #survival analysis
library(eha) #used for data 
data(oldmort) #create the data
#create surv data set
mort <-Surv(oldmort$exit-oldmort$enter,event=oldmort$event)
covariates="sex"
fullCalibrationList <- crossValidation(data=oldmort,fu=oldmort$exit-oldmort$exit,event=oldmort$event,k=2,nrRuns=10,ValidationFunction = calibration,FittingFunction=evalCox,by=0.1,time=1,covariates=covariates)
fullCalibrationData <-  lapply(fullCalibrationList,bind_rows) %>% bind_rows()
fullData <- fullCalibrationData %>% group_by(deciles) %>% summarise(proportion=mean(proportion,na.rm=T),lower=mean(lower,na.rm=T),upper=mean(upper,na.rm=T))%>%
    filter(!is.na(proportion))

linearModel=with(fullData,lm(proportion~deciles))
textDat= data.frame(x=rep(min(fullData$deciles),2),y=c(max(fullData$upper)/2 + 0.05 ,max(fullData$upper)/2 - 0.05),
                    text=paste(c("A","B"),format(linearModel$coefficients,digits = 2),sep=" = "))
ggplot(data=fullData,aes(deciles,proportion)) + geom_point(na.rm=T) + geom_smooth(method = "lm",se=F,na.rm=T) + geom_text(data=textDat,aes(x,y,label=text))+
    geom_errorbar(aes(ymin = lower,ymax = upper,width=0.01,color="red"),na.rm=T) + theme_bw() +
    theme(legend.position="none")

library(survivalValidation)
library(survival)
library(epitools)
data("wcgs")
CovariateData <- wcgs %>% select(-id,-chd69,-time169)
Covariates <- names(CovariateData)
formula <- paste("surv_object",paste(Covariates,collapse=" + "),sep=" ~ ") 
surv_object <- with(wcgs,Surv(time169,event=chd69))
Model <- coxph(as.formula(formula),data=CovariateData)
#newFormula <- paste("surv_object",paste(Covariates[summary(Model)$coefficients[,5]<0.05],collapse=" + "),sep=" ~ ")

#######-Calibration-#########
fullCalibrationList <- crossValidation(data=CovariateData,fu=wcgs$time169,event=wcgs$chd69,k=2,nrRuns=10,ValidationFunction = calibration,FittingFunction=evalCox,by=0.1,time=3000,covariates=Covariates)
fullCalibrationData <-  lapply(fullCalibrationList,bind_rows) %>% bind_rows()
fullData <- fullCalibrationData %>% group_by(deciles) %>% summarise(proportion=mean(proportion,na.rm=T),lower=mean(lower,na.rm=T),upper=mean(upper,na.rm=T))%>%
    filter(!is.na(proportion))

linearModel=with(fullData,lm(proportion~deciles))
textDat= data.frame(x=rep(min(fullData$deciles),2),y=c(max(fullData$upper)/2 + 0.05 ,max(fullData$upper)/2 - 0.05),
                    text=paste(c("A","B"),format(linearModel$coefficients,digits = 2),sep=" = "))
ggplot(data=fullData,aes(deciles,proportion)) + geom_point(na.rm=T) + geom_smooth(method = "lm",se=F,na.rm=T) + geom_text(data=textDat,aes(x,y,label=text))+
    geom_errorbar(aes(ymin = lower,ymax = upper,width=0.01,color="red"),na.rm=T) + theme_bw() +
    theme(legend.position="none")
#######-IBS-#########
t <- Sys.time()
IntegratedBS <- crossValidation(data=CovariateData,fu=wcgs$time169,event=wcgs$chd69,k=10,ValidationFunction = IBS,FittingFunction=evalCox,nrRuns = 10,covariates=Covariates)
Sys.time() - t
IBSData <- sapply(IntegratedBS,unlist) %>% rowMeans() %>% data.frame(IntegratedBrierScore = (.))

ModelType <- factor("Cox")
ggplot(data = IBSData,aes(ModelType,IntegratedBrierScore)) + geom_boxplot()

