## ----setup, include = FALSE---------------------------------------------------
#file.edit(normalizePath("~/.Renviron"))
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=TRUE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----pine---------------------------------------------------------------------
library(bootPLS)
library(plsRglm)
data(pine, package = "plsRglm")
Xpine<-pine[,1:10]
ypine<-log(pine[,11])

## ----pinedisplay, cache=FALSE, eval=LOCAL-------------------------------------
pairs(pine)

## ----pinecv, cache=FALSE, eval=LOCAL------------------------------------------
bbb <- plsRglm::cv.plsR(log(x11)~.,data=pine,nt=6,K=nrow(pine),NK=1,verbose=FALSE)
plsRglm::cvtable(summary(bbb))

## ----pinecvrep, cache=FALSE, eval=LOCAL---------------------------------------
bbb2 <- plsRglm::cv.plsR(log(x11)~.,data=pine,nt=6,K=6,NK=100,verbose=FALSE)

## ----pinecvtable, cache=FALSE, eval=LOCAL-------------------------------------
plsRglm::cvtable(summary(bbb2))

## ----pineplotcvtableQ2, cache=FALSE, eval=LOCAL-------------------------------
plot(plsRglm::cvtable(summary(bbb2)),type="CVQ2")

## ----pineplotcvtablePress, cache=FALSE, eval=LOCAL----------------------------
plot(plsRglm::cvtable(summary(bbb2)),type="CVPress")

## ----nbcompbootplsR, cache=FALSE, eval=LOCAL----------------------------------
set.seed(4619)
nbcomp.bootplsR(Y=ypine,X=Xpine,R =500)

## ---- resbootrep, cache=FALSE, eval=LOCAL-------------------------------------
set.seed(4619)
res_boot_rep <- replicate(20,nbcomp.bootplsR(Y=ypine,X=Xpine,R =500,verbose =FALSE,parallel = "multicore",ncpus = 2))

## ----barplotresbootrep, cache=FALSE, eval=LOCAL-------------------------------
barplot(table(res_boot_rep))

## ----nbcompbootspls, cache=FALSE, eval=LOCAL----------------------------------
nbcomp.bootspls(x=Xpine,y=ypine,eta=.5)

## ----spls, cache=FALSE, eval=LOCAL--------------------------------------------
nbcomp.bootspls.para(x=Xpine,y=ypine,eta=.5)
nbcomp.bootspls.para(x=Xpine,y=ypine,eta=c(.2,.5))

## ----pinereload, cache=FALSE, eval=LOCAL--------------------------------------
library(bootPLS)
library(plsRglm)
data(pine, package = "plsRglm")
Xpine<-pine[,1:10]
ypine<-log(pine[,11])
datasetpine <- cbind(ypine,Xpine)

## ----pineadaptyx, cache=FALSE, eval=LOCAL-------------------------------------
coefs.plsR.adapt.ncomp(datasetpine,sample(1:nrow(datasetpine)))

## ----pineadaptyxreplicate, cache=FALSE, eval=LOCAL----------------------------
replicate(20,coefs.plsR.adapt.ncomp(datasetpine,sample(1:nrow(datasetpine))))

## ----pineadaptmc, cache=FALSE, eval=LOCAL-------------------------------------
coefs.plsR.adapt.ncomp(datasetpine,sample(1:nrow(datasetpine)),ncpus=2,parallel="multicore")

## ----pineadaptmcreplicate, cache=FALSE, eval=LOCAL----------------------------
replicate(20,coefs.plsR.adapt.ncomp(datasetpine,sample(1:nrow(datasetpine)),ncpus=2,parallel="multicore"))

## ----azeload, cache=FALSE, eval=LOCAL-----------------------------------------
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
dataset <- cbind(y=yaze_compl,Xaze_compl)

## ----azefit, cache=FALSE, eval=LOCAL------------------------------------------
modplsglm <- plsRglm(y~.,data=dataset,10,modele="pls-glm-family",family="binomial")

## ----azebootcomp, cache=FALSE, eval=LOCAL-------------------------------------
set.seed(4619)
aze_compl.bootYT <- suppressWarnings(nbcomp.bootplsRglm(modplsglm))

## ----azeplotbootcomp, cache=FALSE, eval=LOCAL---------------------------------
plsRglm::boxplots.bootpls(aze_compl.bootYT)

## ----azebootcompCI, cache=FALSE, eval=LOCAL-----------------------------------
plsRglm::confints.bootpls(aze_compl.bootYT)

## ----azebootcompplotCI, cache=FALSE, eval=LOCAL-------------------------------
suppressWarnings(plsRglm::plots.confints.bootpls(plsRglm::confints.bootpls(aze_compl.bootYT)))

## ----prostateload, cache=FALSE, eval=LOCAL------------------------------------
set.seed(4619)
data(prostate, package="spls")
nbcomp.bootsgpls((prostate$x)[,1:30], prostate$y, R=250, eta=0.2, typeBCa = FALSE)

## ----prostatespls, cache=FALSE, eval=LOCAL------------------------------------
nbcomp.bootsgpls.para((prostate$x)[,1:30], prostate$y, R=250, eta=c(.2,.5), maxnt=10, typeBCa = FALSE)

## ----clean, cache=FALSE, eval=LOCAL-------------------------------------------
rm(list=c("Xpine","ypine","bbb","bbb2","aze_compl","Xaze_compl","yaze_compl","aze_compl.bootYT","dataset","modplsglm","pine","datasetpine","prostate","res_boot_rep"))

