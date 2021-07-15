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

