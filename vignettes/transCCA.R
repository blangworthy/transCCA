## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, results="hide", message=FALSE--------------------------------------
library(transCCA)
data("rdef6yr")

## ----est, results="hide", message=FALSE---------------------------------------
ccest <- transCCA(rdef6yr[,2:21],rdef6yr[,22:26],ndir = 1)
ccvar <- transCorVar(rdef6yr[,2:21],rdef6yr[,22:26],i = 1)

## ----ci, results="hide", message=FALSE----------------------------------------
xub <- ccest$xcoef + qnorm(0.975)*sqrt(diag(ccvar$XCoefVar)/210)
xlb <- ccest$xcoef - qnorm(0.975)*sqrt(diag(ccvar$XCoefVar)/210)

yub <- ccest$ycoef + qnorm(0.975)*sqrt(diag(ccvar$YCoefVar)/210)
ylb <- ccest$ycoef - qnorm(0.975)*sqrt(diag(ccvar$YCoefVar)/210)

corub <- min(ccest$cancor + qnorm(0.975)*sqrt(ccvar$corVar/210),1)
corlb <- max(ccest$cancor - qnorm(0.975)*sqrt(ccvar$corVar/210),0)

