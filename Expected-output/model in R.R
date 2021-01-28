library(aod)
library(Rcpp)
library(MASS)

data.cnv <- read.table("sample.wise.input", header=T)

#model 
model.burden <- glm(RESPONSE ~ PREDICTOR + Sex + PC1 +PC2 +PC3, data=data.cnv, family=binomial)
sum.mod <- summary(model.burden)
sum.tab <-sum.mod$coefficients
p.mod<-with(model.burden, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
or.ci.mod<-exp(cbind(OR = coef(model.burden), confint(model.burden)))
main.results<-cbind(p.mod, or.ci.mod)
