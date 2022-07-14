library(haven)
library(ivreg)

dat1 <- read_dta('colonial_origins/maketable1/maketable1.dta')
dat4 <- read_dta('colonial_origins/maketable4/maketable4.dta')
dat5 <- read_dta('colonial_origins/maketable5/maketable5.dta')
dat6 <- read_dta('colonial_origins/maketable6/maketable6.dta')

dat4_clean <- dat4[!is.na(dat4$baseco),]
#lm(logpgp95 ~ avexpr , data = dat)

est_simples <- ivreg(logpgp95 ~ avexpr | logem4 , data = dat4_clean)
# dat <- as.matrix(dat1[,c(4,3,9)])
# dat <- dat[!rowSums(!is.finite(dat)),]
# twoSLS(1,2,3,dat=dat)
est_simples
sqrt(vcov(est_simples)[2,2])

# latitude is correlated with settler mortality and therefore harms
est_lat <- ivreg(logpgp95 ~ avexpr + lat_abst| logem4 +lat_abst, data = dat4_clean)
est_lat
summary(est_lat)
sqrt(vcov(est_lat)[2,2])

summary(ivreg(formula = logpgp95 ~ avexpr | logem4 + lat_abst, data = dat4_clean))

# ethnic fragmentation 
dat6_clean <- dat6[!is.na(dat6$baseco),]
est_fragm <- ivreg(logpgp95 ~ avexpr + avelf| logem4 +avelf, data = dat6_clean)
est_fragm
summary(est_fragm)
sqrt(vcov(est_fragm)[2,2])

# european ancestry is correlated with settler mortality
est_large <- ivreg(logpgp95 ~ avexpr + avelf + edes1975| logem4 +avelf + edes1975, data = dat6_clean)
est_large
summary(est_large)
sqrt(vcov(est_large)[2,2])

#
est_large_ins <- ivreg(logpgp95 ~ avexpr + avelf| logem4 + edes1975 + avelf, data = dat6_clean)
est_large_ins
summary(est_large_ins)
sqrt(vcov(est_large_ins)[2,2])



######## mortlaity + european as instrumnets
est_n <-  ivreg(logpgp95 ~ avexpr| logem4 + edes1975, data = dat6_clean)
est_n 
sqrt(vcov(est_n)[2,2])

######## mortlaity + european as conditioning variable
est_n2 <-  ivreg(logpgp95 ~ avexpr+ edes1975| logem4 + edes1975, data = dat6_clean)
est_n2 
sqrt(vcov(est_n2)[2,2])


