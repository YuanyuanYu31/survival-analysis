library(survival)
library(ggplot2)
library(survminer)
###### read in data #########
heartfailure <- read.csv("heart_failure_clinical_records_dataset.csv",
                         header = TRUE)
head(heartfailure)


dim(heartfailure)
summary(heartfailure)
str(heartfailure)

age<-heartfailure$age
anaemia<-heartfailure$anaemia
creatinine_phosphokinase<-heartfailure$creatinine_phosphokinase
diabetes<-heartfailure$diabetes
ejection_fraction<-heartfailure$ejection_fraction
high_blood_pressure<-heartfailure$high_blood_pressure
platelets<-heartfailure$platelets
serum_creatinine<-heartfailure$serum_creatinine
serum_sodium<-heartfailure$serum_sodium
sex<-heartfailure$sex
smoking<-heartfailure$smoking
time<-heartfailure$time
DEATH_EVENT<-heartfailure$DEATH_EVENT

######### EDA ##########
# descriptive summary
summary(heartfailure)

# histogram of continuous variables
hist(time, breaks = 15, freq = F, main = "Hitogram of Time")
lines(density(time), col = "dark blue", lwd =2)

hist(age, breaks = 15, freq = F, main = "Hitogram of Age")
lines(density(age), col = "dark blue", lwd =2)

hist(creatinine_phosphokinase, breaks = 15, freq = F, 
     main = "Hitogram of Creatinine_phosphokinase")
lines(density(creatinine_phosphokinase), col = "dark blue", lwd =2)
# normal range 10-120

hist(ejection_fraction, breaks = 15, freq = F, 
     main = "Hitogram of Ejection_fraction")
lines(density(ejection_fraction), col = "dark blue", lwd =2)
# normal range 50-70

hist(platelets, breaks = 15, freq = F, 
     main = "Hitogram of Platelets")
lines(density(platelets), col = "dark blue", lwd =2)

hist(serum_creatinine, breaks = 15, freq = F, 
     main = "Hitogram of Serum_creatinine")
lines(density(serum_creatinine), col = "dark blue", lwd =2)
# normal range 0.7-1.2

hist(serum_sodium, breaks = 15, freq = F, 
     main = "Hitogram of Serum_sodium")
lines(density(serum_sodium), col = "dark blue", lwd =2)
# normal range 135-145

hist(log(serum_creatinine), breaks = 15, freq = F)

# categorical variables
table(anaemia)
#  0   1 
# 170 129 
table(diabetes)
# 0  1 
# 174 125
table(high_blood_pressure)
#  0   1 
# 194 105 
table(sex)
#  0   1 
#  105 194 
table(smoking)
# 0   1 
# 203  96

# compare age density for event and non-event
ggplot(data = heartfailure, aes(x=age, fill=factor(DEATH_EVENT))) +
  geom_density(alpha=.3)

# compare ejection_fraction density for event and non-event
ggplot(data = heartfailure, aes(x=ejection_fraction, fill=factor(DEATH_EVENT)))+
  geom_density(alpha=.3)

# compare creatinine_phosphokinase density for event and non-event
ggplot(data = heartfailure, aes(x=creatinine_phosphokinase, 
            fill=factor(DEATH_EVENT))) + geom_density(alpha=.3)

# compare serum_creatinine density for event and non-event
ggplot(data = heartfailure, aes(x=serum_creatinine, fill=factor(DEATH_EVENT)))+
  geom_density(alpha=.3)

# Kaplan-Meier
table(DEATH_EVENT)
delta<-(DEATH_EVENT==1)

fit_sex <- survfit(Surv(time,delta)~sex,  type='kaplan', conf.type='log-log', 
          data = heartfailure) 
plot(fit_sex)

fit_anaemia <- survfit(Surv(time,delta)~anaemia,  type='kaplan', 
                       conf.type='log-log', data = heartfailure) 
plot(fit_anaemia) # seems significant

fit_diabetes <- survfit(Surv(time,delta)~diabetes,  type='kaplan', 
                      conf.type='log-log', data = heartfailure) 
plot(fit_diabetes)

fit_blood_pressure <- survfit(Surv(time,delta)~high_blood_pressure,  
              type='kaplan', conf.type='log-log', data = heartfailure) 
plot(fit_blood_pressure) # seems significant

fit_smoking <- survfit(Surv(time,delta)~smoking,  type='kaplan', 
              conf.type='log-log', data = heartfailure) 
plot(fit_smoking)
# anaemia  blood pressure seems to have effect on survival 

# the general survival curve
fitKM<-survfit(Surv(time,delta)~1, type='kaplan', conf.type='log-log', 
               data = heartfailure)
fitKM
summary(fitKM)
plot(fitKM, xlab='Time (Days)', ylab='Survival Probability', 
     main="The Estimated Survival Curve")


########################
# model building #
########################
###### Cox PH model ########
fit<- coxph(Surv(time,delta)~age+anaemia+creatinine_phosphokinase+diabetes+
              ejection_fraction+high_blood_pressure+platelets+serum_creatinine+
              serum_sodium+sex+smoking, data=heartfailure)
fit
summary(fit)

# Natural history model
fit1<-coxph(Surv(time,delta)~age, data=heartfailure)
fit1#sig

fit2<-coxph(Surv(time,delta)~anaemia, data=heartfailure)
fit2#no

fit3<-coxph(Surv(time,delta)~creatinine_phosphokinase, data=heartfailure)
fit3#no

fit4<-coxph(Surv(time,delta)~diabetes, data=heartfailure)
fit4#no

fit5<-coxph(Surv(time,delta)~ejection_fraction, data=heartfailure)
fit5#sig

fit6<-coxph(Surv(time,delta)~high_blood_pressure, data=heartfailure)
fit6#sig

fit7<-coxph(Surv(time,delta)~platelets, data=heartfailure)
fit7#no

fit8<-coxph(Surv(time,delta)~serum_creatinine, data=heartfailure)
fit8#sig

fit9<-coxph(Surv(time,delta)~serum_sodium, data=heartfailure)
fit9#sig

fit10<-coxph(Surv(time,delta)~sex, data=heartfailure)
fit10#no

fit11<-coxph(Surv(time,delta)~smoking, data=heartfailure)
fit11#no


#backward elimination

fita<- coxph(Surv(time,delta)~age+anaemia+creatinine_phosphokinase+diabetes+
               ejection_fraction+high_blood_pressure+platelets+serum_creatinine+
               serum_sodium+sex+smoking, data=heartfailure)
fita
summary(fita)
#delete platelets with the largest p value 0.6806 


fitb<- coxph(Surv(time,delta)~age+anaemia+creatinine_phosphokinase+diabetes+
               ejection_fraction+high_blood_pressure+serum_creatinine+
               serum_sodium+sex+smoking, data=heartfailure)
fitb
summary(fitb)
#delete smoking with the largest p value 0.6455    


fitc<- coxph(Surv(time,delta)~age+anaemia+creatinine_phosphokinase+diabetes+
               ejection_fraction+high_blood_pressure+serum_creatinine+
               serum_sodium+sex, data=heartfailure)
fitc
summary(fitc)
#delete diabetes with the biggest p value 0.5681    


fitd<- coxph(Surv(time,delta)~age+anaemia+creatinine_phosphokinase+
               ejection_fraction+high_blood_pressure+serum_creatinine+
               serum_sodium+sex, data=heartfailure)
fitd
summary(fitd)
# delete sex with the biggest p value  0.4109

fite<- coxph(Surv(time,delta)~age+anaemia+creatinine_phosphokinase+
               ejection_fraction+high_blood_pressure+serum_creatinine+
               serum_sodium, data=heartfailure)
fite
summary(fite)
# the p value of serum_sodium is 0.0505

###### functional form for age ##########

fit_age <- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                   ejection_fraction+high_blood_pressure+serum_creatinine
                 , method='breslow',data=heartfailure)
# martingale res. vs. age
plot(age,resid(fit_age))
# resid(fitd_age) gives martingale residuals
lines(lowess(age,resid(fit_age),f=.5))
abline(lm(resid(fit_age)~age),lty=3, col="red")
# not linear, seems there are two different slopes.

# cut off at age == 65
old<-(age>=65)
fite_old<- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                   ejection_fraction+high_blood_pressure+serum_creatinine+
                   serum_sodium+age+old, method='breslow',data=heartfailure)
summary(fite_old)
# when both age and old are included in the model, old is not sig.

##interaction between age and old
interac_age_old<-age*old

fite_old_inter<- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                  ejection_fraction+high_blood_pressure+serum_creatinine+
                 serum_sodium+age+old+interac_age_old, method='breslow',
                  data=heartfailure)
summary(fite_old_inter) # all not sig

# age^2
age2<-age^2
fite_age2<- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                    ejection_fraction+high_blood_pressure+serum_creatinine+
                    serum_sodium+age+age2, method='breslow',data=heartfailure)

summary(fite_age2) #both not sig

# only age2 included
fite_age2_only<- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                      ejection_fraction+high_blood_pressure+serum_creatinine+
                      serum_sodium+age2, method='breslow',data=heartfailure)

summary(fite_age2_only) #both not sig

# martingale res. vs. age2
fit_age <- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                   ejection_fraction+high_blood_pressure+serum_creatinine
                 , method='breslow',data=heartfailure)
# martingale res. vs. age
plot(age2,resid(fit_age))
# resid(fitd_age) gives martingale residuals
lines(lowess(age2,resid(fit_age),f=.5))
abline(lm(resid(fite)~age2),lty=3, col="red")
# there are also two slopes, age2 can not solve the problem


#log age

fite_log_age<- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                       ejection_fraction+high_blood_pressure+serum_creatinine+
                       serum_sodium, method='breslow',data=heartfailure)
# martingale res. vs. age
plot(log(age),resid(fite_log_age))
# resid(fitd_age) gives martingale residuals
lines(lowess(log(age),resid(fite_log_age),f=.5))
abline(lm(resid(fite_log_age)~log(age)),lty=3, col="red")

#sqrt age
fite_sqrt_age<- coxph(Surv(time,delta)~anaemia+creatinine_phosphokinase+
                        ejection_fraction+high_blood_pressure+serum_creatinine+
                        serum_sodium, method='breslow',data=heartfailure)
# martingale res. vs. age
plot(sqrt(age),resid(fite_sqrt_age))
# resid(fitd_age) gives martingale residuals
lines(lowess(sqrt(age),resid(fite_sqrt_age),f=.5))
abline(lm(resid(fite_sqrt_age)~sqrt(age)),lty=3, col="red")

## cube age
fite_cub_age<- coxph(Surv(time,delta)~age+I(age^2)+I(age^3)+anaemia+
                creatinine_phosphokinase+ejection_fraction+high_blood_pressure+
            serum_creatinine+serum_sodium, method='breslow',data=heartfailure)
summary(fite_cub_age)
# all are not significant


# factor age into 4 groups now
agegroup<-ifelse(age<65,1,2) # create new categorical variable agegroup, cut off age at 65 years old
agegroup<-ifelse(age>=65 & age<75,2,agegroup)
agegroup<-ifelse(age>=75 & age<85,3,agegroup)
agegroup<-ifelse(age>=85,4,agegroup)
table(agegroup)
agegroup <- factor(agegroup)

agegroup2 <- ifelse(agegroup==2,1,0)
agegroup3 <- ifelse(agegroup==3,1,0)
agegroup4 <- ifelse(agegroup==4,1,0)

###### functional form for creatinine_phosphokinase ########
fite_enz<- coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                   ejection_fraction+high_blood_pressure+serum_creatinine+
                   serum_sodium, method='breslow',data=heartfailure)

# martingale res. vs. creatinine_phosphokinase
plot(creatinine_phosphokinase,resid(fite_enz))
# resid(fitd_enz) gives martingale residuals
lines(lowess(creatinine_phosphokinase,resid(fite_enz),f=.5))
abline(lm(resid(fite_enz)~creatinine_phosphokinase),lty=3, col="red")
# keep it as linear form


###### functional form for ejection_fraction #########
fite_eject<- coxph(Surv(time,delta)~factor(agegroup)+anaemia+
              creatinine_phosphokinase+high_blood_pressure+serum_creatinine+
                     serum_sodium, method='breslow',data=heartfailure)

# martingale res. vs. ejection_fraction
plot(ejection_fraction,resid(fite_eject))
# resid(fitd_eject) gives martingale residuals
lines(lowess(ejection_fraction,resid(fite_eject),f=.5))
abline(lm(resid(fite_eject)~ejection_fraction),lty=3, col="red")
# there is turning point at 35.

##### factor eject_fraction (this is what we decided to use)#########
eject_group<-ifelse(ejection_fraction<35,1,2)
eject_group<-ifelse(ejection_fraction>50,3,eject_group)
eject_group <- factor(eject_group)

eject_group1 <- ifelse(eject_group==1,1,0)
eject_group2 <- ifelse(eject_group==2,1,0)
eject_group3 <- ifelse(eject_group==3,1,0)

fite_eject_factor<- coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                creatinine_phosphokinase +high_blood_pressure+serum_creatinine+
                    serum_sodium+factor(eject_group),
                   method='breslow',data=heartfailure)
summary(fite_eject_factor)


###### functional form for serum_creatinine #############
fite_serum_creatinine<- coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                creatinine_phosphokinase+high_blood_pressure+
                        serum_sodium+factor(eject_group),
                               method='breslow',data=heartfailure)

# martingale res. vs. serum_creatinine
plot(serum_creatinine,resid(fite_serum_creatinine))
# resid(fitd_serum_creatinine) gives martingale residuals
lines(lowess(serum_creatinine,resid(fite_serum_creatinine),f=.5))
abline(lm(resid(fite_serum_creatinine)~serum_creatinine),lty=3, col="red")

#log transformation
fite_serum_creatinine_trans<- coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                creatinine_phosphokinase+high_blood_pressure+serum_creatinine+
                  log(serum_creatinine)+serum_sodium+factor(eject_group),
                               method='breslow',data=heartfailure)

summary(fite_serum_creatinine_trans) #log is sig

# omit log(serum_creatinine)
fite_serum_creatinine_log<-coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                          creatinine_phosphokinase +high_blood_pressure
                                 +serum_sodium+factor(eject_group)
                                 , method='breslow',data=heartfailure)

# martingale res. vs. log serum_creatinine
plot(log(serum_creatinine),resid(fite_serum_creatinine_log))
# resid(fitd_serum_creatinine_log) gives martingale residuals
lines(lowess(log(serum_creatinine),resid(fite_serum_creatinine_log),f=.5))
abline(lm(resid(fite_serum_creatinine_log)~log(serum_creatinine)),lty=3,col="red")

# the log serum_creatinine is almost linear in the model,keep log term.


#factor creatinine

serum_creatinine_g<-ifelse(serum_creatinine>=1.3, 1, 0)

fite_serum_creatinine_factor2<-coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                              creatinine_phosphokinase
                              +high_blood_pressure+factor(serum_creatinine_g)
                                     +serum_sodium+factor(eject_group)
                                     , method='breslow',data=heartfailure)
summary(fite_serum_creatinine_factor2)

# all are sig except serum_sodium, hence we drop the serum_sodium

fite_final_test<-coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                         creatinine_phosphokinase
                       +high_blood_pressure+factor(serum_creatinine_g)
                       +factor(eject_group)
                       , method='breslow',data=heartfailure)
summary(fite_final_test)

temp<-cox.zph(fite_final_test)
temp
ggcoxzph(temp)

# From the temp, we found the p value of eject_group is smaller than 0.05, 
# which means this violates the cox PH assumption, we fixed it by strata
fite_final_strata<-coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                           creatinine_phosphokinase
                         +high_blood_pressure+factor(serum_creatinine_g)
                         +strata(eject_group)
                         , method='breslow',data=heartfailure)
summary(fite_final_strata)

temp<-cox.zph(fite_final_strata)
temp
ggcoxzph(temp)
# now all variables have p values larger than 0.05


# we hope to see which is better between log(serum_creatinine) and
# factor(serum_creatinine) before strata eject_group
fite_final_test_log<-coxph(Surv(time,delta)~factor(agegroup)+anaemia+
                             creatinine_phosphokinase
                           +high_blood_pressure+log(serum_creatinine)
                           +factor(eject_group)
                           , method='breslow',data=heartfailure)

summary(fite_final_test_log)
temp<-cox.zph(fite_final_test_log)
temp
ggcoxzph(temp)

# using AIC to compare which one is better
AIC(fite_final_test, fite_final_test_log)

#                     df      AIC
# fite_final_test      9 953.6133
# fite_final_test_log  9 947.0568

# Hence, we decided to use the model with log(serum_creatinine) 

##################
#Final model#
##################

fit_final<-coxph(Surv(time,delta)~agegroup2+agegroup3+agegroup4+anaemia+
                   creatinine_phosphokinase
                 +high_blood_pressure+log(serum_creatinine)
                 +strata(eject_group)
                 , method='breslow',data=heartfailure)
summary(fit_final)

temp<-cox.zph(fit_final)
temp
ggcoxzph(temp)


########################
# Detecting outliers #
########################
fit_final<-coxph(Surv(time,delta)~agegroup2+agegroup3+agegroup4+anaemia+
                   creatinine_phosphokinase
                 +high_blood_pressure+log(serum_creatinine)
                 +strata(eject_group)
                 , method='breslow',data=heartfailure)
summary(fit_final)

# calculate the risk scores
rs<-rbind(coef(fit_final)) %*% rbind(agegroup==2,agegroup==3, agegroup==4, 
                                     anaemia,creatinine_phosphokinase,
                                     high_blood_pressure,log(serum_creatinine))
rs<-as.vector(rs)
rs
plot(rs,resid(fit_final,'dev'),xlab='Risk Score',ylab='Deviance Residual')
abline(0,0)
identify(rs,resid(fit_final,'dev'))
#[1] 2 and 229


###################################
# Identifying influence points #
###################################

fit_final<-coxph(Surv(time,delta)~agegroup2+agegroup3+agegroup4+anaemia+
                   creatinine_phosphokinase
                 +high_blood_pressure+log(serum_creatinine)
                 +strata(eject_group)
                 , method='breslow',data=heartfailure)
summary(fit_final)
dim(sresid)
# 299 5

plot(factor(agegroup),sresid[,1],xlab='Age<65=1, Age>=65=2',ylab='Residual')
identify(factor(agegroup),sresid[,1])
#[1] 229
plot(anaemia,sresid[,2],xlab='anaemia',ylab='Residual')
identify(anaemia,sresid[,2])
#[1] 49

plot(creatinine_phosphokinase,sresid[,3],xlab='creatinine_phosphokinase',
                           ylab='Residual')
identify(creatinine_phosphokinase,sresid[,3])
# [1] 118

plot(high_blood_pressure,sresid[,4],xlab='high_blood_pressure',ylab='Residual')
identify(high_blood_pressure,sresid[,4])
# [1] 229

plot(log(serum_creatinine),sresid[,5],xlab='log(serum_creatinine)',
                          ylab='Residual')
identify(log(serum_creatinine),sresid[,5])
# [1] 2

mtext(outer=T,line=-1,cex=1.5,"Figure 2: influential cases")

# the identified observations
heartfailure[2,] #creatinine_phosphokinase-highest
heartfailure[229,]# ejection_fraction-low, time 207 days.
heartfailure[49,] 
heartfailure[118,]

#################################################################
# Estimated survival function for a given set of covariate values
#################################################################

## fit the final cox PH model
fit_final<-coxph(Surv(time,delta)~agegroup2+agegroup3+agegroup4+anaemia+
                   creatinine_phosphokinase
                 +high_blood_pressure+log(serum_creatinine)
                 +strata(eject_group1,eject_group2,eject_group3)
                 , method='breslow',data=heartfailure)
summary(fit_final)


## Provide covariate values

zvalues<-data.frame(agegroup2=1,agegroup3=0,agegroup4=0,anaemia=0,
                  creatinine_phosphokinase=146,high_blood_pressure=0,
            serum_creatinine=1.3,eject_group1=1,eject_group2=0,eject_group3=0)


##  Calculate estimated survival function
fit_predict<-survfit(fit_final, newdata = zvalues)
plot(fit_predict, xlab = "Time(Days)", ylab = "Survival Probability",
  main="Survival Curve for Hypothetical Patient(Ejection Group 1)",col = "red")

summary(fit_predict)
