---
title: "Cross-sectional associations between 24-hour time-use composition, grey matter
  volume and cognitive function in healthy older adults"
author: "Maddison Mellow, Dot Dumuid and Ty Stanford"
date: "2023-05-24"
output: 
  html_document:
    keep_md: true
---



## Background

Increasing physical activity (PA) is an effective strategy to slow reductions in cortical volume and maintain cognitive function in older adulthood. However, PA does not exist in isolation, but coexists with sleep and sedentary behaviour to make up the 24-hour day. We investigated how the balance of all three behaviours (24-hour time-use composition) is associated with grey matter volume in healthy older adults, and whether grey matter volume influences the relationship between 24-hour time-use composition and cognitive function. 

The following analysis pipelines were used for our study, "Cross-sectional associations between 24-hour time-use composition, grey matter volume and cognitive function in healthy older adults" (https://www.medrxiv.org/content/10.1101/2023.05.15.23289982v1). 

## Data overview

The code presented here was replicated for each brain volume outcome (total grey matter, temporal lobe, hippocampus, lateral ventricle, frontal lobe). Only the code used for the total grey matter volume outcome is presented here, for simplicity. 


The following variables were included in analyses:

**Grey matter volume**

- `Neurotix.GM.Vol` = total grey matter volume (uncorrected)
- `GM.corrected` = total grey matter volume (corrected)

**Cognitive outcome measure**

- `longtermmem` = long-term memory
- `execfunc` = executive functions
- `procspeed` = processing speed

**Time-use variables**

- `all_days_sleeptime` = total time spent in sleep per day (averaged over recording period)
- `all_days_sedtime` = total time spent in sedentary behaviour per day (averaged over recording period)
- `all_days_mvtime` = total time spent in moderate-vigorous physical activity per day (averaged over recording period)
- `all_days_lighttime` = total time spent in light physical activity per day (averaged over recording period)


**Covariates/variables used in correction**

- `visit_age` (age in years at baseline visit)
- `sex` = 2 levels (male or female)
- `site_ch` = site, 2 levels (Adelaide or Newcastle)
- `edu_yrs` = total years education (primary + secondary + tertiary)
- `dist_correct` = use of distortion correction during imaging (on/off)
- `Neurotix.TIV` = total intracranial volume (ml)

## Code 

#### 1. Load required packages

```r
library("dplyr") 
library("tidyr") 
library("readr")
library("ggplot2")
library("tibble")
library("lubridate") 
library("compositions") 
library("foreach")
library("knitr") 
library("skimr")
library("GGally") 
library("car")
library("rstatix")
library("GGally")
library("robCompositions")
library("boot")
library("mice")
library("performance")
```

#### 2. Load and review dataset


```r
load('data/FinalData/act_dat.RData', verbose = TRUE)
skim(act_dat)
```

#### 3. Refine dataset 
Make sure factor variables are classified correctly


```r
str(act_dat)
act_dat$site_ch <- factor(act_dat$site_ch)
act_dat$Sex <- factor(act_dat$Sex)
act_dat$Dist_correct <- factor(act_dat$Dist_correct)
act_dat$valid_file <- factor(act_dat$valid_file)
str(act_dat)
```


Look at missing data across entire dataframe


```r
md.pattern(act_dat, rotate.names = TRUE)
```

Remove participants who are missing time-use data (in this case, all participants who were missing time-use data were missing sleep variables, so participants were removed based on this rule)


```r
nrow(act_dat)
act_dat <-
  act_dat %>%
  dplyr::filter(!is.na(all_days_sleeptime))
nrow(act_dat)
```

#### 4. Correct volume data for total intracranial volume, site, and distortion correction. 

For the purpose of this script, only the correction for total grey matter volume is displayed here.


```r
fit.GM <- lm(Neurotix.GM.Vol ~ Neurotix.TIV + site_ch + Dist_correct, data = act_dat)
summary(fit.GM)

act_dat[['GM.corrected']] <-  mean(act_dat$Neurotix.GM.Vol) +  act_dat$Neurotix.GM.Vol - predict(fit.GM, newdata=act_dat, type='response')
```

Check dispersion of raw VS. corrected values

```r
ggplot(act_dat, aes(x=Neurotix.GM.Vol, y=GM.corrected, col=site_ch)) + geom_point(size = 1.5)
```

#### 5. Clean and refine dataset

Check that time-use data add up to 1440 minutes (24-hrs)

```r
#rowSums = sums values in specified columns for ALL rows- to check if they add up to roughly 1440 
rowSums(act_dat[, 
                c("all_days_sleeptime", "all_days_sedtime", "all_days_lighttime", "all_days_mvtime")
])
```

Look for any invalid files (coded as NA or 0, as 1=valid)

```r
table(act_dat$valid_file, useNA = 'ifany')
act_dat <- act_dat[act_dat$valid_file %in% 1, ] #remove participants that don't have valid accelerometry files
table(act_dat$valid_file, useNA = 'ifany') #check that they were removed 
```

Create new dataframes: cols_outcome and cols_pred (predictor)

```r
cols_covar <- c("ID", "Age", "Sex", "edu_yrs")
cols_outcome <- c("GM.corrected", "procspeed", "longtermmem", "execfunc")
cols_pred <- c("all_days_sleeptime", "all_days_sedtime", "all_days_lighttime", "all_days_mvtime")
```

Combine cols_outcome and cols_pred in to cols_want dataframe

```r
cols_want <- c(cols_covar, cols_outcome, cols_pred)
```

Create new dataframe 'actdat' which contains all rows from act_dat but only cols_want columns

```r
cols_want[!(cols_want %in% colnames(act_dat))]
actdat <- act_dat[, cols_want]
as_tibble(actdat)
```

Check that time use == 1440 minutes

```r
hist(rowSums(actdat[, cols_pred])) #create histogram of all rows but only predictor variable columns
actdat[rowSums(actdat[, cols_pred]) > 1500, ] #show which rows add up to >1500 mins
actdat[rowSums(actdat[, cols_pred]) < 1400, ] #show which rows add up to <1400 mins
```

Remove participants with more than 1500 minutes of time-use data
***note that this was an arbitrary but pragmatic decision to remove observations with potentially misleading time-use recordings made by researchers for this study.

```r
nrow(actdat)
actdat <- actdat[!(rowSums(actdat[, cols_pred]) > 1500), ]
nrow(actdat)
rowSums(actdat[, cols_pred])
```

Apply closure function (`compositions` package), which rescales time-use data so that all participants have 1440 mins time use

```r
actdat[, cols_pred] <- clo(actdat[, cols_pred], total = 1440)
rowSums(actdat[, cols_pred])
```

#### 6. Correlations between independent time-use behaviours, ROI volumes, cognitive outcomes and continuous covariates

Create `cols_to_use` dataframe containing variables of interest

```r
cols_to_use <- c("Age", "edu_yrs" , "all_days_lighttime",
          "all_days_mvtime", "all_days_sedtime",
          "all_days_sleeptime", "GM.corrected", "Temp.corrected",
          "Hippo.corrected", "Vent.corrected", "Frontal.corrected",
          "procspeed", "execfunc", "longtermmem")
```

The following code was used to investigate Pearson correlations between variables outside of the composition (i.e., all correlations except between time-use behaviours)

```r
nc <- length(cols_to_use)
actdat[, cols_to_use]
cor <- matrix("", nrow = nc, ncol = nc, dimnames = list(cols_to_use, cols_to_use))
for (i in 1:(nc - 1)) { # i <- 1
  for (j in (i + 1):nc) { # j <- 3
    
    print(i)
    print(j)
    
    var_i <- actdat[[cols_to_use[i]]]
    var_j <- actdat[[cols_to_use[j]]]
      cor_ij <- cor.test(var_i, var_j, subset = !(is.na(var_i) | is.na(var_j)))
      cor[i, j] <- 
        sprintf("%1.3f (%1.3f, %1.3f)", cor_ij$estimate, cor_ij$conf.int[1], cor_ij$conf.int[2])
    
  }
}
cor[1:3, 1:3]
```

The following code was then used to investigate correlations between time-use behaviours, using an alternative method (symmetric balanced isometric log-ratios)

Create dataframe containing time-use variables only

```r
corrs <- subset(actdat, select=c('all_days_lighttime', 'all_days_mvtime', 'all_days_sedtime', 'all_days_sleeptime'))
head(corrs)
```

Apply `corCoDa` function (within `robCompositions` package) to time-use variables

```r
corrs <- corCoDa(corrs)
corrs
corrs.v = as.vector(corrs[upper.tri(corrs)]) 
```

Transfer values in the correlation matrix as a vector (goes down the columns)

```r
ndat = actdat %>% dplyr::select('all_days_lighttime', 'all_days_mvtime', 'all_days_sedtime', 'all_days_sleeptime')
ndat=as.data.frame(ndat)
```

Create non-parametric 95% confidence intervals using non-parametric bootstrapping (1000 samples)

```r
compCorr = function(d, i)
{
  lpa= d[i,1]
  mvpa=d[i,2]
  sb=d[i,3]
  sleep=d[i,4]
  
  comp=cbind.data.frame(lpa, mvpa, sb, sleep)
  matrix = corCoDa(comp)
  mat1 = as.vector(matrix[upper.tri(matrix)]) #gives the values in the correlation matrix as a vector (goes down the cols)
  mat1
}
compCorr(ndat)

(xd.boot <- boot(ndat, compCorr, 1000))

limper<-matrix(0,2,6) 
for (var in 1:6)
{
  xd.boot.ci<-boot.ci(xd.boot, index=var,conf=0.95,type = "perc")
  limper[,var]<-xd.boot.ci$per[4:5]}
```

View compositional correlations (first columns) and 95% bootstrap CI (second and third column)

```r
cbind.data.frame(corrs.v, as.data.frame(t(limper)))

round(cbind.data.frame(corrs.v, as.data.frame(t(limper))),2)
```

#### 7. Compositional data analysis setup

Note that the order in which the predictor variables have been arranged will impact the interpretation of the isometric log ratios. In this instance, we are setting up the ilrs so that the first ilr represents sleep:remaining behaviours, the second ilr represents sedentary:active behaviours (light PA + MVPA), and the third ilr represents light PA:MVPA. All three ilrs are entered in to linear regression models to capture the entire 24-hour time-use composition. 

Arrange data

```r
cols_pred <- c("all_days_sleeptime","all_days_sedtime","all_days_lighttime", "all_days_mvtime")
```

Create sequential binary partition matrix

```r
sbp4 = matrix(c( 1, -1, -1,-1, 
                 0, 1, -1, -1,
                 0, 0, 1, -1),
              ncol=4, byrow=TRUE)


psi4 = gsi.buildilrBase(t(sbp4)) 
```

Compute ilrs of all rows, but only for columns in cols_pred

```r
ilrs = ilr(acomp(actdat[, (cols_pred)]), V=psi4) 
head(ilrs)
```
Rename ilrs to 'ilr1, ilr2, and ilr3'

```r
colnames(ilrs) <- paste0("ilr", 1:3)
head(ilrs)
nrow(ilrs)
nrow(actdat)
```

Create new dataframe '`ilr_dat`' that contains the data to be used in regression models

```r
ilr_dat <- cbind(actdat[, c(cols_outcome, cols_covar)], ilrs) 
head(ilr_dat)
```

Add polynomial time-use composition as a variable in `ilr_dat` to assess whether quadratic terms improve model fit

```r
c_ <- with(ilr_dat, cbind(ilr1, ilr2, ilr3))
ilr_dat$poly_c <- poly(c_, degree = 2, raw = TRUE)
```

#### 8. Linear regression models investigating primary research aim (associations between 24-hour time-use composition and ROI volumes) 

Model 1 = total GM volume as the outcome, covariates as predictors (age, sex, education)
Model 2 = total GM volume as the outcome, covariates AND time-use composition as predictors
Model 3 = total GM volume as the outcome, covariates AND time-use composition (polynomial terms) as predictors


```r
ilr_mod_1 <- lm(GM.corrected ~ Age + Sex + edu_yrs, dat = ilr_dat)
ilr_mod_2 <- lm(GM.corrected ~  Age + Sex + edu_yrs + cbind(ilr1, ilr2, ilr3), dat = ilr_dat)
ilr_mod_3 <- lm(GM.corrected ~  Age + Sex + edu_yrs + poly_c, dat = ilr_dat)
```

Check assumptions of model using `performance` package

```r
check_model(ilr_mod_1)
check_model(ilr_mod_2)
```

View outputs of linear regression models (and then overall F test using `car::Anova`)

```r
summary(ilr_mod_1)
summary(ilr_mod_2)
summary(ilr_mod_3)
car::Anova(ilr_mod_1)
car::Anova(ilr_mod_2) 
car::Anova(ilr_mod_3)
```

Use `anova` test to determine if addition of polynomial terms improve model prediction

```r
anova(ilr_mod_1, ilr_mod_2, ilr_mod_3)
```

In this case, model 1 was used as the final model.

Adjust p-values of final F test using false discovery rate adjustment (Benjamini-Hochberg FDR adjustment) 

```r
pvals <- c(2.499e-08, 0.8613, 0.2876)
p.adjust(pvals, method = "BH", n=length(pvals))
```

#### 9. Linear regression models investigating secondary research aim (whether cognitive function is associated with interaction between ROI volume and time-use composition)

Model 1 = long-term memory as the outcome; covariates, main-effects of time-use composition and total GM volume, and interaction between time-use composition and total GM volume as predictors.


```r
ilr_mod_1 <- lm(longtermmem ~ Age + Sex + edu_yrs + cbind(ilr1, ilr2, ilr3) * GM.corrected, dat=ilr_dat)
```

View outputs of linear regression models (and then overall F test using `car::Anova`)

```r
summary(ilr_mod_1)
car::Anova(ilr_mod_1)
```

Adjust p-values of final F test using false discovery rate adjustment (Benjamini-Hochberg FDR adjustment)

```r
pvals <- c(0.74672, 0.08606, 0.24769, 0.03582, 0.37437, 0.01063)
p.adjust(pvals, method = "BH", n=length(pvals))
```

NOTE: This process is then repeated with executive function and processing speed as the outcome variables (and each of the three models then repeated for other ROI volume outcomes)

#### 10. Predictive modelling

In these analyses, we found a significant time-use composition * total GM volume interaction for executive function outcomes. The following code was used to generate the predictive modelling plots, which demonstrate the predicted effects of reallocating time towards or away from each time-use behaviour on executive function z-score. These relationships were plotted separately across participants with higher and lower total GM volume (above and below the sample mean).

NOTE: To run this code, a `GM.grp` variable was created whereby participants were classified into "upper" (above mean GM volume) or "lower" (below GM volume) groups. 

Make one-for-remaining reallocated compositions


```r
act=actdat[,cols_pred]
act=acomp(act)
(m=clo(mean(act),total=1))

# Make one-for-remaining reallocated compositions
#change sleep
#add 15 min Sleep
r=15/1440/m[1]
s=r* m[1]/(1-m[1])
sl.15i=acomp(cbind(m[1]*(1+r),m[2]*(1-s),m[3]*(1-s),m[4]*(1-s)))
clo(sl.15i, total=1440) #sleep has 15 min more than the mean
#this is the mean
clo(m, total=1440)

#add 30 min Sleep
r=30/1440/m[1]
s=r* m[1]/(1-m[1])
sl.30i=acomp(cbind(m[1]*(1+r),m[2]*(1-s),m[3]*(1-s),m[4]*(1-s)))
clo(sl.30i, total=1440)

#take 15 min Sleep
r=-15/1440/m[1]
s=r* m[1]/(1-m[1])
sl.15d=acomp(cbind(m[1]*(1+r),m[2]*(1-s),m[3]*(1-s),m[4]*(1-s)))
clo(sl.15d, total=1440)

#take 30 min Sleep
r=-30/1440/m[1]
s=r* m[1]/(1-m[1])
sl.30d=acomp(cbind(m[1]*(1+r),m[2]*(1-s),m[3]*(1-s),m[4]*(1-s)))
clo(sl.30d, total=1440)

#change sb
#add 15 min Sb
r=15/1440/m[2]
s=r* m[2]/(1-m[2])
sb.15i=acomp(cbind(m[1]*(1-s),m[2]*(1+r),m[3]*(1-s),m[4]*(1-s)))
clo(sb.15i, total=1440)

#add 30 min Sb
r=30/1440/m[2]
s=r* m[2]/(1-m[2])
sb.30i=acomp(cbind(m[1]*(1-s),m[2]*(1+r),m[3]*(1-s),m[4]*(1-s)))
clo(sb.30i, total=1440)

#take 15 min Sb
r=-15/1440/m[2]
s=r* m[2]/(1-m[2])
sb.15d=acomp(cbind(m[1]*(1-s),m[2]*(1+r),m[3]*(1-s),m[4]*(1-s)))
clo(sb.15d, total=1440)

#take 30 min Sb
r=-30/1440/m[2]
s=r* m[2]/(1-m[2])
sb.30d=acomp(cbind(m[1]*(1-s),m[2]*(1+r),m[3]*(1-s),m[4]*(1-s)))
clo(sb.30d, total=1440)

#change lpa
#add 15 min lpa
r=15/1440/m[3]
s=r* m[3]/(1-m[3])
lpa.15i=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1+r),m[4]*(1-s)))
clo(lpa.15i, total=1440)

#add 30 min lpa
r=30/1440/m[3]
s=r* m[3]/(1-m[3])
lpa.30i=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1+r),m[4]*(1-s)))
clo(lpa.30i, total=1440)

#take 15 min lpa
r=-15/1440/m[3]
s=r* m[3]/(1-m[3])
lpa.15d=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1+r),m[4]*(1-s)))
clo(lpa.15d, total=1440)

#take 30 min lpa
r=-30/1440/m[3]
s=r* m[3]/(1-m[3])
lpa.30d=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1+r),m[4]*(1-s)))
clo(lpa.30d, total=1440)

#change mvpa
#add 15 min mvpa
r=15/1440/m[4]
s=r* m[4]/(1-m[4])
mvpa.15i=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1-s),m[4]*(1+r)))
clo(mvpa.15i, total=1440)

#add 30 min mvpa
r=30/1440/m[4]
s=r* m[4]/(1-m[4])
mvpa.30i=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1-s),m[4]*(1+r)))
clo(mvpa.30i, total=1440)

#take 15 min mvpa
r=-15/1440/m[4]
s=r* m[4]/(1-m[4])
mvpa.15d=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1-s),m[4]*(1+r)))
clo(mvpa.15d, total=1440)

#take 30 min mvpa
r=-30/1440/m[4]
s=r* m[4]/(1-m[4])
mvpa.30d=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1-s),m[4]*(1+r)))
clo(mvpa.30d, total=1440)

#gather comps
comps <- rbind(sl.30d, sl.15d, m, sl.15i, sl.30i,
            sb.30d, sb.15d, m, sb.15i, sb.30i,
            lpa.30d, lpa.15d, m, lpa.15i, lpa.30i,
            mvpa.30d, mvpa.15d, m, mvpa.15i, mvpa.30i)

new.ilrs=ilr(comps, V=psi4)#make ilrs with the one-for-remaining reallocated compositions
```


Make linear regression model again using numeric variable type for binary variables

```r
ilr_dat$sexN=as.numeric(ilr_dat$Sex)
table(ilr_dat$sexN)
ilr_mod_1 <- lm(execfunc ~ Age + sexN + edu_yrs +
                  GM.grp*cbind(ilr1, ilr2, ilr3), dat = ilr_dat)
```

Check that the interaction remains significant when total GM volume is classified as a categorical (binary) variable

```r
car::Anova(ilr_mod_1)
mean(actdat$edu_yrs, na.rm=TRUE)
```

Make new dataset to predict for

```r
newdata=expand.grid(new.ilrs[,1], Age=mean(actdat$Age), sexN=1.5, edu_yrs=16.55733, GM.grp=c('Upper', 'Lower'))
newdata$ilr1=new.ilrs[,1]#add ilr1
newdata$ilr2=new.ilrs[,2]#add ilr2
newdata$ilr3=new.ilrs[,3]#add ilr3
head(newdata)
nrow(newdata)
```

Use new data to predict executive function

```r
pred.EF=predict(ilr_mod_1, newdata, interval="confidence")
Activity=clo(comps,total=1440)
plotdata=cbind(newdata, pred.EF, sleep=Activity[,1], sb=Activity[,2], lpa=Activity[,3], mvpa=Activity[,4])
head(plotdata)
plotdata$Activity <- rep(c(rep("Sleep",5), rep("SB",5), rep("LPA",5), rep("MVPA",5)),2) #createS  a variable for plotting 
#it indicates the "one" variable in the "one-to-remaining" reallocations, i.e., which variable we are considering the -30, -15 etc mins for. 
plotdata$Activity=factor(plotdata$Activity, levels=c("Sleep","SB","LPA","MVPA"))
plotdata$Reallocation=c(rep(seq(-30,30,15),8))#create a variable to show what the reallocations are 
head(plotdata)
plotdata$Predicted_EF = plotdata$fit
plotdata$GM.grp=factor(plotdata$GM.grp, levels = c("Lower", "Upper"))#make a factor to help plotting
str(plotdata)
```

Create plot using ggplot

```r
ggplot(plotdata, aes(x=Reallocation, y=Predicted_EF, group=GM.grp))+
  geom_hline(yintercept = 0, col = "black", alpha= 0.25) +
  geom_point(aes(colour=GM.grp))+
  geom_line(aes(colour=GM.grp))+
  scale_color_manual("Total grey matter volume", values=c("orange", "purple"))+
  scale_fill_manual("Total grey matter volume", values=c("orange", "purple"))+
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill=GM.grp),alpha=0.1) +
  facet_wrap(~Activity, scales="free_x")+
  labs(x = "Time reallocation (minutes)", y = "Predicted difference in executive function z-score") +
  scale_x_continuous(breaks = seq(-30,30,15))+
  theme_bw() + 
  theme(legend.position="bottom")
```






