## Background

Increasing physical activity (PA) is an effective strategy to slow
reductions in cortical volume and maintain cognitive function in older
adulthood. However, PA does not exist in isolation, but coexists with
sleep and sedentary behaviour to make up the 24-hour day. We
investigated two main questions/aims: 1) is 24-hour time-use composition
associated with grey matter volume in healthy older adults?, and 2) does
grey matter volume moderate the relationship between 24-hour time-use
composition and cognitive function? Herein, we will refer to these as
aims 1 and 2.

The following analysis pipelines were used for our study,
“Cross-sectional associations between 24-hour time-use composition, grey
matter volume and cognitive function in healthy older adults”, published
in the International Journal of Behavioral Nutrition and Physical
activity
(<https://link.springer.com/article/10.1186/s12966-023-01557-4>).

## Data overview

The code presented here was replicated for each ROI volume in aim 1
(total grey matter, frontal lobe, temporal lobe, hippocampi, lateral
ventricles) and each cognitive outcome in aim 2 (long-term memory,
executive function, processing speed). For simplicity, here we will
present the code used for total grey matter and executive function only.

The following variables were included in analyses:

***Brain volume measures***

-   `GM.vol` = total grey matter volume (ml), uncorrected
-   `GM.corrected` = total grey matter volume (ml), corrected for total
    intracranial volume, site, and distortion correction

***Cognitive outcome measure***

-   `execfunc` = measure of executive function

***Time-use variables***

-   `all_days_sleeptime` = total time spent in sleep per day (averaged
    over recording period)
-   `all_days_sedtime` = total time spend in sedentary behaviour per day
    (averaged over recording period)
-   `all_days_mvtime` = total time spent in moderate-vigorous physical
    activity per day (averaged over recording period)
-   `all_days_lighttime` = total time spend in light physical activity
    per day (averaged over recording period)

***Covariates***

-   `Age` = age in years
-   `Sex` = sex assigned at birth (male/female only in this study)
-   `edu_yrs` = years of education

***Variables used to correct volume values***

-   `TIV` = total intracranial volume
-   `site` = study site (Adelaide or Newcastle)
-   `DC` = distortion correction feature used during MRI scan (yes/no)

## Code

#### 1. Load required packages

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
    library("gtsummary")
    library("mice")
    library("rstatix")
    library("robCompositions")
    library("boot")
    library("performance")

#### 2. Load and refine dataset

    load('data/FinalData/act_dat.RData', verbose = TRUE)
    skim(act_dat)

Check that variables are classified correctly (only a few examples
provided here)

    str(act_dat)
    act_dat$site <- as.factor(act_dat$site)
    act_dat$Age <- as.numeric(act_dat$Age)

Check completeness/missingness in dataset

    mice::md.pattern(act_dat, rotate.names=TRUE)

Correct brain volume variables for TIV, DC, and site (only one example
provided here, but repeated for all ROIs)

    fit.GM <- lm(GM.vol ~ TIV + site + DC, data = act_dat)
    summary(fit.GM)
    act_dat['GM.corrected'] <- act_dat$GM.vol - (predict(fit.GM, newdata=act_dat, type='response') - mean(act_dat$GM.vol))
    #visualise difference between original and corrected vols
    ggplot(act_dat, aes(x=GM.vol, y=GM.corrected)) + geom_point(size = 1.5)

Log-transform data

*Note: Completed only for lateral ventricle volume which won’t be
discussed further in this script. This was added to code once model-fit
statistics were examined later in the code*

    act_dat$vent.corrected <- log(act_dat$vent.corrected)

Look for and remove invalid accelerometry files (coded as NA or 0, as
1=valid)

    table(act_dat$valid_file, useNA = 'ifany')
    act_dat <- act_dat[act_dat$valid_file %in% 1, ] #remove participants that don't have valid accelerometry files
    table(act_dat$valid_file, useNA = 'ifany') #check that they were removed 

Create additional “Frontal.grp” variable (which separates participants
to above and below median total grey matter volume)

*Note: this step was added retrospectively (following significant
interaction) to allow for plotting by total GM volume*

    act_dat <- act_dat %>%
      mutate(GM.grp = ntile(GM.corrected, 2)) %>%
      mutate(GM.grp = if_else(GM.grp==1, 'Lower', 'Upper'))

    #check where upper and lower are separated; check difference between mean and median (point of separation)
    act_dat %>%
      mutate(y=1) %>%
      ggplot(., aes(x = GM.corrected, y = y, col = GM.grp)) +
      geom_point(alpha = 0.25)+
      theme_classic()+
      theme(axis.text.y = element_blank()) +
      geom_vline(xintercept = median(act_dat$GM.corrected)) +
      geom_vline(xintercept = mean(act_dat$GM.corrected), col = "grey50")

Refine dataset to only include required variables

    cols_covar <- c("ID", "Age", "Sex", "edu_yrs")
    cols_outcome <- c("GM.corrected", "GM.grp", "execfunc")
    cols_pred <- c("all_days_sleeptime", "all_days_sedtime", "all_days_lighttime", "all_days_mvtime")

    #combine cols_outcome and cols_pred in to cols_want dataframe
    cols_want <- c(cols_covar, cols_outcome, cols_pred)

    #create new dataframe 'actdat' which contains all rows from act_dat but only cols_want columns
    cols_want[!(cols_want %in% colnames(act_dat))]
    actdat <- act_dat[, cols_want]
    as_tibble(actdat)

Refine dataset based on time-use data

    #1. Visualise the minute values of data for each participant (is it close to 1440?)
    hist(rowSums(actdat[, cols_pred]))
    actdat[rowSums(actdat[, cols_pred]) > 1500, ] #show which rows add up to >1500 mins
    actdat[rowSums(actdat[, cols_pred]) < 1400, ] #show which rows add up to <1400 mins

    # remove participants with more than 1500 minutes of time use data
    #Note that this was an arbitrary but pragmatic decision to remove observations with potentially misleading time-use recordings made by researchers for this study.
    nrow(actdat)
    actdat <- actdat[!(rowSums(actdat[, cols_pred]) > 1500), ]

Review dataset and obtain arithmetic means (before closure to 1440
minutes)

    actdat %>%
      tbl_summary(
        statistic = list(all_continuous() ~ "{mean} ({sd})")
        # other option; statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
      )

Apply clo() function to rescale composition to 1440 minutes

    actdat[, cols_pred] <- clo(actdat[, cols_pred], total = 1440)
    rowSums(actdat[, cols_pred]) #checks that all participants have 1440 mins

Section summary: at this stage of the code, we have 1) loaded the
dataset, 2) corrected brain volume data for TIV, DC and site, 3) refined
the dataset to only include those with valid wear time indicators, 4)
summarised the dataset using tbl\_summary(), and 5) applied the clo()
function so that compositional variables have been re-scaled to sum to
1440 minutes (24h).

For the purpose of this script, it is important to note that data were
frequently inspected for normality and outliers throughout the cleaning
process. Further, the steps outlined in this section were adjusted to
include all brain volume measures and cognitive outcomes.

#### 3. Compute correlations between variables

The following code evaluates correlations between time-use variables,
continuous covariates, brain volume outcomes and cognitive outcomes. As
specified in the manuscript, the correlations between time-use variables
(i.e., within the composition) were computed using separate code based
on Kynclova, Hron and Filzmoser (2017):
<https://link.springer.com/article/10.1007/s11004-016-9669-3>

Define the variables that will be used to compute correlations

    cols_to_use <- c("Age",  "edu_yrs" , "all_days_lighttime",
              "all_days_mvtime", "all_days_sedtime",
              "all_days_sleeptime", "GM.corrected", "execfunc")

Create correlation matrix with 95% CIs Note: within-composition
correlations are not extracted/interpreted from this matrix

    nc <- length(cols_to_use)
    actdat[, cols_to_use]
    m_cor <- matrix("", nrow = nc, ncol = nc, dimnames = list(cols_to_use, cols_to_use))
    for (i in 1:(nc - 1)) { # i <- 1
      for (j in (i + 1):nc) { # j <- 3
        
        print(i)
        print(j)
        
        var_i <- actdat[[cols_to_use[i]]]
        var_j <- actdat[[cols_to_use[j]]]
          cor_ij <- cor.test(var_i, var_j, subset = !(is.na(var_i) | is.na(var_j)))
          m_cor[i, j] <- 
            sprintf("%1.3f (%1.3f, %1.3f)", cor_ij$estimate, cor_ij$conf.int[1], cor_ij$conf.int[2])
        
      }
    }
    #view correlation matrix in console
    m_cor[1:3, 1:3]

Create within-composition correlations using corCoDa() function

    corrs <- subset(actdat, select=c('all_days_lighttime', 'all_days_mvtime', 'all_days_sedtime', 'all_days_sleeptime'))
    corrs <- corCoDa(corrs)
    corrs
    corrs.v = as.vector(corrs[upper.tri(corrs)]) 
    #gives the values in the correlation matrix as a vector (goes down the cols)
    ndat = actdat %>% dplyr::select('all_days_lighttime', 'all_days_mvtime', 'all_days_sedtime', 'all_days_sleeptime')
    ndat=as.data.frame(ndat)
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
    ##these are the compositional correlations (first col) and 95% bootstrap CI (2 and 3 col)
    cbind.data.frame(corrs.v, as.data.frame(t(limper)))
    round(cbind.data.frame(corrs.v, as.data.frame(t(limper))),2)

#### 4. Create ilrs for analyses

Double-check the order of time-use variables in `cols_pred`

    cols_pred <- c("all_days_sleeptime","all_days_sedtime","all_days_lighttime", "all_days_mvtime")

Create sequential binary partition to use as ilr base

    sbp4 = matrix(c( 1, -1, -1,-1, 
                     0, 1, -1, -1,
                     0, 0, 1, -1),
                  ncol=4, byrow=TRUE)
    psi4 = gsi.buildilrBase(t(sbp4)) 

Create ilrs

    ilrs = ilr(acomp(actdat[, (cols_pred)]), V=psi4)
    colnames(ilrs) <- paste0("ilr", 1:3)

Create new dataframe `ilr_dat` to use in analyses

    ilr_dat <- cbind(actdat[, c(cols_outcome, cols_covar)], ilrs)

Simplify presentation of polynomial time-use composition

    c_ <- with(ilr_dat, cbind(ilr1, ilr2, ilr3))
    ilr_dat$poly_c <- poly(c_, degree = 2, raw = TRUE)

#### 5. Aim 1: Associations between 24-hr time-use composition and grey matter volume

The following section outlines the linear regression modelling pipelines
used to investigate aim 1. For simplicity, only the code for total grey
matter volume outcome is presented. Here is an excerpt of the manuscript
which describes the statistical modelling approach:

*First, linear regression models were used to derive the associations
between 24-hour time-use composition (predictor) and each ROI volume
(outcome). Model 1 included covariates only (age, sex, education),
whilst Model 2 included covariates and the predictor of interest
(time-use composition). Model fit was examined using the performance
package in R \[48\]. To improve the model fit, lateral ventricle volumes
were log-transformed, whilst all other model fit diagnostics passed
assumption checks and variables were therefore not transformed. To
account for the possibility that associations between time-use
composition and ROI volume outcomes may be non-linear (e.g., inverted
U-shaped relationships between sleep and ROI volume), an additional
model (Model 3) which was identical to Model 2 but expressed time-use
composition using quadratic (squared) terms was fit. We used an F-test
to explore whether quadratic terms (Model 3) improved model fit compared
to Models 1 and 2 *

Create models

    model1 <- lm(GM.corrected ~ Age + Sex + edu_yrs, dat = ilr_dat) #covariate only model
    model2 <- lm(GM.corrected ~ Age + Sex + edu_yrs + cbind(ilr1, ilr2, ilr3), dat = ilr_dat) #time-use composition + covariates
    model3 <- lm(GM.corrected ~  Age + Sex + edu_yrs + poly_c, dat = ilr_dat) #polynomial terms

Check the model fit of models 1 and 2 using performance() package

    check_model(model1)
    check_model(model2)

Summarise linear regression outputs

    summary(model1)
    summary(model2)
    summary(model3)

Use ANOVA to determine whether inclusion of squared(polynomial) terms
improves model fit

    anova(model2, model3)

*Note: polynomial terms did not improve model fit for any outcomes in
this study and won’t be described/explored further in this code
example.*

Use car() package to determine overall significance of associations
between time-use composition and total GM volume

    car::Anova(model2)

Adjust p-values for false discovery rate using Benjamini-Hochberg method
Note: where ‘0.111’, ‘0.222’, ‘0.333’, and ‘0.444’ are, replace with
p-values from car::Anova output.

    pvals <- c(0.111, 0.222, 0.333, 0.444)
    p.adjust(pvals, method = "BH", n=length(pvals))

#### 6. Aim 2: Does grey matter volume moderate the association between 24-hr time-use composition and cognitive function?

The following section outlines the linear regression modelling pipelines
used to investigate aim 2. For simplicity, only the code for total grey
matter volume outcome and executive function cognitive outcome is
presented. For each ROI volume, we tested the interaction with time-use
composition for all three cognitive outcomes. Here is an excerpt of the
manuscript which describes the statistical modelling approach:

*Next, we investigated whether the associations between time-use
composition and the cognitive function outcomes were moderated by grey
matter volume, frontal lobe volume, temporal lobe volume or hippocampus
volume. For each cognitive outcome (long-term memory, executive
function, processing speed), a series of linear regression models were
fit, with each incorporating main effects of ROI volume and time-use
composition, as well as the interaction between time-use composition and
the respective ROI volume. All models were adjusted for covariates (age,
sex, education).*

Create and summarise model

    model1 <- lm(execfunc ~ Age + Sex + edu_yrs + cbind(ilr1, ilr2, ilr3) * GM.corrected, dat=ilr_dat)
    summary(model1)

Use car() package to determine overall significance of associations

    car::Anova(model1)

Adjust p-values for false discovery rate using Benjamini-Hochberg method
Note: where ‘0.111’, ‘0.222’, ‘0.333’, ‘0.444’, and ‘0.555’ are, replace
with p-values from car::Anova output.

    pvals <- c(0.111, 0.222, 0.333, 0.444, 0.555)
    p.adjust(pvals, method = "BH", n=length(pvals))

#### 7. Plot interaction between total GM volume and time-use composition for executive function outcome

To better understand the signficant interactions, we used compositional
isotemporal substitution modelling. Here is an excerpt from the
manuscript which describes the approach:

*Similarly, where an interaction between 24-hour time-use composition
and a volumetric outcome was significantly associated with a cognitive
outcome, we planned to plot model-generated predictive response curves
to demonstrate how cognitive function was associated with meaningful
reallocations of time across different brain volume levels (dichotomized
to above or below the mean ROI volume) using one-for-remaining and
one-for-one swaps*.

The following code was adapted from the codaredistlm() package by
Stanford, Rasmussen and Dumuid (2022):
<https://cran.r-project.org/web/packages/codaredistlm/index.html>.

    act=actdat[,cols_pred]
    #make predictive compositions
    act=acomp(act)
    (m=clo(mean(act),total=1))

The following section makes one-for-remaining reallocated compositions.
For simplicity, we have presented the code used to create +15min
reallocations for the each time-use variable. This was repeated for all
other reallocation increments (e.g., +30 min, -15 min, -30 min) across
all time-use variables.

    #add 15 min Sleep
    r=15/1440/m[1]
    s=r* m[1]/(1-m[1])
    sl.15i=acomp(cbind(m[1]*(1+r),m[2]*(1-s),m[3]*(1-s),m[4]*(1-s)))
    clo(sl.15i, total=1440) #sleep has 15 min more than the mean
    #this is the mean
    clo(m, total=1440)

    #add 15 min Sb
    r=15/1440/m[2]
    s=r* m[2]/(1-m[2])
    sb.15i=acomp(cbind(m[1]*(1-s),m[2]*(1+r),m[3]*(1-s),m[4]*(1-s)))
    clo(sb.15i, total=1440)

    #add 15 min lpa
    r=15/1440/m[3]
    s=r* m[3]/(1-m[3])
    lpa.15i=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1+r),m[4]*(1-s)))
    clo(lpa.15i, total=1440)

    r=15/1440/m[4]
    s=r* m[4]/(1-m[4])
    mvpa.15i=acomp(cbind(m[1]*(1-s),m[2]*(1-s),m[3]*(1-s),m[4]*(1+r)))
    clo(mvpa.15i, total=1440)

Gather compositions that were created and make new ilrs (note that all
new reallocations will need to be listed here)

    comps <- rbind(sl.30d, sl.15d, m, sl.15i, sl.30i,
                sb.30d, sb.15d, m, sb.15i, sb.30i,
                lpa.30d, lpa.15d, m, lpa.15i, lpa.30i,
                mvpa.30d, mvpa.15d, m, mvpa.15i, mvpa.30i)
    new.ilrs=ilr(comps, V=psi4)

Recreate model for plotting

*Note: This model includes the binary ‘GM.grp’ variable. Here we can
double-check that the interaction remains significant when we express
total GM vol as a categorical variable. We have also expressed sex as a
numeric variable to allow for 1.5 to be used (rather than male or female
as the average case).*

    ilr_dat$sexN=as.numeric(ilr_dat$Sex)
    model <- lm(execfunc ~ Age + sexN + edu_yrs + 
                      GM.grp*cbind(ilr1, ilr2, ilr3), dat = ilr_dat)
    car::Anova(model)

Make new dataset to predict for.

    newdata=expand.grid(new.ilrs[,1], Age=mean(actdat$Age), sexN=1.5, edu_yrs=mean(actdat$edu_yrs), GM.grp=c('Upper', 'Lower'))
    newdata$ilr1=new.ilrs[,1]
    newdata$ilr2=new.ilrs[,2]
    newdata$ilr3=new.ilrs[,3]

Use newdata to predict executive function

    pred.EF=predict(model, newdata, interval="confidence")

Manipulate data, create `plotdata` dataframe

    #add original reallocated compositions to the plotting dataset for later
    Activity=clo(comps,total=1440)
    plotdata=cbind(newdata, pred.EF, sleep=Activity[,1], sb=Activity[,2], lpa=Activity[,3], mvpa=Activity[,4])
    head(plotdata)
    #Create variable for plotting and indicate the "one" variable in the "one-to-remaining" reallocations
    #the '5' denotes that for each time-use variable we are computing 5 reallocations (including zero)
    plotdata$Activity <- rep(c(rep("Sleep",5), rep("SB",5), rep("LPA",5), rep("MVPA",5)),2)
    plotdata$Activity=factor(plotdata$Activity, levels=c("Sleep","SB","LPA","MVPA"))
    plotdata$Reallocation=c(rep(seq(-30,30,15),8))
    head(plotdata)
    plotdata$pred.EF = plotdata$fit
    plotdata$GM.grp=factor(plotdata$GM.grp, levels = c("Lower", "Upper"))#make a factor to help plotting

Plot reallocations separately for ‘upper’ and ‘lower’ total GM volume
groups

    ggplot(plotdata, aes(x=Reallocation, y=pred.EF, group=GM.grp))+
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

#### Please contact Dr Maddison Mellow (<maddison.mellow@unisa.edu.au>) if you have any questions about the code or methodology.
