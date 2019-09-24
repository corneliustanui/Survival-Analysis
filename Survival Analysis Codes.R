#---------------SURVIVAL ANALYSIS--------------------------
library(survival)
Data <- lung
summary(Data)

# The three main components of survival analysis are:
# a) Time to event, in this data set it is "time"
# b) The grouping variable(s)(any categorical IVs, e.g. sex, religion, etc.)
# c) Event, this is binary with 1 = success (outcome observed), 0 = failure
# Hence change "status" to 0|1

Data$status <- ifelse(Data$status == 2, 1, 0)

# 1) DEFINITION OF TERMS----

# In survival package the function "surv" takes agruments:
# "time" = Start time for interval data
# "time2" = End time for interval data
# "event" = the event of interest
# "type" = the type of censoring, i.e. "right", "left", "interval", "counting", "mstate"
# "origin" = origin of hazard function, (s(x)) for counting process data.

# "right censored" means subjects leave the study before time to event
# "left censored" means subjects join the study after it already started.
# "interval" means "time" is recodred in ranges, e.g. (1-20-hours, where t1 = 1 hr, t2 = 20 hrs) etc.
# "counting" means
# "mstate" means the event is a "multi-state factor", i.e several possible outcomes/levels, where the firt level represents censoring.

# 2) DECLARING SURVIVAL DATA----
survive <- with(Data, Surv(time = time, event = status))
summary(survive)
survive # The plus (+) in the results stand for censored cases (no outcome)

# 3) NON-PARAMETRIC ESTIMATION OF S(X)----
# These models have no assumptions on the distribution of factors.
# Do not take in continuous predictors.
# Do not take in more than one predictors.
# They are not regression models.

# 3.a) Kaplan-Maier Estimation----

# Single S(X) with no comparison (no predictor) (use ~1)
KM_1 <- survfit(survive~1)
KM_1
summary(KM_1)

plot(KM_1, main = "Overall S(X)", ylab = "Time to event", xlab = "S(X)", conf.int = "none") # Without 95%CI
plot(KM_1, main = "Overall S(X)", ylab = "Time to event", xlab = "S(X)") 

# NB: Median Survival Time occurs when s(x) ==0.5

abline(h = 0.5, v = 310, col = "red") # Plotting median tome and mean s(x)

# S(X) with comparison/Predictor(use ~x)
Data <- within(Data, {
  sex <- factor(sex,
                levels = 1:2,
                labels = c("Male","Female"))})

KM_2<- survfit(survive~sex, data = Data)
KM_2
summary(KM_2)

plot(KM_2, 
     main = "S(X) for Sex",   # Title
     ylab = "Time to event",  # Y-axis label
     xlab = "S(X)",           # X-axis label
     conf.int = "both",       # To plot 95%CI for both groups, default is "none".
     col = c("red","blue"),   # Color to distinguish groups
     lwd = 2,                 # Plot a line thicker than default
     mark = 3)                # To put marks on censored regions, "3" means "vertical" lines.
legend("topright", 
       legend = c("Male", "Female"), 
       col = c("red", "blue"), 
       lwd = 2)
abline(h = 0.5,              # Mean s(x)
       v = c(270, 426),      # Median times
       col = c("black", "red", "blue"), 
       lty = 2)              # Broken lines
#____________________________________________________________
# Inverse Plot of S(x)
plot(KM_2, 
     fun = "event",         # To flip Y-axis
     main = "S(X) for Sex",   
     ylab = "Time to event",  
     xlab = "S(X)",           
     conf.int = "both",       
     col = c("red","blue"),   
     lwd = 2,                 
     mark = 3)                
legend("topleft", 
       legend = c("Male", "Female"), 
       col = c("red", "blue"), 
       lwd = 2)
#______________________________________________________________
#______________________________________________________________
broom::tidy(KM_2) # Clean the output

ggplot2::autoplot(KM_2) # Requires "ggfortify" package to work.
#_____________________________________________________________

# TESTING WHETHER THE DIFFERENCE IN SURVIVAL TIME BTWN THE TWO GROUPS IS SIGNIFICANT OR DUE TO CHANCE
# A) Log-Rank Test, also called Mantel-Haenszel Test (rho = 0) (More weight on higher values of time)
Diff_1 <- survdiff(survive~sex, data = Data, rho = 0)
Diff_1

# B) Generalised Wilcoxon (rho = 1) (More weight on lower values of time)
Diff_2 <- survdiff(survive~sex, data = Data, rho = 1)
Diff_2

# C) Tarobe-Ware Test (rho = 1.5) (In between Log-rank and wilcoxon)
Diff_3 <- survdiff(survive~sex, data = Data, rho = 1.5)
Diff_3

# NB: p-value < 0.05 implies there is enough evidence that the groups differ significantly. 

# 3.b) Nelson-Aalen Estimation----
N_A <- survfit(coxph(survive~sex, data = Data), type = "aalen")
N_A
summary(N_A)

# 4) SEMI-PARAMETRIC ESTIMATION OF S(X)----
# 4.a) Cox Proportional Hazards----
# Allows inclusion of multiple predictors, both categorical and continuous.
# Propotional Hazard Assumption:
# i) Constant Hazard Ratio over time
# ii) Groups compared therefore diverge at the same rate (or are parallel), but never converge nor cross.
# (Bj) < 0 means increase in the respective X is associated with lower log(HR) and longer survival time.
# (Bj) >0 means increase in the respective X is associated with higher log(HR) and shorter survival time.

# exp(Bj) == HR;

# HR < 1 means increase in the respective X is associated with lower HR and longer survival time.
# HR > 1 means increase in the respective X is associated with higher HR and shorter survival time.

# To test Ho:B = 0, we can use Wald test,  or Likelihood  Ratio Test(LRT) if events are fewer. 

Cox <- coxph(survive~sex, data = Data)                          # Handles ties by deafult by "Efron" method
Cox_Bres <- coxph(survive~sex, data = Data, method = "breslow") # Handles ties by "Breslow" method.

Cox
summary(Cox)    # Reports LRT and Wald.
anova(Cox)      # Reports logLTR.

# 5) PARAMETRIC ESTIMATION OF S(X)----
# 5.a) Exponential Estimation----
Expo <- survreg(survive~sex, data = Data, dist =  "exponential")
Expo
summary(Expo)

# 5.b) Weibull Estimation----
Weibull <- survreg(survive~sex, data = Data, dist =  "weibull")
Weibull
summary(Weibull)

# 5.c) Loglogistic Estimation----
loglogistic <- survreg(survive~sex, data = Data, dist =  "loglogistic")
loglogistic
summary(loglogistic)

# 5.d) Gaussian Estimation----
gaussian <- survreg(survive~sex, data = Data, dist =  "gaussian")
gaussian
summary(gaussian)

# 5.e) Lognormal Estimation----
lognormal <- survreg(survive~sex, data = Data, dist =  "lognormal")
lognormal
summary(lognormal)

#____________________________END_____________________________________
