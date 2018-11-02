# Lesson 15: Logistic, Poisson & Nonlinear Regression
# https://onlinecourses.science.psu.edu/stat501/node/433/
# Leukemia remission (logistic regression)
# 
# Load the leukemia data.
# Fit a logistic regression model of REMISS vs CELL + SMEAR + INFIL + LI + BLAST + TEMP.
# Calculate 95% confidence intervals for the regression parameters based on asymptotic normality and based on profiling the least-squares estimation surface.
# Fit a logistic regression model of REMISS vs LI.
# Create a sctterplot of REMISS vs LI and add a fitted line based on the logistic regression model.
# Calculate the odds ratio for LI and a 95% confidence interval.
# Conduct a likelihood ratio (or deviance) test for LI.
# Calculate the sum of squared deviance residuals and the sum of squared Pearson residuals.
# Use the hoslem.test function in the ResourceSelection package to conduct the Hosmer-Lemeshow goodness-of-fit test.
# Calculate a version of R2 for logistic regression.
# Create residual plots using Pearson and deviance residuals.
# Calculate hat values (leverages), studentized residuals, and Cook's distances.

leukemia <- read.table ("DataL15/leukemia_remission.txt", header=T)
attach(leukemia)

model.1 <- glm(REMISS ~ CELL + SMEAR + INFIL + LI + BLAST + TEMP, family="binomial")
summary(model.1)

confint.default(model.1) # based on asymptotic normality
confint(model.1) # based on profiling the least-squares estimation surface

model.2 <- glm(REMISS ~ LI, family="binomial")
summary(model.2)

plot(x=LI, y=REMISS,
     panel.last = lines(sort(LI), fitted(model.2)[order(LI)]))

exp(coef(model.2)[2]) # odds ratio = 18.12449
exp(confint.default(model.2)[2,]) # 95% CI = (1.770284, 185.561725)

anova(model.2, test="Chisq")

sum(residuals(model.2, type="deviance")^2) # 26.07296
model.2$deviance # 26.07296
sum(residuals(model.2, type="pearson")^2) # 23.93298

library(ResourceSelection)
hoslem.test(model.2$y, fitted(model.2), g=9)

1-model.2$deviance/model.2$null.deviance # "R-squared" = 0.2414424

plot(1:27, residuals(model.2, type="pearson"), type="b")
plot(1:27, residuals(model.2, type="deviance"), type="b")

summary(influence.measures(model.2))

hatvalues(model.2)[8] # 0.1498395
residuals(model.2)[8] # -1.944852
rstudent(model.2)[8] # -2.185013
cooks.distance(model.2)[8] # 0.5833219

detach(leukemia)

# Disease outbreak (logistic regression)
# 
# Load the disease outbreak data.
disease <- read.table("DataL15/DiseaseOutbreak.csv", sep=",", header=T)
attach(disease)
head(disease)

# Create interaction variables.
Age.Middle <- Age*Middle
Age.Lower <- Age*Lower
Age.Sector <- Age*Sector
Middle.Sector <- Middle*Sector
Lower.Sector <- Lower*Sector

# Fit "full" logistic regression model of Disease vs four predictors and five interactions.
model.1 <- glm(Disease ~ Age + Middle + Lower + Sector + Age.Middle + Age.Lower +
                 Age.Sector + Middle.Sector + Lower.Sector, family="binomial")

# Fit "reduced" logistic regression model of Disease vs four predictors.
model.2 <- glm(Disease ~ Age + Middle + Lower + Sector, family="binomial")

# Conduct a likelihood ratio (or deviance) test for the five interactions.
anova(model.2, model.1, test="Chisq")

# Display the analysis of deviance table with sequential deviances.
anova(model.1, test="Chisq")
detach(disease)

###########################################################################
# Toxicity and insects (logistic regression using event/trial data format)#
###########################################################################
# Load the toxicity data.
toxicity <- read.table("DataL15/toxicity.csv", sep=",",  header=T)
attach(toxicity)

# Create a Survivals variable and a matrix with Deaths in one column and Survivals in the other column.
Survivals <- SampSize - Deaths
y <- cbind(Deaths, Survivals)

# Fit a logistic regression model of Deaths vs Dose.
model.1 <- glm(y ~ Dose, family="binomial")
summary(model.1)

# Calculate 95% confidence intervals for the regression parameters based on asymptotic normality and based on profiling the least-squares estimation surface.
confint.default(model.1) # based on asymptotic normality
confint(model.1) # based on profiling the least-squares estimation surface


# Calculate the odds ratio for Dose and a 95% confidence interval.
exp(coef(model.1)[2]) # odds ratio = 1.962056
exp(confint.default(model.1)[2,]) # 95% CI = (1.817279, 2.118366)


# Display the observed and fitted probabilities.
cbind(Dose, SampSize, Deaths, Deaths/SampSize, fitted(model.1))

# Create a sctterplot of observed probabilties vs Dose and add a fitted line based on the logistic regression model.
plot(x=Dose, y=Deaths/SampSize,
     panel.last = lines(sort(Dose), fitted(model.1)[order(Dose)]))

detach(toxicity)

### Poisson example (Poisson regression) ###
# 
# Load the poisson data.
poisson <- read.table("DataL15/poisson_simulated2.csv", sep=",",  header=T)
attach(poisson)
head(poisson)
names(poisson)
str(poisson)
# Create a scatterplot of the data.
with (data=poisson, plot(x=x, y=y)) # somehow attach() does not link well

# Fit a Poisson regression model of y vs x.
model.1 <- glm(y ~ x, family="poisson", data=poisson)
summary(model.1)

# Calculate 95% confidence intervals for the regression parameters based on asymptotic normality and based on profiling the least-squares estimation surface.
confint.default(model.1) # based on asymptotic normality
confint(model.1) # based on profiling the least-squares estimation surface

# Create a sctterplot of y vs x and add a fitted line based on the Poisson regression model.
with(data=poisson, plot(x=x, y=y,
     panel.last = lines(sort(x), fitted(model.1)[order(x)])))

# Conduct a likelihood ratio (or deviance) test for x.
anova(model.1, test="Chisq")

# Calculate the sum of squared deviance residuals and the sum of squared Pearson residuals and calculate p-values based on chi-squared goodness-of-fit tests.
sum(residuals(model.1, type="deviance")^2) # 27.84209
model.1$deviance # 27.84209
pchisq(model.1$deviance, 28, lower.tail=F) # p-value = 0.4728389

sum(residuals(model.1, type="pearson")^2) # 26.09324
pchisq(sum(residuals(model.1, type="pearson")^2), 28, lower.tail=F) # p-value = 0.5679192

# Calculate pseudo R2 for Poisson regression.
1-model.1$deviance/model.1$null.deviance # Pseudo R-squared = 0.423676

# Create residual plots using Pearson and deviance residuals.
plot(fitted(model.1), residuals(model.1, type="pearson"))
plot(fitted(model.1), residuals(model.1, type="deviance"))

summary(influence.measures(model.1))

# Calculate hat values (leverages) and studentized residuals.
residuals(model.1)[8] # 1.974329
rstudent(model.1)[8] # 2.028255

detach(poisson)

### Hospital recovery (exponential regression) ###
# 
# Load the recovery data.
recovery <- read.table("DataL15/recovery.csv", sep=",", header=T)
attach(recovery)
head(recovery)

# Create log(prog) variable.
logprog <- log(prog)
summary(lm(logprog ~ days))

# Obtain starting values for nonlinear model parameters from fitting a simple linear regression model of log(prog) vs days.
exp(4.037159) # 56.66513

# Fit nonlinear regression model to data using these starting values.
model.1 <- nls(prog ~ theta1 * exp(theta2 * days),
               start=list(theta1=56.7, theta2=-0.038))
summary(model.1)

# Create a scatterplot of prog vs days and add a fitted line based on the nonlinear regression model.
plot(x=days, y=prog,
     panel.last = lines(sort(days), fitted(model.1)[order(days)]))

detach(recovery)

### U.S. census population (population growth nonlinear regression) ###
# 
# Load the census data.
census <- read.table("DataL15/us_census.txt", header=T)
attach(census)
head(census)

# Obtain starting values for nonlinear model parameters from observing features of a scatterplot of population vs year.
plot(x=year, y=population)
log(350/3.929-1) # 4.478259
log(350/5.308-1) - log(350/3.929-1) # -0.3048229

# Fit nonlinear regression model to data using these starting values.
model.1 <- nls(population ~ beta1 / (1 + exp(beta2 + beta3 * (year - 1790) / 10)),
               start=list(beta1=350, beta2=4.5, beta3=-0.3))
summary(model.1)

# Create a scatterplot of population vs year and add a fitted line based on the nonlinear regression model.
plot(x=year, y=population,
     panel.last = lines(sort(year), fitted(model.1)[order(year)]))

# Create a residual plot.
plot(x=year, y=residuals(model.1),
     panel.last = abline(h=0, lty=2))

detach(census)
