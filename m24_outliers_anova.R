cat("\014")

# Load in libraries
library(emmeans)
library(nlme)

data <- read.csv("*")

data$Subject <- as.factor(data$Subject)
data$Sequence <- as.factor(data$Sequence)
data$Order_S <- as.factor(data$Order_S)
data$Run <- as.factor(data$Run)


# Set parameters for linear mixed effect model
ctrl <<- lmeControl(maxIter = 1000, msMaxIter = 1000, optimMethod = "optim", msMaxEval = 10000, niterEM = 10000)

# Fit linear mixed effect model
model_out <- lme(data=data, Outliers ~ Sequence + Order_S + Run, random = ~ 1 + Sequence + Order_S + Run | Subject, na.action = na.omit, control = ctrl)
print(summary(model_out))

emm <- emmeans(model_out, ~ Sequence)
print(pairs(emm))

emm <- emmeans(model_out, ~ Order_S)
print(pairs(emm))

emm <- emmeans(model_out, ~ Run)
print(pairs(emm))


# #ANOVA (NOT USED!)
# res.aov <- aov(Outliers ~ Sequence + Order_S + Run + Error(Subject/(Sequence + Order_S + Run)), data=data)
# print(summary(res.aov))
#
# emm <- emmeans(res.aov, ~ Sequence)
# print(pairs(emm))
#
# emm <- emmeans(res.aov, ~ Order_S)
# print(pairs(emm))
#
# emm <- emmeans(res.aov, ~ Run)
# print(pairs(emm))
