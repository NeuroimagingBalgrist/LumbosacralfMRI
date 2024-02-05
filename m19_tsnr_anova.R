cat("\014")

library(emmeans)

data <- read.csv("/*/sublvl_tsnr/tsnr_anova_data.txt")

data$Subject <- as.factor(data$Subject)
data$Sequence <- as.factor(data$Sequence)
data$Slice <- as.factor(data$Slice)


#Two-way ANOVA
res.aov <- aov(tSNR ~ Sequence + Slice + Error(Subject/(Slice +Sequence)) , data=data)
print(summary(res.aov))

emm <- emmeans(res.aov, ~ Sequence)
print(pairs(emm))

emm <- emmeans(res.aov, ~ Slice)
print(pairs(emm))

# res.aov <- aov(tSNR ~ Sequence * Slice + Error(Subject/(Sequence*Slice)) , data=data)
# print(summary(res.aov))
# emm <- emmeans(res.aov, ~ Sequence, at = list(Slice = "3"))
# pairs(emm)
