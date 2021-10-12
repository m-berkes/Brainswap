########### Libraries #######################################

# Check.Packages
# Function to simplify package install
# Requires 'packages' as array of strings using 'c'

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# check dependencies and install if needed. These are just a bunch of packages I have found useful in the past, not necessarily used here.

packages<-c("ggplot2",  "plyr", "dplyr", "tidyverse", "car", "DescTools", "plotrix", "lsr", "e1071", "patchwork", "RColorBrewer",
            "emmeans", "ggpubr", "fastmatch", "BayesianFirstAid", "extrafont", "MatchIt", "FactoMineR", "modelr", "infer", "interactions",
            "factoextra")
check.packages(packages)

# run this after installing new fonts on system
font_import()
loadfonts()

## needed for extrafonts if saving graph to PDF
## Needed only on Windows - run once per R session
## Adjust the path to match your installation of Ghostscript - so obviously you also need to download Ghostscript
Sys.setenv(R_GSCMD = "D:/Program Files/gs9.26/bin/gswin32c.exe")


########### Set up #####################

## load file
brainswap <- read_csv(
  file = "brainswap.csv",         #change this location
  col_names=TRUE, na = "."
  )

brainswap$LangGroup <- factor(brainswap$LangGroup, levels = c("ML", "BL"))
brainswap$CogProfile <- factor(brainswap$CogProfile, levels = c("CN", "MCI", "AD"))
brainswap$Sex <- factor(brainswap$Sex, levels = c("M", "F"))
brainswap$SiteID <- factor(brainswap$SiteID)

# prep dataframes for matching
# MatchIt requires complete cases only
swapcomplete <- brainswap[complete.cases(brainswap), ]

## change language group to a new logical variable - FALSE is untreated (monolingual) and TRUE is treated (bilingual)
swapcomplete$Group <- swapcomplete$LangGroup == "BL"

# assign new numbers based on CogProfile: 0 is healthy, 1 is MCI or AD
swapcomplete$CogProfileNum <- as.numeric(swapcomplete$CogProfile)
swapcomplete$CogProfileNum[ swapcomplete$CogProfileNum == 1 ] <- 0
swapcomplete$CogProfileNum[ swapcomplete$CogProfileNum > 1 ] <- 1


################## Base descriptives ################

ddply(swapcomplete, .(LangGroup), summarise, 
      N = length(LangGroup),
      Age_mean = mean(Age),
      Age_sd = sd(Age),
      FA_mean = mean(FA),
      FA_sd = sd(FA),
      AD_mean = mean(AD), 
      AD_sd = sd(AD),
      RD_mean = mean(RD),
      RD_sd = sd(RD),
      Edu_mean = mean(Edu, na.rm = TRUE),
      Edu_sd = sd(Edu, na.rm = TRUE),
      MMSE_mean = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      CogProfile_mean = mean(CogProfileNum),
      CogProfile_sd = sd(CogProfileNum)
      )



############### PCA ########################

# PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE), requires FactoMiner package

for.pca <- subset(swapcomplete, select = c("FA", "AD", "RD"))
res.pca <- PCA(for.pca, graph = TRUE)

# because the first PC score accounts for almost 80% of variance, we only use PC1 in analyses
scoresPCA <- res.pca$ind$coord[,1]
swapcompletePCA <- cbind(swapcomplete, scoresPCA)

ddply(swapcompletePCA, .(LangGroup), summarise, 
      PCA = mean(scoresPCA),
      PCA_sd = sd(scoresPCA),
      PCA2 = mean(scoresPCA2),
      PCA2_sd = sd(scoresPCA2)
      )

# match based on PCA scores
matchedPCA <- matchit(Group ~ Sex + Edu + Age + scoresPCA,
                     data = swapcompletePCA, method = "nearest", distance = "logit", discard = "treat")
plot(matchedPCA, type = 'jitter', interactive=FALSE)
plot(matchedPCA, type = "hist")
matchedPCA.out <- matchedPCA$match.matrix     # this shows which subjects are matched together
summary(matchedPCA)                           # information about pre and post match
df.matchPCA <- match.data(matchedPCA)[1:ncol(swapcompletePCA)]



# PCA descriptives
ddply(df.matchPCA, .(LangGroup), summarise, 
      N = length(LangGroup),
      Age_m = mean(Age),
      Age_sd = sd(Age),
      FA_mean = mean(FA),
      FA_sd = sd(FA),
      AD_mean = mean(AD), 
      AD_sd = sd(AD),
      RD_mean = mean(RD),
      RD_sd = sd(RD),
      Edu_mean = mean(Edu, na.rm = TRUE),
      Edu_sd = sd(Edu, na.rm = TRUE),
      MMSE_mean = mean(MMSE, na.rm = TRUE),
      MMSE_sd = sd(MMSE, na.rm = TRUE),
      CogProfileNum = mean(CogProfileNum),
      PCA = mean(scoresPCA),
      PCA_sd = sd(scoresPCA)
    )


# list of covariates
matchPCA_cov <- c('Edu', 'Age', 'scoresPCA')

# t tests comparing covariates means
lapply(matchPCA_cov, function(v) {
  t.test(df.matchPCA[, v] ~ df.matchPCA$LangGroup, var.equal=TRUE)
})

# cohensd - requires lsr package
cohensD(Edu ~ LangGroup, data = df.matchPCA, method = "pooled")
cohensD(Age ~ LangGroup, data = df.matchPCA, method = "pooled")
cohensD(scoresPCA ~ LangGroup, data = df.matchPCA, method = "pooled")


# matched group comparisons
mod1 <- lm(MMSE ~ LangGroup, data = df.matchPCA)
summary(mod1)
EtaSq(mod1)
mod2 <- lm(scoresPCA ~ LangGroup, data = df.matchPCA)
summary(mod2)
EtaSq(mod2)
mod3 <- lm(Edu ~ LangGroup, data = df.matchPCA)
summary(mod3)
EtaSq(mod3)
mod4 <- lm(Age ~ LangGroup, data = df.matchPCA)
summary(mod4)
EtaSq(mod4)

chisq_test(df.matchPCA, Sex ~ LangGroup)


# examining site effects using linear models
MLcomplete <- subset(swapcompletePCA, LangGroup =="ML")

mod_site <- lm(scoresPCA ~ SiteID, data = df.matchPCA)
summary(mod_site)
Anova(mod_site, type = 2)
EtaSq(mod_site)
leveneTest(scoresPCA~SiteID, df.matchPCA)
lsmeans(mod_site, pairwise~SiteID, adjust="tukey")


############## Permutations #################

## we want to run 1000 permutations of a random sample of MLs to see the distribution of cognitively healthy to impaired individuals

# create subset with only MLs
MLcomplete <- subset(swapcompletePCA, LangGroup =="ML")


## calculate the proportion of unhealthy matched ML cognitive profiles
# edit matched sample to include string variable of brain health
df.matchPCA <- df.matchPCA %>%
  mutate(Health = case_when(
    CogProfile %in% c("CN") ~ "Healthy",
    CogProfile %in% c("MCI", "AD") ~ "Unhealthy"
  )) %>%
  select(ID:scoresPCA, Health)

# create dataset of ML only matched sample
df.matchML <- subset(df.matchPCA, LangGroup == "ML")

# mean of unhealthy ML matched sample
p_hat <- df.matchML %>%
  summarize(mean(Health == "Unhealthy")) %>%
  pull()

# edit ML only dataset of full ADNI sample
MLcomplete <- MLcomplete %>%
  mutate(Health = case_when(
    CogProfile %in% c("CN") ~ "Healthy",
    CogProfile %in% c("MCI", "AD") ~ "Unhealthy"
  )) %>%
  select(ID:scoresPCA, Health)


# compute 1000 random samples
library(infer)
set.seed(2019)
null_distn <- MLcomplete %>%
  specify(response = Health, success = "Unhealthy") %>%
  hypothesize(null = "point", p = 0.2732919) %>%
  generate(reps = 1000, type = "simulate") %>%
  calculate(stat = "prop")


# bar graph of sample distribution
ggplot(data = null_distn, mapping = aes(x = stat)) +
  geom_bar() +
  geom_vline(xintercept = p_hat, color = "red")

visualise(null_distn) +
  shade_p_value(obs_stat = p_hat, direction = "two_sided")

# get p-value from this randomisation-based test for a single proportion

null_distn %>%
  summarize(p_value = mean(stat >= p_hat) * 2)

# same result, different method
null_distn %>%
  get_p_value(obs_stat = p_hat, direction = "two_sided")



# plot of null distribution to publication specifications

png("Figure_2.png", units="mm", width = 80, height = 55.8, res=300)
ggplot(null_distn, aes(x = stat, colour = stat)) +
  geom_bar(fill = "black", alpha = 0.7) +
  geom_vline(xintercept = p_hat, color = "indianred2", size = 1) +
  geom_vline(xintercept = 0.12, color = "cornflowerblue", size = 1) +
  annotate("text", x = p_hat+.008, y = 45, label = "Matched Monolingual Sample", srt = 270, colour = "grey10", size = 2.25, fontface = 2, family = "Segoe UI") +
  annotate("text", x = 0.12+.008, y = 60, label = "Bilingual Sample", srt = 270, colour = "grey10", size = 2.25, fontface = 2, family = "Segoe UI") +
  labs(x = "\nMean Cognitive Profile", y = "Count\n") +
  theme(
    text = element_text(family = "Segoe UI"),
    axis.title.x = element_text(colour="#08306B", size = 10, face = 2, lineheight = 0.25),
    axis.title.y = element_text(colour="#08306B", size = 10, face = 2, lineheight = 0.25),
    axis.text.x = element_text(colour="#08306B", size = 8),
    axis.text.y = element_text(colour="#08306B", size = 8),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(linetype = 'dashed', colour = "grey90"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()



# original code of null distribution
ggplot(null_distn, aes(x = stat, colour = stat)) +
  geom_bar(fill = "black", alpha = 0.6) +
  geom_vline(xintercept = p_hat, color = "indianred2", size = 1) +
  geom_vline(xintercept = 0.12, color = "cornflowerblue", size = 1) +
  annotate("text", x = p_hat+.007, y = 52, label = "Matched Monolingual Sample", srt = 270, colour = "grey10", size = 5, fontface = 2, family = "Segoe UI") +
  annotate("text", x = 0.12+.007, y = 64, label = "Bilingual Sample", srt = 270, colour = "grey10", size = 5, fontface = 2, family = "Segoe UI") +
  labs(x = "\nMean Cognitive Profile", y = "Count\n") +
  theme(
    text = element_text(family = "Segoe UI"),
    plot.title = element_text(colour="grey10", size = 18, face="bold"),
    axis.title.x = element_text(colour="#08306B", size = 12, face = 2, lineheight = 0.75),
    axis.title.y = element_text(colour="#08306B", size = 12, face = 2, lineheight = 0.75),
    axis.text.x = element_text(colour="#08306B", size = 10),
    axis.text.y = element_text(colour="#08306B", size = 10),
    legend.title = element_blank(),
    legend.text = element_text(colour="grey10",size = 12),
    legend.position = c(.75, .75),
    legend.justification = c("right", "top"),
    legend.background=element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(linetype = 'dashed', colour = "grey80"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )


plot(density(df.test$Mean))






############# PCA plot ##################

png("Figure_1.png", units="mm", width = 80, height = 80, res=300)
fviz_pca_var(res.pca, col.var = "steelblue", label="none", arrowsize = 1) + 
  labs(x = "\nDimension 1 (78.5%)", y = "Dimension 2 (16.2%)\n") +
  annotate("text", x = -0.89, y = 0.44, label = "FA", colour = "steelblue", size = 5, fontface = 1, family = "Segoe UI") +
  annotate("text", x = 0.8, y = 0.67, label = "DA", colour = "steelblue", size = 5, fontface = 1, family = "Segoe UI") +
  annotate("text", x = 0.91, y = -0.05, label = "DR", colour = "steelblue", size = 5, fontface = 1, family = "Segoe UI") +
  theme(text = element_text(family = "Segoe UI"),
        plot.title = element_blank(),
        axis.title.x = element_text(colour="#08306B", size = 10, face = 2, lineheight = 0.25),
        axis.title.y = element_text(colour="#08306B", size = 10, face = 2, lineheight = 0.25),
        axis.text = element_text(colour="#08306B", size = 8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed')
  )
dev.off()



############ Graphing prep ################

## from this point on includes some additional graphing code that was not used in the final result. For curiosity only

# rename levels of language group
levels(swapcompletePCA$LangGroup)[levels(swapcompletePCA$LangGroup)=="ML"] <- "Monolingual"
levels(swapcompletePCA$LangGroup)[levels(swapcompletePCA$LangGroup)=="BL"] <- "Bilingual"
levels(df.matchPCA$LangGroup)[levels(df.matchPCA$LangGroup)=="ML"] <- "Monolingual"
levels(df.matchPCA$LangGroup)[levels(df.matchPCA$LangGroup)=="BL"] <- "Bilingual"

# rename education
names(swapcompletePCA)[names(swapcompletePCA) == 'Edu'] <- 'Education'
names(df.matchPCA)[names(df.matchPCA) == 'Edu'] <- 'Education'

#transformations
complete_trans <- swapcompletePCA
complete_trans$ref_MMSE <- 31-complete_trans$MMSE #reflect MMSE values to make positive skew
complete_trans$ln_r_MMSE <- log1p(complete_trans$ref_MMSE) #ln transform
complete_trans$ln_MMSE <- 4.1-complete_trans$ln_r_MMSE #reflect ln transform to get accurate sign on correlation value

match_trans <- df.matchPCA
match_trans$ref_MMSE <- 31-match_trans$MMSE
match_trans$ln_r_MMSE <- log1p(match_trans$ref_MMSE)
match_trans$ln_MMSE <- 4.091042-match_trans$ln_r_MMSE






######### Interaction plot ##############

# requires Interactions package

# First set up your model below:
# in our case, we'll want to know about the effect of white matter on MMSE scores at different levels of education.


#model using transformed MMSE with matched sample
model_match_trans <- lm(ln_MMSE ~ scoresPCA * Education, data = match_trans)
summary(model_match_trans)

model_plot_match_trans <- interact_plot(model_match_trans, pred = scoresPCA, modx = Education, plot.points = TRUE, jitter = 0.2, data = match_trans) +
  labs(x="\nWhite Matter\nPrincipal Component", y = "MMSE\n(Transformed)\n"
       #title="White Matter and Education Predicting MMSE"
  ) +
  scale_y_continuous(limits = c(0.6, 3.8)) +
  theme(legend.position = 'top') +
  coord_fixed(ratio = 4/3)



model <- lm(MMSE ~ scoresPCA * Education, data = df.matchPCA)
summary(model)

model_plot <- interact_plot(model, pred = scoresPCA, modx = Education, plot.points = TRUE, hitter = 0.2, data = df.matchPCA) +
  labs(title="White Matter and Education Predicting MMSE", x="\nWhite Matter\nPrincipal Component", y = "MMSE\n")

model_full <- lm(MMSE ~ scoresPCA * Education, data = swapcompletePCA)
summary(model_full)

model_plot_full <- interact_plot(model_full, pred = scoresPCA, modx = Education, plot.points = TRUE, jitter = 0.2, data = swapcompletePCA) +
  labs(x="\nWhite Matter\nPrincipal Component", y = "MMSE\n"
       #title="White Matter and Education Predicting MMSE"
      ) +
  scale_y_continuous(limits = c(9, 31)) +
  theme(legend.position = 'top')
  +  coord_fixed(ratio = 1/2)

#model using transformed MMSE dataset
model_full_trans <- lm(ln_MMSE ~ scoresPCA * Education, data = complete_trans)
summary(model_full_trans)
model_plot_full_trans <- interact_plot(model_full_trans, pred = scoresPCA, modx = Education, plot.points = TRUE, data = complete_trans) +
  labs(title="White Matter and Education Predicting MMSE", x="White Matter\nPrincipal Component", y = "MMSE")



##correlations by full sample
cor.test(swapcompletePCA$MMSE, swapcompletePCA$Education, method = c("pearson"))
cor.test(~ MMSE + Education, alternative="two.sided", data=swapcompletePCA)


corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))
}

ddply(swapcompletePCA, .(LangGroup), summarise,
      z=corfun(MMSE, Education)$statistic,
      pval=corfun(MMSE, Education)$p.value,
      tau.est=corfun(MMSE, Education)$estimate,
      alt=corfun(MMSE, Education)$alternative
    )


##correlations by matched sample
cor.test(df.matchPCA$MMSE, df.matchPCA$Education, method = c("pearson"))
cor.test(~ MMSE + Education, alternative="two.sided", data=df.matchPCA)


ddply(df.matchPCA, .(LangGroup), summarise,
      z=corfun(MMSE, Education)$statistic,
      pval=corfun(MMSE, Education)$p.value,
      tau.est=corfun(MMSE, Education)$estimate,
      alt=corfun(MMSE, Education)$alternative
      )


# MMSE and education correlation by laguage group
cor.test(match_trans$ln_MMSE, match_trans$Education, method = c("pearson"))

ddply(match_trans, .(LangGroup), summarise,
      z=corfun(ln_MMSE, Education)$statistic,
      pval=corfun(ln_MMSE, Education)$p.value,
      tau.est=corfun(ln_MMSE, Education)$estimate,
      alt=corfun(ln_MMSE, Education)$alternative
      )


#MMSE and PCA correlation by language group
cor.test(match_trans$ln_MMSE, match_trans$scoresPCA, method = c("pearson"))

ddply(match_trans, .(LangGroup), summarise,
      z=corfun(ln_MMSE, scoresPCA)$statistic,
      pval=corfun(ln_MMSE, scoresPCA)$p.value,
      tau.est=corfun(ln_MMSE, scoresPCA)$estimate,
      alt=corfun(ln_MMSE, scoresPCA)$alternative
      )


############# Correlation plots ############

#scatter plot of correlation for transformed MMSE and Education by group using matched sample
cor_edu_match <- ggscatter(match_trans, x = "Education", y = "ln_MMSE",
                           add = "reg.line", conf.int = TRUE, 
                           cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 19, label.y = 2.4),
                           xlab = "\nEducation", #ylab = "MMSE\n(Log Transformed)\n",
                           color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  scale_x_continuous(breaks=c(12,14,16,18,20,22)) + scale_y_continuous(breaks=c(1,2,3), limits = c(0.6, 3.8)) +
  stat_cor(aes(color = LangGroup), label.x = 13, label.y = c(1.8, 3.6), show.legend = FALSE) +
  labs(x = "\nEducation", y = "") +
  theme(text = element_text(family = "Segoe UI"),
        legend.text = element_text(size = 10), legend.position = 'top',
        legend.background=element_blank(), legend.title = element_blank(),
        axis.title.x = element_text(colour="#08306B", size = 12, face = 2),
        #axis.title.y = element_text(colour="#08306B", size = 12, face = 2),
        axis.text = element_text(colour="#08306B", size = 10),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed')
  )
  +  coord_fixed(ratio = 2.5/1)

# correlation plot of reflected, log transformed, reflected MMSE and PC score using matched sample
cor_PCA_match <- ggscatter(match_trans, x = "scoresPCA", y = "ln_MMSE",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0, label.y = 3.7),
          xlab = "\nPrincipal Component Score", #ylab = "MMSE\n(Log Transformed)\n",
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  scale_x_continuous(breaks=c(0,1,2,3,4)) + scale_y_continuous(breaks=c(1,2,3), limits = c(0.6, 3.8)) +
  stat_cor(aes(color = LangGroup), label.x = 1, label.y = c(2.0, 3.25), show.legend = FALSE) +
  labs(x = "\nPrincipal Component Score", y = "") +
  theme(text = element_text(family = "Segoe UI"),
        legend.text = element_text(size = 10), legend.position = 'top',
        legend.background=element_blank(), legend.title = element_blank(),
        axis.title.x = element_text(colour="#08306B", size = 12, face = 2),
        #axis.title.y = element_text(colour="#08306B", size = 12, face = 2),
        axis.text = element_text(colour="#08306B", size = 10),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed')
  )





#scatter plot of correlation for original MMSE and Education by group
cor_edu_full <- ggscatter(swapcompletePCA, x = "Education", y = "MMSE", #title = "MMSE and Education Correlation by Language Group",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 19, label.y = 31),
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter", alpha = 0.7, size = 2
          ) +
  stat_cor(aes(color = LangGroup), label.x = 13, label.y = c(25, 31), show.legend = FALSE) +
  scale_x_continuous(breaks=c(12,14,16,18,20,22)) +
  scale_y_continuous(limits = c(9, 31)) +
  labs(x = "\nEducation", y = "") +
  theme(text = element_text(family = "Segoe UI"),
        legend.text = element_text(size = 10), legend.position = 'top',
        legend.background=element_blank(), legend.title = element_blank(),
        axis.title.x = element_text(colour="#08306B", size = 12, face = 2),
        axis.title.y = element_text(colour="#08306B", size = 12, face = 2),
        axis.text = element_text(colour="#08306B", size = 10),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed')
        ) + coord_fixed(ratio = 1/2)


#scatter plot of correlation for MMSE and Education by group using matched sample
ggscatter(df.matchPCA, x = "Education", y = "MMSE", title = "MMSE and Education Correlation by Language Group",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 19, label.y = 31),
          xlab = "\nEducation", ylab = "MMSE\n",
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  stat_cor(aes(color = LangGroup), label.x = 13, label.y = c(25, 31), show.legend = FALSE) +
  scale_x_continuous(breaks=c(12,14,16,18,20,22)) +
  theme(legend.position = c(0.90, 0.75), legend.background=element_blank(), legend.title = element_blank())


#scatter plot of correlation for MMSE and PCA using full sample
ggscatter(swapcompletePCA, x = "scoresPCA", y = "MMSE", title = "MMSE and Principal Component Correlation by Language Group",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 2.5, label.y = 28),
          xlab = "\nPrincipal Component Score", ylab = "MMSE\n",
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  stat_cor(aes(color = LangGroup), label.x = 2.5, label.y = c(22, 31), show.legend = FALSE) +
  theme(legend.position = c(0.80, 0.30), legend.background=element_blank(), legend.title = element_blank())

#scatter plot of correlation for MMSE and PCA using matched sample
ggscatter(df.matchPCA, x = "scoresPCA", y = "MMSE", title = "MMSE and Principal Component Correlation by Matched Language Groups",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 3, label.y = 27.5),
          xlab = "\nPrincipal Component Score", ylab = "MMSE\n",
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  stat_cor(aes(color = LangGroup), label.x = 3, label.y = c(20, 31), show.legend = FALSE) +
  theme(legend.position = c(0.80, 0.30), legend.background=element_blank(), legend.title = element_blank())




# correlation plot of reflected, log transformed, reflected MMSE and PC score using full sample
ggscatter(complete_trans, x = "scoresPCA", y = "ln_MMSE", title = "MMSE and Principal Component Correlation by Language Group",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = -3, label.y = 3.7),
          xlab = "\nPrincipal Component Score", ylab = "MMSE\n(Log Transformed)\n",
          color = "LangGroup", palette = c("#FC4E07", "#00AFBB"), position = "jitter") +
  stat_cor(aes(color = LangGroup), label.x = 1, label.y = c(2.24, 3.25), show.legend = FALSE) +
  scale_y_continuous(breaks=c(1.0,1.5,2.0,2.5,3.0,3.5)) +
  theme(legend.position = c(0.80, 0.30), legend.background=element_blank(), legend.title = element_blank())





#scatter plot of correlation for transformed MMSE and FA by group using matched sample
cor_FA_match <- ggscatter(match_trans, x = "FA", y = "ln_MMSE",
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.47, label.y = 3.6),
                          xlab = "\nFractional Anisotropy", ylab = "MMSE\n(Log Transformed)\n",
                          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  scale_x_continuous(breaks=c(0.40,0.42,0.44,0.46, 0.48,0.50)) + scale_y_continuous(breaks=c(1,2,3), limits = c(0.6, 3.8)) +
  stat_cor(aes(color = LangGroup), label.x = 0.42, label.y = c(2.0, 3.6), show.legend = FALSE) +
  labs(x = "\nEducation", y = "") +
  theme(text = element_text(family = "Segoe UI"),
        legend.text = element_text(size = 10), legend.position = 'top',
        legend.background=element_blank(), legend.title = element_blank(),
        axis.title.x = element_text(colour="#08306B", size = 12, face = 2),
        axis.title.y = element_text(colour="#08306B", size = 12, face = 2),
        axis.text = element_text(colour="#08306B", size = 10),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed')
  )


#scatter plot of correlation for MMSE and AD by group using matched sample
ggscatter(match_trans, x = "AD", y = "ln_MMSE", title = "MMSE and Axial Diffusivity Correlation by Language Group",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.00127, label.y = 3.6),
          xlab = "\nAxial Diffusivity", ylab = "MMSE\n(Log Transformed)\n",
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  stat_cor(aes(color = LangGroup), label.x = 0.00127, label.y = c(2.75, 3.2), show.legend = FALSE) +
  scale_y_continuous(breaks=c(1.0,1.5,2.0,2.5,3.0,3.5)) +
  theme(legend.position = c(0.50, 0.25), legend.background=element_blank(), legend.title = element_blank())

#scatter plot of correlation for MMSE and RD by group using matched sample
ggscatter(match_trans, x = "RD", y = "ln_MMSE", title = "MMSE and Radial Diffusivity Correlation by Language Group",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.0007, label.y = 3.6),
          xlab = "\nRadial Diffusivity", ylab = "MMSE\n(Log Transformed)\n",
          color = "LangGroup", palette = c("indianred2", "cornflowerblue"), position = "jitter") +
  stat_cor(aes(color = LangGroup), label.x = 0.0007, label.y = c(2.6, 3.0), show.legend = FALSE) +
  scale_y_continuous(breaks=c(1.0,1.5,2.0,2.5,3.0,3.5)) +
  theme(legend.position = c(0.50, 0.25), legend.background=element_blank(), legend.title = element_blank())





