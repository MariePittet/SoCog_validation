### Psychometric validation of MASC in patients with acquired brain injury. CHUV, Lausanne, 2023 ###
### Main author: Marie Pittet
### PIs: Dr. Peggy d'Honincthun, Prof. Arseny A. Sokolov
### Contributors: Margot Picard, Anaïs Piolet, Magali Descombes

# Environment -------------------------------------------------------------

#install (uncomment) and load necessary packages
#install.packages("readxl")
library("readxl")
#install.packages("car")
library("car")
#install.packages("plyr")
library("plyr")
#install.packages("ggpubr")
library("ggpubr")
#install.packages("outliers")
library("outliers")
#install.packages("pROC")
library("pROC")
#install.packages("ggplot2")
library("ggplot2")
#install.packages("rstatix")
library("rstatix")
#install.packages("psych")
library("psych")
#install.packages("irr")
library("irr")
#install.packages("Hmisc")
library("Hmisc")
#install.packages("ggdist")
library("ggdist")
#install.packages("corrplot")
library("corrplot")
#install.packages("finalfit")
library("finalfit")

# Data wrangling ----------------------------------------------------------

#load data
df <- read_excel("SC_results.xlsx")

#indicate factor variables and name levels if needed.
df$Lesion <- factor(df$Lesion,
                    levels = c(0, 1),
                    labels = c("Controls", "ABI"))
df$SCL <- factor(df$SCL)
df$Sex <- factor(df$Sex,
                 levels = c(0, 1),
                 labels = c("Females", "Males"))
df$Diagnostic <- factor(df$Diagnostic)


# Demographic Table -------------------------------------------------------

explanatory <- c("Age", "Sex", "SCL")
df %>% summary_factorlist("Lesion", explanatory, p=TRUE, na_include=FALSE)

# Custom functions (for efficiency) ---------------------------------------

# Function 1: arranging data frames to perform t-tests or ANOVAS
arrange_df <- function(testname) {
  vec1 <- c()#initialize empty vector to store patients scores
  vec2 <- c()#initialize empty vector to store controls scores
  
  for (i in 1:nrow(df)) {
    #loop to store test scores according to lesion status (healthy vs brain-injured)
    if (df$Lesion[i] == "ABI") {
      vec1 <- append(vec1, testname[i])
    } else {
      vec2 <- append(vec2, testname[i])
    }
  }
  return(list(patients_scores = vec1, controls_scores = vec2))
}

# Function 2: identifying extreme values (Q1-3IQR / Q3 + 3IQR). Removing extremes
flag_extreme_values <-
  function(dataframe,
           independent_var1,
           independent_var2,
           dependent_var) {
    outliers <-  dataframe %>%
      group_by({
        {
          independent_var1
        }
      }, {
        {
          independent_var2
        }
      }) %>%
      identify_outliers({
        {
          dependent_var
        }
      })
    extremes <-  outliers %>%
      filter(outliers$is.extreme == TRUE)
    extremes <- as.data.frame(extremes)
    cleaned_df <- dataframe
    for (i in 1:nrow(extremes)) {
      cleaned_df <-
        cleaned_df[!(cleaned_df$ID == extremes$ID[i] &&
                       clean_df$error_type == extremes$error_type[i]), ]
    }
    
    return(extremes)
  }


#------------------------------------ 1. Difference of MASC total scores between ABI patients and controls ----------------------------------------

# arranging a dataframe with function 1
MASC <- arrange_df(df$MASCtot)

# Verifying t-tests assumptions
hist(MASC$patients_scores)
hist(MASC$controls_scores)

shapiro.test(MASC$patients_scores) #if p <.05, the distribution significantly differs from a normal distribution
shapiro.test(MASC$controls_scores) # same

var.test(MASC$patients_scores, MASC$controls_scores) #if p-value <.05, variances are different between the groups

# t-test for independent samples. If previous var.test SIG --> var.equal = FALSE
MASC_ttest <- t.test(MASC$patients_scores, MASC$controls_scores, var.equal = TRUE)
MASC_ttest

# Cute plots
library(RColorBrewer)
compare_means(MASCtot ~ Lesion, data = df)

plot_MASC1 <- ggplot(df, aes(x = Lesion, y = MASCtot, fill = Lesion))+
              ggdist::stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA, alpha = .5)+
              geom_boxplot (width = .15, outlier.color = NA)+
              ggdist::stat_dots(side = "left", justification = 1.12, binwidth = .25, alpha = .5, size = 2)+
              coord_cartesian(xlim = c(1,2.2))+
              #ggtitle("Main effect of group on MASC total scores") +
              xlab("Group")+
              ylab("MASC total score")+ 
              theme_pubclean()+
              #scale_fill_brewer(palette = "Set2")+
              scale_fill_manual(values = c("lightskyblue3","salmon2")) +
              theme(legend.position = "none")+
              #theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
              theme(axis.title = element_text(face = "bold", size = 18))+
              theme(axis.text = element_text(size = 16))+
              scale_y_continuous(expand = c(0, 0), limits = c(20, 45))+
              stat_compare_means(label =  "p.signif", label.x = 1.5,label.y = 42, method="t.test")
  
plot_MASC1  

#------------------------------------- 2. Difference of MASC control items scores between ABI patients and controls ---------------------------------------

# Again, arranging a dataframe to perform mean comparison with function 1
MASC_ctrl <- arrange_df(df$MASCcontrol)

# Verifying t-tests assumptions
hist(MASC_ctrl$patients_scores)
hist(MASC_ctrl$controls_scores)

shapiro.test(MASC_ctrl$patients_scores) #if p <.05, the distribution significantly differs from a normal distribution
shapiro.test(MASC_ctrl$controls_scores) # same
var.test(MASC_ctrl$patients_scores, MASC_ctrl$controls_scores) #if p-value <.05, variances are different between the groups

#Non parametric test for independent samples (as distributions are unlikely to be normal on control items)
MASC_ctrl_kruskal <- kruskal.test(df$MASCcontrol, df$Lesion)
MASC_ctrl_kruskal

# cute plot
compare_means(MASCcontrol ~ Lesion, data = df)
plot_MASC_ctrl1 <- ggplot(df, aes(x = Lesion, y = MASCcontrol, fill = Lesion))+
  ggdist::stat_halfeye(adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA, alpha = .5)+
  geom_boxplot (width = .12, outlier.color = NA)+
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = .07, alpha = .5)+
  coord_cartesian(xlim = c(1.2,1.8))+
  ggtitle("Main effect of group on MASC control questions") +
  xlab("Group")+
  ylab("MASC control questions score")+ 
  theme_pubclean()+
  scale_fill_manual(values = c("lightskyblue3","salmon2")) +
  #scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "none")+
  #theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
  theme(axis.title = element_text(face = "bold", size = 18))+
  theme(axis.text = element_text(size = 16))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7))+
  stat_compare_means(label.x = 1.4,label.y = 6.5, method="kruskal.test")

plot_MASC_ctrl1  

#------------------------------------------------------- 3. Exploring MASC error type by group ------------------------------------------------------

#Arranging a new data frame to analyze error type by group. This must be the least elegant way to do it, but hey, it works.
data_rows <- nrow(df)
error_type <- c(rep(1:3, times = data_rows))
error_score <- c()

for (i in 1:nrow(df)) {
  error_score <-
    c(error_score, df$MASChyper[i], df$MASChypo[i], df$MASCabsent[i])
}

lesionstatus1 <- c()
lesionstatus0 <- c()

for (i in 1:nrow(df)) {
  if (df$Lesion[i] == "ABI") {
    status = 1
    lesionstatus1 <-
      append(lesionstatus1, status)
  } else {
    status = 0
    lesionstatus0 <-
      append(lesionstatus0, status)
  }
}

lesion_status <-
  c(
    lesionstatus1,
    lesionstatus1,
    lesionstatus1,
    lesionstatus0,
    lesionstatus0,
    lesionstatus0
  )

ID <- c()
for (i in 1:nrow(df)) {
  ID <- append(ID, rep(df$ID[i], 3))
}

df_errortype <- data.frame(ID, lesion_status, error_type, error_score)
df_errortype$lesion_status <-
  factor(
    df_errortype$lesion_status,
    levels = c(0, 1),
    labels = c("Controls", "ABI")
  )
df_errortype$error_type <-
  factor(
    df_errortype$error_type,
    levels = c(1, 2, 3),
    labels = c("HyperToM", "HypoToM", "absent ToM")
  )

#Visualizing data
plot1 <-ggboxplot(df_errortype,
            x = "error_type",
            y = "error_score",
            color = "lesion_status")

#cleaning data
df_errortype_cleaned<- df_errortype[-49,]

# Summary statistics
df_errortype_cleaned %>%
  group_by(lesion_status, error_type) %>%
  get_summary_stats(error_score, type = "mean_sd")

# Assessment of normality (Shapiro-Wilk)
df_errortype_cleaned %>%
  group_by(lesion_status, error_type) %>%
  shapiro_test(error_score)
ggqqplot(df_errortype_cleaned, "error_score", ggtheme = theme_bw()) # qqplot to see if points follow reference line
facet_grid(lesion_status ~ error_type)

# Homogeneity of variances of within-subject variable (if p<.05, variances are not homogeneous across groups)
df_errortype_cleaned %>%
  group_by(error_type) %>%
  levene_test(error_score ~ lesion_status)

# Homogeneity of covariances of between-subject variable (if p<.001, covariances are not homogeneous. In that case, separate analyses and do multiple repeated measures ANOVAs).
box_m(df_errortype_cleaned[, "error_score", drop = FALSE], df_errortype_cleaned$lesion_status)

# Anova
require(rstatix)
fit <- anova_test(data = df_errortype_cleaned,
        dv = error_score,
        wid = ID,
        between = lesion_status,
        within = error_type)

summary(fit)
get_anova_table(fit)

# Contrasts
pwc <- df_errortype_cleaned %>%
  group_by(error_type) %>%
  pairwise_t_test(error_score ~ lesion_status, p.adjust.method = "bonferroni")
pwc

# Visualizing ANOVA results in a cute way
# cute plot
errors_anova_plot1 <- ggplot(df_errortype_cleaned, aes(x = lesion_status, y = error_score, fill = lesion_status))+
  geom_boxplot (width = .50, outlier.color = NA)+
  facet_wrap(~error_type)+
  #ggtitle("Error scores by error type by group") +
  xlab("Error type") +
  ylab("Error score") +
  theme_pubclean()+
  #scale_fill_brewer(palette = "Set2")+
  scale_fill_manual(values = c("lightskyblue3","salmon2")) +
  theme(legend.position = "none")+
  #theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
  theme(axis.title = element_text(face = "bold", size = 18))+
  theme(axis.text = element_text(size = 16))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 13))+
  stat_compare_means(label.y = 12, label = "p.signif", label.x = 1.5)

errors_anova_plot1


# Exploring affective/cognitive differences -------------------------------

#standardizing affective and cognitive scores
df$MASCaff<- df$MASCaff/15
df$MASCcog<- df$MASCcog/18

# Re-arranging data.frame for ANOVA
ID1<- c()
for (i in 1:nrow(df)){
  ID1<- append(ID1, rep(df$ID[i],2))
}

ToM_type<- c(rep(1:2, times=data_rows))

ToM_score<- c()
for (i in 1:nrow(df)){
  ToM_score<- c(ToM_score,df$MASCaff[i],df$MASCcog[i])
}

lesionstatus1<- c()
lesionstatus0<- c()
for (i in 1:nrow(df)){
  if (df$Lesion[i] == "ABI"){
    status = 1
    lesionstatus1<- append(lesionstatus1,status)
  } else {
    status = 0
    lesionstatus0<- append(lesionstatus0, status)
  }
}
lesion_status<- c(lesionstatus1, lesionstatus1,lesionstatus0,lesionstatus0)

df_ToMtype<- data.frame(ID1,lesion_status,ToM_type,ToM_score)
df_ToMtype$lesion_status<- factor(df_ToMtype$lesion_status, levels=c(0,1), labels=c("Controls", "ABI"))
df_ToMtype$lesion_status<- factor(df_ToMtype$lesion_status, levels=c("ABI","Controls"))
df_ToMtype$ToM_type<- factor(df_ToMtype$ToM_type, levels=c(1,2), labels=c("affective ToM", "cognitive ToM"))

fit1 <- anova_test(data = df_ToMtype,
                  dv = ToM_score,
                  wid = ID1,
                  between = lesion_status,
                  within = ToM_type)

summary(fit1)
get_anova_table(fit1)


# Visualizing data
ToMType_plot<-  ggboxplot(df_ToMtype, x = "ToM_type", y = "ToM_score", fill = "lesion_status", facet.by = "lesion_status", short.panel.labs = TRUE, lwd=1)+
  scale_fill_manual(values = c("darkorange2", "royalblue"))+
  xlab("ToM type")+
  ylab("ToM score")+
  ggtitle("ToM score by ToM type by group")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(axis.title = element_text(face="bold"))+
  theme(legend.position ="none")+
  stat_compare_means(label.x = 1.5)
ToMType_plot


# AUC ---------------------------------------------------------------------

#generating an error vector (MASC max points - participant'score)
MASC_errors<- c()
for (i in 1:nrow(df)){
  errors<- 45-df$MASCtot[i]
  MASC_errors<- append(MASC_errors,errors)
}

# generating an error vector (RMET)
RMET_errors<- c()
for (i in 1:nrow(df)){
  errors<- 33-df$RMETtot[i]
  RMET_errors<- append(RMET_errors,errors)
}

# generating an error vector (Attribution of intentions)
Sarfati_errors<- c()
for (i in 1:nrow(df)){
  errors<- 14-df$ToMSarfati[i]
  Sarfati_errors<- append(Sarfati_errors,errors)
}

# generating an error vector (emotion recognition)
EmoRecog_errors<- c()
for (i in 1:nrow(df)){
  errors<- 60-df$EmoRecogScore[i]
  EmoRecog_errors<- append(EmoRecog_errors,errors)
}

# generating an error vector (Faux-Pas)
FP_errors<- c()
for (i in 1:nrow(df)){
  errors<- 30-df$`Faux pasFP`[i]
  FP_errors<- append(FP_errors,errors)
}

# generating an error vector (Strange Stories)
StrangeStories_errors<- c()
for (i in 1:nrow(df)){
  errors<- 16-df$StrangeStoriesToM[i]
  StrangeStories_errors<- append(StrangeStories_errors,errors)
}

#generating a data frame for the ROC analysis
ROC_df<- data.frame(df$Lesion,MASC_errors, RMET_errors, Sarfati_errors, EmoRecog_errors, FP_errors, StrangeStories_errors)
ROC_df$df.Lesion <- factor(ROC_df$df.Lesion , levels=c("Controls", "ABI"))

#ROC analysis with Lesion status as classification criterion
roc1 = roc(ROC_df$df.Lesion,ROC_df$MASC_errors, smooth=FALSE)
roc2 = roc(ROC_df$df.Lesion,ROC_df$RMET_errors, smooth = FALSE)
roc3 = roc(ROC_df$df.Lesion,ROC_df$Sarfati_errors, smooth = FALSE)
roc4 = roc(ROC_df$df.Lesion,ROC_df$EmoRecog_errors, smooth = FALSE)
roc5 = roc(ROC_df$df.Lesion,ROC_df$FP_errors, smooth = FALSE)
roc6 = roc(ROC_df$df.Lesion, ROC_df$StrangeStories_errors, smooth = FALSE)

plot(roc1, col = "salmon", print.auc = TRUE, print.auc.y = 0.6, xlab = "1-Specificity", ylab = "Sensitivity")
plot(roc2, col = "lightsteelblue", add = TRUE, print.auc = TRUE, print.auc.y = 0.3)
plot(roc3, col = "navyblue", add = TRUE, print.auc = TRUE, print.auc.y = 0.1)
plot(roc4, col = "violet", add = TRUE, print.auc = TRUE, print.auc.y = 0.2)
plot(roc5, col = "gold", add = TRUE, print.auc = TRUE, print.auc.y = 0.5)
plot(roc6, col = "grey85", add = TRUE, print.auc = TRUE, print.auc.y = 0.4)
title(main = "Receiver Operator Curves", line = 2.5)
legend("bottomright",
       legend=c("MASC","Faux-Pas", "Strange Stories", "RMET", "Emotion Recognition", "Attribution of intentions"),
       col=c("salmon", "gold", "grey85", "lightsteelblue", "violet", "navyblue"),
       lwd=3, cex=0.6, xpd = TRUE, horiz = FALSE)


optimal_coords <- coords(roc1, "best", ret = c("threshold", "sensitivity", "specificity"))
ci.auc(roc1)

# 7. MASC convergent validity -------------------------------------------------------------
df_corr_all<- data.frame(df$MASCtot,df$MASCaff,df$MASCcog,df$EmoRecogScore,df$RMETtot,df$EmoFluency,df$ToMSarfati,df$`Faux pasFP`,df$LEAStot,df$StrangeStoriesToM)
df_corr_all.cor<- cor(df_corr_all, method = c("spearman"), use="complete.obs")
df_corr_all.rcorr<- rcorr(as.matrix(df_corr_all))
df_corr_all.p<- df_corr_all.rcorr$P
p_all<- df_corr_all.p[1:3,1:10]

colnames(df_corr_all)<- c('MASC total', 'MASC affective ToM', 'MASC cognitive ToM', "Emotion Recognition", "RMET", "Emotional fluency","Cartoon task", "Faux-Pas", "LEAS", "Strange Stories")

corrplot(df_corr_all.rcorr$r,type="upper", method = "number",  # full plot
         p.mat = df_corr_all.rcorr$P, sig.level = 0.05, insig = "blank")

corrplot(df_corr_all.cor[1:3,1:10], method = "number", p.mat = p_all, sig.level = 0.05, insig = "blank") # reduced plot


# Test-retest -------------------------------------------------------------

retest<- rcorr(df$MASCtot,df$REMASCtot, type = c("spearman"))
retest$P
retest$r

# Predictions of questionnaires -------------------------------------------

#OSCARS self ~ MASC total
lm1<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$MASCtot[!is.na(df$OSCARSself)])
sum_lm1<- summary(lm1)
plot(df$MASCtot[!is.na(df$OSCARSself)], df$OSCARSself[!is.na(df$OSCARSself)], xlab = c("MASC total"),ylab = c("OSCARS self"))
abline(lm1)
subtitle<- sprintf("p-value =  %s", sum_lm1$coefficients[2,4])
title(main = "OSCARS self ~ MASC", sub = subtitle)

lm1a<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$MASChypo[!is.na(df$OSCARSself)])
plot(df$MASChyper[!is.na(df$OSCARSself)], df$OSCARSself[!is.na(df$OSCARSself)], xlab = c("MASC HyperToM"),ylab = c("OSCARS self"))
abline(lm1a)
lm1b<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$MASChyper[!is.na(df$OSCARSself)])
plot(df$MASChypo[!is.na(df$OSCARSself)], df$OSCARSself[!is.na(df$OSCARSself)], xlab = c("MASC HypoToM"),ylab = c("OSCARS self"))
abline(lm1b)
lm1c<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$MASCabsent[!is.na(df$OSCARSself)])
plot(df$MASCabsent[!is.na(df$OSCARSself)], df$OSCARSself[!is.na(df$OSCARSself)], xlab = c("MASC absent ToM"),ylab = c("OSCARS self"))
abline(lm1c)


OSCARS_MASC<- ggplot(df, aes(x=MASCtot, y=OSCARSself)) +
  labs(x='MASC total score', y='OSCARS total scores', title='by total scores') +
  geom_point() +
  geom_smooth(method=lm , se=TRUE, col = "red") +
  geom_text(label = paste("Adj R² = ", round(summary(lm1)$r.squared, 2)), x = 35, y = 41) +
  geom_text(label = "p = ***", x = 35, y = 45) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face = "bold")) 

OSCARS_Hyper<- ggplot(df, aes(x=MASChyper, y=OSCARSself)) +
  labs(x='MASC hypermentalization errors', y='OSCARS total scores', title='by hypermentalization errors') +
  geom_point() +
  geom_smooth(method=lm , se=TRUE, col = "red") +
  geom_text(label = paste("Adj R² = ", round(summary(lm1a)$r.squared, 2)), x = 8, y = 41) +
  geom_text(label = "p = *", x = 8, y = 45) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face ="bold")) 


OSCARS_Hypo<- ggplot(df, aes(x=MASChypo, y=OSCARSself)) +
  labs(x='MASC hypomentalization errors', y='OSCARS total scores', title='by hypomentalization errors') +
  geom_point() +
  geom_smooth(method=lm , se=TRUE, col = "red") +
  geom_text(label = paste("Adj R² = ", round(summary(lm1b)$r.squared, 2)), x = 8, y = 41) +
  geom_text(label = "p = *", x = 8, y = 45) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face = "bold")) 


OSCARS_absent<- ggplot(df, aes(x=MASCabsent, y=OSCARSself)) +
  labs(x='MASC absent mentalization errors', y='OSCARS total scores', title='by absent mentalization errors') +
  geom_point() +
  geom_smooth(method=lm , se=TRUE, col = "red") +
  geom_text(label = paste("Adj R² = ", round(summary(lm1c)$r.squared, 2)), x = 4, y = 41) +
  geom_text(label = "p = **", x = 4, y = 45) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face = "bold")) 

OSCARS_errors<- ggarrange(OSCARS_MASC, OSCARS_Hyper, OSCARS_Hypo, OSCARS_absent,
                          ncol = 2, nrow = 2)
annotate_figure(OSCARS_errors, top = text_grob("Prediction of everyday life social cognition by MASC scores", 
                                      color = "red", face = "bold", size = 18))
OSCARS_MASC
OSCARS_Hyper
OSCARS_Hypo
OSCARS_absent


lm2<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$EmoRecogScore[!is.na(df$OSCARSself)])
summary(lm2)

OSCARS_Emorecog<- ggplot(df, aes(x=MASCtot, y=EmoRecogScore)) +
  labs(x=' Emotion recognition scores', y='OSCARS total scores', title='Prediction of OSCARS scores by emotion recognition scores') +
  geom_point() +
  geom_smooth(method=lm , se=TRUE, col = "red") +
  geom_text(label = paste("Adj R² = ", round(summary(lm2)$r.squared, 2)), x = 35, y = 43) +
  geom_text(label = "p = ***", x = 35, y = 45) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold')) 
OSCARS_Emorecog

lm3<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$RMETtot[!is.na(df$OSCARSself)])
summary(lm3)

lm4<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$ToMSarfati[!is.na(df$OSCARSself)])
summary(lm4)

lm5<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$EmoFluency[!is.na(df$OSCARSself)])
summary(lm5)

lm6<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$`Faux pasFP`[!is.na(df$OSCARSself)])
summary(lm6)

lm7<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$LEAStot[!is.na(df$OSCARSself)])
summary(lm7)

lm8<- lm(df$OSCARSself[!is.na(df$OSCARSself)] ~ df$StrangeStoriesToM[!is.na(df$OSCARSself)])
summary(lm8)

test_names<- c("MASC", "Emotion Recognition", "RMET", "Cartoon task", "Emotional fluency", "Faux-pas", "LEAS", "Strange Stories")
p_values<- c(summary(lm1)$coefficients[2,4], summary(lm2)$coefficients[2,4],summary(lm3)$coefficients[2,4],summary(lm4)$coefficients[2,4],summary(lm5)$coefficients[2,4],summary(lm6)$coefficients[2,4],summary(lm7)$coefficients[2,4],summary(lm8)$coefficients[2,4])
r2_values<- c(summary(lm1)$r.squared,summary(lm2)$r.squared,summary(lm3)$r.squared,summary(lm4)$r.squared,summary(lm5)$r.squared,summary(lm6)$r.squared,summary(lm7)$r.squared,summary(lm8)$r.squared)
df_predictive_power<- data.frame(test_names,p_values, r2_values)

# Visualization 

viz<- ggplot(df_predictive_power, aes(x = test_names, y = p_values, colors = test_names)) +
      geom_point(size = 3) +
      scale_y_continuous(trans = "log10") +
      ylab("p value (log10)") +
      ggtitle("Predictive value of social cognition tasks on OSCARS questionnaire") +
      geom_hline(yintercept=0.05, linetype="dashed", color = "red", linewidth=1)+
      geom_vline(xintercept = 1.5, color="gray") +
      geom_vline(xintercept = 2.5, color="gray") +
      geom_vline(xintercept = 3.5, color="gray") +
      geom_vline(xintercept = 4.5, color="gray") +
      geom_vline(xintercept = 5.5, color="gray") +
      geom_vline(xintercept = 6.5, color="gray") +
      geom_vline(xintercept = 7.5, color="gray") +
      theme_test()+
      theme(plot.title = element_text(hjust=0.5, size=15, face='bold')) +
      theme(axis.text.y = element_text(angle = 90))+
      theme(axis.text.x = element_text(size = 12, face = "bold"))+
      scale_x_discrete(guide=guide_axis(n.dodge=2)) +
      theme(axis.title.x = element_blank())

viz2<-  ggplot(df_predictive_power, aes(x = test_names, y = r2_values, colors = test_names)) +
        geom_point(size = 3) +
        ylab("Adjusted R squared") +
        xlab("test name") +
        geom_vline(xintercept = 1.5, color="gray") +
        geom_vline(xintercept = 2.5, color="gray") +
        geom_vline(xintercept = 3.5, color="gray") +
        geom_vline(xintercept = 4.5, color="gray") +
        geom_vline(xintercept = 5.5, color="gray") +
        geom_vline(xintercept = 6.5, color="gray") +
        geom_vline(xintercept = 7.5, color="gray") +
        theme_test() +
        theme(axis.title.x=element_blank()) +
        theme(axis.text.y = element_text(angle = 90)) +
        theme(axis.text.x = element_text(size = 12, face = "bold"))+
        scale_x_discrete(guide=guide_axis(n.dodge=2)) 

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(viz), ggplotGrob(viz2), size = "last"))


#OSCARS hetero ~ MASC total
lm2<- lm(df$OSCARShetero[!is.na(df$OSCARShetero)] ~ df$MASCtot[!is.na(df$OSCARShetero)])
sum_lm2<- summary(lm2)
plot(df$MASCtot[!is.na(df$OSCARShetero)], df$OSCARShetero[!is.na(df$OSCARShetero)], xlab = c("MASC total"),ylab = c("OSCARS hetero"))
abline(lm2)
subtitle<- sprintf("p-value =  %s", sum_lm2$coefficients[2,4])
title(main = "OSCARS hetero ~ MASC", sub = subtitle)


# Adjusting for neuropsychological covariates ----------------------------------------------------------

# correlations between MASC and npsy scores
df_MASC_npsy<- data.frame(df$MASCtot,df$TMTA,df$SDMT,df$LitFluency,df$CatFluency)
colnames(df_MASC_npsy)<- c('MASC total', 'TMT-A', 'SDMT', "Literal Fluency", "Categorial Fluency")
corr_MASC_npsy<- rcorr(as.matrix(df_MASC_npsy))
corrplot(corr_MASC_npsy$r,type="upper", method = "number",  # plot with significant correlation
         p.mat = corr_MASC_npsy$P, sig.level = 0.10, insig = "blank")


# ANCOVA (covariable comes first to suppress its effect)
#install.packages("emmean")
library(emmeans)
ancova_npsy<- anova_test(MASCtot ~ TMTA +  Lesion, data = df)
get_anova_table(ancova_npsy)
reg_mul<- lm(MASCtot ~ TMTA + SDMT + LitFluency + CatFluency + Lesion, data = df)
summary(reg_mul)

marginal_means<- emmeans_test(MASCtot ~ Lesion, covariate = TMTA, data = df)
get_emmeans(marginal_means)


# Item Response Theory ----------------------------------------------------
#install.packages("mirt")
library("mirt")

# data wrangling
irt_data<- read.csv("MASC_items.csv", header = TRUE, sep = ";")
irt_data %>% mutate(MASC_study = factor(MASC_study), Lesion_status = factor(Lesion_status))
irt_data1<- irt_data[c(1:40), -c(1,2,3,49)]

# exploring data dimensionality
library(psych)

efa_result <- fa(irt_data1, nfactors = 2, rotate = "oblimin", fm = "minres")
print(efa_result)

library(lavaan)

cfa_model <- '
F1 =~ Q6 + Q8 + Q11 + Q19 + Q31 + Q32 + Q33 + Q40 + Q41 + Q43 + Q44 + Q45
  F2 =~ Q3 + Q5 + Q10 + Q14 + Q18 + Q20 + Q36 + Q39 + Q40
'
fit <- cfa(cfa_model, data = irt_data1)
summary(fit, fit.measures = TRUE)

# fitting the 1-PL model (only item difficulty)
model <- mirt.model("socialcognition = 1 - 45")
model1PL <- mirt(irt_data1, model = 1, itemtype = "Rasch", SE = TRUE)

M2(model1PL) # model performance
residuals <- residuals(model1PL, type = "Q3") # checking for local independence
plot(residuals)
eigenvalues<- eigen(cor(irt_data1))$values # checking for unidimensionality
plot(eigenvalues)

bidimensional_model<- mirt(irt_data1, model = 2, itemtype = "Rasch")


itemfit(model1PL) # item fit (significant p values indicate poor fit)
residuals(model1PL, df.p = T) # checking for local independence assumption

model1PL_coefs<- coef(model1PL, simplify = TRUE, IRTpar = TRUE)$items # obtaining difficulty coefficients
mean_difficulty<- mean(model1PL_coefs[,2]) 
sd_difficulty<- sd(model1PL_coefs[,2]) # items are on average pretty easy Mean diff= -1.24, SD = 1.23


# Extract the plot object without displaying it
plot_obj <- plot(model1PL, type = "trace", facet_items = FALSE, 
                 auto.key = list(space = "right", cex = 0.7),  # Reduce legend text size
                 return = TRUE)

# Modify the plot object to change axis labels
plot_obj$xlab <- "ToM ability (θ)"
plot_obj$ylab <- "Probability of Correct Response"
plot_obj$main <- "Item Characteristic Curves"
plot_obj$strip <- FALSE
print(plot_obj)

plot(model1PL, type = "trace", facet_items = FALSE)
plot(model1PL, type = "infoSE") # overall test information function. Test is most informative in lower range of ability


# Fitting the 2-PL model (item difficulty and discrimination)
model2PL <- mirt(irt_data1, model, itemtype = "2PL", SE = TRUE)
M2(model2PL)
itemfit(model2PL) # item fit (significant p values indicate poor fit)
residuals(model2PL, df.p = T) # checking for local independence assumption

model2PL_coefs<- coef(model2PL, simplify = TRUE, IRTpar = TRUE)$items # obtaining difficulty and discrimination coefficients
mean_difficulty<- mean(model2PL_coefs[,2]) 
sd_difficulty<- sd(model2PL_coefs[,2]) # mean difficulty = -2.11, SD = 7.29
mean_discrimination<- mean(model2PL_coefs[,1]) 
sd_discrimination<- sd(model2PL_coefs[,1]) # discrimination = 0.73, SD = 1.33

plot(model2PL, type = "trace", facet_items = FALSE) # plotting ICCs
which(model2PL_coefs[,1] < 0) # showing which items have negative discrimination (which is bad news)
plot(model2PL, type = "infoSE") # overall test information function

#comparing model fit
anova(model1PL, model2PL) # most criteria favor model1PL (except for the HQ)


# Palying around (polytomous model) ---------------------------------------
poly_data<- read.csv("MASC_polytomous_data.csv", header = TRUE, sep = ";")
poly_data<- poly_data[-c(7,16),]

model <- mirt.model("socialcognition = 1 - 45")

model_poly <- mirt(poly_data, 1,"nominal")

summary(model_poly)

plot(model_poly, type = "trace", facet = TRUE)
plot(model_poly, type = "info", theta_lim = c(-4, 4))
M2(model_poly) # model performance

anova(model1PL, model_poly)


# Quick visualization of the effect of stroke vs TBI on scores -----------------------------------------------------------
p_total <- ggplot(df, aes(x=Diagnostic, y=MASCtot, fill = Diagnostic)) + 
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("lightskyblue3", "salmon1", "salmon3"))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill = "black")+
  labs(title="MASC total score", x ="Diagnostic", y = "MASC total score (45pts)")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face ="bold"))+
  theme(axis.title.x = element_blank())
p_total

p_hyper <- ggplot(df, aes(x=Diagnostic, y=MASChyper, fill = Diagnostic)) + 
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("lightskyblue3", "salmon1", "salmon3"))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill = "black")+
  labs(title="MASC hypermentalizations", x ="Diagnostic", y = "MASC n hyperToM")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  theme(axis.title.x = element_blank())
p_hyper

p_hypo <- ggplot(df, aes(x=Diagnostic, y=MASChypo, fill = Diagnostic)) + 
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("lightskyblue3", "salmon1", "salmon3"))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill = "black")+
  labs(title="MASC hypomentalizations", x ="Diagnostic", y = "MASC n hypoToM")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  theme(axis.title.x = element_blank()) 
p_hypo

p_noToM <- ggplot(df, aes(x=Diagnostic, y=MASCabsent, fill = Diagnostic)) + 
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("lightskyblue3", "salmon1", "salmon3"))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, fill = "black")+
  labs(title="MASC absences of mentalization", x ="Diagnostic", y = "MASC n noToM")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  theme(axis.title.x = element_blank())
p_noToM

grid.arrange(p_total, p_hyper, p_hypo, p_noToM, ncol=2)

