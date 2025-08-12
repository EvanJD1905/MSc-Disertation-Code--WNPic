rm(list=ls())
gc()
graphics.off()


library(readxl)
library(lme4)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(corrplot)
library(performance)
library(ggeffects)


#------------------------------------------------------------------------------------------------------------------------
Combined_JRC_Data <- read_xlsx("Combined_Colony_JRC_Data.xlsx")

Clean_data <- Combined_JRC_Data[, c(
  "TOTAL_No_Active_Nests", "TOTAL_No_Inactive_Nests", "ColonyID", "UniqueID", 
  "Prop_Active_Nests", "Protected_Area_Category", "Distance_to_Settlement", "Distance_to_Road", 
  "Distance_to_Footpaths_and_Trails", "Distance_to_Small_River", 
  "Survey_Year_yyyy", "Evidence_of_human_activity", "dist_to_large_patch_edge_m", 
  "prop_degraded_500m", "prop_degraded_1km")]


for(colony in unique(Clean_data$ColonyID)){
  toto <- Clean_data[which(Clean_data$ColonyID==colony),]
  if(length(unique(toto$Protected_Area_Category))>1){
    # print(c(colony,unique(toto$Protected_Area_Category)))
    # print("")
    Clean_data[which(Clean_data$ColonyID==colony),"Protected_Area_Category"] <- unique(toto$Protected_Area_Category)[1]
  }
}
# tata <- Clean_data[which(Clean_data$Survey_Year_yyyy>=2023),]
# unique(tata$Protected_Area_Category)


# Remove rows where Prop_Active_Nests is NA
Clean_data$Prop_Active_Nests <- as.numeric(Clean_data$Prop_Active_Nests)
Clean_data_noNA <- Clean_data[!is.na(Clean_data$Prop_Active_Nests), ]

Clean_data_noNA$Survey_Year_yyyy <- Clean_data_noNA$Survey_Year_yyyy - 2006
Clean_data_noNA$ObsID <- factor(1:nrow(Clean_data_noNA))
Clean_data_noNA$Protected_Area_Category <- relevel(
  factor(Clean_data_noNA$Protected_Area_Category),
  ref = "Off Site"
)


#Roads
Clean_data_noNA$Distance_to_Road_s <- scale(Clean_data_noNA$Distance_to_Road)
Clean_data_noNA$Log_Distance_to_Road<- log(Clean_data_noNA$Distance_to_Road)
Clean_data_noNA$Log_Distance_to_Road_s <- scale(Clean_data_noNA$Log_Distance_to_Road)
#Footpaths & Trails
Clean_data_noNA$Distance_to_FandT_s <- scale(Clean_data_noNA$Distance_to_Footpaths_and_Trails)
Clean_data_noNA$Log_Distance_to_FandT<- log(Clean_data_noNA$Distance_to_Footpaths_and_Trails)
Clean_data_noNA$Log_Distance_to_FandT_s <- scale(Clean_data_noNA$Log_Distance_to_FandT)
#Rivers
Clean_data_noNA$Distance_to_SRiver_s <- scale(Clean_data_noNA$Distance_to_Small_River)
Clean_data_noNA$Log_Distance_to_SRiver<- log(Clean_data_noNA$Distance_to_Small_River)
Clean_data_noNA$Log_Distance_to_SRiver_s <- scale(Clean_data_noNA$Log_Distance_to_SRiver)
#Settlements
Clean_data_noNA$Distance_to_Settlement_s <- scale(Clean_data_noNA$Distance_to_Settlement)
Clean_data_noNA$Log_Distance_to_Settlement<- log(Clean_data_noNA$Distance_to_Settlement)
Clean_data_noNA$Log_Distance_to_Settlement_s <- scale(Clean_data_noNA$Log_Distance_to_Settlement)
#Distance to Forest Edge
Clean_data_noNA$Distance_to_Forest_Edge_s <- scale(Clean_data_noNA$dist_to_large_patch_edge_m)
Clean_data_noNA$Log_Distance_to_Forest_Edge<- log(Clean_data_noNA$dist_to_large_patch_edge_m)
Clean_data_noNA$Log_Distance_to_Forest_Edge_s <- scale(Clean_data_noNA$Log_Distance_to_Forest_Edge)
#Degradation 500m
mean_val <- mean(Clean_data_noNA$prop_degraded_1km, na.rm = TRUE)
sd_val <- sd(Clean_data_noNA$prop_degraded_1km, na.rm = TRUE)

Clean_data_noNA$prop_degraded_1km_s <- (Clean_data_noNA$prop_degraded_1km - mean_val) / sd_val

#Degradation 1km
mean_val <- mean(Clean_data_noNA$prop_degraded_500m, na.rm = TRUE)
sd_val <- sd(Clean_data_noNA$prop_degraded_500m, na.rm = TRUE)

Clean_data_noNA$prop_degraded_500m_s <- (Clean_data_noNA$prop_degraded_500m - mean_val) / sd_val


#------------------------------------------------------------------------------------------------------------------------

#Overdispersion
Overdisp <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type= "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail = FALSE)
  c(chisq= Pearson.chisq, ratio= prat, rdf = rdf, p = pval)
}

#------------------------------------------------------------------------------------------------------------------------

Clean_data_noNA$Evidence_of_human_activity <- as.numeric(as.factor(Clean_data_noNA$Evidence_of_human_activity))
Clean_data_noNA$Protected_Area_Category <- as.numeric(as.factor(Clean_data_noNA$Protected_Area_Category))

numeric_data <- Clean_data_noNA[, sapply(Clean_data_noNA, is.numeric)]

selected_predictors <- c("Log_Distance_to_Forest_Edge_s", 
                         "prop_degraded_500m", 
                         "Log_Distance_to_FandT_s",
                         "Distance_to_Road_s", 
                         "Log_Distance_to_SRiver_s", 
                         "Log_Distance_to_Settlement_s", 
                         "Evidence_of_human_activity", 
                         "Survey_Year_yyyy", 
                         "Protected_Area_Category")

selected_predictors <- intersect(selected_predictors, names(numeric_data))

cor_matrix <- cor(numeric_data[, selected_predictors], use = "complete.obs")

corrplot(cor_matrix, method = "color", addCoef.col = "black", tl.cex = 0.8)


#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------

#Final Models
#Proportion
Proportion <- glmer(cbind(TOTAL_No_Active_Nests, TOTAL_No_Inactive_Nests) ~ Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Log_Distance_to_FandT_s  + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity + Survey_Year_yyyy + Protected_Area_Category + (1|ColonyID), family=binomial, data=Clean_data_noNA)

#Active
Active_P <- glmer(TOTAL_No_Active_Nests ~ Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Log_Distance_to_FandT_s + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity  + Survey_Year_yyyy + Protected_Area_Category + (1 | ColonyID) + (1|UniqueID), data=Clean_data_noNA, family=poisson)

summary(Proportion)
summary(Active_P)