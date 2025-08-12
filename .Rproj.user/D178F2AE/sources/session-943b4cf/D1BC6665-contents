rm(list=ls())
gc()
graphics.off()

library(readxl)
library(lme4)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(corrplot)
library(car)

Combined_JRC_Data <- read_xlsx ("Combined_Colony_JRC_Data.xlsx")

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
Clean_data_noNA <- Clean_data[!is.na(Clean_data$Prop_Active_Nests),]

Clean_data_noNA <- Clean_data_noNA %>%
  filter(!Survey_Year_yyyy %in% c(2016, 2017))
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


#Only including colonies that were involved in this programme.
CommunityColonyIDs <- c(
  "FA1", "BE1", "BE2", "BE3", "SK1", "SK2", "SK2a",
  "NB1", "NB2a", "NB3", "MG1a", "MG1b", "NJ2", "NJ1",
  "WA1", "WA2", "BA1", "BA1a"
)

CommunityColonies <- Clean_data_noNA[Clean_data_noNA$ColonyID %in% CommunityColonyIDs, ]

#Splitting into before and after 2016-2017
Before <- CommunityColonies[CommunityColonies$Survey_Year_yyyy <= 9,]
After <- CommunityColonies[CommunityColonies$Survey_Year_yyyy >=12,]
During <- CommunityColonies[CommunityColonies$Survey_Year_yyyy %in% c(10, 11),]


CommunityColonies<- CommunityColonies %>%
  mutate(BeforeOrAfter = case_when(
    Survey_Year_yyyy < 10 ~ "0",
    Survey_Year_yyyy >= 10 ~ "1"
  ))

#------------------------------------------------------------------------------------

##Before models
#Proportion
Before_prop <- glmer(cbind(TOTAL_No_Active_Nests, TOTAL_No_Inactive_Nests) ~ Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity + Survey_Year_yyyy + Protected_Area_Category  + (1 |ColonyID), family=binomial, data=Before)

#Active
Before_Active_P <- glmer(TOTAL_No_Active_Nests ~ Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity + Survey_Year_yyyy + Protected_Area_Category  + (1 | ColonyID), data=Before, family=poisson)


summary(Before_prop)
summary(Before_Active_P)


#-------------------------------------------------------


##After models

#Proportion

After_prop <- glm(cbind(TOTAL_No_Active_Nests, TOTAL_No_Inactive_Nests) ~ Log_Distance_to_Forest_Edge + prop_degraded_500m_s  + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Survey_Year_yyyy + Evidence_of_human_activity + Protected_Area_Category, family=binomial, data=After)


#Active
After_Active_P <- glm(TOTAL_No_Active_Nests ~ Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity + Survey_Year_yyyy + Protected_Area_Category, data=After, family=poisson)


summary(After_prop)
summary(After_Active_P)


#----------------------------------------------------------------------------------

before_after <- rbind(Before,After)
before_after$Intervention <- c(rep("Before",nrow(Before)),rep("After",nrow(After)))
before_after$Intervention <- factor(before_after$Intervention,levels=c("Before","After"))

#----------------------------------------------------------------------------------

BAProp<- glmer(cbind(TOTAL_No_Active_Nests, TOTAL_No_Inactive_Nests) ~ Survey_Year_yyyy * Intervention + Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity + Protected_Area_Category + (1 |ColonyID), family=binomial, data=before_after)

BAActive <- glmer(TOTAL_No_Active_Nests ~ Survey_Year_yyyy * Intervention + Log_Distance_to_Forest_Edge + prop_degraded_500m_s + Distance_to_Road_s + Distance_to_SRiver_s + Log_Distance_to_Settlement_s + Evidence_of_human_activity + Protected_Area_Category  + (1 |ColonyID), data=before_after, family=poisson)


summary(BAProp)
summary(BAActive)


BAPRes<- simulateResiduals(fittedModel = BAProp)
BAARes<- simulateResiduals(fittedModel = BAActive)

plot(BAPRes)
plot(BAARes)

