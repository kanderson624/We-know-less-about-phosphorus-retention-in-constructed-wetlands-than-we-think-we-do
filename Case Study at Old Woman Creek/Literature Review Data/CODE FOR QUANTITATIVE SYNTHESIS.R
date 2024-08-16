set.seed(1680) # for reproducibility
library(dplyr) # for data cleaning
library(ISLR) # for college dataset
library(cluster) # for gower similarity and pam
library(Rtsne) # for t-SNE plot
library(ggplot2) # for visualization
library(viridis)
library(lubridate)
library(dplyr)
library(tidyverse)
library(ggpubr)

########################
## Clustering
########################
# code from https://dpmartin42.github.io/posts/r/cluster-mixed-types
# Version with MY data
setwd("C:/Users/kande120/OneDrive - Kent State University/Meta-Analysis")
data<-read.csv("Final_Dataset_by_wetland.csv")
#data<-read.csv("groupingpapers.csv")
# remove 

data1 = subset(data, select = c(Wetland.name, System,Class,Regime,Artificial,Source,Landuse,Surface.area))
data1$System = as.factor(data1$System)
data1$Class = as.factor(data1$Class)
data1$Regime = as.factor(data1$Regime)
data1$Artificial = as.factor(data1$Artificial)
data1$Source = as.factor(data1$Source)
data1$Landuse = as.factor(data1$Landuse)
#data1$Age = as.numeric(data1$Age)
data1$Surface.area = as.numeric(data1$Surface.area)
#data1$cat_area = as.factor(data1$cat_area)


# Remove wetland name before clustering
gower_dist <- daisy(data1[, -1],
                    metric = "gower",
                    type = list(logratio = 3))
summary(gower_dist)


gower_mat <- as.matrix(gower_dist)
# Output most similar pair
data1[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]



# Output most dissimilar pair
data1[
  which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
        arr.ind = TRUE)[1, ], ]

# Calculate silhouette width for many k using PAM
sil_width <- c(NA)
for(i in 2:10){
  
  pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}
# Plot sihouette width (higher is better)
plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)

pam_fit <- pam(gower_dist, diss = TRUE, k = 2)
pam_results <- data1 %>%
  dplyr::select(-Wetland.name) %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary # this shows average statistics for each of the clusters

data1[pam_fit$medoids, ]
print(pam_fit)

# this makes a vector of the clusters
clusters3 = pam_fit$clustering
# this adds the clusters to the dataframe
dataclus = cbind(data1, clusters3)
# print out the csv with the clusters: 

pam_fit <- pam(gower_dist, diss = TRUE, k = 4)
pam_results <- data1 %>%
  dplyr::select(-Wetland.name) %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
pam_results$the_summary # this shows average statistics for each of the clusters

#data1[pam_fit$medoids, ]
#print(pam_fit)

# this makes a vector of the clusters
clusters8 = pam_fit$clustering
# this adds the clusters to the dataframe
dataclus = cbind(dataclus, clusters8)
# print out the csv with the clusters: 
#write.csv(dataclus,"clusters.csv", row.names = FALSE)



tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)
tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = data1$name)
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster))+theme_classic()+scale_colour_viridis(discrete=TRUE)

## Stats by grouping
setwd("C:/Users/kande120/OneDrive - Kent State University/Meta-Analysis")
data<-read.csv("Final_Dataset_by_wetland.csv")

data$clusters8 = as.factor(data$clusters8)

data %>% 
  group_by(clusters8) %>% 
  summarize(avg=mean(TPretention.g.m2.y,na.rm=TRUE), n=n(), sd=sd(TPretention.g.m2.y,na.rm=TRUE), se=sd/sqrt(n), median = median(TPretention.g.m2.y,na.rm=TRUE), max=max(TPretention.g.m2.y,na.rm=TRUE), min = min(TPretention.g.m2.y,na.rm=TRUE))

a1 = aov(data=data,TPretention.g.m2.y~clusters8 )
summary(a1)
TukeyHSD(a1)

data %>% 
  group_by(clusters8) %>% 
  summarize(avg=mean(TPloading.g.m2.y,na.rm=TRUE), n=n(), sd=sd(TPloading.g.m2.y,na.rm=TRUE), se=sd/sqrt(n), median = median(TPloading.g.m2.y,na.rm=TRUE), max=max(TPloading.g.m2.y,na.rm=TRUE), min = min(TPloading.g.m2.y,na.rm=TRUE))

a2 = aov(data=data,TPloading.g.m2.y~clusters8 )
summary(a2)
TukeyHSD(a2)

# Plot this out
preten <- ggplot(data, aes(x=clusters8, y=TPretention.g.m2.y)) + 
  geom_boxplot()+ylim(-200,1400)
preten + geom_jitter(shape=16, position=position_jitter(0.2))+theme_classic()
ploading <- ggplot(data, aes(x=clusters8, y=TPloading.g.m2.y)) + 
  geom_boxplot()
ploading + geom_jitter(shape=16, position=position_jitter(0.2))+theme_classic()

# log axis version
data$logTP = log(data$TPretention.g.m2.y)
data$logTP[is.infinite(data$logTP)] <- NA
preten <- ggplot(data, aes(x=clusters8, y=logTP)) + 
  geom_boxplot()+ylim(-5,7.5)
preten + geom_jitter(shape=16, position=position_jitter(0.2))+theme_classic()
# log10 axis version
data$log10TP = log10(data$TPretention.g.m2.y)
data$log10TP[is.infinite(data$log10TP)] <- NA
preten <- ggplot(data, aes(x=clusters8, y=log10TP)) + 
  geom_boxplot()+ylim(-2.5,3.5)
preten + geom_jitter(shape=16, position=position_jitter(0.2))+theme_classic()

########################
## Sampling Years Monitored
########################
setwd("C:/Users/kande120/OneDrive - Kent State University/Meta-Analysis")
data<-read.csv("Final_Dataset_by_wetland.csv")
# first take only a single row from each manuscript
data1 = data[!duplicated(data$Citation),]

# make it define years monitored as a character:
data1$Years.monitored = as.character(data1$Years.monitored)

df=data[c("Wetland.name","Years.monitored")]


# group wetlands by site (this takes out once where we're getting higher values for number of years monitored because we have multiple years as rows)
d1=df %>%                                       # Set up the pipe
  subset(complete.cases(df)) %>%               # Removes rows with NA values
  group_by(Wetland.name) %>%                           # Groups by the Name column
  count(Years.monitored)


# tell R to make years monitored back into a number: 
d1$Years.monitored=as.numeric(d1$Years.monitored)
# and tell R to round it to the nearest whole number
d1$Years.monitored=ceiling(d1$Years.monitored)
# make the plot
p1 = ggplot(d1, aes(x=Years.monitored)) + geom_histogram(binwidth=1)+theme_classic()+
  xlab("Years Monitored")+ylab("Number of Studies")

data1$Age..years.=as.numeric(data1$Age..years.)

p2 = ggplot(data1, aes(x=Age..years.)) + geom_histogram(binwidth=1)+theme_classic()+
  xlab("Age of wetland")+ylab("Number of Studies")
ggarrange(p1,p2)


########################
## Sampling Methods
########################
setwd("C:/Users/kande120/OneDrive - Kent State University/Meta-Analysis")
data<-read.csv("Final_Dataset_by_wetland.csv")
# first take only a single row from each manuscript
df = data[!duplicated(data$Citation),]

# order the factors
df$Load_method <- ordered(df$Load_method, levels = c("Constant", "event", "flow-weighted","Linear Interpolation","other", "ISCO"))
df$Frequency <- ordered(df$Frequency, levels = c("daily", "biweekly", "weekly", "bimonthly","monthly","seasonally","event"))


# remove NAs
dfwater = subset(df, Load_method != "NA")

ggplot(dfwater, aes(x = Load_method,fill=Frequency)) +
  geom_bar()+theme_classic()+
  theme(axis.text.x = element_text(angle=45))
# save pdf as portrait and 4.5 x 4.5 inches


dfwater$Load_method <- factor(dfwater$Load_method, levels = c("ISCO", "other", "Constant", "flow-weighted", "Linear Interpolation", "event"))

# alternate version: 
ggplot(dfwater, aes(x = Frequency,fill=Load_method)) +
  geom_bar()+theme_classic()+
  theme(axis.text.x = element_text(angle=45))

dfwater %>% 
  count(Frequency)
dfwater %>% 
  count(Load_method)


########################
## Sampling Frequency Demonstration with OWC data
########################
setwd("C:/Users/kande120/OneDrive - Kent State University/Phosphorus Budgets/OWC/Data")
owcmasterdata<-read.csv("Master_Nutirent_Dataset_3-7-24.csv")
owcmasterdata$Date<-mdy(owcmasterdata$Date) # tell R that the date is a date 

# now lets subset the hell out of this dataset:
owcRL <- subset(owcmasterdata, Project == "RL") # only the retention and load project 
owcRLinflow <- subset(owcRL, Site == "BR") # only the inflow at Berlin Road
owcRLinflow23 <- subset(owcRLinflow, Date >= "2023-01-01" & Date <= "2023-02-10") # only from 2015 on
owc23 = owcRLinflow23 %>% group_by(Site, Date = as.Date(Date,"%Y-%m-%d")) %>% summarise(across(c(TP), mean)) # average any daily TP replicates together 

# make monthly and weekly points
owc23$month<-month(owc23$Date) # Create a variable for each month
owc23$daynum<-day(owc23$Date) # Create a variable for each day of the month
owc23$week<-isoweek(owc23$Date) # Take out the number for each isoweek.
owc23$day <- weekdays(owc23$Date) # Take out the day of the week


owc23month = subset(owc23, daynum == 6)
owc23week = subset(owc23, day == "Tuesday")
owc23event = subset(owc23, TP >= 400)

ggplot(owc23, aes(x=Date, y=TP)) +
  geom_line(color = "#440154FF") + 
  geom_point(data=owc23event, color = "#FDE725FF", size = 3)+
  geom_point(data=owc23week, color = "#31688EFF", size = 3)+
  geom_point(data=owc23month, color = "#35B779FF", size = 3)+
   xlab("") + ylab("TP (ug/L)")+theme_classic()


########################
## OWC Discharge - TP relationship
########################
setwd("C:/Users/kande120/OneDrive - Kent State University/Phosphorus Budgets/OWC/Data")
owcmasterdata<-read.csv("Master_Nutirent_Dataset_3-7-24.csv")
owcmasterdata$Date<-mdy(owcmasterdata$Date) # tell R that the date is a date 

# now lets subset the hell out of this dataset:
owcRL <- subset(owcmasterdata, Project == "RL") # only the retention and load project 
owcRLinflow <- subset(owcRL, Site == "BR") # only the inflow at Berlin Road
owcRLinflowrecent <- subset(owcRLinflow, Date >= "2015-01-01") # only from 2015 on
owcDF = owcRLinflowrecent %>% group_by(Site, Date = as.Date(Date,"%Y-%m-%d")) %>% summarise(across(c(TP), mean)) # average any daily TP replicates together 
library(dataRetrieval)
# Old Woman Creek at Berlin Rd near Huron OH
siteNumber <- "04199155"
owcInfo <- readNWISsite(siteNumber)
parameterCd <- "00060"

# As a note they have sub-daily data, but we don't need it. We're going to be pairing with daily nutrient data so it's actually super convenient that this code makes the data daily for us!
# Raw daily data: 
rawDailyData <- readNWISdv(
  siteNumber, parameterCd,
  "2015-05-05", "2024-01-01")
colnames(rawDailyData)[colnames(rawDailyData) == 'X_00060_00003'] <- 'discharge_cfs'  #rename the discharge column to something meaningful

owcDF <- owcDF %>% 
  left_join(y=rawDailyData, by=c("Date")) # merge out 2 datasets

owcBR = owcDF[c("Date","TP","discharge_cfs")] # clean up our dataset so it's just our 3 variables we're interested in
owcBR$discharge_L = owcBR$discharge_cfs * 86400 * 28.3168 # Convert cubic feet/second to liters/day
owcBR$TPmgL = owcBR$TP / 1000 # convert ug/L to mg/L
owcBR$P_lbs = owcBR$discharge_L * owcBR$TPmgL / 453600 # calculate P load and convert to lbs 


# we're going to gap fill. It'll make everything easier in terms of coding this out
# first we need to fill in all missing dates with NAs so r doesn't get mad at us:
owcBRcomplete = owcBR %>% 
  complete(Date = seq.Date(as.Date("2015-05-04"), max(Date), by = "1 day"))

# Now we need to define the weeks and years in our dataset so r can pick them out
owcBRcomplete$day <- weekdays(owcBRcomplete$Date) # Take out the day of the week
owcBRcomplete$week<-isoweek(owcBRcomplete$Date) # Take out the number for each isoweek.
owcBRcomplete$year<-isoyear(owcBRcomplete$Date) # Create a variable for each isoyear

# then we fill all the gaps with the median value (we're going to use the median for the ENTIRE DATASET)
owcBR_GF = owcBRcomplete %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

QClm = lm(TP~discharge_L, data = owcBR_GF)
summary(QClm)

d1 = owcBR_GF %>% 
  group_by(year)%>%
  summarize(meanP = mean(TPmgL,na.rm = T),medP = median(TPmgL,na.rm = T))
mean(d1$medP)
########################
## OWC Discharge - TP relationship at the outlet
########################
setwd("C:/Users/kande120/OneDrive - Kent State University/Phosphorus Budgets/OWC/Data")
owcmasterdata<-read.csv("Master_Nutirent_Dataset_3-7-24.csv")
owcmasterdata$Date<-mdy(owcmasterdata$Date) # tell R that the date is a date 

# now lets subset the hell out of this dataset:
owcRL <- subset(owcmasterdata, Project == "RL") # only the retention and load project 
owcRLoutflow <- subset(owcRL, Site == "WM") # only the outflow at the wetland mouth
owcRLoutflowrecent <- subset(owcRLoutflow, Date >= "2015-01-01" & Date < "2019-01-01") # only from 2015 on
owcoDF = owcRLoutflowrecent %>% group_by(Site, Date = as.Date(Date,"%m/%d/%Y")) %>% summarise(across(c(TP), mean)) # average any daily TP replicates together 
library(dataRetrieval)
# Old Woman Creek at Berlin Rd near Huron OH
siteNumber <- "04199155"
owcInfo <- readNWISsite(siteNumber)
parameterCd <- "00060"

# As a note they have sub-daily data, but we don't need it. We're going to be pairing with daily nutrient data so it's actually super convenient that this code makes the data daily for us!
# Raw daily data: 
rawDailyData <- readNWISdv(
  siteNumber, parameterCd,
  "2015-05-05", "2024-01-01")
colnames(rawDailyData)[colnames(rawDailyData) == 'X_00060_00003'] <- 'discharge_cfs'  #rename the discharge column to something meaningful

owcoDF <- owcoDF %>% 
  left_join(y=rawDailyData, by=c("Date")) # merge out 2 datasets

owcoBR = owcoDF[c("Date","TP","discharge_cfs")] # clean up our dataset so it's just our 3 variables we're interested in
owcoBR$discharge_L = owcoBR$discharge_cfs * 86400 * 28.3168 # Convert cubic feet/second to liters/day
owcoBR$TPmgL = owcoBR$TP / 1000 # convert ug/L to mg/L
owcoBR$P_lbs = owcoBR$discharge_L * owcoBR$TPmgL / 453600 # calculate P load and convert to lbs 


# we're going to gap fill. It'll make everything easier in terms of coding this out
# first we need to fill in all missing dates with NAs so r doesn't get mad at us:
owcoBRcomplete = owcoBR %>% 
  complete(Date = seq.Date(as.Date("2015-05-04"), max(Date), by = "1 day"))

# Now we need to define the weeks and years in our dataset so r can pick them out
owcoBRcomplete$day <- weekdays(owcoBRcomplete$Date) # Take out the day of the week
owcoBRcomplete$week<-isoweek(owcoBRcomplete$Date) # Take out the number for each isoweek.
owcoBRcomplete$year<-isoyear(owcoBRcomplete$Date) # Create a variable for each isoyear

# then we fill all the gaps with the median value (we're going to use the median for the ENTIRE DATASET)
owcoBR_GF = owcoBRcomplete %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))

QClm = lm(TP~discharge_L, data = owcoBR_GF)
summary(QClm)

d1 = owcoBR_GF %>% 
  group_by(year)%>%
  summarize(meanP = mean(TPmgL,na.rm = T),medP = median(TPmgL,na.rm = T))
mean(d1$medP)
