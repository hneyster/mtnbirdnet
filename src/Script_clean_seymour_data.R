# Mountain Bird Network 2021 analysis
# first draft to clean up seymour data and plot elevational patterns of species richness

library(rgbif)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic()) #set ggplot theme here


# read in csv with raw data from ebird download

mydata <- read.csv("2021_MBN_data.csv")

# add elevation to each row using "elevation" function from rgbif package
# takes ~ 2 minutes to extract elevations for each row
elevation <- elevation(latitude = mydata$Latitude,
          longitude = mydata$Longitude,
          elevation_model = "srtm1",
          username = "bgfreeman")

mydata <- cbind(mydata, elevation[,3])
rename(mydata, elevation = elevation[, 3])
mydata <- rename(mydata, c("elevation" = "elevation[, 3]"))

# this dataframe includes data for several elevational transects
# you can see that quickly by looking at the different counties (note some of the New Mexico sites span multiple counties)
table(mydata$County)

# for now just look at seymour data
seymour <- mydata[mydata$County == "Metro Vancouver",]


########
# clean up Seymour data
# first synonymize names
# second deal with individuals > 50 m using checklist comments
########

table(seymour$Common.Name)
# notes: 
# Dark-eyed Junco coded both as "Dark-eyed Junco" and "Dark-eyed Junco (Oregon)"
# Evening Grosbeak coded both as "Evening Grosbeak" and "Evening Grosbeak (type 1)"
# Red Crossbill coded both as "Red Crossbill" and "Red Crossbill (Western Hemlock or type 3)"
# also a couple "sp's"

# these are reflected in the scientific names as well, so need to synonymize

seymour$Common.Name[seymour$Common.Name == "Dark-eyed Junco (Oregon)"] <- "Dark-eyed Junco"
seymour$Common.Name[seymour$Common.Name == "Evening Grosbeak (type 1)"] <- "Evening Grosbeak"
seymour$Common.Name[seymour$Common.Name == "Red Crossbill (Western Hemlock or type 3)"] <- "Red Crossbill"


table(seymour$Observation.Details)

# quite a few with varying text to indicate all individuals were beyond 50 m
all_beyond_50 <- c("&gt; 50", "&gt;50m", "&gt; 50 (both)", "&gt; 50m", "&gt;50", "&gt;59", "both likely beyond 50m",
                   "both whooping beyond 50m", "more than 50 m", "whooping beyond 50 m", "whooping beyond 50m", "whooping. more than 50 m", "whopping beyond 50m")

# also text to indicate some individuals were beyond 50 m (all 1 individual > 50 m)
some_beyond_50 <- c("1 &gt; 50", "1 &gt;50", "1 &gt;50m", "1 in", "1 seen within 50m. another heard beyond 50 m", "1&gt; 50", "1&gt;50")

seymour$beyond_50 <- ifelse(seymour$Observation.Details %in% all_beyond_50, seymour$Count, 
                                  ifelse(seymour$Observation.Details %in% some_beyond_50, 1, 0))

seymour$within_50 <- seymour$Count - seymour$beyond_50

# write a csv with the cleaned data ready for analysis

write.csv(seymour, "seymour.csv")

#######
# make plot of raw species richness vs. elevation
######


seymour_pointcounts <-
  seymour %>%
  group_by(Submission.ID) %>%
  dplyr::summarise(sprichness = n_distinct(Common.Name), sprichness_within50 = n_distinct(Common.Name[within_50!=0]),
                     abundance_total = sum(Count), abundance_within = sum(within_50))


# add elevation to this

seymour_pointcounts$elevation <- seymour$elevation[match(seymour_pointcounts$Submission.ID, seymour$Submission.ID)]

# make plots of species richness vs. elevation

# for within 50 m
ggplot(seymour_pointcounts, aes(x = elevation, y = sprichness_within50))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = loess)+
  scale_x_continuous("elevation (m)", breaks = c(0, 250, 500, 750, 1000, 1250, 1500))+
  scale_y_continuous("species richness", breaks = c(0, 2, 4, 6, 8, 10, 12, 14))
  
# for overall
ggplot(seymour_pointcounts, aes(x = elevation, y = sprichness))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = loess)+
  scale_x_continuous("elevation (m)", breaks = c(0, 250, 500, 750, 1000, 1250, 1500))+
  scale_y_continuous("species richness", breaks = c(0, 2, 4, 6, 8, 10, 12, 14))
