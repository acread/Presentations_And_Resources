This is a little code block to make a map color coded by count from each country \
I did it in a hacky way by first pulling the full list of countries and adding a \
second column with the count per country (this was really easy in this case because \
there was a relatively low number of countries.

The scale here is not ideal because it is discrete and does not include all values from 0-max (9) \


`````R
library(ggplot2)
library(dplyr)
require(maps)
require(viridis)
library("MetBrewer")

setwd("/Users/read0094/Desktop")
Speaker_Inst_Country=fread("MPMI_speaker_list_simple.txt")
#I madet this .txt file by pulling the list of countries and appending a
# "SpeakerCount" column

#head(Speaker_Inst_Country)
#region SpeakerCount
#1:       Aruba            0
#2: Afghanistan            0
#3:      Angola            0
#4:    Anguilla            0
#5:     Albania            0
#6:     Finland            0

world_map <- map_data("world")
#world_regions=as.data.frame(unique(world_map$region))
#write_tsv(world_regions, "world_regions_ggplot.txt")

Speaker.map=left_join(Speaker_Inst_Country,world_map, by="region")
Speaker.map$SpeakerCount=as.factor(Speaker.map$Speaker)

x=ggplot(Speaker.map, aes(long, lat, group = group))+
  geom_polygon(aes(fill = SpeakerCount ), color = "darkgrey", size=0.1)+
  #scale_fill_viridis_d(direction =1, option = "magma")
  scale_fill_manual(values=met.brewer("VanGogh3"))+
  theme_void()+
  theme(legend.position = "none")
`````

<img width="760" alt="image" src="https://user-images.githubusercontent.com/43852873/220715207-8d574ca0-0f4d-499e-bee0-b8f9983039eb.png">

