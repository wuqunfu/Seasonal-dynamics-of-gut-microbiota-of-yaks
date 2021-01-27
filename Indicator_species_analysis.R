#Indicator species

install.packages("indicspecies")

library(indicspecies)

#In this case, I have one column called “season”, that contains the information for grouping my samples into “spring”,“summer”,"autumn","winter" categories. Next I want to create a data frame with all my abundance data, and a vector that contains the information from the “season” column:

pc = read.csv("C:/Users/Fuyu Shi/Desktop/bar/Indicator species/grass.csv", header= TRUE)
#Because my abundance data starts in the 3rd column of my original “pc” data frame. I put a “3” in the code above. Change this number to match the column where your abundance data starts. Also change “season” to whatever you have called your grouping variable.
abund = pc[,3:ncol(pc)]
#Next we can run the indicator species command:

time = pc$Season
#multipatt is the name of the command from the indicspecies package. The mulitpatt command results in lists of species that are associated with a particular group of samples. If your group has more than 2 categories, multipatt will also identify species that are statistically more abundant in combinations of categories.

#The parameters for multipatt are as follows:
  
#  the community abundance matrix:"abund"
#the vector that contains your sample grouping information: "season"
#the function that multipatt is using to identify indicator species: func = "r.g"
#the number of permutations used in the statistical test: control = how(nperm=9999)


inv = multipatt(abund, time, func = "r.g", control = how(nperm=9999))
summary(inv)




install.packages("tidyverse") 

if(!require(devtools))	install.packages("devtools") devtools::install_github("kassambara/ggpubr") 
install.packages("dplyr")
install.packages("ggpubr") 

library("ggplot2")
library("ggpubr")
library("tidyverse")

library("readr") 
library("dplyr")


data=read.table("C:/Users/Fuyu Shi/Desktop/bar/indicator species/grass.txt",header=TRUE)  
theme_set(theme_pubr())
#df	<-	data	%>%	group_by(plant)	%>%	summarise(counts	=	n())
#rownames(df)=df[1,]
#df=c(df,stat)
df	<-	data	%>%	arrange(desc(plant))	%>%	mutate(plant=plant, Indicator_value=Indicator_value) 

df$plant <- factor(df$plant,levels=c("Rosaceae","Salicaceae","Elaeagnaceae","Scrophulariaceae","Compositae","Gramineae","Leguminosae","Gentianaceae","Polygonaceae","Umbelliferae","Cyperaceae","Ranunculaceae"))
head(df,	4) 

par(mar=c(4,6,2,20))
col1 = rgb(250/255,60/255,60/255)
col2 = rgb(0/255,220/255,0/255)
col3 = rgb(30/255,60/255,255/255)
col4 = rgb(0/255,200/255,200/255)
col5 = rgb(240/255,0/255,130/255)
col6 = rgb(230/255,220/255,50/255)
col7 = rgb(240/255,130/255,40/255)
col8 = rgb(160/255,0/255,200/255)
col9 = rgb(160/255,230/255,50/255)
col10 = rgb(0/255,160/255,255/255)
col11 = rgb(230/255,175/255,45/255)
col12 = rgb(0/255,210/255,140/255)
col13 = rgb(247/255,153/255,209/255)

color <- c(col1,col2,col4,col6,col5,col3,col9,col7,col8,col10,col12,col11,col13)
color <- c(col1,col2,col3,col4,col5,col6,col7,col8,col15,col16,col10)

mydata<-read.csv("C:/Users/Fuyu Shi/Desktop/bar/indicator species/grass_family_vertical.csv",sep=",",na.strings="NA",stringsAsFactors=FALSE)

mydata$family <- factor(mydata$family,levels=c("Gentianaceae","Leguminosae","Gramineae","Compositae","Scrophulariaceae","Elaeagnaceae","Salicaceae","Ranunculaceae","Cyperaceae","Umbelliferae","Rosaceae","Polygonaceae"))
head(mydata,	4) 
color <- c(col1,col2,col4,col6,col5,col3,col9,col7,col8,col10,col12,col11,col13)

ggplot(mydata, aes(Indicator_value, family)) +
  geom_segment(aes(x=0, 
                   xend=Indicator_value,
                   y=family, 
                   yend=family),color	=	"lightgray",	size	=	1	)+
  geom_point(aes(color	=	family),	size	=	19, shape = 19)+
  scale_x_continuous(limits = c(0,1))+theme_classic()+
  ggpubr::color_palette(c(col12,col8,col4,col6,col10,col13,col11,col7,col3,col9,col2,col1))+
  theme(legend.position="none",
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=10,face="plain",color="black"),
        legend.title=element_text(size=14,face="plain",color="black")
  )