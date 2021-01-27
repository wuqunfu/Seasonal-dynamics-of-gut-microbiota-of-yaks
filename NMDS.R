#NMDS
library(ggplot2)
library(vegan)
otus = read.csv(file = "C:/Users/Fuyu Shi/Desktop/bar/NMDS/grass.csv", header = TRUE, check.names = FALSE, row.names = 1)
metadata = read.csv(file = "C:/Users/Fuyu Shi/Desktop/bar/NMDS/EPIPHYTE_METADATA1.csv", header = TRUE, check.names = FALSE, row.names = 1)
good_samples <- colnames(otus[(colSums(decostand(otus,"pa")) >= 1)])     # decostand(x,"pa") counts presence/absence
otus = otus[,good_samples]
metadata = metadata[good_samples,]
t_otus <- as.data.frame(t(otus))
min_depth = min(colSums(otus))
t_otus_rarefied <- as.data.frame(round(rrarefy(t_otus, min_depth)))
sqrt_t_otus_rarefied = sqrt(t_otus_rarefied)
rank.totus <- rankindex(as.matrix(sqrt_t_otus_rarefied), t_otus_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rank was given by the", names(sort(rank.totus, decreasing = TRUE)[1]), "method."))
otus_dist = as.matrix((vegdist(t_otus_rarefied, "bray")))
write.table(otus_dist, file="bray dissimilarity.biom", append = F, sep="\t", quote=F, row.names=T, col.names=T)
#perform NMDS
NMDS = metaMDS(otus_dist)
#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Season = metadata$Season)
NMDS$Season <- factor(NMDS$Season,levels=c("Spring","Summer","Autumn","Winter"))
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

col14 = rgb(109/255,227/255,255/255)
col16 = rgb(175/255,240/255,255/255)
col15 = rgb(153/255,153/255,153/255)



cols <- c(col7, col1,col8,col13)
head(NMDS)
 p1=ggplot(NMDS, aes(x=MDS1, y=MDS2, col=NMDS$Season)) +
  geom_point(data=metadata,size=6,shape=17,alpha=0.8) +
  stat_ellipse(type = "norm", linetype = 2)+
  theme_bw() +
  scale_color_manual(values=cols)+
  guides(colour=guide_legend(override.aes = list(shape = 17)))+
  labs(title = "Trans-humance management",x="NMDS1",y="NMDS2",size=14)+
  theme(axis.text.x = element_text(color="black", size=12), axis.text.y = element_text(color="black", size=12),
   axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),legend.text=element_text(lineheight=.8,size=13), legend.title = element_text(size = 13))+
   theme(plot.title = element_text(size = 12, face = "bold") , legend.title=element_text(size=12) ,legend.text=element_text(size=10))
 p1

#Anosim analysis
anosim_location = anosim(otus_dist, metadata$Season,permutations = 9999)
anosim_location # take a look at results
summary(anosim_location)
plot(anosim_location)

adonis_location = adonis(otus_dist ~ Season, metadata,permutations = 9999)
adonis_location # take a look at results; there are no summary() or plot() methods included
distance_data<-vegdist(otus_dist)
anova(betadisper(distance_data,metadata$Season))