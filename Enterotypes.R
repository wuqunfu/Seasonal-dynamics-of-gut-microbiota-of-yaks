###########Get enterotypes######################
mothur > get.communitytype(shared=genus_abund.shared, method=pam, label=0.03)


###########Enterotypes#####################################################################################################################3

install.packages("clusterSim")
install.packages("flexclust")
install.packages("fpc")
library(cluster)
library(clusterSim)
library(ade4)
library(vegan)
library(flexclust)
library(fpc)

###read data
setwd("C:/Users/Fuyu Shi/Desktop/bar/enterope")
data=read.table("all.txt", header=T, row.names=1, dec=".", sep="\t")
data
#a probability distribution distance metric related to Jensen-Shannon divergence (JSD) to cluster the samples.
#Here is an example complete function to create JSD distance matrix in R:
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
#we can apply this function on our data
data.dist=dist.JSD(data)

###calculate the PS based on PAM cluster using code from Dr.Knights.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Import the code, relicate five runs:
PS1=prediction.strength(data.dist)
PS1
PS2=prediction.strength(data.dist)
PS2
PS3=prediction.strength(data.dist)
PS3
PS4=prediction.strength(data.dist)
PS4
PS5=prediction.strength(data.dist)
PS5

PS <- (PS1+PS2+PS3+PS4+PS5)/5  ###Artifical
PS
#We can get average PS index for supporting cluster division.

1.2. ##calculate the rJSD distance.
data.dist=sqrt(data.dist)####data.dist was generated from the above code.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Import the code, relicate five runs:
PS1=prediction.strength(data.dist)
PS1
PS2=prediction.strength(data.dist)
PS2
PS3=prediction.strength(data.dist)
PS3
PS4=prediction.strength(data.dist)
PS4
PS5=prediction.strength(data.dist)
PS5

PS <- (PS1+PS2+PS3+PS4+PS5)/5###Artifical
PS

#We can get average PS index for supporting cluster division.
1.3. ###calculate the Bray-Curtis distance(BC)Using the "Vegan" package and then calculate the PS using the code of Dr.Knights.
data.transfer<-t(data) ###transform the rows and column
data.dist<-vegdist(data.transfer, method="bray", binary=FALSE)
data.dist
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Import the code, relicate five runs:
PS1=prediction.strength(data.dist)
PS1
PS2=prediction.strength(data.dist)
PS2
PS3=prediction.strength(data.dist)
PS3
PS4=prediction.strength(data.dist)
PS4
PS5=prediction.strength(data.dist)
PS5

PS <- (PS1+PS2+PS3+PS4+PS5)/5###Artifical
PS



2. ###calculating the Calinski-Harabasz (CH) Index based on the "clusterSim" pakage.  Large sample size!!!!!
pam(as.dist(x), k, diss=TRUE) # x is a distance matrix and k the number of clusters
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
nclusters=NULL  #### if we use the BC,rJSD and Euclidean distance, do not use t(data), just use data.transfer.

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

### evaluating the cluster number using CH index.
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")### evaluating the cluster number using CH index.
###determining the cluster number
data.cluster=pam.clustering(data.dist, k=3)
data.cluster
3. ##calculating the aveage silhouette index (SI)based on the "cluster" package
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
obs.silhouette


###Based on the PS result, we determined the best nurmber of cluster is k=3 using JSD.
data.dist=dist.JSD(data)

data.cluster=pam.clustering(data.dist, k=3)
data.cluster
4. ###plot the cluster result based on the JDS(JDS的PS值最高).
data.dist=dist.JSD(data)

data.cluster=pam.clustering(data.dist, k=3)
data.cluster
#Principal coordinates analysis (PCoA)
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
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

color <- c(col1,col3,col7)

s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=c(2,4))  ###obs.pcoa$li represent the dataframe dfxy.
s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=c(1,8))
s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=palette(rainbow(10)))
s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=color)

plot(s.class, fac=as.factor(data.cluster))
labs <- rownames(obs.pcoa$li)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=3)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col="green")
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col="#8B6914")
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=1:8)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=1:16)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=palette(rainbow(10)))
#Anosim analysis
data.group=read.table("simper_JSD_group.txt", header=T, row.names=1, dec=".", sep="\t")

anosim_location = anosim(data.dist, data.group$JSD)
anosim_location # take a look at results
anosim_location$statistic
anosim_location$signif
cdis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="1"]
fdis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="2"]
adis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="3"]

bedis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="Between"]

par(bty='l')
boxplot(bedis,cdis,fdis,adis,
        xaxt='n',ylab='Distance rank',col=c(col4,col1,col3,col7),ylim=c(0,50000),
       width=c(0.2,0.2,0.2,0.2),border=c('turquoise4','red4','navyblue','orange4'))
axis(1,at=c(1,2,3,4),label=c('Between','P1','P2','P3'))
text(3,48000,paste('p value', anosim_location$signif,sep=' = ') ,cex=1.2)
text(2,48000,paste('R ',anosim_location$statistic ,sep=' = ') ,cex=1.2)
data=t(data)  ###### row with the variables (aniaml name), and the column with the taxa name.
data.group=read.table("simper_JSD_group.txt", header=T, row.names=1, dec=".", sep="\t")
data.group
##### this file define the group for the samples, which could be obtained from the "data.cluster" using "as.factor" fuction, like "data.factor=as.factor(data.cluster)". 
sim <- with(data.group, simper(data,permutations = 0, group=c("1","2","3")))
summary(sim,ordered = TRUE) ###ordered "the species be ordered by their average contribution".



#######################################BC DISSIMILARITY#############################################################33

###Based on the PS result, we determined the best nurmber of cluster is k=3 using JSD.
data.transfer<-t(data) ###transform the rows and column
data.dist<-vegdist(data.transfer, method="bray", binary=FALSE)
data.dist

data.cluster=pam.clustering(data.dist, k=3)
data.cluster
4. ###plot the cluster result based on the BC
#Principal coordinates analysis (PCoA)
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
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

color <- c(col1,col3,col7)

s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=c(2,4))  ###obs.pcoa$li represent the dataframe dfxy.
s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=c(1,8))
s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=palette(rainbow(10)))
s.class(obs.pcoa$li, fac=as.factor(data.cluster), cellipse = 1.5,grid=F,col=color)

plot(s.class, fac=as.factor(data.cluster))
labs <- rownames(obs.pcoa$li)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=3)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col="green")
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col="#8B6914")
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=1:8)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=1:16)
text(obs.pcoa$li,labels=labs,adj=c(-.1,-.8),cex=0.6,col=palette(rainbow(10)))
###identifying taxa contributing to groups based on BC distance, the SIMPER method in "vegan" package was used.

data=t(data)  ###### row with the variables (aniaml name), and the column with the taxa name.
data.group=read.table("simper_BC_group.txt", header=T, row.names=1, dec=".", sep="\t")
data.group
##### this file define the group for the samples, which could be obtained from the "data.cluster" using "as.factor" fuction, like "data.factor=as.factor(data.cluster)". 
sim <- with(data.group, simper(data,permutations = 0, group=c("1","2","3")))
summary(sim,ordered = TRUE) ###ordered "the species be ordered by their average contribution".







#Anosim analysis
data.group=read.table("simper_BC_group.txt", header=T, row.names=1, dec=".", sep="\t")

anosim_location = anosim(data.dist, data.group$entrotype)
anosim_location # take a look at results
anosim_location$statistic
anosim_location$signif
cdis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="1"]
fdis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="2"]
adis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="3"]

bedis=anosim_location$dis.rank[as.vector(anosim_location$class.vec)=="Between"]

par(bty='l')
boxplot(bedis,cdis,fdis,adis,
        xaxt='n',ylab='Distance rank',col=c(col4,col1,col3,col7),ylim=c(0,50000),
        width=c(0.2,0.2,0.2,0.2),border=c('turquoise4','red4','navyblue','orange4'))
axis(1,at=c(1,2,3,4),label=c('Between','P1','P2','P3'))
text(3,48000,paste('p value', anosim_location$signif,sep=' = ') ,cex=1.2)
text(2,48000,paste('R ',anosim_location$statistic ,sep=' = ') ,cex=1.2)


#compared with BC and JSD, we concluded that BC is better,so the BC will be used in the next analysis.


###PCoA analysis for the raw data based on BC distance.

data=t(data)###transform the data format.
PCoA.res<-capscale(data~1,distance="bray")
PCoA.res

######VEGDIST: Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao" or "mahalanobis".
######DIST: the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.

summary(PCoA.res) 


###kernel density estimation(kde) using the "ks" package, kde is a popular tool for visualising the distribution of data.
???
  
  # Kernel Density Plot
  data.kes<-kde(x=data,positive=TRUE)

###identifying taxa contributing to groups based on BC distance, the SIMPER method in "vegan" package was used.

data=t(data)  ###### row with the variables (aniaml name), and the column with the taxa name.
data.group=read.table("simper1.group.txt", header=T, row.names=1, dec=".", sep="\t")
data.group
##### this file define the group for the samples, which could be obtained from the "data.cluster" using "as.factor" fuction, like "data.factor=as.factor(data.cluster)". 
sim <- with(data.group, simper(data,permutations = 0, group=c("P1","P2","P3")))
summary(sim,ordered = TRUE) ###ordered "the species be ordered by their average contribution".










#######################################Fisher-test

Genus <-
  matrix(c(15,3,3,15),
         nrow = 2,
         dimnames =
           list(c("P1", "P2"),
                c("Overrepresented", "Not Overrepresented")))
fisher.test(Genus, alternative = "two.sided", conf.level = 0.99, simulate.p.value = FALSE, B = 2000)
p=fisher.test(Genus, alternative = "two.sided", conf.level = 0.99, simulate.p.value = FALSE, B = 10000)[1]
p.adjust(p, method = "fdr",n=length(p))




################################################
###read data
setwd("C:/Users/Fuyu Shi/Desktop/bar/enterope/BC")
library(ggplot2)
data=read.table("BC_grass_frequency_season.txt", header=T, row.names=1, dec=".", sep="\t")
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

color <- c(col1,col7,col3)

data$season <- factor(data$season,levels=c("spring","summer","autumn","winter"))



#	Use	position	=	position_dodge()
p	<-	ggplot(data,	aes(x	=	season,	y	=	value))	+	geom_bar(			
  aes(color	=	entrotype,	fill	=	entrotype),	
  stat	=	"identity",	position	=	position_dodge(0.8),				
  width	=	0.7	)	+	
  coord_cartesian(ylim=c(0,100)) + 
  scale_color_manual(values	=	c(col1,col3,col7))+
  scale_fill_manual(values	=	c(col1,col3,col7)) +
  theme_classic()+
  theme(axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=2), 
        axis.title.y = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) 
p




