setwd("C:/Users/Kexin/Dropbox/CogMaster/S2/SuperSem")
dat = read.csv("call.csv", header = TRUE, stringsAsFactors=FALSE)

# Data preprocessing

## convert NA and "" as 0; percentage as 1
dat[is.na(dat)] = 0
dat[,7:39][dat[,7:39] == ""] = 0
dat[,7:39] = sapply(dat[,7:39], as.character)

## we don't consider % first, convert them to 1
dat[,7:39][apply(dat[,7:39],2,function(x){grepl("%",x)})] = 1

## Method 1: generate a code for each row based on 0/1 values in context
dat$code = ""
dat$code = apply(dat[,7:39],1, paste,collapse = "")


## Method 2: 2 dimensions = [call, context], value/color = species
### create specie index: specielist
nspecie = length(unique(dat$species))
specielist = data.frame(matrix(ncol = 2, nrow = nspecie))
specielist[,1:2] = cbind(unique(dat$species),(1:nspecie) )

### add specie index to original data
dat$specieInd = sapply(dat$species, function(x){specielist[match(x,specielist[,1]),2]})

### create the final matrix
ncall = length(unique(dat$call))
ncontext = 39 - 7 + 1

dat2 = data.frame(matrix(ncol = ncall, nrow = ncontext))
colnames(dat2) <-  unique(dat$call)
rownames(dat2) <-  colnames(dat[,7:39])

### give values to matrix
for (i in 7:39){
row_t = which(dat[,i] == 1)
if (length(row_t) == 0) next
call_t = dat$call[row_t]
specie_t = as.matrix(dat$specieInd[row_t])
for (j in 1:nrow(specie_t )){
dat2[i-6,call_t[j]] = paste(specie_t[j,1], collapse = ",")
}}

### viz
#dat2[is.na(dat2)] = NA
x = (1:ncol(dat2))
y = (1:nrow(dat2))
dat2 = as.matrix(dat2)
dat2 = apply(dat2,2,as.numeric)

par(mar=c(5.1, 5.1, 4.1, 13.1), xpd=TRUE)
image(y = (1:ncol(dat2)), x=(1:nrow(dat2)), z = dat2, col=specielist[,2],xlab = "context",ylab = "calls")

abline(v=c(12.5,19.5,24.5,29.5),lty=1, col = "grey")
grid( nrow(dat2),ncol(dat2), lty = 6, col = "grey")
op <- par(cex = 0.7)
legend(35,150, legend = specielist[,1],pch = 15, col = specielist[,2])


#dat3: in dat 2, convert "species label" to 1
dat3 = dat2
dat3[!is.na(dat3)] = 1
dat3[is.na(dat3)] = 0
v = (1:33)
rownames(dat3) = c(as.character(v))


#------------------------------------------MDS-------------------------------------------------------------
## create distance metric
dis = matrix(nrow = 33, ncol =33)
for (i in 1:nrow(dat3)){
  for (j in 1:nrow(dat3)){
  dis[i,j] = sum(dat3[j,which(dat3[i,]==1)])
  }}

### solution 1: auto dist matrix --> not work
#dis = dist(dat3)
#dis = 1/dis

### solution 2: calc dis matrix by 1/#points
dis = 1/dis
dis[which(dis == Inf)] = 1000
dis=as.dist(dis)

### mds + plot
####2d
fit = cmdscale(dis, k =2)
x = fit[,1]
y = fit[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",  main="Context Metric MDS",	type="n")
text(x, y, labels = row.names(dat3), cex=1.2, col = c(1,2,3))

####3d
fit = cmdscale(dis, k =3)
x = fit[,1]
y = fit[,2]
z = fit[,3]
plot3d(x,y,z)
text3d(x=x,y=y,z=z,texts=c(as.character(v)),col='red') 

####1d

fit = cmdscale(dis, k =2)
x = as.vector(fit[,1])
y = as.numeric(rep(1, 33))
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",  main="Context Metric	MDS",	type="n")
text(x, y, labels = row.names(dat3), cex=1,col = c(1:15))


#--------------------------------------------------------clustering------------------------------------------
## k-means

###Elbow Method for finding the optimal number of clusters
set.seed(123)
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(fit, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

### clustering plot
kc <- kmeans(fit[,1:2],4)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",  main="K-means clustering for Context MDS Metric",	type="n")
text(x, y, labels = row.names(dat3), cex=1,col = kc$cluster)

kc$cluster

#-------------------------------------------------------------------PCA--------------------------------------------------------
## Method 1: dim = calls
#call.pca <- prcomp(dat3,center = TRUE) 
#print(call.pca)
#plot(call.pca, type = "l")
#plot(call.pca$x[,1],call.pca$x[,2])
#text(call.pca$x[,1],call.pca$x[,2], labels = row.names(dat3), cex=1,col = c(1:4))


##Method 2: dim = dist in MDS
####3d MDS
fit = cmdscale(dis, k =3)
x = fit[,1]
y = fit[,2]
z = fit[,3]
#plot3d(x,y,z)
#text3d(x=x,y=y,z=z,texts=c(as.character(v)),col='red') 
mds.pca <- prcomp(fit,
                 center = TRUE,
                 scale. = TRUE) 
print(mds.pca)
plot(mds.pca, type = "l")
biplot(mds.pca,col=c(1:4)) 
summary(mds.pca)
##Method 3: dim = context, data = calls
cont.pca <- prcomp(t(dat3),
                  center = TRUE) 
print(cont.pca)
plot(cont.pca, type = "l")

plot(cont.pca$x[,1],cont.pca$x[,2],type = "p",cex = 1)
text(cont.pca$x[,1],cont.pca$x[,2], labels = row.names(dat3), cex=0.8,col = c(1:4))

#plot3d(cont.pca$rotation[,1],cont.pca$rotation[,2],cont.pca$rotation[,3],type = "p",cex = 0.2)
#text3d(cont.pca$rotation[,1],cont.pca$rotation[,2], cont.pca$rotation[,3],texts = row.names(dat3), cex=0.8,col = c(1:4))

#----------------------------------------------------------------HC------------------------------------------------------
cont.class = c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5)
hc <- hclust(dis)
plot(hc, main = "Hierarchical Clustering for Context")

# install the package:
if (!require('dendextend')) install.packages('dendextend'); library('dendextend')

# Create the dend:
dend <- as.dendrogram(hc)
type = kc$cluster
type <- factor(type)
n_types <- length(unique(type))
cols_5 <- colorspace::rainbow_hcl(n_types, c = 200, l  = 50)
col_type <- cols_5[type]
labels_colors(dend) <- col_type[order.dendrogram(dend)]
dend <- color_branches(dend, k=1)
plot(dend)

summary(cont.pca)
write.csv(x = dat, file = "call-0.csv")

