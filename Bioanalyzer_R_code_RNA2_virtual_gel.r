#set the working directory
setwd("H:\\virtual_gel")
myfiles<-list.files()
which(substr(myfiles,69,79)=="Results.csv")
results<-myfiles[which(substr(myfiles,69,79)=="Results.csv")]
n <- which(substr(myfiles,69,79)=="Results.csv")

targets <- read.delim(myfiles[n[1]], sep=",",header=F,skip=14,as.is=T)
ba.lane <- read.delim("sample_numbering.txt", sep="\t",skip=0,as.is=T,header=T)
targets[1:5,]
tbl<-which(targets[,1]=="Overall Results:")
which(targets[,1]=="Sample Name")
length(which(targets[,1]=="Sample Name"))
samples<-as.character(targets[which(targets[,1]=="Sample Name"),2])


#read in the targets file                               
colnms<-as.character(targets[(tbl[1]+1):(tbl[1]+4),1])
mat1 <- matrix(NA, ncol=length(colnms), nrow=length(samples)*length(results))
colnames(mat1)<-colnms
rownames(mat1)<-rep(samples,length(results))


read.delim("samples02.txt", sep="\t",header=T,as.is=T)
for(i in 1:length(results))

{

targets <- read.delim(results[i], sep=",",header=F,skip=14,as.is=T)
targets[1:5,]
which(targets[,1]=="Sample Name")
length(which(targets[,1]=="Sample Name"))
samples<-as.character(targets[which(targets[,1]=="Sample Name"),2])
tbl<-which(targets[,1]=="Overall Results:")
length(which(targets[,1]=="Overall Results:"))

tbl[1]
tbl[1]
colnms<-as.character(targets[(tbl[1]+1):(tbl[1]+4),1])

for (j in 1:length(samples))
{
mat1[((12*(i-1))+j),]<-as.matrix(targets[(tbl[j]+1):(tbl[j]+4),2])
}
}
mat1
cbind(mat1,ba.lane[,1:4])

#write.table(cbind(mat1,ba.lane[,1:4]),file="S:\\Genomics\\CWS\\Sanchez\\kim_tu\\ktu2\\BioAnalyzer\\conc.txt", sep="\t",col.names=NA)



###STOP###
names <- ba.lane[,2]

require(graphics)



#redo chip1
ramp <- colorRamp(c("darkmagenta","white"))
elec <- read.delim("samples02.txt", sep="\t",header=T,as.is=T)
#elec <- myfiles[which(substr(myfiles[1:100],69,74)=="Sample")]

#write.table(elec,file="S:\\Genomics\\CWS\\Sanchez\\kim_tu\\ktu2\\BioAnalyzer\\samples01.txt", sep="\t",col.names=NA)

plot(1:65,rep(80,65),xlab="BioAnalyzer Lane #",ylab="seconds",pch=NA,ylim=c(0,80))
for(i in 1:57)
{
for(k in 1:1060)
{
mini<-min(as.numeric(read.delim(elec[57,3], sep=",",header=F,skip=18,as.is=T)[,2]),na.rm=T)+1
palette(rgb(ramp(seq(.9,0,len=round(max(as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[,2]),na.rm=T)+1))),max=255))
y<-as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[k,1])
x<-round(as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[k,2])+mini ,0)
segments(i,y,i+1,y,col=x)
text(i+.5,75,format(elec[i,4],digits=4),cex=.7,srt = 90)
#text(i+.5,77,"ng/uL",cex=.7,srt = 60)
text(i+.5,10,elec[i,2],cex=.7,srt = 90)
}
}  



#######refined peakfinder####################



#read bioanalyzer data into a matrix called dta
#since the total RNA and mRNA assays run differently
#skip more lines for the total RNA assay


dta <- matrix(NA, nrow=1060, ncol=nrow(elec))
for(i in 1:nrow(elec))
{
x <- read.csv(as.character(elec[i,3]), header=F, skip=18, nrows=1060)
dta[,i] <- x[,2]
}
time<- x[,1]


#plot ladder and define peak threshold in HD
lad<-57
plot(time,dta[,lad],type="h",col="purple",xlab="seconds",ylab="Fluorescence Units",xlim=c(17,75),ylim=c(0,40))
abline(h=5)
which(dta[,lad] > 5)
#draw some segments and store max fluor in max.peaks matrix
t1<-17
t2<-65
p.row<-length(dta[which(time==t1):which(time==t1+1),lad])*(t2-t1+1)
peaks <- matrix(NA, nrow=p.row, ncol=1)
max.peaks <- matrix(NA, nrow=length(seq(t1,t2,.25)), ncol=1)

for (i in 1:length(seq(t1,t2,.25)))
{ 
ii<-seq(t1,t2,.25)[i]
x0<- ii
y0<- ii*.25
x<- ii+1
y<- ii*.25
segments(x0,y0,x,y)
max.peaks[i]<- max(dta[which(time==x0):which(time==x),lad])
}                
#define the middle of the segments
max.avg<-seq(t1+.5,t2+.5,.25)

#plot the maximum value of each segment
max.peaks
lines(max.avg,max.peaks,type="b",col="red")
#lines(max.avg,min.peaks,type="b",col="blue")

#find time associated with peaks > 5 fu and fill seconds matrix with them
sizes <- which(max.peaks>5)
points(max.avg[sizes],max.peaks[sizes],col="blue")
abline(h=5,col="red")
#define a slopes matrix to find the peaks
slopes <- matrix(NA, nrow=length(max.peaks[sizes])-1, ncol=1)
#fill slopes matrix with slopes
for (i in 1:length(max.peaks[sizes])-1)
{
slopes[i,1]<-max.peaks[sizes][i+1]-max.peaks[sizes][i]
}
#define where the slope is equal to zero
bp <- c("25 nt","200 nt","500 nt","1000 nt","2000 nt","4000 nt","6000 nt")
zilch<-which(slopes[,1]==0)
ladder <- matrix(NA, nrow=length(unique(max.peaks[sizes][zilch])), ncol=1)


for(j in 1:length(unique(max.peaks[sizes][zilch])))
{
abline(v=time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])])
ladder[i]<-time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])]
vert <- time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])]
horz <- unique(max.peaks[sizes][zilch])[j]+5
text(vert,horz,bp[j],cex=1.5,srt = 90)
}
###




 


#####FINAL CODE FOR VIRTUAL GEL###

bp <- c("- 25 nt","- 200 nt","- 500 nt","- 1000 nt","- 2000 nt","- 4000 nt","- 6000 nt")
ramp <- colorRamp(c("red","black"))
plot(1:65,rep(80,65),xlab="BioAnalyzer Lane #",ylab="seconds",pch=NA,ylim=c(0,80))
for(i in 1:57)
{
for(k in 1:1060)
{
mini<-min(as.numeric(read.delim(elec[57,3], sep=",",header=F,skip=18,as.is=T)[,2]),na.rm=T)+1
palette(rgb(ramp(seq(.9,0,len=round(max(as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[,2]),na.rm=T)+1))),max=255))
y<-as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[k,1])
x<-round(as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[k,2])+mini ,0)
segments(i,y,i+1,y,col=x)
text(i+.5,73,format(elec[i,4],digits=4),cex=.8,srt = 90)
text(i+.5,77,"ng/uL",cex=.7,srt = 90)
text(i+.5,12,elec[i,2],cex=.8,srt = 90)
for(j in 1:length(unique(max.peaks[sizes][zilch])))
{
#abline(v=time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])])
ladder[i]<-time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])]
vert <- time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])]
horz <- unique(max.peaks[sizes][zilch])[j]+5
text(60,vert,bp[j],cex=.8)
}
}  
}


##########################
#####FINAL CODE FOR RAW VIRTUAL GEL2##

bp <- c("- 25 nt","- 200 nt","- 500 nt","- 1000 nt","- 2000 nt","- 4000 nt","- 6000 nt")
ramp <- colorRamp(c("red","black"))
plot(1:65,rep(80,65),xlab="BioAnalyzer Lane #",ylab="seconds",pch=NA,ylim=c(0,80))
for(i in 1:57)
{
for(k in 1:1060)
{
mini<-min(as.numeric(read.delim(elec[57,3], sep=",",header=F,skip=18,as.is=T)[,2]),na.rm=T)+1
#read.delim(elec[3,3]...x=3, y=3, x is sample with most fluorescence
palette(rgb(ramp(seq(.9,0,len=round(max(as.numeric(read.delim(elec[3,3], sep=",",header=F,skip=18,as.is=T)[,2]),na.rm=T)+1))),max=255))
y<-as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[k,1])
x<-round(as.numeric(read.delim(elec[i,3], sep=",",header=F,skip=18,as.is=T)[k,2])+mini ,0)
segments(i,y,i+1,y,col=x)
text(i+.5,73,format(elec[i,4],digits=4),cex=.8,srt = 90)
text(i+.5,77,"ng/uL",cex=.7,srt = 90)
text(i+.5,12,elec[i,2],cex=.8,srt = 90)
for(j in 1:length(unique(max.peaks[sizes][zilch])))
{
#abline(v=time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])])
ladder[i]<-time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])]
vert <- time[which(unique(max.peaks[sizes][zilch])[j]==dta[,lad])]
horz <- unique(max.peaks[sizes][zilch])[j]+5
text(60,vert,bp[j],cex=.8)
}
}  
}

#########STOP##########







nano <- read.delim("ba_vs_nd.txt", sep="\t",header=T,,as.is=T)
barplot(nano[,4],col="lightseagreen",main=colnames(nano)[4],ylab="concentration (ng/uL)",names=nano[,2],horiz=FALSE,cex.names=.75,las=3,ylim=c(0,600))
barplot(nano[,3],col="lightgreen",main=colnames(nano)[3],ylab="concentration (ng/uL)",names=nano[,2],horiz=FALSE,cex.names=.75,las=3,ylim=c(0,600))
barplot(nano[,5],col="lightblue",main=colnames(nano)[5],ylab="concentration (ng/uL)",names=nano[,2],horiz=FALSE,cex.names=.75,las=3,ylim=c(0,3))
barplot(nano[,6],col="cadetblue2",main=colnames(nano)[6],ylab="concentration (ng/uL)",names=nano[,2],horiz=FALSE,cex.names=.75,las=3,ylim=c(0,3))


plot(nano[,3],nano[,4], main="BA vs ND", pch=19, col="blue",xlab="BA conc",ylab="ND conc")
