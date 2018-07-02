#!/usr/bin/sh/
a<-read.table("Raw_data.txt",check.names=F)
b<-read.table("clean_data.txt",check.names=F)
table<-merge(a,b,by="V2")
colnames(table)<-c("sample","Raw","clean")
chim<-read.table("otus_table1.txt",header=T,row.names=1,sep="\t",check.names=F)
data1<-data.frame()
k=1
for (i in (1:length(colnames(chim))))
{d<-(sum(chim[,i]))
e<-(sum(chim[,i]>0))
data1[k,1]<-colnames(chim[i])
data1[k,2]<-d
k<-k+1
}
colnames(data1)<-c("sample","Chimerta_check")
table1_3<-merge(table,data1,by="sample")


c<-read.table("otus_table.txt",header=T,row.names=1,sep="\t",check.names=F)
data<-data.frame()
j=1
for (i in (1:length(colnames(c))))
{d<-(sum(c[,i]))
e<-(sum(c[,i]>0))
data[j,1]<-colnames(c[i])
data[j,2]<-d
data[j,3]<-e
j<-j+1
}
colnames(data)<-c("sample","OTU_seq","OTU")
all_table<-merge(table1_3,data,by="sample")
write.table(all_table,"sample2sequence_summary.xls",col.names=T,sep="\t",quote=F,row.names=F)
