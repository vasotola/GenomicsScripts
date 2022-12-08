## Tally across SNP Windows----
setwd("~/Dropbox/Research/Georgia/UGA_Research/Mimulus/RIL_Computing/LinkageMapping/LepMap_pubRun/Windowed_run")

####read in and combine outfiles from vcftools -012----
gtype<-read.table("out.012",header=F) #47851 markers
gtype[1:10,1:5]
dim(gtype)
gtype<-gtype[,-1]
gtype[1:10,1:5]
dim(gtype)

# replace missing (-1) with NA
gtype <- replace(gtype, gtype == -1, NA) # NA
gtype[1:10,1:5]

# combine with inds and loci for full matrix
inds<-read.table("out.012.indv",header=F)
gtype.inds<-cbind(inds,gtype)
gtype.inds[1:10,1:5]
loci<-read.table("out.012.pos",header=F)
loci[1:5,]
loci$V3<-paste(loci$V1, loci$V2, sep="_")
loci[1:5,]
loci_comb<-loci$V3
length(loci_comb)
dim(gtype.inds)
colnames(gtype.inds)[2:47852]<-loci_comb
gtype.inds[1:5,1:5]

# transpose
t.gtype.inds<-t(gtype.inds)

sum(is.na(t.gtype.inds))
# write out if needed
#write.csv(t.gtype.inds,"combined_VCFgtypeMatrix.csv")

gtype[1:10,1:5]
dim(gtype)
t.gtype<-t(gtype)
t.gtype<-as.data.frame(t.gtype)
t.gtype[1:5,1:5]

dim(t.gtype)

# add in chromosome N to separate out
gtypematrix<-cbind(loci$V1,t.gtype)

##Start here if already formatted data above----
# write out
#write.csv(gtypematrix,"gtypematrix.csv")
#str(gtypematrix)
gtypematrix<-read.csv("gtypematrix.csv")
gtypematrix<-gtypematrix[,-1]

# split to window call over each chromosome separately
chromes <- split(gtypematrix,gtypematrix$loci.V1)


####Get started with tallying----
library(tidyverse)
inds<-read.table("out.012.indv",header=F)
str(inds)
inds<-inds$V1

loci<-read.table("out.012.pos",header=F)
loci[1:5,]
loci$V3<-paste(loci$V1, loci$V2, sep="_")
loci[1:5,]
loci_comb<-loci$V3

loci.bp<-split(loci,loci$V1)

######chr1----
chromes$CE10_chr1[1:5,1:5]
chr1<-chromes$CE10_chr1[,-1]
dim(chr1)
colnames(chr1)[1:147]<-as.character(inds)
chr1[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr1[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr1[,-1])))
chr1<-cbind(list.rep,chr1)
names(chr1)[1]<-"list.rep"

str(chr1)
chr1[1:5,1:5]
dim(chr1)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr1 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr1<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr1[,1]
d.sum.chr1<-d.sum.chr1[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr1<-matrix(NA,nrow=nrow(d.sum.chr1),ncol=ncol(d.sum.chr1))

for(i in 1:nrow(call.chr1)){
  for(j in 1:ncol(call.chr1)){
    call.chr1[i,j]<-ifelse(d.sum.chr1[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr1,file="d.sum.chr1.csv")
#write.csv(call.chr1,file="call.chr1.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr1
tally.d_ab_freq<-tally.d_ab/d.sum.chr1
tally.d_bb_freq<-tally.d_bb/d.sum.chr1

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr1 <- matrix(NA,nrow=nrow(call.chr1),ncol=ncol(call.chr1))

for(i in 1:nrow(ave.chr.called.chr1)){
  for(j in 1:ncol(ave.chr.called.chr1)){
    ave.chr.called.chr1[i,j] <- ifelse(call.chr1[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}

write.csv(ave.chr.called.chr1,file="ave.chr.called88.chr1.csv")

######chr2----
chromes$CE10_chr2[1:5,1:5]
chr2<-chromes$CE10_chr2[,-1]
dim(chr2)
colnames(chr2)[1:147]<-as.character(inds)
chr2[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr2[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr2[,-1])))
chr2<-cbind(list.rep,chr2)
names(chr2)[1]<-"list.rep"

str(chr2)
chr2[1:5,1:5]
dim(chr2)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr2 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr2<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr2[,1]
d.sum.chr2<-d.sum.chr2[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr2<-matrix(NA,nrow=nrow(d.sum.chr2),ncol=ncol(d.sum.chr2))

for(i in 1:nrow(call.chr2)){
  for(j in 1:ncol(call.chr2)){
    call.chr2[i,j]<-ifelse(d.sum.chr2[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr2,file="d.sum.chr2.csv")
#write.csv(call.chr2,file="call.chr2.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr2
tally.d_ab_freq<-tally.d_ab/d.sum.chr2
tally.d_bb_freq<-tally.d_bb/d.sum.chr2

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr2 <- matrix(NA,nrow=nrow(call.chr2),ncol=ncol(call.chr2))

for(i in 1:nrow(ave.chr.called.chr2)){
  for(j in 1:ncol(ave.chr.called.chr2)){
    ave.chr.called.chr2[i,j] <- ifelse(call.chr2[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



write.csv(ave.chr.called.chr2,file="ave.chr.called88.chr2.csv")


######chr3----
chromes$CE10_chr3[1:5,1:5]
chr3<-chromes$CE10_chr3[,-1]
dim(chr3)
colnames(chr3)[1:147]<-as.character(inds)
chr3[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr3[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr3[,-1])))
chr3<-cbind(list.rep,chr3)
names(chr3)[1]<-"list.rep"

str(chr3)
chr3[1:5,1:5]
dim(chr3)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr3 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr3<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr3[,1]
d.sum.chr3<-d.sum.chr3[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr3<-matrix(NA,nrow=nrow(d.sum.chr3),ncol=ncol(d.sum.chr3))

for(i in 1:nrow(call.chr3)){
  for(j in 1:ncol(call.chr3)){
    call.chr3[i,j]<-ifelse(d.sum.chr3[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr3,file="d.sum.chr3.csv")
#write.csv(call.chr3,file="call.chr3.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr3
tally.d_ab_freq<-tally.d_ab/d.sum.chr3
tally.d_bb_freq<-tally.d_bb/d.sum.chr3

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr3 <- matrix(NA,nrow=nrow(call.chr3),ncol=ncol(call.chr3))

for(i in 1:nrow(ave.chr.called.chr3)){
  for(j in 1:ncol(ave.chr.called.chr3)){
    ave.chr.called.chr3[i,j] <- ifelse(call.chr3[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



write.csv(ave.chr.called.chr3,file="ave.chr.called88.chr3.csv")


######chr4----
chromes$CE10_chr4[1:5,1:5]
chr4<-chromes$CE10_chr4[,-1]
dim(chr4)
colnames(chr4)[1:147]<-as.character(inds)
chr4[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr4[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr4[,-1])))
chr4<-cbind(list.rep,chr4)
names(chr4)[1]<-"list.rep"

str(chr4)
chr4[1:5,1:5]
dim(chr4)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr4 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr4<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr4[,1]
d.sum.chr4<-d.sum.chr4[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr4<-matrix(NA,nrow=nrow(d.sum.chr4),ncol=ncol(d.sum.chr4))

for(i in 1:nrow(call.chr4)){
  for(j in 1:ncol(call.chr4)){
    call.chr4[i,j]<-ifelse(d.sum.chr4[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr4,file="d.sum.chr4.csv")
#write.csv(call.chr4,file="call.chr4.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr4
tally.d_ab_freq<-tally.d_ab/d.sum.chr4
tally.d_bb_freq<-tally.d_bb/d.sum.chr4

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr4 <- matrix(NA,nrow=nrow(call.chr4),ncol=ncol(call.chr4))

for(i in 1:nrow(ave.chr.called.chr4)){
  for(j in 1:ncol(ave.chr.called.chr4)){
    ave.chr.called.chr4[i,j] <- ifelse(call.chr4[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



write.csv(ave.chr.called.chr4,file="ave.chr.called88.chr4.csv")


######chr5----
chromes$CE10_chr5[1:5,1:5]
chr5<-chromes$CE10_chr5[,-1]
dim(chr5)
colnames(chr5)[1:147]<-as.character(inds)
chr5[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr5[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr5[,-1])))
chr5<-cbind(list.rep,chr5)
names(chr5)[1]<-"list.rep"

str(chr5)
chr5[1:5,1:5]
dim(chr5)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr5 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr5<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr5[,1]
d.sum.chr5<-d.sum.chr5[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr5<-matrix(NA,nrow=nrow(d.sum.chr5),ncol=ncol(d.sum.chr5))

for(i in 1:nrow(call.chr5)){
  for(j in 1:ncol(call.chr5)){
    call.chr5[i,j]<-ifelse(d.sum.chr5[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr5,file="d.sum.chr5.csv")
#write.csv(call.chr5,file="call.chr5.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr5
tally.d_ab_freq<-tally.d_ab/d.sum.chr5
tally.d_bb_freq<-tally.d_bb/d.sum.chr5

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr5 <- matrix(NA,nrow=nrow(call.chr5),ncol=ncol(call.chr5))

for(i in 1:nrow(ave.chr.called.chr5)){
  for(j in 1:ncol(ave.chr.called.chr5)){
    ave.chr.called.chr5[i,j] <- ifelse(call.chr5[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



write.csv(ave.chr.called.chr5,file="ave.chr.called88.chr5.csv")


######chr6----
chromes$CE10_chr6[1:5,1:5]
chr6<-chromes$CE10_chr6[,-1]
dim(chr6)
colnames(chr6)[1:147]<-as.character(inds)
chr6[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr6[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr6[,-1])))
chr6<-cbind(list.rep,chr6)
names(chr6)[1]<-"list.rep"

str(chr6)
chr6[1:5,1:5]
dim(chr6)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr6 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr6<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr6[,1]
d.sum.chr6<-d.sum.chr6[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr6<-matrix(NA,nrow=nrow(d.sum.chr6),ncol=ncol(d.sum.chr6))

for(i in 1:nrow(call.chr6)){
  for(j in 1:ncol(call.chr6)){
    call.chr6[i,j]<-ifelse(d.sum.chr6[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr6,file="d.sum.chr6.csv")
#write.csv(call.chr6,file="call.chr6.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr6
tally.d_ab_freq<-tally.d_ab/d.sum.chr6
tally.d_bb_freq<-tally.d_bb/d.sum.chr6

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr6 <- matrix(NA,nrow=nrow(call.chr6),ncol=ncol(call.chr6))

for(i in 1:nrow(ave.chr.called.chr6)){
  for(j in 1:ncol(ave.chr.called.chr6)){
    ave.chr.called.chr6[i,j] <- ifelse(call.chr6[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



write.csv(ave.chr.called.chr6,file="ave.chr.called88.chr6.csv")

######chr7----
chromes$CE10_chr7[1:5,1:5]
chr7<-chromes$CE10_chr7[,-1]
dim(chr7)
colnames(chr7)[1:147]<-as.character(inds)
chr7[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr7[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr7[,-1])))
chr7<-cbind(list.rep,chr7)
names(chr7)[1]<-"list.rep"

str(chr7)
chr7[1:5,1:5]
dim(chr7)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr7 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr7<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr7[,1]
d.sum.chr7<-d.sum.chr7[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr7<-matrix(NA,nrow=nrow(d.sum.chr7),ncol=ncol(d.sum.chr7))

for(i in 1:nrow(call.chr7)){
  for(j in 1:ncol(call.chr7)){
    call.chr7[i,j]<-ifelse(d.sum.chr7[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr7,file="d.sum.chr7.csv")
#write.csv(call.chr7,file="call.chr7.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr7
tally.d_ab_freq<-tally.d_ab/d.sum.chr7
tally.d_bb_freq<-tally.d_bb/d.sum.chr7

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr7 <- matrix(NA,nrow=nrow(call.chr7),ncol=ncol(call.chr7))

for(i in 1:nrow(ave.chr.called.chr7)){
  for(j in 1:ncol(ave.chr.called.chr7)){
    ave.chr.called.chr7[i,j] <- ifelse(call.chr7[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



write.csv(ave.chr.called.chr7,file="ave.chr.called88.chr7.csv")


######chr8----
chromes$CE10_chr8[1:5,1:5]
chr8<-chromes$CE10_chr8[,-1]
dim(chr8)
colnames(chr8)[1:147]<-as.character(inds)
chr8[1:5,1:5]

#number of SNPs in your window
n=18

#make a list of 18 repeated however many times (number of SNPs)
list.rep<-list(rep(1:(nrow(chromes$CE10_chr8[,-1])%/%n+1),each=n,len=nrow(chromes$CE10_chr8[,-1])))
chr8<-cbind(list.rep,chr8)
names(chr8)[1]<-"list.rep"

str(chr8)
chr8[1:5,1:5]
dim(chr8)

#tally the number of genotype likelihood (0,1,2,NA) in each window for each individual
tally.d<-chr8 %>% # %>% is equivalent to | (pipe) in linux
  group_by(list.rep) %>%
  pivot_longer(names_to="ind",values_to="gl",cols=2:148) %>%
  group_by(list.rep,ind,gl) %>%
  summarise(counts=n())

#change back to wide form from long form above
tally.d2<-tally.d %>%
  pivot_wider(names_from=ind,values_from=counts)

#cleans up and adds total counts in each window and each GL class (0,1,2)
tally.d3<-aggregate(tally.d2,by=list(tally.d2$list.rep,tally.d2$gl),FUN=sum,na.rm=T)
str(tally.d3)

#write.csv(tally.d3,file="tally_d3.csv")
tally.d4<-tally.d3[,-c(1:4)]

#tally.df <- replace(tally.df, tally.df == 0, NA)
tally.d5<-cbind(tally.d3[,3:4],tally.d4)

names(tally.d5)[1:2]<-c("list.rep","gl")

str(tally.d5)

tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 0, "AA")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 1, "AB")
tally.d5$gl <- replace(tally.d5$gl, tally.d5$gl == 2, "BB")

#write.csv(tally.d5,file="tally_d5.csv")

list<-tally.d5$list
d.sum<-cbind(tally.d5$list,tally.d5[,3:149])
colnames(d.sum)[1]<-"list"

d.sum.chr8<-aggregate(.~list,d.sum,FUN='sum')
window<-d.sum.chr8[,1]
d.sum.chr8<-d.sum.chr8[,-1]


tally.d_aa<-filter(tally.d5,gl=="AA")
tally.d_aa<-tally.d_aa[,-c(1:2)]

tally.d_ab<-filter(tally.d5,gl=="AB")
tally.d_ab<-tally.d_ab[,-c(1:2)]

tally.d_bb<-filter(tally.d5,gl=="BB")
tally.d_bb<-tally.d_bb[,-c(1:2)]

call.chr8<-matrix(NA,nrow=nrow(d.sum.chr8),ncol=ncol(d.sum.chr8))

for(i in 1:nrow(call.chr8)){
  for(j in 1:ncol(call.chr8)){
    call.chr8[i,j]<-ifelse(d.sum.chr8[i,j] < 8, "NoCall","call")
  }
}

#write.csv(tally.d_aa,file="tally.d_aa.csv")
#write.csv(tally.d_ab,file="tally.d_ab.csv")
#write.csv(tally.d_bb,file="tally.d_bb.csv")
#write.csv(d.sum.chr8,file="d.sum.chr8.csv")
#write.csv(call.chr8,file="call.chr8.csv")

tally.d_aa_freq<-tally.d_aa/d.sum.chr8
tally.d_ab_freq<-tally.d_ab/d.sum.chr8
tally.d_bb_freq<-tally.d_bb/d.sum.chr8

#write.csv(tally.d_aa_freq,file="tally.d_aa_freq.csv")
#write.csv(tally.d_ab_freq,file="tally.d_ab_freq.csv")
#write.csv(tally.d_bb_freq,file="tally.d_bb_freq.csv")


ave.chr.called.chr8 <- matrix(NA,nrow=nrow(call.chr8),ncol=ncol(call.chr8))

for(i in 1:nrow(ave.chr.called.chr8)){
  for(j in 1:ncol(ave.chr.called.chr8)){
    ave.chr.called.chr8[i,j] <- ifelse(call.chr8[i,j]=="NoCall","NoCall",
                                       ifelse(tally.d_aa_freq[i,j] >= 0.88, "AA",
                                              ifelse(tally.d_bb_freq[i,j] >= 0.88, "BB",
                                                     ifelse(tally.d_ab_freq[i,j] >= 0.88, "AB","NoCall"))))
    
  }
}



#write.csv(ave.chr.called.chr8,file="ave.chr.called88.chr8.csv")

##Combine tallied datasets----
dim(ave.chr.called.chr8)

chrs.1<-as.data.frame(rep("CHR1",nrow(ave.chr.called.chr1)))
chrs.2<-as.data.frame(rep("CHR2",nrow(ave.chr.called.chr2)))
chrs.3<-as.data.frame(rep("CHR3",nrow(ave.chr.called.chr3)))
chrs.4<-as.data.frame(rep("CHR4",nrow(ave.chr.called.chr4)))
chrs.5<-as.data.frame(rep("CHR5",nrow(ave.chr.called.chr5)))
chrs.6<-as.data.frame(rep("CHR6",nrow(ave.chr.called.chr6)))
chrs.7<-as.data.frame(rep("CHR7",nrow(ave.chr.called.chr7)))
chrs.8<-as.data.frame(rep("CHR8",nrow(ave.chr.called.chr8)))

colnames(chrs.1)<-"chr"
colnames(chrs.2)<-"chr"
colnames(chrs.3)<-"chr"
colnames(chrs.4)<-"chr"
colnames(chrs.5)<-"chr"
colnames(chrs.6)<-"chr"
colnames(chrs.7)<-"chr"
colnames(chrs.8)<-"chr"

chrs<-rbind(chrs.1,chrs.2,chrs.3,chrs.4,chrs.5,chrs.6,chrs.7,chrs.8)

## Combine into one dataframe
tally.all.chrs<-rbind(ave.chr.called.chr1,ave.chr.called.chr2,ave.chr.called.chr3,ave.chr.called.chr4,
                      ave.chr.called.chr5,ave.chr.called.chr6,ave.chr.called.chr7,ave.chr.called.chr8)
dim(tally.all.chrs)

######Marker name formatting----
loci<-read.table("out.012.pos",header=F)
loci[1:5,]
loci$V3<-paste(loci$V1, loci$V2, sep="_")
loci[1:5,]
loci_comb<-loci$V3

#number of SNPs in your window
n=18

loci.bp<-split(loci,loci$V1)

#chr1 names
max.bp.chr1<-aggregate(loci.bp$CE10_chr1[,-1],list(rep(1:(nrow(loci.bp$CE10_chr1[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr1[,-1]))),max,na.rm=T)[-1]
min.bp.chr1<-aggregate(loci.bp$CE10_chr1[,-1],list(rep(1:(nrow(loci.bp$CE10_chr1[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr1[,-1]))),min,na.rm=T)[-1]

mean.bp.chr1<-aggregate(loci.bp$CE10_chr1[,-1],list(rep(1:(nrow(loci.bp$CE10_chr1[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr1[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr1<-cbind(min.bp.chr1$V2,max.bp.chr1$V2,mean.bp.chr1$V2)

colnames(minMaxMean.bp.chr1)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr1<-as.data.frame(minMaxMean.bp.chr1)
minMaxMean.bp.chr1$comb<-paste("chr1",minMaxMean.bp.chr1$minBP,minMaxMean.bp.chr1$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr1))

minMaxMean.bp.chr1$marker.n<-paste(minMaxMean.bp.chr1$comb,seq,sep="_")

#chr2 names
max.bp.chr2<-aggregate(loci.bp$CE10_chr2[,-1],list(rep(1:(nrow(loci.bp$CE10_chr2[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr2[,-1]))),max,na.rm=T)[-1]
min.bp.chr2<-aggregate(loci.bp$CE10_chr2[,-1],list(rep(1:(nrow(loci.bp$CE10_chr2[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr2[,-1]))),min,na.rm=T)[-1]

mean.bp.chr2<-aggregate(loci.bp$CE10_chr2[,-1],list(rep(1:(nrow(loci.bp$CE10_chr2[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr2[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr2<-cbind(min.bp.chr2$V2,max.bp.chr2$V2,mean.bp.chr2$V2)

colnames(minMaxMean.bp.chr2)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr2<-as.data.frame(minMaxMean.bp.chr2)
minMaxMean.bp.chr2$comb<-paste("chr2",minMaxMean.bp.chr2$minBP,minMaxMean.bp.chr2$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr2))

minMaxMean.bp.chr2$marker.n<-paste(minMaxMean.bp.chr2$comb,seq,sep="_")


#chr3 names
max.bp.chr3<-aggregate(loci.bp$CE10_chr3[,-1],list(rep(1:(nrow(loci.bp$CE10_chr3[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr3[,-1]))),max,na.rm=T)[-1]
min.bp.chr3<-aggregate(loci.bp$CE10_chr3[,-1],list(rep(1:(nrow(loci.bp$CE10_chr3[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr3[,-1]))),min,na.rm=T)[-1]

mean.bp.chr3<-aggregate(loci.bp$CE10_chr3[,-1],list(rep(1:(nrow(loci.bp$CE10_chr3[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr3[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr3<-cbind(min.bp.chr3$V2,max.bp.chr3$V2,mean.bp.chr3$V2)

colnames(minMaxMean.bp.chr3)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr3<-as.data.frame(minMaxMean.bp.chr3)
minMaxMean.bp.chr3$comb<-paste("chr3",minMaxMean.bp.chr3$minBP,minMaxMean.bp.chr3$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr3))

minMaxMean.bp.chr3$marker.n<-paste(minMaxMean.bp.chr3$comb,seq,sep="_")


#chr4 names
max.bp.chr4<-aggregate(loci.bp$CE10_chr4[,-1],list(rep(1:(nrow(loci.bp$CE10_chr4[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr4[,-1]))),max,na.rm=T)[-1]
min.bp.chr4<-aggregate(loci.bp$CE10_chr4[,-1],list(rep(1:(nrow(loci.bp$CE10_chr4[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr4[,-1]))),min,na.rm=T)[-1]

mean.bp.chr4<-aggregate(loci.bp$CE10_chr4[,-1],list(rep(1:(nrow(loci.bp$CE10_chr4[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr4[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr4<-cbind(min.bp.chr4$V2,max.bp.chr4$V2,mean.bp.chr4$V2)

colnames(minMaxMean.bp.chr4)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr4<-as.data.frame(minMaxMean.bp.chr4)
minMaxMean.bp.chr4$comb<-paste("chr4",minMaxMean.bp.chr4$minBP,minMaxMean.bp.chr4$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr4))

minMaxMean.bp.chr4$marker.n<-paste(minMaxMean.bp.chr4$comb,seq,sep="_")


#chr5 names
max.bp.chr5<-aggregate(loci.bp$CE10_chr5[,-1],list(rep(1:(nrow(loci.bp$CE10_chr5[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr5[,-1]))),max,na.rm=T)[-1]
min.bp.chr5<-aggregate(loci.bp$CE10_chr5[,-1],list(rep(1:(nrow(loci.bp$CE10_chr5[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr5[,-1]))),min,na.rm=T)[-1]

mean.bp.chr5<-aggregate(loci.bp$CE10_chr5[,-1],list(rep(1:(nrow(loci.bp$CE10_chr5[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr5[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr5<-cbind(min.bp.chr5$V2,max.bp.chr5$V2,mean.bp.chr5$V2)

colnames(minMaxMean.bp.chr5)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr5<-as.data.frame(minMaxMean.bp.chr5)
minMaxMean.bp.chr5$comb<-paste("chr5",minMaxMean.bp.chr5$minBP,minMaxMean.bp.chr5$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr5))

minMaxMean.bp.chr5$marker.n<-paste(minMaxMean.bp.chr5$comb,seq,sep="_")


#chr6 names
max.bp.chr6<-aggregate(loci.bp$CE10_chr6[,-1],list(rep(1:(nrow(loci.bp$CE10_chr6[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr6[,-1]))),max,na.rm=T)[-1]
min.bp.chr6<-aggregate(loci.bp$CE10_chr6[,-1],list(rep(1:(nrow(loci.bp$CE10_chr6[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr6[,-1]))),min,na.rm=T)[-1]

mean.bp.chr6<-aggregate(loci.bp$CE10_chr6[,-1],list(rep(1:(nrow(loci.bp$CE10_chr6[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr6[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr6<-cbind(min.bp.chr6$V2,max.bp.chr6$V2,mean.bp.chr6$V2)

colnames(minMaxMean.bp.chr6)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr6<-as.data.frame(minMaxMean.bp.chr6)
minMaxMean.bp.chr6$comb<-paste("chr6",minMaxMean.bp.chr6$minBP,minMaxMean.bp.chr6$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr6))

minMaxMean.bp.chr6$marker.n<-paste(minMaxMean.bp.chr6$comb,seq,sep="_")


#chr7 names
max.bp.chr7<-aggregate(loci.bp$CE10_chr7[,-1],list(rep(1:(nrow(loci.bp$CE10_chr7[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr7[,-1]))),max,na.rm=T)[-1]
min.bp.chr7<-aggregate(loci.bp$CE10_chr7[,-1],list(rep(1:(nrow(loci.bp$CE10_chr7[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr7[,-1]))),min,na.rm=T)[-1]

mean.bp.chr7<-aggregate(loci.bp$CE10_chr7[,-1],list(rep(1:(nrow(loci.bp$CE10_chr7[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr7[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr7<-cbind(min.bp.chr7$V2,max.bp.chr7$V2,mean.bp.chr7$V2)

colnames(minMaxMean.bp.chr7)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr7<-as.data.frame(minMaxMean.bp.chr7)
minMaxMean.bp.chr7$comb<-paste("chr7",minMaxMean.bp.chr7$minBP,minMaxMean.bp.chr7$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr7))

minMaxMean.bp.chr7$marker.n<-paste(minMaxMean.bp.chr7$comb,seq,sep="_")

#chr8 names
max.bp.chr8<-aggregate(loci.bp$CE10_chr8[,-1],list(rep(1:(nrow(loci.bp$CE10_chr8[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr8[,-1]))),max,na.rm=T)[-1]
min.bp.chr8<-aggregate(loci.bp$CE10_chr8[,-1],list(rep(1:(nrow(loci.bp$CE10_chr8[,-1])%/%n+1),
                                                       each=n,len=nrow(loci.bp$CE10_chr8[,-1]))),min,na.rm=T)[-1]

mean.bp.chr8<-aggregate(loci.bp$CE10_chr8[,-1],list(rep(1:(nrow(loci.bp$CE10_chr8[,-1])%/%n+1),
                                                        each=n,len=nrow(loci.bp$CE10_chr8[,-1]))),mean,na.rm=T)[-1]

minMaxMean.bp.chr8<-cbind(min.bp.chr8$V2,max.bp.chr8$V2,mean.bp.chr8$V2)

colnames(minMaxMean.bp.chr8)<-c("minBP","maxBP","meanBP")

minMaxMean.bp.chr8<-as.data.frame(minMaxMean.bp.chr8)
minMaxMean.bp.chr8$comb<-paste("chr8",minMaxMean.bp.chr8$minBP,minMaxMean.bp.chr8$maxBP,sep="_")
seq<-seq(from=1,to=nrow(minMaxMean.bp.chr8))

minMaxMean.bp.chr8$marker.n<-paste(minMaxMean.bp.chr8$comb,seq,sep="_")

## combine into one dataframe
marker.names<-rbind(minMaxMean.bp.chr1,minMaxMean.bp.chr2,minMaxMean.bp.chr3,minMaxMean.bp.chr4,
                    minMaxMean.bp.chr5,minMaxMean.bp.chr6,minMaxMean.bp.chr7,minMaxMean.bp.chr8)
dim(marker.names)

write.csv(marker.names,file="marker_names_info.csv")

tally.chr.called<-cbind(chrs,marker.names$marker.n,tally.all.chrs)
tally.chr.called[1:5,1:5]

tally.chr.called.final<-tally.chr.called[,-c(1)]
tally.chr.called.final[1:5,1:5]

dim(tally.chr.called.final)

inds<-colnames(tally.d2)
inds<-inds[-c(1:2)]

colnames(tally.chr.called.final)[2:148]<-inds
colnames(tally.chr.called.final)[1]<-"marker_name"
tally.chr.called.final[1:10,1:5]

# checking no calls or miscalls in par and card
#YY1.01A is card
which(tally.chr.called.final$YY1.01A == "BB")
which(tally.chr.called.final$YY1.01A == "AB")

#YY1.12H is par
which(tally.chr.called.final$YY1.12H == "AA")
het.call<-which(tally.chr.called.final$YY1.12H == "AB")
tally.chr.called.final$YY1.12H[het.call]<-"BB"


nocall.rep.12h<-which(tally.chr.called.final$YY1.12H == "NoCall")
nocall.rep.01a<-which(tally.chr.called.final$YY1.01A == "NoCall")

nocall.rep.12h
nocall.rep.01a

tally.chr.called.final$YY1.12H[nocall.rep.12h]<-"BB"
tally.chr.called.final$YY1.01A[nocall.rep.01a]<-"AA"

# write out final datafile
write.table(tally.chr.called.final,file="tally_binnedSNPsCalled_88.txt",sep="\t",quote=F,row.names=F)
write.csv(tally.chr.called.final,file="tally_binnedSNPsCalled_88.csv")

##Filter windows and inds----
rm(list=setdiff(ls(), "tally.chr.called.final"))
str(tally.chr.called.final)

#get rid of parents (last two columns)
tally.chr.called.final.noPC<-tally.chr.called.final[,-c(147:148)]


ind.freq<-matrix(NA,nrow=145,ncol=4)
colnames(ind.freq)<-c("freq.AA","freq.AB","freq.BB","freq.NoCall")
inds<-colnames(tally.chr.called.final.noPC)
inds<-inds[-1]

ind.freq<-cbind(inds,ind.freq)

for(i in 2:ncol(tally.chr.called.final.noPC)){
  ind.freq[i-1,2]<-sum(tally.chr.called.final.noPC[,i]=="AA")/(nrow(tally.chr.called.final.noPC))
  ind.freq[i-1,3]<-sum(tally.chr.called.final.noPC[,i]=="AB")/(nrow(tally.chr.called.final.noPC))
  ind.freq[i-1,4]<-sum(tally.chr.called.final.noPC[,i]=="BB")/(nrow(tally.chr.called.final.noPC))
  ind.freq[i-1,5]<-sum(tally.chr.called.final.noPC[,i]=="NoCall")/(nrow(tally.chr.called.final.noPC))
}

ind.freq<-as.data.frame(ind.freq)

filt.n<-which(ind.freq$freq.NoCall > 0.5)

# get rid of those inds w/ >.5 freq of NoCall
tally.chr.called.final.filt<-tally.chr.called.final.noPC[,-(filt.n+1)]

# now onto windows
window.freq<-matrix(NA,nrow=nrow(tally.chr.called.final.filt),ncol=4)
colnames(window.freq)<-c("freq.AA","freq.AB","freq.BB","freq.NoCall")

window.freq<-cbind(tally.chr.called.final.filt$rounded_bp,window.freq)

sum(tally.chr.called.final[1,]=="AA")/145

for(i in 1:nrow(tally.chr.called.final)){
  window.freq[i,2]<-sum(tally.chr.called.final[i,]=="AA")/145
  window.freq[i,3]<-sum(tally.chr.called.final[i,]=="AB")/145
  window.freq[i,4]<-sum(tally.chr.called.final[i,]=="BB")/145
  window.freq[i,5]<-sum(tally.chr.called.final[i,]=="NoCall")/145
}

window.freq<-as.data.frame(window.freq)

filt.win<-which(window.freq$freq.NoCall > 0.5)

# get rid of those inds
tally.chr.called.final.filt<-tally.chr.called.final.filt[-(filt.win),]

# need to add back in the par and card inds
# filter windows from them

parents<-tally.chr.called.final[,c(147:148)]
parents<-parents[-(filt.win),]


tally.chr.called.final.filt<-cbind(tally.chr.called.final.filt,parents)

# write out final filtered datafile
write.table(tally.chr.called.final.filt,file="tally_binnedSNPsCalled_88_filt.txt",sep="\t",quote=F,row.names=F)
write.csv(tally.chr.called.final.filt,file="tally_binnedSNPsCalled_88_filt.csv")

##Window Size Stats----
d<-read.csv("Window_markers_Info.csv",header=T)

str(d)
hist(d$window_size,breaks=100)
max(d$window_size)
min(d$window_size) #end of chr
summary(d$window_size)


#Post LEPMAP filtering----
setwd("~/Dropbox/Research/Georgia/UGA_Research/Mimulus/RIL_Computing/LinkageMapping/LepMap_pubRun/Windowed_run/LepMap/fixed_map")
lep<-read.csv("lepmap_gtype_matrix.csv",header=T,stringsAsFactors = T)
og<-read.csv("tally_binnedSNPsCalled_88_lepmap.csv",header=T,stringsAsFactors = T)

which(lep$marker_id != og$marker_name)

lep.info<-lep[,c(1:3)]
lep<-lep[,-c(1:3)]

og<-og[,-c(1:3)]

lep[1:5,1:5]
og[1:5,1:5]

lep<-as.matrix(lep)
og<-as.matrix(og)

new.matrix<-matrix(NA,nrow=nrow(lep),ncol=ncol(lep))

##little funky w/ first and last marker if they are H, but can easily edit those manually after the fact
##might fill those in w/ NAs

for(i in 1:nrow(new.matrix)){
  for(j in 1:ncol(new.matrix)){
    new.matrix[i,j]<-ifelse(lep[i,j]=="A","A",
                                 ifelse(lep[i,j]=="B","B",
                                        ifelse(lep[i,j] =="H" & og[i,j] =="AB","H",
                                               ifelse(lep[i,j] == "H" & lep[i-1,j] == "H", "H",
                                                      ifelse(lep[i,j] == "H" & lep[i+1,j] == "A" | lep[i+1,j] == "B","NoCall","H")))))
  }
}


yes<-which(lep == new.matrix)
no<-which(lep != new.matrix)


length(yes)/(length(yes)+length(no))
length(no)/(length(yes)+length(no))

final.matrix<-cbind(lep.info,new.matrix)
write.csv(final.matrix,file="finalFixed_gtype_matrix.csv")



