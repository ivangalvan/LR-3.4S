## The simulated data contains 25,000 SNPs and 1000 individuals (100 PO, FS, 3/4S, 2ND and UN pairs)

setwd("C:/Users/ivangalvan/Desktop/three_quarter/simulations_family_relationships/GitHub/")

library(data.table)
library(dplyr)

source("LR_IBD.R")

# make (k0,k1)-plot with PLINK

system("plink2 --bfile data/example --genome --out example")

genome = fread("example.genome")

np= 100
pairs1 = c(paste0("un_",1:np),
           paste0("second_",1:np),
           paste0("3_4_S_",1:np),
           paste0("fs_",1:np),
           paste0("po_",1:np))

pairs2 = paste0(pairs1,"_bis")

pairs = as.data.frame(cbind(pairs1,pairs2))

# select a random subset of 20 FS, 3/4S and 2ND degree pairs

set.seed(123)

pairs = pairs[c(sample(101:200,20),
                sample(201:300,20),
                sample(301:400,20)),]

pairs$Z0 = 0
pairs$Z1 = 0
pairs$Z2 = 0

for(i in 1:nrow(pairs)){
  aux = genome %>% 
    filter((IID1 %in% pairs$pairs1[i] & IID2 %in% pairs$pairs2[i]))
  
  pairs[i,3:5] = aux[,c("Z0","Z1","Z2")]
  
  print(i)
}

my_cols = c(rep("darkorchid2",20),
            rep("black",20),
            rep("dodgerblue2",20))

pairs$cols = my_cols

dir.create("plots")

png("plots/k0_k1_plot.png",width = 6,height = 6,units = "in",res = 300)
plot(pairs$Z0,pairs$Z1,col=my_cols,pch=16,xlab=expression(hat(k)[0]),ylab="",asp=1,xlim=0:1,ylim=0:1)
legend("topright",
       c("2nd","3/4S","FS"), col=c("darkorchid2","black",
                                             "dodgerblue2"),pch=16)
mtext(side=2, text=expression(hat(k)[1]), line=2)
dev.off()

pairs$pairs1 = as.character(pairs$pairs1)
pairs$pairs2 = as.character(pairs$pairs2)


##  Calculate minor allele frequency ########

system("plink2 --bfile data/example --freq --out data/example")

freq = fread("data/example.frq")

system("plink2 --bfile data/example --recode A-transpose --out data/example")

geno_data = fread("data/example.traw")

geno_data[1:10,1:10]

geno_data = geno_data[,-c(1:6)]

fam_file = fread("data/example.fam")

## Calculate likelihood ratios #########

# FS vs UN

pairs$fs_un = 0

for(j in 1:nrow(pairs)){
  
  sample1 = as.numeric(t(geno_data[,which(fam_file$V2 %in% pairs$pairs1[j]),with=F]))
  sample2 = as.numeric(t(geno_data[,which(fam_file$V2 %in% pairs$pairs2[j]),with=F]))
  
  pairs$fs_un[j] = LR_IBD(sample1,sample2,freq$MAF,
                              rel.num = c(0.25,0.5,0.25),rel.den=c(1,0,0),e=0.001)
  
  print(j)
}


# 3/4S vs UN

pairs$three4s_un = 0

for(j in 1:nrow(pairs)){
  
  sample1 = as.numeric(t(geno_data[,which(fam_file$V2 %in% pairs$pairs1[j]),with=F]))
  sample2 = as.numeric(t(geno_data[,which(fam_file$V2 %in% pairs$pairs2[j]),with=F]))
  
  pairs$three4s_un[j] = LR_IBD(sample1,sample2,freq$MAF,
                                   rel.num = c(3/8,1/2,1/8),rel.den=c(1,0,0),e=0.001)
  
  print(j)
}

# 2nd vs UN

pairs$second_un = 0

for(j in 1:nrow(pairs)){
  
  sample1 = as.numeric(t(geno_data[,which(fam_file$V2 %in% pairs$pairs1[j]),with=F]))
  sample2 = as.numeric(t(geno_data[,which(fam_file$V2 %in% pairs$pairs2[j]),with=F]))
  
  pairs$second_un[j] = LR_IBD(sample1,sample2,freq$MAF,
                                  rel.num = c(0.5,0.5,0),rel.den=c(1,0,0),e=0.001)
  
  print(j)
}

### LR plot ###########

n_pairs = nrow(pairs)

png("plots/LR_FS_2nd_34S.png",width = 8,height = 8,units = "in",res = 300)

y_fs = jitter(rep(1,n_pairs))
y_34s = jitter(rep(0,n_pairs))
y_second = jitter(rep(-1,n_pairs))

plot(pairs$fs_un,y_fs + 10,ylim = c(-1,1),
     xlab="log10 LikelihoodRatio",ylab="",yaxt='n',
     xlim = c(min(pairs[,7:9]),max(pairs[,7:9])))

for(i in 21:40){
  lines(rbind(c(pairs$fs_un[i],y_fs[i]),c(pairs$three4s_un[i],y_34s[i])),lty=2,col="gray60")
  lines(rbind(c(pairs$three4s_un[i],y_34s[i]),c(pairs$second_un[i],y_second[i])),lty=2,col="gray60")
}

for(i in 1:20){
  lines(rbind(c(pairs$fs_un[i],y_fs[i]),c(pairs$three4s_un[i],y_34s[i])),lty=2,col=cm.colors(5)[5])
  lines(rbind(c(pairs$three4s_un[i],y_34s[i]),c(pairs$second_un[i],y_second[i])),lty=2,col=cm.colors(5)[5])
}

for(i in 41:60){
  lines(rbind(c(pairs$fs_un[i],y_fs[i]),c(pairs$three4s_un[i],y_34s[i])),lty=2,col="deepskyblue")
  lines(rbind(c(pairs$three4s_un[i],y_34s[i]),c(pairs$second_un[i],y_second[i])),lty=2,col="deepskyblue")
}


labs = c("FS~UN","3/4S~UN","2nd~UN")

points(pairs$fs_un,y_fs,col=pairs$cols)
points(pairs$three4s_un,y_34s,col=pairs$cols)
points(pairs$second_un,y_second,col=pairs$cols)

axis(side = 2,at = c(1:(-1)),labels = c("","",""),las=3,cex.axis=1)

legend("bottomright",c("2nd","3/4S","FS"),col=c("darkorchid2","black","dodgerblue2"),pch=16)

text(cex=1, x=-0.021, y=c(1:(-1)), labs, xpd=TRUE, srt=0)

dev.off()
