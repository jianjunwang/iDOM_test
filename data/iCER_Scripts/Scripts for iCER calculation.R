

###################################################
#### Scripts uploading to Github #####

### an example for calculating iCTR


##### clean the environment----
rm(list=ls())

##### set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
getwd()



##### load library
library(dplyr)  
library(tidyverse)
library(reshape2)
library(plyr)
library(vegan)



##### load DOM data 
data.ra = read.csv("iCER_Scripts/data_table.csv", header = T, row.names = 1, check.names = F)

## load meta data
envi = read.csv("iCER_Scripts/envi.csv", header = T, row.names = 1, check.names = F)

## temperature
vars.list = "Temperature" 
envi.tmp =  envi[, vars.list, drop = FALSE]



######################################################################################################################
######################################################################################################################
####(1) Generating independent data sets for MER and iCER calculation (i.e., "MER dataset" and "iCER dataset") 

## Function for data splitting method to generate independent data sets
random_ID_select = function(dat, prop, ele_num, end_n) {
  
  # dat: meta data including sample ID and temperature
  # prop: split ratio, that is the proportion of samples for MER calculation, such as 0.80
  # ele_num: the number of temperature levels used for MER calculation, such as 5 
  # end_n: the randomization times at least for a given sample for iCER calculation, such as 50 
  
  ############################
  #### 
  set.seed(12345)
  num <-  round(nrow(dat) * prop)            
  ele_gra_ser <- sort(unique(dat$Elevation))  
  
  samp.list <- list()
  unselect.list <- list()
  ID_num <- matrix(nrow = num)
  while (TRUE) {    
    ele_gra <- sort(sample(ele_gra_ser, ele_num, replace = FALSE))    
    random_ID <- sample(which(dat$Elevation %in% ele_gra), num, replace = FALSE)
    df_check <- dat[random_ID, ]
    while (!all(ele_gra %in% df_check$Elevation)) {
      random_ID <- sample(which(dat$Elevation %in% ele_gra), num, replace = FALSE)
      df_check <- dat[random_ID, ]
    }
    random_ID <- sort(random_ID)
    if (any(apply(ID_num, 2, function(dup_det) identical(dup_det, random_ID)))) {
    } else {
      samp.list <- c(samp.list, list(dat[random_ID, ]))
      unselect.list <- c(unselect.list, list(dat[-random_ID, 1]))
    }
    ID_num <- cbind(ID_num, random_ID)
    df <- do.call(cbind, unselect.list)
    if (all(sapply(dat$ID, function(x) sum(x == df)) >= end_n)) {
      break 
    }
  }
  ID_num <- ID_num[,-1]
  
  return(list(samp.list, unselect.list))
  
}  ## end of function


dat = tibble(ID = envi$cDOM.IDs, Elevation = envi$Temperature)
sort(unique(dat$Elevation))
rand.id.tmp = random_ID_select(dat, prop=0.80, ele_num=7, end_n=10)
samp.rand.list = rand.id.tmp[[1]]
unselect.list = rand.id.tmp[[2]]



########################################################################################################
########################################################################################################
####(2) Calculating MER using "MER dataset"

cutoff = 0.3
data.ra.go = data.ra
samp.tmp = samp.rand.list

rho.out.tmp = data.frame()
peak.rand.list = c()
kk=1   
ii=1   
for (kk in 1:length(samp.tmp)) {      # randomization times
  samp.list.go = samp.tmp[[kk]]
  samp.list.go = as.character(samp.list.go$ID)
  
  data.ra.go2 = data.ra.go[samp.list.go, ]
  data.ra.go2 = data.ra.go2[, colSums(data.ra.go2) > 0]
  data.pa.go2 = decostand(data.ra.go2, method = "pa")
  data.ra.go2.keep = data.ra.go2[, colSums(data.pa.go2) >= nrow(data.pa.go2) * cutoff]
  
  ## combining temperature
  dat.tmp = merge(envi.tmp, data.ra.go2.keep,
                  by.x="row.names", by.y="row.names", all.y = TRUE)
  dat.tmp.long = reshape2::melt(dat.tmp, id = c(colnames(dat.tmp)[c(1:(1+length(vars.list)))]))
  
  ## peak
  peak.go = colnames(data.ra.go2.keep)
  
  for (ii in 1:length(peak.go)) {     # peak number
    
    dat.tmp.long.i = subset(dat.tmp.long, variable %in% peak.go[ii])
    
    if (nrow(dat.tmp.long.i) < 5) next
    if (nrow(dat.tmp.long.i) >= 5) {
      
      tmp = dat.tmp.long.i[, c("Temperature", "value")]
      corr.out.tmp = data.frame(spearman.cor = with(tmp, cor(value, Temperature, method="spearman")), 
                                spearman.cor.p = with(tmp, cor.test(value, Temperature, method = "spearman"))$p.value,
                                Permutation = kk,
                                peak = peak.go[ii]
      )
      
      ## collect MERs
      rho.out.tmp = rbind(rho.out.tmp, corr.out.tmp)
    }
  }
  
  ## collect peaks
  peak.go.tmp = list(peak.go); names(peak.go.tmp) = kk
  peak.rand.list = c(peak.rand.list, peak.go.tmp)
}
dim(rho.out.tmp); head(rho.out.tmp,2)
str(peak.rand.list)



########################################################################################################
########################################################################################################
####(3) Calculating iCER using "iCER dataset"

##########################
#### calculating iCER using ALL molecules with MERs

CER.all.out = data.frame()
k=1
for (k in 1:length(peak.rand.list)) {   
  
  rho.peak.go0 = subset(rho.out.tmp, Permutation == k)   # MTRs from "MER dataset"
  rho.peak.go = rho.peak.go0[, c("spearman.cor", "spearman.cor.p", "peak")]  
  
  peak.permu.go = as.character(rho.peak.go$peak)
  
  unselect.permu.go = unselect.list[[k]]
  unselect.permu.go = as.character(unselect.permu.go$ID)
  
  data.ra.go = data.ra[unselect.permu.go, ]  # Relative abundance of molecules in each sample ("iCER dataset")
  data.ra.go = data.ra.go[, colSums(data.ra.go) > 0]
  
  peak.comm.list = intersect(peak.permu.go, colnames(data.ra.go))
  
  rho.go = subset(rho.peak.go, peak %in% peak.comm.list)  
  rownames(rho.go) = rho.go$peak; rho.go = rho.go[, "spearman.cor", drop=FALSE]
  
  ra.go = data.ra.go[, rownames(rho.go)]   
  ra.go = ra.go[rowSums(ra.go) > 0, ]
  CER.tmp = FD::functcomp(rho.go, as.matrix(ra.go))   # iCER calculation
  CER.tmp2 = data.frame(CER.tmp, cDOM.IDs = rownames(CER.tmp))
  
  ## collect iCERs
  CER.all.out.tmp = data.frame(CER.tmp2, Permutation = k)
  CER.all.out = rbind(CER.all.out, CER.all.out.tmp)
}
dim(CER.all.out); head(CER.all.out, 2)



##########################
#### calculating iCER using the molecules with significant MERs

CER.sig.out = data.frame()
k=1
for (k in 1:length(peak.rand.list)) {   
  
  rho.peak.go0 = subset(rho.out.tmp, Permutation == k & spearman.cor.p <= 0.05)   # MTRs from "MER dataset"
  rho.peak.go = rho.peak.go0[, c("spearman.cor", "spearman.cor.p", "peak")]  
  
  peak.permu.go = as.character(rho.peak.go$peak)
  
  unselect.permu.go = unselect.list[[k]]
  unselect.permu.go = as.character(unselect.permu.go$ID)
  
  data.ra.go = data.ra[unselect.permu.go, ]   # Relative abundance of molecules in each sample ("iCER dataset")
  data.ra.go = data.ra.go[, colSums(data.ra.go) > 0]
  
  peak.comm.list = intersect(peak.permu.go, colnames(data.ra.go))
  
  rho.go = subset(rho.peak.go, peak %in% peak.comm.list)  
  rownames(rho.go) = rho.go$peak; rho.go = rho.go[, "spearman.cor", drop=FALSE]
  
  ra.go = data.ra.go[, rownames(rho.go)]   
  ra.go = ra.go[rowSums(ra.go) > 0, ]
  
  CER.tmp = FD::functcomp(rho.go, as.matrix(ra.go))   # iCER calculation
  CER.tmp2 = data.frame(CER.tmp, cDOM.IDs = rownames(CER.tmp))
  
  ## collect iCERs
  CER.sig.out.tmp = data.frame(CER.tmp2, Permutation = k)
  CER.sig.out = rbind(CER.sig.out, CER.sig.out.tmp)
}
dim(CER.sig.out); head(CER.sig.out, 2)

## mean iCERs for each sample
CER.mean = ddply(CER.sig.out, .(cDOM.IDs), summarise, iCER=mean(spearman.cor))
dim(CER.mean); head(CER.mean, 2)



##########################
### plotting iCER against temperature
dat_plot = merge(CER.mean, envi, by="cDOM.IDs", all=TRUE)
dim(dat_plot); head(dat_plot,2)

plot.iCER.temp =
  ggplot()+
  geom_point(data = dat_plot, mapping=aes(x=Temperature, y=iCER), color="black", size=1.5)+
  geom_smooth(data = dat_plot, mapping=aes(x=Temperature, y=iCER), color="red",
              method="lm",formula=y~x,size=0.9)+  
  xlab("Temperature (°C)")+ylab("iCER")+
  theme_bw()+theme(legend.title=element_blank())+
  theme(legend.position="none")+
  theme(panel.grid.major = element_line(size=0.05)) +
  theme(panel.grid.minor = element_line(size=0.05)) +
  theme(axis.title=element_text(size=9, color="black"))+
  theme(axis.text=element_text(size=8, color="black"))
plot.iCER.temp





























