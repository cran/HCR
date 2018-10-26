# HCR

## Overview

This code provides a method to fit the hidden compact representation model as well as to discover the causal direction on discrete data. 

Please cite "Ruichu Cai, Jie Qiao, Kun Zhang, Zhenjie Zhang, Zhifeng Hao. Causal Discovery from Discrete Data using Hidden Compact Representation. NIPS 2018." 
## Installation

```
install.packages("HCR")
```

### Quick Start

This package contains the data synthetic process for HCR model. Here are some examples to make a quick start:

```
# HCR data
library(data.table)
set.seed(1)
data=simuXY(sample_size=2000)
r1<-HCR(data$X,data$Y,score_type = "aic")
r2<-HCR(data$Y,data$X,score_type = "aic")
if(r1$score>r2$score){
  print("X cause Y.")
}else{
  print("Y cause X.")
}
unique(r1$data[,c("X","Yp")]) # The fitted hidden representation
unique(data.frame(data$X,data$Yp)) # The true hidden representation

# A fast implementation of HCR
r1<-HCR.fast(data$X,data$Y)
r2<-HCR.fast(data$Y,data$X)
if(r1$score>r2$score){
  print("X cause Y.")
}else{
  print("Y cause X.")
}
unique(r1$data[,c("X","Yp")]) # The fitted hidden representation
unique(data.frame(data$X,data$Yp)) # The true hidden representation

# ANM data
# Note that the complete generate of the ANM data can be found in http://webdav.tuebingen.mpg.de/causality/
set.seed(0)
X=sample(1:5,2000,replace = T)
N=sample(-1:1,2000,replace = T)
Y=X+N
r1=HCR(X,Y,score_type = "bic",is_anm = T)
r2=HCR(Y,X,score_type = "bic",is_anm = T)
if(r1$score>r2$score){
  print("X cause Y.")
}else{
  print("Y cause X.")
}


```

