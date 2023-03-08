--- 
title: TFXcan preliminary data 
author: Haky Im 
date: '2023-01-19' 
slug: tfxcan-preliminary-data 
categories: 
  - analysis
tags: [] 
---

Based on Temi's (FOXA1) TF prediction from DNAseq

```r
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
PRE = "/Users/haekyungim/Library/CloudStorage/Box-Box/LargeFiles/imlab-data/data-Github/web-data"
##PRE="/Users/margaretperry/Library/CloudStorage/Box-Box/imlab-data/data-Github/web-data "
##PRE="/Users/temi/Library/CloudStorage/Box-Box/imlab-data/data-Github/web-data"
## COPY THE DATE AND SLUG fields FROM THE HEADER
SLUG="tfxcan-preliminary-data" ## copy the slug from the header
bDATE='2023-01-19' ## copy the date from the blog's header here
DATA = glue("{PRE}/{bDATE}-{SLUG}")
if(!file.exists(DATA)) system(glue::glue("mkdir {DATA}"))
WORK=DATA

## move data to DATA
#tempodata=("~/Downloads/tempo/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv")
#system(glue::glue("cp {tempodata} {DATA}/"))
system(glue("open {DATA}")) ## this will open the folder 
```

Read Temi's data

======= \## system(glue(\"open {DATA}\")) \## this will open the folder

Read Temi's data

```r
##source qqunif gist
devtools::source_gist("38431b74c6c0bf90c12f")
```

    ## ℹ Sourcing gist "38431b74c6c0bf90c12f"
    ## ℹ SHA-1 hash of file is "c5aba9ddce06b95125b727d96bffe7bd1557fcc3"

\>\>\>\>\>\>\> ee1d5e03a5e1098ffb92aa0a139a5ee3bceb2afe

```r
tempo <- readRDS(glue("{DATA}/freedman_FOXA1_predicted_vs_observed_matrix_2023-01-19.rds"))
head(tempo$observed)[,1:5]
```

    ##                     region LuCaP_173 LuCaP_167 LuCaP_93 LuCaP_81
    ## 1 chr1_100004103_100004112         0         0        0        0
    ## 2 chr1_100047865_100047874         0         0        0        0
    ## 3 chr1_100057135_100057144         0         0        0        0
    ## 4 chr1_100118800_100118809         0         0        0        0
    ## 5 chr1_100130867_100130876         0         1        0        0
    ## 6 chr1_100162937_100162946         0         0        0        0

```r
head(tempo$predicted)[,1:5]
```

    ##                     region LuCaP_173 LuCaP_167  LuCaP_93  LuCaP_81
    ## 1 chr1_100004103_100004112 0.2428343 0.2438005 0.2439815 0.2422898
    ## 2 chr1_100047865_100047874 0.1331662 0.1348004 0.1343903 0.1313637
    ## 3 chr1_100057135_100057144 0.2555964 0.2728396 0.2716254 0.2644257
    ## 4 chr1_100118800_100118809 0.1309184 0.1311011 0.1311917 0.1304301
    ## 5 chr1_100130867_100130876 0.6093923 0.6089561 0.6101078 0.6093229
    ## 6 chr1_100162937_100162946 0.1385310 0.1386486 0.1376289 0.1395273

```r
identical(tempo$observed$region, tempo$predicted$region)
```

    ## [1] TRUE

```r
obs <- tempo$observed %>% select(-region) %>% as.matrix
rownames(obs) <- tempo$observed$region

pred <- tempo$predicted %>% select(-region) %>% as.matrix
rownames(pred) <- tempo$predicted$region

identical(rownames(obs), rownames(pred))
```

    ## [1] TRUE

```r
## remove sites with no binding
nbou <- apply(obs,1,sum)
ind <- nbou>3 & nbou < 14
obs <- obs[ind,]
pred <- pred[ind,]
identical(rownames(obs), rownames(pred))
```

    ## [1] TRUE

```r
pvec = rep(NA,nrow(obs))
tvec = rep(NA,nrow(obs))
for(cc in 1:nrow(obs) ) 
{
  res = t.test(pred[cc,obs[cc,]==1], pred[cc,obs[cc,]==0])
  pvec[cc] = res$p.value
  tvec[cc] = res$statistic
  }


pred_bou <- pred
pred_unbou <- pred
pred_bou[obs==0] = NA
pred_unbou[obs==1] = NA

tfscore_bou_site <- apply(pred_bou,1,mean,na.rm=TRUE)
tfscore_unbou_site <- apply(pred_unbou,1,mean,na.rm=TRUE)
plot(tfscore_unbou_site,tfscore_bou_site); abline(0,1)
```

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-2-1.png){width="672"}

```r
## can't see much difference here between bound and unbound ind

tfscore_bou_ind <- apply(pred_bou,2,mean,na.rm=TRUE)
tfscore_unbou_ind <- apply(pred_unbou,2,mean,na.rm=TRUE)
rango <- range(c(tfscore_bou_ind,tfscore_unbou_ind))
plot(tfscore_unbou_ind,tfscore_bou_ind,xlim=rango, ylim=rango); abline(0,1); 
title("predictions get var across genome right")
```

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-2-2.png){width="672"}

```r
## bound tf scores higher than unbound tf scores
```

======= hist(pvec)

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-2-1.png){width="672"}

```r
qqunif(pvec)
```

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-2-2.png){width="672"}

::: {#permute-sites-for-observed .section .level2}
## permute sites for observed

```r
for(pp in 1:4)
{
obs_perm <- obs[sample(1:nrow(obs),nrow(obs),replace=FALSE),]
pvec_perm = rep(NA,nrow(obs))
for(cc in 1:nrow(obs) ) 
{
  res = t.test(pred[cc,obs_perm[cc,]==1], pred[cc,obs_perm[cc,]==0])
  pvec_perm[cc] = res$p.value
  }

titulo = glue("perm-{pp}:")
print(titulo)
##hist(pvec_perm,main=titulo)
qqunif(pvec_perm,main=titulo)
qqplot(-log10(pvec_perm), -log10(pvec), 
       main=paste(titulo ))
abline(0,1)
par(mfrow=c(1,1))

}
```

    ## perm-1:

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-1.png){width="672"}![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-2.png){width="672"}

    ## perm-2:

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-3.png){width="672"}![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-4.png){width="672"}

    ## perm-3:

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-5.png){width="672"}![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-6.png){width="672"}

    ## perm-4:

![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-7.png){width="672"}![](%7B%7B%3C%20blogdown/postref%20%3E%7D%7Dindex_files/figure-html/unnamed-chunk-3-8.png){width="672"}
:::

