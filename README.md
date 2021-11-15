## To install the package

## Build from the source
```
devtools::build()
devtools::install()
```
## Pre-requisite
The installation depends on the (RCTD)[https://github.com/dmcable/RCTD] package 
for initial deconvolution and annotation transfer from an existing scRNA-seq
dataset. Install RCTD by following the instructions given in the linked repo.   

## Install from the github repo
```
devtools::github_install('kharchenkolab/slideseqde/')
```

## Introduction
The `slideseqde` package aims to find the differentially expressed genes 
for a target cell-type in a "context" sensitive manner. While the context 
depends on the biological question one sets out to ask and can vary, we 
consider coherant spatial segments as contexts. Below is a graphical
outline of a typical slide-seq DE pipeline we propose. A more detailed
explanation can be obtained from the attached (vignette)[https://github.com/kharchenkolab/slideseqde/blob/master/vignettes/slideseqde.Rmd]

![](https://i.imgur.com/dNolHws.png)
