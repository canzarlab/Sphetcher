# Sphetcher
A software package for sampling massive singe-cell RNA secquencing datasets based on the spherical thresholding algorithm. It selects a small subset of cells referred to as ```sketch``` that evenly cover the transcriptomic space occupied by the original dataset. Such a sketch can accelerate downstream analyses and highlight rare cell types.

<img src=img/overviewv2.png  width="100%" height = "550">

## Installation ##
A compiler that supports C++11 is needed to build sphetcher. You can download and compile the latest code from github as follows:

```
git clone  https://github.com/canzarlab/Sphetcher
cd src
make
```

## Running Sphetcher ##

### To begin ###

First you need to provide a matrix where rows are samples (cells) and columns are features (genes, transcripts, principal components PCs).

Additionally you can provide the prior information (e.g, cell label, collection time point) in case you want to preserve certain number of samples from each category. 

An example of inputs is provided in the directory ```/data```. 

### Usage ###

Once you have compiled Sphetcher it can be run easily with one of the following two options:

```
sphetcher expression_matrix.csv sketch_size sketch_indicator_output.csv
```
or 
```
sphetcher expression_matrix.csv sketch_size class_labels.csv l_min sketch_indicator_output.csv
```
For an example provided in ```/data```
```
sphetcher zeisel_pca.csv 1000 sketch_indicator_output.csv
or 
sphetcher zeisel_pca.csv 1000 zeisel_pca_labels.csv 3 sketch_indicator_output.csv
```

#### Input/Output formats

Input: 

`expression_matrix.csv`
  : expression matrix in comma-separated values (CSV) format: rows are cells, columns are features. <br/>
 `sketch_size` 
  : number of samples to obtain from the data set. <br/>
`class_labels.csv`
  : prior information stored in a column vector, each class is presented by an integer between 1 and K, where K is the number of classes. <br/>
`l_min`
  : minimum number of representatives we want to sample from each class <br/>

Output:

`sketch_indicator_output.csv` : an indicator vector of `n` samples where 1 indicates the sample is in the sketch, 0 otherwise (there are `sketch_size` 1's in the vector).

## Sketches of large single cell datasets ##
#### Adult mouse brain cells from [Saunders et al. (2018)](https://www.cell.com/cell/fulltext/S0092-8674(18)30955-3)  (665,858 cells) ####
[`6646 cells (1%)`](https://syncandshare.lrz.de/dl/fi313Db6B6BeRPnggCuy7vzJ/saunders6646.tar.xz "Saunders")
[`26573 cells (4%)`](https://syncandshare.lrz.de/dl/fiEF5PeaX2yUtnjPwyLapFMJ/saunders26573.tar.xz "Saunders")
[`66649 cells (10%)`](https://syncandshare.lrz.de/dl/fi8MWmP4uVMHK9DejCKaL4yB/saunders66649.tar.xz "Saunders")
#### Mouse organogenesis cell atlas (MOCA) from [Cao et al. (2019)](https://www.nature.com/articles/s41586-019-0969-x) (2,026,641 cells) ####
[20296 cells (1%)](https://syncandshare.lrz.de/dl/fiS6FYf4Fx4HKq2ksuGg7uh6/trapnell20296.tar.xz "Trapnell"),
[81093 cells (4%)](https://syncandshare.lrz.de/dl/fiURXP6azi99M9fn2b9gdVpH/trapnell81093.tar.xz "Trapnell")  
#### Mouse nervous system from [Zeisel et al. (2018)](https://www.cell.com/cell/fulltext/S0092-8674(18)30789-X) (465,281 cells) #### 
[4643 cells (1%)](https://syncandshare.lrz.de/dl/fiE1S8cA9sMf7gMSAJ7hEyPe/zeiselCNS4643.tar.xz "zeiselCNS"),
[18577 cells (4%)](https://syncandshare.lrz.de/dl/fiFc6aCNPp6P8tDeh8WaGoQa/zeiselCNS18577.tar.xz "zeiselCNS"),
[46450 cells (10%)](https://syncandshare.lrz.de/dl/fiHeJebsTd2CfLuU555wvtxd/zeiselCNS46450.tar.xz "zeiselCNS")  
