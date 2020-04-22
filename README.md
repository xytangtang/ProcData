
<!-- README.md is generated from README.Rmd. Please edit that file -->
ProcData: An R Package for Process Data Analysis
================================================

<!-- badges: start -->
<!-- badges: end -->
`ProcData` provides tools for exploratory process data analysis. It contains an example dataset and functions for

-   reading responses from a csv file
-   process manipulation
-   action sequence generators
-   feature extraction methods
-   fitting and making prediction from sequence models

Installation
------------

Download the package from [the download page](http://www.scientifichpc.com/processdata/download.html) and execute the following command in `R`

``` r
install.packages(FILENAME, repos = NULL, dependencies = TRUE)
```

where `FILENAME` should be replaced by the name of the package file downloaded including its path. The development version can be installed from [GitHub](https://github.com/) with:

``` r
devtools::install_github("xytangtang/ProcData")
```

ProcData depends on packages `Rcpp` and [`keras`](https://keras.rstudio.com). A C compiler and python are needed. Some functions in `ProcData` calls functions in `keras` to fit neural networks. To make sure these functions run properly, execute the following command in `R`.

``` r
library(keras)
install_keras(tensorflow = "1.13.1")
```

Note that if this step is skipped, `ProcData` can still be installed and loaded, but calling the functions that depends on `keras` will give an error.

Contents
--------

### Data Structure

`ProcData` organizes response processes as an object of class `proc` which is a list containing the action sequences and the timestamp sequences. Functions are provided to summarize and manipulate `proc` objects.

### Dataset

`ProcData` includes a dataset `cc_data` of the action sequences and binary item responses of 16920 respondents of item CP025Q01 in PISA 2012. The item interface can be found [here](http://www.oecd.org/pisa/test-2012/testquestions/question3/). To load the dataset, run

``` r
data(cc_data)
```

`cc_data` is a list of two elements:

-   `seqs` is a \`proc' object.
-   `responses` is a numeric vector containing the binary responses outcomes.

For data stored in csv files, `read.seqs` can be used to read response processes into R and to organize them into a `proc` object. In the input csv file, each process can be stored in a single line or multiple lines. The sample files for the two styles are example\_single.csv and example\_multiple.csv. The processes in the two files can be read by running

``` r
seqs1 <- read.seqs(file="example_single.csv", style="single", id_var="ID", action_var="Action", time_var="Time", seq_sep=", ")
seqs2 <- read.seqs(file="example_multiple.csv", style="multiple", id_var="ID", action_var="Action", time_var="Time")
```

`write.seqs` can be used to write `proc` objects in csv files.

### Data Generators

`ProcData` also provides three action sequences generators:

-   `seq_gen` generates action sequences of an imaginary simulation-experiment-based item;
-   `seq_gen2` generates action sequences according to a given probability transition matrix;
-   `seq_gen3` generates action sequences from a recurrent neural network. It depends on `keras`.

### Feature Extraction Methods

`ProcData` implements two feature extraction methods that compress varying length response processes into fixed dimension numeric vectors. One of the methods is based on multidimensional scaling (MDS) and the other one is based on sequence-to-sequence autoencoders (seq2seq AE). Details of the two methods can be found [here](http://www.scientifichpc.com/processdata/method.html).

#### MDS

The following functions implement the MDS methods.

-   `seq2feature_mds` extracts `K` features from a given set of response processes or their dissimilarity matrix.
-   `chooseK_mds` selects the number of features to be extracted by cross-validation.

``` r
seqs <- seq_gen(100)
K_res <- chooseK_mds(seqs, K_cand=5:10, return_dist=TRUE)
theta <- seq2feature_mds(K_res$dist_mat, K_res$K)$theta
```

#### seq2seq AE

Similar to MDS, the seq2seq AE method is implemented by two functions. Both functions depend on `keras`.

-   `seq2feature_seq2seq` extracts `K` features from a given set of response processes.
-   `chooseK_seq2seq` selects the number of features to be extracted by cross-validation.

``` r
seqs <- seq_gen(100)
K_res <- chooseK_seq2seq(seqs, K_cand=c(5, 10), valid_prop=0.2)
seq2seq_res <- seq2feature_seq2seq(seqs, K_res$K, samples_train=1:80, samples_valid=81:100)
theta <- seq2seq_res$theta
```

Note that if the number of candidates of `K` is large and a large number of epochs is needed for training the seq2seq AE, `chooseK_seq2seq` can be slow. One can parallel the selection procedure via multiple independent calls of `seq2feature_seq2seq` with properly specified training, validation, and test sets.

### Sequence Models

A sequence model relates response processes and covariates with a response variable. The model combines a recurrent neural network and a fully connected neural network.

-   `seqm` fits a sequence model. It returns an object of class \`seqm'.
-   `predict.seqm` predicts the response variable with a given fitted sequence model. Both `seqm` and `predict.seqm` depends on `keras`.

``` r
n <- 100
seqs <- seq_gen(n)
y1 <- sapply(seqs$action_seqs, function(x) "CHECK_A" %in% x)
y2 <- sapply(seqs$action_seqs, function(x) log10(length(x)))

index_test <- sample(1:n, 10)
index_train <- setdiff(1:n, index_test)
seqs_train <- sub_seqs(seqs, index_train)
seqs_test <- sub_seqs(seqs, index_test)

actions <- unique(unlist(seqs))

# a simple sequence model for a binary response variable
seqm_res1 <- seqm(seqs = seqs_train, response = y1, response_type = "binary",
             actions=actions, K_emb = 5, K_rnn = 5, n_epoch = 5)
pred_res1 <- predict(seqm_res1, new_seqs = seqs_test)

# a simple sequence model for a numeric response variable
seqm_res2 <- seqm(seqs = seqs_test, response = y2, response_type = "scale",
             actions=actions, K_emb = 5, K_rnn = 5, n_epoch = 5)
pred_res2 <- predict(seqm_res2, new_seqs = seqs_test)
```
