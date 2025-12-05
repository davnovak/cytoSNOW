<img src="./logo.png" alt="cytoSNOW" width="250"/>

*cytoSNOW* speeds up computational analyses of flow and mass cytometry data in R.
Using the [`snow`](https://cran.r-project.org/web/packages/snow/index.html) framework, it runs essential components of these workflows using multiple CPU cores.
This makes it possible to run things locally and **fast**.
I have written about an initial experiment with *cytoSNOW* in a blog post ([here](https://davnovak.github.io/docs/cytoSNOW/vignette.html)).

Below are the tasks that *cytoSNOW* accelerates right now.
The scope will be extended continuously.

| Task | Original implementation | Fast implementation |
| :-- | :-- | :-- |
| Pre-processing               | `flowCore` compensation and transformation | `cytoSNOW::ParallelPreprocess`  |
| FCS file aggregation         | `FlowSOM::AggregateFlowFrames`             | `cytoSNOW::ParallelAggregate`   |
| `FlowSOM` feature extraction | `FlowSOM::GetFeatures`                     | `cytoSNOW::ParallelGetFeatures` |

*cytoSNOW* will speed up your analysis as long as you can use more than 1 CPU core on your machine.
This is especially useful if you want to run the analysis locally, on your desktop or laptop.
If you're running your analysis in a cloud, make sure you are using multiple CPU cores (otherwise *cytoSNOW* won't help).
We will also add implementations and guidelines for accelerating workflows run on HPCs.

## Installation

Install using `devtools`:

```R
devtools::install_github('davnovak/cytoSNOW')
```

## Usage

The vignette (`vignette.Rmd`) contains a nice walk-through that quantifies the speed-up brought about by using *cytoSNOW*.

Alternatively, the code snippets below can be used to generate a synthetic large dataset and deploy a generic analytical workflow using *cytoSNOW*.
The workflow is based on a standard *FlowSOM* protocol, as published [here](https://www.nature.com/articles/s41596-021-00550-0).

Each *cytoSNOW* function used below is carefully documented.
You can view the documentation using ```?cytoSNOW:::NameOfFunction```.

<details>
<summary><b>Input data generation</b></summary>
<br>

We begin by generating a synthetic dataset of `N` samples (2000 by default).
(In reality, we only create a single FCS file but reuse it, pretending there are 2000 of them.)

```r
## Simulate synthetic data ----

N  <- 2000          # sample count
nr <- 3e5           # cells per sample
nc <- 30            # number of markers
idcs_type <- 1:20   # markers for cell type
idcs_state <- 21:30 # markers for cell state

markers <- paste0('Marker', seq_len(nc))
set.seed(1); ff <- cytoSNOW::ValidateFCS(
  `colnames<-`( # Gaussian noise
    matrix(rnorm(nr*nc, mean = 10, sd = 5), ncol = nc),
    markers
  )
)

fname_input <- 'InputSample.fcs'
flowCore::write.FCS(ff, fname_input)
fnames <- rep(fname_input, times = N)

Sys.setenv( # exception for cytoSNOW to use same file multiple times
  'DUPLICATE_EXCEPTION' = TRUE 
)
```

To be able to simulate pre-processing, we also generate a spillover matrix for compensation and a `flowCore::transformList` for signal transformation.

```r
## Create a spillover matrix for compensation ----

set.seed(1); spillover <- # Gaussian noise with 1 on the diagonal
  `diag<-`(matrix(abs(rnorm(nc**2, mean = 1e-2, sd = 1e-3)), ncol = nc), 1.)
rownames(spillover) <- colnames(spillover) <- markers

## Create transformation instructions per channel ----

tf_list <- flowCore::transformList(
  from = markers, tfun = flowCore::arcsinhTransform(b = 120)
)
```
<hr>
</details>
<details>
<summary><b>Aggregation</b></summary>
<br>

The next step in the standard protocol is to aggregate expression data from all the pre-processed files, to obtain training data for the clustering model.

```r
agg <- cytoSNOW::ParallelAggregate(fnames = fnames_pre, N = 1e6)
```

This creates a 1-million-cell expression matrix that samples cells from all the files, making sure each sample is represented.
<hr>
</details>
<details>
<summary><b>Clustering</b></summary>
<br>

*FlowSOM* clustering itself is actually very fast.
Despite numerous approaches optimising this process to run faster, this is rarely the real bottleneck.

```r
fsom <- FlowSOM::FlowSOM(agg1, nClus = 40, colsToUse = markers[idcs_type])
```
<hr>
</details>
<details>
<summary><b>Feature extraction</b></summary>
<br>

To be able to compare cell types and cell state across cytometry samples, we would typically use feature extraction as implemented in *FlowSOM*.
Here, we accelerate the feature extraction process.

```r
fe <- cytoSNOW::ParallelGetFeatures(
  fsom          = fsom,
  fnames        = fnames_pre,
  level         = c('clusters', 'metaclusters'),
  type          = c('counts', 'proportions', 'medians'),
  state_markers = markers[idcs_state]
)
```
<hr>
</details>