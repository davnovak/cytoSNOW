<img src="./logo.png" alt="cytoSNOW" width="250"/>

*cytoSNOW* speeds up computational analyses of flow and mass cytometry data in R.
Using the [`snow`](https://cran.r-project.org/web/packages/snow/index.html) framework, it runs essential components of these workflows using multiple CPU cores.
This makes it possible to run things **fast**.
This is especially useful if you're running workflows **locally**, but *cytoSNOW* can power analyses **on high-performance computing (HPC) clusters** as well.

I wrote about an initial experiment with *cytoSNOW* in a blog post ([here](https://davnovak.github.io/docs/cytoSNOW/vignette.html)).

<hr>

Below are the tasks that *cytoSNOW* accelerates right now.
The scope is being extended continuously: for instance, differential abundance/state testing via *diffcyt* will be included.

| Task | Original implementation | Fast implementation |
| :-- | :-- | :-- |
| Preprocessing               | `flowCore` compensation and transformation | `cytoSNOW::ParallelPreprocess`  |
| FCS file aggregation         | `FlowSOM::AggregateFlowFrames`             | `cytoSNOW::ParallelAggregate`   |
| `FlowSOM` feature extraction | `FlowSOM::GetFeatures`                     | `cytoSNOW::ParallelGetFeatures` |
| `CytoNorm` batch effect correction | `CytoNorm::CytoNorm.train` and `CytoNorm::CytoNorm.normalize` | `cytoSNOW::ParallelNormalize.Train` and `cytoSNOW::ParallelNormalize.Apply` **(caution: still undergoing tests)** |

*cytoSNOW* will speed up your analysis as long as you can use more than 1 CPU core.
We'll include guidelines for use in cloud computing and with HPCs eventually.

## Installation

Install using `devtools`:

```R
devtools::install_github('davnovak/cytoSNOW')
```

## Usage

The vignette (`vignette.Rmd`) contains an easy walk-through that quantifies the speed-up brought about by using *cytoSNOW*.
It mirrors the standard *FlowSOM* protocol, as published [here](https://www.nature.com/articles/s41596-021-00550-0).

## Notes

Some general notes on parallelizing computational cytometry workflows.

### Overhead

Multi-core processing has some overhead.
This means that your computer takes a little bit of time to set up your tasks to run in parallel.
The more FCS samples you're processing, the more using *cytoSNOW* pays off.
However, there will be rare cases where *cytoSNOW* doesn't provide a boost in speed, if you're processing very few samples.

### Low-hanging fruit

Many parts of a typical computational cytometry workflow can be optimized.
But some, the 'low-hanging fruit', are worth the effort of optimizing more than others.
For instance, *FlowSOM* clustering often takes up very little time compared to preprocessing, data aggregation, feature extraction, or batch effect correction&horbar;at least in workflows with many FCS files.
That's why we focus preferentially on those.
