# cytoSNOW

This is a package that uses multi-threading for speeding up some operations in cytometry data analysis workflows.
It relies on [`snow`](https://cran.r-project.org/web/packages/snow/index.html) to set up a local cluster and use multiple CPU cores at a time.

## Progress

*cytoSNOW* is work in progress.
I have used parallelisation for a bunch of analytical workflows, but have not yet added all of the implementations to the *cytoSNOW* package.
Below is a task list with implementation status.

| Task | Implementation |
| :-- | :-- |
| FCS file aggregation         | `ParallelAggregate` |
| `FlowSOM` feature extraction |  |
| `diffcyt` testing            |  |
| `CytoNorm` training          |  |
| `CytoNorm` normalisation     |  |

## Installation

See dependencies in `DESCRIPTION`.

```R
devtools::install_github('davnovak/cytoSNOW')
```

All exported functions are documented.