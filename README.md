# cytoSNOW

*cytoSNOW* speeds up computational analyses of flow and mass cytometry data in R.
Using the [`snow`](https://cran.r-project.org/web/packages/snow/index.html) framework, it runs essential components of these workflows using multiple CPU cores.
This makes it possible to run things locally and **fast**.

Below are the tasks that *cytoSNOW* accelerates.

| Task | Original implementation | Fast implementation |
| :-- | :-- | :-- |
| Preprocessing                | `flowCore` compensation and transformation | `ParallelPreprocess`  |
| FCS file aggregation         | `FlowSOM::AggregateFlowFrames`             | `ParallelAggregate`   |
| `FlowSOM` feature extraction | `FlowSOM::GetFeatures`                     | `ParallelGetFeatures` |

## Installation

See dependencies in `DESCRIPTION`.
Install using `devtools`:

```R
devtools::install_github('davnovak/cytoSNOW')
```

All exported functions are documented.


