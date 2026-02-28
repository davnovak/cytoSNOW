<img src="./logo.png" alt="cytoSNOW" width="250"/>

*cytoSNOW* is an R package for **fast** analysis of flow and mass cytometry data.

Using the [`snow`](https://cran.r-project.org/web/packages/snow/index.html) framework, it runs essential workflow components using multiple CPU cores.
It can be used both **locally** or **on high-performance computing (HPC) clusters**.

<i>Check out a proof-of-concept <b>[here](https://davnovak.github.io/docs/cytoSNOW/vignette.html)</b>.</i>

<hr>

See *cytoSNOW*-accelerated tasks below.
We're extending the scope continually.
For example, differential abundance/state testing via [*diffcyt*](https://www.nature.com/articles/s42003-019-0415-5) will be included soon.

| Task | Original implementation | Fast implementation |
| :-- | :-- | :-- |
| Preprocessing               | `flowCore` compensation &amp; transformation | `cytoSNOW::ParallelPreprocess`  |
| Aggregating FCS files         | `FlowSOM::AggregateFlowFrames`             | `cytoSNOW::ParallelAggregate`   |
| `FlowSOM` feature extraction | `FlowSOM::GetFeatures`                     | `cytoSNOW::ParallelGetFeatures` |
| `CytoNorm` batch effect correction | `CytoNorm::CytoNorm.train` and `CytoNorm::CytoNorm.normalize` | `cytoSNOW::ParallelNormalize.Train` and `cytoSNOW::ParallelNormalize.Apply` |

We'll include guidelines for use in cloud computing and with HPCs eventually.

## Installation

Install using `devtools` from your R console:

```R
devtools::install_github('davnovak/cytoSNOW')
```

## Usage

`vignette.Rmd` gives an easy walkthrough.
It quantifies speed-up by *cytoSNOW* in applying the standard [*FlowSOM* protocol](https://www.nature.com/articles/s41596-021-00550-0).

## Is *cytoSNOW* for you?

See the **example situations** below.
Does any of them apply to your work?
Then you will like *cytoSNOW*.

<details>
<summary><i><b>"We need to redo the preprocessing"</b></i></summary>
<br>

Cytometry data preprocessing takes many iterations, as parameters need to be optimized.
Despite automation efforts, **large analyses require 10s or 100s** of re-runs and coordination between wet and dry lab.
*cytoSNOW* dramatically reduces the time spent on this.

</details>

<details>
<summary><i><b>"We need to explore the data"</b></i></summary>
<br>

Many computational cytometry pipelines focus on **exploratory data analysis (EDA)**, as opposed to strictly hypothesis-driven testing.
Cytometry EDA can involve **massive amounts of statistical tests**, to detect compositional and phenotypic differences between samples and groups.
*cytoSNOW* will soon accelerate that as well, running [*diffcyt*](https://www.nature.com/articles/s42003-019-0415-5) under the hood.
We implemented this initially in [*iidx*](https://github.com/saeyslab/iidx).

</details>

<details>
<summary><i><b>"There's too much data!"</b></i></summary>
<br>

Running analyses overnight because you're **looping over hundreds or thousands of FCS files**?
Chances are, *cytoSNOW* will make you much more efficient.

</details>

<details>
<summary><i><b>"Please don't make me use the cluster"</b></i></summary>
<br>

Maybe you don't have **access to a university/company HPC cluster**.
Or it's busy and you can't wait to queue.
Or it's a pain to install all your packages on it.
With *cytoSNOW*, you might be able to just run even big analyses on your work or home computer.

</details>

## Original features

*cytoSNOW* speeds up analyses that use *flowCore*, *FlowSOM*, *CytoNorm*, and other R packages under the hood.
In addition to this, it offers **extra functions that the base packages do not directly implement**:

+ `GetPanels` efficiently extracts names of channels and markers from FCS files without opening them
+ `ValidateFCS` sanitizes metadata when creating new/modified FCS files in R
+ `ParallelAggregate` can use a gating or *FlowSOM* model to only extract a filtered subset of events per file
+ `ParallelNormalize.Train` and `ParallelNormalize.Apply` can use a gating model instead of *FlowSOM* to split samples into compartments for batch effect correction

## Notes on parallelizing

Some general notes on parallelizing computational cytometry workflows.

<details>
<summary><b>Low-hanging fruit</b></summary>
<br>

Many parts of a typical workflow can be optimized.
But some, **the "low-hanging fruit," are worth the effort of optimizing more than others**.
For instance, *FlowSOM* clustering often takes up very little time compared to preprocessing, data aggregation, feature extraction, or batch effect correction&horbar;at least in workflows with many FCS files.
That's why we focus preferentially on those.

</details>

<details>
<summary><b>Computational overhead</b></summary>
<br>

Multi-core processing has some overhead.
This means that your computer takes a little bit of time to set up your tasks for running in parallel.
The more FCS samples you're processing, the more *cytoSNOW* pays off.
However, in rare cases with very few samples *cytoSNOW* will not accelerate your workflow.

</details>
