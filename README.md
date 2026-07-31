# kinetics-from-total-rna

This repository contains a pipeline that uses total RNA-seq to estimate changes in the kinetic
parameters of transcription and RNA processing: the transcription initiation rate, the elongation
speed, the intron splicing speed, and the RNA degradation rate. Because total RNA-seq observes RNA
at steady state, these parameters are identifiable only relative to one another, so every effect
size the pipeline reports is a fold change of a *ratio* of two of them; see
[What the model estimates](#what-the-model-estimates) below.

The code is organized as a Nextflow pipeline. To run it, you will need FASTQ
data of total RNA-seq; please read the documentation below for more details.

The model itself is implemented using PyTorch; the code can be found in
the [rna_kinetics](rna_kinetics) folder
(which can also be installed as a pip package).

> **Note on the preprint.** A preprint is
> [available on bioRxiv](https://doi.org/10.1101/2025.08.24.672013), but it describes an earlier
> version of this method. The code in this repository has since changed substantially — including
> the quantities that are reported and how they are named — and no longer corresponds to what that
> preprint describes. An updated manuscript is in preparation; until then, please treat this README
> as the authoritative description of what the code does.
>
> The current draft is included in this repository as [preprint/preprint.tex](preprint/preprint.tex)
> (with a rendered [PDF](preprint/preprint.pdf)). It is **unfinished** — the results and discussion
> are placeholders — but *Section 2, Dynamical Model of RNA Metabolism*, is essentially complete.
> That section derives the underlying physical model of RNA metabolism, its steady-state solution,
> and the structural non-identifiability described below, and is the best starting point if you
> want to understand what the pipeline is actually fitting.

## What the model estimates

**Every reported effect size is a log fold change (LFC) of a *ratio* of two rates, never of a single
rate.** This is a structural property of the data, not a shortcoming of the fit: total RNA-seq
measures RNA at steady state, and steady-state abundances constrain the underlying kinetic rates
only relative to one another. Scaling all rates by a common factor leaves the predictions
unchanged, so no individual rate is identifiable.

Concretely, the pipeline reports the following parameters (using the same names that appear in the
`parameter_type` column of the output):

| Parameter | Interpretation |
| --- | --- |
| `lfc_init_over_deg` | LFC(initiation rate) − LFC(degradation rate) |
| `lfc_elong_over_deg` | LFC(elongation speed) − LFC(degradation rate) |
| `lfc_splice_over_deg` | LFC(splicing speed) − LFC(degradation rate) |
| `lfc_elong_over_splice` | LFC(elongation speed) − LFC(splicing speed) |

The practical consequence is that a significant result identifies which *pair* of rates has moved
relative to the other, and never by itself which member of the pair is responsible. A positive
`lfc_elong_over_deg`, for instance, is equally consistent with faster elongation and with slower
degradation. Interpreting the parameters jointly helps: if `lfc_elong_over_deg` and
`lfc_splice_over_deg` shift together by a similar amount while `lfc_elong_over_splice` stays near
zero, a change in degradation is the more parsimonious explanation than a coincident change in both
elongation and splicing.

Of the four, only `lfc_init_over_deg` maps onto a familiar quantity: it is the log fold change in
steady-state abundance, i.e. what a conventional differential-expression analysis reports. Further
ratios can be formed by subtracting the parameters from each other — for example
`lfc_init_over_deg` − `lfc_elong_over_deg` = LFC(initiation / elongation speed), a change in
polymerase density along the gene body — but only the four parameters in the table are covered by
the likelihood-ratio tests.

Note also that `lfc_init_over_deg` is always estimated **per gene**, in every model. A global
version of it would not be identifiable: RNA-seq is compositional data, capturing only relative
abundances, so an expression shift shared by *all* genes is indistinguishable from a change in
sequencing depth and is absorbed by the library size factors. The other two parameters do not have
this problem, because they are identified by the intron-to-exon read ratio *within* each gene,
which is unaffected by library size. Accordingly, the global RNA kinetics model tests only
`lfc_elong_over_deg` and `lfc_splice_over_deg`.

Effect sizes are always reported as log2 fold changes. Significance comes from likelihood-ratio
tests against a reduced design (see *LRT contrasts* below).

## Choosing a model

The pipeline provides six models. They differ along two axes, and you can enable any combination of
them in `dataset_params.yaml`; each writes to its own output subfolder.

**Axis 1 — which data are used.**

* **RNA kinetics models** fit exon read counts, intron read counts, and intron coverage jointly.
  Using the read counts as well lets them resolve all three degradation-referenced parameters
  (`lfc_init_over_deg`, `lfc_elong_over_deg`, `lfc_splice_over_deg`).
* **Intron coverage models** fit the intron coverage alone, and report a single parameter,
  `lfc_elong_over_splice`. This is the simpler and more robust model, at the cost of not separating
  elongation and splicing from degradation.

**Axis 2 — granularity**, i.e. how many introns share one estimated effect. Coarser granularity
pools more data and therefore has more statistical power, but assumes the effect is genuinely
shared.

| Granularity | RNA kinetics | Intron coverage |
| --- | --- | --- |
| **Global** — one effect for the whole dataset | `fit_global_rna_kinetics_model` | `fit_global_intron_coverage_model` |
| **Per gene** — one effect per gene, shared by its introns | `fit_rna_kinetics_model` | `fit_gene_specific_intron_coverage_model` |
| **Per intron** — a separate effect for every intron | `fit_intron_specific_rna_kinetics_model` | `fit_intron_specific_intron_coverage_model` |

A good starting point is to run the global and per-gene models of both families, which is what
`dataset_params_template.yaml` enables by default. The global fits tell you whether there is a
dataset-wide shift and stay well powered even in small experiments, while the per-gene fits
identify which genes drive it.

The two per-intron models are off by default. They are the most granular but the least powered,
since each estimate draws on a single intron; the per-intron RNA kinetics model is additionally the
most expensive to compute. They are worth enabling when you specifically expect introns within a
gene to behave differently — otherwise the per-gene models answer the same question with more data
behind each estimate. As with every flag here, you can enable any combination you like.

Granularity also changes what the output contains. Per-gene and per-intron models run an additional
regularization step, shrinking effect sizes towards zero with an empirical Bayes prior fitted by
[ashr](https://github.com/stephens999/ashr) across all genes; their `test_results.tsv` therefore
carries both `l2fc_unregularized` and `l2fc_regularized`. A global model estimates a single effect,
so there is nothing to pool a prior across and no shrinkage is applied — its output has one
estimate, `l2fc_mle`.

## Running the pipeline

Besides the FASTQ data, you will need to prepare two files to run the pipeline:

* ```samplesheet.csv```, containing sample annotation.
* ```dataset_params.yaml```, which specifies dataset metadata and several paths on your filesystem.

### Preparing samplesheet

The file ```samplesheet.csv``` needs to contain the following 4 columns:

* *sample*, specifying sample names (the names can be arbitrary).
* *fastq_1*, name of FASTQ file with reads 1 (compressed by gunzip).
* *fastq_2*, same as above for reads 2.
* *strandedness*, specifying the strandedness (read orientation) of the samples. Currently, only values ```forward```
  and ```reverse``` are supported.
  (The most common Illumina RNA-seq protocols use ```reverse``` orientations).

Besides these 4 mandatory columns, the samplesheet can contain arbitrary explanatory variables (e.g., genotype,
intervention, etc.), that can
be used in the dataset parameter file to specify a design formula (see below).

### Preparing dataset parameters

The file ```dataset_params.yaml``` can be created by copying
the [dataset_params_template.yaml](dataset_params_template.yaml) file
and filling it according to the comments. We will now discuss several parameter details:

* **Reference genome and transcriptome** (*gtf_file*, *genome_fasta*, *transcriptome_fasta*, and *gtf_source*): In
  the pipeline, we are aligning reads both to the genome and transcriptome. Therefore, the reference
  files (*gtf_file*, *genome_fasta*, and *transcriptome_fasta*) need to be all compatible (from the same source and
  version).

  We are using [STAR](https://github.com/alexdobin/STAR) for aligning reads to the genome; please follow the
  recommendations
  from the STAR documentation for details on how to choose the reference genome files
  (see section *2: Generating genome indexes* of
  the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)). Typically, the *primary assembly* version of the annotation is the correct one to use.

  Since Ensembl and Gencode annotations slightly differ, the used source needs to be stated in the *gtf_source*
  parameter
  (supported values are ```"ensembl"``` and ```"gencode"```).

  We are using [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) package to import transcript abundances from Salmon. Tximport requires
  argument *ignoreTxVersion*, which typically should be true for human or mouse reference annotation from Ensemble, and false 
  for most other references. 
* **STAR and Salmon indices**: You can optionally pre-build the STAR index and/or Salmon index yourself and provide the
  path to it. If no path
  is provided (the parameter is left ```null```), the index will be built from the provided genome/transcriptome files.
* **Design formula**: the parameter *design_formula* uses R syntax to generate a design matrix from the formula (
  see [R documentation](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula) for details).
  The explanatory variables used in the formula need to be the columns of the ```samplesheet.csv```.
* **LRT contrasts**: We are using likelihood-ratio tests to assess the significance of LFC parameters. The parameter
  *lrt_contrasts* specifies a list of tests that are going to be performed.

  For categorical variables, please specify (for each test) the *variable* (name of the column in the samplesheet), and 
  comparison groups *group_1* and *group_2* (levels of the variable). The null hypothesis being tested is whether the LFC
  parameters between *group_1* and *group_2* are equal to 0.

  For example, assume that the samplesheet contains column *genotype*, with levels *wild_type*, *knockout_1*, and *knockout_2*.
  To test whether each of the knockout group is different from the wild-type, you can specify:
  ```
  lrt_contrasts:
  - variable: 'genotype'
    group_1: 'knockout_1'
    group_2: 'wild_type'
    
  - variable: 'genotype'
    group_1: 'knockout_2'
    group_2: 'wild_type'
  ```
  You can include as many tests as you wish, e.g., compare also groups in some other explanatory variable, or include also
  the comparison between *knockout_1* and *knockout_2* in the *genotype* variable.

  For a continuous variable, specify only the *variable* parameter, without the comparison groups. 

  We currently support testing only the main-effect term labels from an R terms object, not interaction or other
  higher-order term labels. If you wish to test e.g., an interaction term between two variables, please manually create 
  the interaction variable as a separate column in the *samplesheet*; then, you can add it to the design formula 
  and also include it in the tests.

* **Stage**: Our pipeline consists of 2 workflows: [pre-processing](./workflows/preprocessing.nf)
  and [modeling](./workflows/modeling.nf).
  The pre-processing contains read alignment and related steps, and is supposed to be
  run only once per dataset. The modeling, on the other hand, may be run multiple times, e.g, experimenting with
  different
  design formulas.

  The parameter *stage* specifies whether both or only one workflow should be run. Possible values are
  ```"all"``` (to run both workflows), ```"preprocess"``` and ```"model"```. 
  
  Please note that you can also generally use the ```-resume``` argument for the ```nextflow run``` command. Setting ```stage: 'model'```
  is simply an orthogonal way to re-use preprocessed data, independent of the caching done via  Nextflow work folder.

### Executing the pipeline

Please install [Nextflow](https://www.nextflow.io/) (version >= 25.04) and [Docker](https://www.docker.com/) on your system.
Then, the pipeline can be run by

```commandline
nextflow run main.nf -params-file dataset_params.yaml
```

See the Nextflow documentation for [details](https://www.nextflow.io/docs/latest/executor.html) on how to run the
pipeline on your HPC/cloud system.

**Note on Slurm**: On our HPC system using Slurm, we noticed the following bug: when several of the
model-fitting processes — that is, any of the processes that run PyTorch — are executed on the same
cluster node, their CPU usages somehow collide, resulting in orders of magnitude slower process
execution. We suspect that this may be related to the underlying BLAS settings and/or our cluster
setup, and we are currently trying to resolve this issue. In the case that you experience similar
behavior, please execute at most one model-fitting process per cluster node.

The `beyer_cluster` profile in [nextflow.config](nextflow.config) contains the crude workaround we
currently use for this: each fitting process requests substantially more memory than it actually
needs (130 GB, somewhat over half of a node), so that Slurm cannot place two of them on the same
node. This is deliberately a hack — the requested memory goes unused, and it relies on a side
effect of the scheduler rather than expressing what we actually want.

We kept it because the alternatives we tried were worse. Setting explicit thread limits
(`OMP_NUM_THREADS`, `torch.set_num_threads` and similar) did not resolve the collision in our tests,
and most variants made performance worse rather than better. Requesting nodes exclusively
(`--exclusive`) does prevent the problem, but reserves the entire node, whereas over-requesting
memory still leaves the remaining cores available to other jobs.

## Contact

If you have any questions or experience any problems with the code,
please open an issue on this repository, or reach out at
[jkoubele@uni-koeln.de](mailto:jkoubele@uni-koeln.de).