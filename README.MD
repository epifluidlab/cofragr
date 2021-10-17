# cofragr: cfDNA co-fragmentation analysis in R

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#understanding-cofragr">Understanding cofragr</a></li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#preparation">Preparation</a></li>
        <li><a href="#analysis-with-rscript">Analysis with RScript</a></li>
        <li><a href="#analysis-with-snakemake">Analysis with Snakemake</a></li>
      </ul>
    </li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

## About The Project

cofragr is an R package for cofragmentation patterns of cell-free DNA.

## Getting Started

### Prerequisites

The cofragr package is still under development and hasn’t been released to CRAN, so you must install cofragr by using the **devtools** workflow. You can install **devtools** two different ways:

You can install the library directly from CRAN:

```R
install.packages("devtools")
```

You can also install the development version from GitHub:

```R
devtools::install_github("r-lib/devtools")
```

You can learn more about **devtools** and get help with its installation [here](https://devtools.r-lib.org/).


### Installation

Now that you've downloaded **devtools**, you can use the following command to install the package:

```R
devtools::install_github("epifluidlab/cofragr")
```

Congratulations! You've successfully downloaded and installed **cofragr**!

*Note: If you're having any problems installing cofragr, make sure that you have the latest version of R installed.*

## Understanding cofragr

Understanding how **cofragr** works can be crucial to its successful usage. There are three core functions to this R package:

##### read_fragments

Calling this function will read the fragment data file, which is typically a BED file. Since this functions uses the package **[bedtools](https://haizi-zh.github.io/bedtorch/)** to read the fragment files, the arguments/options for `range` and `genome` will follow the same format. 

###### Options

```
cofragr::read_fragments(file_path, range = NULL, genome = NULL)
```

`range; default: range = NULL`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specifies the ranges of the BED file that read_fragments will read.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note that `range` follows the standard genomic range notation, ex. "chr1:1001-2000". 

`genome; default: genome = NULL`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Specifies the reference genome for the BED file in question.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note that `genome` can take any of the `genome` arguments that [`GenomeInfoDb::Seqinfo()`](https://rdrr.io/bioc/GenomeInfoDb/man/Seqinfo-class.html) can.

Again, you can read more about the options behind `range` and `genome` [here](https://haizi-zh.github.io/bedtorch/reference/read_bed.html).

##### contact_matrix

Calling this function will calculate the contact matrix from the fragment data file, which should have already been called by the `read_fragments` function. 

###### Options

```
cofragr::contact_matrix(frag, bin_size = 500e3L, n_workers = 1L, subsample = 10e3L, min_sample_size = 100L, bootstrap = 1L, seed = NULL)
```
Note that `frag` is already defined as the fragment file as a result of the `read_fragments` function.

`bin_size; default: bin_size = 500e3L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the size of the bins for the contact matrix in base-pair units.

`n_workers; default: n_workers = 1L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the number of workers for the task at hand.

`subsample; default: subsample = 10e3L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the sampling size for the contact matrix calculations. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Note:* In this case, the subsample means for each genomic bin, we only randomly select 1 in every 10,000 fragments for the analysis. This counteracts the fact that the number of fragments are different between different bins.

`bootstrap; default: bootstrap = 1L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the number of bootstrap samples in the contact matrix calculations.

##### write_contact_matrix

Calling this function will write the contact matrix to a specified filepath; the contact matrix should have already been calculated by the `contact_matrix` function. 

###### Options

```
cofragr::write_contact_matrix(cm, file_path, comments = NULL)
```

`comments; default: comments = NULL`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Allows the addition of a character vector to be appended at the top of the written BED file as a header.

## Usage

cofragr provides both R APIs and an R script for the analysis. 

### Preparation

The pipeline requires fragment data files as input. A fragment data file is essentially a BED file, with each row representing a cfDNA fragment. Below is an example:

    14      19000035        19000198        .       27      +
    14      19000044        19000202        .       42      +
    14      19000045        19000202        .       20      -
    14      19000049        19000202        .       12      +

The six columns are: chrom, start, end, feature name (which we do not use), MAPQ (which is the smallest MAPQ of the two paired reads), and strand (which we do not use).

You can easily prepare the fragment data file from a query-sorted BAM file, or refer to our FinaleDB paper for more details (Zheng, Zhu, and Liu 2020).

### Analysis with RScript

The quickest way to start is using the **RScript** for your dataset.

You can find the script in R package installation:

    system.file("extdata/scripts", "cofragr.R", package = "cofragr")       

Or simply clone the source code from GitHub and find the script in the directory `inst/ext/scripts/`.


Then you can simply run some variation of the following, inputting any arguments as you desire.

    Rscript cofragr.R 
    -i examplefile.bed.gz \
    -o output_file_path \
    --min-mapq 30 \
    --chroms 1:2:3 \
    --subsample 10000 \
    --min-fraglen 50 \
    --max-fraglen 350 \
    --sample-id example \

We encourage you to look through the script file, `cofragr.R`, to understand all of the available arguments.

##### Output

You will be outputted two files, the calculated contact matrix and its accompanying index file. 

cofrag_cm.bed.gz is the calculated contact matrix:

    #chrom  start   end     chrom2  start2  end2    score   n_frag1 n_frag2 p_value p_value_sd
    1       0       500000  1       0       500000  15.14129840182        435     435     0.58172039812       .
    1       0       500000  1       500000  1000000 0       435     2539    0       .
    1       0       500000  1       1000000 1500000 0       435     3997    0       .
    1       0       500000  1       1500000 2000000 0       435     3707    0       .
    1       0       500000  1       2000000 2500000 0       435     4361    0       .


### Analysis with Snakemake

For large-scale real applications, we suggest use workflow management tools such as **snakemake**, nextflow, etc. We provide a snakemake file. Similar to cofragr.R, you can find it:

In R package installation:

    system.file("extdata/scripts", "cofragr.smk", package = "cofragr")

Or simply clone the source code from GitHub and find the script in the directory `inst/ext/scripts/`.

*Note: Usually snakemake workflow is more related to your specific computing environment. As a result, the example above is only for your reference. Have a look in both the snakemake file and the command above, and modify it accordingly.*

## License

See `LICENSE` for more information.

## Contact

* Haizi Zheng haizi.zh@gmail.com
* Yaping Liu lyping1986@gmail.com
* Ravi Bandaru ravi14.bandaru@gmail.com
