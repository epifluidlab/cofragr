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
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>

## About The Project

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

## Usage

There are three core functions to this R package:

##### read_fragments

Calling this function will read the fragment data file, which is typically a BED file. Since this functions uses the package **[bedtools](https://haizi-zh.github.io/bedtorch/)** to read the fragment files, the arguments/options for `range` and `genome` will follow the same format. 

To call the function in R, use:

```
cofragr::read_fragments(file_path, [options])
```

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

To call the function in R, use:

```
cofragr::contact_matrix(frag, [options])
```

Note that `frag` is already defined as the fragment file as a result of the `read_fragments` function.

###### Options

```
cofragr::contact_matrix(frag, bin_size = 500e3L, n_workers = 1L, subsample = 10e3L, min_sample_size = 100L, bootstrap = 1L, seed = NULL)
```

`bin_size; default: bin_size = 500e3L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the size of the bins for the contact matrix in base-pair units.

`n_workers; default: n_workers = 1L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the number of workers for the task at hand.

`subsample; default: subsample = 10e3L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the sampling size for the contact matrix calculations.

`bootstrap; default: bootstrap = 1L`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Specifies the number of bootstrap samples in the contact matrix calculations.

##### write_contact_matrix

Calling this function will write the contact matrix to a specified filepath; the contact matrix should have already been calculated by the `contact_matrix` function. 

To call the function in R, use:

```
cofragr::write_contact_matrix(cm, file_path, [OPTIONS])
```

###### Options

```
cofragr::write_contact_matrix(cm, file_path, comments = NULL)
```

`comments; default: comments = NULL`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Allows the addition of a character vector to be appended at the top of the written BED file as a header.

## License

See `LICENSE` for more information.

## Contact

* Haizi Zheng haizi.zh@gmail.com
* Yaping Liu lyping1986@gmail.com

## Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-1" class="csl-entry">

Spatial co-fragmentation pattern of cell-free DNA recapitulates in vivo chromatin organization and identifies tissues-of-origin
Yaping Liu, Tzu-Yu Liu, David E. Weinberg, Brandon W. White, Chris J. De La Torre, Catherine L. Tan, Anthony D. Schmitt, Siddarth Selvaraj, Vy Tran, Louise C. Laurent, Luc Cabel, François-Clément Bidard, Girish Putcha, Imran S. Haque
*bioRxiv* 564773; doi: https://doi.org/10.1101/564773

</div>

</div>
