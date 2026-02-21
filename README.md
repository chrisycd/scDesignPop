## scDesignPop

------------------------------------------------------------------------

scDesignPop is a simulator for population-scale single-cell RNA-sequencing (scRNA-seq) data. By incorporating eQTL effects from genotype data with other covariates, scDesignPop has several key applications: 

- 1) performing eQTL power analysis across cell types, 

- 2) benchmarking single-cell eQTL mapping methods with user-specified eQTLs, 

- 3) protecting genomic privacy by mitigating eQTL-based re-identification of individuals via linking attacks

<span style="color:blue"> **Detailed tutorials that illustrate various functionalities of scDesignPop are available at this [website](https://chrisycd.github.io/scDesignPop/docs/index.html)**</span>. The following illustration figure summarizes the workflow of scDesignPop:

<img src="man/figures/scDesignPop_illustration.jpg" width="600"/>

# Table of contents
1. [Installation](#installation)
2. [Tutorials](#tutorials)
3. [Contact](#contact)
4. [Related Manuscripts](#related-manuscripts)

## Installation<a name="installation"></a>

To install the key dependence of scDesignPop, we recommend:

``` r
if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
}
```

You can install the development version of scDesignPop from
[GitHub](https://github.com) with:

``` r
if(!require("remotes")){
    install.packages("remotes")
}
remotes::install_github("chrisycd/scDesignPop")
```

## Tutorials<a name="tutorials"></a>

For tutorials, please check the [website](https://chrisycd.github.io/scDesignPop/docs/index.html). The tutorials will demonstrate the applications of **scDesignPop** as follows:

-   [Model population-scale scRNA-seq data](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop.html)

-   [Modeling batch effect in population-scale scRNA-seq data](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-batch-effect.html)

-   [Modeling trans-eQTL effect using public data](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-trans-eQTL.html)

-   [Model linear dynamic eQTL effects in continuous cell states](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-dynamic-eQTL.html)

-   [Model nonlinear dynamic eQTL effects in continuous cell states](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-dynamic-eQTL-NL.html)

-   [Model cell type proportions for new individuals](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-celltype-prop-modeling.html)

-   [Modify eQTL effect for eGenes / non-eGenes](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-modify-eQTL-effects.html)
    
-   [Power analysis for selected genes](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-power-analysis-selected.html)

-   [Power analysis based on user specified eQTL effect sizes](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-power-analysis-ES-specification.html)

## Contact<a name="contact"></a>

Any questions or advice on `scDesignPop` are welcomed! Please report it on [issues](https://github.com/chrisycd/scDesignPop/issues), or contact Chris Dong ([cycd\@g.ucla.edu](mailto:cycd@g.ucla.edu){.email}) or Yihui Cen ([yihuicen\@g.ucla.edu](mailto:yihuicen@g.ucla.edu){.email}).


## Related Manuscripts<a name="related-manuscripts"></a>

-   The original **scDesignPop** paper (TBD. Coming soon...)

-   The simulator for general single-cell counts: **scDesign3**
    -   **scDesign3**: [Song, D., Wang, Q., Yan, G. *et al.* scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. *Nat Biotechnol* **42**, 247–252 (2024).](https://www.nature.com/articles/s41587-023-01772-1)
-   The predecessors of **scDesign3**
    -   **scDesign**: [Li, W. V., & Li, J. J. (2019). A statistical simulator scDesign for rational scRNA-seq experimental design. *Bioinformatics*, **35**(14), i41-i50.](https://academic.oup.com/bioinformatics/article/35/14/i41/5529133)
    -   **scDesign2**: [Sun, T., Song, D., Li, W. V., & Li, J. J. (2021). scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured. *Genome biology*, **22**(1), 1-37.](https://link.springer.com/article/10.1186/s13059-021-02367-2)
-   The simulator for single-cell multi-omics reads developed by our lab member Guanao Yan
    -   **scReadSim**: [Yan, G., Song, D. & Li, J.J. scReadSim: a single-cell RNA-seq and ATAC-seq read simulator. *Nat Commun* **14**, 7482 (2023)](https://doi.org/10.1038/s41467-023-43162-w)
