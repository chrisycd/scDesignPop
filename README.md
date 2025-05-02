# scDesignPop

------------------------------------------------------------------------

scDesignPop is a simulator for population-scale single-cell RNA-sequencing (scRNA-seq) data. By incorporating eQTL effects from genotype data with other covariates, scDesignPop has several key applications: 

- 1) performing eQTL power analyses at cell-type resolution, 

- 2) protecting genomic privacy by mitigating eQTL-based re-identification of individuals via linking attacks, 

- 3) simulating scRNA-seq data for new individuals using either simulated or real genotype data, and 

- 4) generating positive- and negative-control data for cell-type-specific eQTLs.

<span style="color:blue"> **Detailed tutorials that illustrate various functionalities of scDesignPop are available at this [website](https://chrisycd.github.io/scDesignPop/docs/index.html)**</span>. The following illustration figure summarizes the workflow of scDesignPop:

<img src="man/figures/scDesignPop_illustration.jpg" width="600"/>

# Table of contents
1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Tutorials](#tutorials)
4. [Contact](#contact)
5. [Related Manuscripts](#related-manuscripts)

## Installation<a name="installation"></a>

You can install the development version of scDesignPop from
[GitHub](https://github.com) with:

``` r
# install.packages("remotes")
remotes::install_github("chrisycd/scDesignPop")
```

## Quick Start<a name="quick-start"></a>

## Tutorials<a name="tutorials"></a>

For all detailed tutorials, please check the [website](https://chrisycd.github.io/scDesignPop/docs/index.html). The tutorials will demonstrate the applications of **scDesignPop** from the following perspectives: preprocessing, data simulation, model alteration, and power analysis.

-   Preprocessing

-   Data simulation

-   Model alteration
    -   [Modify eQTL effect for eGenes / non-eGenes](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-modify-eQTL-effect.html)
    
-   Power Analysis
    -   [Power analysis based on a fitted scDesignPop marginal model](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-power-analysis-fitted.html)
    -   [Power analysis for selected genes](https://chrisycd.github.io/scDesignPop/docs/articles/scDesignPop-power-analysis-selected.html)

## Contact<a name="contact"></a>

Any questions or advice on `scDesignPop` are welcomed! Please report it on [issues](https://github.com/chrisycd/scDesignPop/issues), or contact Chris Dong ([cycd\@g.ucla.edu](mailto:cycd@g.ucla.edu){.email}) or Yihui Cen ([yihuicen\@g.ucla.edu](mailto:yihuicen@g.ucla.edu){.email}).


## Related Manuscripts<a name="related-manuscripts"></a>

-   The original **scDesignPop** paper

-   The simulator for general single-cell counts: **scDesign3**
    -   **scDesign3**: [Song, D., Wang, Q., Yan, G. *et al.* scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. *Nat Biotechnol* **42**, 247–252 (2024).](https://www.nature.com/articles/s41587-023-01772-1)
-   The predecessors of **scDesign3**
    -   **scDesign**: [Li, W. V., & Li, J. J. (2019). A statistical simulator scDesign for rational scRNA-seq experimental design. *Bioinformatics*, **35**(14), i41-i50.](https://academic.oup.com/bioinformatics/article/35/14/i41/5529133)
    -   **scDesign2**: [Sun, T., Song, D., Li, W. V., & Li, J. J. (2021). scDesign2: a transparent simulator that generates high-fidelity single-cell gene expression count data with gene correlations captured. *Genome biology*, **22**(1), 1-37.](https://link.springer.com/article/10.1186/s13059-021-02367-2)
-   The simulator for single-cell multi-omics reads developed by our lab member Guanao Yan
    -   **scReadSim**: [Yan, G., Song, D. & Li, J.J. scReadSim: a single-cell RNA-seq and ATAC-seq read simulator. *Nat Commun* **14**, 7482 (2023)](https://doi.org/10.1038/s41467-023-43162-w)
