# SMR-IVW

SMR-IVW is an extension of the [SMR](https://cnsgenomics.com/software/smr/) software that additionally allows to perform univariable and multivariable inverse-variance weighting (IVW) Mendelian randomization (MR) analyses on omics and GWAS summary data.

## Citations

If you use this software please consider citing the original multi-SNP-based SMR method and software:

Wu Y, Zeng J, Zhang F, Zhu Z, Qi T, Zheng Z, Lloyd-Jones LR, Marioni RE, Martin NG, Montgomery GW, Deary IJ, Wray NR, Visscher PM, McRae AF & Yang J (2018) Integrative analysis of omics summary data reveals putative mechanisms underlying complex traits. [Nature Communications, 9: 918](https://doi.org/10.1038/s41467-018-03371-0)

as well as the extension to univariable and multivariable MR-IVW analyses:

Quantifying mediation between omics layers and complex traits
Sadler MC, Auwerx C, Porcu E & Kutalik Z (2021) Quantifying mediation between omics layers and complex traits. [bioRxiv 2021.09.29.462396]( https://doi.org/10.1101/2021.09.29.462396)

## Installation

The software works on Unix / Linux platforms. Simply download the latest release [here](https://github.com/masadler/smrivw/releases).

Alternatively compile the package making sure all the paths (e.g. Eigen library) point to the right location on your system.

## Usage

Extensive documentation on data management functions (e.g. formatting input files, filtering input files) can be found on the official [SMR website](https://cnsgenomics.com/software/smr/). 

Documentation of univariable and multivariable MR-IVW is provided on this [Github-Wiki](https://github.com/masadler/smrivw/wiki).

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[GNU GPL version 3](https://www.gnu.org/licenses/gpl-3.0.en.html)
