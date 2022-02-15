# SMR-IVW

SMR-IVW is an extension of the [SMR](https://cnsgenomics.com/software/smr/) software that additionally allows to perform univariable and multivariable inverse-variance weighting (IVW) Mendelian randomization (MR) analyses on omics and GWAS summary data.

## Citations

If you use this software please consider citing the original multi-SNP-based SMR method and software:

Wu Y, Zeng J, Zhang F, Zhu Z, Qi T, Zheng Z, Lloyd-Jones LR, Marioni RE, Martin NG, Montgomery GW, Deary IJ, Wray NR, Visscher PM, McRae AF & Yang J (2018) Integrative analysis of omics summary data reveals putative mechanisms underlying complex traits. [Nature Communications, 9: 918](https://doi.org/10.1038/s41467-018-03371-0)

as well as the extension to univariable and multivariable MR-IVW analyses:

Sadler MC, Auwerx C, Porcu E & Kutalik Z (2021) Quantifying mediation between omics layers and complex traits. [bioRxiv 2021.09.29.462396]( https://doi.org/10.1101/2021.09.29.462396)

## Installation

The software works on Unix / Linux platforms (tested on CentOS Linux release 7.9.2009 (Core)). Download the latest release of the executable [here](https://github.com/masadler/smrivw/releases) (about 2MB < 1min download time). 

Then unzip (if necessary), rename the software (optional) and change permissions if necessary as follows:

```bash
gzip -d smrivw_1.0.gz
mv smrivw_1.0 smrivw
chmod +rwx smrivw
```

Alternatively compile the package making sure all the paths (e.g. Eigen library) point to the right location on your system.

## Usage

Extensive documentation on data management functions (e.g. formatting input files, filtering input files) can be found on the official [SMR website](https://cnsgenomics.com/software/smr/). 

Documentation of univariable and multivariable MR-IVW is provided on this [Github-Wiki](https://github.com/masadler/smrivw/wiki).

As a starter run the examples inside the `examples/` folder as follows (if necessary change the path of the executable):

```bash
bash run_examples.sh
```

The expected output is shown in the `examples/expected_output/` folder. The first example computes the univariable and multivariable MR effects of the DNA methylation probe cg07689907 on Triglycerides (TRI) via the transcript *FADS2*. The second examples computes the same MR causal effects but with the output being Apolipoprotein A (APOA). Details on the column headers and used formulas are found in the [Multivariable MR (Mediation analysis) Wiki Page](https://github.com/masadler/smrivw/wiki/Multivariable-MR-(Mediation-analysis)). The example data contains a subset of the 1000 Genomes (European ancestry) reference panel, DNA methylation QTL effects from the [GoDMC consortium](http://mqtldb.godmc.org.uk/downloads), eQTL effects from the [eQTLGen consortium](https://www.eqtlgen.org/cis-eqtls.html) and GWAS effect sizes from the [UK Biobank](http://www.nealelab.is/uk-biobank). Note that in order to reproduce results in the manuscript, it is recommended to use the [UK10K reference panel](https://doi.org/10.1038/nature14962) which better matches allele frequencies of the other datasets. The expected run time should be less than 1 minute.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[GNU GPL version 3](https://www.gnu.org/licenses/gpl-3.0.en.html)
