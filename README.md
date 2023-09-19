# Welcome to EISA-EXPOSOME

EISA-EXPOSOME platform for high throughput suspect chemical screening, which can help reduce the dependence on chemical standards in traditional chemical analysis and significantly enhance the chemical
coverage. This algorithm used a targeted peak extraction strategy to rescue features that cannot
be extracted by traditional peak extraction algorithms, improving the coverage of chemica  substance annotation.
### The workflow of EISA-EXPOSOME is as follows
![enter image description here](https://github.com/Lab-XUE/EISA-EXPOSOME/blob/main/Figure/workflow.png)

# Database for EISA-EXPOSOME

If you are building your own database, your file (.xlsx /.csv) must contain the following columns:|NAME|PrecursorMZ|ProductMZ|Intensity|RT|ID|, **RT** is not essential. 
We also provide the compiled T3DB database file in .xlsx format.
|NAME|PrecursorMZ|ProductMZ|Intensity|RT|ID|
|-|-|-|-|-|-|
|Methamidophos|142.0086|94.0046|100|2.182|1|
|Methamidophos|142.0086|124.9816|27|2.182|1|


# Shiny
We provide a Rshiny program for EISA-EXPOSOME, which runs with the interface shown below, and you can filter the results according to the visualisation interface！
![enter image description here](https://github.com/Lab-XUE/EISA-EXPOSOME/blob/main/Figure/shiny.tif)


# Literature related to EISA technology

**EISA** is a computational tool for exposome targered extraction and annotation based on databases. The following literature can help you better understand EISA technology.
- [Quantitative multiple fragment monitoring with enhanced in-source fragmentation/annotation mass spectrometry](https://www.nature.com/articles/s41596-023-00803-0)
- [Single Quadrupole Multiple Fragment Ion Monitoring Quantitative Mass Spectrometry](https://pubs.acs.org/doi/full/10.1021/acs.analchem.1c01246)
- [Enhanced in-Source Fragmentation Annotation Enables Novel Data Independent Acquisition and Autonomous METLIN Molecular Identification](https://pubs.acs.org/doi/10.1021/acs.analchem.0c00409)
- [Proteomics with Enhanced In-Source Fragmentation/Annotation: Applying XCMS-EISA Informatics and Q-MRM High-Sensitivity Quantification](https://pubs.acs.org/doi/10.1021/jasms.1c00188)
