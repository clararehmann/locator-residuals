# SLiMulation

Running Locator on simulated, randomly-dispersing populations -
*how does the spatial distribution of training samples affect prediction bias?*

Simulations were run using `scripts/biased_migration.slim` with a `BIAS` parameter of 0.0;
output tree sequences were processed using `scripts/process_treeseq.py` and VCF genotype files were then converted to Zarr format.
