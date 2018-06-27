# Metagen-o-matic
## Pipeline for automated analysis of metagenomic data
 
### Whole Metagenome Sequencing (WMS) data

1. ./pipe_wgs_preproc.sh
	- Description:
		Pre-processing of reads (adapter and quality trimming)
	- Parameters:
		- Input directory
		- An identifier for the analysis
		- Output directory
	- Example:
	```
	./pipe_wgs_preproc.sh ./indir example ./procdir > ./pipe_wgs_preproc.out.txt 2> ./pipe_wgs_preproc.err.txt
	```


2. pipe_wgs_asm_pe.sh
	- Description:
		Assemble reference metagenome sequences based on sequencing reads data.
	- Parameters:
		- Input directory
		- Regular expression used to capture sample name patterns and group samples from biological replicates of the same sample groups
		- An identifier for the analysis
	- Example:
	```
	./pipe_wgs_asm_pe.sh ./procdir/example '^(sample.)' example > ./pipe_wgs_asm_pe.out.txt 2> ./pipe_wgs_asm_pe.err.txt
	```

3. pipe_wgs_map.sh
	- Description:
		Taxonomic and functional annotation
	- Parameters:
		- Input directory
		- An identifier for the analysis
	- Example:
	```
	./pipe_wgs_map.sh ./procdir/example example > ./pipe_wgs_map.out.txt 2> ./pipe_wgs_map.err.txt
	```

