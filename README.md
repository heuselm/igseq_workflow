## igseq_workflow
Workflow to analyze targeted immunoglobulin sequencing data for patient-specific immunoproteogenomics

# Dependencies
Python >= 3.8
Singularity >= 3.6
Snakemake >= 5.24.0

# Usage
```
# optional, if that's how you set up the env
source ~/venvs/[my_Python38_sing36_sm524_env]/bin/activate

# clone wf repository
git clone https://github.com/heuselm/igseq_workflow.git

# descend into wf directory
cd igseq_workflow

# stage input data
copy/move fastq_R1.gz and fastq_R2.gz to 00_input/
define sample_id, master_barcode, slave_barcode, species in header-less, tab-separated barcodes.txt in 00_input/

# run
run ./run.sh
```
