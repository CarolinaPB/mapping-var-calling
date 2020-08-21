# WUR_mappingtry


[Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)


Options: create your own files/ copy from my directory /home/WUR/moiti001/snakemake-try or clone my private repository (https://github.com/CarolinaPB/WUR_mappingtry.git) - probably will need to request access (send me a message with your github username and I'll add you to the repo)


## If you copy my dir/ github repo:
### Create conda environment
```
conda create --name <env> --file requirements.txt
```

### Activate environment
```
conda activate env_name
```
(might be missing a few steps, let me know)

### Create hpc config file ([good example](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/))

Necessary for snakemake to prepare and send jobs

#### Start with creating the directory
```
mkdir -p ~/.config/snakemake/slurm
```

#### Add config.yaml to that directory and add the specifications:
```
jobs: 10
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} --job-name={rule} --output=logs_slurm/{rule}.out --error=logs_slurm/{rule}.err --mail-type=[CHOOSE] --mail-user=[EMAIL]"
default-resources: [time_min=180, cpus=16, mem_mb=16000]

use-conda: true
```
(change the options between square brackets)

### Check if the config.yaml in the snakemake directory has the file paths you want

## How to run
First it's good to always make a dry run: shows if there are any problems with the rules and we can use it to look at the commands and verify that all the fields are in the correct place

Dry run (prints execution plan and commands that will be run)
```
snakemake -np 
```
Run in the HPC with conda environments (necessary for some steps)
```
snakemake --profile slurm --use-conda
```

Other flags:
- --forceall : run all the steps, even if it's not needed
- --rerun-incomplete : rerun incomplete steps
- -r [rulename] : run this specific rule
- --max-jobs-per-second \<N> : sometimes there are some problems with the job timings/ many jobs being submitted at once so it's good to choose a low number (I tried with 2 and it worked well)



Directory tree (at the end of the run)
```
snakemake_try
├── config.yaml
├── envs
│   ├── bwa.yaml
│   ├── freebayes.yaml
│   └── qualimap.yaml
├── fixmate
│   └── DTG-SG-188.fixmate.bam
├── logs_slurm
├── mapped_reads
│   └── DTG-SG-188.bam
├── others
│   └── mapping_bwa.sh
├── README.md
├── requirements.txt
├── results
│   └── qualimap
│       ├── genome_results.txt
│       └── report.pdf
├── scripts
│   ├── bwa_index_check.sh
│   └── mapping.sh
├── Snakefile
├── sorted_reads
│   ├── DTG-SG-188.fixmate.sort.bam
│   ├── DTG-SG-188.fixmate.sort.bam.bai
│   └── DTG-SG-188.fixmate.sort_stats
└── variant_calling
    └── var.vcf.gz
```


## If starting from scratch


### Create directory for the snakemake files and move to that directory

```
mkdir XXX
cd XXX
```

### Create conda environment
```
conda create --name <env> --file requirements.txt
```

### Activate environment
```
conda activate env_name
```

### Create hpc config file ([good example](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/))

#### Start with creating the directory
```
mkdir -p ~/.config/snakemake/slurm
```

#### Add config.yaml to that directory and add the specifications:
```
jobs: 10
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} --job-name={rule} --output=logs_slurm/{rule}.out --error=logs_slurm/{rule}.err --mail-type=[CHOOSE] --mail-user=[EMAIL]"
default-resources: [time_min=180, cpus=16, mem_mb=16000]

use-conda: true
```
### Create a directory for the conda environments
First make sure you're in you snakemake directory
```
mkdir envs
```
### Create directory for the slurm log files
```
mkdir logs_slurm
```
### Create a Snakefile 
```
configfile: "config.yaml"

rule all:
    "Target file"

rule abc:
    input:
        "XXX"
    output:
        "Target file"
    shell:
        "bash commands"
```
It's possible to specify the conda environment for each rule. Place the following after the output
```
conda:
    envs/my_env.yaml
```


## How to run

Dry run (prints execution plan and commands that will be run)
```
snakemake -np 
```
Run in the HPC with conda environments (necessary for some steps)
```
snakemake --profile slurm --use-conda
```

Other flags:
- --forceal : run all the steps, even if it's not needed
- --rerun-incomplete : rerun incomplete steps