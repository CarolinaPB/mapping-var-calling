# WUR_mappingtry


[Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)


Options: create your own files/ copy from my directory /home/WUR/moiti001/snakemake-try or clone my private repository (https://github.com/CarolinaPB/WUR_mappingtry.git) - probably will need to request access (send me a message with your github username and I'll add you to the repo)


## If you copy my dir/ github repo:
### Create conda environment
```
conda create --name env_name
```

### Activate environment
```
conda activate env_name
```
(might be missing a few steps, let me know)

### Create hpc config file ([good example](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/))

#### Start with creating the directory
```
mkdir -p ~/.config/snakemake/slurm
```

#### Add config.yaml to that directory and add the specifications:
```
jobs: 10
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} --output=logs_slurm/{rule}.out --error=logs_slurm/{rule}.err --mail-type=[CHOOSE] --mail-user=[EMAIL]"
default-resources: [time_min=180, cpus=16, mem_mb=16000]

use-conda: true
```
(change the options between square brackets)

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

Directory tree
```
snakemake-try/
├── checks
├── config.yaml
├── envs
│   ├── bwa.yaml
│   ├── freebayes.yaml
│   └── qualimap.yaml
├── fixmate
├── logs_slurm
├── mapped_reads
├── others
├── README.md
├── scripts
│   ├── bwa_index_check.sh
│   └── mapping.sh
├── Snakefile
├── sorted_reads
└── variant_calling
```


## If starting from scratch


### Create directory for the snakemake files and move to that directory

```
mkdir XXX
cd XXX
```

### Create conda environment
```
conda create --name env_name
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
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} --output=logs_slurm/{rule}.out --error=logs_slurm/{rule}.err --mail-type=[CHOOSE] --mail-user=[EMAIL]"
default-resources: [time_min=180, cpus=16, mem_mb=16000]

use-conda: true
```
### Create a directory for the conda environments
First make sure you're in you snakemake directory
```
mkdir envs
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