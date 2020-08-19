# WUR_mappingtry

## To use:
[Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

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

### Create hpc config file [good example](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)

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
Add your config file at the top 
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
With 
```
conda:
    envs/my_env.yaml
```
you can specify the conda environment for each rule
