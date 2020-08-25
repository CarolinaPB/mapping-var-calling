# WUR_mappingtry
Install `conda` if you don't have it

### Create conda environment
```
conda create --name <env> --file requirements.txt
```

### Activate environment
```
conda activate <env>
```

### To deactivate the environment (if you want to leave the conda environment)
```
conda deactivate
```

### Create hpc config file ([good example](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/))

Necessary for snakemake to prepare and send jobs

#### Start with creating the directory
```
mkdir -p ~/.config/snakemake/<profile name>
```

#### Add config.yaml to that directory and add the specifications:
```
jobs: 10
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} --job-name={rule} --exclude=fat001,fat002,fat101,fat100 --output=logs_slurm/{rule}.out --error=logs_slurm/{rule}.err --mail-type=[CHOOSE] --mail-user=[EMAIL]"
default-resources: [time_min=180, cpus=16, mem_mb=16000]

use-conda: true
```
(change the options between square brackets)

### Change the paths and choose where the results will be saved (OUTDIR). Snakemake will run in this directory

## How to run
First it's good to always make a dry run: shows if there are any problems with the rules and we can use it to look at the commands and verify that all the fields are in the correct place

Dry run (prints execution plan and commands that will be run)
```
snakemake -np 
```
Run in the HPC with conda environments (necessary for some steps)
```
snakemake --profile <profile name> --use-conda
```

Other flags:
- --forceall : run all the steps, even if it's not needed
- --rerun-incomplete : rerun incomplete steps
- -r [rulename] : run this specific rule
- --max-jobs-per-second \<N> : sometimes there are some problems with the job timings/ many jobs being submitted at once so it's good to choose a low number





