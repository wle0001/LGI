#!/bin/bash
### SLURM batch script

#execute this script as following:
# sbatch this.script.path/this.script.name

# You can call it from cron. Type 'crontab -e' to edit and see example

### Email address
#SBATCH --mail-user=wle0001@uah.edu

### Job name
#SBATCH -J LGI_Batch

### Partition (queue), select shared or standard
#SBATCH -p standard

### TOTAL processors (number of tasks)
#SBATCH --ntasks 1

### total run time estimate (D-HH:MM)
#SBATCH -t 5-01:00

### allocated memory for EACH processor (GB))
### NOTE: TOTAL MEMORY ALLOCATED = mem-per-cpu X ntasks
### Be careful about making this number large, especially when using a large number of processors (ntasks)
#SBATCH --mem-per-cpu=98G

### Mail to user on an event
### common options are FAIL, BEGIN, END, REQUEUE
#SBATCH --mail-type=END,FAIL

### Ouput files
#SBATCH -o slurm1.out # STDOUT
#SBATCH -e slurm1.err # STDERR


echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo

### Add your module statments here

module load python/v3


# Source python environment

conda activate nwm

# run python scripts

#example:
 
python LGI.py
