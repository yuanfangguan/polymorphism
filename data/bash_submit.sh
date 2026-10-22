#!/bin/bash
                                       ## REQUIRED: #!/bin/bash must be on the 1st line
                                       ## and it must be the only string on the line
#SBATCH --job-name=coinflip            ## Name of the job for the scheduler
#SBATCH --account=bioinf585w23_class ## Generally your PI's uniqname will go here
##SBATCH --partition=standard           ## name of the queue to submit the job to.
#SBATCH --partition=gpu ## name of the queue to submit the job to.
                                       ## (Choose from: standard, debug, largemem, gpu,
   ##  spgpu) 
#SBATCH --gpus=1  			   ## if partition=gpu, number of GPUS needed
   ## make the directive = #SBATCH, not ##SBATCH 
#SBATCH --nodes=1                      ## number of nodes you are requesting
#SBATCH --ntasks=1                     ## how many task spaces do you want to reserve
#SBATCH --cpus-per-task=3              ## how many cores do you want to use per task
#SBATCH --time=1:30:00                 ## Maximum length of time you are reserving the 
   ## resources for (bill is based on time used)
#SBATCH --mem=128g ## Memory requested per core
#SBATCH --mail-user=gyuanfan@umich.edu  ## send email notifications to umich email listed
#SBATCH --mail-type=END                ## when to send email (standard values are:
                                       ## NONE,BEGIN,END,FAIL,REQUEUE,ALL.
   ## (See documentation for others)
#SBATCH --output=./%x-%j               ## send output and error info to the file listed

##(optional: different name format than default) 

module load tensorflow/2.7.0

my_job_header

#python3 train_early_stop.py
#python3 pickle_chr_bit.py NC_000001.11
#python3 pickle_chr_bit.py NC_000002.12
#python3 pickle_chr_bit.py NC_000003.12
python3 pickle_chr_text.py NC_000001.11



