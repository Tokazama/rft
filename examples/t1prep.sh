#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8192M
#SBATCH -o /fslhome/zach8769/logfiles/sobik/test_output.txt
#SBATCH -e /fslhome/zach8769/logfiles/sobik/test_error.txt
#SBATCH -J "jlf"
#SBATCH --mail-user=zchristensen7@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

rootdir=/fslhome/zach8769/

ARTHOME=$rootdir/bin/art
export ARTHOME

export ANTSPATH=$rootdir/bin/antsbin/bin/
PATH=${ANTSPATH}:${PATH}

IFS=$'\n'
array=( $(find /fslhome/zach8769/compute/sobic/*/t1/ -type f  -name "acpc.nii") )
for i in 0; do
t1=$(dirname ${array[$i]})
subjDir=$(dirname $t1)
mkdir $subjDir/atlas

N4BiasFieldCorrection -d 3 -i ${array[$i]} -o ${t1}/n4.nii.gz -s 8 -b [200] -c [50x50x50x50,0.000001]
N4BiasFieldCorrection -d 3 -i ${t1}/n4.nii.gz -o ${t1}/n4.nii.gz -s 4 -b [200] -c [50x50x50x50,0.000001]
N4BiasFieldCorrection -d 3 -i ${t1}/n4.nii.gz -o ${t1}/n4.nii.gz -s 2 -b [200] -c [50x50x50x50,0.000001]

/fslhome/zach8769/bin/c3d ${t1}/n4.nii.gz -resample-mm 1x1x1mm -o ${t1}/n4_resliced.nii.gz

antsBrainExtraction.sh -d 3 \
-a ${t1}/n4_resliced.nii.gz
-e $rootdir/compute/NKI10AndUnder/T_template0.nii.gz \
-m $rootdir/compute/NKI10AndUnder/T_template0_BrainExtractionMask.nii.gz
-o $t1/extract

temp=$rootdir/compute/jlf/OASIS-TRT-20
antsJointLabelFusion.sh -d 3 -o ${subjDir}/atlas/ -t ${t1}/extract*.nii.gz -c 5 \
-g $temp-1.nii.gz -l $temp-1_DKT31_CMA_labels.nii.gz \
-g $temp-2.nii.gz -l $temp-2_DKT31_CMA_labels.nii.gz \
-g $temp-3.nii.gz -l $temp-3_DKT31_CMA_labels.nii.gz \
-g $temp-4.nii.gz -l $temp-4_DKT31_CMA_labels.nii.gz \
-g $temp-5.nii.gz -l $temp-5_DKT31_CMA_labels.nii.gz \
-g $temp-6.nii.gz -l $temp-6_DKT31_CMA_labels.nii.gz \
-g $temp-7.nii.gz -l $temp-7_DKT31_CMA_labels.nii.gz \
-g $temp-8.nii.gz -l $temp-8_DKT31_CMA_labels.nii.gz \
-g $temp-9.nii.gz -l $temp-9_DKT31_CMA_labels.nii.gz \
-g $temp-10.nii.gz -l $temp-10_DKT31_CMA_labels.nii.gz \
-g $temp-11.nii.gz -l $temp-11_DKT31_CMA_labels.nii.gz \
-g $temp-12.nii.gz -l $temp-12_DKT31_CMA_labels.nii.gz \
-g $temp-13.nii.gz -l $temp-13_DKT31_CMA_labels.nii.gz \
-g $temp-14.nii.gz -l $temp-14_DKT31_CMA_labels.nii.gz \
-g $temp-15.nii.gz -l $temp-15_DKT31_CMA_labels.nii.gz \
-g $temp-16.nii.gz -l $temp-16_DKT31_CMA_labels.nii.gz \
-g $temp-17.nii.gz -l $temp-17_DKT31_CMA_labels.nii.gz \
-g $temp-18.nii.gz -l $temp-18_DKT31_CMA_labels.nii.gz \
-g $temp-19.nii.gz -l $temp-19_DKT31_CMA_labels.nii.gz \
-g $temp-20.nii.gz -l $temp-20_DKT31_CMA_labels.nii.gz
done

