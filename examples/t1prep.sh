

export ANTSPATH=$rootdir/bin/antsbin/bin/
PATH=${ANTSPATH}:${PATH}

rootdir=/fslhome/zach8769/

for command

subjDir=$(dirname ${array[$i]})
mkdir $subjDir/atlas
N4BiasFieldCorrection -d 3 -i ${array[$i]} -o ${subjDir}/n4.nii.gz -s 8 -b [200] -c [50x50x50x50,0.000001]
N4BiasFieldCorrection -d 3 -i ${subjDir}/n4.nii.gz -o ${subjDir}/n4.nii.gz -s 4 -b [200] -c [50x50x50x50,0.000001]
N4BiasFieldCorrection -d 3 -i ${subjDir}/n4.nii.gz -o ${subjDir}/n4.nii.gz -s 2 -b [200] -c [50x50x50x50,0.000001]

/fslhome/zach8769/bin/c3d ${subjDir}/n4.nii.gz -resample-mm 1x1x1mm -o ${subjDir}/n4_resliced.nii.gz

antsBrainExtraction.sh -d 3 \
-a $rootdir/compute/NKI10AndUnder/T_template0.nii.gz \
-m $rootdir/compute/NKI10AndUnder/T_template0_BrainExtractionMask.nii.gz
-o $subjDir/extract

temp=$rootdir/compute/jlf/OASIS-TRT-20
antsJointLabelFusion.sh -d 3 -o $subjDir/atlas/ -t $subjDir/extract*.nii.gz -c 5 \
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

