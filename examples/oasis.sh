# convert .hdr/.img to .nifti
for i in $(find ~/Desktop/oasis/*/RAW/ -type f -name "*.hdr")
do 
sub=$(dirname $i)
id=$(dirname $sub)
echo $id
dcm2nii -a n -g n -r n -x n -o id/nifti/*.nii $i
done

# Organize oasis brains
for i in $(find ~/Desktop/oasis/*/RAW/ -type d)
do
sub=$(dirname $i)
cd $sub
num=1
for scan in $(find RAW/*.nii -type f)
  do
  filename=$sub/nifti/t1_$num
  mkdir $filename
  echo $num
  cp $scan $filename/t1.nii
  num=$(($num+1))
  done
done

#ACPC detect



antsBrainExtraction.sh -d 3 -a $i
