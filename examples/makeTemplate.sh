inputPath=/fslhome/zach8769/compute/oasistemplate/

${ANTSPATH}/buildtemplateparallel.sh \
  -d 2 \
  -o T_ \
  -i 4 \
  -g 0.25 \
  -j 4 \
  -c 5 \
  -m 100x70x50x10 \
  -n 1 \
  -r 1 \
  -s CC \
  -t GR 
