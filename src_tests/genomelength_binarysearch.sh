./build.sh
genome=$1
safile=$2
saplingfile=$3
min=$4
max=$5
outfile=$6
echo '' > $outfile
for i in `seq $min $max`;
do
    pow=$((2**i))
    echo 'Length:' $pow
    echo 'Length:' $pow >> $outfile
    ./genomelength_binarysearch $1 'genomelengthbinarysearch_'$2'_'$pow $pow >> $outfile
done  
