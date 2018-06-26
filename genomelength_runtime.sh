./build.sh
genome=$1
safile=$2
saplingfile=$3
min=$4
max=$5
outfile=$6
buckets=15
echo '' > $outfile
for i in `seq $min $max`;
do
    pow=$((2**i))
    echo 'Length:' $pow
    echo 'Length:' $pow >> $outfile
    ./genomelength_runtime $1 'genomelengthruntime_'$2'_'$buckets'_'$pow 'genomelengthruntime_'$3'_'$buckets'_'$pow $buckets $pow >> $outfile
done  
