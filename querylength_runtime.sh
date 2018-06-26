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
    echo 'Query Length:' $i
    echo 'Query Length:' $i >> $outfile
    ./querylength_runtime $1 'querylengthruntime_'$2'_'$buckets'_'$i 'querylengthruntime_'$3'_'$buckets'_'$i $buckets $i >> $outfile
done  
