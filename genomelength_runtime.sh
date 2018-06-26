./build.sh
genome=$1
safile=$2
saplingfile=$3
min=$4
max=$5
outfile=$6
buckets=15
echo '' > $outfile!
for i in `seq $min $max`;
do
    rm $2'_'$buckets'_'$pow
    rm $3'_'$buckets'_'$pow
    pow=$((2**i))
    echo 'Length:' $pow
    echo 'Length:' $pow >> $outfile
    ./genomelength_runtime $1 $2'_'$buckets'_'$pow $3'_'$buckets'_'$pow $buckets $pow >> $outfile
done  
