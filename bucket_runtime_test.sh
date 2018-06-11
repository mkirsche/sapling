./build.sh
genome=$1
safile=$2
saplingfile=$3
outfile=$4
echo '' > $outfile
for i in `seq 5 25`;
do
    echo 'Buckets:' $i
    ./sapling $1 $2 $3'_'$i $i >> $outfile
done 

cat 
