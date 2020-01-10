g++ -o count count.cpp
for i in `ls *.fa`
do
  echo $i
  bn=`basename $i .fa`
  ./count $i | tee $bn
done
python plot.py
