if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

WORKINGDIR=`pwd`

fn=$WORKINGDIR/$1
highlightfile=$fn'.highlights'

javac $BINDIR/*.java
java -cp $BINDIR BestAndWorstBins $fn $highlightfile 10
python $BINDIR/hist.py $highlightfile
