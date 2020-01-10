if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

WORKINGDIR=`pwd`

javac $BINDIR/*.java
java -cp $BINDIR SuffixArraySimulatedSequences
python $BINDIR/plotSimSa.py gc.png
python $BINDIR/plotSimSa.py repeats.png
