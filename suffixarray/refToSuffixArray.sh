BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
echo 'BINDIR '$BINDIR

cd $BINDIR
git submodule update --init --recursive
cd libdivsufsort
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_INSTALL_PREFIX="/usr/local" ..
sed -i 's/int32_t/int64_t/g' include/divsufsort.h
make

cd $WORKINGDIR

ref_untrimmed=$1
ref=$ref_untrimmed'.ref'
if [ ! -r ref ]
then
  g++ $BINDIR/trimRef.cpp -o $BINDIR/trimRef
  $BINDIR/trimRef $ref_untrimmed
fi
final_sa=$ref_untrimmed'.sa'
unrefined_sa=$ref'.sa'

echo 'Trimmed reference: ' $ref
$BINDIR/libdivsufsort/build/examples/mksary $ref $unrefined_sa

g++ $BINDIR/addlcp.cpp -o $BINDIR/addlcp
echo 'final_sa: '$final_sa
if [ ! -r $final_sa ]
then
  $BINDIR/addlcp $ref $unrefined_sa $final_sa
fi
