#! /bin/bash

#####################################################
#                                                   #
# Here is a little bash script that runs the code   #
# in an exemplary manner to wring out any changes   #
# we have made. Ideally this is run before a commit #
# to sanity check that the code behaves the same    #
# way between versions. Must have bash.             #
#                                                   #
#####################################################

## counters for the tests
tests=0
failures=0

## print out our name
function corr {
    echo ""
    echo "CORR CORR CORR     CORR CORR CORR     CORR CORR CORR     CORR CORR CORR"
    echo "CORR CORR CORR     CORR CORR CORR     CORR CORR CORR     CORR CORR CORR"
    echo "CORR               CORR      CORR     CORR      CORR     CORR      CORR"
    echo "CORR               CORR      CORR     CORR CORR          CORR CORR"
    echo "CORR               CORR      CORR     CORR   CORR        CORR   CORR"
    echo "CORR CORR CORR     CORR CORR CORR     CORR      CORR     CORR      CORR"
    echo "CORR CORR CORR     CORR CORR CORR     CORR      CORR     CORR      CORR"
    echo ""
}

## checksums
function check_checksums {
    cksum="$( cat output.txt | grep "CHECKSUM" | cut -d" " -f4 )"
    let 'tests=tests+1'
    if [ $cksum != "passed" ] ; then
       echo "[CHECKSUMS] checksums of $1 do not match!"
       let 'failures=failures+1'
    fi
}

function baryon_test {
    ## simple test of first element of matrix
    first="$(./MESONS $1 0,0,0,0,0 | grep "CORR 0" | cut -d" " -f3)"
    let 'tests=tests+1'
    if [ $first != $2 ] ; then
       echo "[BARYONS] $1 not equal to known result $first != $2"
       let 'failures=failures+1'
    fi
    ## project a parity and check xyz degeneracy for the first six digits
    ./BARYONS $1 CHIRAL NONE L5 false 0 4,4,4,8 tmp > output.txt
    check_checksums output.txt
    x="$(./MESONS tmp 0,0,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    y="$(./MESONS tmp 1,1,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    z="$(./MESONS tmp 2,2,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    let 'tests=tests+1'
    if [ $x != $y ] || [ $x != $z ] || [ $y != $z ] ; then
      echo "[BARYONS] $1 non xyz degeneracy in parity projected (x,y,z) :: ($x,$y,$z)"
      let 'failures=failures+1'
    fi
    ## clean up
    rm tmp output.txt
}

## test our baryon output
function test_baryons {
    echo "[BARYONS] Testing uuu type"
    baryon_test baryon.uuu 9.838464389472e-01

    echo "[BARYONS] Testing uud type"
    baryon_test baryon.uud 4.088883194057e-01

    echo "[BARYONS] Testing uds type"
    baryon_test baryon.uds 3.258534193379e-01

    ## clean up
    rm baryon.uuu baryon.uud baryon.uds
}

## bash script to check our output
function test_mesons {
    echo "[MESONS] Testing flavour diagonal meson code"

    ## compute the pseudoscalar at some timeslices
    ./MESONS meson1 5,5,0,0,0 >tmp
    first="$(cat tmp | grep "CORR 0" | cut -d" " -f3)"
    second="$(cat tmp | grep "CORR 1" | cut -d" " -f3)"
    last="$(cat tmp | grep "CORR 7" | cut -d" " -f3)"
    let 'tests=tests+1'
    if [ $first != "8.583031475225e-01" ] ; then
      echo "[MESONS] Error $first != 8.583031475225e-01 "
      let 'failures=failures+1'
    fi
    ## test that t=1 == t=7
    let 'tests=tests+1'
    if [ $second != $last ]; then 
      echo "[MESONS] Error $second is not the same as $last"
      let 'failures=failures+1'
    fi
    ## make sure that the flavour off diagonal code gives
    ## the same result as flavour diagonal for same props
    echo "[MESONS] Testing flavour off diagonal code"
    ./MESONS meson2 5,5,0,0,0 >tmp
    first2="$(cat tmp | grep "CORR 0" | cut -d" " -f3)"
    second2="$(cat tmp | grep "CORR 1" | cut -d" " -f3)"
    last2="$(cat tmp | grep "CORR 7" | cut -d" " -f3)"
    let 'tests=tests+1'
    if [ $first != $first2 ] || [  $second != $second2 ] ; then
       echo "[MESONS] Error flavour off diagonal not flavour diagonal!"
       let 'failures=failures+1'
    fi
    ## finally test that x y z are degenerate
    ./MESONS meson1 0,0,0,0,0 >tmp
    x="$(cat tmp | grep "CORR 4" | cut -d" " -f3)"
    ./MESONS meson1 1,1,0,0,0 >tmp
    y="$(cat tmp | grep "CORR 4" | cut -d" " -f3)"
    ./MESONS meson1 2,2,0,0,0 >tmp
    z="$(cat tmp | grep "CORR 4" | cut -d" " -f3)"
    let 'tests=tests+1'
    if [ $x != $y ] || [ $y != $z ] || [ $x != $z ] ; then
       echo "[MESONS] xyz are not degenerate! (x,y,z) -> ($x,$y,$z)"
       let 'failures=failures+1'
    fi
    ## clean up our temporaries
    rm tmp meson1 meson2
}

## test degeneracy for tetraquark contractions
function tetra_test {
    ## diquark-diquark
    ## project a parity and check xyz degeneracy for the first six digits
    x="$(./MESONS $1 0,0,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    y="$(./MESONS $1 0,1,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    z="$(./MESONS $1 0,2,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    let 'tests=tests+1'
    if [ $x != $y ] || [ $x != $z ] || [ $y != $z ] ; then
      echo "[TETRA] $1 non xyz degeneracy in diquark (x,y,z) :: ($x,$y,$z)"
      let 'failures=failures+1'
    fi
    ## dimeson
    ## project a parity and check xyz degeneracy for the first six digits
    x="$(./MESONS $1 1,0,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    y="$(./MESONS $1 1,1,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    z="$(./MESONS $1 1,2,0,0,0 | grep "CORR 4" | cut -d" " -f3 | cut -b 2,3,4,5,6,7,8,9)"
    let 'tests=tests+1'
    if [ $x != $y ] || [ $x != $z ] || [ $y != $z ] ; then
      echo "[TETRA] $1 non xyz degeneracy in dimeson (x,y,z) :: ($x,$y,$z)"
      let 'failures=failures+1'
    fi
}

## test our tetraquark contraction code
function test_TETRA {
    echo "[TETRA] Testing tetraquark contractions"
    
    echo "[TETRA] Testing udbb degeneracy"
    tetra_test tet.udbb

    echo "[TETRA] Testing udcb degeneracy"
    tetra_test tet.udcb

    echo "[TETRA] Testing usbb degeneracy"
    tetra_test tet.usbb

    echo "[TETRA] Testing uscb degeneracy"
    tetra_test tet.uscb

    ## clean up
    rm tet.udbb tet.udcb tet.usbb tet.uscb
}

## test our vpfs output, test WIs are sane and then look at the actual data
function test_VPF {
    echo "[VPF] testing vacuum polarisation function code"
    echo "[VPF] Testing configuration space WI"
    cspace="$(cat runout.txt | grep "backward config-space violation" | cut -d" " -f5)"
    let 'tests=tests+1'
    c1="$(echo $cspace | cut -d" " -f1)"
    c2="$(echo $cspace | cut -d" " -f2)"
    if [ $c1 != $c2 ] ; then
	echo "[VPF] diagonal and off diagonal conserved-local are not the same $c1,$c2"
	let 'failures=failures+1'
    fi
    ## test p_mu ~0
    echo "[VPF] Testing momentum space WI"
    pmuPimunu="$(cat runout.txt | grep "p_{mu}" | cut -d" " -f6,7,8)"
    p1="$(echo $pmuPimunu | cut -d" " -f1)"
    p2="$(echo $pmuPimunu | cut -d" " -f3)"
    let 'tests=tests+1'
    if [ $p1 != $p2 ] ; then
	echo "[VPF] p_mu Pimunus are different between diagonal and off $p1,$p2"
	let 'failures=failures+1'
    fi
    p1="$(echo $pmuPimunu | cut -d" " -f2)"
    p2="$(echo $pmuPimunu | cut -d" " -f4)"
    let 'tests=tests+1'
    if [ $p1 != $p2 ] ; then
	echo "[VPF] p_mu Pimunus are different between diagonal and off $p1,$p2"
	let 'failures=failures+1'
    fi
    if [ test -a cl1.CVLV.trans.bin ] ; then
	echo "[VPF] Can't test projection for cl1.CVLV.trans.bin"
    else
	## vpfread the stuff
	zero="$(./VPFREAD cl1.CVLV.trans.bin 4,4,4,8 | grep "0.000000" | cut -d" " -f14)"
	let 'tests=tests+1'
	## this thing is a garbage value because of the projection
	if [ $zero != "0.220843" ] ; then
	    echo "[VPF] zero momentum point does not match previous code $zero,0.220843"
	    let 'failures=failures+1'
	fi
    fi
    ## clean up the files it creates
    rm cl1.* cl2.*
}

## little ascii art for our users, who am I kidding? for me
corr

## run the code :: pipe to runout.txt
./CORR -i infile -c lat..400 >runout.txt

## test some select values are correct
test_baryons
test_mesons
test_TETRA
test_VPF

## clean up the output
rm runout.txt

printf "\n[REGRESSION] Tests run :: $tests || tests failed :: $failures\n"
