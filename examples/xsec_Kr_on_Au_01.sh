# calculate excitation cross section and visualize the b-dependent excitation probability
./nuco  \
	--Ap=80 --Zp=36 --At=197 --Zt=79 \
	--E=60 \
	--b=10 --bmax=200 --cross-section-integral-steps=100 \
	--method=relativistic --accuracy=1e-7 \
	--calc-cross-section \
	--levelE=0 --levelI=0 \
	--levelE=.6166 --levelI=2 \
	--MEfrom=0 --MEto=1 --MElambda=2 --MEvalue=61.8

cat << EOF > visualize.gnuplot
set xlabel "b [fm]"
set ylabel "amplitude"
set grid 
set term wxt size 600,1000
set key top right box opaque spacing 1.5 width 1
set multiplot layout 3,1
plot \
	"squared_amplitudes.dat" using 1:2 lt 1 sm cs title "ground state", \
	"squared_amplitudes.dat" using 1:3 lt 2 sm cs title "2^+_1 state",  \
	"squared_amplitudes.dat" using 1:2 lt 1 pt 7  notitle "ground state", \
	"squared_amplitudes.dat" using 1:3 lt 2 pt 7  notitle "2^+_1 state",
set log y
set ylabel "d sigma / d b [mb/fm]"
plot \
	"diff_xsec.dat" using 1:8 lt 2 w l   title "2^+_1 state", \
	"diff_xsec.dat" using 1:8 lt 2 pt 7  notitle ,

set xlabel "theta [rad]"
set ylabel "sin(theta) (d sigma / d Omega) [mb/sr]"
plot \
	"diff_xsec.dat" using 2:9 lt 2 w l   title "2^+_1 state", \
	"diff_xsec.dat" using 2:9 lt 2 pt 7  notitle ,
unset multiplot 
EOF
gnuplot --persist visualize.gnuplot