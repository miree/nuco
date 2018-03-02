# calculate excitation cross section and visualize the b-dependent excitation probability
./nuco  \
	--Ap=80 --Zp=36 --At=197 --Zt=79 \
	--E=60 \
	--b=12 --bmax=200 --cross-section-integral-steps=100 \
	--method=relativistic --accuracy=1e-7 \
	--calc-cross-section \
	--levelE=0 --levelI=0 \
	--levelE=0.6166 --levelI=2 \
	--levelE=1.436  --levelI=4 \
	--MEfrom=0 --MEto=1 --MElambda=2 --MEvalue=61.8 \
	--MEfrom=1 --MEto=2 --MElambda=2 --MEvalue=113.6

cat << EOF > visualize.gnuplot
reset 
set xlabel "b [fm]"
set ylabel "amplitude"
set grid 
set term wxt size 600,1000
set key top right box opaque spacing 1.5 width 1
set multiplot layout 3,1
plot[0:] \
	"squared_amplitudes.dat" using 1:2 w l lt 1 lw 3 sm cs title "ground state", \
	"squared_amplitudes.dat" using 1:3 w l lt 2 lw 3 sm cs title "2^+_1 state",  \
	"squared_amplitudes.dat" using 1:4 w l lt 3 lw 3 sm cs title "4^+_1 state", 
set ylabel "d sigma / d b [mb/fm]"
unset key
set log y
plot \
	"diff_xsec.dat" using 1:6 w l lt 2 lw 3 title "2^+_1 state", \
	"diff_xsec.dat" using 1:8 w l lt 3 lw 3 title "4^+_1 state" ,

set xlabel "theta [fm]"
set ylabel "sin(theta) (d sigma / d Omega) [mb/sr]"
plot \
	"diff_xsec.dat" using 2:7 w l lt 2 lw 3 title "2^+_1 state", \
	"diff_xsec.dat" using 2:9 w l lt 3 lw 3 title "4^+_1 state"  ,
unset multiplot 
EOF
gnuplot --persist visualize.gnuplot