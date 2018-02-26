# calculate and visualize a trajectory
./nuco  \
	--Ap=80 --Zp=36 --At=197 --Zt=79 \
	--E=10 \
	--b=2 \
    --method=relativistic --accuracy=1e-7

cat << EOF > visualize.gnuplot
set xlabel "x [fm]"
set ylabel "y [fm]"
set grid 
set term wxt size 500,1000
set key top left box opaque spacing 1.5 width 1
set multiplot layout 2,1
set title "^{80}Kr on ^{197}Au at b = 2fm and E_{kin} = 10 MeV/u"
plot \
	"steps.dat" using 1:2 lt 1 w l title "projectile", \
	"steps.dat" using 3:4 lt 2 w l title "target", \
	"steps.dat" using 1:2 lt 1 pt 7 title "projectile", \
	"steps.dat" using 3:4 lt 2 pt 7 title "target"
unset title
plot[-40:40][-40:40] \
	"steps.dat" using 1:2 lt 1 w l title "projectile", \
	"steps.dat" using 3:4 lt 2 w l title "target", \
	"steps.dat" using 1:2 lt 1 pt 7 title "projectile", \
	"steps.dat" using 3:4 lt 2 pt 7 title "target"
unset multiplot 
EOF
gnuplot --persist visualize.gnuplot