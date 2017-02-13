GSL_FLAGS = `gsl-config --libs | sed 's/-/-L-/g'`
RDMDFLAGS = --build-only --chatty --compiler=gdc #--force --chatty

all:
	dmd `find -iname "*.d"` `gsl-config --libs | sed 's/-/-L-/g'` -O -release -of=nuco 

# this does not work anymore ? :-/
#	rdmd $(RDMDFLAGS) $(GSL_FLAGS)  source/app.d && mv source/app nuco
	
