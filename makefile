GSL_FLAGS = `gsl-config --libs | sed 's/-/-L-/g'`
RDMDFLAGS = --build-only -O  -release  #--force --chatty

all: 
	rdmd $(RDMDFLAGS) $(GSL_FLAGS)  source/app.d && mv source/app nuco
	
