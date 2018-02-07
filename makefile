GSL_FLAGS = `gsl-config --libs | sed 's/-/-L-/g'`
RDMDFLAGS = --build-only --chatty --compiler=gdc #--force --chatty


SOURCES_NUCO = \
		  source/relat_dynamics.d      \
		  source/types.d               \
		  source/parameters.d          \
		  source/nonrelat_dynamics.d   \
		  source/app.d                 \
		  source/lebedev_quadrature.d  \
		  source/nucd/em.d             \
		  source/nucd/nucleus.d        \
		  source/nucd/geometry.d       \
		  source/nucd/kinematics.d     \
		  source/integrate.d           \
		  source/output_result.d      \


SOURCES_COULEX = \
		  source/relat_dynamics.d      \
		  source/types.d               \
		  source/parameters.d          \
		  source/nonrelat_dynamics.d   \
		  source/coulex.d              \
		  source/lebedev_quadrature.d  \
		  source/nucd/em.d             \
		  source/nucd/nucleus.d        \
		  source/nucd/geometry.d       \
		  source/nucd/kinematics.d     \
		  source/integrate.d           \

SOURCES_SLCOULEX = \
		  source/slcoulex.d            \
		  source/types.d               \
		  source/nucd/em.d             \
		  source/nucd/nucleus.d        \
		  source/nucd/geometry.d       \
		  source/nucd/kinematics.d     \

SOURCES_GSL = \
		  source/gsl/gsl_odeiv2.d      \
		  source/gsl/gsl_types.d       \
		  source/gsl/gsl_errno.d       \
		  source/gsl/gsl_spline.d      \
		  source/gsl/gsl_interp.d      \
		  source/gsl/gsl_sf_legendre.d \
		  source/gsl/gsl_sf_result.d   \
		  source/gsl/gsl_sf_bessel.d   \
		  source/gsl/gsl_mode.d        \
		  source/gsl/gsl_precision.d   \
		  source/gsl/gsl_sf_coupling.d \

SOURCES_TEST = \
		  source/relat_dynamics.d      \
		  source/types.d               \
		  source/parameters.d          \
		  source/nonrelat_dynamics.d   \
		  source/test_integration.d    \
		  source/lebedev_quadrature.d  \
		  source/nucd/em.d             \
		  source/nucd/nucleus.d        \
		  source/nucd/geometry.d       \
		  source/nucd/kinematics.d     \
		  source/integrate.d           \

all: nuco test_integration #slcoulex

nuco: $(SOURCES_NUCO)
	ldc2 $(SOURCES_NUCO) $(SOURCES_GSL) `gsl-config --libs | sed 's/-/-L-/g'` -O  -of=nuco

test_integration: $(SOURCES_TEST)	
	dmd $(SOURCES_TEST) $(SOURCES_GSL) `gsl-config --libs | sed 's/-/-L-/g'` -O -release -of=test_integration 

coulex: $(SOUCES_COULEX)
	dmd $(SOURCES_COULEX) $(SOURCES_GSL) `gsl-config --libs | sed 's/-/-L-/g'` -O -release -of=coulex 

slcoulex: source/slcoulex.d $(SOUCES_SLCOULEX)
	dmd $(SOURCES_SLCOULEX) $(SOURCES_GSL) `gsl-config --libs | sed 's/-/-L-/g'` -O -release -of=slcoulex 


# this does not work anymore ? :-/
#	rdmd $(RDMDFLAGS) $(GSL_FLAGS)  source/app.d && mv source/app nuco
	
