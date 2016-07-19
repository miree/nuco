import types;
import parameters;

extern(C) int ode_classsical(double t, double* ys, double *dydts, void *params)
{
	// get parameters
	auto pars = cast(Parameters*)params;
	real  m1 = pars.Ap*u;
	real  m2 = pars.At*u;
	real  q1 = pars.Zp;
	real  q2 = pars.Zt;

	// map the ys[]-array to meaningful names
	auto x1 = Vec2([ys[ 0],ys[ 1]]);
	auto x2 = Vec2([ys[ 2],ys[ 3]]);
	auto v1 = Vec2([ys[ 4],ys[ 5]]);
	auto v2 = Vec2([ys[ 6],ys[ 7]]);
	
	// do the calculation for particle 2
	auto r12   = x2-x1;                           // vector from particle 1 to particle 2
	auto E12   = r12 * (ahc*q1 / r12.length^^3); // E-field created by particle 1 at position of particle 2
	auto dx2dt = v2;
	auto dv2dt = ( E12  ) * (q2/m2) * c^^2;

	// do the calculation for particle 1
	auto r21   = x1-x2;                           // vector from particle 2 to particle 1
	auto E21   = r21 * (ahc*q2 / r21.length^^3); // E-field created by particle 2 at position of particle 1
	auto dx1dt = v1;
	auto dv1dt = ( E21  ) * (q1/m1) * c^^2;

	
	// map the calculated quantities to output array
	dydts[0] = dx1dt[0];
	dydts[1] = dx1dt[1];
	dydts[2] = dx2dt[0];
	dydts[3] = dx2dt[1];
	
	dydts[4] = dv1dt[0];
	dydts[5] = dv1dt[1];
	dydts[6] = dv2dt[0];
	dydts[7] = dv2dt[1];
	
	pars.a1 = dv1dt; // These values need to be stored in the history. We use
	pars.a2 = dv2dt; // the parameter as a transport vehicle to the outside world
	
	return 0;
}

// same as classical, but with magnetic fields added
extern(C) int ode_magnetic(double t, double* ys, double *dydts, void *params)
{
	// get parameters
	auto pars = cast(Parameters*)params;
	real  m1 = pars.Ap*u;
	real  m2 = pars.At*u;
	real  q1 = pars.Zp;
	real  q2 = pars.Zt;

	// map the ys[]-array to meaningful names
	auto x1 = Vec2([ys[ 0],ys[ 1]]);
	auto x2 = Vec2([ys[ 2],ys[ 3]]);
	auto v1 = Vec2([ys[ 4],ys[ 5]]);
	auto v2 = Vec2([ys[ 6],ys[ 7]]);
	
	// do the calculation for particle 2
	auto r12   = x2-x1;                           // vector from particle 1 to particle 2
	auto E12   = r12 * (ahc*q1 / r12.length^^3); // E-field created by particle 1 at position of particle 2
	auto Bz12  = (v1[0]*r12[1] - v1[1]*r12[0])/c * 1.44*q1 / r12.length^^3; 
	auto dx2dt = v2;
	auto dv2dt = ( E12 + Vec2([v2[1]*Bz12,-v2[0]*Bz12])/c - v2*(v2*E12)/c^^2 ) * (q2/m2) * c^^2;

	// do the calculation for particle 1
	auto r21   = x1-x2;                           // vector from particle 2 to particle 1
	auto E21   = r21 * (ahc*q2 / r21.length^^3); // E-field created by particle 2 at position of particle 1
	auto Bz21  = (v2[0]*r21[1] - v2[1]*r21[0])/c * 1.44*q2 / r21.length^^3;
	auto dx1dt = v1;
	auto dv1dt = ( E21 + Vec2([v1[1]*Bz21,-v1[0]*Bz21])/c - v1*(v1*E21)/c^^2 ) * (q1/m1) * c^^2;

	
	// map the calculated quantities to output array
	dydts[0] = dx1dt[0];
	dydts[1] = dx1dt[1];
	dydts[2] = dx2dt[0];
	dydts[3] = dx2dt[1];
	
	dydts[4] = dv1dt[0];
	dydts[5] = dv1dt[1];
	dydts[6] = dv2dt[0];
	dydts[7] = dv2dt[1];
	
	pars.a1 = dv1dt; // These values need to be stored in the history. We use
	pars.a2 = dv2dt; // the parameter as a transport vehicle to the outside world
	
	return 0;
}
