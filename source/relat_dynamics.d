import types;
import parameters;
import nucd.kinematics;


void transform_fields(Vec2 beta, ref Vec2 E, ref real B)
{
	auto B_old = B;
	auto E_old = E;
	auto g = gamma(beta.length);
	B = g*(B_old - (beta[0]*E_old[1] - beta[1]*E_old[0]));
	
	auto e_beta = beta/beta.length;
	auto E_para = e_beta*(E_old*e_beta); // Efield parallel to velocity
	auto E_perp = E_old - E_para;        // Efield perpendiculear to velocity
	
	E_perp = (E_perp + Vec2([beta[1]*B_old, 
                            -beta[0]*B_old]))*g;
	E = E_para + E_perp;
	import std.stdio;
	//writeln(E_old.length, "   ", (E_para+E_perp).length);
	//writeln(E_old, " , ", B_old, "  -->  ", E, " , ", B);
}

auto transform_potentials(Vec2 beta, Vec3 pot)
{
	import nucd.geometry;
	Mat3 boost;
	boost.boost_direction(beta);
	return boost*pot;
}

real get_retarded_time(ref History h, ref Vec2 x, real t, real eps = 5e-8)
{
	// get retarded times, positions and velocity of particle 1
	real delta_t_ret_min = 1; // lower bound for t-t_ret
	real delta_t_ret_max = 0; // upper bound for t-t_ret
	for (;;)
	{
		auto x_ret = h.get(t-delta_t_ret_min).x;
		if (c*delta_t_ret_min > (x-x_ret).length)
			break;
		delta_t_ret_max = delta_t_ret_min;
		delta_t_ret_min *= 2.;
	}
	real t_ret_min = t-delta_t_ret_min;
	real t_ret_max = t-delta_t_ret_max;
	// Now we know that t_ret is inside the interval [t_ret_min, t_ret_max].
	// Look for t_ret by means of bisecting that interval
	while (t_ret_max-t_ret_min > eps)
	{
		real t_ret_med = 0.5*(t_ret_min+t_ret_max);
		auto x_ret = h.get(t_ret_med).x;
		if (c*(t-t_ret_med) < (x-x_ret).length)
		{
			t_ret_max = t_ret_med;
		}
		else
		{
			t_ret_min = t_ret_med;
		}
	}
	return 0.5*(t_ret_min+t_ret_max);	
}

// computes the potential 4-vector (which is a 3-vector because 
// only 2 spatial dimensions are relevant)
auto potential(Vec2 r_ret, Vec2 beta_ret, real q_ret)
{
	auto R   = r_ret.length;
	auto phi = q_ret/(R - r_ret*beta_ret);
	auto A   = beta_ret * q_ret / ( R - r_ret*beta_ret);
	return Vec3([phi,A[0],A[1]]); 
}
auto potential3D(Vec3 r_ret, Vec3 beta_ret, real q_ret)
{
	auto R   = r_ret.length;
	auto phi = q_ret/(R - r_ret*beta_ret);
	auto A   = beta_ret * q_ret / ( R - r_ret*beta_ret);
	return Vec3([phi,A[0],A[1]]); 
}
//auto phi_potential( Vec2 r_ret, Vec2 beta_ret, real q_ret)
//{
//	auto R  = r_ret.length;
//	return q_ret/(R - r_ret*beta_ret);
//}
//auto A_potential( Vec2 r_ret, Vec2 beta_ret, real q_ret)
//{
//	auto R  = r_ret.length;
//	return beta_ret * q_ret / ( R - r_ret*beta_ret);
//}

auto E_field_uniform( Vec2 r_ret,  Vec2 beta_ret, real q_ret)
{
	auto R  = r_ret.length;
	auto er = r_ret/R;
	auto g  = gamma(beta_ret.length);
	return (er-beta_ret) * q_ret / ( (R*g)^^2 * (1 - beta_ret*er)^^3);
}
auto E_field_radiative( Vec2 r_ret,  Vec2 beta_ret,  Vec2 dbetadt_ret, real q_ret)
{
	auto R  = r_ret.length;
	auto er = r_ret/R;
	auto g  = gamma(beta_ret.length);
	return ((er-beta_ret)*(er*dbetadt_ret) - dbetadt_ret*(er*(er-beta_ret))) * q_ret 
	    /  (c * R * (1 - beta_ret*er)^^3);
}
auto B_field(Vec2 x,  Vec2 x_ret,  Vec2 E)
{
	auto r  = x-x_ret;
	auto R  = r.length;
	auto er = r/R;
	return (er[0]*E[1] - er[1]*E[0]);
}

auto acceleration(Vec2 E, real Bz, Vec2 v, real q, real m)
{
	auto g = gamma(v.length/c);
	return (E + Vec2([v[1]*Bz,-v[0]*Bz])/c - v*(v*E)/c^^2) * ahc*q/(m*g) * c^^2;
}

extern(C) int ode_relativistic(double t, double* ys, double *dydts, void *params)
{
	// get parameters
	auto pars = cast(Parameters*)params;
	real  m1 = pars.Ap*u;
	real  m2 = pars.At*u;
	real  q1 = pars.Zp;
	real  q2 = pars.Zt;

	// map the ys[]-array to meaningful names
	auto x1 = Vec2([ys[0],ys[1]]);
	auto x2 = Vec2([ys[2],ys[3]]);
	auto v1 = Vec2([ys[4],ys[5]]);
	auto v2 = Vec2([ys[6],ys[7]]);
	auto tau1 = ys[8];
	auto tau2 = ys[9];

	if (pars.h1.points.length == 0)
	{
		assert(pars.h1.points.length == pars.h2.points.length);
		pars.h1.add(t, x1, v1, pars.a1, tau1, Vec2([0,0]), 0, Vec3([0,0,0]));
		pars.h2.add(t, x2, v2, pars.a2, tau2, Vec2([0,0]), 0, Vec3([0,0,0]));

		// get retarded time, position and velocity of particle 2
		real t_ret2 = get_retarded_time(pars.h2, x1, t);
		pars.p2_ret = pars.h2.get(t_ret2);
		// get retarded time, position and velocity of particle 1
		real t_ret1 = get_retarded_time(pars.h1, x2, t);
		pars.p1_ret = pars.h1.get(t_ret1);

		pars.h1.points_partner_ret ~= pars.p2_ret;
		pars.h2.points_partner_ret ~= pars.p1_ret;
	}

	// get retarded time, position and velocity of particle 2
	real t_ret2 = get_retarded_time(pars.h2, x1, t);
	pars.p2_ret = pars.h2.get(t_ret2);
	// compute dynamics of particle 1 in EM-field created by particle 2
	auto E21_u  = E_field_uniform  (x1-pars.p2_ret.x, pars.p2_ret.v/c, q2);
	auto E21_r  = E_field_radiative(x1-pars.p2_ret.x, pars.p2_ret.v/c, pars.p2_ret.a/c, q2);
	auto E21    = E21_u + E21_r;
	auto B21    = B_field(x1, pars.p2_ret.x, E21);
	auto dx1dt = v1;
	auto dv1dt = acceleration(E21, B21, v1, q1, m1);
	auto dtau1dt = 1./gamma(v1.length/c);

	// get retarded time, position and velocity of particle 1
	real t_ret1 = get_retarded_time(pars.h1, x2, t);
	pars.p1_ret = pars.h1.get(t_ret1);
	// compute dynamics of particle 2 in EM-field created by particle 1
	auto E12_u  = E_field_uniform  (x2-pars.p1_ret.x, pars.p1_ret.v/c, q1);
	auto E12_r  = E_field_radiative(x2-pars.p1_ret.x, pars.p1_ret.v/c, pars.p1_ret.a/c, q1);
	auto E12    = E12_u + E12_r;
	auto B12    = B_field(x2, pars.p1_ret.x, E12);
	auto dx2dt  = v2;
	auto dv2dt  = acceleration(E12, B12, v2, q2, m2);
	auto dtau2dt = 1./gamma(v2.length/c);

	// map the calculated quantities to output array
	dydts[0] = dx1dt[0];
	dydts[1] = dx1dt[1];
	dydts[2] = dx2dt[0];
	dydts[3] = dx2dt[1];

	dydts[4] = dv1dt[0];
	dydts[5] = dv1dt[1];
	dydts[6] = dv2dt[0];
	dydts[7] = dv2dt[1];
	
	dydts[8] = dtau1dt;
	dydts[9] = dtau2dt;

	pars.a1 = dv1dt; // These values need to be stored in the history. We use
	pars.a2 = dv2dt; // the parameter as a transport vehicle to the outside world
	
	transform_fields(v1/c, E21, B21); // transform into rest frame of nucleus 1
	transform_fields(v2/c, E12, B12); // transform into rest frame of nucleus 2
	
	pars.E21 = E21;  // save the local E- and B-fields and potentials for both particles
	pars.B21 = B21;
	pars.E12 = E12;
	pars.B12 = B12;
	
	pars.pot21 = transform_potentials(v1/c ,potential(x1-pars.p2_ret.x, pars.p2_ret.v/c, q2));
	pars.pot12 = transform_potentials(v2/c ,potential(x2-pars.p1_ret.x, pars.p1_ret.v/c, q1));

	return 0;
}
