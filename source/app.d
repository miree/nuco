import std.stdio;
import nucd.geometry;
import nucd.kinematics;
import gsl.gsl_odeiv2;
import gsl.gsl_errno;
import std.math;
import std.conv;
import std.getopt;
import std.complex;

// local imports
import types;
import nonrelat_dynamics;
import relat_dynamics;
import integrate;
import parameters;

// this is the same as get_retarded time in "relat_dynamics" but with a 3D test point
real get_retarded_time_3D(ref History h, double x, double y, double z, real t, real eps = 5e-8)
{
	// get retarded times, positions and velocity of particle 1
	real delta_t_ret_min = 1; // lower bound for t-t_ret
	real delta_t_ret_max = 0; // upper bound for t-t_ret
	for (;;)
	{
		auto x_ret = h.get(t-delta_t_ret_min).x;
		if (c*delta_t_ret_min > sqrt((x-x_ret[0])^^2 + (y-x_ret[1])^^2 + z^^2))
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
		if (c*(t-t_ret_med) < sqrt((x-x_ret[0])^^2 + (y-x_ret[1])^^2 + z^^2))
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

/// returns b2 as seen from an observer moved with b1. 
/// Arguments are in units of c
Vec2 velocity_addition(Vec2 b1, Vec2 b2)
{
	import nucd.geometry;
	Mat3 boost;
	boost.boost_direction(b1);
	//boost.print();
	
	auto g2 = gamma(b2.length);
	Vec3 bb2 = Vec3([g2, g2*b2[0], g2*b2[1]]);
	import std.stdio;
	bb2 = boost*bb2;
	return Vec2([bb2[1]/bb2[0],bb2[2]/bb2[0]]);
}

// transfrom an angle theta (measured in lab frame)
// into the CM frame. 
real transform_theta(real betaCM, real beta, real theta)
{
	auto beta_CM = Vec2([betaCM, 0]);
	auto v = Vec2([cos(theta),sin(theta)])*beta;
	auto v_CM = velocity_addition(beta_CM, v);
	return  atan2(v_CM[1],v_CM[0]);
}

Vec2 transform_direction(Vec2 r, Vec2 beta)
{
//	write(r, "  ");
	Vec2 e_beta = beta/beta.length;
	Vec2 r_para = e_beta*(r*e_beta);
	Vec2 r_perp = r-r_para;
	r_para *= gamma(beta.length);
	Vec2 retval = r_para + r_perp;
//	write("    ", retval);
	return retval;
	//import geometry;
	//Mat3 boost;
	//boost.boost_direction(beta);
	//auto rr = Vec3([0,r[0],r[1]]);
	//rr = boost*rr;
	//writeln("      ", rr);
	//return Vec2([rr[1],rr[2]]);
}

void help(string arg0)
{
	stderr.writeln("call for example: " ~ arg0 ~ " --Ap=85 --Zp=35 --At=197 --Zt=79 --E=10  --method=relativistic --accuracy=1e-7 --b=10");
}

void main(string[] args)
{
	// deal with input parameters from the command line
	Parameters params;
	getopt(args,
			"Ap",        &params.Ap,
			"Zp",        &params.Zp,
			"b",         &params.bp,
			"beta",      &params.betap,
			"E",         &params.Ep,
			"At",        &params.At,
			"Zt",        &params.Zt,
			"CM",        &params.CM,
			"rotate",    &params.rotate,
			"method",    &params.method,
			"accuracy",  &params.accuracy,
			"timeframe", &params.timeframe,
			"compare-rutherford", &params.compare_rutherford,
			"compare-rutherford-w", &params.compare_rutherford_w,
			"compare-rutherford-N", &params.compare_rutherford_N,
			"compare-SL-field", &params.compare_SL_field
			);
	
	if (isNaN(params.betap) && !isNaN(params.Ep))
	{
		params.betap = beta(gamma(params.Ep*params.Ap, params.Ap*u));
	}
	else if (!isNaN(params.betap) && isNaN(params.Ep))
	{
		params.Ep = kinEnergy(gamma(params.betap), params.Ap*u)/params.Ap;
	}
	else
	{
		stderr.writeln("ERROR: Either beta or E has to be defined. Defining none or both is impossible");
		help(args[0]);
		return;
	}
	
	// check if all needed parameters are set (i.e. not NaN)
	if (isNaN( params.Ap
	          *params.Zp 
	          *params.bp 
	          *params.betap 
	          *params.Ep 
	          *params.At 
		      *params.Zt))
	{
		stderr.writeln("ERROR: Not all parameters defined");
		help(args[0]);
		return;
	}
	
	if (params.CM == false && params.compare_rutherford == true)
	{
		stderr.writeln("ERROR: --compare-rutherford works only together with --CM");
		return;
	}
	
	// summarize all parameters for the user
	writeln("parameter summary:");
	writeln("   Ap    = ", params.Ap, " ");
	writeln("   Zp    = ", params.Zp, " ");
	writeln("   bp    = ", params.bp, " fm");
	writeln("   betap = ", params.betap, " ");
	writeln("   Ep    = ", params.Ep, " MeV/u");
	writeln("   At    = ", params.At, " ");
	writeln("   Zt    = ", params.Zt, " ");
	if (params.CM)
	{
		writeln("calculation is done in center-of-mass system");
	}
	else
	{
		writeln("calculation is done in laboratory system");
	}
	
	
	final switch(params.method)
	{
		case IntegrationMethod.classical:
			integrate.integrate(&ode_classsical, gsl_odeiv2_step_rkf45, params);
		break;
		case IntegrationMethod.magnetic:
			integrate.integrate(&ode_magnetic, gsl_odeiv2_step_rkf45, params);
		break;
		case IntegrationMethod.relativistic:
			integrate.integrate(&ode_relativistic, gsl_odeiv2_step_rkf45, params);
		break;
	}

	// calculate distance of closest approach
	real t0 = 0;
	real dt = 100;
	real dmin;
	do
	{
		real ti = t0+dt; // test time
		auto x1  = params.h1.get(t0).x;
		auto x2  = params.h2.get(t0).x;
		auto x1i = params.h1.get(ti).x;
		auto x2i = params.h2.get(ti).x;

		auto d  = (x2 - x1).length;
		auto di = (x2i-x1i).length;
		dmin = d;
		//writeln("distance ", d, " ", di, " ", dt);

		if (di < d)
		{
			t0 = ti;
		}
		else 
		{
			dt /= -2;
		}
	}
	while (dt*dt > 1e-12);
	//writeln("hallo");
	writeln("distance of closest approach: d_min=" , dmin , " fm    at   t_min=", t0 , " zs");
	
	if (params.compare_SL_field)
	{
		double tau01;
		double tau02;
		// determine the proper time of highest field strength
		{
			double[] taus;
			double[] Bs;
			foreach(i ; 0..params.h1.fields.length)
			{
				taus ~= params.h1.fields[i].tau;
				Bs   ~= params.h1.fields[i].B;
			}
			auto spl = Spline(taus, Bs);
			double tau0 = 0;
			double dtau = 100;
			do
			{
				auto taui = tau0+dtau;
				auto B  = spl.get(tau0);
				auto Bi = spl.get(taui);
				if (Bi*Bi > B*B) // find maximum
				{
					tau0 = taui;
				}
				else
				{
					dtau /= -2;
				}
			}
			while(dtau*dtau > 1e-12);
			tau01 = tau0;
			auto compareout  = File( "field-compare1.dat", "w+");
			foreach (field ; params.h1.fields)
			{
				real tau = field.tau - tau0;
				real q = params.Zt; // this is the field on the position of the projectile (caused by the target charge)
				real g = gamma(params.betap);
				real v = params.betap*c;
				real b = params.bp;
				real Ex =  q*g*v*tau  /(b^^2+(g*v*tau)^^2)^^(3./2.);
				real Ey =  q*g*b      /(b^^2+(g*v*tau)^^2)^^(3./2.);
				real B  = -q*g*b*(v/c)/(b^^2+(g*v*tau)^^2)^^(3./2.);
				compareout.writefln("%.10f %.10f %.10f %.10f %.10f %.10f %.10f", tau, field.E[0], field.E[1], field.B, Ex, Ey, B);
			}
			
		}
		{
			double[] taus;
			double[] Bs;
			foreach(i ; 0..params.h2.fields.length)
			{
				taus ~= params.h2.fields[i].tau;
				Bs   ~= params.h2.fields[i].B;
			}
			auto spl = Spline(taus, Bs);
			double tau0 = 0;
			double dtau = 100;
			do
			{
				auto taui = tau0+dtau;
				auto B  = spl.get(tau0);
				auto Bi = spl.get(taui);
				if (Bi*Bi > B*B) // find maximum
				{
					tau0 = taui;
				}
				else
				{
					dtau /= -2;
				}
			}
			while(dtau*dtau > 1e-12);
			tau02 = tau0;
		}
		writeln("tau01 = ", tau01);
		writeln("tau02 = ", tau02);
	}
	
	if (params.compare_rutherford)
	{
		auto compareout  = File( "rutherford-compare.dat", "w+");
		assert(params.h1.points.length == params.h1.points.length);
	
		Matrix!(2,2,0,real) rot; 
		rot.identity();
		double 	angle = atan2(params.h1.points[$-1].v[1], params.h1.points[$-1].v[0]);
		rot.rotation(0,1,(PI-angle)/2);
		
		real mu = u * params.Ap*params.At/(params.Ap+params.At);
		real a  = params.Zp*params.Zt*ahc/mu/params.betap^^2;
		real theta_half = atan2(a,cast(real)params.bp);
		foreach(n; 0..params.compare_rutherford_N)
		{
			real e = 1/sin(theta_half);
			real w = params.compare_rutherford_w*(2.0*n/params.compare_rutherford_N - 1.0);
			real xx = a * (cosh(w) + e);
			real yy = a * sqrt(e^^2-1)*sinh(w);
			real tt = a * (e * sinh(w) + w) / params.betap / c;
			
			auto pp1 = params.h1.get(tt+t0).x;
			auto pp2 = params.h2.get(tt+t0).x;
			auto p1 = rot*pp1;
			auto p2 = rot*pp2;

			compareout.writefln("%.10f %.10f %.10f    %.10f        %.10f", 
			                     tt,   -xx,   yy,   p1[0]-p2[0], p1[1]-p2[1]);
		}
	}


	auto stepout  = File( "steps.dat", "w+");
	assert(params.h1.points.length == params.h2.points.length);
	assert(params.h1.points_partner_ret.length == params.h2.points_partner_ret.length);
	assert(params.h1.points.length == params.h1.points_partner_ret.length);
	foreach(i;0..params.h1.points.length)
	{
		Matrix!(2,2,0,real) rot; 
		rot.identity();
		if (params.rotate)
		{
			double 	angle = atan2(params.h1.points[$-1].v[1], params.h1.points[$-1].v[0]);
			rot.rotation(0,1,(PI-angle)/2);
		}
		
		auto pp1 = params.h1.points[i].x;
		auto pp2 = params.h2.points[i].x;
		auto p1 = rot*(pp1);
		auto p2 = rot*pp2;
		auto rp1 = params.h1.points_partner_ret[i].x;
		auto rp2 = params.h2.points_partner_ret[i].x;
		auto r1 = rot*(rp1);
		auto r2 = rot*rp2;
		stepout.writefln("%.10f  %.10f  %.10f  %.10f    %.10f %.10f %.10f %.10f", 
					      p1[0], p1[1], p2[0], p2[1],
					      r1[0], r1[1], r2[0], r2[1]);
	}
	auto fieldout  = File( "field.dat", "w+");
	assert(params.h1.fields.length == params.h1.fields.length);
	foreach(i;0..params.h1.points.length)
	{
		auto E1   = params.h1.fields[i].E.length;
		auto B1   = params.h1.fields[i].B;
		auto phi1 = params.h1.fields[i].pot[0];
		auto A1   = Vec2([params.h1.fields[i].pot[1],params.h1.fields[i].pot[2]]).length;
		auto tau1 = params.h1.fields[i].tau;
		
		auto E2   = params.h2.fields[i].E.length;
		auto B2   = params.h2.fields[i].B;
		auto phi2 = params.h2.fields[i].pot[0];
		auto A2   = Vec2([params.h2.fields[i].pot[1],params.h2.fields[i].pot[2]]).length;
		auto tau2 = params.h2.fields[i].tau;
		//stepout.writefln("%.10f %.10f", p1[0], p1[1]);
		fieldout.writefln("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f", tau1, E1, B1, phi1, A1, tau2, E2, B2, phi2, A2);
	}

	double[9] theta;

	if (params.CM)
	{
		writeln("-------------------------CM system--------------------------------");
		writeln("beta_proj_in   = ", params.h1.points[0  ].v[0]/c, " ", params.h1.points[0].v[1]/c,   "   ", params.h1.points[0].v.length/c);
		writeln("E_proj_in      = ", kinEnergy(gamma(params.h1.points[0  ].v[0]/c), params.Ap*u)/params.Ap);
		writeln("beta_proj_out  = ", params.h1.points[$-1].v[0]/c, " ", params.h1.points[$-1].v[1]/c, "   ", params.h1.points[$-1].v.length/c);
		auto beta_rel = addVelocity(params.h1.points[$-1].v/c, params.h2.points[$-1].v*(-1/c));
		writeln("theta_proj_out = ", atan2(beta_rel[1],beta_rel[0]));
		theta[0] = atan2(beta_rel[1],beta_rel[0]);
		writeln("beta_targ_in   = ", params.h2.points[0  ].v[0]/c, " ", params.h2.points[0].v[1]/c,   "   ", params.h2.points[0].v.length/c);
		writeln("E_targ_in      = ", kinEnergy(gamma(params.h2.points[0  ].v[0]/c), params.At*u)/params.At);
		writeln("beta_targ_out  = ", params.h2.points[$-1].v[0]/c, " ", params.h2.points[$-1].v[1]/c, "   ", params.h2.points[$-1].v.length/c);
		
		
		// transfrom the velocities into the lab system:
		auto beta_CM = Vec2([-params.beta_CM, 0]);
		//writeln(beta_CM);
		//writeln(params.h1.points[0].v/c);
		auto beta_proj_in_lab  = velocity_addition(beta_CM, params.h1.points[0  ].v/c);
		auto beta_proj_out_lab = velocity_addition(beta_CM, params.h1.points[$-1].v/c);
		auto beta_targ_in_lab  = velocity_addition(beta_CM, params.h2.points[0  ].v/c);
		auto beta_targ_out_lab = velocity_addition(beta_CM, params.h2.points[$-1].v/c);
		
		writeln("-------------------------lab system--------------------------------");
		writeln("beta_proj_in   = ", beta_proj_in_lab[0] , " ", beta_proj_in_lab[1] , "   ", beta_proj_in_lab.length);
		writeln("beta_proj_out  = ", beta_proj_out_lab[0], " ", beta_proj_out_lab[1], "   ", beta_proj_out_lab.length);
		writeln("E_proj_in      = ", kinEnergy( gamma(beta_proj_in_lab.length), u), " MeV/u");
		writeln("E_proj_out     = ", kinEnergy( gamma(beta_proj_out_lab.length), u), " MeV/u");
		beta_rel = addVelocity(beta_proj_out_lab, beta_targ_out_lab*(-1));
		writeln("theta_proj_out (relative velocity) = ", atan2(beta_rel[1],beta_rel[0]));
		writeln("theta_proj_out (lab velocity)      = ", atan2(beta_proj_out_lab[1],beta_proj_out_lab[0]));
		theta[1] = atan2(beta_rel[1],beta_rel[0]);
		theta[2] = atan2(beta_proj_out_lab[1],beta_proj_out_lab[0]);
		writeln("beta_targ_in   = ", beta_targ_in_lab[0] , " ", beta_targ_in_lab[1] , "   ", beta_targ_in_lab.length);
		writeln("beta_targ_out  = ", beta_targ_out_lab[0], " ", beta_targ_out_lab[1], "   ", beta_targ_out_lab.length);
		writeln("E_targ_in      = ", kinEnergy( gamma(beta_targ_in_lab.length), u), " MeV/u");
		writeln("E_targ_out     = ", kinEnergy( gamma(beta_targ_out_lab.length), u), " MeV/u");
	}
	else
	{
		writeln("-------------------------lab system--------------------------------");
		writeln("beta_proj_in   = ", params.h1.points[0  ].v[0]/c, " ", params.h1.points[0].v[1]/c,   "   ", params.h1.points[0].v.length/c);
		writeln("beta_proj_out  = ", params.h1.points[$-1].v[0]/c, " ", params.h1.points[$-1].v[1]/c, "   ", params.h1.points[$-1].v.length/c);
		writeln("E_proj_in      = ", kinEnergy( gamma(params.h1.points[0].v.length/c), u), " MeV/u");
		writeln("E_proj_out     = ", kinEnergy( gamma(params.h1.points[$-1].v.length/c), u), " MeV/u");
		auto beta_rel = addVelocity(params.h1.points[$-1].v/c, params.h2.points[$-1].v*(-1/c));
		writeln("theta_proj_out (relative velocity) = ", atan2(beta_rel[1],beta_rel[0]));
		writeln("theta_proj_out (lab velocity)      = ", atan2(params.h1.points[$-1].v[1],params.h1.points[$-1].v[0]));
		theta[1] = atan2(beta_rel[1],beta_rel[0]);
		theta[2] = atan2(params.h1.points[$-1].v[1],params.h1.points[$-1].v[0]);
		writeln("beta_targ_in   = ", params.h2.points[0  ].v[0]/c, " ", params.h2.points[0].v[1]/c,   "   ", params.h2.points[0].v.length/c);
		writeln("beta_targ_out  = ", params.h2.points[$-1].v[0]/c, " ", params.h2.points[$-1].v[1]/c, "   ", params.h2.points[$-1].v.length/c);
		writeln("E_targ_in      = ", kinEnergy( gamma(params.h2.points[0].v.length/c), u), " MeV/u");
		writeln("E_targ_out     = ", kinEnergy( gamma(params.h2.points[$-1].v.length/c), u), " MeV/u");
		
		// transfrom the velocities into the CM system:
		auto beta_CM = Vec2([params.beta_CM, 0]);
		//writeln(beta_CM);
		//writeln(params.h1.points[0].v/c);
		auto beta_proj_in_cm  = velocity_addition(beta_CM, params.h1.points[0  ].v/c);
		auto beta_proj_out_cm = velocity_addition(beta_CM, params.h1.points[$-1].v/c);
		auto beta_targ_in_cm  = velocity_addition(beta_CM, params.h2.points[0  ].v/c);
		auto beta_targ_out_cm = velocity_addition(beta_CM, params.h2.points[$-1].v/c);
		
		writeln("-------------------------CM system--------------------------------");
		writeln("beta_proj_in   = ", beta_proj_in_cm[0] , " ", beta_proj_in_cm[1] , "   ", beta_proj_in_cm.length);
		writeln("beta_proj_out  = ", beta_proj_out_cm[0], " ", beta_proj_out_cm[1], "   ", beta_proj_out_cm.length);
		beta_rel = addVelocity(beta_proj_out_cm, beta_targ_out_cm*(-1));
		writeln("theta_proj_out = ", atan2(beta_rel[1],beta_rel[0]));
		theta[0] = atan2(beta_rel[1],beta_rel[0]);
		writeln("beta_targ_in   = ", beta_targ_in_cm[0] , " ", beta_targ_in_cm[1] , "   ", beta_targ_in_cm.length);
		writeln("beta_targ_out  = ", beta_targ_out_cm[0], " ", beta_targ_out_cm[1], "   ", beta_targ_out_cm.length);
	}

		
	writeln("-------------------------approximative results----------------------");
	{
		real modified_impact_parameter = params.bp + (PI/2) * params.Zt*params.Zp*ahc / (params.Ap*u*params.betap^^2*gamma(params.betap));
		writeln("impact parameter: ", params.bp, "   modified impact parameter: ", modified_impact_parameter);
		real theta_straight_line = 2*params.Zt*params.Zp*ahc/(modified_impact_parameter*params.Ep*params.Ap);
		
		//real theta_straight_line = 2*params.Zt*params.Zp*ahc/(params.bp*params.Ep*params.Ap);
		writeln("theta straight-line approximation (compare with theta in labframe with lab velocity) = ", theta_straight_line);
		theta[3] = theta_straight_line;
		auto theta_straight_line_CM = transform_theta(params.beta_CM, params.betap, theta_straight_line);
		theta[4] = theta_straight_line_CM;
		writeln("theta straight-line approximation CM (compare with theta in CM frame) = ", theta_straight_line_CM);
	}
	{
		real mu = u * params.Ap*params.At/(params.Ap+params.At);
		real a  = params.Zp*params.Zt*ahc/mu/params.betap^^2;
		real theta_half = atan2(a,cast(real)params.bp);
		writeln("theta classical (compare with theta in CM system) = ", theta_half*2);
		theta[5] = theta_half*2;

		auto beta_CM = Vec2([-params.beta_CM, 0]);
		auto v = Vec2([cos(2*theta_half),sin(2*theta_half)])*(params.betap-params.beta_CM);
		auto v_lab = velocity_addition(beta_CM, v);
		auto theta_lab = atan2(v_lab[1],v_lab[0]);
		theta[6] = theta_lab;
		writeln("theta classical (compare with theta in lab system with relative velocity?) = ", theta_lab);
	}
	{
		real mu = u * params.Ap*params.At/(params.Ap+params.At);
		real a  = params.Zp*params.Zt*ahc/mu/params.betap^^2/gamma(params.betap);
		real theta_half = atan2(a,cast(real)params.bp);
		writeln("theta classical (Bertulani's modification, compare with theta in CM system) = ", theta_half*2);
		theta[7] = theta_half*2;

		auto beta_CM = Vec2([-params.beta_CM, 0]);
		auto v = Vec2([cos(2*theta_half),sin(2*theta_half)])*(params.betap-params.beta_CM);
		auto v_lab = velocity_addition(beta_CM, v);
		auto theta_lab = atan2(v_lab[1],v_lab[0]);
		theta[8] = theta_lab;
		writeln("theta classical (Bertulani's modification, compare with theta in lab system relative velocity) = ", theta_lab);
	}


	auto theta_file = File("theta_file.dat", "a+");
	theta_file.write(params.Ep, "   ");
	theta_file.write(params.bp, "   ");
	foreach(t;theta)
	{
		theta_file.writef("%20.20f \t ", t);
	}
	theta_file.writeln();
	
	
//	writeln("multipoles at 1st");
//	auto multipole_file = File("multipole_at_1.dat", "w+");
//	for (double t = -0.5; t <= 0.5; t += 0.0002)
//	{
//		//writeln("multipole moment calculation at t = ", t);
//		auto center = params.h1.get(t);  // the position of particle 1 at time t;
//		//writeln(center.v);
//		import lebedev_quadrature;
//		Vec3 pot;
//		import std.complex;
//		Complex!double[int] mE1;
//		foreach(m;-1..2) mE1[m] = complex(0,0);
//		//auto lq = LQ(0,0,0,0,0,1);
//		foreach(lq;lq0110)
//		{
//			real r = 1;
//			auto dr = transform_direction(Vec2([lq.x,lq.y])*r, center.v/c);
//			//multipole_file.writeln(dr[0]," ", dr[1]);
//			auto center3D = Vec3([center.x[0]+dr[0], center.x[1]+dr[1],lq.z*r]);
//			auto t_ret = get_retarded_time_3D(params.h2, center3D[0], center3D[1], center3D[2], t);
//			auto point_ret = params.h2.get(t_ret);
//			auto retx3D   = Vec3([point_ret.x[0], point_ret.x[1], 0]);
//			auto retv3D   = Vec3([point_ret.v[0], point_ret.v[1], 0]); 
//			pot = transform_potentials(center.v/c, 
//						potential3D(center3D - retx3D, retv3D/c, params.Zt));
//			foreach(m;-1..2)
//			{
//				import nucd.em;
//				mE1[m] += pot[0] * lq.w ;//* Ylm(1, m, lq.theta, lq.phi);
//			}	
//		}
//		auto t_ret = get_retarded_time(params.h2, center.x, t);
//		auto point_ret = params.h2.get(t_ret);
//		pot = transform_potentials(center.v/c, 
//						potential(center.x - point_ret.x, point_ret.v/c, params.Zp));
//		multipole_file.writeln(center.tau, " ", pot[0], " ", mE1[-1].abs, " ", mE1[0].abs, " ", center.v.length );
//	}
//
//	writeln("multipoles at 2nd");
//	auto multipole_file2 = File("multipole_at_2.dat", "w+");
//	for (double t = -0.5; t <= 0.5; t += 0.0002)
//	{
//		//writeln("multipole moment calculation at t = ", t);
//		auto center = params.h2.get(t);  // the position of particle 1 at time t;
//		//writeln(center.v);
//		import lebedev_quadrature;
//		Vec3 pot;
//		import std.complex;
//		Complex!double[int] mE1;
//		foreach(m;-1..2) mE1[m] = complex(0,0);
//		//auto lq = LQ(0,0,0,0,0,1);
//		foreach(lq;lq0110)
//		{
//			real r = 1;
//			auto dr = transform_direction(Vec2([lq.x,lq.y])*r, center.v/c);
//			//multipole_file2.writeln(dr[0]," ", dr[1]);
//			auto center3D = Vec3([center.x[0]+dr[0], center.x[1]+dr[1],lq.z*r]);
//			auto t_ret = get_retarded_time_3D(params.h1, center3D[0], center3D[1], center3D[2], t);
//			auto point_ret = params.h1.get(t_ret);
//			auto retx3D   = Vec3([point_ret.x[0], point_ret.x[1], 0]);
//			auto retv3D   = Vec3([point_ret.v[0], point_ret.v[1], 0]); 
//			pot = transform_potentials(center.v/c, 
//						potential3D(center3D - retx3D, retv3D/c, params.Zp));
//			foreach(m;-1..2)
//			{
//				import nucd.em;
//				mE1[m] += pot[0] * lq.w ;//* Ylm(1, m, lq.theta, lq.phi);
//			}	
//		}
//		auto t_ret = get_retarded_time(params.h1, center.x, t);
//		auto point_ret = params.h1.get(t_ret);
//		pot = transform_potentials(center.v/c, 
//						potential(center.x - point_ret.x, point_ret.v/c, params.Zt));
//						
//		multipole_file2.writeln(center.tau, " ", pot[0], " ", mE1[-1].abs, " ", mE1[0].abs, " ", center.v.length);
//	}

}