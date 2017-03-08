/* 
 * Copyright (C) 2015,2016 Michael Reese
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */



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
import nucd.kinematics;



void help(string arg0)
{
	stderr.writeln("call for example: " ~ arg0 ~ " --Ap=85 --Zp=35 --At=197 --Zt=79 --E=10  --method=relativistic --accuracy=1e-7 --b=10");
}

void main(string[] args)
{
	// deal with input parameters from the command line
	Parameters params;

	// hard-code nuclear structure information



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

	writeln("distance of closest approach: d_min=" , params.dmin , " fm    at   t_min=", params.t0 , " zs");
	
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
			tau01 = params.p1_tau0;
			auto compareout  = File( "field-compare1.dat", "w+");
			foreach (field ; params.h1.fields)
			{
				real tau = field.tau - tau01;
				real q = params.Zt; // this is the field on the position of the projectile (caused by the target charge)
				real g = gamma(params.betap);
				real v = params.betap*c;
				real mu = u * params.Ap*params.At/(params.Ap+params.At);
				real b = params.bp + (PI/2) * params.Zt*params.Zp*ahc / (mu*params.betap^^2*gamma(params.betap));
				//real b = params.bp;
				real Ex =  q*g*v*tau  /(b^^2+(g*v*tau)^^2)^^(3./2.);
				real Ey =  q*g*b      /(b^^2+(g*v*tau)^^2)^^(3./2.);
				real B  = -q*g*b*(v/c)/(b^^2+(g*v*tau)^^2)^^(3./2.);
				compareout.writefln("%.10f %.10f %.10f %.10f %.10f %.10f %.10f", tau, field.E[0], field.E[1], field.B, Ex, Ey, B);
			}
		}
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
			
			auto pp1 = params.h1.get(tt+params.t0).x;
			auto pp2 = params.h2.get(tt+params.t0).x;
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
		writeln("beta_CM       = ", params.beta_CM);
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
		real mu = u * params.Ap*params.At/(params.Ap+params.At);
		real modified_impact_parameter = params.bp + (PI/2) * params.Zt*params.Zp*ahc / (mu*params.betap^^2*gamma(params.betap));
		writeln("impact parameter: ", params.bp, "   modified impact parameter: ", modified_impact_parameter);
		real theta_straight_line = 2*params.Zt*params.Zp*ahc/(params.Ap*u * params.betap^^2 * gamma(params.betap) * modified_impact_parameter);
		
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
		writeln("theta classical (infinite mass limit, compare with theta in CM system) = ", theta_half*2);
		theta[7] = theta_half*2;

		auto beta_CM = Vec2([-params.beta_CM, 0]);
		auto v = Vec2([cos(2*theta_half),sin(2*theta_half)])*(params.betap-params.beta_CM);
		auto v_lab = velocity_addition(beta_CM, v);
		auto theta_lab = atan2(v_lab[1],v_lab[0]);
		theta[8] = theta_lab;
		writeln("theta classical (infinite mass limit, compare with theta in lab system relative velocity) = ", theta_lab);
	}


	auto theta_file = File("theta_file.dat", "a+");
	theta_file.write(params.Ep, "   ");
	theta_file.write(params.bp, "   ");
	foreach(t;theta)
	{
		theta_file.writef("%20.20f \t ", t);
	}
	theta_file.writeln();


	writeln("multipols at projectile position");
	auto infofile = File("infofile.dat", "w+");	
	foreach(p; params.h1.points)
	{
		real t = p.t;//params.t0;
		auto projectile_center = params.h1.get(t);
		auto projectile_pos_2d = projectile_center.x;
		auto t_ret             = get_retarded_time(params.h2, projectile_pos_2d, t);
		auto target_center     = params.h2.get(t_ret);
		real r  = 0.001; // radius of the sphere
		import lebedev_quadrature;
		int l = 3;
		Complex!double[int] mEl;
		foreach(m;-l..l+1) mEl[m] = complex(0,0);

		////writeln((target_center.x - projectile_center.x).length, " ", c*(t_ret-t));
		auto R    = direction_from_mooved_system(projectile_center.v/c,  projectile_center.x - target_center.x);
		auto beta = velocity_addition(projectile_center.v/c, target_center.v/c);
		auto potential = potential(R, beta, params.Zt);


		//foreach(lq;lq0110)
		//{
		//	import types;
		//	auto dr = sim_frame_vector(lq.x,lq.y,lq.z)*r;//transform_direction(sim_frame_vector(lq.x,lq.y,lq.z)*r, Vec3(projectile_center.v/c));
		//	//infofile.writeln(dr[0]," ",dr[1]," ",dr[2]);
		//	auto R = projectile_center.x-target_center.x;
		//	auto R3d = Vec3([R[0],R[1],0])+dr;
		//	auto potential = transform_potentials(projectile_center.v/c,
		//										  potential3D(R3d, Vec3(target_center.v/c), params.Zt));
		//	foreach(m;-l..l+1)
		//	{
		//		import nucd.em;
		//		mEl[m] += potential[0] * lq.w * Ylm(l, m, lq.theta, lq.phi);
		//	}		
		//}
		//foreach(m;-l..l+1) mEl[m] *= 4*PI/r^^l;
		//writeln("-----");
		infofile.write(t, " ", p.tau, " ");//, p.fields.E[0], " " , p.fields.E[1], " " , p.fields.B, " ");
		infofile.writeln(potential[0], " ", potential[1], " ", potential[2], " ", R.length);
		//foreach(m;-l..l+1)
		//{
		//	//infofile.write(mEl[m].abs, " " );
		//	//writeln(m,":", (mEl[m].re!=0)?(((-1)^^m)*mEl[m].re/mEl[-m].re):0, " ", (mEl[m].im!=0)?(((-1)^^m)*mEl[m].im/mEl[-m].im):0 );
		//}		
		//infofile.writeln();

	}
	foreach(field;params.h1.fields)
	{
		//infofile.writeln(field.tau, " ", field.pot[0], " ", field.pot[1], " ", field.pot[2]);
	}

	//writeln("multipoles at 1st");
	//foreach(p; params.h1.points)
	//{
	//	double t = p.t;
	//	//writeln("multipole moment calculation at t = ", t);
	//	auto center = params.h1.get(t);  // the position of particle 1 at time t;
	//	//writeln(center.v);
	//	import lebedev_quadrature;
	//	Vec3 pot;
	//	import std.complex;
	//	Complex!double[int] mE1;
	//	foreach(m;-1..2) mE1[m] = complex(0,0);
	//	//auto lq = LQ(0,0,0,0,0,1);
	//	foreach(lq;lq0110)
	//	{
	//		real r = 1;
	//		auto dr = transform_direction(Vec2([lq.x,lq.y])*r, center.v/c);
	//		//multipole_file.writeln(dr[0]," ", dr[1]);
	//		auto center3D = Vec3([center.x[0]+dr[0], center.x[1]+dr[1],lq.z*r]);
	//		auto t_ret = get_retarded_time(params.h2, center3D, t);
	//		auto point_ret = params.h2.get(t_ret);
	//		auto retx3D   = Vec3([point_ret.x[0], point_ret.x[1], 0]);
	//		auto retv3D   = Vec3([point_ret.v[0], point_ret.v[1], 0]); 
	//		pot = transform_potentials(center.v/c, 
	//					potential3D(center3D - retx3D, retv3D/c, params.Zt));
	//		foreach(m;-1..2)
	//		{
	//			import nucd.em;
	//			mE1[m] += pot[0] * lq.w * Ylm(1, m, lq.theta, lq.phi);
	//		}	
	//	}
	//	auto t_ret = get_retarded_time(params.h2, center.x, t);
	//	auto point_ret = params.h2.get(t_ret);
	//	pot = transform_potentials(center.v/c, 
	//					potential(center.x - point_ret.x, point_ret.v/c, params.Zp));
	//	multipole_file.writeln(center.tau, " ", pot[0], " ", mE1[-1].abs, " ", mE1[0].abs, " ", center.v.length );
	//}
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

	writeln("hits1 = ", params.h1.n_hit, "   miss1 = ", params.h1.n_miss, "   all lookups1 = ", params.h1.n_all_lookups, "   avoided lookups1 = ", params.h1.n_direct_shortcut, "   shortcuts1 = ", params.h1.n_shortcut);
	writeln("hits2 = ", params.h2.n_hit, "   miss2 = ", params.h2.n_miss, "   all lookups2 = ", params.h2.n_all_lookups, "   avoided lookups2 = ", params.h2.n_direct_shortcut, "   shortcuts2 = ", params.h2.n_shortcut);
}
