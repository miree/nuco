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


import parameters;
import types;
import nucd.kinematics;
import nucd.geometry;
import gsl.gsl_odeiv2;
import gsl.gsl_errno;
import std.math;
import std.conv;
import std.stdio;

void integrate(ode, type)(ode func, 
							type T, 
							ref Parameters params)
{
	// setup the system
	auto system = gsl_odeiv2_system(func,   null, 10, cast(void*)&params);

	// setup the integrator
	//auto T       = gsl_odeiv2_step_rkf45;
	auto step    = gsl_odeiv2_step_alloc(T, 10);
	auto control = gsl_odeiv2_control_y_new(params.accuracy, 0.0);
	auto evolve  = gsl_odeiv2_evolve_alloc(10);

	double t1 = -params.timeframe; // [zs] (zeptoseconds)
	double t2 =  params.timeframe; // [zs] 

	double beta1 = params.betap;
	double beta2 = 0;
	
	params.beta_CM = nucd.kinematics.betaCM(params.Ap, params.betap, params.At, 0);
	if (params.CM)
	{
		double beta_CM = nucd.kinematics.betaCM(params.Ap, params.betap, params.At, 0);
		writeln("betaCM = ", beta_CM);

		//beta1 = (beta1-beta_CM)/(1-beta1*beta_CM);
		//beta2 = (beta2-beta_CM)/(1-beta2*beta_CM);
		beta1 = velocity_addition(Vec2([beta_CM,0]),Vec2([beta1,0]) ) [0];
		beta2 = velocity_addition(Vec2([beta_CM,0]),Vec2([beta2,0]) ) [0];
		beta_CM = nucd.kinematics.betaCM(params.Ap, beta1, params.At, beta2);
		writeln("betaCM = ", beta_CM);
		auto Ekin1 = kinEnergy(gamma(beta1),params.Ap*u);
		auto Ekin2 = kinEnergy(gamma(beta2),params.At*u);
		// kinetic energy of the CM moovement is missing in the CM system.
		// This is the kinetic energy of the total system mass (gamma1*m1 + gamma2*m2) mooving with velocity beta_CM
		auto EkinCM = kinEnergy(gamma(params.beta_CM),(gamma(beta1)*params.Ap+gamma(beta2)*params.At)*u);
		writeln("Ekin_CM[MeV]  = ", Ekin1,"+",Ekin2,"+",EkinCM," = ", Ekin1+Ekin2+EkinCM);
		writeln("Ekin_lab[MeV] = ", params.Ep*params.Ap);
	}
	
	
	double[10] y = [t1*beta1*c, params.bp/2,  //  x1,y1
									t1*beta2*c,-params.bp/2,  //  x2,y2
									
									beta1*c,          0.0*c,  // vx1,vy1
									beta2*c,          0.0*c,  // vx2,vy2
									 
									t1/gamma(beta1),          // proper time tau1
									t1/gamma(beta2),          // proper time tau2
								];
	
	//auto stepout  = File( "steps.dat", "w+");
	real angle;
	
	// integrate the ode
	double t = t1;
	double h = 2; // [zs] (zeptoseconds)
	while (t < t2)
	{
		int status = gsl_odeiv2_evolve_apply(evolve, control, step,
											 &system, 
											 &t, t2,
											 &h, y.ptr);
		if (status != GSL_SUCCESS)
			break;

		// call the function once more to get the acceleration and the EM-Fields at the advanced position.
		// Not doing this implies that the last function call is at the desired position, which is (I believe)
		// not guaranteed by the higher order integration procedures.
		double[10] dydt;
		func(t, y.ptr, dydt.ptr, cast(void*)&params);
		
		// add the new step to the history 
		auto x1 = Vec2([y[0],y[1]]);  // projectile (if calculated in lab frame)
		auto v1 = Vec2([y[4],y[5]]);  // projectile (if calculated in lab frame)
		auto x2 = Vec2([y[2],y[3]]);
		auto v2 = Vec2([y[6],y[7]]);
		auto tau1 = y[8];
		auto tau2 = y[9];
		
		params.h1.add(t, x1, v1, params.a1, tau1, params.E21, params.B21, params.pot21);
		params.h2.add(t, x2, v2, params.a2, tau2, params.E12, params.B12, params.pot12);
		
		params.h1.points_partner_ret ~= params.p2_ret;
		params.h2.points_partner_ret ~= params.p1_ret;
		
		//auto beta_CM = betaCM(params.Ap, v1/c, params.At, v2/c);
		//writeln(beta_CM[0], " ", beta_CM[1]);
	}
	
	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(control);
	gsl_odeiv2_step_free(step);

	// do some post-processing

	// search distance of closest approach
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
	while (dt*dt > 1e-16);
	params.t0 = t0;
	params.dmin = dmin;
	params.p1_tau0 = params.h1.get(t0).tau;
	params.p2_tau0 = params.h2.get(t0).tau;


}


// calculate Coulomb excitation amplitudes
void excite(ode, type)(ode func, 
											 type T, 
											 ref Parameters params)
{
	// setup the system
	uint N = cast(uint)(2*params.amplitudes.length); // number of independent components: two components for each complex amplitude
	auto system = gsl_odeiv2_system(func,   null, N, cast(void*)&params);

	// setup the integrator
	//auto T       = gsl_odeiv2_step_rkf45;
	auto step    = gsl_odeiv2_step_alloc(T, N);
	auto control = gsl_odeiv2_control_y_new(params.excitation_accuracy, 0.0);
	auto evolve  = gsl_odeiv2_evolve_alloc(N);

	double t1 = params.h1.points[1].t;   // time of first trajectory point
	double t2 = params.h1.points[$-1].t; // time of last trajectory point

	// extract the step size from trajectory integration
	// to resue them for excitation calculation
	double[] hs;
	for(uint i = 2; i < params.h1.points.length; ++i)
	{
		hs ~= params.h1.points[i].t - params.h1.points[i-1].t;
	}
	
	double[] dydts = new double[N];
	double[] ys    = new double[N];
	foreach(idx;0..N) 
	{
		ys[idx]    = 0;
		dydts[idx] = 0;
	}
	ys[0] = 1.0;  // initialy only ground state is occupied
	
	auto stepout  = File( "amp.dat", "w+");
	
	// integrate the ode
	double t = t1;
	double myt = t1;
	auto Nsteps = hs.length;
	//double stp = 0.0001;
	bool force_h = true;
	double stp;
	foreach(cnt,h; hs)
	{
		if (force_h) stp = h;
		double told = t;
//		write("t=",t, "  -> dt=", stp);
		int status = gsl_odeiv2_evolve_apply(evolve, control, step,
											 &system, 
											 &t, t2,
											 &stp, ys.ptr);
		if (stp < h) force_h = false;
//		writeln("   -> t=",t, "  -> expected t=",told+h);
		if (status != GSL_SUCCESS || t > 50) // nothing will happen 50zs after the collision
			break;

		myt += h;
		// call the function once more to get the acceleration and the EM-Fields at the advanced position.
		// Not doing this implies that the last function call is at the desired position, which is (I believe)
		// not guaranteed by the higher order integration procedures.
		func(t, ys.ptr, dydts.ptr, cast(void*)&params);
		
		double sum = 0;
		foreach(ampl; params.amplitudes)
		{
			import std.complex;
			sum += abs(ampl.a)^^2;
		}
		writeln(force_h?'*':'-'," dt=",stp,"   ",t, " zs:  ", " sum = ", sum-1 , "   a[gs]=", params.amplitudes[0].a, " dadt[gs]=", params.amplitudes[0].dadtau);


		stepout.write(t, " ");
		foreach(y; ys)
		{
			stepout.write(y, " ");
		}
		stepout.writeln();

		//break;
	}
	
	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(control);
	gsl_odeiv2_step_free(step);

	// do some post-processing

}



// calculate Coulomb excitation amplitudes with classical analytic trajectories
void excite_classical_coulex(ode, type)(ode func, 
											 type T, 
											 ref Parameters params)
{
	// setup the system
	uint N = cast(uint)(2*params.amplitudes.length); // number of independent components: two components for each complex amplitude
	auto system = gsl_odeiv2_system(func,   null, N, cast(void*)&params);

	// setup the integrator
	//auto T       = gsl_odeiv2_step_rkf45;
	auto step    = gsl_odeiv2_step_alloc(T, N);
	auto control = gsl_odeiv2_control_y_new(params.excitation_accuracy, 0.0);
	auto evolve  = gsl_odeiv2_evolve_alloc(N);

	double w1 = -50;
	double w2 =  50;

	double[] dydts = new double[N];
	double[] ys    = new double[N];
	foreach(idx;0..N) 
	{
		ys[idx]    = 0;
		dydts[idx] = 0;
	}
	ys[0] = 1.0;  // initialy only ground state is occupied
	
	auto stepout  = File( "amp.dat", "w+");
	
	// integrate the ode
	double w = w1;
	//double stp = 0.0001;
	double h = 0.01;
	while (w < w2)
	{
		int status = gsl_odeiv2_evolve_apply(evolve, control, step,
											 &system, 
											 &t, t2,
											 &h, ys.ptr);
		if (status != GSL_SUCCESS || t > 50) // nothing will happen 50zs after the collision
			break;

		// call the function once more to get the acceleration and the EM-Fields at the advanced position.
		// Not doing this implies that the last function call is at the desired position, which is (I believe)
		// not guaranteed by the higher order integration procedures.
		func(w, ys.ptr, dydts.ptr, cast(void*)&params);
		
		double sum = 0;
		foreach(ampl; params.amplitudes)
		{
			import std.complex;
			sum += abs(ampl.a)^^2;
		}
		writeln("h=",h,"   w=",w, "  ", " sum = ", sum-1 , "   a[gs]=", params.amplitudes[0].a, " dadt[gs]=", params.amplitudes[0].dadtau);


		stepout.write(t, " ");
		foreach(y; ys)
		{
			stepout.write(y, " ");
		}
		stepout.writeln();

		//break;
	}
	
	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(control);
	gsl_odeiv2_step_free(step);

	// do some post-processing

}

