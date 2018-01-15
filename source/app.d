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
import output_result;


void help(string arg0)
{
	stderr.writeln("call for example: " ~ arg0 ~ " --Ap=85 --Zp=35 --At=197 --Zt=79 --E=10 --levelE=0 --levelI=0 --levelE=.6166 --levelI=2 --MEfrom=0 --MEto=1 --MElambda=2 --MEvalue=61.8 --method=relativistic --accuracy=1e-7 --b=10");
}

void main(string[] args)
{
	// deal with input parameters from the command line
	Parameters params;

	// hard-code nuclear structure information
	//params.levels ~= Parameters.Level(0.0,    0); // 0: ground state
	//params.levels ~= Parameters.Level(0.6166, 2); // 1: first 2+ 
	//params.levels ~= Parameters.Level(1.4361, 4); // 2: first 4+ 
	//params.levels ~= Parameters.Level(15.2, 1); // 2: first 4+ at 1.0 MeV (8 is half spin: 8/2 = l = 4)
	// transitions (only electrical ones so far)
	//params.matrix_elements ~= Parameters.MatrixElement(0,1, 2, complex(61.8)); // E2 transition from ground to 2+1 state 
	//params.matrix_elements ~= Parameters.MatrixElement(1,2, 2, complex(113.6)); // E2 transition from 4+1 to 2+1 state 

	//params.matrix_elements ~= Parameters.MatrixElement(1,2, 2, complex(61.8)); // E2 transition from ground to 2+1 state with 100 e fm^2
	//params.matrix_elements ~= Parameters.MatrixElement(1,1, 1, 17.5); // E2 transition from ground to 2+1 state with 100 e fm^2
	//params.matrix_elements ~= Parameters.MatrixElement(0,2, 2, 50);  // E2 transition from ground to 4+1 state with 100 e fm^2
	//params.matrix_elements ~= Parameters.MatrixElement(1,2, 2, 55.3);  // E2 transition from    2+1 to 4+1 state with 100 e fm^2

	double[] levelsE;
	double[] levelsI;

	int[]    ME_from_indices;
	int[]    ME_to_indices;
	int[]    ME_lambdas;
	double[] ME_values;


	getopt(args,
			"Ap",        &params.Ap,                // mass number of projectile
			"Zp",        &params.Zp,                // atomic number of projectile
			"b",         &params.bp,                // impact parameter
			"beta",      &params.betap,             // velocity (specify this OR the kinetic energy)
			"dmin",		 &params.d_min,             // distance of closest approach
			"E",         &params.Ep,                // kinetic energy (specify this OR the velocity)
			"At",        &params.At,                // mass number of target
			"Zt",        &params.Zt,                // atomic number of target
			"levelE",    &levelsE,                  // add a level with given energy
			"levelI",    &levelsI,                  // add a level with given spin
			"MEfrom",    &ME_from_indices,          // add a reduced matrix element with given from index (into the level array)
			"MEto"  ,    &ME_to_indices,            // add a reduced matrix element with given to   index (into the level array)
			"MElambda",  &ME_lambdas,               // add a reduced matrix element with given lambda
			"MEvalue",   &ME_values,                // add a reduced matrix element with given value in units of [e fm^lambda]
			"CM",        &params.CM,                // if this argument is given, initial velocities will be transformed into center of mass sytem before the simulation
			"rotate",    &params.rotate,            // rotate the resulting trajectories
			"method",    &params.method,            // methods are "relativistic" "classical" "magnetic" (only relativistic makes sense)
			"accuracy",  &params.accuracy,          // accuracy in the trajectory integration
			"timeframe", &params.timeframe,         // simulation is started at -timeframe (in units of [zs] = zeptoseconds)
			"compare-rutherford", &params.compare_rutherford,  
			"compare-rutherford-w", &params.compare_rutherford_w,
			"compare-rutherford-N", &params.compare_rutherford_N,
			"compare-SL-field", &params.compare_SL_field,
			"calc-cross-section", &params.calc_cross_section
			);

	// get the level information
	if (levelsE.length != levelsI.length)
	{
		writeln("ERROR: need same number of --levelE=... and levelI=... arguments");
		return;
	}
	foreach (n;0..levelsE.length)
	{
		params.levels ~= Parameters.Level(levelsE[n],levelsI[n]);
	}

	// get the transitions
	if (ME_from_indices.length != ME_to_indices .length ||
		ME_from_indices.length != ME_lambdas .length    ||
		ME_from_indices.length != ME_values.length  )
	{
		writeln("ERROR: need same number of --MEfrom=... , --MEto=... , --MElambda=... , --MEvalue=... arguments");
		return;
	}
	foreach (n; 0..ME_from_indices.length)
	{
		params.matrix_elements ~= Parameters.MatrixElement(ME_from_indices[n], ME_to_indices[n], ME_lambdas[n], complex(ME_values[n]));
	}

	// preprocess transitions to be in a more useful form for the computation
	params.preprocess_transitions();
	
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
	          *params.betap 
	          *params.Ep 
	          *params.At 
		      *params.Zt))
	{
		stderr.writeln("ERROR: Not all parameters defined 1");
		help(args[0]);
		return;
	}
	if (isNaN(params.bp) && isNaN(params.d_min))
	{
		stderr.writeln("ERROR: Not all parameters defined 2");
		help(args[0]);
		return;
	}
	
	if (params.CM == false && params.compare_rutherford == true)
	{
		stderr.writeln("ERROR: --compare-rutherford works only together with --CM");
		return;
	}
	


	if (!isNaN(params.d_min))
	{
		writeln("calculating impact parameter for distance of closest approach of ", params.d_min, " fm");
		params.bp = 0;
		params.integrate_trajectory();
		if (params.dmin < params.d_min)
		{
			double b_min = 0;
			double b_max = params.d_min;
			do 
			{
				params.bp = 0.5*(b_min+b_max);
				params.integrate_trajectory();
				//write("bp = " , params.bp , "  -> dmin = ", params.dmin, " (d_min requested) = " , params.d_min, "        ::: ");
				if (params.dmin > params.d_min)
				{
					b_max = params.bp;
				}
				else
				{
					b_min = params.bp;
				}

				//writeln("b_min = ", b_min, "     b_max = " , b_max);
			}
			while (b_max-b_min > 1e-8);
			params.bp=0.5*(b_min+b_max);

			// here, we know that we are too fast and have to look for the correct impact parameter

		}
	}
	//else
	{
		params.integrate_trajectory();
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
	

	writeln("distance of closest approach: d_min=" , params.dmin , " fm    at   t_min=", params.t0 , " zs");

	output_result.output_result(params);	

	// calculate the coulex cross section integral (equation 3.48 in thesis)
	writeln("params.calc_cross_section = ", params.calc_cross_section);
	if (params.calc_cross_section)
	{
		double b_min = params.bp;
		params.debug_on = false;

		double sum_old = 0;
		double b_old;

		double integral = 0;
		//double integral_lin = 0;
		int n_step = 20;
		for (int n = 0; n < n_step; ++n)
		{
			double t = 1.0*(n_step-n)/n_step;
			double b = b_min+(1-t^^2)/t^^2;
			write(b);
			params.bp = b;
			params.integrate_trajectory();

			double sum = 0;

			integrate.excite(&ode_excitation, gsl_odeiv2_step_rkf45, params);
			foreach(idx,amp;params.amplitudes) if (idx > 0) sum += abs(amp.a)^^2;

			writeln("  ", sum*b);

			if (t != 1) // approximat last two points by an analytic function f(x)=Ca/x^Cb, and integrate that analytically
			{
				double u0 = sum_old*b_old;
				double u1 = sum*b;
				double Cb = log(u0/u1) / log(b/b_old);
				double Ca = u0*b_old^^Cb;
				double Cc = 1-Cb;

				//writeln("a=",Ca, "b=",Cb);
				integral     += (Ca/Cc)*(b^^Cc - b_old^^Cc);
				//integral_lin += 0.5*(u0+u1)*(b-b_old);
			}
			b_old   = b;
			sum_old = sum;
		}
		writeln("integral     = " , 2*PI*integral, " fm^2 = ", 2*PI*integral*10, " mbarn");
		//writeln("integral_lin = " , integral_lin);
	}


	//import relat_dynamics;
	//writeln("multipols at projectile position");
	//auto infofile = File("infofile.dat", "w+");	
	//foreach(p; params.h1.points)
	//{
	//	real t = p.t;//params.t0;
	//	int  l = 2;
	//	Complex!double[int] Slm = projectile_S_lm(l, params, t);

	//	infofile.write(t, " ", p.tau, " ");
	//	foreach(m;-l..l+1)
	//	{
	//		infofile.write(Slm[m].abs, " " );
	//	}		
	//	infofile.writeln();

	//}

	//foreach(i;-5..6)
	//{
	//	writeln("(-1)^^",i,"=",(-1)^^i, " complex(0,1)^^",i,"=",complex(0,1)^^i);
	//}

	//import integrate;
	//uint N = cast(uint)(2*params.amplitudes.length); // number of independent components: two components for each complex amplitude
	//double[] dydts = new double[N];
	//double[] ys    = new double[N];
	//foreach(idx;0..N) 
	//{
	//	ys[idx]    = 0;
	//	dydts[idx] = 0;
	//}
	//ys[0] = 1.0;  // initialy only ground state is occupied	
	//params.debug_on = true;
	//ode_excitation(0, ys.ptr, dydts.ptr, cast(void*)&params);



	params.debug_on = false;
	integrate.excite(&ode_excitation, gsl_odeiv2_step_rkf45, params);

	//import nucd.em;
	//foreach(Mr;-2..3)
	//{
	//	auto aif = relativistic_coulex_excitation_amplitude(
	//				E2,
	//				0,0,
	//				2,Mr,
	//				61.8, // e*fm^lambda
	//				params.betap,
	//				params.bp,
	//				0.6166, // MeV
	//				params.Zp,
	//				params.Ap,
	//				params.Zt,
	//				params.At
	//				) ;
	//	writeln("Alder/Winther amplitudes: ", aif, "     P_if = ", abs(aif)^^2);
	//}
	//foreach(idx,amp;params.amplitudes)
	//{
	//	writeln(idx,":",amp.a,"    P_if=",abs(amp.a)^^2);
	//}

	//// output for automatic data capture
	write("xxx:");
	write(params.Ep, " ");
	foreach(idx,amp;params.amplitudes) write(abs(amp.a)^^2, " ");
	writeln();
}
