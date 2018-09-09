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
import nucd.em;
import gsl.gsl_odeiv2;
import gsl.gsl_errno;
import std.math;
import std.conv;
import std.getopt;
import std.complex;
import std.range;
import std.file;

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

	auto arglength = args.length;

	auto rslt = getopt(args,
			"Ap|a",                "mass number of projectile "                                                                                        ,   &params.Ap,                // mass number of projectile
			"Zp|z",                "atomic number of projectile "                                                                                      ,   &params.Zp,                // atomic number of projectile
			"impact_parameter|b",  "impact parameter "                                                                                                 ,   &params.bp,                // impact parameter
			"beta|v",              "velocity (specify this OR the kinetic energy) "                                                                    ,   &params.betap,             // velocity (specify this OR the kinetic energy)
			"dmin|d",		       "distance of closest approach "                                                                                     ,   &params.d_min,             // distance of closest approach
			"bmax|B" ,             "when calculating cross sections, trajectories will be calculated up to an impact parameter of bmax"                ,   &params.b_max,
			"energy|e",            "kinetic energy (specify this OR the velocity) "                                                                    ,   &params.Ep,                // kinetic energy (specify this OR the velocity)
			"At|A",                "mass number of target "                                                                                            ,   &params.At,                // mass number of target
			"Zt|Z",                "atomic number of target "                                                                                          ,   &params.Zt,                // atomic number of target
			"levelE|E",            "add a level with given energy "                                                                                    ,   &levelsE,                  // add a level with given energy
			"levelI|I",            "add a level with given spin "                                                                                      ,   &levelsI,                  // add a level with given spin
			"MEfrom|m",            "add a reduced matrix element with given from index (into the level array) "                                        ,   &ME_from_indices,          // add a reduced matrix element with given from index (into the level array)
			"MEto|M"  ,            "add a reduced matrix element with given to   index (into the level array) "                                        ,   &ME_to_indices,            // add a reduced matrix element with given to   index (into the level array)
			"MElambda|l",          "add a reduced matrix element with given lambda "                                                                   ,   &ME_lambdas,               // add a reduced matrix element with given lambda
			"MEvalue|x",           "add a reduced matrix element with given value in units of [e fm^lambda] "                                          ,   &ME_values,                // add a reduced matrix element with given value in units of [e fm^lambda]
			"CM",                  "if this argument is given, initial velocities will be transformed into center of mass sytem before the simulation ",   &params.CM,                // if this argument is given, initial velocities will be transformed into center of mass sytem before the simulation
			"rotate",              "rotate the resulting trajectories "                                                                                ,   &params.rotate,            // rotate the resulting trajectories
			"method",              "methods are \"relativistic\" \"classical\" \"magnetic\" (only relativistic makes sense) "                          ,   &params.method,            // methods are "relativistic" "classical" "magnetic" (only relativistic makes sense)
			"accuracy",            "accuracy in the trajectory integration "                                                                           ,   &params.accuracy,          // accuracy in the trajectory integration
			"timeframe",           "simulation is started at -timeframe (in units of [zs] = zeptoseconds) "                                            ,   &params.timeframe,         // simulation is started at -timeframe (in units of [zs] = zeptoseconds)
			"compare-rutherford",  " ",        &params.compare_rutherford,  
			"compare-rutherford-w"," ",        &params.compare_rutherford_w,
			"compare-rutherford-N"," ",        &params.compare_rutherford_N,
			"compare-SL-field",    " ",        &params.compare_SL_field,
			"calc-cross-section",  "giving this will calculate the cross section based on trajectories between b and bmax"                             ,  &params.calc_cross_section,
			"cross-section-integral-steps", " number of different trajectories between bmin and bmax" , &params.cross_section_integral_steps
			);

	if (rslt.helpWanted || (arglength == 1))
	{
		defaultGetoptPrinter("nuco - nuclear collision. \n\nexample usage: " ~ args[0] ~ "  --Ap=85 --Zp=35 --At=197 --Zt=79 --E=10 --levelE=0 --levelI=0 --levelE=.6166 --levelI=2 --MEfrom=0 --MEto=1 --MElambda=2 --MEvalue=61.8 --method=relativistic --accuracy=1e-7 --b=10\n\n arguments:"  , rslt.options);
		return;
	}

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
	

	double grazing_b;
	// calculate grazing angle
	{
		double dmin = params.d_min; // safe dmin parameter to override with sum of radii
		double bp = params.bp;
		params.d_min = R0(cast(int)params.Ap)+R0(cast(int)params.At);
		calculate_b_from_dmin(params);
		grazing_b = params.bp;
		writeln("grazing impact parameter: ", grazing_b, " fm => dmin = ", params.dmin, " fm (= ", R0(cast(int)params.Ap), " fm + ", R0(cast(int)params.At), " fm )");
		params.d_min = R0(cast(int)params.Ap)+R0(cast(int)params.At)+5;
		calculate_b_from_dmin(params);
		writeln("min. safe Coulex impact parameter: ", params.bp, " fm => dmin = ", params.dmin, " fm (= ", R0(cast(int)params.Ap), " fm + ", R0(cast(int)params.At), " fm + 5 fm )");

		params.d_min = dmin;
		params.bp    = bp;
	}


	if (!isNaN(params.d_min))
	{
		writeln("calculating impact parameter for distance of closest approach of ", params.d_min, " fm");
		calculate_b_from_dmin(params);
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

	// output some results of trajectory calculation
	output_result.output_result(params);	

	// calculate the coulex cross section integral (equation 3.48 in thesis)
	writeln("params.calc_cross_section = ", params.calc_cross_section);
	if (params.calc_cross_section)
	{
		int n_step = params.cross_section_integral_steps;
		double[][] sum_level   = new double[][](n_step, params.levels.length);
		double[]   bs          = new double[]  (n_step);
		double[]   thetalabs   = new double[]  (n_step);
		double[] dthetalab_dbs = new double[]  (n_step);

		double b_min = params.bp;
		params.debug_on = false;

		double b_old;

		import std.parallelism;
		Parameters[] params_parallel = new Parameters[](n_step);
		writeln("b_grazing = ", grazing_b, "    b_min = ", b_min);
		writeln("coulex integral in ", params.cross_section_integral_steps, " steps");
		foreach (n, ref params_p; taskPool.parallel(params_parallel))
		{
			params_p = params;

			// compute reasonalble step size for impact parameter
			double b;
			{
				import std.math;
				double mu = u * params_p.Ap*params_p.At/(params_p.Ap+params_p.At);
				double aa = params_p.Zp*params_p.Zt*alpha*hbar*c/(mu*params_p.betap*params_p.betap*c*c);
				double theta_min = 2*atan2(aa,params_p.b_max);
				double theta_max = 2*atan2(aa,b_min);

				double theta = theta_max - n*(theta_max-theta_min)/(n_step-1);
				b = aa/tan(theta/2);
			}

			//b = b_min+n*(params.b_max-b_min)/n_step;
			double b_plus_delta = b+1e-5;
			double b_minus_delta = b-1e-5;

			params_p.bp = b_plus_delta;
			params_p.integrate_trajectory();
			double thetalab_plus_delta = params_p.thetalab;
			params_p.bp = b_minus_delta;
			params_p.integrate_trajectory();
			double thetalab_minus_delta = params_p.thetalab;
			params_p.bp = b;
			params_p.integrate_trajectory();

			// save the scattering angle and the derivative of the scattering angle
			thetalabs[n] = params_p.thetalab;
			dthetalab_dbs[n] = (thetalab_plus_delta-thetalab_minus_delta)/(b_plus_delta-b_minus_delta);

			double sum = 0;
			sum_level[n][] = 0;
			integrate.excite(&ode_excitation, gsl_odeiv2_step_rkf45, params_p);
			foreach(idx,amp;params_p.amplitudes) 
			{
				//writeln("sum_level[", amp.level_index, "] = ", sum_level[amp.level_index]);
				sum_level[n][amp.level_index] += amp.a.sqAbs;
				if (idx > 0) sum += amp.a.sqAbs;
			}
			bs[n] = b;

			write(".");
			stdout.flush();

		}
		writeln();

		// write sum of amplitudes to a file
		auto f = File("squared_amplitudes.dat","w+");
		import std.algorithm;
		foreach(zs ; zip(bs, sum_level)) {
			f.write(zs[0], " ");
			foreach(sum; zs[1]) f.write(sum, " "); // this has to be multiplied by b before integrating
			f.writeln();
		}

		f = File("diff_xsec.dat","w+");
		import std.algorithm;
		f.writef("#%15s %15s %15s", "b", "thetalab", "dthetalab_db");
		foreach(n; 0..params.levels.length) f.writef(" %20s %20s %30s", "amplitude^2[" ~ to!string(n) ~ "]", "dsigma_db[" ~ to!string(n) ~ "]", "sin(theta)*dsigma_dOmega[" ~ to!string(n) ~ "]", );
		f.writeln("    AW_amplitude^2[1] AW_dsigma_db[1] AW_sin(theta)*dsigma_dOmega[1]");
		foreach(zs ; zip(bs, thetalabs, dthetalab_dbs, sum_level)) {
			double b = zs[0];
			double thetalab = zs[1];
			double dthetalab_db = zs[2];
			f.writef("%15s", b);
			f.writef("%15s", thetalab);
			f.writef("%15s", dthetalab_db);
			foreach(sum; zs[3]) f.writef("%20s %20s %30s ", sum, 2*PI*b*sum, -b*sum*(1./dthetalab_db)); // this has to be multiplied by b before integrating


			double AlderWintherSum = 0;
			for(int M = -cast(int)params.levels[1].L; M <= cast(int)params.levels[1].L; ++M)
			{
				auto aif = relativistic_coulex_excitation_amplitude(
							Multipolarity(Multipolarity.Mode.E, params.matrix_elements[0].lambda),
							0,0,
							cast(int)params.levels[1].L, M,
							sqrt(params.matrix_elements[0].ME.sqAbs), // e*fm^lambda
							params.betap,
							b,
							params.levels[1].E, // MeV
							params.Zp,
							params.Ap,
							params.Zt,
							params.At);
				AlderWintherSum += aif.sqAbs;
			}
			f.writef("%20s %20s %20s ", AlderWintherSum, 2*PI*b*AlderWintherSum, -b*AlderWintherSum*(1./dthetalab_db));
			f.writeln();
		}

		// accurately calculate the cross section integral
		foreach(i ; 0..params.levels.length) {
			auto as = new double[](n_step);
			foreach(j; 0..sum_level.length)  as[j] = sum_level[j][i];
			write("level ", i, " cross section = ", cross_section_integrate_curve(bs,as)*2*PI, " fm^2  = ", cross_section_integrate_curve(bs,as)*10*2*PI, " mbarn");
			if (i == 1) // give Alder/Winther Cross-section for the first excited state
			{
				double AWxsec = relativistic_coulex_cross_section(Multipolarity(Multipolarity.Mode.E, params.matrix_elements[0].lambda),  // multipolarity
												 reducedTransitionStrength(sqrt(params.matrix_elements[0].ME.sqAbs), 2*cast(int)params.levels[0].L),            // reduced transition strength
												 params.betap,         // velocity in units of c, i.e. v/c
												 params.bp,            // impact parameter is approx. distance of closest approach
												 params.levels[1].E,       // energy difference of initial and final level
												 params.Zp,           // charge of excited nucleus
												 params.Ap,           // mass number of excited nucleus
												 params.Zt,           // charge of field providing nucleus
												 params.At);           // mass number of field providing nucleus
				write(" ;  (Alder/Winther straigh line) = ", AWxsec, " fm^2 = ", AWxsec*10, " mbarn");
			}
			writeln;
		}
	}

	params.debug_on = false;
	if (params.levels.length > 0)
	{
		integrate.excite(&ode_excitation, gsl_odeiv2_step_rkf45, params);
		write("xxx:");
		write(params.Ep, " ");
		foreach(idx,amp;params.amplitudes) write(amp.a.sqAbs, " ");
		writeln();	
	}

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


}

double cross_section_integrate_curve(double[] bs1, double[] as1)
{
	auto bs = bs1.dup;
	auto as = as1.dup;

	import gsl.gsl_errno;
	import gsl.gsl_spline;

	assert(bs.length == as.length);
	assert(bs.length > 1);

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	scope(exit) gsl_interp_accel_free(acc);
    gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, bs.length);
    scope(exit) gsl_spline_free(spline);

    for(int i =0; i < as.length; ++i) { as[i] = log(as[i]); }
    gsl_spline_init(spline, bs.ptr, as.ptr, bs.length);

	import std.algorithm;
	double sum = 0;

	//auto f = File("curve_int.txt","w+");


	double db = 0.001;
	for (double b = bs[0]; b <= bs[bs.length-1]-1.5*db; b += db)
	{
		double b0 = b;
		double b1 = b+db;
		double a0 = exp(gsl_spline_eval(spline, b, acc));
		double a1 = exp(gsl_spline_eval(spline, b+db, acc));
		sum += (b1-b0)*(b1*a1+b0*a0)/2;
	}

	//for (int i = 1; i < bs.length; ++i)
	//{
	//	sum += (bs[i]-bs[i-1])*(bs[i]*as[i]+bs[i-1]*as[i-1])/2;
	//}
	//writeln("-------");
	//zip(bs,as).each!(a => writeln(a[0]," ",a[0]*a[1]));
	//zip(bs[1..$],as[1..$]).each!(a => writeln(a[0]," ",a[0]*a[1]));
	return sum;
}

void calculate_b_from_dmin(ref Parameters params, double epsilon = 1e-8)
{
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
			if (params.dmin > params.d_min)
			{
				b_max = params.bp;
			}
			else
			{
				b_min = params.bp;
			}
		}
		while (b_max-b_min > epsilon);
		params.bp=0.5*(b_min+b_max);
	}	
}

