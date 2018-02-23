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
			"calc-cross-section",  "giving this will calculate the cross section based on trajectories between b and bmax"                             ,  &params.calc_cross_section
			);

	//if (rslt.helpWanted || args.length == 1)
	//{
	//	defaultGetoptPrinter("nuco - nuclear collision. \n\nexample usage: " ~ args[0] ~ "  --Ap=85 --Zp=35 --At=197 --Zt=79 --E=10 --levelE=0 --levelI=0 --levelE=.6166 --levelI=2 --MEfrom=0 --MEto=1 --MElambda=2 --MEvalue=61.8 --method=relativistic --accuracy=1e-7 --b=10\n\n arguments:"  , rslt.options);
	//	return;
	//}

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

	output_result.output_result(params);	


	//auto bs_ = [1.,2,3,4,5,6,7,8,9];
	//auto as_ = [10.,20,30,40,50,60,70,80,90];
	//cross_section_integrate_curve(bs_,as_);


	// calculate the coulex cross section integral (equation 3.48 in thesis)
	writeln("params.calc_cross_section = ", params.calc_cross_section);
	if (params.calc_cross_section)
	{
		int n_step = 40;
		double[][] sum_level = new double[][](n_step, params.levels.length);
		double[]   bs        = new double[]  (n_step);
		//double[] sum_level_old = new double[](params.levels.length); 
		//double[][] integral_level = new double[](n_steps, params.levels.length); 

		double b_min = params.bp;
		params.debug_on = false;

		double b_old;

		import std.parallelism;
		Parameters[] params_parallel = new Parameters[](n_step);
		//double b_exp  = 1.2;
		//double b_base = params.b_max / b_exp^^(n_step/2-1);
		//b_exp = (params.b_max/b_base)^^(1.0/(n_step-1));
		//if (b_min > 1e-6)
		//{
		//	writeln("b_min>1e-6");
		//	b_base = b_min;
		//	b_exp = (params.b_max/b_base)^^(1.0/(n_step-1));
		//}
		//writeln("b_min = ", b_min);
		//writeln("b_base = ", b_base);
		//writeln("b_exp  = ", b_exp);
		writeln("b_grazing = ", grazing_b, "    b_min = ", b_min, "  .....  ", );
		foreach (n, ref params_p; taskPool.parallel(params_parallel))
		{
			params_p = params;
			//double t = 1.0*(n_step-n)/n_step; // t is a variable in the interval ]0,1]
			//double b = (b_min>0.1)?(b_min*1.2^^n):(0.1*1.2^^n);   // b is in the interval [bmin,\infinity]

			double b;
			//if (grazing_b < 1e6)
			//{
			//	double b_mid = params.b_max/8;
			//	if (n < n_step/2)
			//	{
			//		b = b_min+2*n*(b_mid-b_min)/n_step;
			//	}
			//	else
			//	{
			//		double B = (params.b_max/b_mid)^^(1./(n_step-n_step/2));
			//		b = b_mid*B^^(n-n_step/2);
			//	}			
			//}
			//else 
			//{
			//	if (n < n_step/2)
			//	{
			//		b = b_min+2*n*(2*grazing_b-b_min)/n_step;
			//	}
			//	else
			//	{
			//		double b_left = 2*grazing_b;
			//		double B = (params.b_max/b_left)^^(1./(n_step-n_step/2));
			//		b = b_left*B^^(n-n_step/2);
			//	}
			//}
			//write(b);
			b = b_min+n*(params.b_max-b_min)/n_step;
			params_p.bp = b;
			params_p.integrate_trajectory();

			double sum = 0;
			sum_level[n][] = 0;
			integrate.excite(&ode_excitation, gsl_odeiv2_step_rkf45, params_p);
			foreach(idx,amp;params_p.amplitudes) 
			{
				//writeln("sum_level[", amp.level_index, "] = ", sum_level[amp.level_index]);
				sum_level[n][amp.level_index] += abs(amp.a)^^2;
				if (idx > 0) sum += abs(amp.a)^^2;
			}
			bs[n] = b;

			//writeln(b, "  ", sum*b);
			//if (n != 0) // approximate last two points by an analytic function f(x)=Ca/x^Cb, and integrate that analytically
			//{
			//	foreach(ulong level; 0..params_parallel[n].levels.length)
			//	{

			//		double u0_ = sum_level_old[level]*b_old;
			//		double u1_ = sum_level[level]*b;
			//		double Cb_ = log(u0_/u1_) / log(b/b_old);
			//		double Ca_ = u0_*b_old^^Cb_;
			//		double Cc_ = 1-Cb_;

			//		integral_level[level] += (Ca_/Cc_)*(b^^Cc_ - b_old^^Cc_);
			//		assert(integral_level[level] > 0);
			//	}
			//}
			//b_old   = b;
			//sum_old = sum;
			//sum_level_old[] = sum_level[];

		}

		// write sum of amplitudes to a file
		auto f = File("curve.txt","w+");
	
		import std.algorithm;
		//zip(bs,sum_level).each!(a => f.writeln(a[0], " ", a[1]));//, a[1].each!(b => f.write(b, " ")), f.writeln());
		foreach(zs ; zip(bs, sum_level)) {
			f.write(zs[0], " ");
			foreach(sum; zs[1]) f.write(sum, " "); // this has to be multiplied by b before integrating
			f.writeln();
		}


		foreach(i ; 0..params.levels.length) {
			auto as = new double[](n_step);
			foreach(j; 0..sum_level.length)  as[j] = sum_level[j][i];
			writeln("level ", i, " cross section = ", cross_section_integrate_curve(bs,as)*2*PI, " fm^2  = ", cross_section_integrate_curve(bs,as)*10*2*PI, " mbarn");
		}

		//writeln("integral     = " , 2*PI*integral, " fm^2 = ", 2*PI*integral*10, " mbarn");
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

	auto f = File("curve_int.txt","w+");


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

