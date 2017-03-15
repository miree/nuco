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


import types;

enum IntegrationMethod { classical, magnetic, relativistic };

struct Parameters
{
	double Ap;
	double Zp;
	double At;
	double Zt;
	double Ep;
	double betap;
	double bp;
	
	real timeframe = 1000000; // integrate the equations of motion from -timeframe to +timeframe
							  // this is given in units of zeptoseconds (the default is 1 ps)
	
	double beta_CM;    // center of mass velocity in the beginning of the simulation
	bool CM = false;  // true if the velocties should be
					  // transformed into center-of-mass system
					  // before integrating the differential equation
	bool rotate = false; // rotate the trajectory such that it is symmetric
						 // to the x-axis
	bool compare_rutherford = false; // tells the program to output a file containing
									 // simulated trajectories together with the 
									 // non-relativistic Rutherford trajectory
									 // for quantitative comparison
	double compare_rutherford_w = 10; // this will calculate analytic solution for w in [-10,10]
									  // and compare to the simulation
	int compare_rutherford_N = 1000; // determines the number of points calculated in the analytic solution
	
	bool compare_SL_field = false; // compare the EM field with the prediction of the straight-line approximation

	IntegrationMethod method;
	double accuracy = 1e-5; // integration accuracy
	double excitation_accuracy = 1e-6;

	bool debug_on = false;

	Vec2   a1 = Vec2([0,0]);
	Vec2   a2 = Vec2([0,0]);
	Vec2  E21 = Vec2([0,0]);
	real  B21 = 0;
	Vec2  E12 = Vec2([0,0]);
	real  B12 = 0;

	Vec3  pot21   = Vec3([0,0,0]);
	Vec3  pot12   = Vec3([0,0,0]);
	
	History h1;  // store the history of particle 1 (projectile)
	History h2;  // store the history of particle 2 (target)
	
	History.Point p1_ret; // temporary storage for the retarded position of particle 1 as seen from particle 2
	History.Point p2_ret; // temporary storage for the retarded position of particle 2 as seen from particle 1

	real dmin;    // distance of closest approach
	real t0;      // simultation time of closest approach
	real p1_tau0; // proper time of particle 1 of closest approach
	real p2_tau0; // proper time of particle 2 of closest approach

	import std.complex;
	struct Amplitude
	{
		Complex!double a;
		Complex!double dadtau;
		uint ys_index_re;
		uint ys_index_im;
		uint level_index; // index into the input level array
		int L;
		int M; 
		double E; // energy of the corresponding level
		Transition[] transitions; // all transitions that connect this level
	}
	Amplitude[] amplitudes;
	struct Transition
	{
		int L; // half-spin of connected level
		int M; // half-spin of connected level
		double E; // energy of connected level
		int lambda; // lambda=1 for E1, lambda=2 for E2, ... 
		Complex!double ME; // matrix element
		uint idx; // index into the amplitudes array
	}
	bool[] lambdas; // flags that are set to true if a value of labmda (e.g. 2 for E2 transition) is needed
	                   // if only E1 transitions are there, lambdas will contain only 1 value: lambdas[0] = 1.
	                   // if E1 and E3 transitions are present: lambas[0] = 1, lambdas[1] = 3

	// for input only
	struct Level
	{
		double E; // excitation energy
		double L; // integer or half integer
	}
	Level[]      levels;
	struct MatrixElement
	{
		int from;  // index into the levels array, pointing to the 
		int to;    // index into the levels array 
		int lambda; // oder of multipolarity: E2 transitions would have lambda=2, E1 transitions would have lambda=1
		Complex!double ME; // matrix element in units of e*fm^lambda
	}
	MatrixElement[] matrix_elements;
	void preprocess_transitions()
	{
		uint ys_idx = 0;
		foreach(level_idx,level; levels)
		{
			import std.math;
			int twoL = cast(int)lround(level.L*2);
			for(int twoM = -twoL; twoM <= twoL; twoM += 2)
			{
				amplitudes ~= Amplitude(complex(0,0), complex(0,0),
										cast(uint)(ys_idx++), cast(uint)(ys_idx++), cast(uint)(level_idx), 
										twoL, twoM, level.E);
			}
		}
		// now that all amplitudes are present... 
		//   ... include the transitions
		foreach(index_to,   ref ampl_to  ; amplitudes)
		foreach(index_from, ref ampl_from; amplitudes)
		{
			foreach(me; matrix_elements)
			{
				if (me.to == ampl_to.level_index && me.from == ampl_from.level_index )
				{
					if (me.lambda >= lambdas.length) lambdas.length = me.lambda+1;
					lambdas[me.lambda] = true;

					import std.math;
					ampl_to.transitions ~= Transition( ampl_from.L, ampl_from.M, ampl_from.E,
													   me.lambda*2, me.ME, 
													   cast(uint)(index_from));
					ampl_from.transitions ~= Transition( ampl_to.L, ampl_to.M, ampl_to.E,
														 me.lambda*2,  (-1)^^((me.lambda*2+ampl_to.L-ampl_from.L)/2)* me.ME, 
														 cast(uint)(index_to));
														
				}
			}
		}
		foreach(amplitude; amplitudes)
		{
			import std.stdio;
			writeln(amplitude);
		}

		foreach(i,lambda;lambdas)
		{
			import std.stdio;
			writeln("lambda",i,"=",lambda);
		}
	}
}

struct Spline
{
	import gsl.gsl_spline;
	gsl_interp_accel *acc;
	gsl_spline       *spline;
	this(double[] ts, double[] ys)
	{
		assert(ts.length == ys.length);
		auto N = ts.length;
		acc    = gsl_interp_accel_alloc();
		spline = gsl_spline_alloc(gsl_interp_cspline, N);
		gsl_spline_init(spline, ts.ptr, ys.ptr, N);
	}
	double get(real t)
	{
		return gsl_spline_eval(spline, t, acc);
	}
	~this()
	{
		gsl_spline_free(spline);
		gsl_interp_accel_free(acc);
	}
}


struct History // This could also be called "Trajectory", but since it is used to do position lookups at earlier times, the name "History" seems appropriate.
{
	// helper value to accelerate the lookup of retarded time
	// (it stores the result of the last calculated retarded time)
	real last_t_ret;
	
	// helper value to remember the index of the where the last history lookup.
	long last_index_low;
	long last_index_high;
	
	int n_hit, n_miss, n_all_lookups, n_direct_shortcut, n_shortcut;

	struct Point
	{
		real t;   // time in the simulation frame
		Vec2 x;   // position
		Vec2 v;   // velocity
		Vec2 a;   // acceleration
		real tau; // proper time (Eigenzeit), and its first and second derivative
		real dtaudt;
		real d2taudt2;
		
		// the spline coefficients are buffered inside the supporting points
		// to make computation faster. the left(earlier) point stores the spline coefficients
		// for the interpolation between left and right supports.
		real a0,a1,a2;
		real b0,b1,b2;
		real c0,c1,c2;
		real d0,d1,d2;
		real e0,e1,e2;
		real f0,f1,f2;	
		bool coefficients_available = false;
				
	}
	Point[] points;
	Point[] points_partner_ret; // the retarded time/position/.../ of the collision partner
	
	// For the quantities tau, E, B no first or second 
	// derivative is known at the time of computation.
	// Therefore natural splines have to be used for interpolation
	// and GSL requires arrays for interpolation.
	// (that is why the following is not put in a struct)
	//real ts[]
	//real taus[];
	//real Exs[];
	//real Eys[];
	//real Bs[];
	
	struct Field
	{
		real tau;
		
		Vec2   E;
		real   B;
		
		Vec3 pot;
	}
	Field[] fields;
	
	// add a point into the history
	void add(real t, Vec2 x, Vec2 v, Vec2 a, 
			 real tau, Vec2 e, real b, Vec3 pot)
	in
	{
		// it is not allowed to add points that are earlyer than the
		// latest history entry 
		if (points.length > 0)
		{
			assert(t > points[$-1].t);
		}
	}
	body
	{
		import nucd.kinematics;
		// add position, velocity, acceleration, tau, first and second derivative dtaudt, d2taudt2
		points ~= Point(t,x,v,a, tau,1.0/gamma(v.length/c),-gamma(v.length/c)*v.length*a.length/c^^2);
		fields ~= Field(tau,e,b,pot);
		//ts   ~= t;
		//taus ~= tau;
		//Exs  ~= e[0];
		//Eys  ~= e[1];
		//Bs   ~= b;
	}
	
	// get interpolated values of position and field from the history
	Point get(real t)
	in
	{
		assert(points.length > 0);
	}
	body
	{
		if (t <= points[0].t) // project into the past (before the first recorded point)
		{
			return Point(t, 
						 points[0].x - points[0].v*(points[0].t-t),
						 points[0].v,
						 Vec2([0,0]));
		}
		if (t >= points[$-1].t) // project into the future
		{
			return Point(t,
						 points[$-1].x + points[$-1].v*(t-points[$-1].t),
						 points[$-1].v,
						 Vec2([0,0]));
		}
		
		
		auto pa = &points[last_index_low];  // the two supporting points between which 
		auto pb = &points[last_index_high]; // we have to do a spline interpolation
		
		// see if we are still at the same point in the history table
		// doing this takes 70% of the time than not 
		//if (last_index_high - last_index_low == 1)
		if (!(pa.t < t && pb.t > t))
		{
			++n_all_lookups;
			long index_low  = 0;
			long index_high = points.length-1;
			// try to see if expanding the interval by 1 works (a small optimization)
			if (last_index_low != 0 && last_index_high != points.length)
			{
				pa = &points[last_index_low-1];
				pb = &points[last_index_high+1];
				if (pa.t < t && pb.t > t)
				{
					index_low = last_index_low-1;
					index_high = last_index_high+1;
					++n_shortcut;
				}			
			}

			// find two nearest neigbors to the requested time 
			// by means of bisection
			int n_lookup = 0;
			while(index_high - index_low > 1)
			{
				++n_lookup;
				long index_mid = (index_high + index_low)/2;
				real t_mid   = points[index_mid].t;
				
				if (t_mid > t)
				{
					index_high = index_mid;
				}
				else
				{
					index_low = index_mid;
				}
			}
			// remember the indices for the next time
			last_index_low  = index_low;
			last_index_high = index_high;
			
			pa = &points[index_low];   
			pb = &points[index_high];
		}
		else
		{
			++n_direct_shortcut;
		}

		if (!pa.coefficients_available)
		{
			// 5-th order spline interpolation between the two closest points pa and pb.
			real t2 = pb.t-pa.t;
			real t2_p2 = t2*t2;    // = t2^^2
			real t2_p3 = t2_p2*t2; // = t2^^3
			real t2_p4 = t2_p3*t2; // = t2^^4
			real t2_p5 = t2_p4*t2; // = t2^^5
			// coefficients for the x-coordinate
			pa.a0 = pa.x[0];
			pa.b0 = pa.v[0];
			pa.c0 = pa.a[0]/2.;
			pa.d0 =  (20*pb.x[0]-20*pa.x[0]- 8*t2*pb.v[0]-12*t2*pa.v[0]+(  pb.a[0]-3*pa.a[0])*t2_p2)/(2*t2_p3);
			pa.e0 = -(30*pb.x[0]-30*pa.x[0]-14*t2*pb.v[0]-16*t2*pa.v[0]+(2*pb.a[0]-3*pa.a[0])*t2_p2)/(2*t2_p4);
			pa.f0 =  (12*pb.x[0]-12*pa.x[0]- 6*t2*pb.v[0]- 6*t2*pa.v[0]+(  pb.a[0]-  pa.a[0])*t2_p2)/(2*t2_p5);
	
			// coefficients for the y-coordinate
			pa.a1 = pa.x[1];
			pa.b1 = pa.v[1];
			pa.c1 = pa.a[1]/2.;
			pa.d1 =  (20*pb.x[1]-20*pa.x[1]- 8*t2*pb.v[1]-12*t2*pa.v[1]+(  pb.a[1]-3*pa.a[1])*t2_p2)/(2*t2_p3);
			pa.e1 = -(30*pb.x[1]-30*pa.x[1]-14*t2*pb.v[1]-16*t2*pa.v[1]+(2*pb.a[1]-3*pa.a[1])*t2_p2)/(2*t2_p4);
			pa.f1 =  (12*pb.x[1]-12*pa.x[1]- 6*t2*pb.v[1]- 6*t2*pa.v[1]+(  pb.a[1]-  pa.a[1])*t2_p2)/(2*t2_p5);
			
			// coefficients for tau
			pa.a2 = pa.tau;
			pa.b2 = pa.dtaudt;
			pa.c2 = pa.d2taudt2/2.;
			pa.d2 =  (20*pb.tau-20*pa.tau- 8*t2*pb.dtaudt-12*t2*pa.dtaudt+(  pb.d2taudt2-3*pa.d2taudt2)*t2_p2)/(2*t2_p3);
			pa.e2 = -(30*pb.tau-30*pa.tau-14*t2*pb.dtaudt-16*t2*pa.dtaudt+(2*pb.d2taudt2-3*pa.d2taudt2)*t2_p2)/(2*t2_p4);
			pa.f2 =  (12*pb.tau-12*pa.tau- 6*t2*pb.dtaudt- 6*t2*pa.dtaudt+(  pb.d2taudt2-  pa.d2taudt2)*t2_p2)/(2*t2_p5);
			
			pa.coefficients_available = true;
			++n_miss;
		}
		else
		{
			++n_hit;
		}
		
		// calculate x and y results
		real dt = t-pa.t;
		real dt_p2 = dt*dt;     // = dt^^2
		real dt_p3 = dt_p2*dt;  // = dt^^3
		real dt_p4 = dt_p3*dt;  // = dt^^4
		real dt_p5 = dt_p4*dt;  // = dt^^5
		real x  = pa.a0 + pa.b0*dt +   pa.c0*dt_p2 +   pa.d0*dt_p3 +    pa.e0*dt_p4 +    pa.f0*dt_p5;
		real y  = pa.a1 + pa.b1*dt +   pa.c1*dt_p2 +   pa.d1*dt_p3 +    pa.e1*dt_p4 +    pa.f1*dt_p5;
		real vx =         pa.b0    + 2*pa.c0*dt    + 3*pa.d0*dt_p2 +  4*pa.e0*dt_p3 +  5*pa.f0*dt_p4;
		real vy =         pa.b1    + 2*pa.c1*dt    + 3*pa.d1*dt_p2 +  4*pa.e1*dt_p3 +  5*pa.f1*dt_p4;
		real ax =                    2*pa.c0       + 6*pa.d0*dt    + 12*pa.e0*dt_p2 + 20*pa.f0*dt_p3;
		real ay =                    2*pa.c1       + 6*pa.d1*dt    + 12*pa.e1*dt_p2 + 20*pa.f1*dt_p3;
		real tau= pa.a2 + pa.b2*dt +   pa.c2*dt_p2 +   pa.d2*dt_p3 +    pa.e2*dt_p4 +    pa.f2*dt_p5;

		//import std.stdio;
		//writeln("n_lookup = ", n_lookup);

		return Point(t,
					 Vec2([x,y]),
					 Vec2([vx,vy]),
					 Vec2([ax,ay]),
					 tau);
	}
	
	//Field get_field(double tau)
	//{
	//	
	//}
	
};

