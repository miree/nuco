/* 
 * Copyright (C) 2015,2016 Michael reese
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

/*
 * Author:  M. Reese
 */

module nucd.em;

import std.stdio;
import std.complex;
import std.math;
import std.algorithm;


//Complex!double iPow(int i)
//{
//	writeln("iPow(",i,")");
//	switch(i%4)
//	{
//		case 0: return Complex!double( 1, 0);
//		case 1: return Complex!double( 0, 1);
//		case 2: return Complex!double(-1, 0);
//		case 3: return Complex!double( 0,-1);
//		default : {}
//	}
//	return Complex!double.init;
//}
//
double fac(int n)
{
	if (n < 0) return 0;
	static real[] table = [1,1,2];
	if (n >= table.length)
		for (auto i = table.length; i <= n; ++i)
			table ~= table[$-1]*i;
	return table[n];		
}
double facc(int n)
{
	if (n < 0) return 0;
	static real[] table = [1,1,2];
	if (n >= table.length)
		for (auto i = table.length; i <= n; ++i)
			table ~= table[$-2]*i;
	return table[n];		
}


double Plm(int lambda, int mu, double x)
{
	import gsl.gsl_sf_legendre;

	if (mu > lambda)
		return 0;
	if (lambda < 0)
		return 0;
	if (mu < 0)
		return 0;
	
	writeln(lambda-1," ",mu," ",x);
	return gsl_sf_legendre_Plm(lambda,mu,x);
}

/// Type of multipole radiation, containing mode (electric or magnetic,
/// i.e. E or M) and the order of the multipole (1 = dipole, 
/// 2 = quadrupole, 3 = octupole, 4 = hexadecupole, ...)
/// Mode and order are typically labeled \sigma \lambda in nuclear 
/// physics text books.
struct Multipolarity
{
private:
	Mode mode;
	int order;

public:
	enum Mode
	{
		E  = 1,
		M  = 2
	}
	this(Mode md, int or)
	{
		mode  = md;
		order = or;
	}
}

/// Predefined lowest order multipoles
immutable
{
	auto E0 = Multipolarity(Multipolarity.Mode.E,0);
	auto E1 = Multipolarity(Multipolarity.Mode.E,1);
	auto E2 = Multipolarity(Multipolarity.Mode.E,2);
	auto E3 = Multipolarity(Multipolarity.Mode.E,3);
	auto E4 = Multipolarity(Multipolarity.Mode.E,4);
	auto E5 = Multipolarity(Multipolarity.Mode.E,5);
	auto E6 = Multipolarity(Multipolarity.Mode.E,6);
	auto E7 = Multipolarity(Multipolarity.Mode.E,7);

	auto M0 = Multipolarity(Multipolarity.Mode.M,0);
	auto M1 = Multipolarity(Multipolarity.Mode.M,1);
	auto M2 = Multipolarity(Multipolarity.Mode.M,2);
	auto M3 = Multipolarity(Multipolarity.Mode.M,3);
	auto M4 = Multipolarity(Multipolarity.Mode.M,4);
	auto M5 = Multipolarity(Multipolarity.Mode.M,5);
	auto M6 = Multipolarity(Multipolarity.Mode.M,6);
	auto M7 = Multipolarity(Multipolarity.Mode.M,7);
}

/// allow to use foreach with multipolarities and iterate over all
/// magnetic substates
struct Orientation
{
private:
	Multipolarity mult;
	int           m;
public:
	this(Multipolarity mlt)
	{
		mult =  mlt;
		m    = -mlt.order;
	}
	@property bool    empty() { return m > mult.order; }
	@property ref int front() { return m; }
	void           popFront() { ++m; }
}

/// The G function as defined in [1]. This is just hard coded.
/// A generic implementation using compile-time evaluation would be cool.
/// [1] A.Winter, K.Alder, Relativistic Coulomb Excitation, Nuclear Physics A319 (1979) 518 - 532
Complex!real G(Multipolarity mult, int mu, double x)
{
	switch(mult.mode)
	{
		case Multipolarity.Mode.E:
		if (mu < 0)	return (-1)^^mu * G(mult, -mu, x);
		switch(mult.order)
		{
			case 1:
				switch(mu)
				{
					case 0:  return complex(0,-4./3.*sqrt(PI)*sqrt(x^^2-1));
					case 1:  return complex(1./3.*sqrt(8*PI)*x);
					default: return complex!real(0.);
				}
			case 2:
				switch(mu)
				{
					case 0:  return complex(2./5.*sqrt(PI)*x*sqrt(x^^2-1));
					case 1:  return complex(0, 2./5.*sqrt(1./6.*PI)*(2*x^^2-1));
					case 2:  return complex(-2./5.*sqrt(1./6.*PI)*x*sqrt(x^^2-1));
					default: return complex!real(0.);
				}
			case 3:
				switch(mu)
				{
					case 0:  return complex(0,2./105.*sqrt(PI)*(5*x^^2-1)*sqrt(x^^2-1));
					case 1:  return complex(-1./105*sqrt(1./3.*PI)*x*(15*x^^2-11));
					case 2:  return complex(0,-1./21.*sqrt(2./15.*PI)*(3*x^^2-1)*sqrt(x^^2-1));
					case 3:  return complex(1./21.*sqrt(1./5.*PI)*x*(x^^2-1));
					default: return complex!real(0.);
				}
			default: return complex!real(0.);
		}
		case Multipolarity.Mode.M:
		if (mu < 0)	return -(-1)^^mu * G(mult, -mu, x);
		switch(mult.order)
		{
			case 1:
				switch(mu)
				{
					case 0:  return complex!real(0.);
					case 1:  return complex(0.,-1./3.*sqrt(8*PI));
					default: return complex!real(0.);
				}
			case 2:
				switch(mu)
				{
					case 0:  return complex!real(0.);
					case 1:  return complex(2./5.*sqrt(1./6.*PI)*x);
					case 2:  return complex(0., 2./5.*sqrt(1./6.*PI)*sqrt(x^^2-1));
					default: return complex!real(0.);
				}
			default: return complex!real(0.);
		}		
		default:
		return Complex!real(0.);
	}
}

double g_mu_xi(int mu, double xi)
{
	if (mu < 0) return g_mu_xi(-mu,xi);
	
	import gsl.gsl_sf_bessel : gsl_sf_bessel_Kn;
	import std.math : PI;
	
	double K_mu        = gsl_sf_bessel_Kn(mu  , xi);
	double K_mu_plus_1 = gsl_sf_bessel_Kn(mu+1, xi);
	
	return PI*xi^^2 * (
				K_mu_plus_1^^2  -  K_mu^^2
			  - (2*mu/xi) * K_mu_plus_1 * K_mu 
		              );
}

// calculates cross section in fm^2
double relativistic_coulex_cross_section(Multipolarity mult,  // multipolarity
										 double B,            // reduced transition strength
										 double beta,         // velocity in units of c, i.e. v/c
										 double R,            // impact parameter is approx. distance of closest approach
										 double deltaE,       // energy difference of initial and final level
										 double Ze,           // charge of excited nucleus
										 double Ae,           // mass number of excited nucleus
										 double Zf,           // charge of field providing nucleus
										 double Af)           // mass number of field providing nucleus
{
	import std.math;
	import std.complex;

	double k     = deltaE/197.33;      // 1/fm
	double gamma = 1./sqrt(1-beta^^2); // 1
	//double xi    = k*R/gamma/beta;     // 1

	double A0    = Ae*Af/(Ae+Af);
	double xi    = k/gamma/beta * (R + PI/2. * Ze*Zf*(197.33/137.04)/(A0*940*gamma*beta*beta));
	
	double sigma = k^^(2*(mult.order-1)) * B * (Zf/137.04)^^2;
	
	double sum = 0;
	foreach(mu;Orientation(mult)) // -order,-order+1,...,order-1,order
	{
		double G = abs(G(mult,mu,1/beta))^^2;
		double g = g_mu_xi(mu,xi);
		//writefln("%s   %s   %s", mu, GE2, g);
		sum += G*g;
	}
	sigma *= sum;
	
	return sigma;
}	

double reducedMatrixElement(double reducedTransitionStrength, int twoIi)
{
	import std.math;
	return sqrt((twoIi + 1) * reducedTransitionStrength);
}	

double reducedTransitionStrength(double reducedMatrixElement, int twoIi)
{
	import std.math;
	return reducedMatrixElement^^2 / (twoIi + 1);
}

//Complex!double relativistic_coulex_transition_amplitude(
//	Multipolarity mult,  // multipolarity
//	double B,            // reduced transition strength
//	double beta,         // velocity in units of c, i.e. v/c
//	double R,            // impact parameter is approx. distance of closest approach
//	double deltaE,       // energy difference of initial and final level
//	double Ze,           // charge of excited nucleus
//	double Ae,           // mass number of excited nucleus
//	double Zf,           // charge of field providing nucleus
//	double Af)           // mass number of field providing nucleus
//{
//	import std.math;
//	import std.complex;
//	import gsl.gsl_sf_bessel : gsl_sf_bessel_Kn;
//
//	double k     = deltaE/197.33;      // 1/fm
//	double gamma = 1./sqrt(1-beta^^2); // 1
//
//	double A0    = Ae*Af/(Ae+Af);
//	double xi    = k/gamma/beta * (R + PI/2. * Ze*Zf*(197.33/137.04)/(A0*940*gamma*beta*beta));
//
//	double a = k^^(2*(mult.order-1)) * B * (Zf/137.04)^^2;
//	
//	double sum = 0;
//	foreach(mu;Orientation(mult))
//	{
//		double G = abs(G(mult,mu,1/beta))^^2;
//		double g = g_mu_xi(mu,xi);
//		//writefln("%s   %s   %s", mu, GE2, g);
//		sum += G*g;
//	}
//	sigma *= sum;
//	
//	return sigma;
//	
//}

auto clebsch_gordan(real j1, real m1, real j2, real m2, real j3, real m3)
{
	import gsl.gsl_sf_coupling;
	
	int _j1 = cast(int)(round(2*j1));
	int _j2 = cast(int)(round(2*j2));
	int _j3 = cast(int)(round(2*j3));

	int _m1 = cast(int)(round(2*m1));
	int _m2 = cast(int)(round(2*m2));
	int _m3 = cast(int)(round(2*m3));
	
	return ((((_j1-_j2+_m3)/2)%2)?(-1):(1)) * sqrt(_j3 + 1.)
			* gsl_sf_coupling_3j(_j1, _j2, _j3, _m1, _m2, _m3);
}


double coefficientF(int k, int twoL1, int twoL2, int twoI1, int twoI2)
{
	import std.math;
	import gsl.gsl_sf_coupling;
	import gsl.gsl_errno;
	
	if (k == 0)
	{
		return (twoL1==twoL2)?1.0:0.0;
	}
	
	double A = (-1)^^((twoI1+twoI2)/2 - 1);
	double B = sqrt( cast(real)((twoL1 + 1)*(twoL2 + 1)*(twoI2 + 1)*(2*k + 1)) );
	double C = gsl_sf_coupling_3j(twoL1, twoL2, 2*k,
									   2,    -2,   0);
	double D = gsl_sf_coupling_6j(twoL1, twoL2, 2*k,
									twoI2, twoI2, twoI1);
	return A*B*C*D;							  
}
unittest
{
	import std.stdio;
 
 	//writeln("F-coefficient unittest");
	assert(approxEqual( coefficientF(0, 0,0,0,0),  1.0    ,1e-3));

	assert(approxEqual( coefficientF(2, 4,4,8,4), -0.1707 ,1e-3));
	assert(approxEqual( coefficientF(2, 4,4,0,4), -0.5976 ,1e-3));

	assert(approxEqual( coefficientF(2, 10,8,14,14), -0.209937 ,1e-3));


	assert(approxEqual( coefficientF(4, 4,4,8,4), -0.0085 ,1e-2));
	assert(approxEqual( coefficientF(4, 4,4,0,4), -1.0690 ,1e-2));
	//writeln("passed");
}


double relativisticCoulexAngularDistributionW(double theta, 
											  Multipolarity mult,
											  int twoIi,
											  int twoIf,
											  int twoIff,
											  double deltaE,
											  double beta,
											  double Ae,
											  double Ze,
											  double Af,
											  double Zf,
											  double R) // R = minimum distance
{
	import gsl.gsl_sf_coupling;
	import gsl.gsl_sf_legendre;
	
	double kk     = deltaE/197.33;      // 1/fm
	double gamma  = 1./sqrt(1-beta^^2); // 1
	double A0     = Ae*Af/(Ae+Af);
	double xi     = kk/gamma/beta * (R + PI/2. * Ze*Zf*(197.33/137.04)/(A0*940*gamma*beta*beta));
	
	int  lambda = mult.order;
	int      l1 = 2;
	int 	 l2 = 2;
	
	double W = 0;
	for(int k = 0; k <= 4; k += 2)
	{
		foreach(mu;Orientation(mult))
		{
			double G      = abs(G(mult,mu,1/beta))^^2;
			double g      = g_mu_xi(mu,xi);
			double phase  = (-1)^^mu;
			//writeln(2*lambda, "   ", 2*mu    , "   ", 2*lambda, "   ", -2*mu    , "   ", 2*k, "   ", 0    );
			double w3j    = gsl_sf_coupling_3j(2*lambda, 2*lambda, 2*k, 2*mu    , -2*mu    , 0    );
			double w6j    = gsl_sf_coupling_6j(  twoIf ,   twoIf , 2*k, 2*lambda,  2*lambda, twoIi);
			double F      = coefficientF(k, 2*l1,2*l2, twoIff, twoIf);
			W += G*g*phase*w3j*w6j*F*sqrt(2.*k+1)*gsl_sf_legendre_Pl(k,cos(theta));
		}
	}
	return W;
}											  



double R0(int A)
{
	return 1.25 * A^^(1./3.);
}

/// spherical harmonics. Angles theta and phi are in rad.
Complex!real Ylm(int l, int m, double theta, double phi)
{
	import gsl.gsl_sf_legendre;
	
	int factor = 1;
	bool cc = false;
	if (m < 0) 
	{
		m = -m;
		cc = true;
		if (m%2 == 1)
			factor *= -1;
	}
	
	if (m > l) 
	{
		return Complex!real();
	}
	//writeln("factor ", factor, " " , m%2);
	auto Plm = gsl_sf_legendre_sphPlm(l,m,cos(theta));
	auto result = factor * Plm * std.complex.expi(m*phi);
	if (cc)
		return result.conj();
	return result;
}

/// Spherical harmonic squared. Angle theta is  in rad.
double Ylm_2(int l, int m, double theta)
{
	import gsl.gsl_sf_legendre;

	if (m < 0) 
	{
		m = -m;
	}
	if (m > l) 
	{
		return 0;
	}
	auto Plm = gsl_sf_legendre_sphPlm(l,m,cos(theta));
	return Plm*Plm;
}


unittest
{
	//writeln("Sperical harmonics squared unittest");
	int nTheta = 20;
	int nL = 20;
	foreach(t;0..nTheta)
	{
		double theta = t*PI/nTheta;
		foreach(l;0..nL)
		{
			double sum = 0;
			foreach(m;-l..l+1)
			{
				sum += Ylm_2(l,m, 1);
			}
			assert(approxEqual((2*l+1)/(4*PI), sum));
		}
	}
	//writeln("passed");
}

/// Angular distribution of a pure multipole transition of order l and 
/// orientation m. Angle theta is in rad.
double Flm(int l, int m, double theta)
{
	double sum = 0;
	sum += (l*(l+1) - m*(m+1))*Ylm_2(l,m+1,theta);
	sum += (l*(l+1) - m*(m-1))*Ylm_2(l,m-1,theta);
	sum +=  2*m*m * Ylm_2(l,m,theta);
	
	return sum/(2*l*(l+1));
}


unittest
{
	//writeln("F_lm function unittest");
	immutable nTheta = 1000;
	immutable lMax = 20;
	double[lMax] max;
	double[lMax] min;
	double[lMax] norm = 0;
	foreach(l;0..lMax)
	{
		min[l] = double.max;
		max[l] = 0;
	}
	foreach(t;0..nTheta)
	{
		double theta = t*PI/nTheta;
		assert(approxEqual(Flm(1, 0,theta),  3*sin(theta)^^2 / (8*PI)));
		assert(approxEqual(Flm(1, 1,theta),  3*(1+cos(theta)^^2)/16/PI));
		assert(approxEqual(Flm(1,-1,theta),  3*(1+cos(theta)^^2)/16/PI));
		assert(approxEqual(Flm(2, 0,theta), 15*sin(theta)^^2*cos(theta)^^2 / ( 8*PI)));
		assert(approxEqual(Flm(2, 1,theta),  5*(1-3*cos(theta)^^2+4*cos(theta)^^4) / (16*PI)));
		assert(approxEqual(Flm(2,-1,theta),  5*(1-3*cos(theta)^^2+4*cos(theta)^^4) / (16*PI)));
		assert(approxEqual(Flm(2, 2,theta),  5*(1-cos(theta)^^4) / (16*PI)));
		assert(approxEqual(Flm(2,-2,theta),  5*(1-cos(theta)^^4) / (16*PI)));
		
		// checking for isotropy of the sum over all m
		foreach(l;1..lMax) 
		{
			double sum = Flm(l,0,theta);
			foreach(m;1..l+1)
			{
				sum += 2*Flm(l,m,theta);
			}

			min[l] = std.algorithm.min(sum,min[l]);
			max[l] = std.algorithm.max(sum,max[l]);

			// integral over the full sphere
			norm[l] += 2*PI * sin(theta) * PI/nTheta * Flm(l,1,theta); 
		}
	}
	foreach(l;1..lMax)
	{
		assert(approxEqual(norm[l],1, 1e-3));
		assert(approxEqual(max[l]-min[l], 0));
	}
	//writeln("passed");
} 


/// decay from level I1 to I2 with mixed radiation of L1 and L2
/// with the mixing ratio delta = <L2>/<L1>
double coefficientA(int k, double delta, int twoL1, int twoL2, int twoI1, int twoI2)
{
	return (  coefficientF(k, twoL1, twoL1, twoI1, twoI2)
			+ coefficientF(k, twoL1, twoL2, twoI1, twoI2) * 2*delta
			+ coefficientF(k, twoL2, twoL2, twoI1, twoI2) * delta^^2 
			)  /  (1 + delta^^2);
}

/// uses program Dweiko to compute multistep Coulex cross section
string dweikoCrossSection(	double E_kin,
							double b_min,
							double E1   ,
							double M1   ,
							double E2   ,
							double M2   ,
							double Ii   ,
							double If   )
{
	import std.stdio;
	import std.process;
	import std.array;
	import std.algorithm;
	import std.array;
	import std.conv;
	import std.file;
	import std.uni;

	// produce dweiko input file from template
	File("dweiko.in","w+")
		.writef(readText("dweiko.in.template"), 
			E_kin, b_min, Ii, If, E1, E2, M1, M2);
	
	// run dweiko
	auto dweiko = execute("dweiko"); 

	// extract cross section
	return File("dweiko.out")
			.byLine(KeepTerminator.yes)
			.map!(a=>a.idup)
			.array[$-2]			        // take the second last line of dweiko's output
			.split('=')[$-1][0..$-1]    // take the last string after the equals sign and remove the newline
			.strip(' ')
			.toLower;
}
