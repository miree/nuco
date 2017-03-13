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


module nucd.kinematics;

import std.math;
import std.traits;

///   beta = v/c
///   gamma = 1/sqrt(1-(v/c)^2)
T gamma(T)(T beta) if (isFloatingPoint!T)
{
	return 1/sqrt(1-beta^^2);
}
///   beta = v/c
///   gamma = 1/sqrt(1-(v/c)^2)
T beta(T)(T gamma) if (isFloatingPoint!T)
{
	return sqrt(1-1/gamma^^2);
}
unittest
{
	assert(gamma(0.) == 1.);
	assert(beta(1.) == 0.);

	assert(approxEqual(gamma(beta(1.8)), 1.8, 1e-6));
	assert(approxEqual(gamma(beta(1.5)), 1.5, 1e-6));
	assert(approxEqual(gamma(beta(1.2)), 1.2, 1e-6));
    
	assert(approxEqual(beta(gamma(0.8)), 0.8, 1e-6));
	assert(approxEqual(beta(gamma(0.5)), 0.5, 1e-6));
	assert(approxEqual(beta(gamma(0.2)), 0.2, 1e-6));
}

/// m0c2 = m0*c^2 in units of MeV
/// gamma = 1/sqrt(1-(v/c)^2)
T kinEnergy(T)(T gamma, T m0c2) if (isFloatingPoint!T)
{
	return (gamma - 1)*m0c2;
}
/// m0c2 = m0*c^2 in units of MeV
/// gamma = 1/sqrt(1-(v/c)^2)
T gamma(T)(T kinEnergy, T m0c2) if (isFloatingPoint!T)
{
	return kinEnergy/m0c2 + 1;
}
unittest
{
	double m0c2 = 50.0;
	double x = 30.0;
	assert(approxEqual(gamma(x, m0c2).kinEnergy(m0c2), x, 1e-6));
}

/// center of mass velocity
Vector betaCM(Vector)(double m1, Vector beta1, double m2, Vector beta2)
{
	real gamma1 = gamma(sqrt(beta1*beta1));
	real gamma2 = gamma(sqrt(beta2*beta2));
	return  ( beta1*m1*gamma1 + beta2*m2*gamma2 ) / (m1*gamma1 + m2*gamma2);
}

/// addition of velocities
/// returns beta1 boosted by beta
Vector addVelocity_old(Vector)(Vector beta1, Vector beta)
{
	// unit vectors, parallel and perpendicular to beta
	auto para = beta*(beta1*beta/(beta.length*beta.length));
	auto perp = beta1 - para;
	import std.stdio;
	
	if (para.length > 1e-8) para /= para.length;
	if (perp.length > 1e-8) perp /= perp.length;


	auto b = beta * para;
	auto g = gamma(b);

	
	
	auto b1_para = beta1*para;
	auto b1_perp = beta1*perp;
	

	auto b1_para_boosted = (b1_para + b)/ (    1 + b1_para*b  );
	auto b1_perp_boosted = (b1_perp)    / ( g*(1 + b1_para*b) );
	
	return para*b1_para_boosted + perp*b1_perp_boosted;
}


Vector addVelocity(Vector)(Vector beta1, Vector beta)
{
	if (beta.length==0) return beta1;
	if (beta1.length==0) return beta;

	// unit vectors, parallel and perpendicular to beta
	auto para = beta*(beta1*beta/(beta.length*beta.length));
	auto perp = beta1 - para;
	
	auto para_boosted = (para + beta)/ (                     1 + beta1*beta  );
	auto perp_boosted = (perp)       / ( gamma(beta.length)*(1 + beta1*beta) );
	
	return para_boosted + perp_boosted;
}


import nucd.geometry;
/// returns b2 as seen from an observer moved with b1. 
/// Arguments are in units of c
Vector!(2,T) velocity_addition(T)(Vector!(2,T) b1, Vector!(2,T) b2)
{
	Matrix!(3,3,rowMajor,T) boost;
	boost.boost_direction(b1);//Vector!(2,T)([b1[0],b1[1]]));
	//boost.print();
	
	auto g2 = gamma(b2.length);
	Vector!(3,T) bb2 = Vector!(3,T)([g2, g2*b2[0], g2*b2[1]]);
	import std.stdio;
	auto bb3 = boost*bb2;
//	writeln(boost);
//	writeln("*");
//	writeln(bb2);
//	writeln("=");
//	writeln(bb3);
	return Vector!(2,T)([bb3[1]/bb3[0],bb3[2]/bb3[0]]);
}

// transfrom an angle theta (measured in lab frame)
// into the CM frame. 
real transform_theta(real betaCM, real beta, real theta)
{
	import nucd.geometry;
	auto beta_CM = Vector!(2,real)([betaCM, 0]);
	auto v = Vector!(2,real)([cos(theta),sin(theta)])*beta;
	auto v_CM = velocity_addition(beta_CM, v);
	return  atan2(v_CM[1],v_CM[0]);
}

Vector!(2,T) direction_from_mooved_system(T)(Vector!(2,T) beta, Vector!(2,T) dx)
{
	Matrix!(3,3,rowMajor,T) boost;
	boost.boost_direction(beta);

	Vector!(3,T) x0 = Vector!(3,T)([dx.length,dx[0],dx[1]]);
	auto x1 = boost*x0;

	return Vector!(2,T)([x1[1],x1[2]]);
}
Vector!(3,T) direction3d_from_mooved_system(T)(Vector!(3,T) beta, Vector!(3,T) dx)
{
	Matrix!(4,4,rowMajor,T) boost;
	boost.boost_direction(beta);

	Vector!(4,T) x0 = Vector!(4,T)([dx.length,dx[0],dx[1],dx[2]]);
	auto x1 = boost*x0;

	return Vector!(3,T)([x1[1],x1[2],x1[3]]);
}

Vector!(2,T) transform_direction(T)(Vector!(2,T) r, Vector!(2,T) beta)
{
	import nucd.geometry;
//	write(r, "  ");
	Vector!(2,T) e_beta = beta/beta.length;
	Vector!(2,T) r_para = e_beta*(r*e_beta);
	Vector!(2,T) r_perp = r-r_para;
	r_para *= gamma(beta.length);
	Vector!(2,T) retval = r_para + r_perp;
//	write("    ", retval);
	return retval;
	//import geometry;
	//Mat3 boost;
	//boost.boost_direction(beta);
	//auto rr = Vec3([0,r[0],r[1]]);
	//rr = boost*rr;
	//writeln("      ", rr);
	//return Vector!(2,T)([rr[1],rr[2]]);
}

Vector!(3,T) transform_direction(T)(Vector!(3,T) r, Vector!(3,T) beta)
{
	import nucd.geometry;
	Vector!(3,T) e_beta = beta/beta.length;
	Vector!(3,T) r_para = e_beta*(r*e_beta);
	Vector!(3,T) r_perp = r-r_para;
	r_para *= gamma(beta.length);
	Vector!(3,T) retval = r_para + r_perp;
	return retval;
}

///// returns b2 as seen from an observer moved with b1
//Vector!2 velocity_addition(Vector!2 b1, Vector!2 b2)
//{
//	import geometry;
//	Matrix!(3,3) boost;
//	boost.boost_direction(b1);
//	boost.print();
//	
//	auto g2 = gamma(b2.length);
//	Vector!3 bb2 = Vector!3([g2, g2*b2[0], g2*b2[1]]);
//	import std.stdio;
//	bb2 = boost*bb2;
//	return Vector!2([bb2[1]/bb2[0],bb2[2]/bb2[0]]);
//}

