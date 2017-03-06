// This program is used to verify the numerialc expansion into spherical harmonics with lebedev quadrature. 
// It computes for many different points (r_prime) in 3D space one expansion coefficent of the function 1/|r_prime| on a unit sphere around zero.
// This expansion coefficient is calculated analytically and numerically. 
// Both results are printed.


import std.stdio;
import std.math;
import std.random;
import std.complex;
import lebedev_quadrature;
import nucd.geometry;
import types;
import nucd.em;

double f(Vec3 r)
{
	return 1.0/sqrt(r*r);
}

void main()
{

	foreach(i;0..100)
	{	
		double theta = uniform(0,  PI);
		double phi   = uniform(0,2*PI);
		// center of the field-creating charge 
		auto r_prime = eulerVector!real(6,theta,phi);
		writeln(r_prime);

		auto orientation = Orientation(E1);
		auto l = orientation.l;
		foreach(m;orientation)
		{
			writeln(l," : ",m);

			//int l = 2, m = 2;
			// numerical calculation of the Integral
			Complex!real S_numeric = 0;
			double r_len = uniform(0.01,3.0);
			foreach(lq;lq0110)
			{
				auto r = Vec3([lq.x,lq.y,lq.z])*r_len;
				S_numeric += 4*PI* lq.w * f(r - r_prime) * Ylm(l,m, lq.theta, lq.phi);
			}

			Complex!real S_analytic = 4*PI/(2*l+1) * (r_len^^l/r_prime.length^^(l+1)) * Ylm(l,m, theta, phi);


			writeln("numerical result = ", S_numeric);
			writeln("analytic  result = ", S_analytic);

		}
	}
}