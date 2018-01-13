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

void output_result(ref Parameters params)
{


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


	writeln("hits1 = ", params.h1.n_hit, "   miss1 = ", params.h1.n_miss, "   all lookups1 = ", params.h1.n_all_lookups, "   avoided lookups1 = ", params.h1.n_direct_shortcut, "   shortcuts1 = ", params.h1.n_shortcut);
	writeln("hits2 = ", params.h2.n_hit, "   miss2 = ", params.h2.n_miss, "   all lookups2 = ", params.h2.n_all_lookups, "   avoided lookups2 = ", params.h2.n_direct_shortcut, "   shortcuts2 = ", params.h2.n_shortcut);



}