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


import nucd.geometry;

// some types
alias Vector!(2,real) Vec2;
alias Vector!(3,real) Vec3;
alias Vector!(4,real) Vec4;
alias Matrix!(3,3,rowMajor,real) Mat3;



// transfrom the points in the history into a system, where the initial 
// flight direction of the projectile is the z-axis.
Vector!(3,T) to_z_polar(T)(Vector!(3,T) r)
{
	return Vector!(3,T)([r[1],r[2],r[0]]);
}
Vector!(3,T) to_z_polar(T)(Vector!(2,T) r)
{
	return Vector!(3,T)([r[1],0,r[0]]);
}



// some constants
immutable real c   = 299.792458; // speed of light   [fm/zs]
immutable real ahc = 1.43996442; // alpha*hbar*c     [MeV*fm]
immutable real u   = 931.49406;  // atomic mass unit [MeV/c^2]
immutable real mp  = 938.272013; // proton  mass     [MeV/c^2]
immutable real mn  = 939.565346; // neutron mass     [MeV/c^2]
immutable real me  = 0.51099893; // electron mass    [MeV/c^2]
