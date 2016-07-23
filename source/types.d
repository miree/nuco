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

import nucd.geometry;

// some types
alias Vector!(2,real) Vec2;
alias Vector!(3,real) Vec3;
alias Matrix!(3,3,rowMajor,real) Mat3;



// some constants
immutable real c   = 299.792458; // speed of light   [fm/zs]
immutable real ahc = 1.43996442; // alpha*hbar*c     [MeV*fm]
immutable real u   = 931.49406;  // atomic mass unit [MeV/c^2]
immutable real mp  = 938.272013; // proton  mass     [MeV/c^2]
immutable real mn  = 939.565346; // neutron mass     [MeV/c^2]
immutable real me  = 0.51099893; // electron mass    [MeV/c^2]
