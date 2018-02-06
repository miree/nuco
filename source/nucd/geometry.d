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


module nucd.geometry;

import std.math;
import std.stdio;
import std.random;
import std.conv;

struct Vector(uint dim, T = double)
{
	private:
		T[dim] content = 0;
	public:
		this(in T[dim] init)
		{
			content[] = init[];
		}
		
		//this(this)
		//{
		//	content = content.dup;
		//}
		
		this(uint dim2)(in Vector!(dim2,T) vec)
		if (dim2 < dim)
		{
			content[0..dim2] = vec.content[0..dim2];
			content[dim2..$] = 0;
		}
		
		ref Vector unity(uint i)
		{
			assert (i < dim);
			content[] = 0;
			content[i] = 1;
			return this;
		}
		
		ref T opIndex(uint i)
		{
			assert(i < dim);
			return content[i];
		}
		void print()
		{
			foreach(x; content)
				write(x, " ");
			writeln();
		}
		

		ref Vector opUnary(string op)()
		if (op == "+") // +v for Vectors 
		{
			return this; 
		}
		//Vector opUnary(string op)()
		//if (op == "-") // -v for Vectors 
		//{
		//	Vector!dim tmp;
		//	tmp.content[] = -content[];
		//	return tmp; 
		//}

		ref Vector opOpAssign(string op)(double rhs)
		if (op == "*" || op == "/") // u*=s; u/=s; for Vectors 
		{
			static if (op == "*")
				content[] *= rhs;
			else if (op == "/")
				content[] /= rhs;
			
			return this; 
		}
		ref Vector opOpAssign(string op)(Vector rhs)
		if (op == "+" || op == "-") // u+=v; u-=v; for Vectors 
		{
			static if (op == "+")
				content[] += rhs.content[];
			else if (op == "-")
				content[] -= rhs.content[];
			
			return this; 
		}
		ref Vector!(dim2,T) opOpAssign(uint dim2, string op)(Vector!(dim2,T) rhs)
		if ((dim2 < dim) && (op == "+" || op == "-")) // u+=v; u-=v; for Vectors 
		{
			static if (op == "+")
				content[0..dim2] += rhs.content[];
			else if (op == "-")
				content[0..dim2] -= rhs.content[];
			
			return this; 
		}
		ref Vector opOpAssign(string op)(T[dim] rhs)
		{
			return opAssign!op(Vector!(dim,T)(rhs));
		}


		T opBinary(string op)(in Vector rhs)
		if (op == "*") // scalar dot product for Vectors
		{
			T result = 0;
			foreach(uint i ; 0 .. dim)
				result += content[i] * rhs.content[i];
			return result;	
		}
		T opBinary(string op)(in T[dim] rhs)
		if (op == "*") // scalar dot product for arrays
		{
			return opBinary!op(Vector!(dim,T)(rhs));
		}

		Vector opBinary(string op)(in double rhs)
		if (op == "*" || op == "/") // u*s 
		{
			Vector!(dim,T) result = this;
			result.opOpAssign!op(rhs);
			return result;	
		}
		Vector opBinary(string op)(in Vector rhs)
		if (op == "+" || op == "-") // u+v for Vectors
		{
			Vector!(dim,T) result = this;
			result.opOpAssign!op(rhs);
			return result;	
		}
		
		Vector opBinary(uint dim2, string op)(in Vector!(dim2,T) rhs)
		if ((dim2 < dim) && (op == "+" || op == "-")) // u+v for Vectors
		{
			Vector!(dim,T) result = this;
			result.opOpAssign!op(rhs);
			return result;	
		}
		Vector opBinary(string op)(in T[dim] rhs)
		if (op != "*") // only the scalar dot product is excluded because it doesn't return a Vector but a T
		{
			return opBinary!op(Vector!(dim,T)(rhs));
		}

		Vector opBinary(string op)(in Vector rhs)
		if (op == "^^" && dim == 3)
		{
			Vector!(dim,T) result = [content[1]*rhs.content[2] - content[2]*rhs.content[1],
									 content[2]*rhs.content[0] - content[0]*rhs.content[2],
									 content[0]*rhs.content[1] - content[1]*rhs.content[0]];
			return result;
		}
		
		// this is implicitely defined 
		//bool opEquals(in Vector rhs)
		//{
		//	foreach (uint i ; 0 .. dim)
		//		if (content[i] != rhs.content[i])
		//			return false;
		//	return true;
		//}
		
		ref Vector normalize()
		{
			content[] /= this.length;
			return this;
		}
		
		@property double length() 
		{
			return sqrt(squaredLength);
		}

		@property double squaredLength() 
		{
			double ll = 0;
			foreach(x ; content)
				ll += x*x;
			return ll;
		}
		
		@property int dimension()
		{
			return dim;
		}
		
		@property T* ptr()
		{
			return &content[0];
		}
};


Vector!(3,T) eulerAngles(T = double)(Vector!(3,T) r)
{
	T rxy   = sqrt(r[0]^^2+r[1]^^2);
	T theta = atan2(rxy,r[2]);
	T phi   = atan2(r[1],r[0]); 
	return Vector!(3,T)([sqrt(rxy^^2+r[2]^^2),theta,phi]);
}

Vector!(3,T) eulerVector(T = double)(in double r, in double theta, in double phi)
{
	return Vector!(3,T)([r * cos(phi) * sin(theta),
						 r * sin(phi) * sin(theta),
						 r *            cos(theta)]);
}


Vector!(3,T) polarVector(T = double)(in double r, in double phi, in double theta)
{
	return Vector!(3,T)([r * sin(phi) * sin(theta),
						 r * cos(phi) * sin(theta),
						 r *            cos(theta)]);
}
Vector!(4,T) polarVector(T = double)(in double r, in double phi, in double theta, double omega)
{
	return Vector!(3,T)([r * sin(phi) * sin(theta) * sin(omega),
						 r * cos(phi) * sin(theta) * sin(omega),
						 r *            cos(theta) * sin(omega),
						 r *                       * cos(omega)]);
}

// n-dimenstional polar vector

Vector!(dim,T) polarVector(uint dim, T = double)(in T r, in T[] angles...)
{
	assert (angles.length == dim-1);
	
	Vector!(dim,T) result;
	for (uint i = 0; i < dim; ++i)
	{
		result[i] = r;
		if (i == 0)
		{
			for (uint j = 0; j < dim-1; ++j)
				result[i] *= sin(angles[j]);
		}
		else
		{
			for (uint j = i-1; j < dim-1; ++j)
			{
				if (j == i-1)
					result[i] *= cos(angles[j]);
				else
					result[i] *= sin(angles[j]);
			}
		}
	}
	return result;
}


//Vector!(dim, T) orthogonalVector(uint dim, T)(in Vector!(dim,T) v)
//{
//	auto a = v;
//	for (uint i = 0; i < dim; ++i)
//	{
//		a = a ^^ v;
//		auto alen = a*a;
//		auto vlen = v*v;
//		writeln(v, " ", a, " " , vlen , " " , alen);
//	}
//	return a;	
//}
//
//Vector!(dim, T) karthesian_polar(dim, T)(in Vector!(dim,T), in T[] polarAngles)
//{
//	assert(polarAngles.length == dim-1);
//	
//	
//	
//}



unittest
{
	Vector!3 u = [ 1,0,0 ];
	Vector!3 v = [ 0,1,0 ];
	Vector!3 w = [ 0,0,1 ];
	Vector!3 o = [ 1,1,1 ];

	Vector!3 u2 = u*2;
	u2 *= 2;
	assert(u2 == Vector!3([4,0,0]));
	
	assert(*u.ptr == 1);
	
	
	assert(u+v+w == o);
	
	assert(u^^v == w);
	assert(w^^u == v);
	assert(v^^w == u);

	assert((v^^w + w^^v).length < 0.0000001);
	assert((Vector!3([1,1,1]).normalize() - Vector!3([3,3,3]).normalize()).length < 0.00001);

	assert((u*v) == 0);
	assert((w*u) == 0);
	assert((v*w) == 0);

	assert((u*[0,1,0]) == 0);
	assert((w*[1,0,0]) == 0);
	assert((v*[0,0,1]) == 0);

	assert(Vector!3([0,1,0])*u == 0);
	assert(Vector!3([1,0,0])*w == 0);
	assert(Vector!3([0,0,1])*v == 0);
	
	assert((u*u) == 1);
	assert((v*v) == 1);
	assert((w*w) == 1);
	
	assert((polarVector(1,0,0) - w).length < 0.000001);
	assert((polarVector(1,0,PI) + w).length < 0.000001);
	assert((polarVector(1,0,PI/2)    - v).length < 0.000001);
	assert((polarVector(1,PI/2,PI/2) - u).length < 0.000001);
	
	
	writeln( polarVector(1,1,0) );
	writeln( polarVector!5(1.0, [1., 0., 1.0, 0.0]) );
	
	//Vector!3 a = [1,1,1];
	//writeln(a);
	//writeln(orthogonalVector!(3,double)(a));
	
	writeln("Vector unittest passed");
}


// Row major:
// 
//       col
//        | 
//row-- [ -------- ]
//      [ -------- ]
//      [ -------- ]
//
//
//
// Column major:
// 
//       col
//        | 
//row-- [ | | | | ]
//      [ | | | | ]
//      [ | | | | ]
//
enum
{
	rowMajor = 0,
	columnMajor = 1,
}
// orientation == 0 means row    major
// orientation == 1 means column major
struct Matrix(uint nrows, uint ncolumns, uint orientation = rowMajor, T = double)
{
	private: 
		static if (orientation == rowMajor)
		{
			Vector!(ncolumns,T)[nrows] contentRM;
		}
		static if (orientation == columnMajor)
		{
			Vector!(nrows,T)[ncolumns] contentCM;
		}
		
	public:	
		this(this)
		{
			static if (orientation == rowMajor)
			{
				contentRM = contentRM.dup ;
			}
			static if (orientation == columnMajor)
			{
				contentCM = contentCM.dup ;
			}
		}
		
		this(Matrix!(nrows, ncolumns, (orientation+1)%2, T) rhs)
		{
			foreach(uint row; 0 .. nrows)
				foreach(uint column; 0 .. ncolumns)
				{
					this.access(row,column) = rhs.access(row,column);
				}
		}
		
		this(in T[nrows*ncolumns] content)
		{
			foreach(uint row; 0 .. nrows)
				foreach(uint column; 0 .. ncolumns)
				{
					access(row,column) = content[column + ncolumns*row];
				}
		}

		static if (nrows == ncolumns)
		{
			ref Matrix identity()
			{
				static if (orientation == rowMajor)
				{
					foreach(uint i; 0 .. contentRM.length)
						contentRM[i].unity(i);
				}		
				static if (orientation == columnMajor)
				{
					foreach(uint i; 0 .. contentCM.length)
						contentCM[i].unity(i);
				}		
				return this;
			}
			
			
			ref Matrix rotation(Tangle)(uint axis1, uint axis2, Tangle angle)
			{
				identity();
				T cos_angle = cos(1.0*angle);
				T sin_angle = sin(1.0*angle);
				
				access(axis1,axis1) = cos_angle;	access(axis1,axis2) = -sin_angle;
				access(axis2,axis1) = sin_angle;	access(axis2,axis2) = cos_angle;   
				
				return this;
			}
			
			static if (nrows == 3)
			{
				ref Matrix rotation_axis(Vector!(3,T) axis)
				{
					T phi = sqrt(axis*axis);
					if (phi > 1e-10)
					{
						axis /= phi;
						T cp = cos(phi);
						T c = 1 - cp;
						T s = sin(phi);
						T n0 = axis[0];
						T n1 = axis[1];
						T n2 = axis[2];
						access(0,0) = n0*n0*c + cp;		access(0,1) = n0*n1*c - n2*s;	access(0,2) = n0*n2*c + n1*s;
						access(1,0) = n1*n0*c + n2*s;	access(1,1) = n1*n1*c + cp;		access(1,2) = n1*n2*c - n0*s;
						access(2,0) = n2*n0*c - n1*s;	access(2,1) = n2*n1*c + n0*s;	access(2,2) = n2*n2*c + cp;
					}
					else
					{
						identity();
					}
					return this;
				}
			}

			ref Matrix boost(Tbeta)(uint axis1, uint axis2, Tbeta beta)
			{
				identity();
				T gamma = 1/sqrt(1-beta^^2);
				
				access(axis1,axis1) = gamma;		access(axis1,axis2) = beta*gamma;
				access(axis2,axis1) = beta*gamma;	access(axis2,axis2) = gamma;   
				
				return this;
			}
			static if (nrows == 3)
			{
				ref Matrix boost_direction(Vector!(2,T) beta)
				{
					auto b = beta;
					T bb = b*b;
					T g  = 1/sqrt(1-bb);
					T g1 = g - 1;
					if (bb > 1e-10)
					{
						access(0,0) =  g;		access(0,1) = -g*b[0];			access(0,2) = -g*b[1];
						access(1,0) = -g*b[0];	access(1,1) = 1+g1*b[0]^^2/bb;	access(1,2) = g1*b[0]*b[1]/bb;
						access(2,0) = -g*b[1];	access(2,1) = g1*b[0]*b[1]/bb;	access(2,2) = 1+g1*b[1]^^2/bb;
					}
					else
					{
						identity();
					}
					return this;
				}
			}
			static if (nrows == 4)
			{
				ref Matrix boost_direction(Vector!(3,T) beta)
				{
					auto b = beta;
					T bb = b*b;
					T g  = 1/sqrt(1-bb);
					T g1 = g - 1;
					if (bb > 1e-10)
					{
						access(0,0) =  g;		access(0,1) = -g*b[0];			access(0,2) = -g*b[1];			access(0,3) = -g*b[2];
						access(1,0) = -g*b[0];	access(1,1) = 1+g1*b[0]^^2/bb;	access(1,2) = g1*b[0]*b[1]/bb;  access(1,3) = g1*b[0]*b[2]/bb;
						access(2,0) = -g*b[1];	access(2,1) = g1*b[0]*b[1]/bb;	access(2,2) = 1+g1*b[1]^^2/bb;  access(2,3) = g1*b[1]*b[2]/bb;
						access(3,0) = -g*b[2];	access(3,1) = g1*b[2]*b[0]/bb;	access(3,2) = g1*b[2]*b[1]/bb;  access(3,3) = 1+g1*b[2]^^2/bb;
					}
					else
					{
						identity();
					}
					return this;
				}
			}


		}



		@property uint rows()
		{
			return nrows;
		}

		@property uint columns()
		{
			return ncolumns;
		}
		
		@property T* ptr()
		{
			static if (orientation == rowMajor)
			{
				return &(contentRM[0][0]);
			}		
			static if (orientation == columnMajor)
			{
				return &(contentCM[0][0]);
			}
		}
		
		uint vectorOrientation() const
		{
			return orientation;
		}
		
		ref T access(uint row, uint column) 
		{
			//writeln("orientation = ", orientation , "    nrows = ", nrows, "    ncolumns = ", ncolumns, "   row = ", row , "   column = ", column);
			assert(row < nrows && column < ncolumns);
			static if (orientation == rowMajor)
			{
				return contentRM[row][column];
			}
			static if (orientation == columnMajor)
			{
				return contentCM[column][row];
			}
		}
		const T access(uint row, uint column) 
		{
			return access(row,column);
		}
		
		ref T opIndex(uint row, uint column)
		{
			return access(row,column);
		}
		
		// compare to a different orientation
		bool opEquals(Matrix!(nrows,ncolumns,(orientation+1)%2,T) rhs)
		{
			//writeln("opEquals: orientation = " , orientation, "   rhs.orientation = " , rhs.vectorOrientation());
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns)
				{
					if (this.access(row,column) != rhs.access(row,column))
						return false;
				}
			return true;
		}
		bool opEquals(in Matrix!(nrows,ncolumns,orientation,T) rhs)
		{
			static if (orientation == rowMajor)
			{
				return contentRM == rhs.contentRM;
			}
			static if (orientation == columnMajor)
			{
				return contentCM == rhs.contentCM;
			}
		}

		Matrix!(ncolumns, nrows, orientation, T) transpose()
		{
			Matrix!(ncolumns, nrows, orientation, T) result;
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns)
					result.access(column,row) = this.access(row,column);
			return result;
		}
		Matrix!(ncolumns, nrows, (orientation+1)%2, T) transpose2()
		{
			Matrix!(ncolumns, nrows, (orientation+1)%2, T) result;
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns)
					result.access(column,row) = this.access(row,column);
			return result;
		}


		ref Matrix opOpAssign(string op)(T rhs)
		if (op == "*" || op == "/")
		{
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns)
				{
					static if (op == "*")
					{
						this.access(row,column) *= rhs;
					}
					else if (op == "/")
					{
						this.access(row,column) /= rhs;
					}
				}
			return this;
		}
		
		ref Matrix opOpAssign(string op)(Matrix rhs)
		if (op == "+" || op == "-")
		{
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns)
				{
					static if (op == "+")
					{
						this.access(row,column) += rhs.access(row,column);
					}
					else if (op == "-")
					{
						this.access(row,column) -= rhs.access(row,column);
					}
				}
			return this;
		}
		ref Matrix opOpAssign(string op)(Matrix!(nrows, ncolumns, (orientation+1)%2, T) rhs)
		if (op == "+" || op == "-")
		{
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns)
				{
					static if (op == "+")
					{
						this.access(row,column) += rhs.access(row,column);
					}
					else if (op == "-")
					{
						this.access(row,column) -= rhs.access(row,column);
					}
				}
			return this;
		}

		Matrix!(nrows,ncolumns2,orientation,T) opBinary(string op, uint ncolumns2,uint orientation2)(Matrix!(ncolumns,ncolumns2,orientation2,T) rhs)
		if (op == "*")
		{
			Matrix!(nrows,ncolumns2,orientation,T) result;
			foreach (uint row; 0 .. nrows)
				foreach (uint column; 0 .. ncolumns2)
				{
					result.access(row,column) = 0;
					foreach(uint i; 0 .. ncolumns)
					{
						result.access(row,column) += this.access(row,i)*rhs.access(i,column);
					}
				}
			return result;
		}
		Vector!(nrows,T) opBinary(string op)(Vector!(ncolumns,T) rhs)
		if (op == "*")
		{
			Vector!(nrows,T) result;
			static if (orientation == rowMajor)
			{
				foreach ( int i ; 0 .. nrows)
					result[i] = rhs * this.contentRM[i];
			}
			else if (orientation == columnMajor)
			{
				foreach (int i ; 0 .. ncolumns)
					result += this.contentCM[i] * rhs[i];
			}
			return result;
		}

		Matrix opBinary(string op)(in double rhs)
		if (op == "*" || op == "/") // u*s 
		{
			Matrix result = this;
			result.opOpAssign!op(rhs);
			return result;	
		}
		Matrix opBinary(string op)(Matrix rhs)
		if (op == "+" || op == "-") // u+v for Matrix
		{
			Matrix result = this;
			result.opOpAssign!op(rhs);
			return result;	
		}
		Matrix opBinary(string op)(Matrix!(nrows,ncolumns,(orientation+1)%2,T) rhs)
		if (op == "+" || op == "-") // u+v for Matrix of different orientation
		{
			Matrix result = this;
			result.opOpAssign!op(rhs);
			return result;	
		}
		
		void print()
		{
			foreach (uint row; 0 .. nrows)
			{
				foreach (uint column; 0 .. ncolumns)
					writef("%6.3f",access(row,column));
				writeln();	
			}
		}


		// Vector-wise access through row vectors
		// the single index operator always returns the row-vector
		static if (orientation == rowMajor)
		{
			ref Vector!(ncolumns,T) rowVector(uint row)
			{
				assert(row < nrows);
				return contentRM[row];
			}
			ref Vector!(ncolumns,T) opIndex(uint row)
			{
				assert(row < nrows);
				return contentRM[row];
			}
			const ref Vector!(ncolumns,T) opIndex(uint row)
			{
				return opIndex(row);
			}
		}
		static if (orientation == columnMajor)
		{
			Vector!(ncolumns,T) rowVector(uint row)
			{
				assert(row < nrows);
				Vector!(ncolumns,T) result;
				foreach(uint column; 0 .. ncolumns)
					result[column] = access(row, column);
				return result;
			}
			Vector!(ncolumns,T) opIndex(uint row)
			{
				assert(row < nrows);
				return rowVector(row);
			}
		}

		// Vector-wise access through column vectors
		static if (orientation == columnMajor)
		{
			ref Vector!(nrows,T) columnVector(uint column)
			{
				assert(column < ncolumns);
				return contentCM[column];
			}
			const ref Vector!(nrows,T) columnVector(uint column)
			{
				return columnVector(column);
			}
		}
		static if (orientation == rowMajor)
		Vector!(nrows,T) columnVector(uint column)
		{
			assert(column < ncolumns);
			Vector!(nrows,T) result;
			foreach(uint row; 0 .. nrows)
				result[row] = access(row, column);
			return result;
		}
		
}

Matrix!(N,N) kardanMatrix(uint N, uint orientation = rowMajor, T = double)(in T[] phi)
{
	assert (phi.length >= N*(N-1)/2);
	
	Matrix!(N,N) m;
	m.identity();
	uint idx = 0;
	Matrix!(N,N) d;
	
	for (uint i = 0; i < N-1; ++i)
		for (uint j = i+1; j < N; ++j)
			m = m * d.rotation(i,j, phi[idx++]);
			
	return m;
}

Matrix!(N,N) kardanMatrixInverse(uint N, uint orientation = rowMajor, T = double)(in T[]  phi) 
{
	assert (phi.length >= N*(N-1)/2);
	
	Matrix!(N,N) m;
	m.identity();
	uint idx = 0;
	Matrix!(N,N) d;
	
	for (uint i = 0; i < N-1; ++i)
		for (uint j = i+1; j < N; ++j)
			m = d.rotation(i,j, -phi[idx++]) * m ;
			
	return m;
}


T[N*(N-1)/2] kardanAngles(uint N, uint orientation = rowMajor, T = double)(Matrix!(N,N,orientation,T) m)
{
	T[N*(N-1)/2] phi;
	
	int idx = N*(N-1)/2 - 1;
	Matrix!(N,N,orientation,T) d;	
	
	for (int i = N-2; i >= 0; --i)
		for (int j = N-1; j >= i+1; --j)
		{
			phi[idx] =  atan2(m[j][i], m[j][j]);
			m = m * d.rotation(j,i, phi[idx--]) ;
		}
	return phi;
}


unittest
{
	Matrix!(3,4,rowMajor) mRM = Matrix!(3,4,rowMajor)([ 1, 2, 3, 4,
														5, 6, 7, 8,
														9,10,11,12]);
	Matrix!(3,4,columnMajor) mCM = Matrix!(3,4,columnMajor)([ 1, 2, 3, 4,
															  5, 6, 7, 8,
															  9,10,11,12]);

	auto m = mRM;
	auto t = m.transpose();
	//auto mmm = m.multiply(t);
	//mmm.print();
	
	assert(mRM == mCM);
	assert(mRM == m);

	m[0][0] = 10;
	assert(mRM != m);
	
	Matrix!(3,4,columnMajor) mCM2 = mRM;
	mCM2 += mRM;
	//mCM *= 2;
	auto mm = mCM2 - mCM*2;
	mm.print();
	
	
	Matrix!(2,4,columnMajor) M1 = ([ 1, 1, 1, 1,
								  2, 2, 2, 2,]);
	Matrix!(4,2,rowMajor) M2 = ([ 1, 2, 
								 1, 2,
								 1, 2, 
								 1, 2,]);
	auto M1M2 = M1*M2;
	M1M2.identity();
	
	writeln("");
	M1.print();
	writeln("");
	M2.print();
	writeln("");
	M1M2.print();
	writeln("");
	
	Matrix!(3,3) rot1; rot1.rotation(1,2,0.5);
	Matrix!(3,3) rot2; rot2.rotation(1,2,0.5);
	Matrix!(3,3) rot3; rot3.rotation(1,2,1.0);
	//Matrix!(3,3) rot4 = rot1*rot2;
	
	rot3.rotation_axis(Vector!3([0,0,1]));
	//rot4.rotation(0,1,1);
	writeln("---");
	rot3.print();
	writeln("---");
	//rot4.print();
	writeln("---");
		
	double[6] kardan = [0.1,0.2,0.3,0.4,0.5,0.6];
	auto km = kardanMatrix!(4)(kardan) * kardanMatrixInverse!(4)(kardan);
	km.print();
	
	
	auto kardan2 = kardanAngles!(4)( kardanMatrix!(4)(kardan) );
	auto kardan3 = kardanAngles!(4)( kardanMatrix!(4)(kardan2) );
	
	
	writeln(kardan2);
	writeln(kardan3);

	writeln("Matrix unittest passed");
}

