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


module nucd.nucleus;

import std.stdio;
import std.json;
import std.conv;
import std.algorithm;
import std.array;
import std.math;

double radius(int A)
{
	return 1.25 * A^^(1./3.);
}


class NuclearData
{
	struct State
	{
		struct Decay
		{
			enum Mode {M = 'M', E = 'E'} 
			int    destIndex;          // index of the destination state
			Mode   multMode;           // electric or magnetic transition
			int    multLambda;         // order of multipolarity (0: monopole, 1: diploe, 2: quadrupole, ...)
			double relativeIntensity;  // relative intensity of the decay
			double strength;           // decay strength in 1/ps
		}
		string  label;                 // a copy of the the label
		double  E;                     // excitation energy [keV]
		int     twiceSpin = -1;        // two times the spin (to make it integer)
		int     parity;                // valid values are +1 or -1
		double  tau = float.infinity;  // lifetime [ps]
		Decay[] decays;                // list of indices to other states
	}                                  
	string      isotope;               // name of the isotope
	int         N = -1;                // number of neutrons
	int         Z = -1;                // number of protons
	int         A = -1;                // mass number
	State[]     states;                // list of all states
	double 		mass;			       // mass of the nucleus
	int[string] stateIndex;            // maps state names to array indices
	
	void print()
	{
		foreach(state; states)
		{
			writeln(state.label);
			foreach(decay; state.decays)
				writeln("   -> "
					, states[decay.destIndex].label, " "
					, decay.multMode, decay.multLambda, " "
					, decay.strength
					);
		}
	}
	
	this(const JSONValue ncl)
	{
		if (ncl.type != JSONType.object)
		{
			writeln("error: toplevel json is no object");
			return;
		}

		if ("isotope" in ncl.object) isotope = ncl.object["isotope"].str;
		if ("Z" in ncl.object)       Z = to!int(ncl.object["Z"].integer);
		if ("N" in ncl.object)       N = to!int(ncl.object["N"].integer);
		if ("A" in ncl.object)       A = to!int(ncl.object["A"].integer);
		if ("mass" in ncl.object)       A = to!int(ncl.object["mass"].floating);

		if ("states" in ncl.object && ncl.object["states"].type == JSONType.array)
		{
			foreach(state; ncl.object["states"].array.map!(a => a.object))
			{ // state is of type JSONValue[string]
				if ("label" in state) 
				{	
					string label = state["label"].str;
					//writeln("adding state ", label);
					stateIndex[label] = to!int(states.length);
					states ~= NuclearData.State();
					states[$-1].label = label;

					// extract spin and parity from label
					string labelSpin = label[0 .. $-2];
					char labelParity = label[$-2];
					if      (labelParity == '+') 
					{
						states[$-1].parity = +1;
					}
					else if (labelParity == '-') 
					{
						states[$-1].parity = -1;
					}
					else 
					{
						writeln("error: cannot extract parity from label");
					}
					string[] splittedLabel = split(labelSpin,'/');
					if (splittedLabel.length == 1)
					{
						states[$-1].twiceSpin = 2*to!int(splittedLabel[0]);
					}
					else if (splittedLabel.length == 2) // assumes that the label has the form x/2
					{
						states[$-1].twiceSpin = to!int(splittedLabel[0]);
					}
				}
				if ("E"  in state) // get energy
				{
					if (state["E"].type == JSONType.float_)
						states[$-1].E = state["E"].floating;
					if (state["E"].type == JSONType.integer)
						states[$-1].E = state["E"].integer;
				}
				if ("T1/2"  in state) // get half-life
				{
					if (state["T1/2"].type == JSONType.float_)
						states[$-1].tau = state["T1/2"].floating/log(2);
					if (state["T1/2"].type == JSONType.integer)
						states[$-1].tau = state["T1/2"].integer/log(2);
				}
				// get decays
				if ("decays" in state) // if the state decays
				{
					foreach(decay; state["decays"].array.map!(a => a.object))
					{ // decay is of type JSONValue[string]
						if ("dest" in decay) 
						{
							states[$-1].decays ~= NuclearData.State.Decay();
							states[$-1].decays[$-1].destIndex = stateIndex[decay["dest"].str];
						}
						if ("mult" in decay)
						{
							string multString = decay["mult"].str;
							string multLambdaString = multString[1 .. $];
							if (multString[0] == 'E')
							{
								states[$-1].decays[$-1].multMode = NuclearData.State.Decay.Mode.E;
							}
							if (multString[0] == 'M')
							{
								states[$-1].decays[$-1].multMode = NuclearData.State.Decay.Mode.M;
							}
							states[$-1].decays[$-1].multLambda = to!int(multLambdaString);
						}
						if ("intensity" in decay)
						{
							auto intensityJSON = decay["intensity"];
							switch(intensityJSON.type)
							{
								case JSONType.integer:
									states[$-1].decays[$-1].relativeIntensity = intensityJSON.integer;
								break;
								case JSONType.float_:
									states[$-1].decays[$-1].relativeIntensity = intensityJSON.floating;
								break;
								default:
									writeln("error: decay intensity is neither integer nor float");
							}
						}
						//writeln(states[$-1].decays);
					}
				}
				else if (states[$-1].tau < float.infinity)
				{
					writeln("error: states with finite lifetime need decays");
				}
				double sumIntensity = 0;
				foreach(decay; states[$-1].decays)
					sumIntensity += decay.relativeIntensity;
				foreach(ref decay; states[$-1].decays)
					decay.strength = decay.relativeIntensity / ( sumIntensity * states[$-1].tau);
				//writeln(states[$-1].decays);	
			}
		}
	}
}

struct Nucleus
{
	NuclearData levelScheme;  // describes the nuclear structure
	int         levelIdx;     // describes the current state of the nucleus
	int         twiceM;       // two times the magnetic substate (orientation)
}

unittest
{
	//auto n = new NuclearLevelScheme;
	//n.levels ~= NuclearLevelScheme.Level(null,1,100.0,10.0);
	
	// todo: save n to into a json file...
	
	JSONValue[string]  obj;
	obj["name"]       = JSONValue();
	obj["name"].str   = "Charlie";
	obj["id"]         = JSONValue();
	obj["id"].integer = 3;

	JSONValue          newPerson;
	newPerson.object  = obj;
	

	auto ncl = parseJSON
	(`{
	"isotope": "80Kr",
	"Z": 80,
	"N": 36,
	"mass": 79.916378965,
	"states" : [
		{	"label": "0+1"
			,"E": 0
		}
		,{	"label": "2+1"
			,"E": 616.6
			,"T1/2": 8.6
			,"decays": [
				{ "dest":"0+1", "mult":"E2", "intensity":30 }
			]
		}
		,{	"label": "2+2"
			,"E": 1256.24
			,"T1/2": 7.6
			,"decays": [
				 { "dest":"2+1", "mult":"E2", "intensity":90 }
				,{ "dest":"2+1", "mult":"M1", "intensity":10 }
				,{ "dest":"0+1", "mult":"E2", "intensity":30 }
			]
		}
	]
	}`);
	//writeln(ncl.toPrettyString());
	NuclearData nd = new NuclearData(ncl);
	//nd.print();
}

