#pragma once

//This class represents one edge, which connects two particles with a spring
struct Edge
{
	int start;					//Start index for edge (index is to a particle entry)
	int end;					//End index for edge
	int restLengthIndex;		//Index to rest length for this edge (because straight and diagonal lengths may differ)
	int isStartUnconstrained;	//If the 1st particle in the edge is contrained this will be 0; otherwise 1.
	int isEndUnconstrained;		//Same as above but for the 2nd particle
};
