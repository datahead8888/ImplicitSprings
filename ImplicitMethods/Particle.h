#pragma once

const int DIMENSION = 3;	//Number of dimensions - 3 means 3D

//This struct represents one particle
//A particle corresponds to both a vertex in the edge system and to a vertex used in rendering triangles
struct Particle
{
	double position [DIMENSION];	//Position vector in 3d space
	double velocity [DIMENSION];	//Velocity vector in 3d space	
	double normal [DIMENSION];		//Normal vector in 3d space (used in lighting)
	int triangleCount;				//Number of triangles in which the vertex for this particle is contained (helps calculate normals)
};
