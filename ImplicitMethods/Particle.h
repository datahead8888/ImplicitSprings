#pragma once

const int DIMENSION = 3;	//Number of dimensions - 3 means 3D

//This struct represents one particle
//A particle corresponds to both a vertex in the edge system and to a vertex used in rendering triangles
struct Particle
{
	Particle()
	{
		position[DIMENSION] = 1;
		normal[DIMENSION] = 1;
		color[0] = 0;
		color[1] = 0;
		color[2] = 1;
		color[3] = 1;
	}

	float position [DIMENSION+1];	//Position vector in 3d space
	float normal [DIMENSION+1];		//Normal vector in 3d space (used in lighting)
	float color [DIMENSION+1];		//Color of point
	float velocity [DIMENSION+1];	//Velocity vector in 3d space
	int triangleCount;				//Number of triangles in which the vertex for this particle is contained (helps calculate normals)
};
