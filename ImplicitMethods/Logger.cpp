#include "Logger.h"
#include <list>
#include <iostream>
#include <iomanip>

using namespace std;

extern const int DIMENSION;

//Constructor for Logger class
Logger::Logger()
{
	loggingLevel = MEDIUM;
	isLogging = false;
	printConstraintDeltaV = false;
}

//This method prints the A matrix fully - ie it prints all zeroe's in the matrix
//Note: the current implementation allocates a massive doubly allocated array to stores all the zeroes, taking up quite a bit of memory each time its called for large numbers of particles
//Also, the number of couts uses is excessive when a large number of particles are present
//Thus this method is only useful for evaluating systems with a small number of particles
//Parameters:
//ADiagonalEntries - the diagonal entries in A
//numParticles - the number of particles in the particle system
//ANonDiagonalEntries - the non-zero non-diagonal entries in A
//numEdges - the number of edges in the particles system
//edgeList - the array constaining a list of all edges
//name - a string representing the name of the matrix so that it can be printed (currently A)
//requestLevel - the logging level at which to print this
void Logger::printAFull(double * ADiagonalEntries, int numParticles, double * ANonDiagonalEntries, int numEdges, Edge * edgeList, char * name, LoggingLevel requestLevel)
{
	if (loggingLevel < requestLevel)
	{
		return;
	}
	//Note: this will be slow and take a lot of memory, so it should only be used while in debug mode and for small matrices
	//We may want to pursue an stl list solution as a better solution long term, but it also is only for debug mode
	double ** bigA = new double * [numParticles * DIMENSION];	//One matrix representing A with all zeroes present
	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		bigA[i] = new double[numParticles * DIMENSION];
		for (int j = 0; j < numParticles * DIMENSION; j++)
		{
			bigA[i][j] = 0;
		}
	}

	for (int e = 0; e < numEdges; e++)
	{
		int part1 = edgeList[e].start;
		int part2 = edgeList[e].end;
		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				bigA[part1 * DIMENSION + i][part2 * DIMENSION + j] = ANonDiagonalEntries[e * DIMENSION * DIMENSION + i * DIMENSION + j]  * edgeList[e].isStartUnconstrained ;
				bigA[part2 * DIMENSION + i][part1 * DIMENSION + j] = ANonDiagonalEntries[e * DIMENSION * DIMENSION + i * DIMENSION + j]  * edgeList[e].isEndUnconstrained ;
			}
		}
	}

	for (int k = 0; k < numParticles; k++)
	{
		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				bigA[k * DIMENSION + i][k * DIMENSION + j] = ADiagonalEntries[k * DIMENSION * DIMENSION + i * DIMENSION + j];
			}
		}
	}

	cout << name << ":" << endl;
	
	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		for (int j = 0; j < numParticles * DIMENSION; j++)
		{
			cout << setw(7) << bigA[i][j] << " ";
    
		}
		cout << endl;
	}
	
	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		delete [] bigA[i];
	}

	delete [] bigA;


	
		
}

//This method prints one vector
//Parameters:
//vector - the vector to be printed
//size - the number of elements in the vector
//message - the text to be printed for the vector
//requestLevel - the logging level at which to print this
void Logger::printVector(double * vector, int size, char * message, LoggingLevel requestLevel)
{
	if (loggingLevel < requestLevel)
	{
		return;
	}
	cout << message << ":" << endl;
	for (int i = 0; i < size; i++)
	{
		cout << setw(8) << vector[i] << endl;
	}
}

//This method prints the contents of the edge array
//Parameters:
//edges - the array containing the edges
//numEdges - the number of edges
//requestLevel - the logging level at which to print this
void Logger::printEdges(Edge * edges, int numEdges, LoggingLevel requestLevel)
{
	if (loggingLevel < requestLevel)
	{
		return;
	}
	for (int i = 0; i < numEdges; i++)
	{
		Edge * currentEdge = &edges[i];
		cout << "Edge " << i << ": ";
		cout << currentEdge -> start << " " << currentEdge -> end << " " << currentEdge -> restLengthIndex << " " << currentEdge -> isStartUnconstrained << " " << currentEdge -> isEndUnconstrained << endl;
	}

}

//This method prints the contents of the particles list
//It prints positions and velocities
//Parameters:
//particles - the array of particles
//numParticles - the number of particles
//requestLevel - the logging level at which to print this
void Logger::printParticles(Particle * particles, int numParticles, LoggingLevel requestLevel)
{
	if (loggingLevel < requestLevel)
	{
		return;
	}
	cout << "Positions:" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		Particle * particle = &particles[i];
		cout << "Particle " << i << ": ";
		cout << particle -> position[0] << " " << particle -> position[1] << " " << particle -> position[2] << endl;
	}

	cout << "Velocities:" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		Particle * particle = &particles[i];
		cout << "Particle " << i << ": ";
		cout << particle -> velocity[0] << " " << particle -> velocity[1] << " " << particle -> velocity[2] << endl;
	}
}

//This method prints the normals for the particles
//Parameters:
//particles - the array of particles
//numParticles - the number of particles
//requestLevel - the logging level at which to print this
void Logger::printNormals(Particle * particles, int numParticles, LoggingLevel requestLevel)
{
	if (loggingLevel < requestLevel)
	{
		return;
	}
	cout << "Normals:" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		Particle * particle = &particles[i];
		cout << "Particle " << i << ": ";
		cout << particle -> normal[0] << " " << particle -> normal[1] << " " << particle -> normal[2] << endl;
	}
	
}

//This method prints the mass matrix.  It prints all zeroes in the matrix
//Parameters:
//massMatrix - the mass matrix
//numParticles - number of particles in the particle system
//requestLevel - level at which to log this
void Logger::printMassMatrix(double * massMatrix, int numParticles, LoggingLevel requestLevel)
{
	if (loggingLevel < requestLevel)
	{
		return;
	}
	cout << "Mass Matrix:" << endl;
	for (int i = 0 ; i < numParticles * DIMENSION; i++)
	{
		for (int j = 0; j < numParticles * DIMENSION; j++)
		{
			if (i == j)
			{
				cout << setw(7) << massMatrix[i] << " ";
			}
			else
			{
				cout << setw(7) << "0 ";
			}
		}
		cout << endl;
	}
}