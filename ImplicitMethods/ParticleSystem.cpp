#include <gl/glew.h>
#include "ParticleSystem.h"
#include <gl/glut.h>
#include <gl/GLU.h>
#include <cmath>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include "targa.h"

using namespace std;


////////////////////////////////////////
//Main ParticleSystem class implementation
////It builds a grid of particles (represented by Particle objects)
//These are then held together by springs using hooke's laws.  One edge corresponds to each spring.
//Hooke's law is used to represent the springs, with both a spring component and a damping comopnent
//Implicit methods are used for the integration.  This requires a solving system but allows much larger spring constants without instability
//and potentially allows more efficient implementation
//A gravity component is also present
////////////////////////////////////////

//Based on:
//http://graphics.snu.ac.kr/~kjchoi/publication/cloth.pdf - Explicit and implicit formulas for hooke's law - page 3
//http://www.amath.unc.edu/Faculty/mucha/Reprints/SCAclothcontrolpreprint.pdf - Formulas for A and b in Ax = b system - page 5
//http://en.wikipedia.org/wiki/Conjugate_gradient - Algorithm for conjugate gradient (solves Ax = b)

const double epsilon = 1e-12;	//Used to check approximate equality to 0
extern const int DIMENSION;		//DIMENSION of system (3 for 3D)

//Constructor - initializes particles and settings
ParticleSystem::ParticleSystem(Logger * logger)
{
	this -> logger = logger;
	dimensionSquared = DIMENSION * DIMENSION;

	//Declare Constants
	kd = 15;					//Damping constant
	ks = 500;					//Spring constant - 1000 seems to give reasonable stiffness, but much larger values work, too
	earthGravityValue = 9.8;	//meters per second for earth gravity force
	//Set number of rows and columns for grid here
	rows = 40;
	cols = 40;
	sheetWidth = 5;
	sheetHeight = 5;
	double density = 3;	//Density of object - used to find total object mass

	//Constraints - these particles will not move due to any force
	numConstraints = 4;
	constraintParticles = new int[numConstraints];
	constraintParticles[0] = 0;
	constraintParticles[1] = cols - 1;
	constraintParticles[2] = rows * cols - cols;
	constraintParticles[3] = rows * cols - 1;
			
	//Calculate all secondary variables
	double objectMass = density * sheetWidth * sheetHeight;		//Total object mass - density times area (instead of volume since 2D)
	numParticles = rows * cols;
	double perParticleMass = objectMass / (numParticles);		//Mass of one particle - we are evenly distributing object mass across each particle
	
	//Normally the number of edges horizontally / vertically is one less than the number of vertices in that dimension, unless you have 1 vertex.
	restLengthValues[0] = (cols > 1)?sheetWidth / (cols - 1):0;		//Horizontal rest length of particles - for hooke's law
	restLengthValues[1] = (rows > 1)?sheetHeight / (rows - 1):0;    //Vertical rest length of particles
	restLengthValues[2] = sqrt(restLengthValues[0] * restLengthValues[0] + restLengthValues[1] * restLengthValues[1]);  //Use pythagorean theorem to calculate diagonal restLength
	
	//Compute half of width/height - allows grid to be relatively centered while rendering
	halfWidth = sheetWidth / 2;
	halfHeight = sheetHeight / 2;
	
	particles = new Particle [numParticles];
	
	reset();  //Set up the particle positions / velocities

	//Create edge objects
	numEdges = (rows - 1) * (cols - 1) * 4 + (rows - 1) + (cols - 1); //Calculated from logic of loop below
	edgeList = new Edge [numEdges];
	for (int i = 0; i < numEdges; i++)
	{
		edgeList[i].isStartUnconstrained = 1;
		edgeList[i].isEndUnconstrained = 1;
	}
	
	//Set up edges
	int edgeCounter = 0;	//Index to edge currently being processed
	int i = 0, j = 0;		//For loop indexes
	for (i = 0; i < rows; i++) //row
	{
		for (j = 0; j < cols; j++) //column
		{
			//Vertical connection - going "down"
			if (i <= rows - 2)
			{
				edgeList[edgeCounter].start = i * cols + j;
				edgeList[edgeCounter].end = (i + 1) * cols + j;
				edgeList[edgeCounter].restLengthIndex = 1;
				edgeCounter++;
			}
			//Diagonal connection - going "downright"
			if (i <= rows - 2 && j <= cols - 2)
			{
				edgeList[edgeCounter].start = i * cols + j;
				edgeList[edgeCounter].end = (i + 1) * cols + j + 1;
				edgeList[edgeCounter].restLengthIndex = 2;
				edgeCounter++;
			}
			//Horizontal connection - going "right"
			if (j <= cols - 2)
			{
				edgeList[edgeCounter].start = i * cols + j;
				edgeList[edgeCounter].end = i * cols + j + 1;
				edgeList[edgeCounter].restLengthIndex = 0;
				edgeCounter++;
			}
			
			//Diagonal connection - going "upright" FROM the particle that is "one down"
			if (j <= cols - 2 && i <= rows - 2)
			{
				edgeList[edgeCounter].start = (i + 1) * cols + j;  //One down from current position
				edgeList[edgeCounter].end = i * cols + j + 1;    //One right from current position (up right from one down)
				edgeList[edgeCounter].restLengthIndex = 2;
				edgeCounter++;
			}
			

		}

		

	}

	#ifdef DEBUGGING
	if(numEdges != edgeCounter)
	{
		cout << "Difference between numEdges and edgeCounter" << endl;
		cout << "numEdges: " << numEdges << endl;
		cout << "edgeCounter: " << edgeCounter << endl;
	}
	assert(numEdges == edgeCounter);  //If number edges processed does not match calculated number of edges used for size of list - there is a problem
	#endif

	//Mark constrained particles in the edge list with a factor of 0
	//Note this means the matrix A is not actually symmetric
	for (int e = 0; e < numEdges; e++)
	{
		int part1 = edgeList[e].start; //Particle 1
		int part2 = edgeList[e].end;   //Particle 2

		for (int k = 0; k < numConstraints; k++)
		{
			if (part1 == constraintParticles[k])
			{
				edgeList[e].isStartUnconstrained = 0;
			}
			if (part2 == constraintParticles[k])
			{
				edgeList[e].isEndUnconstrained = 0;
			}

		}
	}

	//Allocate memory for all types of dynamic arrays
	ADiagonalEntries = new double [numParticles * dimensionSquared]; //For 3D - sparse matrix representation - 9 entries (3X3 sum of jacobians) per particle
	ANonDiagonalEntries = new double [numEdges * dimensionSquared];  //For 3D - sparse matrix representation - 9 entries (one jacobian) per edge
	
	massMatrix = new double [numParticles * DIMENSION];
	b = new double[DIMENSION * numParticles];

	currentForce = new double[numParticles * DIMENSION];
	zeroVector = new double [DIMENSION * numParticles];
	deltaV = NULL;

	nextX = new double[DIMENSION*numParticles];
	r = new double [numParticles * DIMENSION];
	nextP = new double[DIMENSION*numParticles];
	temp = new double[DIMENSION*numParticles];
	nextR = new double[DIMENSION*numParticles];
	p = new double[numParticles*DIMENSION];

	//Allocate buffers to hold portions of the screen grabbed / swapped for video generation
	//We allocate enough space to work with the entire screen because the window size must be less than the screen
	int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	screenBuffer = new uint8_t[screenWidth * screenHeight * 4];
	screenRowTemp = new uint8_t [screenWidth * 4];
	
	//Create mass matrix.  Represent it with a vector, since it's a diagonal matrix.
	//Referred to http://en.wikipedia.org/wiki/Mass_matrix for overview of mass matrices
	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		massMatrix[i] = perParticleMass;				
	}
	
	for (int i = 0; i < DIMENSION * numParticles; i++)
	{
		zeroVector[i] = 0;
	}
	
	//Set values for boolean variables
	showInfoText = true;
	isAnimating = true;
	renderMode = 1;
	renderToImage = false;
	frameNumber = 1;

	for (int i = 0; i < TEXT_SIZE; i++)
	{
		text[i] = '\0';
	}

	windowWidth = windowHeight = 0;

	system("del images\\*.tga");  //DOS command to remove all tga image files from previous executions of this program

	#ifdef DEBUGGING
	logger -> loggingLevel = logger -> LIGHT;
	if (logger -> isLogging)
	{
		cout << "Version 49 - video generation by set of images" << endl;
		logger -> printParticles(particles, numParticles, logger -> MEDIUM);
		logger -> printEdges(edgeList, numEdges, logger -> MEDIUM);
	}
	#endif

	////////////////////////////LIGHTING////////////////////////////////////////////////
	//If lighting is calculated in eye space, the eye position basically is the origin - use this for the default
	eyePos[0] = 0.0f;
	eyePos[1] = 0.0f;
	eyePos[2] = -5.0f;
 
	lightAmbient[0] = 0.05; 
	lightAmbient[1] = 0.05;
	lightAmbient[2] = 0.05;

	lightFullAmbient[0] = 0.9; //Lots of ambient - useful for debugging in wireframe mode
	lightFullAmbient[1] = 0.9;
	lightFullAmbient[2] = 0.9;
	lightFullAmbient[3] = 1;

	lightDiffuse[0] = 1.0;
	lightDiffuse[1] = 1.0;
	lightDiffuse[2] = 1.0;
	lightDiffuse[3] = 1;

	//Turning off specular for now unless it's needed (and verified as good looking)
	//lightSpecular[0] = 1;
	//lightSpecular[1] = 1;
	//lightSpecular[2] = 1;
	lightSpecular[0] = 0;
	lightSpecular[1] = 0;
	lightSpecular[2] = 0;
	lightSpecular[3] = 1;

	lightPosition[0] = 0;
	lightPosition[1] = 0;
	lightPosition[2] = 1;
	lightPosition[3] = 1;

	matAmbient[0] = 1.0;
	matAmbient[1] = 1.0;
	matAmbient[2] = 1.0;
	matAmbient[3] = 1;

	matAmbientBack[0] = 1.0;
	matAmbientBack[1] = 1.0;
	matAmbientBack[2] = 1.0;
	matAmbientBack[3] = 1;

	matDiffuse[0] = 0;
	matDiffuse[1] = 0.7;
	matDiffuse[2] = 0;
	matDiffuse[3] = 1;

	matDiffuseBack[0] = 0;
	matDiffuseBack[1] = 0;
	matDiffuseBack[2] = 0.7;
	matDiffuseBack[3] = 1;

	matSpecular[0] = 0.9;
	matSpecular[1] = 0.9;
	matSpecular[2] = 0.9;
	matSpecular[3] = 1;

	

	matShininess[0] = 10000;

	ambientMode = false;

}

//Destructor - free all memory for dynamically allocated arrays
ParticleSystem::~ParticleSystem()
{
	delete [] edgeList;

	delete [] particles;
	delete [] constraintParticles;
	delete [] ADiagonalEntries;
	delete [] ANonDiagonalEntries;
	delete [] b;
	delete [] massMatrix;
	delete [] zeroVector;
	delete [] currentForce;
	delete [] nextX;
	delete [] r;
	delete [] nextP;
	delete [] temp;
	delete [] nextR;
	delete [] p;
	delete [] screenBuffer;
	delete [] screenRowTemp;
}

//This method resets the simulation by returning all particles to their original positions and velocities
//It is also used by the constructor to initialize all particles' positions and velocities, forming a sheet (grid) of cloth
//The particles essentially correspond to vertices in a graph
//This method does NOT change view control, changes to the gravity setting, or changes to the rest length setting
void ParticleSystem::reset()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			particles[i * cols + j].position[0] = 0 + j * restLengthValues[0];
			particles[i * cols + j].position[2] = 0 + -i * restLengthValues[1];
			particles[i * cols + j].position[1] = 0;
			//particles[i * cols + j].position.y = 0 + -i * padding;		//Useful if a vertical sheet is needed
			//particles[i * cols + j].position.z = 0;						//Useful if a vertical sheet is needed
			for (int k = 0; k < DIMENSION; k++)
			{
				particles[i * cols + j].velocity[k] = 0;
			}

		}
	}

}

//Update Method - Implements one time step for the animation
//Uses hooke's law with implicit methods to construct matrix A and vector b then calls the conjugate gradient method to find deltaV
//This is then used to update the particle velocities and in turn the positions
//Parameter - deltaT - Amount of time elapsed to use in integrating.  Type double. 
void ParticleSystem::doUpdate(double deltaT)
{
	//Initialize anything involving any type of summation to 0
	for (int i = 0; i < numParticles * dimensionSquared; i++)
	{
		ADiagonalEntries[i] = 0;
	}

	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		b[i] = 0;
		currentForce[i] = 0;
	}

	//Apply gravity
	for (int i = 1; i < numParticles * DIMENSION; i += DIMENSION)
	{
		currentForce[i] = -earthGravityValue * massMatrix[i];
	}
	
	double deltaTSquared = deltaT * deltaT;	//It is slightly more efficient to store this than recompute it

	for (int e = 0; e < numEdges; e++)
	{
		//Get the first and second particles from the edge, as well as the rest length
		int particle1 = edgeList[e].start;
		int particle2 = edgeList[e].end;
		double restLength = restLengthValues[edgeList[e].restLengthIndex];
		
		if (particle1 != particle2)
		{											
			//Do particle2 - particle1 -> the 1st direction of the spring

			double xij[DIMENSION];  //Position difference between particles
			double vij[DIMENSION];  //Velocity difference between particles
							
			double xijDotProduct = 0; //Dot product of xij vector
			for (int ii = 0; ii < DIMENSION; ii++)
			{
				xij[ii] = particles[particle2].position[ii] - particles[particle1].position[ii];
				vij[ii] = particles[particle2].velocity[ii] - particles[particle1].velocity[ii];
				xijDotProduct += xij[ii] * xij[ii];
			}
			double magnitude = sqrt(xijDotProduct); //Magnitude of xij vector
			
			/////////////////////////////
			//Do initial work for implicit components
			/////////////////////////////
			for (int ii = 0; ii < DIMENSION; ii++)	//Current row in Jacobian submatrix
			{
				for (int jj = 0; jj < DIMENSION; jj++) //Current column in Jacobian submatrix
				{
					//Current position and velocity jacobian submatrix elements
					double forcePositionJacobianSubMatrixElement = 0;
					double forceVelocityJacobianSubMatrixElement = 0;

					if (magnitude >= restLength)
					{
						//Choi's formula for the force-position Jacobian submatrix for spring forces - for magnitude >= restLength (' means transpose):
						//ks * (xij * xij') / (xij' * xij) + ks * (1 - restLength / magnitude) * (I - (xij * xij') / (xij' * xij));
						forcePositionJacobianSubMatrixElement = ks * xij[ii] * xij[jj] / xijDotProduct +
							ks * (1 - restLength / magnitude) * ((ii == jj) - (xij[ii] * xij[jj]) / xijDotProduct);
					}
					else
					{
						//Instead of following Choi here, we just cut off the second term
						//You cannot leave the second term here, or you will cause stability issues
						//ks * (xij * xij') / (xij' * xij)
						forcePositionJacobianSubMatrixElement = ks * xij[ii] * xij[jj] / xijDotProduct;
					}

					//Choi's formula for the force-velocity Jacobian submatrix for damping:
					//kd * I
					forceVelocityJacobianSubMatrixElement = (int)(ii == jj) * kd;

					//Complete formula for A (From Mucha's paper):
					//A = M - deltaT * forceVelocityJacobian - deltaT * deltaT * forcePositionJacobian
					//Diagonal and non-digonal components of A are computed separately because of the sparse matrix representation

					//Non-Diagonal-A elements
					//Implements complete formula here:
					//A = M - deltaT * forceVelocityJacobian - deltaT * deltaT * forcePositionJacobian
					//Note 1: There is no corresponding element in M for a non-diagonal element - thus we simply have negative signs
					//Note 2: This is only done once for both particles because our sparse matrix representation stores it once for both
					//Note 3: Non-Diagonal elements are indexed by edge, not particle id
					//For constrained particles, the edge list constrained attributes will zero out entries if needed
					ANonDiagonalEntries[e * dimensionSquared + ii*DIMENSION+jj] = -forcePositionJacobianSubMatrixElement * deltaTSquared
						- forceVelocityJacobianSubMatrixElement * deltaT;
					
					//Diagonal elements - we must decrement the stored value
					//These use the arrays for diagonal "negative summations" that correspond to each particle index
					//Diagonal A elements - particle 1
					//Implements this part of formula for A:
					//deltaT * forceVelocityJacobian + deltaT * deltaT * forcePositionJacobian (writing it with plus to emphasize subtraction takes place later)
					//The purpose of the -= (negative summation) is that diagonal elements of both forePositionJacobian and forceVelocityJacobian
					//are 1) the sum of all other forces acted on the particle 2) The negative result of those forces -- this implements Newton's 3rd Law for opposite forces
					ADiagonalEntries[particle1 * dimensionSquared + ii * DIMENSION + jj] -= (forcePositionJacobianSubMatrixElement * deltaTSquared
					 + forceVelocityJacobianSubMatrixElement * deltaT);
					
					//Diagonal A elements - particle 2
					//This implements the exact same logic for A as the statement above except that this is for particle2
					ADiagonalEntries[particle2 * dimensionSquared + ii * DIMENSION + jj] -= (forcePositionJacobianSubMatrixElement * deltaTSquared
					 + forceVelocityJacobianSubMatrixElement * deltaT);
					
					//Complete formula for b (from Mucha's paper):
					//b = deltaT * (currentExplicitForce + deltaT * forcePositionJacobian * currentVelocity)
					//Here we ONLY compute this part of b:
					//forcePositionJacobian * currentVelocity
					//Note 1: The diagonal forcePositionJacobian part is negated again to enforce Newton's 3rd Law
					//Note 2: the is unconstrained attribute will 0 out constrained particles in the Jacobian for non-diagonal entries

					//Particle 1
					b[particle1 * DIMENSION + ii] += (forcePositionJacobianSubMatrixElement * particles[particle2].velocity[jj] * edgeList[e].isStartUnconstrained  //Non diagonal element of position jacobian
						- forcePositionJacobianSubMatrixElement * particles[particle1].velocity[jj]); //Diagonal element of position Jacobian
					
					//Particle 2
					b[particle2 * DIMENSION + ii] += (forcePositionJacobianSubMatrixElement * particles[particle1].velocity[jj] * edgeList[e].isEndUnconstrained	//Non diagonal element of position jacobian
						- forcePositionJacobianSubMatrixElement * particles[particle2].velocity[jj]); //Diagonal element of position Jacobian
				}
			}

			//Calculate explicit spring force for this edge (used for finding b) based on Choi's formula
			//We are calculating the force for time n based on hooke's law - this will be used to help us approximate time n + 1 using backwards Euler integration
			//Note: xij, vij, and magnitude were already calculated for this edge above
			for (int i = 0; i < DIMENSION; i++)
			{
				//Do particle2 -> particle1 -> the 1st direction of the spring
				double forceExerted = (ks * (magnitude - restLength)) * (xij[i] / magnitude) + kd * vij[i];  //Change in force
				currentForce[particle1*DIMENSION+i] += forceExerted;
				currentForce[particle2*DIMENSION+i] -= forceExerted; //particle 1 -> particle 2 -- opposite direction
			}		
		} //if (particle1 != particle2)
	} //for (int e = 0; e < numEdges; e++)

	//Apply Constraints to Jacobians by setting rows in the matrix for the particles to 0
	//Note that here we only apply them to the diagonal entries
	//The isStartUnconstrained and isEndUnconstrained fields of the Edge struct handle all non-diagonal entries
	for (int p = 0; p < numParticles; p++)
	{
		bool zeroParticle = false;	//True if this particle is constrained and requires zeroing in the Jacobian
		for (int k = 0; k < numConstraints; k++)
		{
			if (p == constraintParticles[k])
			{
				zeroParticle = true;
				break;
			}
		}

		if (zeroParticle)
		{
			for (int ii = 0; ii < DIMENSION; ii++)
			{
				for (int jj = 0; jj < DIMENSION; jj++)
				{
					ADiagonalEntries[p * dimensionSquared + ii * DIMENSION + jj] = 0;
				}
			}
		}
	}

	//Zero out currentForce vector entries corresponding to constraints
	for (int k = 0; k < numConstraints; k++)
	{
		int constrained = constraintParticles[k];
		for (int ii = 0; ii < DIMENSION; ii++)
		{
			currentForce[constrained*DIMENSION+ii] = 0;
		}
		
	}

	//Finish calculating A for diagonal entries
	//Complete formula for A (from Mucha's paper):
	//A = M - deltaT * forceVelocityJacobian - deltaT * deltaT * forcePositionJacobian
	//Here we are adding the M and the subtraction to the formula
	//Note that we already calculated A for non-diagonal entries - they had no corresponding massMatrix element
	for (int ii = 0; ii < numParticles * dimensionSquared; ii++)
	{
		//Step 1 - Determine if this is a diagonal element - only diagonal elements have corresponding mass matrix entries, and the ADiagonalEntries array stores zeroes for nondiagonal elements
		//ii % (DIMENSION * DIMENSION) - For 3 dimensions this gets at the 9 elements for a given particle
		//% (DIMENSION + 1) - This checks whether or not entries are actually diagonal
		//Step 2 - Find the mass matrix index:
		//i / DIMENSION / DIMENSION (to get current particle #) * DIMENSION (because mass matrix is represented with a vector of size (DIMENSION * numParticles))
		//		--simplifies to i / DIMENSION
		ADiagonalEntries[ii] = (((ii % dimensionSquared % (DIMENSION + 1))==0)?(massMatrix[ii / DIMENSION]):0) - ADiagonalEntries[ii];
	}
	
	//Complete formula for b (from Mucha's paper):
	//b = deltaT * (currentExplicitForce + deltaT * forcePositionJacobian * currentVelocity)
	//At this point, we've already computed forcePositionJacobian * currentVelocity and stored this temp result in b
	//Now we finish the rest of the result
	for (int ii = 0; ii < numParticles * DIMENSION; ii++)
	{
		b[ii] = deltaT * (currentForce[ii] + deltaT * b[ii]);
	}

	//Log everything computed above
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		cout << "-------------------------------------------------------------------------------" << endl;

		logger -> printMassMatrix(massMatrix, numParticles, logger -> MEDIUM);
		
		logger -> printAFull(ADiagonalEntries,numParticles,ANonDiagonalEntries,numEdges,edgeList, "A", logger -> MEDIUM);
		
		logger -> printVector(currentForce,numParticles * DIMENSION, "currentForce (force found using explict methods used in finding b)", logger -> MEDIUM);
		
		logger -> printVector(b, numParticles * DIMENSION, "b", logger -> MEDIUM);

	}
	#endif

	//Find delta velocity through the conjugate gradient algorithm
	deltaV = findConjugateGradient(ADiagonalEntries, ANonDiagonalEntries, b, zeroVector);

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printVector(deltaV,numParticles*DIMENSION,"Delta V", logger -> MEDIUM);
	}
	#endif

	//If desired, only print constrained particles
	//Useful for verifying they are truly constrained
	#ifdef DEBUGGING
	if (logger -> printConstraintDeltaV && !logger -> isLogging) //No need to print if already printing all delta V's as part of larger debugging analsysis
	{
		cout << "DeltaV (only constrained particles shown):" << endl;
		for (int i = 0; i < numConstraints; i++)
		{
			int particleId = constraintParticles[i];
			cout << "Constrained Particle #" << particleId << ":" << endl;
			for (int j = particleId * DIMENSION; j < particleId * DIMENSION + DIMENSION; j++)
			{
				cout << setw(8) << deltaV[j] << " ";
			}
			cout << endl;
		}
	}
	#endif

		
	if (isAnimating)
	{
		for (int i = 0; i < numParticles; i++)
		{
			bool processParticle = true;				//True if current particle will move, false if not
			for (int k = 0; k < numConstraints; k++)
			{
				if (i == constraintParticles[k])
				{
					processParticle = false;
					break;
				}
			}

			if (processParticle)
			{
				for (int j = 0; j < DIMENSION; j++)
				{
					//Add delta velocity to each particle's velocity
					particles[i].velocity[j] += deltaV[i*DIMENSION + j];

					//Use explicit integration with velocity to update each particle's position
					particles[i].position[j] += particles[i].velocity[j] * deltaT;
				}
			}
		}
	}
	
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printParticles(particles, numParticles, logger -> MEDIUM);
	}
	#endif
}

//Method to find the solution x of Ax = b using the conjugate gradient algorithm
//Here x represents delta velocity
//Parameter ADiagonalEntries - Diagonal elements of A matrix
//Parameter ANonDiagonalEntries - Non-Diagonal entries of A matrix
//Parameter b - b vector
//x0 - Start point for algorithm (we just start at the zero vector)
//Based on the algorithm at http://en.wikipedia.org/wiki/Conjugate_gradient
double * ParticleSystem::findConjugateGradient(double * ADiagonalEntries, double * ANonDiagonalEntries, double * b, double * x0)
{
	const double ANSWER_BOUNDS = 1e-3;  //Acceptable bounds of answer to system of equations

	//You cannot use this algorithm to solve a system of equations when b is the zero vector
	//If this occurs, return the zero vector as the solution
	bool isBZero = true;
	
	for (int i = 0; i < DIMENSION * numParticles; i++)
	{
		if (abs(b[i]) > epsilon)
		{
			isBZero = false;
		}
	}

	if (isBZero)
	{
		return zeroVector;
	}
	
	//x = x0
	double * x = x0; 

	//r - start sums at 0
	for (int i = 0; i < DIMENSION * numParticles; i++)
	{
		r[i] = 0;
	}
		
	//r: A * x0 - part I
	for (int edgeNumber = 0; edgeNumber < numEdges; edgeNumber++)
	{
		int particle1 = edgeList[edgeNumber].start;
		int particle2 = edgeList[edgeNumber].end;
		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				r[particle1 * DIMENSION + i] += ANonDiagonalEntries[edgeNumber * dimensionSquared + i * DIMENSION + j] * x0[particle2 * DIMENSION + j] * edgeList[edgeNumber].isStartUnconstrained;
				r[particle2 * DIMENSION + i] += ANonDiagonalEntries[edgeNumber * dimensionSquared + i * DIMENSION + j] * x0[particle1 * DIMENSION + j] * edgeList[edgeNumber].isEndUnconstrained;
			}
		}

	}

	//r: A * x0 - part II
	for (int particleNumber = 0; particleNumber < numParticles; particleNumber++)
	{
		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				r[particleNumber * DIMENSION + i] += ADiagonalEntries[particleNumber * dimensionSquared + i * DIMENSION + j] * x0[particleNumber * DIMENSION + j];
			}
		}
	}
				
	//r: plug the result of A * x0 into the rest of the equation:
	//r = b - A * x0;
	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		r[i] = b[i] - r[i];
	}
	
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printVector(r, numParticles * DIMENSION, "r", logger -> FULL);
	}
	#endif
	
	
	//p = r
	for (int i = 0; i < numParticles * DIMENSION; i++)
	{
		p[i] = r[i];
	}

	while (true)
	{
		//Begin work on scalar alpha:
		//alpha = (r' * r) / (p' * A * p);

		//alpha step 1: (r' * r)
		//r.transpose is a (1, 3 * numParticles) size row vector
		//Thus p.transpose is also a (1,3 * numParticles) size row vector
		double dotProduct1 = 0;
		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			dotProduct1 += r[i] * r[i];
		}

		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			temp[i] = 0;
		}

		//alpha step 2: Multiply p' (row vector) * A (matrix) - Part I
		//It multiplies one column element in p' by one row element in A each time.
		//Thus we iterate through columns first and rows second
		for (int edgeNumber = 0; edgeNumber < numEdges; edgeNumber++)
		{
			int particle1 = edgeList[edgeNumber].start;
			int particle2 = edgeList[edgeNumber].end;
			for (int j = 0; j < DIMENSION; j++) //column
			{
				for (int i = 0; i < DIMENSION; i++) //row
				{
					temp[particle1 * DIMENSION + j] += p[particle2 * DIMENSION + i] * ANonDiagonalEntries[edgeNumber * dimensionSquared + i * DIMENSION + j] * edgeList[edgeNumber].isStartUnconstrained;
					temp[particle2 * DIMENSION + j] += p[particle1 * DIMENSION + i] * ANonDiagonalEntries[edgeNumber * dimensionSquared + i * DIMENSION + j] * edgeList[edgeNumber].isEndUnconstrained;
				}
			}

		}

		//alpha step 2: Muliply p' (row vector) * A (matrix) - Part II
		for (int particleNumber = 0; particleNumber < numParticles; particleNumber++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				for (int i = 0; i < DIMENSION; i++)
				{
					temp[particleNumber * DIMENSION + j] += p[particleNumber * DIMENSION + i] * ADiagonalEntries[particleNumber * dimensionSquared + i * DIMENSION + j];
				}
			}

		}
		
		//Alpha step 3: Multiply result of p' * A (it is temp and is a row vector) * p (a column vector)
		//The result is a scalar - dot product #2
		double dotProduct2 = 0;
		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			dotProduct2 += temp[i] * p[i];
		}

		//alpha step 4: Do final division to find alpha
		double alpha = dotProduct1 / dotProduct2;
		
		//nextX = x + alpha * p;
		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			nextX[i] = x[i] + alpha * p[i];
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger -> loggingLevel >= logger -> FULL)
			{
				cout << "Alpha:" << endl;
				cout << alpha << endl;
			}
		}
		#endif
				

		//NextR: start sums at 0
		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			nextR[i] = 0;
		}
		
		//Begin work on nextR
		//nextR = r - alpha * A * p;

		//NextR: (A * p) - part I
		for (int edgeNumber = 0; edgeNumber < numEdges; edgeNumber++)
		{
			int particle1 = edgeList[edgeNumber].start;
			int particle2 = edgeList[edgeNumber].end;
			for (int i = 0; i < DIMENSION; i++)
			{
				for (int j = 0; j < DIMENSION; j++)
				{
					nextR[particle1 * DIMENSION + i] += ANonDiagonalEntries[edgeNumber * dimensionSquared + i * DIMENSION + j] * p[particle2 * DIMENSION + j] * edgeList[edgeNumber].isStartUnconstrained;
					nextR[particle2 * DIMENSION + i] += ANonDiagonalEntries[edgeNumber * dimensionSquared + i * DIMENSION + j] * p[particle1 * DIMENSION + j] * edgeList[edgeNumber].isEndUnconstrained;
				}
			}

		}
		
		//NextR - (A * p) - part II
		for (int particleNumber = 0; particleNumber < numParticles; particleNumber++)
		{
			for (int i = 0; i < DIMENSION; i++)
			{
				for (int j = 0; j < DIMENSION; j++)
				{
					nextR[particleNumber * DIMENSION + i] += ADiagonalEntries[particleNumber * dimensionSquared + i * DIMENSION + j] * p[particleNumber * DIMENSION + j];
				}
			}

		}

		//NextR: plug the result of A * p into the rest of the equation
		//nextR = r - alpha * A * p;
		for (int i = 0; i < numParticles * DIMENSION; i++)
		{
			nextR[i] = r[i] - alpha * nextR[i];
		}
		
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> printVector(nextR, numParticles*DIMENSION,"nextR", logger -> FULL);
		}
		#endif

		//Begin work on beta:
		//beta = (nextR' * nextR) / (r' * r);

		//Beta step 1: (nextR' * nextR) - dot product of nextR
		double beta1 = 0;
		for (int i = 0; i < numParticles * DIMENSION; i++)
		{
			beta1 += nextR[i] * nextR[i];
		}

		//If within acceptable bounds of answer, exit this function
		if (sqrt(beta1) < ANSWER_BOUNDS) //Based directly on wikipedia matlab code
		{
			break;
		}

		//Beta step 2: (r' * r)
		double beta2 = 0;
		for (int i = 0; i < DIMENSION * numParticles; i++) //Can be done more compactly in a single for loop
		{
			beta2 += (r[i] * r[i]);
		}
		//Beta step 3: final division
		double beta = beta1 / beta2;

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger -> loggingLevel >= logger -> FULL)
			{
				cout << "Beta:" << endl;
				cout << beta << endl;
			}
		}
		#endif
		
		//nextP = nextR + beta * p;
		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			nextP[i] = nextR[i] + beta * p[i];
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> printVector(nextP, numParticles * DIMENSION, "nextP", logger -> FULL);
		}
		#endif
		
		//Move to next iteration
		for (int i = 0; i < DIMENSION * numParticles; i++)
		{
			x[i] = nextX[i];
			r[i] = nextR[i];
			p[i] = nextP[i];
		}
	}

	return nextX;

}

//Method to calculate the normals used for lighting
//Note that since the cross product is only defined in 3 dimensions, this method only works properly for 3 dimensions
void ParticleSystem::calculateNormals()
{
	//Cross product and this function only work if DIMENSION == 3
	if (DIMENSION != 3)
		return;

	double vectorDifferenceA[DIMENSION];  //1st vector formed for each cross product
	double vectorDifferenceB[DIMENSION];  //2nd vector formed for each cross product

	//Initialize items to 0
	for (int i = 0; i < numParticles; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			particles[i].normal[j] = 0;
		}
		particles[i].triangleCount = 0;
	}

	int edgeCounter = 0;	//Represents current edge number
	int i = 0, j = 0;		//For loop indexes
	for (i = 0; i < rows - 1; i++) //row
	{
		for (j = 0; j < cols - 1; j++) //column
		{
			Particle * topLeft = &particles[edgeList[edgeCounter].start];			//because this was the start of the first connection
			Particle * bottomLeft = &particles[edgeList[edgeCounter].end];			//because the first connection connected to the bottom
			Particle * bottomRight = &particles[edgeList[edgeCounter + 1].end];		//because the second connection connected to the bottom right
			Particle * topRight = &particles[edgeList[edgeCounter + 2].end];		//because the third connection connected to the right
	
			edgeCounter += 4;

			//Corresponds to this triangle:
			//glVertex3f(bottomLeft->position[0], bottomLeft->position[1], bottomLeft->position[2]); 
			//glVertex3f(topLeft->position[0], topLeft->position[1], topLeft->position[2]);
			//glVertex3f(topRight->position[0], topRight->position[1], topRight->position[2]);
			for (int i = 0; i < DIMENSION; i++)
			{
				vectorDifferenceA[i] = bottomLeft -> position[i] - topLeft -> position[i];
				vectorDifferenceB[i] = topRight -> position[i] - topLeft -> position[i];
			}

			double x = -(vectorDifferenceA[1] * vectorDifferenceB[2] - vectorDifferenceA[2] * vectorDifferenceB[1]);
			double y = -(vectorDifferenceA[2] * vectorDifferenceB[0] - vectorDifferenceA[0] * vectorDifferenceB[2]);
			double z = -(vectorDifferenceA[0] * vectorDifferenceB[1] - vectorDifferenceA[1] * vectorDifferenceB[0]);

			bottomLeft -> normal[0] += x;
			bottomLeft -> normal[1] += y;
			bottomLeft -> normal[2] += z;
			bottomLeft -> triangleCount++;
			topLeft -> normal[0] += x;
			topLeft -> normal[1] += y;
			topLeft -> normal[2] += z;
			topLeft -> triangleCount++;
			topRight -> normal[0] += x;
			topRight -> normal[1] += y;
			topRight -> normal[2] += z;
			topRight -> triangleCount++;
					
			//Corresponds to this triangle:
			//glVertex3f(bottomLeft->position[0], bottomLeft->position[1], bottomLeft->position[2]); 
			//glVertex3f(bottomRight->position[0], bottomRight->position[1], bottomRight->position[2]);
			//glVertex3f(topRight->position[0], topRight->position[1], topRight->position[2]);
			for (int i = 0; i < DIMENSION; i++)
			{
				vectorDifferenceA[i] = topRight -> position[i] - bottomRight -> position[i];
				vectorDifferenceB[i] = bottomLeft -> position[i] - bottomRight -> position[i];
			}

			x = -(vectorDifferenceA[1] * vectorDifferenceB[2] - vectorDifferenceA[2] * vectorDifferenceB[1]);
			y = -(vectorDifferenceA[2] * vectorDifferenceB[0] - vectorDifferenceA[0] * vectorDifferenceB[2]);
			z = -(vectorDifferenceA[0] * vectorDifferenceB[1] - vectorDifferenceA[1] * vectorDifferenceB[0]);

			bottomLeft -> normal[0] += x;
			bottomLeft -> normal[1] += y;
			bottomLeft -> normal[2] += z;
			bottomLeft -> triangleCount++;
			bottomRight -> normal[0] += x;
			bottomRight -> normal[1] += y;
			bottomRight -> normal[2] += z;
			bottomRight -> triangleCount++;
			topRight -> normal[0] += x;
			topRight -> normal[1] += y;
			topRight -> normal[2] += z;
			topRight -> triangleCount++;
			
		}

		edgeCounter++;
	}

	//For each vertex, we need a normal.  We are looking at the normal of each of the adjacent triangles and calculating an average normal.
	//We already did the summing for this average above.  Now we just need to divide by the number of triangles.
	//Once we find an average normal, we will then normalize it (make it unit lenth) for OpenGL's sake
	for (int i = 0; i < numParticles; i++)
	{
		double magnitude = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			particles[i].normal[j] /= particles[i].triangleCount;  //Find average direction of all normals for adjacent triangles
			magnitude += particles[i].normal[j] * particles[i].normal[j];
		}

		magnitude = sqrt(magnitude);

		for (int j = 0; j < DIMENSION; j++)
		{
			particles[i].normal[j] /= magnitude;  //Normalize the normal to be unit length
		}
		
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printNormals(particles,numParticles, logger -> MEDIUM);
	}
	#endif
}

//Render the particles
void ParticleSystem::doRender()
{
	glTranslatef(-halfWidth, 0.0f, halfHeight);

	if (renderMode == 1)
	{
		//Render surfaces
		int edgeCounter = 0;	//Current edge number
		int i = 0, j = 0;		//For loop indices
		for (i = 0; i < rows - 1; i++) //row
		{
			for (j = 0; j < cols - 1; j++) //column
			{
				Particle * topLeft = &particles[edgeList[edgeCounter].start];			//because this was the start of the first connection
				Particle * bottomLeft = &particles[edgeList[edgeCounter].end];			//because the first connection connected to the bottom
				Particle * bottomRight = &particles[edgeList[edgeCounter + 1].end];		//because the second connection connected to the bottom right
				Particle * topRight = &particles[edgeList[edgeCounter + 2].end];		//because the third connection connected to the right
	
				edgeCounter += 4;

				glBegin(GL_TRIANGLES);
				
				glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
				glNormal3f(bottomLeft->normal[0], bottomLeft->normal[1], bottomLeft->normal[2]);
				glVertex3f(bottomLeft->position[0], bottomLeft->position[1], bottomLeft->position[2]); 
				glNormal3f(topLeft->normal[0], topLeft->normal[1], topLeft->normal[2]);
				glVertex3f(topLeft->position[0], topLeft->position[1], topLeft->position[2]);
				glNormal3f(topRight->normal[0], topRight->normal[1], topRight->normal[2]);
				glVertex3f(topRight->position[0], topRight->position[1], topRight->position[2]);
		
				glNormal3f(bottomRight->normal[0], bottomRight->normal[1], bottomRight->normal[2]);
				glVertex3f(bottomRight->position[0], bottomRight->position[1], bottomRight->position[2]);
				glNormal3f(bottomLeft->normal[0], bottomLeft->normal[1], bottomLeft->normal[2]); 
				glVertex3f(bottomLeft->position[0], bottomLeft->position[1], bottomLeft->position[2]); 
				glNormal3f(topRight->normal[0], topRight->normal[1], topRight->normal[2]);
				glVertex3f(topRight->position[0], topRight->position[1], topRight->position[2]);
				
				glEnd();
			}

			edgeCounter++;
		}
	}
	else
	{
		//Lines Rendering logic (wireframe)
		glColor3f(0.0f, 1.0f, 0.0f);
		for (int i = 0; i < numEdges; i++)
		{
			Particle * one = &particles[edgeList[i].start];
			Particle * two = &particles[edgeList[i].end];
				
			glBegin(GL_LINES);
		
			glVertex3f(one->position[0], one->position[1], one->position[2]); 
			glVertex3f(two->position[0], two->position[1], two->position[2]);

			glEnd();
		}
	}

	if (renderToImage)
	{
		//Read the rendered image into a buffer
		glFlush();
		glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, screenBuffer);

		//Vertically flip the image in memory - it flips two rows from the screen in memory at a time
		for (int i = 0; i < windowHeight / 2; i++)
		{
			if (i != windowHeight - i)
			{
				memcpy(screenRowTemp,(uint8_t *)&screenBuffer[i * windowWidth * 4], windowWidth * 4 * sizeof(uint8_t));
				memcpy((uint8_t *)&screenBuffer[i * windowWidth * 4], (uint8_t *)&screenBuffer[(windowHeight - i) * windowWidth * 4], windowWidth * 4 * sizeof(uint8_t));
				memcpy((uint8_t *)&screenBuffer[(windowHeight - i) * windowWidth * 4], screenRowTemp, windowWidth * 4 * sizeof(uint8_t));
			}
		}
		
		//Write a new targa file, with a new number
		sprintf(imageFileName, "images\\ImplicitMethods%d.tga", frameNumber);
		tga_result result = tga_write_rgb(imageFileName, screenBuffer, windowWidth, windowHeight, 32);
		#ifdef DEBUGGING
		if (logger -> isLogging && logger -> loggingLevel >= logger -> LIGHT)
		{
			if (result != TGA_NOERR)
			{
				cout << tga_error(result) << " at frame number " << frameNumber << endl;
			}
		}
		#endif

		frameNumber++;
	}
		
	//Render onscreen text with informational messages
	if (showInfoText)
	{
		//Reference: http://programming-technique.blogspot.com/2012/05/font-rendering-in-glut-using-bitmap.html
		glDisable(GL_LIGHTING);
		glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, windowWidth, 0, windowHeight);
		
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glRasterPos2i(20, 20);
		for (int i = 0; i < strlen(text); i++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
		}
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glEnable(GL_LIGHTING);
	}

}


//This method stores the current eye position (from the camera) for this mesh
//parameter eyePos - the eye position to store
void ParticleSystem::setEyePos(glm::vec3 & eyePos)
{
	this -> eyePos[0] = eyePos[0];
	this -> eyePos[1] = eyePos[1];
	this -> eyePos[2] = eyePos[2];
}

void ParticleSystem::initVBOs()
{
	//Spring mesh (normal rendering)
	glGenBuffers(1, vboHandle);
	glGenBuffers(1, indexVboHandle);

	//Spring mesh (special edge rendering)
	glGenBuffers(1, gridVboHandle);
	glGenBuffers(1, gridIndexVboHandle);
}

void ParticleSystem::sendVBOs()
{
	if (renderMode)
	{
		//Spring Mesh - full rendering
		glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Particle) * numParticles, particles, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		indices.clear();

		//Render surfaces
		int edgeCounter = 0;	//Current edge number
		int i = 0, j = 0;		//For loop indices
		for (i = 0; i < rows - 1; i++) //row
		{
			for (j = 0; j < cols - 1; j++) //column
			{
				int topLeft = edgeList[edgeCounter].start;			//because this was the start of the first connection
				int bottomLeft = edgeList[edgeCounter].end;			//because the first connection connected to the bottom
				int bottomRight = edgeList[edgeCounter + 1].end;		//because the second connection connected to the bottom right
				int topRight = edgeList[edgeCounter + 2].end;		//because the third connection connected to the right
	
				edgeCounter += 4;

				indices.push_back(bottomLeft);
				indices.push_back(topLeft);
				indices.push_back(topRight);
			
				indices.push_back(bottomRight);
				indices.push_back(bottomLeft);
				indices.push_back(topRight);
			}

			edgeCounter++;
		}

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVboHandle[0]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * indices.size(), &indices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	else
	{
		//Clear vertex and index lists for lines rendering logic (wireframe)
		gridVertices.clear();
		gridIndices.clear();

		//Lines Rendering logic (wireframe)
		for (int i = 0; i < numEdges; i++)
		{
			Particle * one = &particles[edgeList[i].start];
			Particle * two = &particles[edgeList[i].end];

			gridVertices.push_back(*one);
			gridVertices.push_back(*two);

			gridIndices.push_back(i * 2);
			gridIndices.push_back(i * 2 + 1);
		}

		glBindBuffer(GL_ARRAY_BUFFER, gridVboHandle[0]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Particle) * gridVertices.size(), &gridVertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gridIndexVboHandle[0]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * gridIndices.size(), &gridIndices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
 
}

void ParticleSystem::doRender(glm::mat4 & projMatrix, glm::mat4 & modelViewMatrix)
{
	sendVBOs();

	modelViewMatrix = glm::translate(modelViewMatrix, glm::vec3(-halfWidth, 0.0f, halfHeight)); 

	glm::mat4 totalMatrix = projMatrix * modelViewMatrix;
	glm::mat4 normalMatrix = glm::inverse(modelViewMatrix);
	normalMatrix = glm::transpose(normalMatrix);

	//Spring mesh rendering
	glUseProgram(programObject);

	//Parameter setup
	GLuint c0 = glGetAttribLocation(programObject, "position");
	GLuint c1 = glGetAttribLocation(programObject, "normal");
	GLuint c2 = glGetAttribLocation(programObject, "color1");

	GLuint l1 = glGetUniformLocation(programObject, "lightAmbient");
	GLuint l2 = glGetUniformLocation(programObject, "lightDiffuse");
	GLuint l3 = glGetUniformLocation(programObject, "lightSpecular");
	GLuint lp = glGetUniformLocation(programObject, "lightPosition");
	GLuint ep = glGetUniformLocation(programObject, "eyePosition");

	GLuint d1 = glGetUniformLocation(programObject, "ambient_coef");
	GLuint d1b = glGetUniformLocation(programObject, "ambient_back_coef");
	GLuint d2 = glGetUniformLocation(programObject, "diffuse_coef");
	GLuint d2b = glGetUniformLocation(programObject, "diffuse_back_coef");
	GLuint d3 = glGetUniformLocation(programObject, "specular_coef");
	GLuint d4 = glGetUniformLocation(programObject, "mat_shininess");

	GLuint m1 = glGetUniformLocation(programObject, "local2clip");
	GLuint m2 = glGetUniformLocation(programObject, "local2eye");
	GLuint m3 = glGetUniformLocation(programObject, "normalMatrix");

	glEnableVertexAttribArray(c0); 
	glEnableVertexAttribArray(c1);
    glEnableVertexAttribArray(c2); 

	if(renderMode)
	{
		glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVboHandle[0]);
	}
	else
	{
		glBindBuffer(GL_ARRAY_BUFFER, gridVboHandle[0]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gridIndexVboHandle[0]);
	}

	glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(Particle),(char*) NULL+0); 
	glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(Particle),(char*) NULL+16); 
    glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(Particle),(char*) NULL+32); 

	//If in ambient mode, add extra ambience to make things really bright
	//Otherwise use normal ambience
	if (ambientMode)
	{
		glUniform4f(l1, lightFullAmbient[0], lightFullAmbient[1], lightFullAmbient[2], 1.0);
	}
	else
	{
		glUniform4f(l1, lightAmbient[0], lightAmbient[1], lightAmbient[2], 1.0);
	}
	glUniform4f(l2, lightDiffuse[0], lightDiffuse[1], lightDiffuse[2], 1.0);
	glUniform4f(l3, lightSpecular[0], lightSpecular[1], lightSpecular[2],1.0);
	glUniform4f(lp, lightPosition[0], lightPosition[1], lightPosition[2], lightPosition[3]);
	glUniform4f(ep, eyePos[0], eyePos[1], eyePos[2], 1); 

	glUniform4f(d1, matAmbient[0], matAmbient[1], matAmbient[2], 1.0);
	glUniform4f(d1b, matAmbientBack[0], matAmbientBack[1], matAmbientBack[2], 1.0);
	glUniform4f(d2, matDiffuse[0], matDiffuse[1], matDiffuse[2], 1.0);
	glUniform4f(d2b, matDiffuseBack[0], matDiffuseBack[1], matDiffuseBack[2], 1.0);
	glUniform4f(d3, matSpecular[0], matSpecular[1], matSpecular[2],1.0);
	glUniform1f(d4, matShininess[0]);

	glUniformMatrix4fv(m1, 1, GL_FALSE, &totalMatrix[0][0]);
	glUniformMatrix4fv(m2, 1, GL_FALSE, &modelViewMatrix[0][0]);
	glUniformMatrix4fv(m3, 1, GL_FALSE, &normalMatrix[0][0]);

	//Draw
	if (renderMode)
	{
		//Normal cloth rendering
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, (char *) NULL + 0);
	}
	else
	{
		//Special grid wireframe mode
		glDrawElements(GL_LINES, gridIndices.size(), GL_UNSIGNED_INT, (char *) NULL + 0);
	}

	glUseProgram(0);

	if (renderToImage)
	{
		//Read the rendered image into a buffer
		glFlush();
		glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, screenBuffer);

		//Vertically flip the image in memory - it flips two rows from the screen in memory at a time
		for (int i = 0; i < windowHeight / 2; i++)
		{
			if (i != windowHeight - i)
			{
				memcpy(screenRowTemp,(uint8_t *)&screenBuffer[i * windowWidth * 4], windowWidth * 4 * sizeof(uint8_t));
				memcpy((uint8_t *)&screenBuffer[i * windowWidth * 4], (uint8_t *)&screenBuffer[(windowHeight - i) * windowWidth * 4], windowWidth * 4 * sizeof(uint8_t));
				memcpy((uint8_t *)&screenBuffer[(windowHeight - i) * windowWidth * 4], screenRowTemp, windowWidth * 4 * sizeof(uint8_t));
			}
		}
		
		//Write a new targa file, with a new number
		sprintf(imageFileName, "images\\ImplicitMethods%d.tga", frameNumber);
		tga_result result = tga_write_rgb(imageFileName, screenBuffer, windowWidth, windowHeight, 32);
		#ifdef DEBUGGING
		if (logger -> isLogging && logger -> loggingLevel >= logger -> LIGHT)
		{
			if (result != TGA_NOERR)
			{
				cout << tga_error(result) << " at frame number " << frameNumber << endl;
			}
		}
		#endif

		frameNumber++;
	}

	//Render onscreen text with informational messages
	if (showInfoText)
	{
		//Reference: http://programming-technique.blogspot.com/2012/05/font-rendering-in-glut-using-bitmap.html
		glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, windowWidth, 0, windowHeight);
		
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glRasterPos2i(20, 20);
		for (int i = 0; i < strlen(text); i++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
		}
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}


 			
}

//This method iteratively keeps applying forces until it detects it has reached equilibrium
//The alternative methode is a large system of equations on paper, since we would need to compute
//restLength of edge1, restLength of edge 2, restLength of edge 3, etc.
//The method will keep running until the amount of change for all forces is less than threshold or until the number of iterations exceeds constant MAX_ITERATIONS
//Note: Excessively large spring constants and some damping constants will mathematically break an explicit method such as this, which can cause the grid to disappear
//Parameter threshold is the threshold described above
void ParticleSystem::goToRestLengthsForGravity(double threshold)
{
	double deltaT = 0.005;				//Amount of time to simulate processing at one time while doing our explict estimate
	const int MAX_EXECUTIONS = 50000;	//Max iterations for this method.  Some systems will never come to rest.  This prevents infinite loops in those cases.  

	int totalExecutions = 0;			//Number of iterations that have actually occurred
	bool keepGoing = true;				//Boolean flag checked each loop iteration - loop terminates when it's false
	while(keepGoing)
	{
		keepGoing = false;
		
		for (int i = 0; i < numParticles * DIMENSION; i++)
		{
			currentForce[i] = 0;
		}

		double earthGravity = -earthGravityValue;
		for (int i = 0; i < numParticles; i++)
		{
			currentForce[i * DIMENSION + 1] += earthGravity * massMatrix[i*DIMENSION];
		}
		
		for (int e = 0; e < numEdges; e++)
		{
			int particle1 = edgeList[e].start;
			int particle2 = edgeList[e].end;
			double restLength = restLengthValues[edgeList[e].restLengthIndex];

		
			if (particle1 != particle2)
			{
				double xij[DIMENSION];	//Change in position between the two particles
				double vij[DIMENSION];	//Change in velocity
				double magnitude = 0;	//Calculated magnitude of xij

				for (int i = 0; i < DIMENSION; i++)
				{
					xij[i] = particles[particle2].position[i] - particles[particle1].position[i];
					vij[i] = particles[particle2].velocity[i] - particles[particle1].velocity[i];
					magnitude += xij[i] * xij[i];
				}

				magnitude = sqrt(magnitude);
				
				//We're dealing with the force another particle exerts on our current particle
				//The force that the other exerts must be negated in order to enforce Newton's Law
				//Thus it's positive rather than negative.
				for (int i = 0; i < DIMENSION; i++)
				{
					double forceChange = (ks * (magnitude - restLength)) * (xij[i] / magnitude) + kd * vij[i];
					currentForce[particle1*DIMENSION+i] += forceChange;
					currentForce[particle2*DIMENSION+i] -= forceChange;
				}
			} //if (particle1 != particle2)

		} //for (int e = 0; e < numEdges; e++)

		for (int i = 0; i < numParticles; i++)
		{	
			bool processParticle = true;
			for (int j = 0; j < numConstraints; j++)
			{
			if (i == constraintParticles[j])
				{
					processParticle = false;
					continue;
				}
			}
			
			if (processParticle)
			{
				for (int j = 0; j < DIMENSION; j++)
				{
					particles[i].velocity[j] += currentForce[i*DIMENSION + j] / massMatrix[i * DIMENSION] * deltaT;
					particles[i].position[j] += particles[i].velocity[j] * deltaT;
				}
			}
		}
		
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger -> loggingLevel >= logger -> MEDIUM)
			{
				for (int i = 0; i < numParticles; i++)
				{
					cout << "Force for Particle " << i << ": " << currentForce[i * DIMENSION] << ", " << currentForce[i * DIMENSION + 1] << ", " << currentForce[i * DIMENSION + 2] << endl;
				}
			}
		}
		#endif
		
				
		for (int i = 0; i < numParticles * DIMENSION; i++)
		{
			bool checkParticle = true;
			for (int j = 0; j < numConstraints; j++)
			{
				if (constraintParticles[j] == i / DIMENSION) //Integer division
				{
					checkParticle = false;
					break;
				}

			}
			
			if (checkParticle)
			{
				if (fabs(currentForce[i]) > threshold)
				{
					keepGoing = true;
					break;
				}
			}
		}
	
		totalExecutions++;
		if (totalExecutions > MAX_EXECUTIONS) //If an issue exists where the object "breaks free" (ie having too small a spring constant), we need a way to keep this method from infinitely looping.
		{
			cout << "Was not able to reach equilibrium in the system within " << MAX_EXECUTIONS << " iterations" << endl;
			break;
		}
		
	} //While (keepGoing)
	
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		if (logger -> loggingLevel >= logger -> LIGHT)
		{
			cout << "Completed execution of the explicit equilibrium finder" << endl;
			cout << "Total executions:" << totalExecutions << endl;
		}
	}
	#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////User Interface methods to alter attributes of the particle system///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Method to increase the earth gravity used in the particle system
//Parameter amount - the amount to increase.  Use a negative value to decrease
void ParticleSystem::increaseEarthGravity(double amount)
{
	earthGravityValue += amount;
	if (true)
	{
		cout << "Earth Gravity value is now " << earthGravityValue << endl;
	}

	sprintf(text, "Earth Gravity: %f", earthGravityValue);
}

//Method to increase the straight rest length
//Parameter amount - the amount by which to increase it (a negative value decreases)
void ParticleSystem::increaseStraightRestLength(double amount)
{
	restLengthValues[0] += amount;
	restLengthValues[1] += amount;
	if (restLengthValues[0] < epsilon || restLengthValues[1] < epsilon)	//Rest length of 0 results in bad calculation of the b vector.  Negative value doesn't make sense.
    {
		restLengthValues[0] -= amount;
		restLengthValues[1] -= amount;
		return;
    }

	restLengthValues[2] = sqrt(restLengthValues[0] * restLengthValues[0] + restLengthValues[1] * restLengthValues[1]);  //Use pythagorean theorem to calculate diagonal restLength
	

	//Half of the width and half of the height have now changed and must be calculated again
	sheetWidth += amount * (cols - 1);
	sheetHeight += amount * (rows - 1);
	halfWidth = sheetWidth / 2;
	halfHeight = sheetHeight / 2;

	//Recompute all constrained particles in the grid.
	//Non constrained particles will be allowed to "fall in line" due to their spring forces
			
	int currentParticle = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			currentParticle = i * cols + j;
			for (int k = 0; k < numConstraints; k++)
			{
				if (currentParticle == constraintParticles[k])
				{
					particles[currentParticle].position[0] = 0 + j * restLengthValues[0];
					particles[currentParticle].position[2] = 0 + -i * restLengthValues[1];
					particles[currentParticle].position[1] = 0;
					//particles[i * cols + j].position.y = 0 + -i * padding;		//Useful for a vertical sheet
					//particles[i * cols + j].position.z = 0;						//Useful for a vertical sheet
					for (int m = 0; m < DIMENSION; m++)
					{
						particles[currentParticle].velocity[m] = 0;
					}
					break;
				}
			}

		}
	}

	if (true)
	{
		cout << "Straight Rest Length is now  " << restLengthValues[0] << endl;
		cout << "Diagonal Rest Length is now " << restLengthValues[1] << endl;
	}

	sprintf(text, "Straight Rest Length: %f", restLengthValues[0]);
}

//Method to toggle whether or not informational messages are shown on screen
void ParticleSystem::toggleInfoText()
{
	showInfoText = !showInfoText;
}

//Method to set the window width and height
//Parameters width and height - the window dimensions
void ParticleSystem::setWindowDimensions(int width, int height)
{
	windowWidth = width;
	windowHeight = height;
}

//Method to toggle whether or not animation occurs
//Particles will not have their positions or velocities updated if animation is turned off
void ParticleSystem::toggleAnimation()
{
	isAnimating = !isAnimating;
}

//Method to toggle between rendering surfaces and wireframe mode
void ParticleSystem::toggleRenderMode()
{
	renderMode = !renderMode;
	//if (renderMode)
	//{
	//	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//}
	//else
	//{
	//	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//}
}

void ParticleSystem::toggleAmbientMode()
{
	ambientMode = !ambientMode;
}

//Method to toggle between rendering to a series of numbered images and not rendering to them
//It ALWAYS renders to the screen regardless
void ParticleSystem::toggleImageRendering()
{
	renderToImage = !renderToImage;
	if (renderToImage)
	{
		sprintf(text, "Image rendering on");
	}
	else
	{
		sprintf(text, "Image rendering off");
	}
}