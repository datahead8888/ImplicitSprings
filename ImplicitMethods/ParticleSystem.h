#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <gl/glut.h>
#include <gl/GLU.H>
#include <vector>
#include "Particle.h"
#include "Edge.h"
#include "Logger.h"
#include "targa.h"

#define TEXT_SIZE 256	//Maximum length of the on screen message text

//Particle System class
//It builds a grid of particles (represented by Particle objects)
//These are then held together by springs using hooke's laws.  One edge corresponds to each spring.
//Hooke's law is used to represent the springs, with both a spring component and a damping comopnent
//Implicit methods are used for the integration.  This requires a solving system but allows much larger spring constants without instability
//and potentially allows more efficient implementation
//A gravity component is also present
class ParticleSystem
{
	public:
	ParticleSystem(Logger * logger);
	~ParticleSystem();
	void initVBOs();
	void sendVBOs();
	void reset();
	void doUpdate(double elapsedSeconds);
	double * findConjugateGradient(double * ADiagonalEntries, double * ANonDiagonalEntries, double * b, double * x0);
	void calculateNormals();
	void doRender();
	void doRender(glm::mat4 & projMatrix, glm::mat4 & modelViewMatrix);
	//UI Methods
	void increaseEarthGravity(double amount);
	void increaseStraightRestLength(double amount);
	void toggleInfoText();
	void setWindowDimensions(int width, int height);
	void toggleAnimation();
	void toggleRenderMode();
	void toggleImageRendering();
	void goToRestLengthsForGravity(double threshold);
	void setProgramObject(GLuint programObject) {this->programObject = programObject;}
	void setEyePos(glm::vec3 & eyePos);
		
	private:
	double halfWidth;					//Half the width of the original grid.  Used to make the grid initially be centered.
	double halfHeight;					//Half the hight of the original grid.  Used to make the grid initially be centered.
	Edge * edgeList;					//Represents the connections between particles (one spring per connection)
	int numEdges;						//Number of edges in the grid
	Particle * particles;				//Collection of Particle struct objects representing the actual particles (vertices)
	std::vector<int> indices;				//Mesh indices for GLSL rendering
	int numParticles;					//Number of particles in the system
	
	double * b;							//b vector (used to implement Ax = b solver)
	double * ADiagonalEntries;			//Sparse matrix representation of diagonal portion of A (for Ax = b)
	double * ANonDiagonalEntries;		//Sparse matrix representation of non-diagonal portion of A (for Ax = b)
	double * massMatrix;				//Mass Matrix (for 3 dimensions, each set of 3 entries should be the same - mass for one particle)
	double * zeroVector;				//Vector containing all 0's
	double * deltaV;					//Computed change in velocity (Result of Ax = b system -- x here represents deltaV)
	double * currentForce;				//Current force, calculated explicitly but used in implicit calculations

	//The following are temporary variables used in the conjugate gradient implementation
	double * nextX;
	double * r;
	double * nextP;
	double * temp;
	double * nextR;
	double * p;

	int dimensionSquared;				//DIMENSION * DIMENSION occurs so frequently that a lot of computation can be saved by storing this

	bool isAnimating;					//True if particles should move; false if not
	bool showInfoText;					//True if informational messages should be rendered on screen
	double kd;							//Damping constant (Following Choi's name of kd)
	double ks;							//Spring constant
	double earthGravityValue;			//Acceleration rate for gravity
	int * constraintParticles;			//Array containing indexes of all constrained particles)
	int numConstraints;					//Number of constrained particles
	double restLengthValues[3];			//Array holding all differrent types of restLength values, which can be referenced by enumerated values
	int rows;							//Number of rows in the grid
	int cols;							//Number of columns in the grid
	double sheetWidth;					//Width of sheet of cloth (in OpenGL coordinates)
	double sheetHeight;					//Height of sheet of cloth

	int renderMode;						//If one will render surfaces.  If not one, will render a wireframe representation.
	char text[TEXT_SIZE];				//Used for on screen informational messages
	int windowWidth, windowHeight;		//Window width and height used for ortho mode when displaying informational text
	
	//Video generation variables (generates a series of numbered images that can be combined into a video with a tool)
	bool renderToImage;					//If true will begin rendering to series of numbered images
	uint8_t * screenBuffer;				//Buffer to hold screen capture
	uint8_t * screenRowTemp;			//Buffer to temporarily hold one row from the screen while flipping image verticaly
	int frameNumber;					//Numbered frame for series of images to be video
	char imageFileName[100];			//Name of output image file

	Logger * logger;					//Reference to Logger class to perform all logging

	GLuint vboHandle[1];	  //handle to vertex buffer object for vertices
	GLuint indexVboHandle[1]; //handle to vertex buffer object for indices

	GLuint programObject;				//Program object needed for shaders (notably lighting)
	GLfloat eyePos[3];		  //Position of the eye (for the camera)

	GLfloat lightAmbient[4];  //Normal ambient light setting
	GLfloat lightDiffuse[4];  //Diffuse light setting
	GLfloat lightSpecular[4]; //Specular light setting
	GLfloat lightPosition[4]; //Light position setting

	GLfloat matAmbient[4];    //Material ambient setting
	GLfloat matAmbientBack[4];    //Material ambient setting
	GLfloat matDiffuse[4];    //Material diffuse setting
	GLfloat matDiffuseBack[4];    //Material diffuse setting
	GLfloat matSpecular[4];   //Material specular setting
	GLfloat matShininess[1];  //Material shininess setting

	GLfloat lightFullAmbient[4];  //Special "full" ambient setting for debugging

	bool ambientMode;
};
