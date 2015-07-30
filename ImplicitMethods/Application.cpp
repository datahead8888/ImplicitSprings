//Cloth simulation using implicit methods to implement springs for a grid of particles
//Features:
//	*Implements springs between particles defined in an edge list using hooke's law with a damping component
//	*Has an earth gravity component to drag the cloth downwards
//	*Uses implicit methods, allowing for a higher spring constant (variable ks in class ParticleSystem)
//	*Uses a conjugate gradient solver to implement implicit methods
//	*Gourad shading for lighting
//Controls:
//Keyboard:
//	SPACE: Begin animating
//	A: Decrease rest length between particles
//	D: Increase rest length between particles
//	W: Increase earth gravity (starts at 9.8 m/sec)
//	S: Decrease earth gravity
//	Z: toggle wireframe mode
//	X: toggle informational text display
//  E: run an explicit implementation of the simulation (useful for comparison; most obvious if you turn automatic implicit animation off with space bar)
//  R: reset the simulation
//  I: render to a series of numbered images so that they can be combined into a video
//	O: toggle logging of positions/velocities of constrained particles (only if DEBUGGING macro is #defined in Logger.h)
//	P: toggle complete logging (only if DEBUGGING macro is #defined in Logger.h)
//Mouse:
//	Left button - hold this whie dragging the mouse to change the rotation angle of the piece of cloth shown
//	Middle button - zoom in
//	Right button - zoom out
//Written by Chris Jacobsen with advisement from Professor Huamin Wang

//Based on:
//http://graphics.snu.ac.kr/~kjchoi/publication/cloth.pdf - Explicit and implicit formulas for hooke's law - page 3
//http://www.amath.unc.edu/Faculty/mucha/Reprints/SCAclothcontrolpreprint.pdf - Formulas for A and b in Ax = b system - page 5
//http://en.wikipedia.org/wiki/Conjugate_gradient - Algorithm for conjugate gradient (solves Ax = b)

#include <gl/glew.h>
#include <gl/glut.h>
#include <gl/GLU.H>

#include <string>
#include <sstream>
#include <iostream>
#include "ParticleSystem.h"
#include "ViewManager.h"
#include "Keyboard.h"
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

using namespace std;

GLuint SetupGLSL(char *fileName);

//Note: the reason these were declared globally is to accommodate Glut's function calling system
ParticleSystem * particleSystem;	//The main particle system
ViewManager viewManager;			//Instance of the view manager to allow user view control
Keyboard * keyboard;				//Instance of the Keyboard class to process key presses
Logger * logger;					//Instance of Logger class to perform all logging
double ar = 0;

int summedTime = 0;
int frameCount = 0;
const int FRAME_COUNT_LIMIT = 500;

//This function is called for rendering by GLUT
void render()
{
	//Timing mechanism for performance evaluation
	int startTime, endTime;
	
	startTime = glutGet(GLUT_ELAPSED_TIME);
	
	//Update Logic
	double timeElapsed = 0.01;  //Size of time step
	particleSystem -> doUpdate(timeElapsed);
	particleSystem -> calculateNormals();

	//Render Logic
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glLoadIdentity();

	//viewManager.doTransform();
	//particleSystem -> doRender();

	glm::mat4 projMatrix = glm::perspective(45.0f, (float)ar, 0.001f, 100.0f); //Projection Matrix
	double cameraHeight = 0.5;  //Max height of camera if looking straight down at character
	
	glm::mat4 modelViewMatrix = viewManager.doTransform();

	particleSystem -> doRender(projMatrix, modelViewMatrix);

	glFlush();

	glutSwapBuffers();
		
	//Calculate time this frame took to compute and render
	endTime = glutGet(GLUT_ELAPSED_TIME);
	//cout << "Total time for frame was: " << (double)(endTime - startTime) / 1000 << " seconds" << endl;
	frameCount++;
	summedTime += (endTime - startTime);

	if (frameCount == FRAME_COUNT_LIMIT) {
		double fps = frameCount / (summedTime / double(1000));
		cout << "Average FPS for " << FRAME_COUNT_LIMIT << " frames is: " << fps << endl;
		frameCount = 0;
		summedTime = 0;
	}
	

	

}

//This function processes window resize events, ensuring the aspect ratio and anything else dependent on window dimensions is correct
//Parameters width and height are the dimensions of the window
void resize(int width, int height)
{
	//Minimizing the window results in sizes of 0, which causes problems with video generation.
	//Do not allow the width or height to be set to 0.
	if (width == 0 || height == 0)
	{
		return;
	}

	

	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();

	//Set the Open GL viewport to the window's dimensions
	glViewport(0, 0, width, height);
	//Set up perspective.
	//45 degree viewing angle, use the above calculated aspect ratio,
	//set up near clipping plane to 0.1, the far clipping pane to 50.0
	//(Clipping plane values must be positive)
	//gluPerspective(45.0, ar, 0.1, 50.0);

	//Calculate the aspect ratio using the x and y dimensions
	ar = (double) width / height;

	//Go back to MODEL VIEW matrix mode
	//glMatrixMode(GL_MODELVIEW);
	particleSystem -> setWindowDimensions(width,height);
}

//This function processes mouse clicks (when the button is pressed)
//Parameters:
//button - which button was pressed
//state - was it up or down?
//x - x coordinate of click
//y - y coordinte of click
void mouseClick(int button, int state, int x, int y)
{
	viewManager.mouseClick(button,state,x,y);
}

//This function processes mouse motion.
//It only takes action if a button is being held
//Parameters:
//x - x coordinate of movement
//y - y coordinate of movement
void mouseMove(int x, int y)
{
	viewManager.mouseMove(x,y);
}

//This function processes keyPress events (a keyboard button going down)
//Parameter key - the button being pressed
//The other parameters are not used
void keyPressed (unsigned char key, int mystery, int mystery2)
{
	keyboard -> keyPressed(key);
}

//This function processes keyRelease events (a keyboard button going up)
//Parameter key - the button being released
//The other parameters are not used
void keyReleased (unsigned char key, int mystery, int mystery2)
{
	keyboard -> keyReleased(key);
}

//Main function
int main(int argCount, char **argValue)
{
	logger = new Logger();
	particleSystem = new ParticleSystem(logger);
	keyboard = new Keyboard(particleSystem, logger);

	glutInit(&argCount,argValue);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE | GLUT_RGBA);

	glutInitWindowPosition(0,0);
	int windowWidth = 1000;
	int windowHeight = 700;
	glutInitWindowSize(windowWidth, windowHeight);
	resize(windowWidth, windowHeight);
	particleSystem -> setWindowDimensions(windowWidth, windowHeight);
	glutCreateWindow("Implicit Methods");

	//Note: this must be AFTER the create window function call to work properly
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);

	//Set up lighting attributes
	//GLfloat lightPosition[] = {0.0, 0,0, 1.0, 0.0};
	//GLfloat lightAmbient[] = {1.0f, 1.0f, 1.0f, 1.0f};
	//GLfloat lightDiffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
	//GLfloat lightSpecular[] = {1.0f, 1.0f, 1.0f, 1.0f};
	//glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
	//glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	//glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);

	//GLfloat dummy[2];  //The glLightModelfv requires a pointer.  Since it's not zero here, this will result in two sided lighting
	//glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, dummy);

	//GLfloat materialAmbientGreen[] = {0.0, 0.0, 0.0, 1.0}; //No ambient color
	//GLfloat materialDiffuseGreen[] = {0.0, 0.7, 0.0, 1.0};
	//GLfloat materialAmbientBlue[] = {0.0, 0.0, 0.0, 1.0};  //No ambient color
	//GLfloat materialDiffuseBlue[] = {0.0, 0.0, 0.7, 1.0};
	//glMaterialfv(GL_BACK, GL_AMBIENT, materialAmbientBlue);
	//glMaterialfv(GL_FRONT, GL_AMBIENT, materialAmbientGreen);
	//glMaterialfv(GL_BACK, GL_DIFFUSE, materialDiffuseBlue);
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, materialDiffuseGreen);

	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	glClearColor(0.0f,0.0f,0.0f,0.0f);
	
	GLenum errorCode = glewInit();
	GLuint programObject = SetupGLSL("maze");
	particleSystem->initVBOs();
	particleSystem->setProgramObject(programObject);

	glutDisplayFunc(render);
	glutIdleFunc(render);
	glutReshapeFunc(resize);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);
	glutKeyboardFunc(keyPressed);
	glutKeyboardUpFunc(keyReleased);
	glutMainLoop();

	delete keyboard;
	delete particleSystem;
	delete logger;

	return 0;
}