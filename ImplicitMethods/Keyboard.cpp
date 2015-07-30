#include "Keyboard.h"
#include <gl/glut.h>

//Constructor for Keyboard class
//Parameter theSystem - a pointer reference to the ParticleSystem class
//Parameter logger - a pointer reference to the Logger class
Keyboard::Keyboard(ParticleSystem * theSystem, Logger * logger)
{
	particleSystem = theSystem;
	this -> logger = logger;
	for (int i = 0; i < NUMBER_KEYS; i++)
	{
		keys[i] = false;
	}
}

//This method responds to keyPressed events, taking the appropriate action for each key
//Parameter key - which key was pressed
void Keyboard::keyPressed(unsigned char key)
{
	if (!keys[key])  //Only respond once when a user presses a key - don't respond if it was already down beforehand
	{
		switch(key)
		{
		case 'W':
		case 'w':
			particleSystem -> increaseEarthGravity(1);
			break;
		case 's':
		case 'S':
			particleSystem -> increaseEarthGravity(-1);
			break;
		case 'a':
		case 'A':
			particleSystem -> increaseStraightRestLength(-0.1);
			break;
		case 'D':
		case 'd':
			particleSystem -> increaseStraightRestLength(0.1);
			break;
		case 'x':
		case 'X':
			particleSystem -> toggleInfoText();
			break;
		case 'z':
		case 'Z':
			particleSystem -> toggleRenderMode();
			break;
		case 'c':
		case 'C':
			particleSystem -> toggleAmbientMode();
			break;
		case ' ':
			particleSystem -> toggleAnimation();
			break;
		case 'e':
		case 'E':
			particleSystem -> goToRestLengthsForGravity(0.01); 
			break;
		case 'r':
		case 'R':
			particleSystem -> reset();
			break;
		case 'i':
		case 'I':
			particleSystem -> toggleImageRendering();
			break;
		case 'p':
		case 'P':
			logger -> isLogging = !logger -> isLogging;
			break;
		case 'O':
		case 'o':
			logger -> printConstraintDeltaV = !logger -> printConstraintDeltaV;
			break;
		}
	}
	keys[key] = true;

}

//This method responds to key released events
//It marks the respective keyboard key as not being pressed
//Parameter key - the key that was released
void Keyboard::keyReleased(unsigned char key)
{
	keys[key] = false;

}
