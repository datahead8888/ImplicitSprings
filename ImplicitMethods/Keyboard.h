#pragma once

#include "ParticleSystem.h"
#include "ViewManager.h"
#define NUMBER_KEYS 256						//Number of keys for which this class tracks state

//This class both tracks whether or not each key is being held down and responds to keypresses, calling
//methods of the appropriate classes
class Keyboard
{
private:
	bool keys[NUMBER_KEYS];					//Boolean array for all possible keyboard keys.  True means a key is held; false means it is not.
	ParticleSystem * particleSystem;		//Reference to the particle system so that control methods can be called
	Logger * logger;						//Reference to the logger class to allow logging to the console
	ViewManager * viewManager;				//Reference to the view manager
public:
	Keyboard(ParticleSystem * theSystem, ViewManager * viewManager, Logger * logger);
	void keyPressed(unsigned char key);
	void keyReleased(unsigned char key);
	void setViewManager(ViewManager * viewManager) {this->viewManager = viewManager;}
};

