#include <gl/glut.h>
#include <cmath>
#include <iostream>
#include "ViewManager.h"

//Constructor for ViewManager class
ViewManager::ViewManager()
{
	isTracking = false;
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		currentPosition[i] = 0;
		previousPosition[i] = 0;
	}
	
	xAngle = yAngle = 0;
	zoomLevel = 0;

}

//This method responds to mouse click events
//Parameters:
//button - which button had an event
//state - was it down or up
//x - x coordinate of click
//y - y coordinate of click
void ViewManager::mouseClick(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)
	{
		if (!isTracking)
		{
			//For brand new clicks the previous and current position must be the same so that we don't end up with a giant delta position
			previousPosition[0] = currentPosition[0] = x;
			previousPosition[1] = currentPosition[1] = y;
		}
		isTracking = true;
	}
	else if (state == GLUT_UP && button == GLUT_LEFT_BUTTON)
	{
		isTracking = false;
	}
	else if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)
	{
		zoomLevel -= 0.5;

	}
	else if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
	{
		zoomLevel += 0.5;

	}
}

//This method responds to mouse move events
//Parameters:
//x - x coordinate of movement
//y - y coordinate of movement
void ViewManager::mouseMove(int x, int y)
{
  //If a button is not being held, do nothing
  if (!isTracking)
  {
	  return;
  }

  currentPosition[0] = x;
  currentPosition[1]= y;

  double deltaPosition[MOUSE_DIMENSION];  //Represents amount of change in position in x and y directions
  for (int i = 0; i < MOUSE_DIMENSION; i++)
  {
	  deltaPosition[i] = currentPosition[i] - previousPosition[i];
  }

  xAngle += deltaPosition[0];
  yAngle += deltaPosition[1];
  
  for (int i = 0; i < MOUSE_DIMENSION; i++)
  {
	  previousPosition[i] = currentPosition[i];
  }
}

//This method performs the transformation for the current view
//Applies the look at angle, the zoom, and the mouse rotation
void ViewManager::doTransform()
{ 
  gluLookAt(0.0f, 0.0f, 10.0f, 0.0f, 0.0f, 9.0f, 0.0f, 1.0f, 0.0f);
  glTranslatef(0.0f, 0.0f, -zoomLevel);
  glRotatef(xAngle, 0.0f, 1.0f, 0.0f);
  glRotatef(yAngle, 1.0f, 0.0f, 0.0f);
}