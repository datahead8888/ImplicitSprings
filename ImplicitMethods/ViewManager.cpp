#include <gl/glut.h>
#include <cmath>
#include <iostream>
#include "ViewManager.h"

//Constructor for ViewManager class
ViewManager::ViewManager()
{
	isTrackingLeft = false;
	isTrackingMiddle = false;
	isTrackingRight = false;
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		currentPosition[i] = 0;
		previousLeftPosition[i] = 0;
		previousMiddlePosition[i] = 0;
		previousRightPosition[i] = 0;
	}
	
	xAngle = yAngle = 0;
	xTranslation = yTranslation = 0;
	zoomLevel = 0;
	autoRotate = false;
	rotationSpeed = 100;

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
		if (!isTrackingLeft)
		{
			//For brand new clicks the previous and current position must be the same so that we don't end up with a giant delta position
			previousLeftPosition[0] = currentPosition[0] = x;
			previousLeftPosition[1] = currentPosition[1] = y;
		}
		isTrackingLeft = true;
	}
	else if (state == GLUT_UP && button == GLUT_LEFT_BUTTON)
	{
		isTrackingLeft = false;
	}
	
	if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)
	{
		//zoomLevel -= 0.5;
		if (!isTrackingMiddle)
		{
			//For brand new clicks the previous and current position must be the same so that we don't end up with a giant delta position
			previousMiddlePosition[0] = currentPosition[0] = x;
			previousMiddlePosition[1] = currentPosition[1] = y;
			isTrackingMiddle = true;
		}


	}
	else if (state == GLUT_UP && button == GLUT_MIDDLE_BUTTON)
	{
		isTrackingMiddle = false;
	}

	if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
	{
		if (!isTrackingRight)
		{
			//For brand new clicks the previous and current position must be the same so that we don't end up with a giant delta position
			previousRightPosition[0] = currentPosition[0] = x;
			previousRightPosition[1] = currentPosition[1] = y;
			isTrackingRight = true;
		}

	}
	else if (state == GLUT_UP && button == GLUT_RIGHT_BUTTON)
	{
		isTrackingRight = false;
	}
}

//This method responds to mouse move events
//Parameters:
//x - x coordinate of movement
//y - y coordinate of movement
void ViewManager::mouseMove(int x, int y)
{
  //If a button is not being held, do nothing
  //if (!isTrackingLeft)
  //{
  //	  return;
  //}

  currentPosition[0] = x;
  currentPosition[1] = y;

  if (isTrackingLeft)
  {
	double deltaPosition[MOUSE_DIMENSION];  //Represents amount of change in position in x and y directions
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		deltaPosition[i] = currentPosition[i] - previousLeftPosition[i];
	}

	xAngle += deltaPosition[0];
	yAngle += deltaPosition[1];
  
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		  previousLeftPosition[i] = currentPosition[i];
	}
  }

  if (isTrackingMiddle)
  {
	double deltaPosition[MOUSE_DIMENSION];  //Represents amount of change in position in x and y directions
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		deltaPosition[i] = currentPosition[i] - previousMiddlePosition[i];
	}

	xTranslation -= deltaPosition[0] * 0.5;
	yTranslation += deltaPosition[1] * 0.5;
  
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		  previousMiddlePosition[i] = currentPosition[i];
	}
  }

  if (isTrackingRight)
  {
	double deltaPosition[MOUSE_DIMENSION];  //Represents amount of change in position in x and y directions
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		deltaPosition[i] = currentPosition[i] - previousRightPosition[i];
	}

	zoomLevel += deltaPosition[0];
	zoomLevel += deltaPosition[1];
  
	for (int i = 0; i < MOUSE_DIMENSION; i++)
	{
		  previousRightPosition[i] = currentPosition[i];
	}
  }


}

//This method performs the transformation for the current view
//Applies the look at angle, the zoom, and the mouse rotation
glm::mat4 ViewManager::doTransform()
{ 
  //gluLookAt(0.0f, 0.0f, 19.0f, 0.0f, 0.0f, 18.0f, 0.0f, 1.0f, 0.0f);
  //glTranslatef(0.0f, 0.0f, -zoomLevel);
  //glRotatef(xAngle, 0.0f, 1.0f, 0.0f);
  //glRotatef(yAngle, 1.0f, 0.0f, 0.0f);

	//Have the camera look directly at the mesh
	//Rotate around him, or zoom in / out - based on the mouse camera settings
	glm::vec3 eyePos(0.0f, 0.0f, 19.0f);
	glm::vec3 interestPoint(0.0f, 0.0f, 18.0f);
	glm::mat4 viewingMatrix = glm::lookAt(eyePos, interestPoint, glm::vec3(0.0, 1.0, 0.0));

	glm::mat4 transformMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -zoomLevel));
    transformMatrix = glm::translate(transformMatrix, glm::vec3(xTranslation, yTranslation, 0.0f));
	transformMatrix = glm::rotate(transformMatrix, xAngle, glm::vec3(0.0f, 1.0f, 0.0f));
	transformMatrix = glm::rotate(transformMatrix, yAngle, glm::vec3(1.0f, 0.0f, 0.0f));
	//transformMatrix = glm::scale(transformMatrix, glmvec3(1.0f, scale_size, 1.0f));

	glm::mat4 modelViewMatrix = viewingMatrix * transformMatrix;
	return modelViewMatrix;
	
}

//Do update logic for view manager
//Currently makes scene auto rotate if this option is turned on
void ViewManager::doUpdate(double timeElapsed)
{
	if (autoRotate)
	{
		xAngle += rotationSpeed * timeElapsed;
	}

}