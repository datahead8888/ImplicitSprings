#pragma once

const int MOUSE_DIMENSION = 2;	//Number of dimensions for processing mouse clicks - 2 for x and y

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

//This class manages the current view
//Currently it manages
//1) Holding the left button and dragging left, right, up, and down - this rotates the grid object
//2) Pressing middle button (zoom in) or right button (zoom out)
class ViewManager
{
public:
	ViewManager();
	void mouseClick(int button, int state, int x, int y);
	void mouseMove(int x, int y);
	glm::mat4 doTransform();
	void toggleAutoRotate() {autoRotate = !autoRotate;}
	void doUpdate(double timeElapsed);
	
private:
	bool isTrackingLeft;							//True if left button is held so that movement turns the object
	bool isTrackingMiddle;							//True if middle button is held so that movement zooms the object
	bool isTrackingRight;							//True if right button is held so that movement zooms the object
	double previousLeftPosition[MOUSE_DIMENSION];	//Coordinates of last stored left mouse button location
	double previousMiddlePosition[MOUSE_DIMENSION];	//Coordinates of last stored middle mouse button location
	double previousRightPosition[MOUSE_DIMENSION];	//Coordinates of last stored right mouse button location
	double currentPosition[MOUSE_DIMENSION];	//Coordinates of latest mouse click
	float yAngle;								//Amount of rotation due to vertical mouse movements
	float xAngle;								//Amount of rotation due to horizontal mouse movements
	float xTranslation;							//Amount of translation due to horizontal mouse movements
	float yTranslation;							//Amount of translation due to vertical mouse movements
	float zoomLevel;							//Amount of zoom that occurred
	bool autoRotate;							//True if automatically spins, false if not
	double rotationSpeed;						//Speed at which scene rotates when auto rotate is turned on
};