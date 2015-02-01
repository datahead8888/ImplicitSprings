#pragma once

const int MOUSE_DIMENSION = 2;	//Number of dimensions for processing mouse clicks - 2 for x and y

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
	void doTransform();
	
private:
	bool isTracking;							//True if left button is held so that movement turns the object
	double previousPosition[MOUSE_DIMENSION];	//Coordinates of last stored mouse location
	double currentPosition[MOUSE_DIMENSION];	//Coordinates of latest mouse click
	double yAngle;								//Amount of rotation due to vertical mouse movements
	double xAngle;								//Amount of rotation due to horizontal mouse movements
	double zoomLevel;							//Amount of zoom that occurred
};