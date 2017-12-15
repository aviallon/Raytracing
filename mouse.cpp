#include "mouse.h"
#include <allegro5/mouse.h>

Mouse::Mouse(){
	state = new struct ALLEGRO_MOUSE_STATE;
	state->x = 0;
	state->y = 0;
}

Mouse::~Mouse(){
	delete state;
}

int Mouse::getX(){
	return state->x;
}

int Mouse::getY(){
	return state->y;
}

int Mouse::getDZ(){
	return mouseDZ;
}

void Mouse::setDZ(int dz){
	mouseDZ = dz;
}

int Mouse::getDX(){
	return mouseDX;
}

int Mouse::getDY(){
	return mouseDY;
}

void Mouse::setDX(int dx){
	mouseDX = dx;
}

void Mouse::setDY(int dy){
	mouseDY = dy;
}

int Mouse::getBtn(){
	return btn;
}

void Mouse::setBtn(int btn){
	this->btn = btn;
}

struct ALLEGRO_MOUSE_STATE* Mouse::getStatePtr(){
	return state;
}
