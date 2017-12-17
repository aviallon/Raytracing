#pragma once
#ifndef ALLEGRO_H_
#define ALLEGRO_H_

#include <allegro5/allegro.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_primitives.h>
#include <allegro5/allegro_font.h>
#include <allegro5/allegro_ttf.h>
#include <allegro5/allegro_memfile.h>
#include <cmath>
#include <vector>
#define MAX(a, b) (((a > b))?(a):(b))
#define MIN(a, b) (((a < b))?(a):(b))

#include <exception>
#include <memory>

#include "mouse.h"

class Mouse;

typedef unsigned char uchar;

class Allegro
{
private:
    ALLEGRO_DISPLAY *display;
	ALLEGRO_BITMAP *display_bitmap;
    ALLEGRO_TIMER *timer;
    ALLEGRO_EVENT_QUEUE *event_queue;
	ALLEGRO_FONT *default_font;
	ALLEGRO_FILE *arial_file;
	
	void (*mouse_clicked_func_ptr)(Allegro*, void*, unsigned char, int, int);
	void (*mouse_moved_func_ptr)(Allegro*, void*, unsigned char, int, int);
	
	void (*key_down_func_ptr)(Allegro*, void*, unsigned char, uint8_t);
	void (*key_up_func_ptr)(Allegro*, void*, unsigned char, uint8_t);
	
	void (*redraw_func_ptr)(Allegro*, float);
	void (*animate_func_ptr)(Allegro*, float);
	
	void* context;
	
	Mouse mouse;
	
	std::vector<bool> keys;

    bool looping, redraw, redraw_paused;
	
	void _exec_mouse_clicked_function();
	void _exec_mouse_moved_function(uchar ev);
	
	void _exec_key_down_function(uint8_t keycode);
	void _exec_key_up_function(uint8_t keycode);
	void _exec_key_repeat_function();
	
	
	float m_FPS;
	bool cursorSticked = false;

public:
    Allegro();
    ~Allegro();

    int init();
    int createWindow(float FPS, int w, int h);
    void gameLoop();
	
	void bindMouseClick(void (*fptr)(Allegro*, void*, unsigned char, int, int));
	void bindMouseMove(void (*fptr)(Allegro*, void*, unsigned char, int, int));
	
	void bindKeyDown(void (*fptr)(Allegro*, void*, unsigned char, uint8_t));
	void bindKeyUp(void (*fptr)(Allegro*, void*, unsigned char, uint8_t));
	
	void setRedrawFunction(void (*fptr)(Allegro*, float));
	void setAnimateFunction(void (*fptr)(Allegro*, float));
	
	static const uchar MOUSE_R_CLICKED = 1<<0;
	static const uchar MOUSE_L_CLICKED = 1<<7;
	static const uchar MOUSE_MOVED = 1<<1;
	static const uchar KEY_DOWN = 1<<2;
	static const uchar KEY_UP = 1<<3;
	static const uchar KEY_REPEAT = 1<<4;
	static const uchar MOUSE_WHEELED = 1<<5;
	static const uchar MOUSE_MOVED_DELTA = 1<<6;
	
	struct ALLEGRO_COLOR rgb(int r, int g, int b);
	
	/* Les fonctions pour dessiner ! Enfin ! */
	
	void set_pixel(int x, int y, ALLEGRO_COLOR color);
	
	void draw_line(int x1, int y1, int x2, int y2, ALLEGRO_COLOR color, int width);
	void draw_line(int x1, int y1, int x2, int y2);
	
	void draw_ellipse(int x1, int y1, int x2, int y2, ALLEGRO_COLOR color, int width, bool filled);
	void draw_ellipse(int x1, int y1, int x2, int y2);
	
	void draw_ellipse_r(int cx, int cy, int rx, int ry, ALLEGRO_COLOR color, int width, bool filled);
	
	void draw_rectangle(int x1, int y1, int x2, int y2, ALLEGRO_COLOR color, int width, bool filled);
	void draw_rectangle(int x1, int y1, int x2, int y2);
	
	void draw_text(int x, int y, std::string text, ALLEGRO_COLOR color, ALLEGRO_FONT* font);
	void draw_text(int x, int y, std::string text, ALLEGRO_COLOR color);
	void draw_text(int x, int y, const char* text, ALLEGRO_COLOR color);
	
	
	// Locks the buffer. Screen is only updated after it as been unlocked
	void lockScreen();
	void unlockScreen();
	
	void setCursorVisibility(bool visible);
	void setStickCursorToCenter(bool stick);
	
	void stopRedraw();
	void resumeRedraw();
	
	void flipDisplay();
	
	long int getTime();
	
	int getMouseX();
	int getMouseY();
	
	bool isKeyDown(int keycode);
	
	void setContext(void* cont);
	void* getContext();
	
	int getDisplayWidth();
	int getDisplayHeight();
	
	void toggleFullscreen(bool activate);
	bool isInFullscreen();
	
	static void _undefined_(Allegro* master, void* context, unsigned char event, int x, int y);
	static void _undefined_(Allegro* master, void* context, unsigned char event, uint8_t keycode);
	static void _undefined_(Allegro* master, float FPS);
	static struct ALLEGRO_COLOR rgbS(int r, int g, int b){
		return al_map_rgb(r, g, b);
	}
};

#endif
