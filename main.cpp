#include <iostream>
#include <allegro5/allegro.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_primitives.h>
#include <cmath>

#ifndef STRING_H_
#define STRING_H_
#include <string.h>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#endif

//#include "avlib.h"
#include <typeinfo>
#include <vector>
#include <chrono>
#include <thread>
#include <mutex>
#include <future>
#define MAX(a, b) (((a > b))?(a):(b))
#define MIN(a, b) (((a < b))?(a):(b))
#define PI 3.14159265359

#define WIDTH 500
#define HEIGHT 500

#include "allegro.h"
#include "raytracer.h"


#define TRANSLATE 0

using namespace std;

std::string toString(size_t entier);

//std::string toString2(size_t entier){
//	return boost::lexical_cast<std::string>(entier);
//}

std::string toStringF(float f){
	return boost::lexical_cast<std::string>(f);
}


ALLEGRO_COLOR colorToAllegro(Color color){
	return Allegro::rgbS(color._r, color._g, color._b);
}

int8_t sens = 1;

void animate(Allegro* allegro, float FPS){
	World* world_ptr = (World*)allegro->getContext();
	Sphere* sp1 = ((Sphere*)world_ptr->getObject(0));
	
	/* Animation */
	
	if(sp1->ct._x > WIDTH-4*sp1->r && sens != -1)
		sens = -1;
	if(sp1->ct._x < sp1->r && sens != 1)
		sens = 1;
	sp1->ct._x += sens;
	
	/* Fin animation */
}

void redraw(Allegro* allegro, float FPS)
{
	const int t = getms();
	
	allegro->lockScreen();
	
	const Color black(0,0,0);
	
	World* world_ptr = (World*)allegro->getContext();
	
	World world2 = *world_ptr;
	
	Vec camera((double)(WIDTH/2+world2.offset_x), (double)(HEIGHT/2+world2.offset_y), -(WIDTH/2)/tan(15*PI/180)+world2.offset_z);

	Sphere light = world2.light;
	
	//thread t1(raytrace, 0, WIDTH, allegro, world_ptr, camera);
	//t1.join();
	//this_thread::sleep_for(10);
//	vector<future<PixelGrid>> rets(4);
//	
//	for(int i=0;i<rets.size();i++){
//		rets[i] = std::async(&raytrace, i*WIDTH/rets.size(), (i+1)*WIDTH/rets.size(), world_ptr, camera);
//	}
//	
//	//std::future<PixelGrid> ret1 = std::async(&raytrace, 0, WIDTH/2, world_ptr, camera);
//	
//	//std::future<PixelGrid> ret2 = std::async(&raytrace, WIDTH/2, WIDTH, world_ptr, camera);
//	
//	vector<PixelGrid> pixes(4);
//	
//	for(int i=0;i<rets.size();i++){
//		pixes[i] = rets[i].get();
//	}
//	
//	PixelGrid pix = pixes[0];
//	for(int x=0;x<WIDTH;x++){
//		for(int y=0; y<HEIGHT; y++){
//			allegro->set_pixel(y, x, allegro->rgb(pix.get(x, y)._r, pix.get(x, y)._g, pix.get(x, y)._b));
//		}
//	}
	Color pixel;
	for(int x=0;x<WIDTH;x++){
		for(int y=0;y<HEIGHT;y++){
			//pixel = getPixelColor(x, y, world_ptr, camera);
			pixel = Color(10,10,10);
			Vec pI(INFINITY, INFINITY, INFINITY);
			Vec ptemp = light.intersect(camera, Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+world2.lacet, 0+world2.offset_z+world2.roulis));
			Impact impact(Vec(), 0);
			if(ptemp.nonVec != true){
				pixel = light.color*(1000/(ptemp-camera).len());
			}else{
				for(unsigned int i=0;i<world2.size();i++){
					double dt = 0;
					bool shadow = false;
					ptemp = world2.intersect(i, camera, Vec(x+world2.offset_x, y+world2.offset_y, 0+world2.offset_z));
					if((ptemp.nonVec != true && (ptemp-camera).len() < (pI-camera).len())){
						pI = ptemp;
						Vec L = light.ct - pI;
						
						shadow = shadowRay(pI, L, world2, light, i);
						//const bool shadow = false;
						
						if(shadow){
							dt = 0.01f;
						} else {
							Vec N = world2.getNormale(i, pI);
							dt = (L.normalize().dot(N));
						}
						
						pixel = world2.getColor(i, pI).mix(light.color)*dt;
					}
				}
			}
			allegro->set_pixel(y, x, allegro->rgb(pixel._r, pixel._g, pixel._b));
		}
	}
	
	allegro->unlockScreen();
	
	stringstream fps_disp;
	fps_disp << fps(t, getms()) << " FPS\0";
	
	allegro->draw_text(30, 10, fps_disp.str(), allegro->rgb(255, 255, 255));
	
}

void mouseMove(Allegro* allegro, void* context, unsigned char event, int x, int y){
	World* world = (World*)context;
	if(event == Allegro::MOUSE_WHEELED)
		world->offset_x += 2*x;
	else if(event == Allegro::MOUSE_MOVED_DELTA){
		//world->lacet += y/2;
		world->offset_z -= y/2;
		world->offset_y += x/2;
	}
}

void mouseClick(Allegro* allegro, void* context, unsigned char event, int x, int y){
	if(event == Allegro::MOUSE_L_CLICKED){
		allegro->setCursorVisibility(false);
		allegro->setStickCursorToCenter(true);
	}
}

void move(Allegro* allegro, void* context, unsigned char event, uint8_t keycode){
	World* world = (World*)context;
	if(event == Allegro::KEY_DOWN){
		switch(keycode){
			case ALLEGRO_KEY_ESCAPE:
				allegro->setStickCursorToCenter(false);
				allegro->setCursorVisibility(true);
				break;
			case ALLEGRO_KEY_UP:
				world->offset_z += 5;
				break;
			case ALLEGRO_KEY_DOWN:
				world->offset_z -= 5;
				break;
			case ALLEGRO_KEY_RIGHT:
				world->offset_y += 5;
				break;
			case ALLEGRO_KEY_LEFT:
				world->offset_y -= 5;
				break;
			case ALLEGRO_KEY_D:
				world->lacet += 10;
				break;
			case ALLEGRO_KEY_Q:
				world->lacet -= 10;
				break;
			case ALLEGRO_KEY_S:
				world->roulis += 10;
				break;
			case ALLEGRO_KEY_Z:
				world->roulis -= 10;
				break;
				
		}
	} else {
		if(allegro->isKeyDown(ALLEGRO_KEY_UP))
			world->offset_z += 5;
		
		if(allegro->isKeyDown(ALLEGRO_KEY_DOWN))
			world->offset_z -= 5;
			
		if(allegro->isKeyDown(ALLEGRO_KEY_RIGHT))
			world->offset_y += 5;
			
		if(allegro->isKeyDown(ALLEGRO_KEY_LEFT))
			world->offset_y -= 5;
	}
}


int main(int argc, char **argv)
{
	Allegro allegro_obj = Allegro();
	Allegro* allegro = &allegro_obj;
    allegro->init();
    allegro->createWindow(30, HEIGHT, WIDTH);
	
	allegro->setStickCursorToCenter(true);
	allegro->setCursorVisibility(false);
	
	World world;
	world.light = Sphere(Vec(150, HEIGHT/2-80, 20), 5, Color(255, 255, 255));
	allegro->setContext(&world);

	world.addObject(Sphere(Vec(WIDTH/2, HEIGHT/2, 10), 20, Color(255, 255, 5)));
	world.addObject(Sphere(Vec(WIDTH/2-30, HEIGHT/2-100, 50), 30, Color(0, 255, 100)));
	world.addObject(Sphere(Vec(210, 220, 10), 6, Color(255, 255, 255)));
	world.addObject(Sphere(Vec(300, 150, 80), 20, Color(0, 200, 255)));
	world.addObject(Sphere(Vec(230, 180, -300), 5, Color(200, 10, 140)));
	//for(int bla = 0;bla<=10;bla++)
	//	world.addObject(Sphere(Vec(20, 600-10*bla, 150*bla), 7, Color(100, 20*bla, 200)));
	world.addObject(Plan(Vec(WIDTH+50, 0, 0), Vec(1, 0, 0), Color(255, 255, 255), true));
	
	allegro->bindKeyDown(move);
	allegro->bindMouseMove(mouseMove);
	allegro->bindMouseClick(mouseClick);

	allegro->setRedrawFunction(redraw);
	allegro->setAnimateFunction(animate);
	
    allegro->gameLoop();
	
	allegro->setStickCursorToCenter(false);
	allegro->setCursorVisibility(true);
	
	delete allegro;

    return 0;
}
