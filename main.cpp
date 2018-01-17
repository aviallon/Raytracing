#include <iostream>
#include <allegro5/allegro.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_primitives.h>
#include <cmath>
#include <numeric>

#ifndef STRING_H_
#define STRING_H_
#include <string.h>
#include <sstream>
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

//#define WIDTH 500
//#define HEIGHT 500

const signed int WIDTH = 500;
const signed int HEIGHT = 500;

signed int width_offset = 0;

bool high_fps_mode = false;

std::vector<int> renderTimeAverage = std::vector<int>(5);

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

//int8_t sens = 1;
//int8_t sens2 = 2;
void rotateSphere(Sphere* sphere, Vec axis, float angle_deg, int r){
	float angle = PI*float(angle_deg)/180;
	//Vec Er = Vec(cos(angle), sin(angle), 0.0);
	//Vec Etheta = Vec(-sin(angle), cos(angle), 0.0);
	//Vec Ez = Vec(0, 0, 1);
	//const int x = sphere->ct._x;
	//const int y = sphere->ct._y;
	const float z = axis._x;
	
	float x = r*sin(angle)+cos(angle)+axis._z;
	float y = r*cos(angle)-sin(angle)+axis._y;
	
	Vec nouv_ct = Vec(z, y,x );
	
	sphere->ct = nouv_ct;
}

void angleInc(float* angle, float i){
	if((*angle)>=360){
		*angle = 0;
	}
	*angle += i;
}

float angleMercure = 0;
float angleVenus = 0;
float angleTerre = 0;
float angleLune = 0;
float angleMars = 0;

float jour = 0;

void animate(Allegro* allegro, float FPS){
	World* world_ptr = (World*)allegro->getContext();
	
	
	Sphere* mercure = ((Sphere*)world_ptr->getObject(0));
	Sphere* venus = ((Sphere*)world_ptr->getObject(1));
	Sphere* terre = ((Sphere*)world_ptr->getObject(2));
	Sphere* lune = ((Sphere*)world_ptr->getObject(3));
	Sphere* mars = ((Sphere*)world_ptr->getObject(4));
	
	/* Animation */
	
	/*if(sp1->ct._x > WIDTH-4*sp1->r && sens != -1)
		sens = -1;
	if(sp1->ct._x < sp1->r && sens != 1)
		sens = 1;
	sp1->ct._x += sens;*/
	const int vRef = 2; // coefficient pour accelerer le mouvement des planetes
	const int rRef = 83.23*world_ptr->light.r / 10; // on divise par 10 les rayons réels afin de pouvoir voir toutes les planètes
	//const int rayon_mercure = 57909176; // juste pour infon
	
	rotateSphere(mercure, world_ptr->light.ct, angleMercure, rRef);
	rotateSphere(venus, world_ptr->light.ct, angleVenus, rRef*1.868597301);
	rotateSphere(terre, world_ptr->light.ct, angleTerre, rRef*2.583319214);
	rotateSphere(lune, terre->ct, angleLune, rRef*0.552537013);
	rotateSphere(mars, world_ptr->light.ct, angleMars, rRef*3.936105687);
	
	angleInc(&angleMercure, vRef*26.275761202/FPS); // vitesse angulaire en degré par jour
	angleInc(&angleVenus, vRef*10.066286396/FPS);
	angleInc(&angleTerre, vRef*0.985452303/FPS);
	angleInc(&angleLune, vRef*82.647771429/FPS);
	angleInc(&angleMars, vRef*3.284567544/FPS);
	
	if(jour >= 365){
		jour = 0;
	}
	jour += vRef/FPS;
	//cout << angle << endl;
	
	/*if(sp2->ct._y > 100 && sens2 != -1)
		sens2 = -1;
	if(sp2->ct._y < 0 && sens != 1)
		sens2 = 1;
	sp2->ct._y += sens2;*/
	
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
	//if(!high_fps_mode && accumulate(renderTimeAverage.begin(), renderTimeAverage.end(), 0.0)/renderTimeAverage.size() > 100)
	//	high_fps_mode = true;
	//#pragma omp target teams distribute parallel for map(from:pixelgrid[0:100])
	for(int x=0;x<allegro->getDisplayHeight();x++){ // oui, il y a une inversion des axes...
		for(int y=0;y<allegro->getDisplayWidth();y++){
			//pixel = getPixelColor(x, y, world_ptr, camera);
			pixel = Color(10,10,10);
			Vec pI(INFINITY, INFINITY, INFINITY);
			Vec ptemp = light.intersect(camera, Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+world2.lacet+width_offset, 0+world2.offset_z+world2.roulis));
			Impact impact(Vec(), 0);
			if(ptemp.nonVec != true){
				pI = ptemp;
				pixel = light.color*(100000/(ptemp-camera).len());
			}
			//}else{
			bool shadow = false;
			for(unsigned int i=0;i<world2.size();i++){
				double dt = 0;
				shadow = false;
				ptemp = world2.intersect(i, camera, Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+world2.lacet+width_offset, 0+world2.offset_z+world2.roulis));
				if((ptemp.nonVec != true && (ptemp-camera).len() < (pI-camera).len())){
					pI = ptemp;
					Vec L = light.ct - pI;
					
					if(!high_fps_mode)
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
			allegro->set_pixel(y, x, allegro->rgb(pixel._r, pixel._g, pixel._b));
		}
	}
	
	allegro->unlockScreen();
	
	stringstream fps_disp;
	
	renderTimeAverage.push_back((getms() - t));
	renderTimeAverage.erase(renderTimeAverage.begin());
	
	fps_disp << fps(renderTimeAverage) << " FPS\0";
	
	allegro->draw_text(30, 10, fps_disp.str(), allegro->rgb(255, 255, 255));
	
	stringstream jourstr;
	jourstr << "Jour " << (int)jour;
	
	allegro->draw_text(30, 30, jourstr.str(), allegro->rgb(255, 255, 255));
	
	
	width_offset = (WIDTH - allegro->getDisplayWidth())/2;
	//cout << allegro->getDisplayHeight() << ", " << allegro->getDisplayWidth() << endl;
	//WIDTH = allegro->getDisplayWidth();
	//HEIGHT = allegro->getDisplayHeight();
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
				allegro->toggleFullscreen(false);
				break;
			case ALLEGRO_KEY_F:
				//cout << (allegro->isKeyDown(ALLEGRO_KEY_ALT) &&  allegro->isKeyDown(ALLEGRO_KEY_LCTRL)) << endl;
				//flush(cout);
				if(allegro->isKeyDown(ALLEGRO_KEY_ALT) &&  allegro->isKeyDown(ALLEGRO_KEY_LCTRL)){
					allegro->toggleFullscreen(true);
				}
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
    allegro->createWindow(60, HEIGHT, WIDTH);
	
	allegro->setStickCursorToCenter(true);
	allegro->setCursorVisibility(false);
	
	World world;
	width_offset = (WIDTH - allegro->getDisplayWidth())/2;
	world.light = Sphere(Vec(WIDTH/2, HEIGHT/2, 50), 50, Color(200, 255, 255));
	allegro->setContext(&world);

	const int rayon_planetes = 30*world.light.r;

	//Mercure
	world.addObject(Sphere(Vec(), 0.00701308*rayon_planetes, Color(100, 100, 100)));
	
	//Vénus
	world.addObject(Sphere(Vec(), 0.017396866*rayon_planetes, Color(200, 250, 200)));
	
	//Terre
	world.addObject(Sphere(Vec(), 0.018335489*rayon_planetes, Color(0, 50, 255)));
	
	//Lune
	world.addObject(Sphere(Vec(), 0.002496766*rayon_planetes, Color(100, 100, 100)));
	
	//Mars
	world.addObject(Sphere(Vec(), 0.009762829*rayon_planetes, Color(255, 100, 100)));
	//world.addObject(Sphere(Vec(300, 150, 80), 20, Color(0, 200, 255)));
	//world.addObject(Sphere(Vec(230, 180, -300), 5, Color(200, 10, 140)));
	//for(int bla = 0;bla<=10;bla++)
	//	world.addObject(Sphere(Vec(20, 600-10*bla, 150*bla), 7, Color(100, 20*bla, 200)));
	//world.addObject(Plan(Vec(WIDTH+50, 0, 0), Vec(1, 0, 0), Color(255, 255, 255), true));
	
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
