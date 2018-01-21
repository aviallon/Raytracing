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
void rotateSphere(Sphere* sphere, Vec point, Vec axis, double angle_deg, int r){

	double angle = PI*angle_deg/180;
	
	Vec er(0, 1, 1);
	
	Vec temp = er*r; // On crée un point situé à r de l'origine
	
	Vec nouv_ct = temp.rotate(angle, axis) + point; // On le fait ensuite tourner de "angle" atour de l'axe (sans oublier l'offset du point).
	
	sphere->ct = nouv_ct;
}

void rotateSphere2Axis(Sphere* sphere, Vec point, double theta, double phi, double r){
	theta = PI*theta/180;
	phi = PI*phi/180;
	
	
	Vec er(cos(theta), sin(theta)*sin(phi), sin(theta)*cos(phi)); // Vecteur unitaire en coordonnées sphériques
	
	sphere->ct = er*r + point;
}

void angleInc(double* angle, double i){
	if((*angle)>=360){
		*angle = 0;
	}
	*angle += i;
}

double angleMercure = 0;
double angleVenus = 0;
double angleTerre = 0;
double angleLune = 0;
double angleMars = 0;
double angleJupiter = 0;
double angleSaturne = 0;
double angleUranus = 0;
double angleNeptune = 0;

double angleTest1 = 0;
double angleTest2 = 0;

//float rotationCamera = 270;

float jour = 0;

const double fov = (WIDTH/2)/tan(15*PI/180);

Vec camera((double)(WIDTH/2), (double)(HEIGHT/2), -fov);

//Sphere cameraSphere = Sphere(camera, 5, Color(0, 0, 0));

void animate(Allegro* allegro, float FPS){
	World* world_ptr = (World*)allegro->getContext();
	
	
	Sphere* mercure = ((Sphere*)world_ptr->getObject(0));
	Sphere* venus = ((Sphere*)world_ptr->getObject(1));
	Sphere* terre = ((Sphere*)world_ptr->getObject(2));
	Sphere* lune = ((Sphere*)world_ptr->getObject(3));
	Sphere* mars = ((Sphere*)world_ptr->getObject(4));
	Sphere* jupiter = ((Sphere*)world_ptr->getObject(5));
	Sphere* saturne = ((Sphere*)world_ptr->getObject(6));
	Sphere* uranus = ((Sphere*)world_ptr->getObject(7));
	Sphere* neptune = ((Sphere*)world_ptr->getObject(8));
	
	//Sphere* test = ((Sphere*)world_ptr->getObject(9));
	
	/* Animation */
	
	const int vRef = 2; // coefficient pour accelerer le mouvement des planetes
	const int rRef = 83.23*world_ptr->light.r / 10; // on divise par 10 les rayons réels afin de pouvoir voir toutes les planètes
	//const int rayon_mercure = 57909176; // juste pour info
	Vec axeRotation = Vec(1, 0, 0);
	
	rotateSphere2Axis(mercure, world_ptr->light.ct, 90, 90-angleMercure, rRef); // Rotation sphérique :D
	rotateSphere(venus, world_ptr->light.ct, axeRotation, angleVenus, rRef*1.868597301);
	rotateSphere(terre, world_ptr->light.ct, axeRotation, angleTerre, rRef*2.583319214);
	rotateSphere(lune, terre->ct, axeRotation, angleLune, rRef*0.552537013);
	rotateSphere(mars, world_ptr->light.ct, axeRotation, angleMars, rRef*3.936105687);
	rotateSphere(jupiter, world_ptr->light.ct, axeRotation, angleJupiter, rRef*pow(13.441946178, 0.8)); // Sinon elle est trop loin
	rotateSphere(saturne, world_ptr->light.ct, axeRotation, angleSaturne, rRef*pow(24.54152986, 0.8));
	rotateSphere(uranus, world_ptr->light.ct, axeRotation, angleUranus, rRef*pow(49.675703933, 0.8));
	rotateSphere(neptune, world_ptr->light.ct, axeRotation, angleNeptune, rRef*pow(77.767358682, 0.8));
	
	//rotateSphere2Axis(test, world_ptr->light.ct, angleTest1, angleTest2, rRef*2);
	
	//angleInc(&angleTest1, 5/FPS);
	//angleInc(&angleTest2, 30/FPS);
//	rotateSphere(&cameraSphere, world_ptr->light.ct, rotationCamera, 100000);
//	camera = cameraSphere.ct;
	
	
	angleInc(&angleMercure, vRef*4.181917279/FPS); // vitesse angulaire en degré par jour
	angleInc(&angleVenus, vRef*1.602099239/FPS);
	angleInc(&angleTerre, vRef*0.985452303/FPS);
	angleInc(&angleLune, vRef*13.153801358/FPS);
	angleInc(&angleMars, vRef*0.522755161/FPS);
	angleInc(&angleJupiter, vRef*0.082992214/FPS);
	angleInc(&angleSaturne, vRef*0.03357874/FPS);
	angleInc(&angleUranus, vRef*0.011719041/FPS);
	angleInc(&angleNeptune, vRef*0.005968861/FPS);
//	
//	angleInc(&rotationCamera, 1/FPS);
	
	if(jour > 365){
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

Vec getVirtualScreenPixel(int x_dec, int y_dec, Vec camera, Allegro* allegro, int width_offset){
	World* world = (World*)allegro->getContext();
	
	const int lacet = world->lacet;
	const float angle_camera = PI*float(lacet)/180;
	
	//const float x = camera._x; // Hauteur de la caméra
	
	//const float z = fov*sin(angle_camera)+cos(angle_camera)+camera._z;
	//const float y = fov*cos(angle_camera)-sin(angle_camera)+camera._y;
	
	//Vec virtScreenCenter = Vec(x+world->offset_x, y+world->offset_y+width_offset, z+world->offset_z);
	
	const float distance_x = (allegro->getDisplayWidth()/2-x_dec);
	const float distance_y = (allegro->getDisplayHeight()/2-y_dec);
	
	const float angleGD = atan2(distance_x, fov)+angle_camera;
	
	const float angleHB = atan2(distance_y, fov);
	
	const float r_x = fov;//distance_x/(sin(angleGD))+fov; // distance de la camera au pixel concerné
	const float r_y = fov;//distance_y/(sin(angleHB))+fov;
	//const float r = sqrt(r_x*r_x + r_y*r_y);
	
	const float x_pixel = -r_x*sin(angleHB) + r_y*cos(angleHB) /*- camera._x*/; // Oui, c'est contre-intuitif
	
	const float z_pixel = r_x*sin(angleGD)*cos(angleHB) + r_y*sin(angleGD)*sin(angleHB) - allegro->getDisplayHeight()/2/*- camera._z*/;
	
	const float y_pixel = r_x*cos(angleGD)*cos(angleHB) + r_y*cos(angleGD)*cos(angleHB) /*- camera._y*/;
	
	return Vec(x_pixel, y_pixel, z_pixel);
}

void redraw(Allegro* allegro, float FPS)
{
	const int t = getms();
	
	allegro->lockScreen();
	
	const Color black(0,0,0);
	
	World* world_ptr = (World*)allegro->getContext();
	
	World world2 = *world_ptr;

	Sphere light = world2.light;
	
	camera = Vec((double)(WIDTH/2+world2.offset_x), (double)(HEIGHT/2+world2.offset_y), -fov+world2.offset_z);
	
	// Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+world2.lacet+width_offset, 0+world2.offset_z+world2.roulis)
	
//	PixelGrid pix = pixes[0];
//	for(int x=0;x<WIDTH;x++){
//		for(int y=0; y<HEIGHT; y++){
//			allegro->set_pixel(y, x, allegro->rgb(pix.get(x, y)._r, pix.get(x, y)._g, pix.get(x, y)._b));
//		}
//	}
	//cout << (string)getVirtualScreenPixel(0, 0, camera, allegro, width_offset) << ", ";
	//cout << (string)Vec(0+world2.offset_x, 0+world2.offset_y+width_offset, 0+world2.offset_z) << endl;
	//cout << fov << endl;
	Color pixel;
	if(!high_fps_mode && accumulate(renderTimeAverage.begin(), renderTimeAverage.end(), 0.0)/renderTimeAverage.size() > 100)
		high_fps_mode = true;
	//#pragma omp target teams distribute parallel for map(from:pixelgrid[0:100])
	for(int x=0;x<allegro->getDisplayHeight();x++){ // oui, il y a une inversion des axes...
		for(int y=0;y<allegro->getDisplayWidth();y++){
			//pixel = getPixelColor(x, y, world_ptr, camera);
			pixel = Color(100,100,100);
			Vec pI(INFINITY, INFINITY, INFINITY);
			Vec ptemp = light.intersect(camera, Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+width_offset, 0+world2.offset_z+world2.roulis));
			Impact impact(Vec(), 0);
			if(ptemp.nonVec != true){
				pI = ptemp;
				pixel = light.color/**(100000/(ptemp-camera).len())*/;
			}
			//}else{
			bool shadow = false;
			for(unsigned int i=0;i<world2.size();i++){
				double dt = 0;
				shadow = false;
				ptemp = world2.intersect(i, camera, Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+width_offset, 0+world2.offset_z+world2.roulis));
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
						dt = (L.normalize().dot(N)); // cos(L, N)
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
	camera = Vec((double)(WIDTH/2+world->offset_x), (double)(HEIGHT/2+world->offset_y), -fov+world->offset_z);
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
	world.light = Sphere(Vec(WIDTH/2, HEIGHT/2, 10000), 60, Color(255, 255, 150));
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
	
	//Jupiter
	world.addObject(Sphere(Vec(), 0.211828148/5*rayon_planetes, Color(255, 50, 50))); // Sinon elle est VRAIMENT trop grosse
	
	//Saturne
	world.addObject(Sphere(Vec(), 0.178571852/5*rayon_planetes, Color(80, 255, 255))); // Sinon elle est VRAIMENT trop grosse
	
	//Uranus
	world.addObject(Sphere(Vec(), 0.07573037/5*rayon_planetes, Color(0, 50, 220))); // Sinon elle est VRAIMENT trop grosse
	
	//Neptune
	world.addObject(Sphere(Vec(), 0.073374815/5*rayon_planetes, Color(0, 0, 255))); // Sinon elle est VRAIMENT trop grosse
	
	// Test Sphere
	//world.addObject(Sphere(Vec(), 20, Color(255, 255, 255)));
	
	//world.addObject(Plan(Vec(allegro->getDisplayHeight()+10, 0, 0), Vec(1, 0, 0), Color(255, 255, 255), true));
	
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
