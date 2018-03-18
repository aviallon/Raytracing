#include <iostream>
#include <allegro5/allegro.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_primitives.h>
#include <cmath>
#include <ctime>
#include <numeric>

#ifndef STRING_H_
#define STRING_H_
#include <string.h>
#include <sstream>
#include <boost/lexical_cast.hpp>
#endif

#include <typeinfo>
#include <vector>
#include <chrono>
#include <thread>
#include <mutex>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <future>

#define MAX(a, b) (((a > b))?(a):(b))
#define MIN(a, b) (((a < b))?(a):(b))
#define PI 3.14159265359

const signed int WIDTH = 500; // Taille **initiale** de la fenêtre
const signed int HEIGHT = 500;

/* 07/03/2009
 * M	V	T	Ma	J	S	U	N
 * 281	153	207	307	308	168	352	324
 */

double vRef = 2; // coefficient pour accelerer le mouvement des planetes

short selector = -1;

bool high_fps_mode = false;

std::vector<int> renderTimeAverage = std::vector<int>(5);

#include "allegro/allegro.h"
#include "raytracer.h"
#include "phys.hpp"
#include "math.hpp"


using namespace std;

std::string toString(size_t entier);

std::string toStringF(float f){
	return boost::lexical_cast<std::string>(f);
}


void rotateSphere(Sphere* sphere, Vec point, Vec axis, double angle_deg, double p, double excentrisme = 0){

	double angle = PI*angle_deg/180;
	
	double r = p/(1+ excentrisme*cos(angle));
	
	Vec temp = (Vec(1, 1, 1) - axis) * r; // On crée un point situé à r de l'origine
	
	Vec nouv_ct = temp.rotate(angle, axis) + point; // On le fait ensuite tourner de "angle" atour de l'axe (sans oublier l'offset du point).
	
	sphere->ct = nouv_ct;
}

void rotateSphere2Axis(Sphere* sphere, Vec point, double theta, double phi, double r){
	theta = PI/2 -PI*theta/180;
	phi = PI/2 -PI*phi/180;
	
	
	Vec er(cos(theta), sin(theta)*sin(phi), sin(theta)*cos(phi)); // Vecteur unitaire en coordonnées sphériques
	
	sphere->ct = er*r + point;
}

void angleInc(double* angle, double i){
	if((*angle)>=360){
		*angle = 0;
	} else if((*angle) < 0){
		*angle = 360;
	}
	*angle += i;
}

double angleMercure = 281;
double angleVenus = 156.37214484;
double angleTerre = 166;
double angleLune = 142.7187;
double angleMars = 307;
double angleJupiter = 308;
double angleSaturne = 168;
double angleUranus = 352;
double angleNeptune = 324;

double angleTest1 = 0;
double angleTest2 = 0;

//float rotationCamera = 270;

double temps = 0;

const double fov = (WIDTH/2)/tan(15*PI/180);


int reflectTest;
int turnTest;

double tTestProfondeur = 0;
std::vector<int> theTest;
//Sphere cameraSphere = Sphere(camera, 5, Color(0, 0, 0));

int physObjectSphereIndex = 0;
int physObjectIndex = 0;
int droneTagIndex = 0;

void animate(Allegro* allegro, float FPS){
	World* world_ptr = (World*)allegro->getContext();
	
	Sphere* physObjectSphere = ((Sphere*)world_ptr->getObject(physObjectSphereIndex));
	
	PhysicObject* physObject = world_ptr->getPhysicObject(physObjectIndex);
	
	Text3D* droneTag = world_ptr->getTextTag(droneTagIndex);
	
	
	//Sphere soleil = world_ptr->lights[0];
	
	//Sphere* second_light = &(world_ptr->lights[1]);
	
	/* Animation */
	
	chrono::time_point<high_resolution_clock> t = chrono::high_resolution_clock::now();
	duration<double, std::ratio<1> > dt = t - physObject->t0;
	
	temps = dt.count();
	

//Color(127*sin(angleTest1 + 0) + 128, 127*sin(angleTest1 + 2*PI/3) + 128, 127*sin(angleTest1 + 4*PI/3) + 128);
	Vec dronePos = physObject->getPos();
	droneTag->setPos(dronePos+Vec(-1.8, 0, 0));
	stringstream tag;
	tag << physObject->rpm << " RPM";
	droneTag->setText(string(tag.str()));
	
	physObjectSphere->ct = dronePos;
	
}

Vec screenPixRotate(int x, int y, double angle, Vec axis, Vec camera, Allegro* allegro){
	
	double lacetRad = PI*angle / 180 + PI/2;
	
	double dist = fov;
	Vec temp = Vec(0, 1, 1)*dist;
	
	Vec rotation = temp.rotate2(lacetRad, 0, Vec(0, 0, 1));
	
	Vec pos = Vec(x-allegro->getDisplayWidth()/2, y-allegro->getDisplayHeight()/2, 0);
	Vec new_pos = pos.rotate(lacetRad, Vec(0, 0, 1)) + rotation + camera;
	
//	angle = PI*angle/180;
//	return Vec(x, y, 0).rotate(angle, axis) + camera;
	return new_pos;
}

void disableHPerf(Allegro* allegro, Button* btn){
	allegro->getGUI()->displayMessage("High performance mode deactivated", 3000);
	high_fps_mode = false;
	
	allegro->getGUI()->eraseBtn(btn);
}

string getDate(double temps){
	using namespace boost::gregorian;
	stringstream dtstream;
	date t0 = date(2009, 03, 7);
	date now = (t0 + days((int)temps));
	dtstream << now.day() << "/" << now.month() << "/" << now.year();
	return dtstream.str();
}

string getDuration(double duration){
	stringstream dtstream;
	if(abs(duration) < 1.0){
		double hours = (24*duration);
		double minutes = abs(hours - (int)hours)*60;
		int seconds = (minutes - (int)minutes)*60;
		dtstream << (int)hours << ":" << (int)minutes << ":" << seconds;
	} else {
		int days = (int)duration;
		int hours = abs(duration - (int)duration)*24;
		dtstream << days << "j " << hours << "h";
	}
	return dtstream.str();
}

void redraw(Allegro* allegro, float FPS)
{
	const int t = getms();
	
	const Color black(0,0,0);
	
	World* world_ptr = (World*)allegro->getContext();
	
	World world2 = *world_ptr;

	world_ptr->camera = Vec((double)(WIDTH/2+world2.offset_x), (double)(HEIGHT/2+world2.offset_y), -fov+world2.offset_z);
	
	//world_ptr->lights[2].ct = world_ptr->camera;
	
	
	//Vec pix = screenPixRotate(x, y, world2.lacet, Vec(0, 1, 0), camera, allegro);
	// Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+world2.lacet+width_offset, 0+world2.offset_z+world2.roulis)
	
	if(!high_fps_mode && accumulate(renderTimeAverage.begin(), renderTimeAverage.end(), 0.0) / renderTimeAverage.size() > 100){ // If FPS are under 1000/_100_ = 10
		high_fps_mode = true;
		allegro->getGUI()->displayMessage("High performance mode activated", 3000);
		allegro->getGUI()->newBtn(string("Disable HPerf"), allegro->getDisplayWidth()/2, allegro->getDisplayHeight()-35, 30, 100, &disableHPerf);
	}
	
	
	vector<vector<Color> > screen(allegro->getDisplayHeight(), vector<Color>(allegro->getDisplayWidth()));
	
	#pragma omp parallel for collapse(2)
	for(int x=0;x<allegro->getDisplayHeight();x++){ // oui, il y a une inversion des axes...
		for(int y=0;y<allegro->getDisplayWidth();y++){
			Vec pix = Vec(x+world2.offset_x, y+world2.offset_y+world2.width_offset, 0+world2.offset_z);
			screen[x][y] = raytrace(world_ptr, world_ptr->camera, pix);
		}
	}
	
	allegro->lockScreen();
	
	for(int x=0;x<allegro->getDisplayHeight();x++){ // oui, il y a une inversion des axes...
		for(int y=0;y<allegro->getDisplayWidth();y++){
			allegro->set_pixel(y, x, allegro->rgb(screen[x][y]._r, screen[x][y]._g, screen[x][y]._b));
		}
	}
	
	allegro->unlockScreen();
	
	
	world_ptr->drawTextTag(droneTagIndex);
	
	stringstream fps_disp;
	
	renderTimeAverage.push_back((getms() - t));
	renderTimeAverage.erase(renderTimeAverage.begin());
	
	fps_disp << fps(renderTimeAverage) << " FPS\0";
	
	allegro->draw_text(5, 10, fps_disp.str(), allegro->rgb(255, 255, 255), ALLEGRO_ALIGN_LEFT);
	
	stringstream tmpStr;
	tmpStr << round(temps*1000)/1000 << " s";
	
	allegro->draw_text(5, 30, tmpStr.str(), allegro->rgb(255, 255, 255), ALLEGRO_ALIGN_LEFT);
	
	
	PhysicObject* drone = world_ptr->getPhysicObject(physObjectIndex);
	
	stringstream alt;
	alt << round(drone->pos._y*100)/100 << " m";
	
	allegro->draw_text(5, 50, alt.str(), allegro->rgb(255, 255, 255), ALLEGRO_ALIGN_LEFT);
	
	
	stringstream speed;
	speed << round(drone->speed._y*100)/100 << " m/s";
	
	allegro->draw_text(5, 70, speed.str(), allegro->rgb(255, 255, 255), ALLEGRO_ALIGN_LEFT);
	
	stringstream corr;
	corr << "Kc1 : " << round(drone->Kc1*10)/10;
	
	allegro->draw_text(5, 90, corr.str(), allegro->rgb(255, 255, 255), ALLEGRO_ALIGN_LEFT);
	
	world_ptr->width_offset = (WIDTH - allegro->getDisplayWidth())/2; // Required when resizing the display in order to move the camera accordingly
}

void mouseMove(Allegro* allegro, void* context, uint16_t event, int x, int y){
	World* world = (World*)context;
	if(event & Allegro::MOUSE_WHEELED)
		world->offset_x += 3*x;
	else if(event & Allegro::MOUSE_MOVED_DELTA){
		//world->lacet += y/2;
		world->offset_z -= y/2;
		world->offset_y += x/2;
	}
}

void mouseClick(Allegro* allegro, void* context, uint16_t event, int x, int y){
	if(event & Allegro::MOUSE_L_CLICKED){
		allegro->setCursorVisibility(false);
		allegro->setStickCursorToCenter(true);
	}
}

void move(Allegro* allegro, void* context, uint16_t event, uint8_t keycode){
	World* world = (World*)context;
	if(event == Allegro::KEY_DOWN){
		switch(keycode){
			case ALLEGRO_KEY_ESCAPE:
				allegro->setStickCursorToCenter(false);
				allegro->setCursorVisibility(true);
				allegro->toggleFullscreen(false);
				break;
			case ALLEGRO_KEY_F:
				if(allegro->isKeyDown(ALLEGRO_KEY_ALT) &&  allegro->isKeyDown(ALLEGRO_KEY_LCTRL)){
					allegro->toggleFullscreen(true);
				}
				break;

			if(allegro->isKeyDown(ALLEGRO_KEY_LSHIFT)){
				case ALLEGRO_KEY_RIGHT:
				{
					if(selector <= 2){
						selector += 1;
					} else {
						selector = 0;
					}
				}
				
				case ALLEGRO_KEY_LEFT:
				{
					if(selector > 0){
						selector -= 1;
					} else {
						selector = 2;
					}
				}
			} else {
				selector = -1;
			}

			case ALLEGRO_KEY_D:
				angleInc(&(world->lacet), 1);
				break;
			case ALLEGRO_KEY_Q:
				angleInc(&(world->lacet), -1);
				break;
			case ALLEGRO_KEY_S:
				angleInc(&(world->roulis), 1);
				break;
			case ALLEGRO_KEY_Z:
				angleInc(&(world->roulis), -1);
				break;
			case ALLEGRO_KEY_PAD_ASTERISK:
				vRef *= 2;
				break;
			case ALLEGRO_KEY_PAD_SLASH:
				vRef /= 2;
				break;
			case ALLEGRO_KEY_PAD_PLUS:
			{
				if(abs(vRef) > 1){
					vRef += 1;
				} else {
					vRef += 1.0/24;
				}
				break;
			}
			case ALLEGRO_KEY_PAD_MINUS:
			{
				if(abs(vRef) > 1){
					vRef -= 1;
				} else {
					vRef -= 1.0/24;
				}
				break;
			}
			case ALLEGRO_KEY_C:
			{
				PhysicObject* drone = world->getPhysicObject(physObjectIndex);
				drone->Kc1-=10;
			}
				break;
			case ALLEGRO_KEY_V:
			{
				PhysicObject* drone = world->getPhysicObject(physObjectIndex);
				drone->Kc1+=10;
			}
				break;
				
		}
	} else {
		//if(selector == -1){
			if(allegro->isKeyDown(ALLEGRO_KEY_RIGHT))
				world->offset_y += 5;
			
			if(allegro->isKeyDown(ALLEGRO_KEY_LEFT))
				world->offset_y -= 5;
		//}
		
		if(allegro->isKeyDown(ALLEGRO_KEY_UP))
			world->offset_z += 5;
		
		if(allegro->isKeyDown(ALLEGRO_KEY_DOWN))
			world->offset_z -= 5;
		
		if(allegro->isKeyDown(ALLEGRO_KEY_D))
			angleInc(&(world->lacet), 1);
		
		if(allegro->isKeyDown(ALLEGRO_KEY_Q))
			angleInc(&(world->lacet), -1);
			
		if(allegro->isKeyDown(ALLEGRO_KEY_S))
			angleInc(&(world->roulis), 1);
			
		if(allegro->isKeyDown(ALLEGRO_KEY_Z))
			angleInc(&(world->roulis), -1);
			
		if(allegro->isKeyDown(ALLEGRO_KEY_X))
			tTestProfondeur += 1;
			
		if(allegro->isKeyDown(ALLEGRO_KEY_W))
			tTestProfondeur -= 1;
	}
	world->camera = Vec((double)(WIDTH/2+world->offset_x), (double)(HEIGHT/2+world->offset_y), -fov+world->offset_z);
}


int main(int argc, char **argv)
{
	using namespace boost::gregorian;
	
	Allegro allegro_obj = Allegro();
	Allegro* allegro = &allegro_obj;
    allegro->init();
    allegro->createWindow(60, HEIGHT, WIDTH);
	
	allegro->setStickCursorToCenter(true);
	allegro->setCursorVisibility(false);
	
	World world;
	world.allegro = allegro;
	world.width_offset = (WIDTH - allegro->getDisplayWidth())/2;
	
	world.addLight(Sphere(physToRt(Vec(0, 150, 0)), 10, Color(255, 255, 170)));
	
	//world.addLight(Sphere(Vec(WIDTH/2, 0, 0), 10, Color(0, 0, 255)));
	// Moving light
	
	
	world.camera = Vec((double)(WIDTH/2), (double)(HEIGHT/2), -fov);
	
	//world.addLight(Sphere(world.camera - Vec(0, 0, -15), 5, Color(255, 255, 255)*0.2));
	
	allegro->setContext(&world);
	
	date aujourdhui = day_clock::local_day();
	
	time_t curr_time = time(NULL);
	struct tm *aTime = localtime(&curr_time);
	int year = aujourdhui.year();
	int month = aujourdhui.month();
	int day = aujourdhui.day();
	int hour = aTime->tm_hour;
	int min = aTime->tm_min;
	int sec = aTime->tm_sec;
	
	stringstream filename;
	filename << "drone-" << year << "-" << month << "-" << day << "-" << hour << "." << min << "." << sec << ".csv";
	
	PhysicObject drone(rtToPhys(Vec(0, 0, 10)), Vec(0, 0, 0), 1.5, 1, 0.01);
	drone.recordData(filename.str());
	
	physObjectIndex = world.addPhysicalObject(&drone);
	
	Sphere droneSphere(drone.getPos(), 0.5, Color(255, 0, 0), 2);
	
	droneSphere.opacity = 0.3;
	
	physObjectSphereIndex = world.addObject(droneSphere);
	
	Text3D droneTag("Drone", drone.getPos()+Vec(1, 0, 0));
	
	droneTagIndex = world.addTextTag(droneTag);
	
	//Vec a = Vec(allegro->getDisplayHeight(),0,0);
	
	//Plan p = Plan(a, (a+Vec(0, 0, 100)), (a+Vec(0, 100, 0)), Color(150, 150, 150), true);
	
	//world.addObject(p);
	
	// The Earth
	Sphere sol(physToRt(Vec(0, -6000*1000, 0)), 6000*1000, Color(58, 157, 35));
	world.addObject(sol);
	
	world.offset_y = -world.camera._y;
	world.offset_x = -10 - world.camera._x;
	world.offset_z = -80 - world.camera._z;
	
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
