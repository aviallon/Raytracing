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
#include <boost/date_time/gregorian/gregorian_types.hpp>
//#include <date_time/gregorian/gregorian_types.hpp>
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

bool high_fps_mode = false;

std::vector<int> renderTimeAverage = std::vector<int>(5);

#include "allegro/allegro.h"
#include "raytracer.h"

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

//void rotatePlanet(World* world, Sphere* sp, Vec point, double incl, double angle_deg, int r){
//	
//	double angle = PI*angle_deg/180;
//	
//	double inclinaison = PI*incl/180;
//	
//	double theta = world->correction - inclinaison;
//	double phi = 0;
//	
//	//Vec axis = Vec(1000.0, 0.1, 0.1).rotate(inclinaison, Vec(0, 1, 0));
//	
//	//Vec er(cos(theta), sin(theta)*sin(phi), sin(theta)*cos(phi));
//	Vec er(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
//	
//	//std::cout << (string)er << endl;
//	
//	Vec temp = er * r;
//	
//	Vec nouv_ct = temp.rotate(angle, er ^ Vec(1, 1, 1)) + point;
//	
//	sp->ct = nouv_ct;
//}

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
	
	Sphere soleil = world_ptr->lights[0];
	
	//Sphere* second_light = &(world_ptr->lights[1]);
	
	/* Animation */
	
	const int rRef = 83.23*soleil.r / 10; // on divise par 10 les rayons réels afin de pouvoir voir toutes les planètes
	//const int rayon_mercure = 57909176; // juste pour info
	Vec axeRotation = Vec(1, 0, 0);
	
	rotateSphere(mercure, soleil.ct, axeRotation, angleMercure, rRef, 0.205630690);
	//rotateSphere2Axis(mercure, world_ptr->light.ct, 0, angleMercure, rRef); // Rotation sphérique :D
	
	//rotatePlanet(world_ptr, venus, world_ptr->light.ct, 1, angleVenus, rRef*1.868597301);
	rotateSphere(venus, soleil.ct, axeRotation, angleVenus, rRef*1.868597301, 0.006800000);
	rotateSphere(terre, soleil.ct, axeRotation, angleTerre, rRef*2.583319214, 0.016710220);
	rotateSphere(lune, terre->ct, axeRotation, angleLune, rRef*0.552537013, 0.0549);
	rotateSphere(mars, soleil.ct, axeRotation, angleMars, rRef*3.936105687, 0.093412330);
	rotateSphere(jupiter, soleil.ct, axeRotation, angleJupiter, rRef*pow(13.441946178, 0.8), 0.048392660); // Sinon elle est trop loin
	rotateSphere(saturne, soleil.ct, axeRotation, angleSaturne, rRef*pow(24.54152986, 0.8), 0.054150600);
	rotateSphere(uranus, soleil.ct, axeRotation, angleUranus, rRef*pow(49.675703933, 0.8), 0.044405586);
	rotateSphere(neptune, soleil.ct, axeRotation, angleNeptune, rRef*pow(77.767358682, 0.8), 0.008585870);
	
	//rotateSphere2Axis(test, world_ptr->light.ct, angleTest1, angleTest2, rRef*2);
	
	angleInc(&angleTest1, 2/FPS);
	angleInc(&angleTest2, vRef*8/FPS);
	
	//rotateSphere(second_light, venus->ct, axeRotation, angleTest2, rRef*0.3, 0);
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
	
	temps += vRef/FPS;
	
	
	double lacetRad = PI*world_ptr->lacet / 180 + PI/2;
	//double roulisRad = PI*world_ptr->roulis / 180;
	
	double dist = 300;
	Vec temp = Vec(0, 1, 1)*dist;
	
	Vec rotation = temp.rotate2(-lacetRad, 0, Vec(1, 0, 0));
	//Vec rotationLacet = temp.rotate(lacetRad, Vec(1, 0, 0));
	//Vec rotationRoulis = (Vec(1, 0, 1)*dist).rotate(roulisRad, Vec(0, 1, 0));
	
	for (unsigned j = 0; j < 5; j++){
		for(unsigned i = 0; i < 5; i++){
			Sphere* sp = ((Sphere*)world_ptr->getObject(9+j+i*5));
			Vec pos = Vec((2*sp->r)*j, (2*sp->r)*i, 0.0);
			Vec new_ct = pos.rotate(-lacetRad + PI/2, Vec(1, 0, 0))/*.rotate(roulisRad, Vec(0, 1, 0))*/ + rotation + soleil.ct;
			//Vec new_ct = sp->ct.rotate(angleTest1, Vec(1, 0, 0)) + (Vec(1, 0, 0)-sp->ct);
			//sp->color = Color(127*sin(angleTest1 + 0) + 128, 127*sin(angleTest1 + 2*PI/3) + 128, 127*sin(angleTest1 + 4*PI/3) + 128);
			sp->ct = new_ct;
		}
	}
	
//	Sphere* refTest = ((Sphere*)world_ptr->getObject(reflectTest));
//	refTest->ct = Vec(world_ptr->roulis, world_ptr->lacet, tTestProfondeur);
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
	//chrono::time_point<chrono::system_clock> 7mars2009;
	date t0 = date(2009, 03, 7);
	date now = (t0 + days((int)temps));
	dtstream << now.day() << "/" << now.month() << "/" << now.year();
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
	
	stringstream fps_disp;
	
	renderTimeAverage.push_back((getms() - t));
	renderTimeAverage.erase(renderTimeAverage.begin());
	
	fps_disp << fps(renderTimeAverage) << " FPS\0";
	
	allegro->draw_text(30, 10, fps_disp.str(), allegro->rgb(255, 255, 255));
	
	allegro->draw_text(30, 30, getDate(temps), allegro->rgb(255, 255, 255));
	
	stringstream corr;
	corr << "N : " << round(world_ptr->correction*100)/100;
	
	allegro->draw_text(30, 50, corr.str(), allegro->rgb(255, 255, 255));
	
	world_ptr->width_offset = (WIDTH - allegro->getDisplayWidth())/2;
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

void rotateTest(World* world, Allegro* allegro){
	allegro->getGUI()->displayMessage("Disabled", 4000);
//	Sphere* tTest = ((Sphere*)world->getObject(turnTest));
//	Sphere* refTest = ((Sphere*)world->getObject(reflectTest)); // Sphere autour de laquelle on tourne
//	//Vec rotate = refTest->getNormale(tTest->ct);
//	
//	
//	Vec new_ct = tTest->ct - refTest->ct;
//	tTest->ct = new_ct.rotate(PI, refTest->getNormale(tTest->ct) ^ Vec(1, 1, 1)) + refTest->ct;
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
			case ALLEGRO_KEY_PAD_PLUS:
				vRef *= 2;
				break;
			case ALLEGRO_KEY_PAD_MINUS:
				vRef /= 2;
				break;
			case ALLEGRO_KEY_C:
			{
				for(unsigned i = 0; i<theTest.size(); i++){
					Sphere* sp = (Sphere*)world->getObject(theTest[i]);
					sp->n += 0.1;
				}
				Sphere* sp1 = (Sphere*)world->getObject(theTest[0]);
				world->correction = sp1->n;
				
			}
				break;
			case ALLEGRO_KEY_V:
			{
				for(unsigned i = 0; i<theTest.size(); i++){
					Sphere* sp = (Sphere*)world->getObject(theTest[i]);
					sp->n -= 0.1;
				}
				Sphere* sp1 = (Sphere*)world->getObject(theTest[0]);
				world->correction = sp1->n;
			}
				break;
			case ALLEGRO_KEY_T:
				stringstream message;
				message << "Test réflection ";
				Sphere* sp = (Sphere*)world->getObject(theTest[0]);
				const bool hidden = sp->hidden;
				if(hidden){
					message << "activé";
				} else {
					message << "désactivé";
				}
				allegro->getGUI()->displayMessage(message.str(), 3500);
				for(unsigned i = 0; i<theTest.size(); i++){
					Sphere* sp = (Sphere*)world->getObject(theTest[i]);
					sp->hidden = (hidden)?false:true;
				}
				//rotateTest(world, allegro);
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
	Allegro allegro_obj = Allegro();
	Allegro* allegro = &allegro_obj;
    allegro->init();
    allegro->createWindow(60, HEIGHT, WIDTH);
	
	allegro->setStickCursorToCenter(true);
	allegro->setCursorVisibility(false);
	
	World world;
	world.allegro = allegro;
	world.width_offset = (WIDTH - allegro->getDisplayWidth())/2;
	
	world.addLight(Sphere(Vec(WIDTH/2, HEIGHT/2, 0), 60, Color(255, 255, 150)));
	//world.lights.push_back();
	
	//world.addLight(Sphere(Vec(WIDTH/2, 0, 0), 10, Color(0, 0, 255)));
	// Moving light
	
	
	//world.light = Sphere(Vec(WIDTH/2, HEIGHT/2, 0), 60, Color(255, 255, 150));
	world.camera = Vec((double)(WIDTH/2), (double)(HEIGHT/2), -fov);
	
	//world.addLight(Sphere(world.camera - Vec(0, 0, -15), 5, Color(255, 255, 255)*0.2));
	
	allegro->setContext(&world);

	const int rayon_planetes = 30*world.lights[0].r;

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
	
	
//	for (unsigned i = 0; i < 9; i++){
//		Sphere* sp = (Sphere*)world.getObject(i);
//		sp->reflectiveness = 1;
//		sp->opacity = 0.4;
//	}
//	
	// test rotation
	
	
	for (unsigned j = 0; j < 5; j++){
		for (unsigned i = 0; i<5; i++){
			Vec ct;
			//Vec ct = Vec(20.0*j, 20.0*i + 30, -world.lights[0].ct._z).rotate(180, Vec(1, 0, 0));
			Sphere sp = Sphere(ct, 10, Color(255, 255, 255));
			sp.hidden = true;
			sp.opacity = 0.1;
			sp.reflectiveness = 1;
			sp.n = world.correction;
			theTest.push_back(world.addObject(sp));
		}
	}
	
	// Test reflexion
	
//	Sphere sp = Sphere(world.light.ct + Vec(100, 0, 0), 20, Color(0, 0, 0));
//	sp.reflectiveness = 0.5;
//	sp.hidden = true; // Sphere masquée (et même désactivée)
//	reflectTest = world.addObject(sp);
//	
//	Sphere sp2 = Sphere(Vec(), 20, Color(255, 255, 255));
//	sp2.hidden = true; // Sphere masquée (et même désactivée)
//	turnTest = world.addObject(sp2);
	
	// Test overload
	
	//for (unsigned i = 0; i<= 50; i++){
	//	world.addObject(Sphere(Vec(), 10, Color(0, 0, 0)));
	//}
	
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
