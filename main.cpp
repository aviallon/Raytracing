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

const signed int WIDTH = 500; // Taille **initiale** de la fenêtre
const signed int HEIGHT = 500;

double vRef = 2; // coefficient pour accelerer le mouvement des planetes

bool high_fps_mode = false;

std::vector<int> renderTimeAverage = std::vector<int>(5);

#include "allegro/allegro.h"
#include "raytracer.h"

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
	
	Vec temp = Vec(0, 1, 1)*r; // On crée un point situé à r de l'origine
	
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
	} else if((*angle) < 0){
		*angle = 360;
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


int reflectTest;
int turnTest;

double tTestProfondeur = 0;
//std::vector<std::vector<int>> theTest;
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
	
	angleInc(&angleTest1, 2/FPS);
	angleInc(&angleTest2, 30/FPS);
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
	
	
	double lacetRad = PI*world_ptr->lacet / 180 + PI/2;
	double roulisRad = PI*world_ptr->roulis / 180;
	
	double dist = 300;
	Vec temp = Vec(0, 1, 1)*dist;
	
	Vec rotation = temp.rotate2(lacetRad, 0, Vec(1, 0, 0));
	//Vec rotationLacet = temp.rotate(lacetRad, Vec(1, 0, 0));
	//Vec rotationRoulis = (Vec(1, 0, 1)*dist).rotate(roulisRad, Vec(0, 1, 0));
	
	for (unsigned j = 0; j < 5; j++){
		for(unsigned i = 0; i<5; i++){
			Vec pos = Vec(20*(j-2), 20*(i-2), 0);
			Sphere* sp = ((Sphere*)world_ptr->getObject(9+j+i*5));
			Vec new_ct = pos.rotate(lacetRad + PI/2, Vec(1, 0, 0))/*.rotate(roulisRad, Vec(0, 1, 0))*/ + rotation + world_ptr->light.ct;
			//Vec new_ct = sp->ct.rotate(angleTest1, Vec(1, 0, 0)) + (Vec(1, 0, 0)-sp->ct);
			sp->color = Color(127*sin(angleTest1 + 0) + 128, 127*sin(angleTest1 + 2*PI/3) + 128, 127*sin(angleTest1 + 4*PI/3) + 128);
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

void redraw(Allegro* allegro, float FPS)
{
	const int t = getms();
	
	allegro->lockScreen();
	
	const Color black(0,0,0);
	
	World* world_ptr = (World*)allegro->getContext();
	
	World world2 = *world_ptr;

	
	world_ptr->camera = Vec((double)(WIDTH/2+world2.offset_x), (double)(HEIGHT/2+world2.offset_y), -fov+world2.offset_z);
	
	
	//Vec pix = screenPixRotate(x, y, world2.lacet, Vec(0, 1, 0), camera, allegro);
	// Vec(x+world2.offset_x+world2.tangage, y+world2.offset_y+world2.lacet+width_offset, 0+world2.offset_z+world2.roulis)
	
	if(!high_fps_mode && accumulate(renderTimeAverage.begin(), renderTimeAverage.end(), 0.0) / renderTimeAverage.size() > 100){ // If FPS are under 1000/_100_ = 10
		high_fps_mode = true;
		allegro->getGUI()->displayMessage("High performance mode activated", 3000);
	}
	
	
	vector<vector<Color> > screen(allegro->getDisplayHeight(), vector<Color>(allegro->getDisplayWidth()));
	
	#pragma omp parallel for collapse(2)
	for(int x=0;x<allegro->getDisplayHeight();x++){ // oui, il y a une inversion des axes...
		for(int y=0;y<allegro->getDisplayWidth();y++){
			Vec pix = Vec(x+world2.offset_x, y+world2.offset_y+world2.width_offset, 0+world2.offset_z);
			screen[x][y] = raytrace(world_ptr, world_ptr->camera, pix);
		}
	}
	
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
	
	stringstream jourstr;
	jourstr << "Jour " << (int)jour;
	
	allegro->draw_text(30, 30, jourstr.str(), allegro->rgb(255, 255, 255));
	
//	stringstream corr;
//	corr << round(world_ptr->correction/PI*100)/100 << "*PI";
//	
//	allegro->draw_text(30, 50, corr.str(), allegro->rgb(255, 255, 255));
	
	world_ptr->width_offset = (WIDTH - allegro->getDisplayWidth())/2;
}

void mouseMove(Allegro* allegro, void* context, uint16_t event, int x, int y){
	World* world = (World*)context;
	if(event == Allegro::MOUSE_WHEELED)
		world->offset_x += 3*x;
	else if(event == Allegro::MOUSE_MOVED_DELTA){
		//world->lacet += y/2;
		world->offset_z -= y/2;
		world->offset_y += x/2;
	}
}

void mouseClick(Allegro* allegro, void* context, uint16_t event, int x, int y){
	if(event == Allegro::MOUSE_L_CLICKED){
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
				world->correction += PI/12;
				break;
			case ALLEGRO_KEY_V:
				world->correction -= PI/12;
				break;
			case ALLEGRO_KEY_T:
				rotateTest(world, allegro);
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
	world.width_offset = (WIDTH - allegro->getDisplayWidth())/2;
	world.light = Sphere(Vec(WIDTH/2, HEIGHT/2, 0), 60, Color(255, 255, 150));
	world.camera = Vec((double)(WIDTH/2), (double)(HEIGHT/2), -fov);
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
	
	
	// test rotation
	
	
	for (unsigned j = 0; j < 5; j++){
		for (unsigned i = 0; i<5; i++){
			Vec ct = Vec(20.0*j, 20.0*i + 30, -world.light.ct._z).rotate(180, Vec(1, 0, 0));
			Sphere sp = Sphere(ct, 10, Color(0, 0, 0));
			//sp.hidden = true;
			sp.reflectiveness = 0.4;
			world.addObject(sp);
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
	
	//Vec test = (Vec(1, 0, 0) ^ Vec(0, 1, 1)) ^ (Vec(1, 0, 0) ^ Vec(0, 1, 1));
	//cout << (string)test << endl;
	
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
