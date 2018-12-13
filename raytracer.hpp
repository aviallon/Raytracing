#pragma once

#ifndef RAYTRACER_H_
#define RAYTRACER_H_ 1


#include "includes.h"
//
//#include <utility>
//#include "phys.hpp"
//#include "math.hpp"
//#include "allegro/allegro.h"

using std::pair;

template <typename T> int sgn(T val);

int getms();

const float fps(std::vector<int> renderTimeAverage);

class UVTexture{
private:
	png::image<png::rgb_pixel> texture;
	bool initialized = false;
	float _w, _h;
public:
	UVTexture();
	UVTexture(std::string filename);
	
	Color getColorUV(float u, float v);
	bool isInitialized();
};

class Plan{
public:
	Plan(Vec a, Vec b, Vec c, Color color, float rn, float opacity, bool damier);
	
	Vec getNormale(Vec pI);
	
	pair<Vec, Vec> intersect(Vec o, Vec d, Vec camera);
	
	Color getColor(Vec p);
	
	Vec a;
	Vec b;
	Vec c;
	Vec n;
	Color color;
	bool damier;
	float opacity;
	float rn;
};


class Sphere{
public:
	Sphere(Vec c, double r, Color color, double n=1.33, double opacity=1);
	
	Vec getNormale(Vec pI);
	
	pair<Vec, Vec> intersect(Vec o, Vec d, bool getAll = false);
	
	Color getColor(Vec pI);
	
	Vec ct;
	double r, n;
	double reflectiveness = 0;
	double opacity = 1;
	Color color;
	bool hidden = false;
	
	UVTexture texture;
};


enum Obj{
		SPHERE	= (1<<0),
		PLAN	= (1<<1),
		LIGHT	= (1<<2)
};

struct Indices{
	Indices(int i=0, int obj_type = Obj::SPHERE){
		this->i = i;
		this->obj_type = obj_type;
	}
	int i;
	int obj_type;
};

class Text3D;

class World{
protected:

	std::vector<Indices> indices;
	
	std::vector<Sphere> spheres;
	std::vector<Plan> plans;
	
	std::vector<Text3D> textTags;
	
	std::vector<PhysicObject*> physicObjects;
public:

	World();
	
	int addObject(Sphere sphere);
	
	int addObject(Plan plan);
	
	int addPhysicalObject(PhysicObject* obj);
	
	int addTextTag(Text3D txtTag);
	
	int addLight(Sphere light);
	
	void* getObject(int i);
	
	PhysicObject* getPhysicObject(int i);
	
	//unsigned getPhysicObjectsNumber(int )
	
	Text3D* getTextTag(int i);
	
	void drawTextTag(int i);
	
	Vec getNormale(int i, Vec pI);
	
	pair<Vec, Vec> intersect(int i, Vec o, Vec d, bool getAll = false, bool nolight = false);
	
	Color getColor(int i, Vec p);
	
	const int getType(int i);
	
	const double getReflectiveness(int i);
	
	const double getOpacity(int i);
	
	const float getRefractionIndice(int i);
	
	const unsigned int size();
	const unsigned int size_no_lights();
	
	Vec getLightCt(int i);
	
	std::vector<Sphere> lights;
	Vec camera;
	
//	int offset_x = 0;
//	int offset_y = 0;
//	int offset_z = 0;
	
	double correction = 1.5;
	
	int width_offset = 0;
	
	double tangage = 0;
	double lacet = 0;
	double roulis = 0;
	
	bool enableAntialiasing = false;
	
	bool high_fps_mode = false;
	bool disableAutoHighFPS = false;
	
	bool screenRotatedResized = false;
	
	Matrice rotation;
	
	Vec direction;
	vector<vector<Vec>> screen;
	
	Allegro* allegro;
};


class Text3D{
public:
	Text3D(std::string text, Vec pos);
	
	void setText(std::string text);
	
	std::string getText();
	
	void setPos(Vec pos);
	
	Vec getPos();
	
	void draw(Allegro* allegro, Vec camera);
	
	Vec pos;
	std::string text;
};

struct Impact{
	Impact(Vec impct, uint objId);
	
	Impact(bool noImpact);
	
	Vec impact;
	uint objectId;
	bool noImpact = false;
};

Impact getNearestImpact(Vec o, Vec d, World& world);

bool shadowRay(Vec pI, Vec L, World& world, unsigned lampe, int i);

// Renvoie Rtm puis Ttm
std::pair<double, double> fresnelCoefficient(double i, double r, float n1, float n2);

Color raytrace(World* world, Vec origine, Vec direction, int depth = 0);

#endif
