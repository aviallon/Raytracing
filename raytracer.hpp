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

class Plan{
public:
	Plan(Vec a, Vec b, Vec c, Color color, bool damier);
	
	Vec getNormale(Vec pI);
	
	pair<Vec, Vec> intersect(Vec o, Vec d, Vec camera);
	
	Color getColor(Vec p);
	
	Vec a;
	Vec b;
	Vec c;
	Vec n;
	Color color;
	bool damier;
};


class Sphere{
public:
	Sphere(Vec c, double r, Color color, double n=1.33){
		this->ct = c;
		this->r = r;
		this->color = color;
		this->n = n;
	}
	
	Vec getNormale(Vec pI){
		return (pI - ct)/r;
	}
	
	pair<Vec, Vec> intersect(Vec o, Vec d, bool getAll = false){
		pair<Vec, Vec> intersections(Vec(true), Vec(true));
		if(hidden){
			return intersections;
		}
		
		Vec d_o = d-o;
		Vec o_ct = o-ct;
		double det = -1;
		double a=0;
		double b=0;
		double c=0;
		try{
			
			a = pow(d_o,2);
			a = (a==0)?1e-8:a;
			
			b = (d_o).dot(o_ct)*2;
			c = pow(o_ct,2) - r*r;
			
			det = b*b - 4*a*c;
		} catch (...){
			std::cout << "An error happened" << std::endl;
			det = -1;
		}
		

		
		if (det < 1e-6){
			return intersections;
		} else {
			double sq_det = sqrt(det);
			double k1 = MIN((-b - sq_det)/(2*a), (-b + sq_det)/(2*a));
			if(k1>=0){
				intersections.first=d_o*k1+o;
			}
			if(getAll == false)
				return intersections;
			double k2 = MAX((-b - sq_det)/(2*a), (-b + sq_det)/(2*a));
			if(k2>=0){
				intersections.second=d_o*k2+o;
			}
		}
		
		return intersections;
	}
	
	Color getColor(Vec pI){
		
		if(textured){
//			Vec M = (pI - ct).normalize();
//			double theta = acos(M.dot(Vec(1, 0, 0)));
//			double phi = acos(M.dot(Vec(0, 1, 0)));
		}
		
		return color;
	}
	
	Vec ct;
	double r, n;
	double reflectiveness = 0;
	double opacity = 1;
	Color color;
	bool hidden = false;
	bool textured = false;
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
public:

	World();
	
	int addObject(Sphere sphere);
	
	int addObject(Plan plan);
	
	int addPhysicalObject(PhysicObject* obj);
	
	int addTextTag(Text3D txtTag);
	
	int addLight(Sphere light);
	
	void* getObject(int i);
	
	PhysicObject* getPhysicObject(int i);
	
	Text3D* getTextTag(int i);
	
	void drawTextTag(int i);
	
	Vec getNormale(int i, Vec pI){
		if((unsigned)i < indices.size()){
			if(indices[i].obj_type == Obj::SPHERE){
				return spheres[indices[i].i].getNormale(pI);
			} else {
				return plans[indices[i].i].getNormale(pI);
			}
		} else {
			return lights[i - indices.size()].getNormale(pI);
		}
	}
	
	pair<Vec, Vec> intersect(int i, Vec o, Vec d, bool getAll = false, bool nolight = false){
		if((unsigned)i < indices.size()){
			if(indices[i].obj_type == Obj::SPHERE){
				return spheres[indices[i].i].intersect(o, d, getAll);
			} else {
				return plans[indices[i].i].intersect(o, d, camera);
			}
		} else if(!nolight){
			return lights[i-indices.size()].intersect(o, d, false);
		}
		
		return pair<Vec, Vec>(Vec(true), Vec(true));
	}
	
	Color getColor(int i, Vec p){
		if((unsigned)i < indices.size()){
			if(indices[i].obj_type & Obj::SPHERE){
				return spheres[indices[i].i].getColor(p);
			} else {
				return plans[indices[i].i].getColor(p);
			}
		} else {
			return lights[i - indices.size()].getColor(p);
		}
	}
	
	const int getType(int i){
		if((unsigned)i < indices.size()){
			return indices[i].obj_type;
		} else {
			return Obj::LIGHT;
		}
	}
	
	const double getReflectiveness(int i){
		if(indices[i].obj_type == Obj::SPHERE && (unsigned)i < indices.size()){
			return spheres[indices[i].i].reflectiveness;
		} else {
			return 0;
		}
	}
	
	const double getOpacity(int i){
		if(indices[i].obj_type == Obj::SPHERE && (unsigned)i < indices.size()){
			return spheres[indices[i].i].opacity;
		} else {
			return 1;
		}
	}
	
	const float getRefractionIndice(int i){
		if(indices[i].obj_type == Obj::SPHERE && (unsigned)i < indices.size()){
			return spheres[indices[i].i].n;
		} else {
			return 1;
		}
	}
	
	const unsigned int size(bool withLights = true){
		if(withLights)
			return indices.size() + lights.size();
		else
			return indices.size();
	}
	
	Vec getLightCt(int i){
		return lights[i].ct;
	}
	
	std::vector<Sphere> lights;
	Vec camera;
	
	int offset_x = 0;
	int offset_y = 0;
	int offset_z = 0;
	
	double correction = 1.5;
	
	int width_offset = 0;
	
	double tangage = 0;
	double lacet = 0;
	double roulis = 0;
	
	bool high_fps_mode = false;
	
	Matrice rotation;
	
	Vec direction;
	
	Allegro* allegro;
	
private:

	std::vector<Indices> indices;
	
	std::vector<Sphere> spheres;
	std::vector<Plan> plans;
	
	std::vector<Text3D> textTags;
	
	std::vector<PhysicObject*> physicObjects;
};


class Text3D{
public:
	Text3D(std::string text, Vec pos){
		this->pos = pos;
		this->text = text;
	}
	
	void setText(std::string text){
		this->text = text;
	}
	
	std::string getText(){
		return text;
	}
	
	void setPos(Vec pos){
		this->pos = pos;
	}
	
	Vec getPos(){
		return pos;
	}
	
	void draw(Allegro* allegro, Vec camera){
		World* world = (World*)(allegro->getContext());
		
		const double fov = (allegro->getDisplayWidth()/2)/tan(15*PI/180);
		Vec d = (pos - camera).normalize()*fov;
		d = world->rotation*d;
		//cout << (string)d << endl;
		Vec ex(1, 0, 0);
		Vec ey(0, 1, 0);
		double x = (d|ex)+allegro->getDisplayHeight()/2;
		double y = (d|ey)+allegro->getDisplayWidth()/2;
		allegro->draw_text(y, x, text, allegro->rgb(255, 255, 255), ALLEGRO_ALIGN_CENTER);
	}
	
	Vec pos;
	std::string text;
};

struct Impact{
	Impact(Vec impct, uint objId){
		impact = impct;
		objectId = objId;
	}
	
	Impact(bool noImpact){
		this->noImpact = noImpact;
	}
	
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
