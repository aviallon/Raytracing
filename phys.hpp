#pragma once 

#ifndef PHYS_HPP_
#define PHYS_HPP_ 1

#include "includes.h"
//#include <iostream>
//#include <chrono>
//#include <cmath>
//#include <sstream>
//#include <cstdio>
//#include <string>
//#include <fstream>
//#include <thread>
//#include <algorithm>
//
//#include "math.hpp"
//#include "allegro/allegro.h"

using namespace std;
using namespace std::chrono;

/**
 * @brief Waits for the specified time in microsecs.
 * 
 * If threadSleep is set to true, the OS will put the thread in sleep during the waiting, so that it doesn't consume much ressources. But it is far less precise than default mode. Use it only for durations above 100ms.
 * @param microsecs Number of microseconds of sleeping.
 * @param threadSleep Uses the OS to pause the thread for the specified time. Less precise.
 */
void waitFor(int microsecs, bool threadSleep=false);

/**
 * @brief Converts a Raytracing vector to a Physics vector
 * @param v Vector to convert
 * @return Returns converted vector
 */
Vec rtToPhys(Vec v);

/**
 * @brief Converts a Physics vector to a Raytracing vector
 * @param v Vector to convert
 * @return Returns converted vector
 */
Vec physToRt(Vec v);

enum Forces{
	None = 0,
	Gravity = 1 << 0,
	FluidFriction = 1 << 1,
	SolidFriction = 1 << 2,
	Electromagnetism = 1 << 3,
	DroneThrust = 1 << 4,
	NormalReaction = 1 << 5
};

inline Forces operator|(Forces a, Forces b)
{return static_cast<Forces>(static_cast<int>(a) | static_cast<int>(b));}

/**
 * @class PhysicObject
 * @author Antoine Viallon
 * @date 17/03/18
 * @file phys.hpp
 * @brief Physics object, obeying laws of physics : they are submitted to weight, friction, and other things. Currently not really modular.
 */
class PhysicObject {
protected:
	static vector<PhysicObject*> physWorld;
	static int id_iterator;
public:
	int id;
	Vec precAccel;
	Vec speed;
	Vec pos;
	double mass; //kg
	double charge;
	double fluidFrictionCoeff;
	
	double radius;
	
	time_point<high_resolution_clock> tSpeed;
	time_point<high_resolution_clock> tPos;
	time_point<high_resolution_clock> tSave;
	time_point<high_resolution_clock> t0;

	double roulis = 0, lacet = 0, tanguage = 0;
	
	Vec positionConsigne;
	
	double Kc1 = 600;
	double Kc2 = 1;
	double Kcapm = 1;
	double KcapM = 10;
	
	bool referential = false;
	
	Forces _forces;
	
	ofstream courbes;
	
	PhysicObject(Vec p0, Vec v0 = Vec(0,0,0), Forces forces = None, double m = 1, double q = 0, double fluidK = 0.01, double radius = 1, bool isreferential = false){
		tSpeed = high_resolution_clock::now();
		tPos = high_resolution_clock::now();
		tSave = high_resolution_clock::now();
		t0 = high_resolution_clock::now();
		mass = m;
		speed = v0;
		pos = p0;
		positionConsigne = Vec(0, 10, 0);
		charge = q;
		fluidFrictionCoeff = fluidK;
		referential = isreferential;
		_forces = forces;
		
		physWorld.push_back(this);
		this->id = id_iterator;
		id_iterator++;
	}

	~PhysicObject(){
		if(courbes.is_open())
			courbes.close();
			
		for(unsigned i=0; i<physWorld.size(); i++){
			if(physWorld.at(i) == this){
				physWorld.erase(physWorld.begin() + i);
				break;
			}
		}
	}
	
	bool recordData(string filename){
		courbes.open(filename, ios::out);
		if(courbes.bad())
			return false;
		courbes << "t;y;v_y;a_y;rpm" << endl;
		if(courbes.good())
			return true;
			
		return false;
	}
	
	void saveData(){
		if(!courbes.is_open())
			return;
		time_point<high_resolution_clock> t = high_resolution_clock::now();
		duration<double, std::ratio<1> > dtSave = t-tSave;
		tSave = high_resolution_clock::now();
		
		duration<double, std::ratio<1> > dt = t-t0;
		
		if(dtSave.count() > 0.05){
			courbes << dt.count() << ";" << pos._y << ";" << speed._y << ";" << precAccel._y << ";" << rpm << endl;
		}
		//cout << "Save !" << endl;
	}
	
	void setConsignePos(Vec posConsigne){
		posConsigne = rtToPhys(posConsigne);
		if(posConsigne._y > 1)
			positionConsigne = positionConsigne;
	}
	
	Vec getConsignePos(){
		return physToRt(positionConsigne);
	}

//	Vec getEarthAxis(){
//		return Vec(sin(roulis)+cos(roulis), cos(roulis)-sin(roulis), 1); // Not really working.
//	}

	bool getSurfaceContact(Vec surfacePos, Vec pos, Vec speed, Vec surfaceNormal){
		//cout << "\tSurfContact : " << (surfaceNormal|speed) << endl;
		return ((surfaceNormal|speed) < 0);
	}

	Vec collisions(Vec force){
		Vec norm(0, 0, 0);
		for(PhysicObject* other1: physWorld){
			//cout << other1->id << endl;
			if(other1 != this){
				Vec N = (this->pos - other1->pos);
				if(N.len() < this->radius + other1->radius){
					norm += N.normalize()*(force.len());
				}
			}
		}
		return norm;
	}

	Vec weight(double mass){
		return Vec(0.0, -9.81*mass, 0.0);
	}

	Vec fluidFriction(double coeff, Vec speed, int degree=1){
		return speed*pow(speed.len(), degree-1)*(-coeff);
	}

	Vec solidFriction(double coeff, Vec speed, Vec normReaction, Vec surface){
		Vec u = (surface*(speed|surface)).normalize();
		if((speed|surface) < 10e-1)
			u = u*0.01;
		return u*coeff*normReaction.len()*-1;
	}
	
	Vec thrust(double rpm, Vec direction = Vec(0, 1, 0), double coeff = 8.3e-5, double power = 1.8512){
		return direction*coeff*pow(rpm, power)*0.001*9.81; // The last coeffs are here to convert gF (gramm Force) to Newtons.
	}
	
	double rpm = 0;
	double asservissementV(Vec speed){
		if(isnan(speed._y)){
			return 0;
		}
		double step = 2*abs(speed._y);
		if(speed._y > 0 && rpm >= step){
			rpm -= step;
		} else if (rpm <= 10000-step){
			rpm += step;
		}
		return between(rpm, 0, 10000);
	}

	double bestEquilibre(Vec accel, double mass, double coeff = 8.3e-5, double power = 1.8512){
		return pow((1000/4)*(mass/coeff), 1/power);
	}

	double asservAlt(double alt, Vec pos, Vec speed, Vec accel, double coeff = 8.3e-5, double power = 1.8512){
		double equilibre = bestEquilibre(accel, mass, coeff, power);
		cout << "\t eq : " << equilibre << endl;
		//cout << "\t Kc1 : " << Kc1 << endl;
		double res = 0;
		if(isnan(pos._y)){
			return 0;
		}
		double vcoeff = 1;
		if(speed._y > -1){
			vcoeff = 0.5;
		}
		double delta = between(vcoeff*Kc1*(alt-pos._y), -equilibre*0.1, equilibre*0.4);
		//double delta = between(Kc1*(alt-pos._y)/between(Kc2*abs(speed._y), Kcapm, KcapM), -equilibre*0.8, equilibre*0.8);//1.2*abs(speed._y)*abs(pos._y-alt)/alt;
		res = equilibre + delta;
		rpm = between(res, 0, 10000);
		return rpm;
	}
	

	Vec getAccel(){
		Vec normReac(0,0,0);
		
		Vec force(0,0,0);
		if(_forces & None){
			
		}
		if(_forces & Gravity){
			force += weight(mass);
		}
		if(_forces & FluidFriction){
			force += fluidFriction(fluidFrictionCoeff, speed, 2); // The 2 is the magnitude
		}
		if(_forces & SolidFriction){
			force += solidFriction(0.6, speed, Vec(0, 9.81*mass, 0), Vec(1, 0, 0));
		}
		if(_forces & Electromagnetism){
			
		}
		if(_forces & DroneThrust){
			rpm = asservAlt(positionConsigne._y, pos, speed, force/mass);
			force += thrust(rpm, force.normalize()*-1)*4;
		}
		if(_forces & NormalReaction){
			force += collisions(force);
		}
		
		//cout << "RPM : " << rpm << ", a_y=" << force._y/mass << ", v_y= " << speed._y << ", y= " << pos._y << endl;
		
		precAccel = force/mass;
		return force/mass;
	}

	void updateSpeed(){
		time_point<high_resolution_clock> t2 = high_resolution_clock::now();
		duration<double, std::milli> dt = t2-tSpeed;
		tSpeed = t2;
		
		double dtS = (dt.count()/1000);
		
		speed += (precAccel + getAccel())/2 * dtS;
	}

	Vec getSpeed(){
		updateSpeed();
		
		return speed;
	}

	void updatePos(){
		time_point<high_resolution_clock> t2 = high_resolution_clock::now();
		duration<double, std::milli> dt = t2-tPos;
		tPos = t2;
		
		pos += getSpeed() * (dt.count()/1000);
	}

	Vec getPos(){
		updatePos();
		
		if(courbes.is_open())
			saveData();
		
		Vec raytracingCompatibility(-pos._y, pos._x, pos._z);
		
		return raytracingCompatibility;
	}
};

// 2638 lines of code at 21:21 03/13/2018 :D

#endif
