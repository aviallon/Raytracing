#ifndef PHYS_HPP_
#define PHYS_HPP_

#include <iostream>
#include <chrono>
#include <cmath>
#include <sstream>
#include <cstdio>
#include <string>
#include <fstream>
#include <thread>
#include <algorithm>

#include "math.hpp"

using namespace std;
using namespace std::chrono;

/**
 * @brief Waits for the specified time in microsecs.
 * 
 * If threadSleep is set to true, the OS will put the thread in sleep during the waiting, so that it doesn't consume much ressources. But it is far less precise than default mode. Use it only for durations above 100ms.
 * @param microsecs Number of microseconds of sleeping.
 * @param threadSleep Uses the OS to pause the thread for the specified time. Less precise.
 */
void waitFor(int microsecs, bool threadSleep=false){
	using namespace std::chrono;
	if(threadSleep){
		this_thread::sleep_for(microseconds(microsecs));
	} else {
		auto t1 = high_resolution_clock::now();
		auto t2 = high_resolution_clock::now();
		duration<double, std::micro> dt = t2-t1;
		while(dt.count() < microsecs){
			t2 = high_resolution_clock::now();
			dt = t2-t1;
		}
	}
}

Vec rtToPhys(Vec v){
	return Vec(v._y, -v._x, v._z);
}

Vec physToRt(Vec v){
	return Vec(-v._y, v._x, v._z);
}

class PhysicObject {
public:
	Vec precAccel;
	Vec speed;
	Vec pos;
	double mass; //kg
	double charge;
	double fluidFrictionCoeff;
	
	time_point<high_resolution_clock> tSpeed;
	time_point<high_resolution_clock> tPos;
	time_point<high_resolution_clock> t0;

	double roulis = 0, lacet = 0, tanguage = 0;
	
	ofstream courbes;
	
	PhysicObject(Vec p0, Vec v0 = Vec(0,0,0), double m = 1, double q = 0, double fluidK = 0.01){
		tSpeed = high_resolution_clock::now();
		tPos = high_resolution_clock::now();
		t0 = high_resolution_clock::now();
		mass = m;
		speed = v0;
		pos = p0;
		charge = q;
		fluidFrictionCoeff = fluidK;
	}

	~PhysicObject(){
		if(courbes.is_open())
			courbes.close();
	}
	
	bool recordData(string filename){
		courbes.open(filename, ios::out);
		if(courbes.bad())
			return false;
		courbes << "t;y;v_y;a_y;rpm" << endl;
		if(courbes.good())
			return true;
	}
	
	void saveData(){
		time_point<high_resolution_clock> t = high_resolution_clock::now();
		duration<double, std::ratio<1> > dt = t-t0;
		
		courbes << dt.count() << ";" << pos._y << ";" << speed._y << ";" << precAccel._y << ";" << rpm << endl;
	}

//	Vec getEarthAxis(){
//		return Vec(sin(roulis)+cos(roulis), cos(roulis)-sin(roulis), 1); // Not really working.
//	}

	bool getSurfaceContact(Vec surfacePos, Vec pos, Vec speed, Vec surfaceNormal){
		//cout << "\tSurfContact : " << (surfaceNormal|speed) << endl;
		return ((surfaceNormal|speed) < 0);
	}

	Vec getNormalReaction(Vec surfaceNormal, Vec force){
		return surfaceNormal*(force|surfaceNormal)*(-1);
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
		double res = 0;
		if(isnan(pos._y)){
			return 0;
		}
		double delta = between(100*(alt-pos._y)/between(abs(speed._y), 1, 10), -equilibre*0.25, equilibre*0.25);//1.2*abs(speed._y)*abs(pos._y-alt)/alt;
		res = equilibre + delta;
		rpm = between(res, 0, 10000);
		return rpm;
	}
	
	
//	double delta_rpm = 0;
//	double asserv(double alt, Vec pos, Vec speed, Vec accel, double mass, double coeff = 8.3e-5, double power = 1.8512){
//		double base_rpm = pow((1000/4*mass/coeff), 1/power);
//		
//		if(isnan(pos._y)){
//			delta_rpm = 0;
//		}
//		
//		double step = 1.2*abs(speed._y)*abs(pos._y-alt)/alt;
//		if(pos._y > alt){
//			delta_rpm -= step;
//		} else{
//			delta_rpm += step;
//		}
//		if(pos._y < 1){
//			delta_rpm = 10000;
//		}
//		delta_rpm = between(delta_rpm, -base_rpm/2, (10000-base_rpm/2));
//		
//		return between(base_rpm + delta_rpm, 0, 10000);
//	}

	Vec getAccel(){
		Vec normReac(0,0,0);
		Vec force = weight(mass) + fluidFriction(fluidFrictionCoeff, speed, 2) /*+ solidFriction(0.6, speed, Vec(0, 9.81*mass, 0), Vec(1, 0, 0))*/;
		
		//asservAccel(force/mass, mass);
		//asservissementV(speed);
		rpm = asservAlt(15, pos, speed, force/mass);
		//double rp = asserv(80, pos, speed, force/mass, mass);
		force += thrust(rpm, force.normalize()*-1)*4;
		
		cout << "RPM : " << rpm << ", a_y=" << force._y/mass << ", v_y= " << speed._y << ", y= " << pos._y << endl;
		
		precAccel = force/mass;
		return force/mass;
		//return Vec(0,0,0);
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
