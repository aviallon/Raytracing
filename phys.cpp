#include "includes.h"

int PhysicObject::id_iterator = 0;
vector<PhysicObject*> PhysicObject::physWorld;

void waitFor(int microsecs, bool threadSleep){
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

