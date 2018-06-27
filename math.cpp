#include "includes.h"

void pVec(Vec v){
	std::cout << (std::string)v << std::endl;
}

double pow(const Vec& v, int i){
	Vec v2 = v;
	if(i==2){
		return v2.dot(v2);
	} else {
		return 0;
	}
}

inline int clamp(double x, double min, double max){
	return floor((x>max) ? max : ((x<min) ? min : x));
}