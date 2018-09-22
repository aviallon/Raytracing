#pragma once

#ifndef MATH_HPP_
#define MATH_HPP_ 1


#include "includes.h"
//#include <vector>
//#include <iostream>
//#include <sstream>
//#include <armadillo>

#define START_BITMASK_SWITCH(x) \
    for (uint64_t bit = 1; x >= bit; bit *= 2) if (x & bit) switch (bit)

class Vec {
public:
	Vec(double x, double y, double z){
		_x = x;
		_y = y;
		_z = z;
	}
	
	Vec(){
		_x = 0;
		_y = 0;
		_z = 0;
	}
	
	Vec(bool pasVec){
		Vec(0,0,0);
		nonVec = pasVec;
	}
	
	Vec operator+(const Vec& v2){
		return Vec(this->_x+v2._x, this->_y+v2._y, this->_z+v2._z);
	}
	
	Vec& operator+=(const Vec& v){
		this->_x += v._x;
		this->_y += v._y;
		this->_z += v._z;
		return *this;
	}
	
	Vec operator-(const Vec& v2){
		return Vec(this->_x-v2._x, this->_y-v2._y, this->_z-v2._z);
	}
	
	Vec operator*(const float k){
		return Vec(this->_x*k, this->_y*k, this->_z*k);
	}
	
	Vec operator/(const float k){
		return Vec(this->_x/k, this->_y/k, this->_z/k);
	}
	
	Vec operator%(const int k){
		return Vec((int)(this->_x)%k, (int)(this->_y)%k, (int)(this->_z)%k);
	}
	
	double operator|(const Vec& v){
		return _x*v._x + _y*v._y + _z*v._z;
	}
	
	Vec operator^(const Vec& v2){
		return Vec(this->_y*v2._z - this->_z*v2._y, this->_z*v2._x - this->_z*v2._z, this->_x*v2._y - this->_y*v2._z);
	}
	
	bool operator==(const Vec& v2){
		if(this->_x==v2._x && this->_y==v2._y && this->_z==v2._z){
			return true;
		} else{
			return false;
		}
	}
	
	bool operator!=(const Vec& v2){
		return not((*this)==v2);
	}
	
	
	double dot(const Vec& v2){
		return (this->_x*v2._x + this->_y*v2._y + this->_z*v2._z);
	}
	
	double len(){
		return sqrt(this->_x*this->_x+this->_y*this->_y+this->_z*this->_z);
	}
	
	/* Angle is in radians */
	Vec rotate(double angle, Vec axis){
		Vec moi = *this;
		
		axis = axis.normalize();
		
		Vec v = moi*cos(angle) + (axis ^ moi)*sin(angle) + axis * (1 - cos(angle))*(moi.dot(axis));
		
		return v;
	}
	
	Vec rotate2(double theta, double phi, Vec axis){
		Vec moi = *this;
		
		axis = axis.normalize();
		
		Vec axis2 = axis ^ Vec(1, 1, 1);
		
		Vec v1 = moi*cos(theta) + (axis ^ moi)*sin(theta) + axis * (1 - cos(theta))*(moi.dot(axis));
		
		Vec v2 = v1*cos(phi) + (axis2 ^ v1)*sin(phi) + axis2 * (1 - cos(phi))*(v1.dot(axis2));
		
		return v2;
	}
	
	Vec normalize(){
		double l = (this->len());
		//float l = sqrt((*this).dot(*this));
		return Vec(this->_x/l, this->_y/l, this->_z/l);
	}
	
	operator std::string() {
		std::string repr;
		repr = repr + typeid(this).name() + "(" + std::to_string(_x) + "," + std::to_string(_y) + "," + std::to_string(_z) + ")";
		return repr;
	}
	
	double _x=0, _y=0, _z=0;
	bool nonVec = false;
	
};

void pVec(Vec v);

double pow(const Vec& v, int i);

using namespace arma;
class Matrice{
	mat matrix;
	
	void setColumn(double j, double x, double y, double z){
		matrix(0, j) = x;
		matrix(1, j) = y;
		matrix(2, j) = z;
	}
	
public:
	const static uint16_t X = 0;
	const static uint16_t Y = 1;
	const static uint16_t Z = 2;
	
	Matrice(unsigned n = 3){
		matrix = mat(n, n, fill::eye);
	}
	
	Matrice(mat mtx){
		matrix = mtx;
	}
	
	Matrice(Vec v){
		matrix = vec(3);
		matrix(0) = v._x;
		matrix(1) = v._y;
		matrix(2) = v._z;
	}
	
	Matrice(uint16_t axis, double angle){
		matrix = mat(3, 3, fill::zeros);
		switch(axis){
			case 0:
				setColumn(0, 1, 0, 0);
				setColumn(1, 0, cos(angle), sin(angle));
				setColumn(2, 0, -sin(angle), cos(angle));
				break;
			case 1:
				setColumn(0., cos(angle), 0, -sin(angle));
				setColumn(1, 0, 1, 0);
				setColumn(2, sin(angle), 0, cos(angle));
				break;
			case 2:
				setColumn(0, cos(angle), sin(angle), 0);
				setColumn(1, -sin(angle), cos(angle), 0);
				setColumn(2, 0, 0, 1);
		}
	}
	
	Vec toVec(){
		if(matrix.is_vec())
			return Vec(matrix(0), matrix(1), matrix(2));
		else
			return Vec(0, 0, 0);
	}
	
	Matrice operator*(Matrice& m2){
		return Matrice(matrix*m2.matrix);
	}
	
	Matrice operator*(Vec v){
		return Matrice(matrix*Matrice(v).matrix);
	}
	
	Matrice operator+(Matrice m2){
		return Matrice(matrix+m2.matrix);
	}
	
	Matrice operator*(double k){
		return Matrice(matrix*k);
	}
	
	operator Vec(){
		return this->toVec();
	}
};

inline int clamp(double x, double min, double max);

// Pre-computed sin values - 8 bit resolution
float st_sin(float x);

float st_cos(float x);

// Pre-computed tan values - 9 bit resolution (512)
float st_tan(float x);

#endif
