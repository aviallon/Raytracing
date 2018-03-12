#ifndef MATH_HPP_
#define MATH_HPP_

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

#endif
