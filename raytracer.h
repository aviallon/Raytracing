#ifndef RAYTRACER_H_
#define RAYTRACER_H_

#include <utility>

using std::pair;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int getms(){
	using namespace std::chrono;
	return duration_cast<milliseconds>(
		system_clock::now().time_since_epoch()
	).count();
}

const float fps(std::vector<int> renderTimeAverage){
	float average = std::accumulate(renderTimeAverage.begin(), renderTimeAverage.end(), 0.0)/renderTimeAverage.size();
	return floor(10000/average)/10;
}

class Vec {
public:
	Vec(double x, double y, double z){
		_x = x;
		_y = y;
		_z = z;
	}
	Vec(int x, int y, int z){
		_x = double(x);
		_y = double(y);
		_z = double(z);
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

//inline int between(double x, int min, int max){
//	return floor((x>max) ? max : ((x<min) ? min : x));
//}

class Plan{
public:
	Plan(Vec a, Vec n, Color color, bool damier = false){
		this->a = a;
		this->n = n.normalize();
		this->damier = damier;
		//std::cout << n.len() << ", " << this->n.len() << std::endl;
		this->color = color;
	}
	
	Vec getNormale(Vec pI){
		return n;
	}
	
	pair<Vec, Vec> intersect(Vec o, Vec d){
		pair<Vec, Vec> intersections(Vec(true), Vec(true));
		
		Vec d_o = d-o;
		
		/** A OBSERVER ET TRAVAILLER !!! */
		d_o = Vec(-d_o._x, d_o._y, d_o._z-WIDTH/2);
		// --> Généralisation non triviale...
		// Voir sur internet Changement de base 3d -> PDF d'une université informatique intéressant.
		
		
		//std::cout << (o-a).len() << std::endl;
		double _d = -(n.dot(a));
		if (n.dot(d_o.normalize()) != 0){
			float k = (n.dot(o)+_d)/(n.dot(d_o));;
			if(k<0){ // WTF il se passe des trucs étranges ici
				return intersections;
			}
			intersections.first = d_o*k+o;
			return intersections;
		}
		return intersections;
	}
	
	Color getColor(Vec p){

		if(!damier)
			return color;
		
		const int motifSize = 100;
		Vec p_inMotif = p%motifSize;
		
		const double sgn_y = (p_inMotif._y < 0)?-1:1;
		const double sgn_z = (p_inMotif._z < 0)?-1:1;
			
		if((p_inMotif._y) <= sgn_y*motifSize/2 && (p_inMotif._z) <= sgn_z*motifSize/2) {
			return Color(50, 50, 50);
		} else if((p_inMotif._y) > sgn_y*motifSize/2 && (p_inMotif._z) > sgn_z*motifSize/2){
			return Color(50, 50, 50);
		}
		
		return color;
	}
	
	Vec a;
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
		double a, b, c;
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

class World{
public:

	World(){
		//sz = 0;
	}
	
	
	int addObject(Sphere sphere){
		spheres.push_back(sphere);
		indices.push_back(Indices(spheres.size()-1, Obj::SPHERE));
		return indices.size()-1;
	}
	
	int addObject(Plan plan){
		plans.push_back(plan);
		indices.push_back(Indices(plans.size()-1, Obj::PLAN));
		return indices.size()-1;
	}
	
	int addLight(Sphere light){
		lights.push_back(light);
		return lights.size()-1;
		//indices.push_back(Indices(lights.size()-1, Obj::LIGHT | Obj::SPHERE));
	}
	
	void* getObject(int i){
		if((unsigned)i < indices.size()){
			if(indices[i].obj_type == Obj::SPHERE){
				return &(spheres[indices[i].i]);
			} else {
				return &(plans[indices[i].i]);
			}
		} else {
			return &(lights[i - indices.size()]);
		}
	}
	
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
				return plans[indices[i].i].intersect(o, d);
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
	//Sphere light = Sphere(Vec(0, 0, 0), 5, Color(255, 255, 255));
	Vec camera;
	
	int offset_x = 0;
	int offset_y = 0;
	int offset_z = 0;
	
	double correction = 1.5;
	
	int width_offset = 0;
	
	double tangage = 0;
	double lacet = 0;
	double roulis = 0;
	
	Allegro* allegro;
	
private:

	std::vector<Indices> indices;
	
	std::vector<Sphere> spheres;
	std::vector<Plan> plans;
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

Impact getNearestImpact(Vec o, Vec d, World& world){
	Vec nearestImpact(INFINITY, INFINITY, INFINITY);
	int objId = -1;
	#pragma omp parallel for
	for(uint j=0;j<world.size();j++){
		Vec impact = world.intersect(j, o, d).first;
		Vec d_o = d - o;
		if(impact.nonVec || (impact-o).len() > d_o.len())
			continue;
			
		if((impact-o).len() < (nearestImpact-o).len()){
			nearestImpact = impact;
			objId = j;
		}
	}
	
	if(objId != -1)
		return Impact(nearestImpact, objId);
	else
		return Impact(true);
	
}

bool shadowRay(Vec pI, Vec L, World& world, unsigned lampe, int i){
	Sphere light = world.lights[lampe];
	for(unsigned obstacle = 0; obstacle < world.size(); obstacle++){
		if((int)obstacle == i || obstacle == lampe + world.size(false)){
			continue;
		}
		
		Vec point_intersection_obstacle = world.intersect(obstacle, pI, light.ct).first;
		
		if (point_intersection_obstacle.nonVec){
			continue;
		}
		
		double angle = (pI - light.ct).normalize().dot((point_intersection_obstacle - light.ct));
		
		if(angle < 0){
			continue;
		}
		
		if((pI - light.ct).len() > (point_intersection_obstacle - pI).len()){
			return true;
		}
	}
	
	return false;
}

// Renvoie Rtm puis Ttm
std::pair<double, double> fresnelCoefficient(double i, double r, float n1, float n2){
	const double n1cosi = n1*cos(i); // Permet d'éviter de recalculer plusieurs fois les mêmes cos, les calculs trigonométriques étant très couteux à l'ordinateur
	const double n2cosr = n2*cos(r);
	
	const double reflect = (n1cosi - n2cosr)/(n1cosi + n2cosr);
	
	const double transmission = 2*n1cosi/(n1cosi + n2cosr);
	
	return std::pair<double, double>(reflect, transmission);
}

Color raytrace(World* world, Vec origine, Vec direction, int depth = 0){
	double nobgreflect = 1;
	int max_depth = 4;
	
	if(depth > 0)
		nobgreflect = 1;
	
	Color bg(50, 50, 50);
	//Color pixel(0, 0, 0);
	Color pixel = bg*nobgreflect;// =  bg*nobgreflect;
	bool light_here = false;
	
	for(unsigned lampe = 0; lampe < world->lights.size(); lampe++){
	
		Sphere light = world->lights[lampe];
		
//		if(lampe == 0)
//			pixel = bg*nobgreflect;
		Color temp_pixel;
		bool ok = false;
			
		Vec pI(INFINITY, INFINITY, INFINITY);
		Vec ptemp = light.intersect(origine, direction).first;
		if(ptemp.nonVec != true){
			pI = ptemp;
			light_here = true;
		}
		bool shadow = false;
		for(unsigned int i=0;i<world->size();i++){
			double dt = 0;
			shadow = false;
			ptemp = world->intersect(i, origine, direction).first;
			
			if((ptemp.nonVec != true && (ptemp-origine).len() < (pI-origine).len())){
				ok = true;
				light_here = false;
				pI = ptemp;
				Vec L = light.ct - pI;
				
				if(!high_fps_mode){
					shadow = shadowRay(pI, L, *world, lampe, i);
				} else {
					max_depth = 1;
				}
				
				Vec N = world->getNormale(i, pI);
				dt = (L.normalize().dot(N));
				
				if(shadow){
					dt = 0.01f;
				}
				
				Color reflection(0, 0, 0);
				Color transmission(0, 0, 0);
				std::pair<double, double> coeffs;
				
				if(depth < max_depth){
				
					Vec r = (origine - pI);
					Vec B = N*r.dot(N);
					Vec A = r - B;
					
					if(world->getReflectiveness(i) > 0){
						if((origine - pI).len() > 10){
							Vec reflect = (r - A);
							
							reflection = raytrace(world, pI, reflect*1000, depth+1)*r.normalize().dot(N);
						}
					}
					
					
					if(world->getOpacity(i) < 1){
						
						const float n = world->getRefractionIndice(i);
						
						const double angle_incidence1 = acos(r.normalize().dot(N));
						
						const double r_angle = asin((1/n) * sin(  angle_incidence1  ) );
						
						
						coeffs = fresnelCoefficient(angle_incidence1, r_angle, 1, n);
						
						Vec Aprime = A.normalize()*tan(r_angle)*B.len();

						Vec refract1 = (Aprime+B)*-1;
						
						const Vec pI2 = world->intersect(i, pI, refract1*1000, true).second;
						
						if(true || !pI2.nonVec){
						
							Vec N2 = world->getNormale(i, pI2)*-1;
							
							Vec r2 = (pI - pI2);
							Vec B_2 = N2*r2.dot(N2);
							Vec A2 = r2 - B_2;
							
							const double angle_incidence2 = acos(r2.normalize().dot(N2));
							
							const double r_angle2 = asin((n) * sin(  angle_incidence2  ) );
							
							const std::pair<double, double> coeffs2 = fresnelCoefficient(angle_incidence2, r_angle2, n, 1);
							const double transmis2 = MAX(MIN(coeffs.second * coeffs2.second, 1), 0);

							
							Vec Aprime2 = A2.normalize()*tan(r_angle2)*B_2.len();
							Vec refract2 = (Aprime2+B_2)*-1;
							
							transmission = raytrace(world, pI2, refract2*1000, depth+1)*transmis2;
						}
					}
					
				}
				
				if(world->getType(i) == Obj::LIGHT){
					temp_pixel = world->getColor(i, pI);
				} else {
					temp_pixel = (world->getColor(i, pI)*world->getOpacity(i)).mix(light.color)*dt + reflection*coeffs.first/*world->getReflectiveness(i)*/ + transmission*(1-world->getOpacity(i));
				}
			}
		}
		
		if(light_here){
			pixel = light.color;
		}
		
		if(ok){
			if(lampe == 0 || light_here){
				pixel = temp_pixel;
			} else {
				pixel = pixel.blend(temp_pixel);
			}
		}
		
		// Fin boucle des lampes
	}

	return pixel;
}

#endif
