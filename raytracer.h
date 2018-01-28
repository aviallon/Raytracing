#ifndef RAYTRACER_H_
#define RAYTRACER_H_

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
		Vec(0,0,0);
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
	
	double _x, _y, _z;
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

inline int between(double x, int min, int max){
	return floor((x>max) ? max : ((x<min) ? min : x));
}

class Color {
public:
	Color(int r, int g, int b){
		_r = between(r, 0, 255);
		_g = between(g, 0, 255);
		_b = between(b, 0, 255);
	}
	
	Color(){
		Color(0, 0, 0);
	}
	
	Color(bool notColor){
		Color(0, 0, 0);
		this->notColor = notColor;
	}
	
	Color operator*(double k){
		return Color(between(this->_r*k, 0, 255), between(this->_g*k, 0, 255), between(this->_b*k, 0, 255));
	}
	
	Color operator+(const Color& c){
		return Color(between(this->_r+c._r, 0, 255), between(this->_g+c._g, 0, 255), between(this->_b+c._b, 0, 255));
	}
	
	Color mix(const Color& c){
		return Color(between(this->_r, 0, c._r), between(this->_g, 0, c._g), between(this->_b, 0, c._b));
	}
	
	int _r, _g, _b;
	bool notColor = false;
};

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
	
	Vec intersect(Vec o, Vec d){
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
				return Vec(true);
			}
			return d_o*k+o;
		}
		return Vec(true);
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
	
	Vec intersect(Vec o, Vec d){
		if(hidden){
			return Vec(true);
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
			return Vec(true);
		} else {
			double sq_det = sqrt(det);
			double k = MIN((-b - sq_det)/(2*a), (-b + sq_det)/(2*a));
			if(k<0){
				return Vec(true);
			}
			return d_o*k+o;
		}
		
		return Vec(true);
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
	Color color;
	bool hidden = false;
	bool textured = false;
};


enum Obj{
		SPHERE	= (1<<0),
		PLAN	= (1<<1)
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
	
	void* getObject(int i){
		if(indices[i].obj_type == Obj::SPHERE){
			return &(spheres[indices[i].i]);
		} else {
			return &(plans[indices[i].i]);
		}
	}
	
	Vec getNormale(int i, Vec pI){
		if(indices[i].obj_type == Obj::SPHERE){
			return spheres[indices[i].i].getNormale(pI);
		} else {
			return plans[indices[i].i].getNormale(pI);
		}
	}
	
	Vec intersect(int i, Vec o, Vec d){
		if(indices[i].obj_type == Obj::SPHERE){
			return spheres[indices[i].i].intersect(o, d);
		} else {
			return plans[indices[i].i].intersect(o, d);
		}
	}
	
	Color getColor(int i, Vec p){
		if(indices[i].obj_type == Obj::SPHERE){
			return spheres[indices[i].i].getColor(p);
		} else {
			return plans[indices[i].i].getColor(p);
		}
	}
	
	const int getType(int i){
		return indices[i].obj_type;
	}
	
	const double getReflectiveness(int i){
		if(indices[i].obj_type == Obj::SPHERE){
			return spheres[indices[i].i].reflectiveness;
		} else {
			return 0;
		}
	}
	
	const unsigned int size(){
		return indices.size();
	}
	
	Sphere light = Sphere(Vec(0, 0, 0), 5, Color(255, 255, 255));
	Vec camera;
	
	int offset_x = 0;
	int offset_y = 0;
	int offset_z = 0;
	
	double correction = PI/2;
	
	int width_offset = 0;
	
	double tangage = 0;
	double lacet = 0;
	double roulis = 0;
	
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
		Vec impact = world.intersect(j, o, d);
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

bool shadowRay(Vec pI, Vec L, World& world, Sphere& light, int i){
	for(unsigned obstacle = 0; obstacle < world.size(); obstacle++){
		if((int)obstacle == i){
			continue;
		}
		
		Vec point_intersection_obstacle = world.intersect(obstacle, pI, light.ct);
		
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

Color raytrace(World* world, Vec origine, Vec direction, int depth = 0){
	double nobgreflect = 1;
	int max_depth = 3;
	
	if(depth > 0)
		nobgreflect = 0;
	
	Color bg(50, 50, 50);
	
	Sphere light = world->light;
	
	Color pixel = bg*nobgreflect;
	Vec pI(INFINITY, INFINITY, INFINITY);
	Vec ptemp = light.intersect(origine, direction);
	if(ptemp.nonVec != true){
		pI = ptemp;
		pixel = light.color;
	}
	bool shadow = false;
	for(unsigned int i=0;i<world->size();i++){
		double dt = 0;
		shadow = false;
		ptemp = world->intersect(i, origine, direction);
		
		if((ptemp.nonVec != true && (ptemp-origine).len() < (pI-origine).len())){
			pI = ptemp;
			Vec L = light.ct - pI;
			
			if(!high_fps_mode){
				shadow = shadowRay(pI, L, *world, light, i);
				max_depth = 1;
			}
			
			Vec N = world->getNormale(i, pI);
			dt = (L.normalize().dot(N));
			
			if(shadow){
				dt = 0.01f;
			}
			
			Color reflection(0, 0, 0);
			
			if(depth < max_depth && world->getReflectiveness(i) > 0){
				if((origine - pI).len() > 10){
					Vec r = (origine - pI);
					Vec B = N*r.dot(N);
					Vec A = r - B;
					Vec reflect = (r - A) *  50;
					reflection = raytrace(world, pI, reflect, depth+1)*r.normalize().dot(N);
				}
			}
			
			pixel = world->getColor(i, pI).mix(light.color)*dt + reflection*world->getReflectiveness(i);
		}
	}

	return pixel;
}


//void displayTextMessage(Allegro* allegro, std::string message){
//	std::cout << "[DISPLAYED] " << message << std::endl;
//	
//	allegro->getGUI()->displayMessage(message, 5000);
//}

#endif
