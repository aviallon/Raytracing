#include "includes.h"

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

UVTexture::UVTexture(){
	
}

bool UVTexture::isInitialized(){
	return initialized;
}

UVTexture::UVTexture(std::string filename){
	texture = png::image<png::rgb_pixel>(filename);
	initialized = true;
	_w = texture.get_width()-1;
	_h = texture.get_height()-1;
}

Color UVTexture::getColorUV(float u, float v){
	const int w = clamp(u*_w, 0.0f, _w);
	const int h = clamp(v*_h, 0.0f, _h);
	const png::rgb_pixel pix = texture.get_pixel(w, h);
	return Color(pix.red, pix.green, pix.blue);
}

Plan::Plan(Vec a, Vec b, Vec c, Color color, bool damier = false){
	this->a = a;
	this->b = b;
	this->c = c;
	
	this->n = ((b-a)^(c-a)).normalize();
	this->damier = damier;
	//std::cout << n.len() << ", " << this->n.len() << std::endl;
	this->color = color;
}

Vec Plan::getNormale(Vec pI){
	return n;
}

pair<Vec, Vec> Plan::intersect(Vec o, Vec d, Vec camera){
	pair<Vec, Vec> intersections(Vec(true), Vec(true));
	
	Vec d_o = o-d;
	
	return intersections;
	
	if((d_o|n) < 0)
		return intersections;
	
	Vec b_a = (a-b);
	Vec c_a = (a-c);
	
	Vec h = a + b_a*(d_o|b_a) + c_a*(d_o|c_a);
	
//		Vec a_h = (a-h).normalize();
//		Vec b_h = (b-h).normalize();
//		Vec c_h = (c-h).normalize();
	
	//float angle = acos(a_h|c_h) + acos(b_h|c_h) + acos(c_h|a_h);
	
	intersections.first = h*-1;
	
	return intersections;
}

Color Plan::getColor(Vec p){

//		if(!damier)
//			return color;
//		
//		const int motifSize = 100;
//		Vec p_inMotif = p%motifSize;
//		
//		const double sgn_y = (p_inMotif._y < 0)?-1:1;
//		const double sgn_z = (p_inMotif._z < 0)?-1:1;
//			
//		if((p_inMotif._y) <= sgn_y*motifSize/2 && (p_inMotif._z) <= sgn_z*motifSize/2) {
//			return Color(50, 50, 50);
//		} else if((p_inMotif._y) > sgn_y*motifSize/2 && (p_inMotif._z) > sgn_z*motifSize/2){
//			return Color(50, 50, 50);
//		}
	
	return Color(255, 255, 255);
	//return color;
}

Sphere::Sphere(Vec c, double r, Color color, double n, double opacity){
	this->ct = c;
	this->r = r;
	this->color = color;
	this->n = n;
	this->opacity = opacity;
}

Vec Sphere::getNormale(Vec pI){
	return (pI - ct)/r;
}

pair<Vec, Vec> Sphere::intersect(Vec o, Vec d, bool getAll){
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

Color Sphere::getColor(Vec pI){
	
	if(texture.isInitialized()){
			Vec M = (pI - ct).normalize();
			//const double dx = M.dot(Vec(1, 0, 0));
			//const double dy = M.dot(Vec(0, 1, 0));
			//const double dz = M.dot(Vec(0, 0, 1));
			//const double u = 0.5 + atan2(M.dot(Vec(0, 0, 1)),  M.dot(Vec(1, 0, 0))) / (2*PI);
			//const double v = 0.5 - asin(M.dot(Vec(0, 1, 0))) / PI;
			return texture.getColorUV(0.5 + atan2(M.dot(Vec(0, 0, 1)),  M.dot(Vec(1, 0, 0))) / (2*PI), 0.5 - asin(M.dot(Vec(0, 1, 0))) / PI);
	}
	
	return color;
}

World::World(){
	
}

int World::addObject(Sphere sphere){
	spheres.push_back(sphere);
	indices.push_back(Indices(spheres.size()-1, Obj::SPHERE));
	return indices.size()-1;
}

int World::addLight(Sphere light){
	lights.push_back(light);
	return lights.size()-1;
}

int World::addObject(Plan plan){
	plans.push_back(plan);
	indices.push_back(Indices(plans.size()-1, Obj::PLAN));
	return indices.size()-1;
}

int World::addPhysicalObject(PhysicObject* obj){
	physicObjects.push_back(obj);
	return physicObjects.size()-1;
}

int World::addTextTag(Text3D txtTag){
	textTags.push_back(txtTag);
	return textTags.size()-1;
}

void* World::getObject(int i){
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

PhysicObject* World::getPhysicObject(int i){
	return physicObjects.at(i);
}

Text3D* World::getTextTag(int i){
	return &(textTags.at(i));
}

void World::drawTextTag(int i){
	if((unsigned)i < textTags.size()){
		textTags[i].draw(allegro, camera);
	}
}

Vec World::getNormale(int i, Vec pI){
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

pair<Vec, Vec> World::intersect(int i, Vec o, Vec d, bool getAll, bool nolight){
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

Color World::getColor(int i, Vec p){
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

const int World::getType(int i){
	if((unsigned)i < indices.size()){
		return indices[i].obj_type;
	} else {
		return Obj::LIGHT;
	}
}

const double World::getReflectiveness(int i){
	if(indices[i].obj_type == Obj::SPHERE && (unsigned)i < indices.size()){
		return spheres[indices[i].i].reflectiveness;
	} else {
		return 0;
	}
}

const double World::getOpacity(int i){
	if(indices[i].obj_type == Obj::SPHERE && (unsigned)i < indices.size()){
		return spheres[indices[i].i].opacity;
	} else {
		return 1;
	}
}

const float World::getRefractionIndice(int i){
	if(indices.at(i).obj_type == Obj::SPHERE){
		return spheres[indices[i].i].n;
	} else {
		return 1;
	}
}

const unsigned int World::size(bool withLights){
	if(withLights)
		return indices.size() + lights.size();
	else
		return indices.size();
}

Vec World::getLightCt(int i){
	return lights[i].ct;
}

Text3D::Text3D(std::string text, Vec pos){
	this->pos = pos;
	this->text = text;
}

void Text3D::setText(std::string text){
	this->text = text;
}

std::string Text3D::getText(){
	return text;
}

void Text3D::setPos(Vec pos){
	this->pos = pos;
}

Vec Text3D::getPos(){
	return this->pos;
}

void Text3D::draw(Allegro* allegro, Vec camera){
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

Impact::Impact(Vec impct, uint objId){
	impact = impct;
	objectId = objId;
}

Impact::Impact(bool noImpact){
	this->noImpact = noImpact;
}

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
	if(world.getType(i) == Obj::PLAN)
		return false;

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


// MOST IMPORTANT PART OF THE WHOLE PROGRAMM
Color raytrace(World* world, Vec origine, Vec direction, int depth){
	static double nobgreflect = 1;
	int max_depth = 4;
	static float (*_sin)(float) = &sin;
	
	if(world->high_fps_mode){
		max_depth = 2;
		//_sin = &st_sin; // Has little positive impact on performance.
		//_cos = &st_cos; // As below
		//_tan = &st_tan; // Has negative impact on performances
	}
	
	if(depth > 0)
		nobgreflect = 1;
	
	static Color bg(135, 206, 235);
	//static Color bg(50, 50, 50);
	//Color pixel(0, 0, 0);
	Color pixel =  bg*nobgreflect;
	bool light_here = false;
	
	for(unsigned lampe = 0; lampe < world->lights.size(); lampe++){
	
		Sphere light = world->lights[lampe];
		
//		if(lampe == 0)
//			pixel = bg*nobgreflect;
		Color temp_pixel(0, 0, 0);
		bool got_real_intersection = false;
			
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
				got_real_intersection = true;
				light_here = false;
				pI = ptemp;
				Vec L = light.ct - pI;
				
				if(!world->high_fps_mode)
					shadow = shadowRay(pI, L, *world, lampe, i);
				
				Vec N = world->getNormale(i, pI);
				dt = MIN((L.normalize().dot(N)), 1.0f);
				
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
						
						const double r_angle = asin((1/n) * _sin(  angle_incidence1  ) );
						
						
						coeffs = fresnelCoefficient(angle_incidence1, r_angle, 1, n);
						
						Vec Aprime = A.normalize()*tan(r_angle)*B.len();

						Vec refract1 = (Aprime+B)*-1;
						
						// Can be changed for stuff while INSIDE ball
						const Vec pI2 = world->intersect(i, pI, refract1*1000, true).second;
						
						if(!pI2.nonVec){
						
							Vec N2 = world->getNormale(i, pI2)*-1;
							
							Vec r2 = (pI - pI2);
							Vec B_2 = N2*r2.dot(N2);
							Vec A2 = r2 - B_2;
							
							const double angle_incidence2 = acos(r2.normalize().dot(N2));
							
							const double r_angle2 = asin((n) * _sin(  angle_incidence2  ) );
							
							const std::pair<double, double> coeffs2 = fresnelCoefficient(angle_incidence2, r_angle2, n, 1);
							const double transmis2 = MAX(MIN(coeffs.second * coeffs2.second, 1), 0);

							
							Vec Aprime2 = A2.normalize()*tan(r_angle2)*B_2.len();
							Vec refract2 = (Aprime2+B_2)*-1;
							
							transmission = raytrace(world, pI2, refract2*1000, depth+1)*transmis2;
							
							transmission = transmission.mix(world->getColor(i, pI));
						}
					}
					
				}
				
				if(world->getType(i) == Obj::LIGHT){
					temp_pixel = world->getColor(i, pI);
				} else {
					//temp_pixel = world->getColor(i, pI)*dt + reflection*coeffs.first/*world->getReflectiveness(i)*/ + transmission*(1-world->getOpacity(i));
					temp_pixel = (world->getColor(i, pI)*world->getOpacity(i)).mix(light.color)*dt + reflection*coeffs.first/*world->getReflectiveness(i)*/ + transmission*(1-world->getOpacity(i));
				}
			}
		}
		
		if(light_here){
			pixel = light.color;
		}
		
		if(got_real_intersection){
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