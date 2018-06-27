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

Color raytrace(World* world, Vec origine, Vec direction, int depth){
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
				
				if(!world->high_fps_mode){
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
							
							transmission = transmission.mix(world->getColor(i, pI));
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