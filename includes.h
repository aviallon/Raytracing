#pragma once

#include <iostream>
#include <allegro5/allegro.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_primitives.h>
#include <cmath>
#include <ctime>
#include <numeric>
#include <fstream>
#include <cstdio>

#include <string.h>
#include <sstream>
#include <boost/lexical_cast.hpp>

#include <typeinfo>
#include <vector>
#include <chrono>
#include <thread>
#include <utility>
#include <algorithm>
#include <chrono>
#include <mutex>
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <future>
#include <armadillo>
#include <png++/png.hpp>

#define MAX(a, b) (((a > b))?(a):(b))
#define MIN(a, b) (((a < b))?(a):(b))
#define PI 3.14159265359

const signed int WIDTH = 500; // Taille **initiale** de la fenêtre
const signed int HEIGHT = 500;

#include "allegro/allegro.h"

#include "math.hpp"
#include "phys.hpp"
#include "raytracer.hpp"