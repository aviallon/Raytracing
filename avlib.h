#pragma once
#ifndef AVLIB_H_
#define AVLIB_H_

#ifndef STRING_H_
#define STRING_H_
#include <string.h>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#endif

std::string toString(size_t entier){
	return boost::lexical_cast<std::string>(entier);
}

inline void d(std::string s){
	#ifndef DEBUG_DISABLED
	std::cout << "DEBUG : " << s << std::endl;
	#endif
}

inline void d(int s){
	#ifndef DEBUG_DISABLED
	std::cout << "DEBUG : " << s << std::endl;
	#endif
}

inline void d(float s){
	#ifndef DEBUG_DISABLED
	std::cout << "DEBUG : " << s << std::endl;
	#endif
}

inline void d(char* s){
	#ifndef DEBUG_DISABLED
	std::cout << "DEBUG : " << s << std::endl;
	#endif
}

#endif
