#ifndef KEY_H_
#define KEY_H_

//#include <string>
#include <vector>
//#include <iostream>

struct Key{
	int keycode = 0;
	bool down = false;
};

class Keys{
	
public:
	Key();
	~Key();

	std::vector<Key> keys()
};

#endif
