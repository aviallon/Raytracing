# Raytracing
This is a Raytracing made from scratch in C++, using a C++ wrapper of Allegro I did.


# How to build :
Compile using : -std=c++11 -lpthread;allegro;allegro_image;allegro_primitives;allegro_memfile;allegro_ttf;allegro_font
-Wl,--format=binary -Wl,arial.ttf -Wl,--format=default

The last line is required to include the arial.ttf file __inside__ the executable, at link time. This avoids having to share your program and the font file as two separate files.

The optimisation level -O2 is highly recommended as it's 5 times as performant as lesser ones. Also -O3 and higher causes weird shadows glitches I didn't identify yet.

