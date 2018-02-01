# Raytracing
This is a Raytracing made from scratch in C++, using a C++ wrapper of Allegro I did.


# How to build :
Compile using :
```bash
-std=c++11 -lpthread;allegro;allegro_image;allegro_primitives;allegro_memfile;allegro_ttf;allegro_font;allegro_dialog -Wl,--format=binary -Wl,allegro/arial.ttf -Wl,--format=default
```

The last line is required to include the arial.ttf file _inside_ the executable, at link time. This avoids having to share your program and the font file as two separate files.

The optimisation level `-O3` is highly recommended as it's 5 times as performant as lesser ones.

