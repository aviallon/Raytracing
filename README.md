# Physics simulator using Raytracing
I'm making a physics simulator based on my already-made Raytracing. The goal is to be able to simulate my drone's flight with good accuracy.


# How to build :
Compile using :
```bash
-std=c++11 -lpthread;allegro;allegro_image;allegro_primitives;allegro_memfile;allegro_ttf;allegro_font;allegro_dialog
-Wl,--format=binary -Wl,allegro/fonts/Arimo-Regular.ttf -Wl,--format=default
-fopenmp
```

The last line is required to include the font file _inside_ the executable, at link time. This avoids having to share your program and the font file as two separate files.

The optimisation level `-O3 -fexpensive-optimizations` is highly recommend.

