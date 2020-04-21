
CLANG_FLAGS = -std=c++17 -Wall -g -framework OpenCL
PNG_FLAGS = `libpng-config --cflags --ldflags`
GTEST_FLAGS = -lpthread -lgtest_main -lgtest  -lpthread

default: test

all: fast slow

headers: gfxalgebra.hpp gfximage.hpp gfxnumeric.hpp gfxpng.hpp gfxraytrace.hpp rayson.hpp CLHelpers.hpp CLStructs.hpp

libraries: /usr/local/Cellar/png++/0.2.10/include/png++/png.hpp

mrraytracer: headers libraries mrraytracer.cpp
	clang++ ${CLANG_FLAGS} ${PNG_FLAGS} mrraytracer.cpp -o mrraytracer

%.png : %.json mrraytracer
	./mrraytracer $< $@

fast: scene_gtri_ortho_flat.png \
	    scene_gtri_persp_flat.png \
			scene_gtri_ortho_phong.png \
			scene_gtri_persp_phong.png \
			scene_2spheres_ortho_flat.png \
			scene_2spheres_persp_flat.png \
			scene_2spheres_ortho_phong.png \
			scene_2spheres_persp_phong.png

slow: teatime.png

test: scene_2spheres_ortho_flat.png

/usr/local/Cellar/png++/0.2.10/include/png++/png.hpp:
	@echo -e "libpng++ library not installed\n"
	@echo -e "Installing libpng++-dev. Please provide the password when asked\n"
	@sudo brew install png++

clean:
		rm -f mrraytracer *png
