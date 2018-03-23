mkdir -p build
cd build

cmake -D CMAKE_BUILD_TYPE=Release \
	  -D CMAKE_INSTALL_PREFIX=$PREFIX \
      -D CMAKE_PREFIX_PATH=$PREFIX \
      ~/projects/poly-tri

make install