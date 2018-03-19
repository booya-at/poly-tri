mkdir -p build
cd build

cmake \
	  -D CMAKE_INSTALL_PREFIX=$PREFIX \
      -D CMAKE_PREFIX_PATH=$PREFIX \
      /home/kati/projects/poly-tri

make
exit 1