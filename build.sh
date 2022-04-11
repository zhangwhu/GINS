echo "Configuring and building RTKLIB ..."
mkdir build
cd build 
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j4
