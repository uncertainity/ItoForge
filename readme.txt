// Test to check if QuantLib is installed.
g++ quant_lib_test.cpp -std=c++17 -I$HOME/local/include -L$HOME/local/lib -lQuantLib -o test
g++ Pricing_Path_Surface.cpp -std=c++17 -I$HOME/local/include -L$HOME/local/lib -lQuantLib -o AllInOne

// Generate Path --> Convert to Q --> Generate Surface
g++ -std=c++17 app/main.cpp src/*.cpp -Iinclude -I$HOME/local/include -I$HOME/arka/eigen-3.4.0 -I$HOME/arka/LBFGSpp/include -L$HOME/local/lib -lQuantLib -o demo

// Run this before running the ./demo (any executable). This is because lQuantLib is not installed system wide.
export LD_LIBRARY_PATH="$HOME/local/lib:$LD_LIBRARY_PATH"

// eigen
export EIGEN_INCLUDE_PATH=/data1/sandesh/arka/eigen-3.4.0

// when you want to compile with both eigen and armadillo
g++ -std=c++17 quant_lib_test.cpp -I$HOME/local/include -I$HOME/arka/eigen-3.4.0 -L$HOME/local/lib -larmadillo -lQuantLib -o allversion

g++ -std=c++17 quant_lib_test.cpp -I$HOME/local/include -I$HOME/arka/eigen-3.4.0 -I$HOME/arka/LBFGSpp/include -L$HOME/local/lib -larmadillo -lQuantLib -o demo
