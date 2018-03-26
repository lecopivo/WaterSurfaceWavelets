#include <iostream>


#include "../math/Integration.h"
#include "../math/ArrayAlgebra.h"


int main(){

  auto bar = [](double x){
	       return x;
	     };


  auto r = integrate(100,0,1,bar);

  std::cout << "Result of integration should be: 0.5" << std::endl;
  std::cout << "Result of integration is: " << r << std::endl;;


  auto foo = [](double x){
	       return std::array<double,4>{0.0,1.0,x,x*x};
	     };

  auto result = integrate(100,0,1,foo);

  std::cout << "Result of integration should be: ";
  std::cout << "0 1 0.5 0.33333" << std::endl;
  std::cout << "Result of integration is: ";
  std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << std::endl;
    
}
