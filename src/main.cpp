#include "WaveGrid.hpp"

int main(){

  WaveGrid::Settings s;
  s.xmin = 0;
  s.xmax = 100;
  s.ymin = 0;
  s.ymax = 100;

  s.n_x = 100;
  s.n_y = 100;
  s.n_theta = 8;
  s.n_k = 1;

  s.spectrumType = WaveGrid::Settings::PiersonMoskowitz;

  WaveGrid grid(s);

  double dt = grid.cflTimeStep();

  grid.timeStep(dt,0.0);
  
  return 0;
}
