#include "WaveGrid.h"

using namespace WaterWavelets;

int main(){

  WaveGrid::Settings s;
  s.size = 50;

  s.n_x = 100;
  s.n_theta = 8;
  s.n_zeta = 1;

  s.spectrumType = WaveGrid::Settings::PiersonMoskowitz;

  WaveGrid grid(s);

  double dt = grid.cflTimeStep();

  grid.timeStep(dt);
  
  return 0;
}
