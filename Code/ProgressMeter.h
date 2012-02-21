// [====================] 100.0 % (xx.x sec/voxel)
// [===============>    ]  75.x % (xx.x sec/voxel)
// [==========>         ]  50.x % (xx.x sec/voxel)
// [=====>              ]  25.x % (xx.x sec/voxel)
// [>                   ]   0.x % (xx.x sec/voxel)

#ifndef _ProgressMeter_h_
#define _ProgressMeter_h_

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <ctime>

class ProgressMeter
{
public:

  ProgressMeter();
  ProgressMeter(unsigned int maxIterations, unsigned int numTicks);

  void tick();
  void tock();
  void render();

private:

  unsigned int m_numTicks;
  double m_iteration;
  double m_maxIteration;
  double m_percentage;
  clock_t m_time0, m_time1;
  double m_elapsed;
};

#endif
