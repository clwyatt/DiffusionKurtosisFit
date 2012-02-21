// [====================] 100.0 % (xx.x sec/voxel)
// [===============>    ]  75.x % (xx.x sec/voxel)
// [==========>         ]  50.x % (xx.x sec/voxel)
// [=====>              ]  25.x % (xx.x sec/voxel)
// [>                   ]   0.x % (xx.x sec/voxel)

#include <sstream>
#include <iomanip>
#include <iostream>
#include <ctime>

#include "ProgressMeter.h"

ProgressMeter::ProgressMeter()
{
  m_iteration = 0;
  m_maxIteration = 100;
  m_numTicks = 20;
}

ProgressMeter::ProgressMeter(unsigned int maxIterations, unsigned int numTicks)
{
  m_iteration = 0;
  m_maxIteration = maxIterations;
  m_numTicks = numTicks;
}

void ProgressMeter::tick()
{
  m_time0 = clock();
}

void ProgressMeter::tock()
{
  m_time1 = clock();

  m_elapsed =  static_cast<double>(CLOCKS_PER_SEC)/static_cast<double>(m_time1 - m_time0);
  m_iteration += 1;
  m_percentage = m_iteration/m_maxIteration;
}

void ProgressMeter::render()
{
  std::ostringstream oss;
  oss << "[";
  unsigned int i;
  for(i = 0; i < (m_numTicks-1)*m_percentage; ++i)
    {
    oss << '=';
    }
  oss << ((m_percentage == 1) ? '=' : '>');
  i += 1;
  while(i < m_numTicks)
    {
    oss << ' ';
    i+= 1;
    }
  oss << "] "
      << std::setw(5)
      << std::setiosflags( std::ios::right )
      << std::setiosflags( std::ios::fixed )
      << std::setprecision(1)
      << m_percentage*100
      << " % " 
      << std::setw(5)
      << std::setiosflags( std::ios::right )
      << std::setiosflags( std::ios::fixed )
      << std::setprecision(1)
      << m_elapsed << " voxel/sec";

  std::cout << oss.str().c_str() << "\r";
  std::flush(std::cout);
}

