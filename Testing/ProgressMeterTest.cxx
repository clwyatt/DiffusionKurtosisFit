#include <unistd.h>
#include "ProgressMeter.h"

int main()
{
  ProgressMeter meter;
  for(unsigned int i = 0; i < 5; ++i)
    {
    meter.tick();
    sleep(1);
    meter.tock();
    meter.render();
    }
  return 0;
}
