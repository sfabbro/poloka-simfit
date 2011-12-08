#ifndef VIGNETSERVER__H
#define VIGNETSERVER__H

#include <string>
#include <dimage.h>


void reserve_vignet_in_server(const std::string& fitsfilename, const Window& window);
void get_vignet_from_server(const std::string& fitsfilename, const Window& window, Kernel& kern, double value_when_outside_fits=0);



#endif
