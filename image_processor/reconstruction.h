#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <iostream>
#include <string>

#include "../config.h"

bool App2OnImagefile(const std::string& imgfile, 
                   const std::string& swcfile, 
                   const std::string& markerfile,
                   const std::string vaa3dpath = __V3DSO__,
                   const std::string app2path = __APP2SO__
                   );

bool GsdtOnImagefile(const std::string& imgfile, 
                   const std::string& newfile,
                   const std::string& vaa3dpath = __V3DSO__,
                   const std::string& gsdtpath = __GSDTSO__);

#endif  // RECONSTRUCTION_H
