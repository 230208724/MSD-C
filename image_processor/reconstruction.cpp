#include <cstdlib>
// <cstdlib> 是 C++ 标准库中的一个头文件，提供了与C语言标准库（stdlib.h）相关的函数和常量。
//它包含了一些常用的函数，用于处理与程序执行环境、内存分配、类型转换等相关的操作。
// 以下是一些常见的在 <cstdlib> 头文件中定义的函数和常量：
// std::abort：异常终止程序的函数。
// std::atof、std::atoi、std::atol、std::atoll：将字符串转换为浮点数或整数。
// std::malloc、std::calloc、std::realloc、std::free：用于动态内存分配和释放。
// std::rand、std::srand：生成伪随机数的函数。
// std::system：执行系统命令的函数。
// std::exit、std::_Exit、std::quick_exit、std::atexit：控制程序的终止和退出。
// 此外，<cstdlib> 头文件还定义了一些常量，如 NULL、EXIT_SUCCESS 和 EXIT_FAILURE。
#include <iostream> 
#include <string>
#include <vector>
#include <array>
#include <fstream>

#include "../config.h"

#include "reconstruction.h"

bool App2OnImagefile(const std::string& imgfile, 
                   const std::string& swcfile, 
                   const std::string& markerfile,
                   const std::string vaa3dpath,
                   const std::string app2path) {
    // 执行系统命令
    std::string run = vaa3dpath+" -x "+app2path+" -f app2 -i "+imgfile+" -o "+swcfile+" -p "+markerfile+" 0 AUTO 0 0 0 1 5 0 0 0";
    int result = std::system(run.c_str());
    //std::exec函数用于执行系统命令ls -l。
    //与std::system函数不同的是，std::exec函数只在命令执行失败时返回，否则它会直接替换当前进程，
    //因此如果程序执行到std::exec函数调用后的代码，说明命令执行失败。

    // 检查执行结果
    if (result == 0) {
        return true;
    } else {
        return false;
    }
}

bool GsdtOnImagefile(const std::string& imgfile, 
                   const std::string& newfile,
                   const std::string& vaa3dpath,
                   const std::string& gsdtpath) {
    int bg = 0;
    int cnn = 1;
    int chan = 0;
    double zthick = 1.0;
    int u8 = 1;

    std::string gsdt_str = vaa3dpath + " -x " + gsdtpath + " -f gsdt -i " + imgfile + " -o " + newfile +
                            " -p " + std::to_string(bg) + " " + std::to_string(cnn) + " " + std::to_string(chan) +
                            " " + std::to_string(zthick) + " " + std::to_string(u8);

    int result = std::system(gsdt_str.c_str());
    if (result == 0) {
        return true;
    } else {
        return false;
    }
}





