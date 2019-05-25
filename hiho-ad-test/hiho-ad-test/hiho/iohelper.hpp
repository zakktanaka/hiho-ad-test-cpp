#pragma once

#include <iostream>
#include <iomanip>
#include <limits>

#define HIHO_IO_MAX_LEN_DOUBLE_LSHOW std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)

#define HIHO_IO_LEFT_COUT  std::cout.setf(std::ios::left);  std::cout
#define HIHO_IO_RIGHT_COUT std::cout.setf(std::ios::right); std::cout

#define HIHO_IO_WIDTH(x) std::setw(x)

#define HIHO_IO_FUNC_WIDTH HIHO_IO_WIDTH(35)
#define HIHO_IO_TIME_WIDTH HIHO_IO_WIDTH(6)

#define HIHO_IO_VALUE(vv) #vv " : " << vv
#define HIHO_IO_TIME(vv)  #vv " : " << HIHO_IO_TIME_WIDTH << vv << " msec"
#define HIHO_IO_VALUE_TIME(vv) #vv " : " << vv.value << "(" << HIHO_IO_TIME_WIDTH << vv.duration() << " msec)"