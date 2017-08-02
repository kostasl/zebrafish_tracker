/*
* This file is part of fastms.
*
* Copyright 2014 Evgeny Strekalovskiy <evgeny dot strekalovskiy at in dot tum dot de> (Technical University of Munich)
*
* fastms is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* fastms is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with fastms. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef UTIL_CUDA_CHECK_H
#define UTIL_CUDA_CHECK_H

#if !defined(DISABLE_CUDA) && defined(__CUDACC__)

#include <string>



#define CUDA_CHECK cuda_check(__FILE__, __LINE__)
void cuda_check(std::string file, int line);



#endif // !defined(DISABLE_CUDA) && defined(__CUDACC__)

#endif // UTIL_CUDA_CHECK_H
