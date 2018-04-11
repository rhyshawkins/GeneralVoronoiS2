//
//    GeneralVoronoiS2 : A general Trans-dimensional Voronoic Cell program
//    for surface spherical problems.
//
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#pragma once
#ifndef generalvoronois2exception_hpp
#define generalvoronois2exception_hpp

#include <exception>

#define GENERALVORONOIS2EXCEPTION(fmt, ...) generalvoronois2exception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class generalvoronois2exception : public std::exception {
public:

  
  generalvoronois2exception(const char *srcfile,
		       const char *function,
		       int lineno,
		       const char *fmt, ...);
  ~generalvoronois2exception();
  
};

#endif // generalvoronois2exception_hpp
