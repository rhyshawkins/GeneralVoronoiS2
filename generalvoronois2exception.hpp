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
