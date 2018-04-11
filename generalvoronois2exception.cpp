
#include <stdio.h>
#include <stdarg.h>

#include "generalvoronois2exception.hpp"

generalvoronois2exception::generalvoronois2exception(const char *srcfile,
						     const char *function,
						     int lineno,
						     const char *fmt, ...)
{
  va_list ap;
  
  fprintf(stderr, "GeneralVoronoiS2 Exception: %s: %s: %d:", srcfile, function, lineno);

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

  fprintf(stderr, "\n");
}

generalvoronois2exception::~generalvoronois2exception()
{
}
