#ifndef utilities_H				// (prevent multiple inclusion)
#define utilities_H

#include <inttypes.h>
typedef int32_t  s32;
typedef uint32_t u32;
typedef int64_t  s64;
typedef uint64_t u64;

// macro to convince gnu c compiler not to complain about unusued function
// arguments

#ifdef __GNUC__
#define arg_dont_complain(arg) arg __attribute__ ((unused))
#else
#define arg_dont_complain(arg) arg
#endif // __GNUC__


// functions in this module

char*  copy_string            (const char* s);
int    strcmp_prefix          (const char* str1, const char* str2);
int    string_to_int          (const char* s);
int    string_to_u32          (const char* s);
int    string_to_unitized_int (const char* s, int byThousands);
double string_to_double       (const char* s);
char*  skip_whitespace        (char* s);
char*  skip_darkspace         (char* s);
char*  duration_to_string     (float seconds);

#endif // utilities_H
