// utilities.c-- miscellaneous utility functions.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <inttypes.h>
#include "utilities.h"

//----------
//
// copy_string--
//	Create (in the heap) a copy of a string or a prefix of a string.
//
//----------
//
// Arguments:
//	const char*	s:	The string to copy.
//	int			n:	(copy_prefix only) the number of characters to copy.
//
// Returns:
//	A pointer to new string;  failures result in program termination.
//
//----------

char* copy_string
   (const char*	s)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc (strlen(s) + 1);
	if (ss == NULL)
		{
		fprintf (stderr, "failed to allocate %lld bytes to copy \"%s\"\n",
		                 (long long) (strlen(s)+1), s);
		exit (EXIT_FAILURE);
		}

	return strcpy (/*to*/ ss, /*from*/ s);
	}

//----------
//
// strcmp_prefix--
//	Determine if a string contains another as a prefix.
//
//----------
//
// Arguments:
//	const char*	str1:	The string.
//	const char*	str2:	The prefix string.
//
// Returns:
//	The same as strcmp(prefix1,str2) would, where prefix1 is str1 truncated
//	to be no longer than str2.
//
//----------

int strcmp_prefix
   (const char*	str1,
	const char*	str2)
	{
	return strncmp (str1, str2, strlen(str2));
	}

//----------
//
// string_to_int, string_to_u32--
//	Parse a string for the integer value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in program
//	termination.
//
//----------

int string_to_int
   (const char*	s)
	{
	char*		ss;
	int			v;
	char		extra;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// convert to number

	if (sscanf (ss, "%d%c", &v, &extra) != 1) goto not_an_integer;

	// make sure signs match

	if ((v < 0) && (*ss != '-')) goto out_of_range;
	if ((v > 0) && (*ss == '-')) goto out_of_range;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fprintf (stderr, "an empty string is not an integer\n");
	exit (EXIT_FAILURE);

not_an_integer:
	fprintf (stderr, "\"%s\" is not an integer\n", s);
	exit (EXIT_FAILURE);

out_of_range:
	fprintf (stderr, "\"%s\" is outside the range of a signed integer\n", s);
	exit (EXIT_FAILURE);

	return 0;
	}


int string_to_u32
   (const char*	s)
	{
	char*		ss;
	u32			v;
	char		extra;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// convert to number

	if (*ss == '-') goto not_an_integer;
	if (sscanf (ss, "%u%c", &v, &extra) != 1) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fprintf (stderr, "an empty string is not an unsigned integer\n");
	exit (EXIT_FAILURE);

not_an_integer:
	fprintf (stderr, "\"%s\" is not an unsigned integer\n", s);
	exit (EXIT_FAILURE);

	return 0;
	}

//----------
//
// string_to_unitized_int, string_to_unitized_int64--
//	Parse a string for the integer value it contains, allowing K, M, and G
//	suffixes.
//
//----------
//
// Arguments:
//	const char*	s:		The string to parse.
//	int byThousands:	true  => K means one thousand
//						false => K means 1,024.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything (except for an opptional suffix) other than a valid integer--
//	failures result in fatality.
//
//----------

int string_to_unitized_int
   (const char*	s,
	int			byThousands)
	{
	char		ss[20];
	int			len = strlen (s);
	char*		parseMe;
	int			v;
	float		vf;
	char		extra;
	int			mult;
	int			isFloat;

	mult = 1;

	if (len >= (int) sizeof (ss))
		parseMe = (char*) s;
	else
		{
		parseMe = ss;
		strcpy (ss, s);

		if (len > 0)
			{
			switch (ss[len-1])
				{
				case 'K': case 'k':
					mult = (byThousands)? 1000 : 1024;
					break;
				case 'M': case 'm':
					mult = (byThousands)? 1000000 : 1024L * 1024L;
					break;
				case 'G': case 'g':
					mult = (byThousands)? 1000000000 : 1024L * 1024L * 1024L;
					break;
				}

			if (mult != 1)
				ss[len-1] = 0;
			}
		}

	isFloat = false;
	if (sscanf (parseMe, "%d%c", &v, &extra) != 1)
		{
		if (sscanf (parseMe, "%f%c", &vf, &extra) != 1) goto bad;
		isFloat = true;
		}

	if (isFloat)
		{
		if ((vf > 0) && ( vf*mult > INT_MAX)) goto overflow;
		if ((vf < 0) && (-vf*mult > INT_MAX)) goto overflow;
		v = (vf * mult) + .5;
		}
	else if (mult != 1)
		{
		if ((v > 0) && ( v > INT_MAX / mult)) goto overflow;
		if ((v < 0) && (-v > INT_MAX / mult)) goto overflow;
		v *= mult;
		}

	return v;

bad:
	fprintf (stderr, "\"%s\" is not an integer\n", s);
	exit (EXIT_FAILURE);

overflow:
	fprintf (stderr, "\"%s\" is out of range for an integer\n", s);
	exit (EXIT_FAILURE);

	return 0;
	}

//----------
//
// string_to_double--
//	Parse a string for the double floating point value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The value of the string.  Note that the string *must not* contain anything
//	other than a valid number-- failures result in fatality.
//
//----------

double string_to_double
   (const char*	s)
	{
	char*		ss;
	double		v;
	char		extra;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// convert to number

	if (sscanf (s, "%lf%c", &v, &extra) != 1) goto not_a_number;

	return v;

empty_string:
	fprintf (stderr, "an empty string is not a number\n");
	exit (EXIT_FAILURE);

not_a_number:
	fprintf (stderr, "\"%s\" is not a number\n", s);
	exit (EXIT_FAILURE);
	}

//----------
//
// skip_whitespace--
//	Skip characters until we get something that ain't whitespace.
// skip_darkspace--
//	Skip characters until we get something that ain't darkspace.
//
//----------
//
// Arguments:
//	char*	s:	The string to read.
//
// Returns:
//	Pointer to the first character at or beyond s that meets the stopping
//	criteria.  Note that we never scan beyond the end of the string.
//
//----------

char* skip_whitespace (char* s)
	{ while ((*s != 0) && (isspace (*s))) s++;  return s; }

char* skip_darkspace (char* s)
	{ while ((*s != 0) && (!isspace (*s))) s++;  return s; }

//----------
//
// duration_to_string--
//	Convert a time (duration) to a string.
//
//----------
//
// Arguments:
//	float	seconds:	The duration to convert.
//
// Returns:
//	a pointer to a string representation of the duration;  note that the
//	buffer holding this string is private to this function.
//
//----------

char* duration_to_string
   (float	seconds)
	{
	static	char buffer[1000];
	int		hours, minutes;

	if (seconds < 60)
		sprintf (buffer, "%.3fs", seconds);
	else if (seconds < 3600)
		{
		minutes =  seconds / 60;
		seconds -= 60 * minutes;
		sprintf (buffer, "%dm%06.3fs", minutes, seconds);
		}
	else
		{
		minutes =  seconds / 60;
		seconds -= 60 * minutes;
		hours   =  minutes / 60;
		minutes -= 60 * hours;
		sprintf (buffer, "%dh%02dm%06.3fs", hours, minutes, seconds);
		}

	return buffer;
	}
