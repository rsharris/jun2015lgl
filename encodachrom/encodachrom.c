//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: encodachrom.c								Author: Bob Harris
//
// Generic C application program
//
//----------

//----------
//
// encodachrom--
//	Convert a file containing lines tagged with chromosomes (such as a typical
//	genomic intervals file) to one in which the chromosomes are represented by
//	numbers (integers and/or reals).  This is useful as a precursor to sorting.
//	The companion program decodachrom restores the original lines.
//
//  Typical use:
//		encodachrom farf.ranges | env LC_ALL=C sort -n -k1 -k2 | decodachrom > farf.sorted
//
//	The input file is expected to have lines beginning with (usually) "chr"
//	followed by a chromosome number or letter.  The chromosomes are encoded
//	with a number (in a way that should work for all species).  The rest of the
//	line is simply copied.  The input may also contain comment lines, which
//	begin with "#".
//
//	Chromosomes are encoded as follows (this causes the sex chromosomes to be
//	ordered X,Y,W,Z, and, along with mitochondrial M, to appear before the
//	other letters):
//		#              =>  0          (this is a comment)
//		chr0 .. chr99  =>  1..100
//		chrX           =>  101
//		chrY           =>  102
//		chrW           =>  103
//		chrZ           =>  104
//		chrM           =>  105
//		chrA .. chrV   =>  106..127   (with a hole at 118 where chrM is absent)
//		other          =>  200 or 300 (see non-chr discussion below)
//
//	If the chromosome is not among those above, the extra stuff is encoded
//	following a decimal point, encoded in a way that will cause sort to put
//	numbers ahead of alphabetics.  Specifically, each character is encoded
//	according to this table:
//		0..9      => 00..09
//		A..Z,a..z => 100..151  (e.g. A=>100, a=>101, B=>102, b=103, ...)
//		others    => 200..455  (with holes where letters and digits are absent)
//  Note that we aren't trying to achieve an efficient coding, just one that
//  the unix sort command will sort correctly and which is invertible.
//
//	If the chromsome does not begin with "chr" or if the chromsome name doesn't
//	start with a letter or digit (e.g. "chr@"), it is encoded as 200 or 300
//	followed by a post-decimal part encoding its entirety (with the same code
//	as above).  200 is used for the non-letter non-digit case.  300 is used for
//	the non-chr case.
//
//	One 'shortcoming' of the above encoding is if there are numbered chromosomes
//	with more than two digits.  For example, chr800 through chr809 will appear,
//	in sorted output, between chr80 and chr81.  This could be corrected with a
//	simple alteration to the encoding, but until I see three digit chromsomes
//	I'm not going to worry about it.
//
//----------

static const char rcsid[]=
"";

//----------
//
// other files--
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include <ctype.h>				// standard C upper/lower stuff

//----------
//
// prototypes for local functions--
//
//----------

int main (int argc, char** argv);
static void ProcessFile (char* fName, FILE* f);
static char* skip_whitespace (const char* s);
static int strncmplax (const char* s1, const char* s2, size_t n);
static unsigned long get_line(char **linep, unsigned long *np, FILE *fp);
static int need(unsigned long n, char **linep, unsigned long *np);
static void ckfree(void *p);

//----------
//
// encodachrom--
//  Main program
//
//----------

void usage (void);
void usage (void)
	{
	fprintf (stderr, "encodachrom-- numerically encode chromosomes in a range list\n");
	fprintf (stderr, "usage: encodachrom [<text file> .. <text file>]\n");
	fprintf (stderr, "input is from file(s) or stdin\n");
	fprintf (stderr, "output is to stdout\n");
	exit(1);
	}

int main
   (int		argc,
	char**	argv)
	{
	FILE*	f = NULL;
	int		argn;

	// process input file(s)

	if (argc == 0)
		usage ();

	if (argc == 1)	
		ProcessFile (NULL, stdin);
	else
		{
		for (argn=1 ; argn<argc ; argn++)
			{
		    f = fopen (argv[argn], "r");
		    if (f == NULL)
		    	{
		    	fprintf (stderr, "Can't open \"%s\"\n", argv[argn]);
		    	return 0;
		    	}
			ProcessFile (argv[argn], f);
			fclose (f);
		    }
	    }

	return 0;
	}

//----------
//
// ProcessFile--
//	Copy a file to stdout, encoding any chromosomes that appear at the start
//	of a line.
//
//----------

static void ProcessFile
   (char*	fName,
	FILE*	f)
	{
	char*			buffer    = NULL;
	unsigned long	bufferLen = 0;
	unsigned long	lineLen   = 0;
	char*			s;
	int				chrNumber;
	char			chrCh;

	// convert the file, line-by-line

	while (true)
		{
		// read the next line;  note that get_line does not allow us to
		// distinguish an error from an end-of-file

		lineLen = get_line (&buffer, &bufferLen, f);
		if (lineLen == (unsigned long) -1) break;
		if (lineLen == 0) continue;
		if (buffer[lineLen-1] == '\n')
			buffer[--lineLen] = 0;

		// if the line begins with "#" (skipping leading whitespace), encode it
		// all as post-decimal point (with a prefix of 0)

		s = skip_whitespace (buffer);

		if (*s == '#')
			{ s++;  printf ("0");  goto fractional_part; }

		// if the line doesn't begin with "chr" encode it all as post-decimal
		// point (with a prefix of 300)

		if (strncmplax (s, "chr", 3) != 0)
			{ printf ("300");  goto fractional_part; }

		// similarly, if the chromsome name doesn't begin with a digit or a
		// letter, encode it all as post-decimal point (with a prefix of 200)

		s += 3;

		if ((!isdigit(*s)) && (!isalpha(*s)))
			{ printf ("200");  goto fractional_part; }

		// otherwise, encode the chromosome as an integer

		if (isalpha(*s))
			{
			chrCh = toupper(*(s++));
			switch (chrCh)
				{
				case 'X': chrNumber = 101; break;
				case 'Y': chrNumber = 102; break;
				case 'W': chrNumber = 103; break;
				case 'Z': chrNumber = 104; break;
				case 'M': chrNumber = 105; break;
				default:  chrNumber = 106 + chrCh-'A';
				}
			}
		else
			{
			chrNumber = *(s++) - '0';
			if (isdigit(*s))
				chrNumber = 10*chrNumber + *(s++)-'0';
			chrNumber += 1;
			}

		// print the encoded simple chromosome name (the integer)

		printf ("%d", chrNumber);

		// if the chromosome name has extra stuff, encode it as post-decimal
		// point

		if ((*s != 0) && (!isspace(*s)))
			{
		fractional_part:
			printf (".");
			for ( ; (*s!=0)&&(!isspace(*s)) ; s++)
				{
				if (isdigit(*s))
					printf ("0%c", *s);
				else  if (isupper(*s))
					printf ("%03d", 100+2*(*s-'A'));
				else  if (islower(*s))
					printf ("%03d", 101+2*(*s-'a'));
				else
					printf ("%03d", 200+(*s));
				}
			}

		// print the remainder of the line

		printf ("%s\n", s);
		}

	ckfree (buffer);
	}

//----------
//
// skip_whitespace--
//  Skip to the next non-whitespace character in a string.
//
//----------
//
// Arguments:
//  const char*	s: The string.
//
// Returns:
//	A pointer to the next such character in a string (or to the end of the
//	string.
//
//----------

static char* skip_whitespace
   (const char*	s)
	{
	for ( ; (*s!=0)&&(isspace(*s)) ; s++)
		;

	return (char*) s;
	}

//----------
//
// strncmplax--
//	Compare the prefixes of two strings, using lax case matching, and
//	determine which precedes the other alphabetically.
//
//----------
//
// Arguments:
//	const char*	s1:	One string.
//	const char*	s2:	The other string.
//	size_t		n:	The length of the prefixes to compare.
//
// Returns:
//	anything < 0 => s1's n-character prefix precedes s2's alphabetically.
//	0            => s1 and s2 have identical n-character prefixes, except for
//	                .. possible case mismatches.
//	anything > 0 => s2's n-character prefix precedes s1's alphabetically.
//
//----------

static int strncmplax
   (const char*	s1,
	const char*	s2,
	size_t		n)
	{
	char		c1, c2;

	if (n <= 0)
		return 0;

	do
		{
		c1 = *(s1++); c1 = tolower (c1);
		c2 = *(s2++); c2 = tolower (c2);
		} while ((c1 == c2) && (c1 != 0) && (--n > 0));

	return c1 - c2;
	}

//----------
//
// Support routines borrowed from multipipmaker, available at
//	http://www.bx.psu.edu/miller_lab
//
//----------

// get_line -- read a line into malloc()ed memory,
// realloc()ing more if necessary.
// plug compatable with gnu getline().

static unsigned long get_line(char **linep, unsigned long *np, FILE *fp)
{
    int ch;
    unsigned long n = 0;
 
    while ((ch = fgetc(fp)) != -1) {
	if (need(n+1,linep,np)) return -1;
	(*linep)[n++] = ch;
	if (ch == '\n') break;
    }
    if (need(n+1,linep,np)) return -1;
    (*linep)[n] = 0;
    if (n == 0 && ch == -1) return -1;
    return n;
}

static int need(unsigned long n, char **linep, unsigned long *np)
{
    char *p;
    if (*np <= n) {
	*np += (*np >> 5) + 16;
	if ((p = realloc(*linep, *np)) == 0)
	    return -1;
	*linep = p;
    }
    return 0;
}

static void ckfree(void *p)
{
	if (p) free(p);
}

