//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: decodachrom.c								Author: Bob Harris
//
// Generic C application program
//
//----------

//----------
//
// decodachrom--
//  Reverse the process of encodachrom.
//
//  Encodachrom copies a file line by line, encoding chromsomes at the
//  beginning of line with a numeric code.  See encodachrom.c for more details
//  on the specific encoding.
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
#include <stdarg.h>				// standard C variable argument list stuff

//----------
//
// prototypes for local functions--
//
//----------

int main (int argc, char** argv);
static void ProcessFile (char* fName, FILE* f);
static unsigned long get_line(char **linep, unsigned long *np, FILE *fp);
static int need(unsigned long n, char **linep, unsigned long *np);
static void print_argv0(void);
static void fatalf(const char *fmt, ...);
static void ckfree(void *p);

//----------
//
// decodachrom--
//  Main program
//
//----------

void usage (void);
void usage (void)
	{
	fprintf (stderr, "decodachrom-- recover a file that has been encoded by encodachrom\n");
	fprintf (stderr, "usage: decodachrom [<text file> .. <text file>]\n");
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
//	Copy a file to stdout, decoding any chromosomes that appear at the start
//	of a line.
//
//----------

static char declet (int number);
static char declet (int number)
	{
	if ((number % 2) == 0) return 'A' + (number-100)/2;
	                  else return 'a' + (number-101)/2;
	}

static void ProcessFile
   (char*	fName,
	FILE*	f)
	{
	char*			buffer    = NULL;
	unsigned long	bufferLen = 0;
	unsigned long	lineLen   = 0;
	int				lineNumber;
	char*			s;
	int				chrNumber, num;

	// convert the file, line-by-line

	lineNumber = 0;

	while (true)
		{
		// read the next line;  note that get_line does not allow us to
		// distinguish an error from an end-of-file

		lineLen = get_line (&buffer, &bufferLen, f);
		if (lineLen == (unsigned long) -1) break;
		lineNumber++;
		if (lineLen == 0) continue;
		if (buffer[lineLen-1] == '\n')
			buffer[--lineLen] = 0;

		// if the line is empty, discard it

		if (buffer[0] == 0)
			{ continue; }

		// if the line begins with "#", remove that and copy the line as is

		if (buffer[0] == '#')
			{ printf ("%s\n", buffer+1);  continue; }

		// parse the chromosome number;  if the first field is not simply a
		// number, the file is errant

		s = buffer;

		chrNumber = 0;
		while (isdigit(*s))
			chrNumber = 10*chrNumber + *(s++)-'0';

		if ((*s != 0) && (!isspace (*s)) && (*s != '.'))
			goto failure;

		// print the decoded simple chromosome

		if      (chrNumber <    0) goto failure;
		else if (chrNumber ==   0) printf ("#");
		else if (chrNumber <= 100) printf ("chr%d", chrNumber-1);
		else if (chrNumber == 101) printf ("chrX");
		else if (chrNumber == 102) printf ("chrY");
		else if (chrNumber == 103) printf ("chrW");
		else if (chrNumber == 104) printf ("chrZ");
		else if (chrNumber == 105) printf ("chrM");
		else if (chrNumber <= 127) printf ("chr%c", 'A'+chrNumber-106);
		else if (chrNumber == 200) printf ("chr");
		else if (chrNumber == 300) ;
		else                       goto failure;

		// if the chromosome had extra stuff, decode it;  characters were
		// encoded as follows:
		//   0..9      => 00..09
		//   A..Z,a..z => 100..151  (e.g. A=>100, a=>101, B=>102, b=103, ...)
		//   others    => 200..455

		if (*s == '.')
			{
			s++;
			while (isdigit(*s))
				{
				num = *(s++)-'0';
				if (!isdigit(*s)) goto failure;
				num = 10*num + *(s++)-'0';
				if (num < 10)  { printf ("%d", num);  continue; }
				if (!isdigit(*s)) goto failure;
				num = 10*num + *(s++)-'0';
				if (num < 100) goto failure; // (can't happen)
				if (num < 152) { printf ("%c", declet (num));  continue; }
				if (num < 200) goto failure;
				if (num < 456) { printf ("%c", num-200);  continue; }
				goto failure;
				}

			if ((*s != 0) && (!isspace (*s)))
				goto failure;
			}

		// copy the remainder of the line

		printf ("%s\n", s);
		}

	// success

	ckfree (buffer);
	return;

	// failure

failure:
	if (fName == NULL)
		fatalf ("*** ERROR: improper input (line %d.%d). ***\n\n",
				lineNumber, s+1-buffer);
	else
		fatalf ("*** ERROR: improper input (%s, line %d.%d). ***\n\n",
				fName, lineNumber, s+1-buffer);
	}

//----------
//
// Support routines borrowed from multipipmaker, available at
//	http://www.bx.psu.edu/miller_lab
//
//----------

static char *argv0;

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

/* print_argv0 ---------------------------------------- print name of program */
static void print_argv0(void)
{
	if (argv0) {
	char *p = strrchr(argv0, '/');
	(void)fprintf(stderr, "%s: ", p ? p+1 : argv0);
	}
}

/* fatalf --------------------------------- format message, print it, and die */
void fatalf(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fflush(stdout);
	print_argv0();
	(void)vfprintf(stderr, fmt, ap);
	(void)fputc('\n', stderr);
	va_end(ap);
	exit(1);
}

static void ckfree(void *p)
{
	if (p) free(p);
}

