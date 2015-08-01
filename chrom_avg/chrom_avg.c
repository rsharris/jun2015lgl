// chrom_avg.c-- (chromosome coverage) read a list of genomic intervals, with
//               values, and report the average value for each position in
//               chromosomes of interest

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "utilities.h"

// program revision vitals (not the best way to do this!))

#define programVersionMajor    "0"
#define programVersionMinor    "0"
#define programVersionSubMinor "1"
#define programRevisionDate    "20140415"

//----------
//
// global data and types--
//
//----------

// linked list for chromosomes-of-interest

typedef struct spec
	{
	struct spec* next;			// next spec in a linked list
	char*	chrom;				// chromosome name
	u32		start;				// number of uninteresting bases at the start
								// .. of the chromosome
	u32		length;				// number of interesting bases in the chromosome
	void*	countVector;		// vector in which to accumulate count
	double*	sumVector;			// vector in which to accumulate value sum
	int		batchNumber;		// number of batches of this chromosome we've
								// .. seen
	} spec;

// command line options

spec* chromsOfInterest = NULL;
int   valColumn        = 4-1;
int   precision        = 0;
int   originOne        = false;
int   reportBatches    = false;
int   reportChroms     = false;

// counting parameters
// $$$ eventually I'd like to let the user set these to 1, 2 or 4 bytes
// $$$ to allow us to use less memory for certain cases

int   countBytes      = sizeof(u32);
u32   maxCount        = (u32) -1;

int   sumBytes        = sizeof(double);

//----------
//
// stuff for crude profiling
//
//----------

//#define useStandardClock  (define at build time)

#ifdef useStandardClock
#define read_clock() clock()
#define clocksPerSec CLOCKS_PER_SEC
#endif // useStandardClock

#ifndef useStandardClock
#include <sys/time.h>
#define read_clock() microsec_clock()
#define clocksPerSec 1000000
u64 microsec_clock (void);
u64 microsec_clock (void)
	{
	static int		failed = false;
	struct timeval	time1;
	int				err;

	if (failed)	return 0;	// (previous call to gettimeofday has failed)
	err = gettimeofday (&time1, NULL);
	if (err != 0) { failed = true;  return 0; }

	return (((u64) time1.tv_sec) * 1000000) + time1.tv_usec;
	}
#endif // not useStandardClock

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);

// private functions

static void  parse_options        (int _argc, char** _argv);
static spec* find_chromosome_spec (char* chrom);
static int   read_interval        (FILE* f,
                                   char* buffer, int bufferLen, int valCol,
                                   char** chrom, u32* start, u32* end,
                                   double* val);

//----------
//
// option parsing--
//
//----------

static void  usage    (char* message);
static void  chastise (const char* format, ...);

char* programName = "chrom_avg";


static void usage
   (char*	message)
	{
	if (message != NULL) fprintf (stderr, "%s\n", message);
	fprintf (stderr, "usage: %s <chromosome>[:<length>] [options]\n", programName);
	fprintf (stderr, "\n");
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  <chromosome>[:<length>]  accumulate coverage for this chromosome;\n");
	fprintf (stderr, "                           this is cumulative, many chromosomes can be given\n");
	fprintf (stderr, "  L=<length>               set length of all chromosomes (in bases) for which\n");
	fprintf (stderr, "                           length is not otherwise specified;  by default we\n");
	fprintf (stderr, "                           assume L=250M\n");
	fprintf (stderr, "                           L=0 means length *must* be set specifically for each\n");
	fprintf (stderr, "                           chromosome\n");
	fprintf (stderr, "  --value=<col>            input intervals contain a value in the specified\n");
	fprintf (stderr, "                           column; by default we assume this is in column 4\n");
	fprintf (stderr, "  --precision=<number>     number of digits to round average values to\n");
	fprintf (stderr, "  --origin=one             input/output intervals are origin-one, closed\n");
	fprintf (stderr, "  --origin=zero            input/output intervals are origin-zero, half-open\n");
	fprintf (stderr, "                           (this is the default)\n");
	fprintf (stderr, "  --progress               report each batch of the chromosome encountered\n");
	fprintf (stderr, "  --progress=chromosome    report each chromosome as we encounter it\n");
	fprintf (stderr, "  --version                report the program version and quit\n");
	exit (EXIT_FAILURE);
	}


static void chastise (const char* format, ...)
	{
	va_list	args;

	va_start (args, format);
	if (format != NULL)
		vfprintf (stderr, format, args);
	va_end (args);

	usage (NULL);
	}


static void parse_options
   (int			_argc,
	char**		_argv)
	{
	int			argc;
	char**		argv;
	char*		arg, *argVal, *argScan, *argScan2;
	int			tempInt;
	u32			chromStart = 0;
	u32			allChromLength = 0, chromLength;
	spec*		scanSpec, *tailSpec, *newSpec;

	// skip program name

	//programName = _argv[0];
	argv = _argv+1;  argc = _argc - 1;

	//////////
	// scan arguments
	//////////

	while (argc > 0)
		{
		arg    = argv[0];
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// L=<length>

		if ((strcmp_prefix (arg, "L=")   == 0)
		 || (strcmp_prefix (arg, "--L=") == 0))
			{
			tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
			if (tempInt < 0)
				chastise ("chromosome length can't be negative (\"%s\")\n", arg);
			allChromLength = (u32) tempInt;
			if (allChromLength == 0) allChromLength = (u32) -1;
			goto next_arg;
			}

		// --value=<col>

		if (strcmp_prefix (arg, "--value=") == 0)
			{
			valColumn = string_to_int (argVal) - 1;
			if (valColumn == -1)
				chastise ("value column can't be 0 (\"%s\")\n", arg);
			if (valColumn < 0)
				chastise ("value column can't be negative (\"%s\")\n", arg);
			if (valColumn < 3)
				chastise ("value column can't be 1, 2 or 3 (\"%s\")\n", arg);
			goto next_arg;
			}

		// --precision=<col>

		if (strcmp_prefix (arg, "--precision=") == 0)
			{
			precision = string_to_int (argVal);
			if (precision < 0)
				chastise ("precision can't be negative (\"%s\")\n", arg);
			goto next_arg;
			}

		// --origin=one, --origin=zero

		if ((strcmp (arg, "--origin=one") == 0)
		 || (strcmp (arg, "--origin=1")   == 0))
			{ originOne = true;  goto next_arg; }

		if ((strcmp (arg, "--origin=zero") == 0)
		 || (strcmp (arg, "--origin=0")   == 0))
			{ originOne = false;  goto next_arg; }

		// --progress and --progress=chromosome

		if (strcmp (arg, "--progress")  == 0)
			{ reportBatches = true;  goto next_arg; }

		if ((strcmp (arg, "--progress=chromosome")  == 0)
		 || (strcmp (arg, "--progress=chromosomes") == 0))
			{ reportChroms = true;  goto next_arg; }

		// --version

		if (strcmp (arg, "--version")  == 0)
			{
			fprintf (stderr, "%s (version %s.%s.%s released %s)\n",
			                 programName,
			                 programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);
			exit (EXIT_SUCCESS);
			}

		// unknown -- argument

		if (strcmp_prefix (arg, "--") == 0)
			chastise ("Can't understand \"%s\"\n", arg);

		// <chromosome>[:<length>] or (undocumented) <chromosome>[:<start>:<end>]
		// .. in this case <start> and <end> are origin-zero half-open, regardless
		// .. of any user setting

		argScan = strchr(arg,':');
		if (argScan == NULL) argScan2 = NULL;
		                else argScan2 = strchr(argScan+1,':');

		if (argScan == NULL)
			{
			chromStart  = 0;
			chromLength = 0;
			}
		else if (argScan2 == NULL)
			{
			chromStart   = 0;
			*(argScan++) = 0;
			chromLength  = string_to_u32 (argScan);
			}
		else
			{
			*(argScan++)  = 0;
			*(argScan2++) = 0;
			chromStart    = string_to_u32 (argScan);
			chromLength   = string_to_u32 (argScan2) - chromStart;
			}

		tailSpec = NULL;
		for (scanSpec=chromsOfInterest ; scanSpec!=NULL ; scanSpec=scanSpec->next)
			{
			tailSpec = scanSpec;
			if (strcmp (arg, scanSpec->chrom) == 0)
				chastise ("can't specify %s more than once\n", arg);
			}

		newSpec = (spec*) malloc (sizeof(spec));
		if (tailSpec == NULL) chromsOfInterest = newSpec;
						 else tailSpec->next   = newSpec;
		newSpec->next        = NULL;
		newSpec->chrom       = copy_string (arg);
		newSpec->start       = chromStart;
		newSpec->length      = chromLength;
		newSpec->countVector = NULL;
		newSpec->sumVector   = NULL;
		newSpec->batchNumber = 0;

	next_arg:
		argv++;  argc--;
		continue;
		}

	//////////
	// sanity checks
	//////////

	// make sure we got a chromosome-of-interest

	if (chromsOfInterest == NULL)
		chastise ("gotta give me some chromosome names\n");

	// assign default chromosome lengths

	if (allChromLength == 0)
		allChromLength = 250*1000*1000;
	else if (allChromLength == (u32) -1)
		allChromLength = 0;

	for (scanSpec=chromsOfInterest ; scanSpec!=NULL ; scanSpec=scanSpec->next)
		{ if (scanSpec->length == 0) scanSpec->length = allChromLength; }
	}

//----------
//
// main program--
//
//----------

int main
   (int		argc,
	char**	argv)
	{
	char	lineBuffer[1000];
	char	prevChrom[1000];
	u32*	cv32 = NULL;
	double*	sv = NULL;
	char*	chrom;
	spec*	chromSpec, *nextSpec;
	u32		start, end, o, adjStart, adjEnd;
	double	val, newVal;
	u32		ix;
	int		active, ok;
	s64		progressClock;
	float	secs;

	progressClock = -((s64) read_clock());

	parse_options (argc, argv);

	//////////
	// allocate memory
	//////////

	for (chromSpec=chromsOfInterest ; chromSpec!=NULL ; chromSpec=chromSpec->next)
		{
		if (chromSpec->length == 0) goto no_length_specified;
		chromSpec->countVector = calloc (chromSpec->length, countBytes);
		if (chromSpec->countVector == NULL) goto cant_allocate_count;

		chromSpec->sumVector = (double*) calloc (chromSpec->length, sumBytes);
		if (chromSpec->sumVector == NULL) goto cant_allocate_sum;
		sv = chromSpec->sumVector;
		for (ix=0 ; ix<chromSpec->length ; ix++)
			sv[ix] = 0.0;
		}

	//////////
	// process intervals
	//////////

	if (originOne) o = 1;
	          else o = 0;

	// read intervals and accumulate counts and values

	prevChrom[0] = 0;
	chromSpec    = NULL;

	while (true)
		{
		ok = read_interval (stdin, lineBuffer, sizeof(lineBuffer), valColumn,
	                        &chrom, &start, &end, &val);
		if (!ok) break;

		//fprintf (stderr, "%s %u %u %f\n", chrom, start, end, val);

		if (strcmp (chrom, prevChrom) != 0)
			{
			cv32 = NULL;
			sv   = NULL;
			chromSpec = find_chromosome_spec (chrom);
			if (chromSpec != NULL)
				{
				cv32 = (u32*) chromSpec->countVector;	// (this assumes countBytes==4)
				sv   = chromSpec->sumVector;
				chromSpec->batchNumber++;
				if (reportBatches || reportChroms)
					{
					progressClock += (s64) read_clock();
					secs = ((float) progressClock) / clocksPerSec;
					progressClock = -((s64) read_clock());

					fprintf (stderr, "(%s) progress: reading %s batch %d\n",
					                 duration_to_string (secs), chrom, chromSpec->batchNumber);
					}
				}
			else if (reportChroms)
				fprintf (stderr, "progress: ignoring %s\n", chrom);
			strncpy (prevChrom, chrom, sizeof(prevChrom));
			}

		if (chromSpec == NULL) continue;

		start -= o;
		adjStart = start;
		adjEnd   = end;

		if (chromSpec->start == 0)
			{
			// if only length has been specified, we *reject* intervals beyond
			// the end

			if (end > chromSpec->length) goto chrom_too_short;
			}
		else
			{
			// if start and end have been specified, we *ignore* intervals, or
			// protions of intervals, beyond the end

			if (end <= chromSpec->start) continue;
			adjEnd = end - chromSpec->start;
			if (start <= chromSpec->start) adjStart = 0;
			                          else adjStart = start - chromSpec->start;
			if (adjStart >= chromSpec->length) continue;
			if (adjEnd   >= chromSpec->length) adjEnd = chromSpec->length;
			}

		for (ix=adjStart ; ix<adjEnd ; ix++)
			{
			// nota bene: when maxCount is reached, we stop accumulating
			//            .. values for that position
			if (cv32[ix] < maxCount)
				{
				cv32[ix] += 1;
				sv[ix] += val;
				}
			}
		}

	// report intervals

	for (chromSpec=chromsOfInterest ; chromSpec!=NULL ; chromSpec=chromSpec->next)
		{
		cv32 = (u32*) chromSpec->countVector;	// (this assumes countBytes==4)
		sv   = chromSpec->sumVector;

		if (reportBatches || reportChroms)
			{
			progressClock += (s64) read_clock();
			secs = ((float) progressClock) / clocksPerSec;
			progressClock = -((s64) read_clock());

			fprintf (stderr, "(%s) progress: processing %s\n",
							 duration_to_string (secs), chromSpec->chrom);
			}

		active = false;
		start = 0;
		val = 0.0;
		for (ix=0 ; ix<chromSpec->length ; ix++)
			{
			if (cv32[ix] == 0)
				{
				if ((active) && (ix != start))
					printf ("%s\t%d\t%d\t%.*f\n",
					        chromSpec->chrom, chromSpec->start+start+o, chromSpec->start+ix,
					        precision, val);
				active = false;  start = 0;  val = 0.0;
				continue;
				}

			newVal = sv[ix] / cv32[ix];
			if (!active)
				{
				active = true;  start = ix;  val = newVal;
				continue;
				}

			if (newVal != val)
				{
				if (ix != start)
					printf ("%s\t%d\t%d\t%.*f\n",
					        chromSpec->chrom, chromSpec->start+start+o, chromSpec->start+ix,
					        precision, val);
				active = true;  start = ix;  val = newVal;
				continue;
				}
			}

		if ((active) && (chromSpec->length != start))
			printf ("%s\t%d\t%d\t%.*f\n",
			        chromSpec->chrom, chromSpec->start+start+o, chromSpec->start+chromSpec->length,
			        precision, val);
		}

	//////////
	// success
	//////////

	for (chromSpec=chromsOfInterest ; chromSpec!=NULL ; chromSpec=nextSpec)
		{
		nextSpec = chromSpec->next;
		if (chromSpec->chrom       != NULL) free (chromSpec->chrom);
		if (chromSpec->countVector != NULL) free (chromSpec->countVector);
		if (chromSpec->sumVector   != NULL) free (chromSpec->sumVector);
		free (chromSpec);
		}
	chromsOfInterest = NULL;

	return EXIT_SUCCESS;

	//////////
	// failure exits
	//////////

no_length_specified:
	fprintf (stderr, "no length was specified for %s\n",
	                 chromSpec->chrom);
	return EXIT_FAILURE;

cant_allocate_count:
	fprintf (stderr, "failed to allocate %d-entry counting vector for %s, %d bytes per entry\n",
	                 chromSpec->length, chromSpec->chrom, countBytes);
	return EXIT_FAILURE;

cant_allocate_sum:
	fprintf (stderr, "failed to allocate %d-entry summing vector for %s, %d bytes per entry\n",
	                 chromSpec->length, chromSpec->chrom, sumBytes);
	return EXIT_FAILURE;

chrom_too_short:
	fprintf (stderr, "%s %d %d is beyond the end of the chromosome (L=%d)\n",
	                 chrom, start, end, chromSpec->length);
	return EXIT_FAILURE;
	}

//----------
//
// find_chromosome_spec--
//	Locate the spec for a specific chromosome name.
//
//----------
//
// Arguments:
//	char*	chrom:	name of the chromosome to look for.
//
// Returns:
//	a pointer to the spec for the chromosome;  NULL if the chromosome is not
//	in our list.
//
//----------

static spec* find_chromosome_spec
   (char*	chrom)
	{
	spec*	scanSpec;

	for (scanSpec=chromsOfInterest ; scanSpec!=NULL ; scanSpec=scanSpec->next)
		{ if (strcmp (chrom, scanSpec->chrom) == 0) return scanSpec; }

	return NULL;
	}

//----------
//
// read_interval--
//	Read the next interval from a file.
//
//----------
//
// Arguments:
//	FILE*	f:			File to read from.
//	char*	buffer:		Buffer to read the line into.  Note that the caller
//						.. should not expect anything about the contents of
//						.. this buffer upon return.
//	int		bufferLen:	Number of bytes allocated for the buffer.
//	int		valCol:		The column that contains interval value.
//	char**	chrom:		Place to return a pointer to the chromosome.  The
//						.. returned value will point into the line buffer, and
//						.. to a zero-terminated string.
//	u32*	start:		Place to return the start.
//	u32*	end:		Place to return the end.
//	double*	val:		Place to return the value.
//
// Returns:
//	true if we were successful;  false if there are no more lines in the file.
//
//----------

static int read_interval
   (FILE*		f,
	char*		buffer,
	int			bufferLen,
	int			valCol,
	char**		_chrom,
	u32*		_start,
	u32*		_end,
	double*		_val)
	{
	static u32	lineNumber = 0;
	static int	missingEol = false;
	int			lineLen;
	char*		scan, *mark, *field;
	int			col;
	char*		chrom;
	u32			start, end;
	double		val;

	// read the next line

try_again:

	if (fgets (buffer, bufferLen, f) == NULL)
		return false;

	lineNumber++;

	// check for lines getting split by fgets (the final line in the file might
	// not have a newline, but no internal lines can be that way)

	if (missingEol) goto missing_eol;

	lineLen = strlen(buffer);
	if (lineLen != 0)
		missingEol = (buffer[lineLen-1] != '\n');

	// parse the line

	scan = skip_whitespace(buffer);
	if (*scan == 0)   goto try_again;  // empty line
	if (*scan == '#') goto try_again;  // comment line

	chrom = scan = buffer;
	if (*scan == ' ') goto no_chrom;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;

	if (*scan == 0) goto no_start;
	field = scan;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;
	start = string_to_u32 (field);

	if (*scan == 0) goto no_end;
	field = scan;
	mark = skip_darkspace(scan);
	scan = skip_whitespace(mark);
	if (*mark != 0) *mark = 0;
	end = string_to_u32 (field);

	for (col=3 ; col<=valCol ; col++)
		{
		if (*scan == 0) goto no_value;
		field = scan;
		mark = skip_darkspace(scan);
		scan = skip_whitespace(mark);
		}
	if (*mark != 0) *mark = 0;
	val = string_to_double (field);

	//////////
	// success
	//////////

	if (_chrom != NULL) *_chrom = chrom;
	if (_start != NULL) *_start = start;
	if (_end   != NULL) *_end   = end;
	if (_val   != NULL) *_val   = val;

	return true;

	//////////
	// failure exits
	//////////

missing_eol:
	fprintf (stderr, "problem at line %u, line is longer than internal buffer\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_chrom:
	fprintf (stderr, "problem at line %u, line contains no chromosome or begins with whitespace\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_start:
	fprintf (stderr, "problem at line %u, line contains no interval start\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_end:
	fprintf (stderr, "problem at line %u, line contains no interval end\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);

no_value:
	fprintf (stderr, "problem at line %u, line contains no interval value\n",
			 lineNumber-1);
	exit (EXIT_FAILURE);
	}

