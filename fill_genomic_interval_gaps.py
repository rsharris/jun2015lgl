#!/usr/bin/env python

from sys import argv,stdin,stderr,exit


def usage(s=None):
	message = """
usage: cat intervals | fill_genomic_interval_gaps [options] > intervals
  --chromosomes=<filename>  read chromosome names and lengths from a file
  --origin=0                intervals are origin-zero, half-open (default)
  --origin=1                intervals are origin-one, closed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	chromsFilename = None
	origin         = "zero"
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1].strip()

		if (arg.startswith("--chromosomes=")) or (arg.startswith("--chroms=")):
			chromsFilename = argVal
		elif (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			assert (origin in ["zero","one"]), "can't understand %s" % arg
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (chromsFilename == None):
		usage("you have to give me a chromosome lengths file")

	# read the chromosome lengths file

	chromToLength = {}
	chroms = []

	f = file(chromsFilename,"rt")
	for (chrom,length) in name_and_length(f):
		assert (chrom not in chromToLength), \
		       "%s occurs twice in %s" % (chrom,lengthsFile)
		chromToLength[chrom] = length
		chroms += [chrom]

	f.close()

	# process the intervals

	chromToIntervals = {}

	for (chrom,start,end,val) in read_intervals(stdin,origin=origin):
		if (chrom not in chromToIntervals):
			assert (chrom in chromToLength), \
			       "%s is in input but not in %s" % (chrom,lengthsFile)
			chromToIntervals[chrom] = []
		chromToIntervals[chrom] += [(start,end,val)]

	for chrom in chroms:
		if (chrom not in chromToIntervals):
			start = 0
			if (origin == "one"): start += 1
			print "%s\t%d\t%d\t%s" % (chrom,start,chromToLength[chrom],"0")
			continue

		intervals = chromToIntervals[chrom]
		intervals.sort()

		prevEnd = 0
		for (start,end,val) in intervals:
			if (start < prevEnd):
				if (origin == "one"): start += 1
				assert (False), "overlapping intervals on %s: ?-%d and %d-%d" \
				              % (chrom,prevEnd,start,end)

			if (prevEnd < start):
				if (origin == "one"): prevEnd += 1
				print "%s\t%d\t%d\t%s" % (chrom,prevEnd,start,"0")

			if (origin == "one"): start += 1
			print "%s\t%d\t%d\t%s" % (chrom,start,end,val)
			prevEnd = end

		if (prevEnd < chromToLength[chrom]):
			if (origin == "one"): prevEnd += 1
			print "%s\t%d\t%d\t%s" % (chrom,prevEnd,chromToLength[chrom],"0")


# returns the next interval as (chrom,start,end)

def read_intervals(f,origin="zero"):
	columnsNeeded = 3
	numFields = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= columnsNeeded), \
		       "not enough fields at line %d (%d, expected at least %d)" \
		     % (lineNumber,len(fields),columnsNeeded)

		if (numFields == None):
			numFields = len(fields)
		else:
			assert (len(fields) == numFields), \
			       "inconsistent number of fields at line %d (%d, expected %d)" \
			     % (lineNumber,len(fields),numFields)

		try:
			chrom = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (numFields < 4): val = 1
			else:               val = float(fields[3])
			if (end < start): raise ValueError
			if (origin == "one"): start -= 1
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		yield (chrom,start,end,val)


# yields the next name,length pair

def name_and_length(f):

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()

		fields = line.split()

		assert (len(fields) == 2), \
		      "inconsistent number of fields at line %d (%d, expected %d)" \
		    % (lineNumber,len(fields),2)

		try:
			name   =     fields[0]
			length = int(fields[1])
		except ValueError:
			assert (False), "bad length at line %d\n%s" % (lineNumber,line)

		yield (name,length)


if __name__ == "__main__": main()
