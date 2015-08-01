#!/usr/bin/env python
"""
"Close" genomic intervals (in the sense of Minkowski/morphology set operations).

see en.wikipedia.org/wiki/Closing_(morphology)
"""

from sys  import argv,stdin,stderr,exit
from math import ceil


def usage(s=None):
	message = """
usage: cat intervals_file | close_intervals <length> [options] > intervals_file
  <length>       the number of bases to 'close' between intervals
  --origin=one   intervals are origin-one, closed
  --origin=zero  intervals are origin-zero, half-open
                 (this is the default)

  Note that we allow incoming intervals to extend beyond the end of a
  chromosome (and thus output intervals might also)."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse args

	closingLength = None
	origin        = "zero"
	debug         = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--origin=")):
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
		elif (closingLength == None):
			closingLength = int_with_unit(arg)
			if (closingLength < 0): usage("length must be non-negative")
		else:
			usage("unrecognized option: %s" % arg)

	if (closingLength == None):
		usage("you must provide the length")

	# collect the intervals

	chromToIntervals = {}
	chroms = []

	for (chrom,start,end) in read_intervals(stdin,origin=origin):
		if (chrom not in chromToIntervals):
			chromToIntervals[chrom] = []
			chroms += [chrom]
		chromToIntervals[chrom] += [(start,end)]

	# perform dilation, then erosion (and merge)

	dilationLength = (closingLength+1) / 2
	erosionLength  = dilationLength

	for chrom in chroms:
		intervals = [(start-dilationLength,end+dilationLength)
		                  for (start,end) in chromToIntervals[chrom]]

		for (start,end) in non_overlapping_intervals(intervals):
			start += erosionLength
			end   -= erosionLength
			if (origin == "one"): start += 1
			print "%s %d %d" % (chrom,start,end)


def read_intervals(f,origin="zero"):
	numFields = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (numFields == None):
			assert (len(fields) >= 3), \
			       "not enough fields at line %d (%d, expected at least %d)" \
			     % (lineNumber,len(fields),3)
			numFields = len(fields)
		else:
			assert (len(fields) == numFields), \
			       "inconsistent number of fields at line %d (%d, expected %d)" \
			     % (lineNumber,len(fields),numFields)

		try:
			chrom = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (end < start): raise ValueError
			if (origin == "one"): start -= 1
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		yield (chrom,start,end)


# merge into non-overlapping intervals

def non_overlapping_intervals(intervals):
	intervals.sort()
	overlaps = []
	start = end = None
	for (s,e) in intervals:
		if (start == None):
			(start,end) = (s,e)
		elif (s < end):
			end = max(end,e)
		else:
			overlaps += [(start,end)]
			(start,end) = (s,e)

	if (start != None):
		overlaps += [(start,end)]

	return overlaps


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
