#!/usr/bin/python
"""
Remove duplicate lines in a file, keeping only the first occurence of a line.
"""

from sys  import argv,stdin,stderr,exit
from math import ceil


def usage(s=None):
	message = """
usage: cat items_file | keep_first [options] > items_file
  --head=<number>        limit the number of input lines
  --progress=<number>    periodically report how many lines we've read"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

    # parse the command line

	headLimit      = None
	reportProgress = None
	debug          = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# read the items

	lineSeen = {}

	lineNum = 0
	for line in stdin:
		lineNum += 1
		if (headLimit != None) and (lineNum > headLimit):
			print >>stderr, "limit of %s lines reached" % (commatize(headLimit))
			break
		if (reportProgress != None) and (lineNum % reportProgress == 0):
			print >>stderr, "progress: line %s" % (commatize(lineNum))

		line = line.rstrip()
		if (line in lineSeen): continue
		lineSeen[line] = True
		print line


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


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):
		(val,suffix) = val.split(".",1)
		suffix = "." + suffix

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix


if __name__ == "__main__": main()
