#!/usr/bin/python
"""
Collect named tags.
"""

from sys  import argv,stdin,stderr,exit
from math import ceil


def usage(s=None):
	message = """
usage: cat items_file | collect_tags [options]
  --separator=<separator>  separator for tags
                           (default is comma)
  --head=<number>          limit the number of input lines

The input consists of lines of (group,tag) pairs.  Additional columns are
ignored.

	E88BQJZ01 A 8894
	E88BQJZ01 B 140509
	E88BQJZ01 C 4638
	E88BQJZ02 B 135134
	E88BQJZ02 C 4274
	FH0VK6D01 A 40470
	FH0VK6D01 C 22830
	FH0VK6D02 A 47624

Output looks like this:

	E88BQJZ01 A,B,C
	E88BQJZ02 B,C
	FH0VK6D01 A,C
	FH0VK6D02 A"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

    # parse the command line

	separator = ","
	headLimit = None
	debug     = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--separator=")) or (arg.startswith("--sep=")):
			separator = argVal
			if   (separator == "tab"):   separator = "\t"
			elif (separator == "space"): separator = " "
			elif (separator == "none"):  separator = ""
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# read the items

	keys = []
	keyToTags = {}

	lineNum = 0
	for line in stdin:
		lineNum += 1
		if (headLimit != None) and (lineNum > headLimit):
			print >>stderr, "limit of %s lines reached" % (commatize(headLimit))
			break

		fields = line.split()
		assert (len(fields) >= 2), \
		       "not enough fields in line %d (expected at least 2)\n%s" \
		     % (lineNumber,line)
		key = fields[0]
		tag = fields[1]

		if (key not in keyToTags):
			keys += [key]
			keyToTags[key] = [tag]
		elif (tag not in keyToTags[key]):
			keyToTags[key] += [tag]

	# report the items

	for key in keys:
		print "%s\t%s" % (key,separator.join([tag for tag in keyToTags[key]]))


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
