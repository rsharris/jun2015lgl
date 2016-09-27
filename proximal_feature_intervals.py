#!/usr/bin/env python
"""
Given two files of interval "features", find pairs of features close to each
other.
"""

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil


def usage(s=None):
	message = """

usage: cat <intervals> | proximal_feature_intervals <intervals filename>
  --proximity=<number>  (mandatory) distance threshold, in base pairs
  --mutual              only keep pairs that are mutually nearest
                        (by default we keep all pairs for which either is the
                        nearest of the other)
  --cutoff=<value>      discard any features with value no larger than <cutoff>
                        (by default, features needn't have values, and no
                        features are discarded)
  --positive            same as cutoff=0
  --origin=one          input intervals are origin-one, closed
  --origin=zero         input intervals are origin-zero, half-open
                        (this is the default)
                        (output intervals are always origin-zero, half-open)
  --head=<number>       limit the number of input lines

Input intervals are of the form <chrom> <start> <end>, and can be in random
order.  However, there can be no overlaps.

Note that we expect the number of intervals to be relatively small, so that we
can hold them in memory."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global headLimit,origin,mutuallyClosest
	global debug

	# parse args

	features2Filename = None
	maxDistance       = None
	mutuallyClosest   = False
	valueCutoff       = None
	origin            = "zero"
	headLimit         = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--proximity=")) or (arg.startswith("D=")):
			maxDistance = int_with_units(argVal)
		elif (arg == "--mutual"):
			mutuallyClosest = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_units(argVal)
		elif (arg.startswith("--cutoff=")):
			valueCutoff = float(argVal)
		elif (arg == "--positive"):
			valueCutoff = 0.0
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
		elif (features2Filename == None):
			features2Filename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (maxDistance == None):
		usage ("you need to give me a proximity threshold")

	if (features2Filename == None):
		usage ("you must give me a second set of features")

	# load the features

	(chromToFeatures1,order1) = read_features(stdin,cutoff=valueCutoff)

	f = file(features2Filename,"rt")
	(chromToFeatures2,order2) = read_features(f,filename=features2Filename,cutoff=valueCutoff)
	f.close()

	chromOrder = order1 + [chrom for chrom in order2 if (chrom not in order1)]

	# process the features

	for chrom in chromOrder:
		if (chrom not in chromToFeatures1): continue
		if (chrom not in chromToFeatures2): continue

		features1 = chromToFeatures1[chrom]
		features2 = chromToFeatures2[chrom]

		# collect all proximal pairs on this chromosome

		f1ToProximal = {}
		f2ToProximal = {}

		for pair in proximal_pairs(maxDistance,features1,features2):
			(start1,end1,start2,end2) = pair

			if   (start1 >= end2): d = 1 + start1 - end2
			elif (start2 >= end1): d = 1 + start2 - end1
			else:                  d = 0  # (they overlap)

			if ((start1,end1) not in f1ToProximal): f1ToProximal[(start1,end1)] = []
			f1ToProximal[(start1,end1)] += [(d,start2,end2)]

			if ((start2,end2) not in f2ToProximal): f2ToProximal[(start2,end2)] = []
			f2ToProximal[(start2,end2)] += [(d,start1,end1)]

		# for each feature in set 1, choose the closest match in set 2
		# $$$ this ignores ties!

		for (start1,end1) in f1ToProximal:
			f1ToProximal[(start1,end1)].sort()
			f1ToProximal[(start1,end1)] = f1ToProximal[(start1,end1)][0]

		# for each feature in set 2, choose the closest match in set 1

		for (start2,end2) in f2ToProximal:
			f2ToProximal[(start2,end2)].sort()
			f2ToProximal[(start2,end2)] = f2ToProximal[(start2,end2)][0]

		# report all features that are mutually closest

		if (mutuallyClosest):
			# report all features that are mutually closest
			for (start1,end1) in f1ToProximal:
				(d,start2,end2) = f1ToProximal[(start1,end1)]
				if (f2ToProximal[(start2,end2)] != (d,start1,end1)): continue
				print "%s\t%d\t%d\t%d\t%d" % (chrom,start1,end1,start2,end2)
		else:
			# report all features with its closest mate
			pairs = set()
			for (start1,end1) in f1ToProximal:
				(d,start2,end2) = f1ToProximal[(start1,end1)]
				pairs.add((start1,end1,start2,end2))

			for (start2,end2) in f2ToProximal:
				(d,start1,end1) = f2ToProximal[(start2,end2)]
				pairs.add((start1,end1,start2,end2))

			pairs = list(pairs)
			pairs.sort()
			for (start1,end1,start2,end2) in pairs:
				print "%s\t%d\t%d\t%d\t%d" % (chrom,start1,end1,start2,end2)


def proximal_pairs(maxDistance,features1,features2):
	if (features1 == []) or (features2 == []): return

	# sort the features into a common list, dilating set 2

	d = maxDistance
	features = [(start  ,end  ,1) for (start,end) in features1] \
	         + [(start-d,end+d,2) for (start,end) in features2]
	features.sort()

	# report any overlaps that came from different sets

	end1 = 0
	end2 = -d
	for (ix,(start,end,which)) in enumerate(features):
		if (which == 1):
			if ("scan" in debug):
				print >>stderr, "[%d] %d %d %d end2=%d" % (ix,which,start,end,end2)
			if (start < end2):
				for iy in xrange(ix-1,-1,-1):
					(s,e,w) = features[iy]
					if (w == 1): continue
					if ("scan" in debug):
						print >>stderr, "  scan [%d] %d %d %d" % (iy,w,s,e)
					if (e <= start): break
					if ("scan" in debug):
						print >>stderr, "  yield"
					yield (start,end,s+d,e-d)
			end1 = max(end,end1)
		else: # if (which == 2):
			if ("scan" in debug):
				print >>stderr, "[%d] %d %d %d end1=%d" % (ix,which,start,end,end1)
			if (start < end1):
				for iy in xrange(ix-1,-1,-1):
					(s,e,w) = features[iy]
					if (w == 2): continue
					if ("scan" in debug):
						print >>stderr, "  scan [%d] %d %d %d" % (iy,w,s,e)
					if (e <= start): break
					if ("scan" in debug):
						print >>stderr, "  yield"
					yield (s,e,start+d,end-d)
			end2 = max(end,end2)


def read_features(f,filename="(unnamed)",cutoff=None):
	chromOrder      = []
	chromToFeatures = {}

	if (cutoff == None):
		intervalNum = 0
		for (lineNumber,chrom,start,end) in read_intervals(f):
			intervalNum += 1
			if (headLimit != None) and (intervalNum > headLimit):
				print >>stderr, "limit of %s intervals reached in %s" \
							  % (commatize(headLimit),filename)
				break

			if (chrom not in chromToFeatures):
				chromOrder += [chrom]
				chromToFeatures[chrom] =  [(start,end,lineNumber)]
			else:
				chromToFeatures[chrom] += [(start,end,lineNumber)]
	else:
		intervalNum = 0
		for (lineNumber,chrom,start,end,val) in read_intervals(f,withValues=True):
			intervalNum += 1
			if (headLimit != None) and (intervalNum > headLimit):
				print >>stderr, "limit of %s intervals reached in %s" \
							  % (commatize(headLimit),filename)
				break

			if (val <= cutoff): continue

			if (chrom not in chromToFeatures):
				chromOrder += [chrom]
				chromToFeatures[chrom] =  [(start,end,lineNumber)]
			else:
				chromToFeatures[chrom] += [(start,end,lineNumber)]

	# positionally sort them, and verify that there are no overlaps

	for chrom in chromOrder:
		features = chromToFeatures[chrom]
		features.sort()

		(prevStart,prevEnd,prevLine) = features[0]

		for (start,end,lineNumber) in features[1:]:
			if (start < prevEnd):
				assert (False), "intervals overlap on \"%s\"\n" \
				                "(%d,%d) overlaps (%d,%d) at lines %d and %d" \
				              % (chrom,prevStart,prevEnd,start,end,prevLine,lineNumber)
			(prevStart,prevEnd,prevLine) = (start,end,lineNumber)

		chromToFeatures[chrom] = [(start,end) for (start,end,_) in features]

	return (chromToFeatures,chromOrder)


def read_intervals(f,withValues=False):
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
			if (end < start): raise ValueError
			if (origin == "one"): start -= 1
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		if (not withValues):
			yield (lineNumber,chrom,start,end)
			continue

		try:
			if (numFields < 4): val = 1
			else:               val = float(fields[3])
			yield (lineNumber,chrom,start,end,val)
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)


# parse a string as an integer, allowing units (e.g. "3.2M")

def int_with_units(s):
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


# int_or_float--

def int_or_float(s):
	try:               return int(s)
	except ValueError: return float(s)


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

