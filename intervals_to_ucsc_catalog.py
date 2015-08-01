#!/usr/bin/env python
"""
Convert genomic intervals to an html "catalog" for controlling the UCSC browser.
"""

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil
try:                from hashlib import md5 as md5_new
except ImportError: from md5     import new as md5_new


def usage(s=None):
	message = """

usage: cat <intervals> | intervals_to_ucsc_catalog [options]
  <filename1>              (usually required) file to write frame holder to
  <filename2>              (required) file to write catalog frame to
  --catalogonly            don't create a frame holder
  --url=<text>             base URL the catalog frame will reside at
                           (default is "current directory")
  --title=<text>           title for catalog frame
  --browser=<text>         base URL the genome browser
                           (default is http://genome.ucsc.edu/cgi-bin)
  --genome=<name>          reference genome
                           (default is hg19)
  --expand=<ratio>         (cumulative) add link to view intervals, expanded
  --center=<width>         (cumulative) add link to view intervals, centered
  --names[=<length>]       assign hash-based names to the intervals
  --show:numbers           number the intervals
  --show:comments          show any comments for each interval
  --framesplit=<fraction>  catalog/content frame split
  --origin=one             input intervals are origin-one, closed
  --origin=zero            input intervals are origin-zero, half-open
                           (this is the default)
  --head=<number>          limit the number of input intervals"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global frameTitle,catalogUrl
	global browserUrlBase,refGenome
	global nameLength
	global frameSplit
	global commentsSeparator
	global debug

	# parse args

	holderFilename    = None
	catalogFilename   = None
	catalogOnly       = False
	catalogUrlBase    = None
	frameTitle        = None
	browserUrlBase    = "http://genome.ucsc.edu/cgi-bin"
	refGenome         = "hg19"
	expansionSpecs    = []
	nameLength        = None
	showNumbers       = False
	showComments      = False
	commentsSeparator = None
	frameSplit        = (20,80)
	origin            = "zero"
	headLimit         = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--title=")):
			frameTitle = argVal
		elif (arg.startswith("--url=")):
			catalogUrlBase = argVal
			if (catalogUrlBase.endswith("/")): catalogUrlBase = catalogUrlBase[:-1]
		elif (arg == "--catalogonly"):
			catalogOnly = True
		elif (arg.startswith("--browser=")):
			browserUrlBase = argVal
			if (browserUrlBase.endswith("/")): browserUrlBase = browserUrlBase[:-1]
		elif (arg.startswith("--genome=")):
			refGenome = argVal
		elif (arg.startswith("--expand=")) or (arg.startswith("--expansion=")):
			try:               expansionRatio = int(argVal)
			except ValueError: expansionRatio = float(argVal)
			expansionSpecs += [("ratio",expansionRatio,argVal)]
		elif (arg.startswith("--center=")):
			expansionWidth = int_with_units(argVal)
			expansionSpecs += [("center",expansionWidth,argVal)]
		elif (arg == "--names"):
			nameLength = 3
		elif (arg == "--show:numbers"):
			showNumbers = True
		elif (arg == "--show:comments"):
			showComments      = True
			commentsSeparator = None
		elif (arg.startswith("--show:comments=")):
			showComments      = True
			commentsSeparator = argVal
		elif (arg.startswith("--framesplit=")) or (arg.startswith("--split=")):
			if (argVal.endswith("%")):
				f = float(argVal[:-1]) / 100
			elif ("," in argVal):
				(n,d) = argVal.split(",")
				(n,d) = (float(n),float(d))
				f = n / (n+d)
			elif ("/" in argVal):
				(n,d) = argVal.split("/")
				(n,d) = (float(n),float(d))
				f = n / d
			else:
				f = float(argVal)
			f = int(round(100*f))
			f = min(max(f,10),90)
			frameSplit = (f,100-f)
		elif (arg.startswith("--names=")):
			nameLength = int(argVal)
			nameLength = max( 3,nameLength)
			nameLength = min(10,nameLength)
		elif (arg.startswith("--head=")):
			headLimit = int_with_units(argVal)
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
		elif (holderFilename == None):
			holderFilename = arg
		elif (catalogFilename == None):
			catalogFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (catalogOnly):
		if (catalogFilename != None):
			usage("you can't give me the name of the frame holder file with --catalogonly")
		(holderFilename,catalogFilename) = (None,holderFilename)

	if (not catalogOnly) and (holderFilename == None):
		usage("you have to give me the names of the frame holder and catalog frame files")

	if (catalogFilename == None):
		usage("you have to give me the name of the catalog frame file")

	if (frameTitle == None): frameTitle = "(no title)"

	# process the intervals

	if (holderFilename != None): holderF = file(holderFilename,"wt")
	else:                        holderF = None

	catalogF = file(catalogFilename,"wt")
	if (catalogUrlBase == None): catalogUrl = catalogFilename.split("/")[-1]
	else:                        catalogUrl = catalogUrlBase + "/" + catalogFilename.split("/")[-1]

	intervalsPrinted = False
	iChrom = None

	intervalNum = 0
	for (chrom,start,end,comment) in read_intervals(stdin,origin=origin):
		intervalNum += 1
		if (headLimit != None) and (intervalNum > headLimit):
			print >>stderr, "limit of %s intervals reached" % (commatize(headLimit))
			break

		if (not intervalsPrinted):
 			intervalsPrinted = True
 			print >>catalogF, interval_frame_head()

		if (iChrom == None):
 			(iChrom,iStart,iEnd) = (chrom,start,end)

		if   (showNumbers):      number  = str(intervalNum)
		else:                    number  = None
		if   (not showComments): comment = None
		elif (comment == None):  comment = ""

		print >>catalogF, link_to_interval(chrom,start,end,number,comment,
		                                   expansions=expansionSpecs)

	if (intervalsPrinted):
		print >>catalogF, interval_frame_tail()
		if (holderF != None): print >>holderF, frame_holder(iChrom,iStart,iEnd)

	if (holderF != None): holderF.close()
	catalogF.close()


def frame_holder(iChrom,iStart,iEnd):
	intervalName = "%s:%d-%d" % (iChrom,iStart+1,iEnd)

	url = browserUrlBase \
	    + "/hgTracks?" \
	    + "db=" + refGenome \
	    + "&position=%s" % intervalName

	lines = []
	lines += ["<html><head><title>%s</title></head>" % frameTitle]
	lines += ["<frameset cols=\"%d%%,%d%%\">" % (frameSplit[0],frameSplit[1])]
	lines += ["  <frame name=\"catalog\" src=\"%s\">" % catalogUrl]
	lines += ["  <frame name=\"browser\" src=\"%s\">" % url]
	lines += ["</frameset>"]
	lines += ["</html>"]
	return "\n".join(lines)

def interval_frame_head():
	lines = []
	lines += ["<html>"]
	lines += ["<head>"]
	lines += ["<title>"]
	lines += [frameTitle]
	lines += ["</title>"]
	lines += ["</head>"]
	lines += [""]
	lines += ["<table border cellpadding=\"1\">"]
	return "\n".join(lines)

def interval_frame_tail():
	lines = []
	lines += ["</table>"]
	lines += [""]
	lines += ["</body>"]
	lines += ["</html>"]
	return "\n".join(lines)


def link_to_interval(chrom,start,end,number,comment,expansions=None):
	if (expansions == None): expansions = []
	intervalString = "%s:%s-%s" % (chrom,commatize(start+1),commatize(end))

	url = browserUrlBase \
	    + "/hgTracks?" \
	    + "db=" + refGenome \
	    + "&amp;position=%s%%3A%d-%d" % (chrom,start+1,end)

	fields = []

	# add number

	if (number != None):
		fields += ["<a href=\"%s\" target=browser>%s</a>" % (url,number)]

	# add name

	if (nameLength != None):
		hashVal = md5_new()
		hashVal.update(intervalString)
		hashVal  = md5_to_value(hashVal.hexdigest()[:25],26**nameLength)
		fakeName = value_to_name(hashVal,nameLength)
		if ("md5" in debug):
			print >>stderr, "%s -> %d -> %s" % (intervalString,hashVal,fakeName)

		fields += ["<a href=\"%s\" target=browser>%s</a>" % (url,fakeName)]

	# add primary interval link

	fields += ["<a href=\"%s\" target=browser>%s</a>" % (url,intervalString)]

	# add expanded interval link(s)

	for expansion in expansions:
		if (expansion[0] == "ratio"):
			(_,ratio,ratioStr) = expansion
			mid = (start+end) / 2
			exStart = int(mid + ((start - mid) * ratio))
			exEnd   = int(mid + ((end   - mid) * ratio))
			if (exStart < 0): exStart = 0

			exUrl = browserUrlBase \
				+ "/hgTracks?" \
				+ "db=" + refGenome \
				+ "&amp;position=%s%%3A%d-%d" % (chrom,exStart+1,exEnd)

			fields += ["<a href=\"%s\" target=browser>zoom out %sx</a>" % (exUrl,ratioStr)]

		if (expansion[0] == "center"):
			(_,width,widthStr) = expansion
			mid     = (start+end) / 2
			length  = end - start
			excess  = length - width
			exEnd   = end - excess/2
			exStart = max(0,exEnd-width)
			if (exStart < 0): exStart = 0

			exUrl = browserUrlBase \
				+ "/hgTracks?" \
				+ "db=" + refGenome \
				+ "&amp;position=%s%%3A%d-%d" % (chrom,exStart+1,exEnd)

			fields += ["<a href=\"%s\" target=browser>%s @ %s:%s</a>" \
					 % (exUrl,widthStr,chrom,mid)]

	# add comment

	if (comment != None):
		if (commentsSeparator != None):
			comment = comment.split(commentsSeparator)
			comment = "\n".join([c.strip() for c in comment])
		fields += [comment]

	return  "<tr>" + "".join(["<td>%s</td>" % s for s in fields]) + "</tr>"


def md5_to_value(s,modulus=None):
	v = int(s,16)
	if (modulus != None):
		v %= modulus
	return v


def value_to_name(val,numLetters):
	name = []
	for ix in xrange(numLetters):
		letter = val % 26
		val    = val / 26
		name += [chr(ord("a") + letter)]

	name.reverse()
	return "".join(name)


def read_intervals(f,origin="zero"):
	columnsNeeded = 3

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

		try:
			chrom = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (end < start): raise ValueError
			if (origin == "one"): start -= 1
			comment = None
			if (len(fields) > 3) and (fields[3].startswith("#")):
				comment = " ".join(fields[3:])[1:]
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		yield (chrom,start,end,comment)


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


# commatize--
#	Convert an integer to a string, with commas

def commatize(val):
	if (val >= 0):  sign      =  ""
	else:          (sign,val) = ("-",-val)

	val = str(val)

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return sign + val


if __name__ == "__main__": main()

