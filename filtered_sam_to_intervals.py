#!/usr/bin/env python
"""
Filter and convert a SAM file to a list of intervals/reads

References:
  [1] The SAM Format Specification (samtools.github.io/hts-specs/SAMv1.pdf)
  [2] Using eval() safely in python (lybniz2.sourceforge.net/safeeval.html)
"""

from sys        import argv,stdin,stderr,exit
from math       import *
from re         import compile
from sam_reader import samFieldToColumn,SAM_MIN_COLUMNS,SAM_QNAME_COLUMN
try:                from hashlib import md5 as md5_new
except ImportError: from md5     import new as md5_new


def usage(s=None):
	message = """
usage: cat sam_file | filtered_sam_to_intervals [options]
  --namesorted             the sam file has been sorted by read names
  --mergemates[=[<min>..<max>] merge the intervals that have the same name,
                           and discard any singletons, multi-chromosomal, or
                           that are outside the expected insert length
  --requiremates[=[<min>..<max>] same as --mergemates, in that mated pairs are
                           required, but intervals are still output separately
  --require:<criterion>    (cumulative) ignore any input lines that don't
                           satisfy the specified criterion;  this is a
                           statement that evaluates to true if the line should
                           be kept, or false if the line should be discarded
  --prohibit:<criterion>   (cumulative) ignore any input lines that satisfy
                           the specified criterion
  --subset:qname=<K>/<N>   conceptually split the input file into <N> groups,
                           and only process the <K>th group; <K> ranges from 1
                           to <N>;  the split utilizes a hash code of the query
                           name, so all records with the same name are in the
                           same subset (this means pairs will remain paired)
  --chromosome[s]=<names>  (cumulative) only output intervals on the specified
                           chromosomes;  <names> is a comma-separated list
                           (default is to report intervals on all chromosomes)
  --origin=one             output intervals as origin-one, closed
  --origin=zero            output intervals as origin-zero, half-open
                           (this is the default)
  --nonames                don't output <read_name> as 4th column
  --samrecords             output entire sam record starting at our 4th column
  --justsamrecords         only output sam records
  --report:<variable>      (cumulative) report additional fields
  --head=<number>          limit the number of input records;  note that
                           records are counted *before* filtering is performed
  --progress=<number>      periodically report how many records we've read
  --progress=output:<number> periodically report how many records we've written

  By default, the output file is a list of <chrom> <start> <end> <read_name>,
  but if --nonames is used, it is just a list of <chrom> <start> <end>

  Criteria for requirements and prohibitions are something like python
  expressions.  Some examples are
    (CIGAR == *)
    (RNEXT == =)
    (FLAGS & 0x9CF == 0x049) or (FLAGS & 0x9CF == 0x089)
    (0x004 in FLAGS)
    (UNCLIP > RLEN*0.10)

  Note that this program has to parse the expressions and convert them to
  actual python statements.  This saves the user from having to deal with the
  shell handling of quote marks-- imagine the effort to get the shell to pass
  the program quote marks in (CIGAR == "*") for example.  Further, the parsing
  allows us to figure out which variables are actually needed in the expression
  context, which allows us to skip the computation of values that aren't
  needed.  The caveat is that this parsing is simplistic.  General expressions
  aren't supported.  It is hoped that those that are supported are useful.

  Names in expresssions come in three flavors.  First, they can be field names
  from the SAM spec:
    QNAME  FLAG   RNAME  POS    MAPQ   CIGAR
    RNEXT  PNEXT  TLEN   SEQ    QUAL
    POS1    POS for first read in a pair   (when --mergemates or
    POS2    POS for second read in a pair   .. --requiremates is used) 
  Second, they can be aliases for tags:
	SCORE   alias for AS tag
	SUBOPT  alias for XS tag
  Third, they can be names we compute from the SAM record:
    RLEN    length of the the read
    BESTBY  difference of alignment score minus suboptimal score
    MAPCLIP number of bases remaining in the alignment after any clipping
    UNCLIP  number of bases clipped from both ends of alignment (not mapped)
    LUNCLIP number of bases clipped from left end of read  (not mapped)
    RUNCLIP number of bases clipped from right end of read (not mapped)
    MINCLIP minimum number of bases clipped from either end of read (not mapped)
    CLIPBRK position of breakpoint suggested by clipped read
    MORIENT orientation of mate within pair (1F, 1R, 2F, or 2R)
    PORIENT orientation of pair (H2H or T2T)
    FLAGS   (same as FLAG)

  There are also a few canned criteria:
    broken mates   stands for  (FLAGS&0x9CF == 0x049) or (FLAGS&0x9CF == 0x089)"""


	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global subsetN,subsetK
	global knownCriteria,computedVariables,tagToVariable,knownVariables,pairVariables
	global variablesNeeded,tagsNeeded,requirements,prohibitions
	global headLimit,reportProgress,progressId,writtenProgress
	global outputWhat,mergeEm,mergeDistanceMin,mergeDistanceMax
	global origin
	global debug

	knownCriteria = { \
		"head to head"               : (None,", use (PORIENT==H2H)"),
		"head-to-head"               : (None,", use (PORIENT==H2H)"),
		"tail to tail"               : (None,", use (PORIENT==T2T)"),
		"tail-to-tail"               : (None,", use (PORIENT==T2T)"),
		"broken mates"               : "(FLAGS & 0x9CF == 0x049) or (FLAGS & 0x9CF == 0x089)",
		"broken mates with unmapped" : "(FLAGS & 0x9CF == 0x049) or (FLAGS & 0x9CF == 0x089) or (FLAGS & 0x9CF == 0x085)"
		}

	tagToVariable = { \
		"AS" : "SCORE",		# primary alignment score
		"XS" : "SUBOPT"		# next best alignment score
		}

	computedVariables = { \
		"RLEN"    : (read_length,       []),
		"MAPCLIP" : (after_clip_length, []),
		"UNCLIP"  : (clip_length,       []),
		"LUNCLIP" : (left_clip_length,  []),
		"RUNCLIP" : (right_clip_length, []),
		"MINCLIP" : (min_clip_length,   []),
		"CLIPBRK" : (clip_breakpoint,   []),
		"BESTBY"  : (score_best_by,     ["SCORE","SUBOPT"]),
		"MORIENT" : (mate_orientation,  []),
		"PORIENT" : (pair_orientation,  [])
		}

	knownVariables = ["QNAME", "FLAG" , "RNAME", "POS"  , "MAPQ" , "CIGAR",
	                  "RNEXT", "PNEXT", "TLEN" , "SEQ"  , "QUAL" , "FLAGS"]
	pairVariables  = ["POS1","POS2"]

	for tag in tagToVariable:
		variable = tagToVariable[tag]
		knownVariables += [variable]

	for variable in computedVariables:
		knownVariables += [variable]

	# parse the command line

	isNameSorted     = False
	mergeEm          = False
	mergeDistanceMin = None
	mergeDistanceMax = None
	mergeButSeparate = False
	requirements     = []
	prohibitions     = []
	subsetN          = None
	subsetK          = None
	chromsOfInterest = None
	origin           = "zero"
	outputWhat       = ["interval","name"]
	headLimit        = None
	reportProgress   = None
	writtenProgress  = None
	progressId       = None
	debug            = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--namesorted"):
			isNameSorted = True
		elif (arg == "--mergemates") \
		  or (arg == "--requiremates"):
			mergeEm          = True
			mergeButSeparate = (arg == "--requiremates")
			mergeDistanceMin = None
			mergeDistanceMax = None
		elif (arg.startswith("--mergemates=")) \
		  or (arg.startswith("--requiremates=")):
			mergeEm          = True
			mergeButSeparate = (arg.startswith("--requiremates="))
			if (".." not in argVal):
				mergeDistanceMin = None
				mergeDistanceMax = int_with_unit(argVal)
			else:
				(argMin,argMax) = argVal.split("..",1)
				if (argMax == ""):
					mergeDistanceMin = int_with_unit(argMin)
					mergeDistanceMax = None
				elif (argMin == ""):
					mergeDistanceMin = None
					mergeDistanceMax = int_with_unit(argMax)
				else:
					mergeDistanceMin = int_with_unit(argMin)
					mergeDistanceMax = int_with_unit(argMax)
		elif (arg.startswith("--require:")):
			argVal = arg.split(":",1)[1]
			requirements += [argVal.strip()]
		elif (arg.startswith("--prohibit:")):
			argVal = arg.split(":",1)[1]
			prohibitions += [argVal.strip()]
		elif (arg.startswith("--subset:qname=")) or (arg.startswith("--subset:name=")):
			assert ("/" in argVal)
			(subsetK,subsetN) = argVal.split("/",1)
			subsetN = int(subsetN)
			subsetK = int(subsetK)
			assert (0 < subsetK <= subsetN)
		elif (arg.startswith("--chromosome=")) or (arg.startswith("--chromosomes=")) \
		  or (arg.startswith("--chrom="))      or (arg.startswith("--chroms=")):
			if (chromsOfInterest == None): chromsOfInterest = []
			chromsOfInterest += argVal.split(",")
		elif (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			assert (origin in ["zero","one"]), "can't understand %s" % arg
		elif (arg == "--nonames"):
			outputWhat = [x for x in outputWhat if (x != "name")]
		elif (arg == "--samrecords"):
			outputWhat =  [x for x in outputWhat if (x not in ["name","sam record"])]
			outputWhat += ["sam record"]
		elif (arg == "--justsamrecords") or (arg == "--justsam"):
			outputWhat = ["sam record"]
		elif (arg.startswith("--report:")):
			argVal = arg.split(":",1)[1]
			for variable in argVal.split(","):
				outputWhat += [variable.strip()]
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=output:")) or  (arg.startswith("--progress=written:")):
			writtenProgress = int_with_unit(argVal.split(":",1)[1])
		elif (arg.startswith("--progress=")):
			if (":" in argVal):
				(progressId,reportProgress) = argVal.split(":",1)
				reportProgress = int_with_unit(reportProgress)
			else:
				reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (mergeEm) and (not mergeButSeparate):
		extras = [x for x in outputWhat if (x not in ["interval","name","sam record"])]
		if (extras != []):
			usage("--report with --mergemates is not implemented yet")

	# preprocess any requirements, changing them into python statements

	variablesNeeded = set()

	incoming = requirements
	requirements = {}
	for criterion in incoming:
		alias = criterion
		if (criterion in knownCriteria):
			alias = knownCriteria[criterion]
			if (alias == None):
				usage("requirement no longer supported: \"%s\"" % criterion)
			elif (type(alias) != str):
				assert (alias[0] == None)
				usage("requirement no longer supported: \"%s\"%s" % (criterion,alias[1]))
		try:
			(expression,variables) = criterion_to_python(alias)
			requirements[criterion] = (expression,variables)
			for name in variables: variablesNeeded.add(name)
		except ValueError:
			usage("uninterpretable requirement: \"%s\"" % criterion)

	incoming = prohibitions
	prohibitions = {}
	for criterion in incoming:
		alias = criterion
		if (criterion in knownCriteria):
			alias = knownCriteria[criterion]
			if (alias == None):
				usage("prohibition no longer supported: \"%s\"" % criterion)
			elif (type(alias) != str):
				assert (alias[0] == None)
				usage("prohibition no longer supported: \"%s\"%s" % (criterion,alias[1]))
		try:
			(expression,variables) = criterion_to_python(alias)
			prohibitions[criterion] = (expression,variables)
			for name in variables: variablesNeeded.add(name)
		except ValueError:
			usage("uninterpretable prohibition: \"%s\"" % criterion)

	for variable in outputWhat:
		if (variable in ["interval","name","sam record"]): continue
		variablesNeeded.add(variable)

	for variable in computedVariables:
		if (variable not in variablesNeeded): continue
		(_,dependencies) = computedVariables[variable]
		for name in dependencies:
			variablesNeeded.add(name)

	variablesNeeded = list(variablesNeeded)
	variablesNeeded.sort()

	if (not mergeEm) and (not mergeButSeparate):
		pairVariablesNeeded = [variable for variable in variablesNeeded
		                                if (variable in pairVariables)]
		assert (pairVariablesNeeded == []), \
		       "pair variables can't be used without --mergemates or --requiremates\n" \
		     + "  pair variables: %s" \
		     % ",".join(pairVariablesNeeded)
	else:
		pairVariablesNeeded = [variable for variable in variablesNeeded
		                                if (variable in pairVariables)]
		assert (pairVariablesNeeded == []), \
		       "sorry, pair variables aren't implemented yet\n" \
		     + "  pair variables: %s" \
		     % ",".join(pairVariablesNeeded)

	tagsNeeded = []
	for tag in tagToVariable:
		variable = tagToVariable[tag]
		if (variable in variablesNeeded):
			tagsNeeded += [tag]
	tagsNeeded.sort()

	if ("evaluation" in debug) and (requirements != {}):
		print >>stderr, "=== requirements ==="
		for criterionStr in requirements:
			(criterion,vars) = requirements[criterionStr]
			if ([var for var in vars if (var in pairVariables)] == []):
				print >>stderr, "  \"%s\" evaluated as \"%s\"" % (criterionStr,criterion)
			else:
				print >>stderr, "  \"%s\" pair-evaluated as \"%s\"" % (criterionStr,criterion)

	if ("evaluation" in debug) and (prohibitions != {}):
		print >>stderr, "=== prohibitions ==="
		for criterionStr in prohibitions:
			(criterion,vars) = prohibitions[criterionStr]
			if ([var for var in vars if (var in pairVariables)] == []):
				print >>stderr, "  \"%s\" evaluated as \"%s\"" % (criterionStr,criterion)
			else:
				print >>stderr, "  \"%s\" pair-evaluated as \"%s\"" % (criterionStr,criterion)

	if ("context" in debug):
		print >>stderr, "=== variables needed ==="
		for name in variablesNeeded:
			print >>stderr, "  %s" % name

	if ("context" in debug) and (tagsNeeded != []):
		print >>stderr, "=== tags needed ==="
		for tag in tagsNeeded:
			print >>stderr, "  %s" % tag

	# process the SAM file

	readIntervals   = None
	readToIntervals = None
	if (mergeEm):
		if (isNameSorted): prevQName = None
		else:              readToIntervals = {}

	for sam in read_sam_simple(stdin):
		(lineNumber,samRecord,context) = sam

		if ("input" in debug):
			print >>stderr, lineNumber,context["QNAME"],context["RNAME"]

		rName = context["RNAME"]
		if (chromsOfInterest != None) and (rName not in chromsOfInterest):
			continue

		cigar = context["CIGAR"]
		if (cigar == "*"):
			rName = rPos = start = end = "*"
		else:
			extent = cigar_to_extent(cigar,lineNumber=lineNumber)
			if (extent == None): continue
			rPos  = context["POS"]
			(left,right) = extent
			start = rPos - left
			end   = rPos + right
			if (start < 0): start = 0
			assert (end >= start), "start > end (%d > %d) at line %d" \
								 % (start,end,lineNumber)
			if (origin == "one"): start += 1

		qName = context["QNAME"]
		if (not mergeEm):
			write_interval(rName,start,end,qName,context,samRecord,rPos,cigar)
			continue

		if (isNameSorted):
			if (prevQName == None):
				readIntervals = []
				prevQName = qName
			elif (qName != prevQName):
				write_merged_interval(prevQName,readIntervals,reportSeparate=mergeButSeparate)
				readIntervals = []
				prevQName = qName
			readIntervals += [(rName,start,end,context,samRecord)]
			continue

		if (qName not in readToIntervals):
			readToIntervals[qName] =  [(rName,start,end,context,samRecord)]
		else:
			readToIntervals[qName] += [(rName,start,end,context,samRecord)]

	# output merged intervals

	if (readIntervals != None):
		write_merged_interval(prevQName,readIntervals,reportSeparate=mergeButSeparate)
	elif (readToIntervals != None):
		for qName in readToIntervals:
			write_merged_interval(qName,readToIntervals[qName],reportSeparate=mergeButSeparate)


def write_merged_interval(qName,intervals,reportSeparate=False):
	if (len(intervals) == 1): return

	if (not reportSeparate):
		rNames = {}
		for (rName,_,_,_,_) in intervals:
			rNames[rName] = True
		if (len(list(rNames)) != 1): return

	intervals.sort()
	(rName1,start1,end1,context1,samRecord1) = intervals[0]
	end = max([e for (_,_,e,_,_) in intervals])
	if (reportSeparate):
		intervals = [(e,s,rName,context,samRecord) for (rName,s,e,context,samRecord) in intervals[1:]]
		intervals.sort()
		(end2,start2,rName2,context2,samRecord2) = intervals[-1]

	if (mergeDistanceMin != None) and (end-start1 < mergeDistanceMin):
		return
	if (mergeDistanceMax != None) and (end-start1 > mergeDistanceMax):
		return

	if (reportSeparate):
		write_interval(rName1,start1,end1,qName,context1,samRecord1)
		write_interval(rName2,start2,end2,qName,context2,samRecord2)
	else:
		write_interval(rName,start1,end,qName)


numberWritten = 0

def write_interval(rName,start,end,qName,
                   context=None,samRecord=None,rPos=None,cigar=None):
	global numberWritten

	line = []

	if ("interval" in outputWhat):
		line += ["%s\t%s\t%s" % (rName,start,end)]

	if ("name" in outputWhat):
		if ("cigar" not in debug):
			line += [qName]
		else:
	 		if (rPos  == None): rPos  = -1
	 		if (cigar == None): cigar = "(none)"
			line += ["%s\t%d\t%s" % (qName,rPos,cigar)]

	for variable in outputWhat:
		if (variable in ["interval","name","sam record"]): continue
		if (context == None) or (variable not in context):
			line += ["(%s)" % variable]
		elif (variable == "FLAGS"):
			line += ["0x%03X" % context[variable]]
		else:
			line += [str(context[variable])]

	if ("sam record" in outputWhat):
 		if (samRecord == None): samRecord = "(sam)"
		line += [samRecord]

	print "\t".join(line)

	numberWritten += 1
	if (writtenProgress != None) and (numberWritten % writtenProgress == 0):
		progressCount = commatize(numberWritten)
		if (outputWhat == ["sam record"]): writingWhat = "sam records"
		else:                              writingWhat = "intervals"
		if (progressId == None):
			print >>stderr, "progress: %s %s written" % (progressCount,writingWhat)
		else:
			print >>stderr, "progress: %s / %s %s written" % (progressId,progressCount,writingWhat)


def read_sam_simple(f):
	lineNumber = recordNumber = 0
	for line in stdin:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("@")):
			if (outputWhat == ["sam record"]): # (nothing but sam is being output)
				if (line.startswith("@SQ")): print line
			continue

		recordNumber += 1
		if (reportProgress != None) and (recordNumber % reportProgress == 0):
			progressCount = commatize(recordNumber)
			if (progressId == None):
				print >>stderr, "progress: %s sam records read" % progressCount
			else:
				print >>stderr, "progress: %s / %s sam records read" % (progressId,progressCount)

		if (headLimit != None) and (recordNumber > headLimit):
			print >>stderr, "limit of %d sam records reached" % headLimit
			break

		fields = line.split()
		numFields = len(fields)
		assert (numFields >= SAM_MIN_COLUMNS), \
		      "not enough columns at line %d (%d, expected %d)" \
		    % (lineNumber,numFields,SAM_MIN_COLUMNS)

		if ("evaluation" in debug):
			print >>stderr
			print >>stderr, "line %d: \"%s\"" % (lineNumber," ".join(fields[:4]))

		# if we are only to process a named-based subset, filter out any reads
		# not in that subset

		if (subsetN != None):
			hashVal = md5_new()
			hashVal.update(fields[SAM_QNAME_COLUMN])
			hashK = 1 + (int(hashVal.hexdigest()[:25],16) % subsetN)
			if (hashK != subsetK): continue

		# create the context in which we will evaluate the criteria

		context = dict(safeDict)
		context["LINENUMBER"] = lineNumber

		for name in samFieldToColumn:
			col = samFieldToColumn[name]
			context[name] = fields[col]

		try:
			context["FLAG" ] = int_or_string(context["FLAG" ])
			context["POS"  ] = int_or_string(context["POS"  ])
			context["MAPQ" ] = int_or_string(context["MAPQ" ])
			context["PNEXT"] = int_or_string(context["PNEXT"])
			context["TLEN" ] = int_or_string(context["TLEN" ])
		except ValueError:
			assert (False), "bad SAM record at line %d\n%s" % (lineNumber,line)

		context["FLAGS"] = context["FLAG"]
		if (type(context["POS"]) != str): context["POS"] -= 1

		for field in fields[SAM_MIN_COLUMNS:]:
			(tag,typeCode,val) = field.split(":",2)
			if (tag not in tagsNeeded):
				continue
			variable = tagToVariable[tag]
			if   (typeCode == "i"): context[variable] = int(val)
			elif (typeCode == "f"): context[variable] = float(val)
			else:                   context[variable] = val

		for tag in tagsNeeded:
			variable = tagToVariable[tag]
			if (variable not in context): context[variable] = None

		for variable in computedVariables:
			if (variable not in variablesNeeded): continue
			(func,_) = computedVariables[variable]
			context[variable] = func(context)

		if ("context" in debug) or ("flags" in debug):
			values = [(name,context[name]) for name in context if (name not in safeDict)]
			if ("flags" in debug):
				values = [(name,context[name]) for name in context if (name == "FLAGS")]
			values.sort()
			print >>stderr, "=== context for line %d ===" % lineNumber
			nameWidth = max([len(name) for (name,_) in values])
			for (name,val) in values:
				if (name in ["SEQ","QUAL"]) and (len(val) > 20):
					val = "\"%s...\"" % val[:20]
				elif (name in ["FLAG","FLAGS"]):
					val = flags_to_string(val)
				elif (type(val) == str): val = "\"%s\"" % val
				print >>stderr, "  %-*s = %s" % (nameWidth,name,val)

		# filter, by evaluating the requirements and prohibitions

		reject = False
		for criterionStr in requirements:
			(criterion,_) = requirements[criterionStr]
			if ("evaluation" in debug):
				print >>stderr, "evaluating requirement \"%s\"" % criterion
			try:
				val = eval(criterion,{"__builtins__":None},context)
			except NameError:
				print >>stderr, "failed to evaluate requirement \"%s\"" % criterion
				raise
			#if ("evaluation" in debug):
			#	print >>stderr, "  (evaluates to %s)" % val
			if (val == False):
				reject = True
				if ("evaluation" in debug): print >>stderr, "  (rejected)"
				break
			if (val != True):
				assert (False), "requirement \"%s\" evaluates to \"%s\" on SAM record at line %d\n%s\nevaluated as \"%s\"" \
				              % (criterionStr,val,lineNumber,line,criterion)
		if (reject): continue

		for criterionStr in prohibitions:
			(criterion,_) = prohibitions[criterionStr]
			if ("evaluation" in debug):
				print >>stderr, "evaluating prohibition \"%s\"" % criterion
			try:
				val = eval(criterion,{"__builtins__":None},context)
			except NameError:
				print >>stderr, "failed to evaluate prohibition \"%s\"" % criterion
				raise
			#if ("evaluation" in debug):
			#	print >>stderr, "  (evaluates to %s)" % val
			if (val == True):
				reject = True
				if ("evaluation" in debug): print >>stderr, "  (rejected)"
				break
			if (val != False):
				assert (False), "prohibition \"%s\" evaluates to \"%s\" on SAM record at line %d\n%s\n(evaluated as \"%s\")" \
				              % (criterionStr,lineNumber,line,criterion)
		if (reject): continue

		# if it made it through all that, keep it

		yield (lineNumber,line,context)


# functions to support "special variables"

def read_length(context):
	return len(context["SEQ"])


def after_clip_length(context):
	cigar = context["CIGAR"]
	cigarInfo = split_cigar(cigar)
	unclip = len(context["SEQ"])
	if (cigarInfo != None):
		(rpt,op) = cigarInfo.operations[0]
		if (op in ["S","H"]): unclip -= rpt
		if (len(cigarInfo.operations) > 1):
			(rpt,op) = cigarInfo.operations[-1]
			if (op in ["S","H"]): unclip -= rpt
	return unclip


def clip_length(context):
	cigar = context["CIGAR"]
	cigarInfo = split_cigar(cigar)
	clip = 0
	if (cigarInfo != None):
		(rpt,op) = cigarInfo.operations[0]
		if (op in ["S","H"]): clip += rpt
		if (len(cigarInfo.operations) > 1):
			(rpt,op) = cigarInfo.operations[-1]
			if (op in ["S","H"]): clip += rpt
	return clip


def left_clip_length(context):
	cigar = context["CIGAR"]
	cigarInfo = split_cigar(cigar)
	if (cigarInfo != None):
		(rpt,op) = cigarInfo.operations[0]
		if (op in ["S","H"]): return rpt
	return 0


def right_clip_length(context):
	cigar = context["CIGAR"]
	cigarInfo = split_cigar(cigar)
	if (cigarInfo != None):
		(rpt,op) = cigarInfo.operations[-1]
		if (op in ["S","H"]): return rpt
	return 0


def min_clip_length(context):
	(lftClip,rgtClip) = cigar_to_clip_lengths(context["CIGAR"])
	return min(lftClip,rgtClip)


def clip_breakpoint(context):
	extent = cigar_to_extent(context["CIGAR"])
	if (extent == None): return "(NO_CLIPBRK)"
	(left,right) = extent
	(lftClip,rgtClip) = cigar_to_clip_lengths(context["CIGAR"])
	if (lftClip == 0) and (rgtClip == 0): return "(NO_CLIPBRK)"

	rPos = context["POS"]
	if (lftClip >= rgtClip):
		start = rPos - left
		if (start < 0): start = 0
		if (origin == "one"): start += 1
		return start
	else:
		end = rPos + right
		return end


def score_best_by(context):
	score  = context["SCORE"]
	subopt = context["SUBOPT"]
	if (score  == None): score  = 0
	if (subopt == None): subopt = 0
	return score - subopt


def mate_orientation(context):
	flags = context["FLAG"]
	if (flags & 0xFF9 == 0x061): return "1F";
	if (flags & 0xFF9 == 0x091): return "2R";
	if (flags & 0xFF9 == 0x051): return "1R";
	if (flags & 0xFF9 == 0x0A1): return "2F";
	return "(MORIENT)"


def pair_orientation(context):
	flags = context["FLAG"]
	tLen  = context["TLEN"]
	if (tLen > 0):
		if (flags & 0xFF9 in [0x061,0x0A1]): return "H2H";
		if (flags & 0xFF9 in [0x091,0x051]): return "T2T";
	elif (tLen < 0):
		if (flags & 0xFF9 in [0x091,0x051]): return "H2H";
		if (flags & 0xFF9 in [0x061,0x0A1]): return "T2T";
	return "(PORIENT)"


# flags_to_string--

nybbleToBits = {  0:"0000",  5:"0101", 10:"1010", 15:"1111",
                  1:"0001",  6:"0110", 11:"1011",
                  2:"0010",  7:"0111", 12:"1100",
                  3:"0011",  8:"1000", 13:"1101",
                  4:"0100",  9:"1001", 14:"1110"}

def flags_to_string(val):
	(v,nybbles) = (val,1)
	while (v > 16):
		(v,nybbles) = (v/16,nybbles+1)
	nybbles = max(nybbles,3)

	s = []
	v = val
	while (nybbles > 0):
		s += [nybbleToBits[v % 16]]
		v /= 16
		nybbles -= 1
	s += ["0x%0*X =" % (nybbles,val)]
	s.reverse()

	return " ".join(s)


# cigar string processing--

def cigar_to_extent(cigar,lineNumber=None):
	cigarInfo = split_cigar(cigar)
	if (cigarInfo == None): return None

	left = 0
	(rpt,op) = cigarInfo.operations[0]
	if (op == "S"): left += rpt
	right = left

	for (rpt,op) in cigarInfo.operations:
		if (op in ["M","X","=","D","N"]):
			right += rpt
		elif (op in ["I","S"]):
			pass
		else:
			if (lineNumber == None):
				assert (False), "unsupported \"%d%s\" in cigar %s" \
				              % (rpt,op,cigar)
			else:
				assert (False), "unsupported \"%d%s\" in cigar %s (line %d)" \
				              % (rpt,op,cigar,lineNumber)

	return (0,right-left)


def cigar_to_clip_lengths(cigar):
	cigarInfo = split_cigar(cigar)
	lftClip = rgtClip = 0
	if (cigarInfo != None):
		(rpt,op) = cigarInfo.operations[0]
		if (op in ["S","H"]): lftClip = rpt
		if (len(cigarInfo.operations) > 1):
			(rpt,op) = cigarInfo.operations[-1]
			if (op in ["S","H"]): rgtClip = rpt
	return (lftClip,rgtClip)


class CigarInfo: pass

def split_cigar(cigar):

	if (cigar == "*"): return None

	# split the cigar into a list of (count,operation)

	operations = []
	rpt = []
	for ch in cigar:
		if (ch.isdigit()):
			rpt += [ch]
		else:
			assert (rpt != []), "bad cigar: \"%s\"" % cigar
			operations += [(int("".join(rpt)),ch)]
			rpt = []
	assert (rpt == []), "bad cigar: \"%s\"" % cigar

	# trim clipping operators from the ends

	startClip = endClip = 0
	if (operations != []):
		(rpt,op) = operations[0]
		if (op == "H"):
			startClip = rpt
			operations = operations[1:]

	if (operations != []):
		(rpt,op) = operations[-1]
		if (op == "H"):
			endClip = rpt
			operations = operations[:-1]

	splitCigar = CigarInfo()
	splitCigar.operations = operations
	splitCigar.startClip  = startClip
	splitCigar.endClip    = endClip
	return splitCigar


# evaluation context stuff--
#	(see reference [2], lybniz2.sourceforge.net/safeeval.html)

def round_int(x): return int(round(x))

def floor_int(x): return int(floor(x))

def ceil_int(x):  return int(ceil(x))

def log2(x):      return log(x) / log(2.0)


safeList = ["e", "exp", "log", "log10", "pi", "pow", "sqrt"]
safeDict = dict([(k,locals().get(k,None)) for k in safeList])
safeDict["abs"]     = abs
safeDict["int"]     = int
safeDict["float"]   = float
safeDict["ceil"]    = ceil_int
safeDict["floor"]   = floor_int
safeDict["round"]   = round_int
safeDict["log2"]    = log2
safeDict["max"]     = max
safeDict["min"]     = min


# criterion to python expression conversion--

eqnVar        = "[A-Z]+"
eqnInt        = "[0-9]+"
eqnFloat      = "[0-9]*\.[0-9]*"
eqnHex        = "0x[0-9A-Fa-f]+"

eqnArith      = "\*|/|\+|\-"
eqnOperator   = "==|!=|>|<|>=|<=|in|not in"
eqnStrConst   = "\*|=|None|1F|1R|2F|2R|H2H|T2T"
eqnObject     = eqnVar + "|" + eqnInt + "|" + eqnHex + "|" + eqnFloat + "|" + eqnStrConst
eqnExpression = eqnVar + "(" + eqnArith + ")" + "(" + eqnInt + "|" + eqnFloat + ")"
eqnFlags      = "FLAGS *& *0x[0-9A-Fa-f]+"

eqnRe = compile("\( *"
              + "(?P<left>" + eqnObject + "|" + eqnExpression + "|" + eqnFlags + ")"
              + " *"
              + "(?P<operator>" + eqnOperator + ")"
              + " *"
              + "(?P<right>" + eqnObject + "|" + eqnExpression + ")"
              + " *\)")
eqnExpressionRe = compile(
                "(?P<left>" + eqnVar + ")"
              + " *"
              + "(?P<operator>" + eqnArith + ")"
              + " *"
              + "(?P<right>" + eqnVar + "|" + eqnInt + "|" + eqnFloat + ")")

eqnFlagsRe    = compile(eqnFlags)


def criterion_to_python(s):
	if ("expressions" in debug):
		print >>stderr, "converting \"%s\" to python" % s

	expression = []
	variables  = set()
	while (True):
		m = eqnRe.search(s)
		if (m == None):
			expression += [s]
			break

		text = m.group(0)
		(prefix,suffix) = s.split(text,1)

		left  = m.group("left")
		op    = m.group("operator")
		right = m.group("right")

		if ("expressions" in debug):
			print >>stderr, "  left  = %s" % left
			print >>stderr, "  op    = %s" % op
			print >>stderr, "  right = %s" % right

		(text,names) = modify_expression(text,left,op,right)
		expression += [prefix,text]
		s = suffix

		for name in names:
			variables.add(name)

	expression = "".join(expression)
	variables  = list(variables)

	if ("expressions" in debug):
		print >>stderr, "  --> \"%s\"" % expression
		print >>stderr, "  --> [%s]" % ",".join(variables)

	return (expression,variables)


def modify_expression(s,left,op,right):
	(preOp,postOp) = s[len(left)+1:-len(right)-1].split(op)

	(left ,leftVarNames)  = sanify(left)
	(right,rightVarNames) = sanify(right)
	varNames = var_names(leftVarNames,rightVarNames)

	if (op == "in") and (right in ["FLAG","FLAGS"]):
		left  = "FLAGS & %s" % left
		op    = "!="
		right = "0"
	elif (op == "not in") and (right in ["FLAG","FLAGS"]):
		left  = "FLAGS & %s" % left
		op    = "=="
		right = "0"

	expression = "".join(["(",left,preOp,op,postOp,right,")"])
	return (expression,varNames)


def sanify(s):
	if (s in knownVariables) or (s in pairVariables):
		return (s,[s])
	if (s == "None"):
		return (s,None)

	if (eqnFlagsRe.match(s) != None):
		return (s,["FLAGS"])

	m = eqnExpressionRe.search(s)
	if (m != None):
		(left ,leftVarNames)  = sanify(m.group("left"))
		(right,rightVarNames) = sanify(m.group("right"))
		if (type(left)  == str) and (left.startswith ("\"")): raise ValueError
		if (type(right) == str) and (right.startswith("\"")): raise ValueError
		varNames = var_names(leftVarNames,rightVarNames)
		return (s,varNames)

	try:
		validate_number(s)
		return (s,None)
	except ValueError:
		return ("\"" + s + "\"",None)


def validate_number(s):
	if (s.upper().startswith("0X")):
		_ = int(s[2:],16)  # just validation, exception raised if it fails
		return

	try:               _ = int(s)
	except ValueError: _ = float(s)


def var_names(names1,names2):
	varNames = []
	if (names1 != None):
		varNames += names1
	if (names2 != None):
		varNames += [name for name in names2 if (name not in varNames)]
	return varNames


# int_or_string--
#	Parse a string as an integer, leaving it as a string if it fails

def int_or_string(s):
	try:               return int(s)
	except ValueError: return s


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
