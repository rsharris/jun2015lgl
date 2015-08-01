#!/usr/bin/env python
"""
Bare-bones SAM reading "class"
"""

import sys
from time import clock

# column indexes for SAM required fields

SAM_QNAME_COLUMN = 0
SAM_FLAG_COLUMN  = 1
SAM_RNAME_COLUMN = 2
SAM_POS_COLUMN   = 3
SAM_MAPQ_COLUMN  = 4
SAM_CIGAR_COLUMN = 5
SAM_MRNM_COLUMN  = 6
SAM_MPOS_COLUMN  = 7
SAM_ISIZE_COLUMN = 8
SAM_SEQ_COLUMN   = 9
SAM_QUAL_COLUMN  = 10

SAM_MIN_COLUMNS  = 11

samFieldToColumn = { "QNAME" : SAM_QNAME_COLUMN,
                     "FLAG"  : SAM_FLAG_COLUMN,
                     "RNAME" : SAM_RNAME_COLUMN,
                     "POS"   : SAM_POS_COLUMN,
                     "MAPQ"  : SAM_MAPQ_COLUMN,
                     "CIGAR" : SAM_CIGAR_COLUMN,
                     "RNEXT" : SAM_MRNM_COLUMN,
                     "PNEXT" : SAM_MPOS_COLUMN,
                     "TLEN"  : SAM_ISIZE_COLUMN,
                     "SEQ"   : SAM_SEQ_COLUMN,
                     "QUAL"  : SAM_QUAL_COLUMN }


# SAM bit-encoded flags

BAM_FPAIRED      =    1	# the read is paired in sequencing, no matter whether it is mapped in a pair
BAM_FPROPER_PAIR =    2	# the read is mapped in a proper pair
BAM_FUNMAP       =    4	# the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
BAM_FMUNMAP      =    8	# the mate is unmapped
BAM_FREVERSE     =   16	# the read is mapped to the reverse strand
BAM_FMREVERSE    =   32	# the mate is mapped to the reverse strand
BAM_FREAD1       =   64	# this is read1
BAM_FREAD2       =  128	# this is read2
BAM_FSECONDARY   =  256	# not primary alignment
BAM_FQCFAIL      =  512	# QC failure
BAM_FDUP         = 1024	# optical or PCR duplicate
BAM_FDUP         = 2048	# suplementary alignment
BAM_NUM_FLAGS    = 12   # total number of flag bits


# read_sam_records--
#	Yields the next sam record, as a methodless class object (see code for
#	object fields).

class SamRecord: pass

def read_sam_records(f,include=None,recordLimit=None,reportProgress=None,progressFmt=None):

	if (include == None):
		include = ["tags"]

	if (reportProgress != None):
		prevTime = clock()
		if (progressFmt == None): progressFmt = "(%.2f) read %d: %s"

	readCount  = 0
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("@")):
			# $$$ we should collect header info and report it, somehow
			continue

		if (recordLimit != None) and (lineNumber > recordLimit):
			print >>sys.stderr, "record limit of %d reached" % recordLimit
			break

		fields = line.split()
		numFields = len(fields)
		assert (numFields >= SAM_MIN_COLUMNS), \
		      "not enough columns at line %d (%d, expected %d)" \
		    % (lineNumber,numFields,SAM_MIN_COLUMNS)

		samrec = SamRecord()

		try:
			samrec.flag = int(fields[SAM_FLAG_COLUMN])
			if (samrec.flag < 0): raise ValueError
		except ValueError:
			assert (False), "bad SAM flag at line %d\n%s" % (lineNumber,line)

		samrec.qName = fields[SAM_QNAME_COLUMN]
		samrec.rName = fields[SAM_RNAME_COLUMN]
		samrec.rPos  = fields[SAM_POS_COLUMN]
		samrec.mapQ  = fields[SAM_MAPQ_COLUMN]
		samrec.cigar = fields[SAM_CIGAR_COLUMN]
		samrec.mrnm  = fields[SAM_MRNM_COLUMN]
		samrec.mPos  = fields[SAM_MPOS_COLUMN]
		samrec.iSize = fields[SAM_ISIZE_COLUMN]
		samrec.seq   = fields[SAM_SEQ_COLUMN]
		samrec.qual  = fields[SAM_QUAL_COLUMN]

		if ("line" in include):
			samrec.line = line
		if ("line number" in include):
			samrec.lineNumber = lineNumber

		if ("tags" in include) and (numFields > SAM_MIN_COLUMNS):
			tags = {}
			for tag in fields[SAM_MIN_COLUMNS:]:
				(tag,val) = tag.split(":",1)
				assert (tag not in tags)
				tags[tag] = val
			samrec.tags = tags

		readCount += 1
		if (reportProgress != None) and (readCount % reportProgress == 0):
			currTime = clock()
			print >>sys.stderr, progressFmt % (currTime-prevTime,readCount,samrec.qName)
			prevTime = currTime

		yield samrec


def sam_flags_to_binary_string(flags):
	s = []
	f = flags
	c = 0
	while (flags > 0) or (c < BAM_NUM_FLAGS):
		if (c != 0) and (c % 4 == 0): s += ["_"]
		if (flags & 1 == 1): s += ["1"]
		else:                s += ["0"]
		flags = flags >> 1
		c += 1
	s.reverse()
	return "".join(s)
