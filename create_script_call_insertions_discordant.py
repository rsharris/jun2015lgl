#!/usr/bin/env python
"""
Create a cluster job file to call insertions and create a corresponding
discriminator track.
"""

from sys import argv,stdin,stdout,stderr,exit


def usage(s=None):
	message = """

usage: create_script_call_insertions [options] > call_insertions.sh
  <sub>_<samp>_<type>        (required) run descriptor; for example, CS_NORM_PE
                             means subject "CS", sample "NORM", and type "PE";
                             one PE and one MP are required; other filenames
                             can use "{perun}" and "{mprun}" to refer to these
                             strings, or just {run}
  --control=<filename>       read control values from a file (see list below)
  --base=<path>              path prefix; other filenames can use "{base}" to
                             refer to this path
  --chromosomes=<filename>   read chromosome names and lengths from a file
                             (default is {base}/data/hg19.chrom_lengths)
  --input=<filename>         (cumulative) track files to process; by default
                             this is these three tracks:
                               {base}/tracks/{mprun}.insert_length_sparse_or_normal_inserts_sparse
                               {base}/tracks/{mprun}.short_or_discordant
                               {base}/tracks/{perun}.clipped_breakpoints.high
  --pipeline=<text>          description of pipeline (for display with track)
  --track=<filename>         track file to create
                             (default is {base}/tracks/{run}.called_insertions)
  --tempinput=<filename>     temporary file to hold input track file(s);  two
                             files will be needed if any of the input tracks
                             are gzipped
                             (default is {base}/tracks/{run}.{num}.called_insertions.scratch)
  --temp=<filename>          temporary file to hold track file and intermediate
                             track file(s)
                             (default is {base}/tracks/{run}.called_insertions.temp)
  --gzip                     compress track file
  --bigwig[=<filename>]      create bigwig file in addition to track file
  --bigwigchroms=<filename>  chromosomes file for bedGraphToBigWig
                             (default is {base}/temp/ucsc.hg19.chrom_lengths)
  --bigwigurl=<url>          url for the bigwig file; this can use {bigwig}
                             for the bigwig filename
  --bigwiglink=<filename>    path at which to create a symbolic link to the
                             bigwig and info files; this can use {bigwig}
                             for the bigwig filename
  --bigwigposition=<interval> intitial UCSC browser interval for bigwig track
  --catalog[=<filename>]     create catalog file in addition to track file
  --cataloglink=<filename>   path at which to create a symbolic link to the
                             catalog file; this can use {catalog} for the
                             catalog filename
  --cataloggenome=<species>  genome assembly catalog will refer to
  --catalogtitle=<species>   title for catalog file
  --catalogspec=<nick>:<track>  (cumulative) spec for catalog comments column
  --initialize=<text>        (cumulative) shell command to add to job beginning
                             "shebang:bash" is mapped "#!/usr/bin/env bash"
                             other commands are copied "as is"

values read from control file:
  call_insertions.closure
  call_insertions.proximity"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global basePath,runId,peRunName,mpRunName
	global debug

	bashShebang = "#!/usr/bin/env bash"

	# parse args

	peRunName            = None
	mpRunName            = None
	runId                = None
	controlFilename      = None
	basePath             = None
	chromsFilename       = None
	pipelineFilenames    = None
	pipelineText         = None
	trackName            = None
	tempInputFilenames   = None
	tempFilename         = None
	gzipOutput           = False
	bigWigFilename       = None
	bigWigChromsFilename = None
	bigWigUrl            = None
	bigWigLink           = None
	bigWigPosition       = None
	catalogFilename      = None
	catalogLink          = None
	catalogGenome        = None
	catalogTitle         = None
	catalogSpecs         = None
	bashInitializers     = ["set -eu"]
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1].strip()

		if (arg.startswith("--control=")):
			controlFilename = argVal
		elif (arg.startswith("--base=")) or (arg.startswith("--basepath=")) or (arg.startswith("--path=")):
			basePath = argVal
		elif (arg.startswith("--chromosomes=")) or (arg.startswith("--chroms=")):
			chromsFilename = argVal
		elif (arg.startswith("--input=")):
			if (pipelineFilenames == None): pipelineFilenames =  [argVal]
			else:                           pipelineFilenames += [argVal]
		elif (arg.startswith("--pipeline=")):
			pipelineText = argVal
		elif (arg.startswith("--track=")):
			trackName = argVal
		elif (arg.startswith("--tempinput=")):
			if (tempInputFilenames == None): tempInputFilenames =  [argVal]
			else:                            tempInputFilenames += [argVal]
		elif (arg.startswith("--temp=")):
			tempFilename = argVal
		elif (arg == "--gzip"):
			gzipOutput = True
		elif (arg == "--bigwig"):
			bigWigFilename = "{track}.bw"
		elif (arg.startswith("--bigwig=")):
			bigWigFilename = argVal
		elif (arg.startswith("--bigwigchromosomes=")) or (arg.startswith("--bigwigchroms=")):
			bigWigChromsFilename = argVal
		elif (arg.startswith("--bigwigurl=")) or (arg.startswith("--url=")):
			bigWigUrl = argVal
		elif (arg.startswith("--bigwiglink=")) or (arg.startswith("--link=")):
			bigWigLink = argVal
		elif (arg.startswith("--bigwigposition=")) or (arg.startswith("--bigwigpos=")):
			bigWigPosition = argVal
		elif (arg == "--catalog"):
			catalogFilename = "{track}.catalog.html"
		elif (arg.startswith("--catalog=")):
			catalogFilename = argVal
		elif (arg.startswith("--cataloglink=")) or (arg.startswith("--link=")):
			catalogLink = argVal
		elif (arg.startswith("--cataloggenome=")):
			catalogGenome = argVal
		elif (arg.startswith("--catalogtitle=")):
			catalogTitle = argVal
		elif (arg.startswith("--catalogspec=")):
			if (":" not in argVal):
				usage("\"%s\" lacks the required \":\"" % arg)
			if (catalogSpecs == None): catalogSpecs = []
			catalogSpecs += [tuple(argVal.split(":",1))]
		elif (arg.startswith("--initialize=")) or (arg.startswith("--init=")):
			if (argVal == "shebang:bash"):
				argVal = bashShebang
			if (argVal == "set -eu"):
				bashInitializers = [x for x in bashInitializers if (x != "set -eu")]
			bashInitializers += [argVal]
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (peRunName == None) or (mpRunName == None):
			fields = arg.split(":",2)
			if (len(fields) != 3):
				fields = arg.split("_")
			if (len(fields) < 3) or (fields[-1] not in ["PE","MP"]):
				usage("\"%s\" is not a valid run descriptor" % arg)
			if (fields[-1] == "PE"):
				if (peRunName != None): usage("can't define more than one PE run: %s" % arg)
				if (runId == None): runId = fields[0]
				elif (runId != fields[0]): usage("run identifiers must match (\"%s\" and \"%s\")" % (runId,fields[0]))
				peRunName = "_".join(fields)
			elif (fields[-1] == "MP"):
				if (mpRunName != None): usage("can't define more than one MP run: %s" % arg)
				if (runId == None): runId = fields[0]
				elif (runId != fields[0]): usage("run identifiers must match (\"%s\" and \"%s\")" % (runId,fields[0]))
				mpRunName = "_".join(fields)
			else:
				usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (peRunName == None):
		usage("you have to give me a PE run descriptor")

	if (mpRunName == None):
		usage("you have to give me a MP run descriptor")

	if (controlFilename == None):
		usage("you have to give me a control filename")

	if (chromsFilename == None):
		chromsFilename = "{base}/data/hg19.chrom_lengths"

	if (pipelineFilenames == None):
		#pipelineFilenames = ["{base}/tracks/{mprun}.insert_length.sparse",
		#                     "{base}/tracks/{mprun}.normal_inserts.depth.sparse",
		#                     "{base}/tracks/{mprun}.short_inserts.depth.dense",
		#                     "{base}/tracks/{perun}.clipped_breakpoints.high"]
		#if (pipelineText == None):
		#	pipelineText = "mp_insert_lengths- mp_normal_inserts- short_inserts+ pe_clipped_breakpoints.high+"
		pipelineFilenames = ["{base}/tracks/{mprun}.insert_length_sparse_or_normal_inserts_sparse",
		                     "{base}/tracks/{mprun}.short_or_discordant",
		                     "{base}/tracks/{perun}.clipped_breakpoints.high"]
		if (pipelineText == None):
			pipelineText = "mp_insert_lengths_OR_mp_normal_inserts- mp_short_OR_discordant+ pe_clipped_breakpoints.high+"
	elif (len(pipelineFilenames) < 2):
		usage("you have to give me at least two input tracks")

	if (pipelineText == None):
		pipelineText = "called insertions"

	if (trackName == None):
		trackName = "{base}/tracks/{run}.called_insertions"

	tempInputFilesNeeded = 1
	for (ix,inputFilename) in enumerate(pipelineFilenames[1:]):
		inputFilename = do_filename_substitutition(inputFilename)
		if (inputFilename.endswith(".gz")) or (inputFilename.endswith(".gzip")):
			tempInputFilesNeeded = 2
			break

	if (tempInputFilenames == None):
		if (tempInputFilesNeeded == 1):
			tempInputFilenames = [trackName + ".scratch"]
		else:
			tempInputFilenames = [trackName + ".{num}.scratch"]
	if (len(tempInputFilenames) == 1) and (tempInputFilesNeeded > 1):
		tempInputFilename = tempInputFilenames[0]
		if ("{num}" in tempInputFilename):
			tempInputFilenames = []
			for num in xrange(1,tempInputFilesNeeded+1):
				tempInputFilenames += [tempInputFilename.replace("{num}",str(num))]
	if (len(tempInputFilenames) < tempInputFilesNeeded):
		usage("you have to give %d temporary input files" % tempInputFilesNeeded)
	elif (len(tempInputFilenames) > tempInputFilesNeeded):
		tempInputFilenames = tempInputFilenames[:tempInputFilesNeeded]

	if (tempFilename == None):
		tempFilename = trackName + ".temp"

	if (bigWigFilename != None):
		if (bigWigChromsFilename == None):
			bigWigChromsFilename = "{base}/temp/ucsc.hg19.chrom_lengths"
		if (bigWigUrl == None):
			usage("you have to give me a url for the bigwig file")

	if (catalogFilename != None):
		if (catalogGenome == None):
			usage("you have to give me a genome/assembly name for the catalog file")

	trackId = "%s.called_insertions" % runId

	##########
	# perform filename substitution
	##########

	if   (basePath == None):       basePath = "."
	elif (basePath.endswith("/")): basePath = basePath[:-1]

	controlFilename = do_filename_substitutition(controlFilename)
	chromsFilename  = do_filename_substitutition(chromsFilename)

	# input track names

	for (ix,inputFilename) in enumerate(pipelineFilenames):
		inputFilename = do_filename_substitutition(inputFilename)
		if   ((not inputFilename.endswith(".gz")) 
		  and (not inputFilename.endswith(".gzip"))
		  and (not inputFilename.endswith(".dat"))):
			inputFilename += ".dat"
		pipelineFilenames[ix] = inputFilename

	# track name

	trackName = do_filename_substitutition(trackName)

	trackFilename = trackName
	if (gzipOutput):
		if (not trackFilename.endswith(".gz")): trackFilename += ".gz"
	else:
		if (not trackFilename.endswith(".dat")): trackFilename += ".dat"

	if (tempInputFilenames != None):
		for (ix,tempInputFilename) in enumerate(tempInputFilenames):
			tempInputFilenames[ix] = do_filename_substitutition(tempInputFilename)

	if (tempFilename != None):
		tempFilename = do_filename_substitutition(tempFilename)

	# big wig name

	if (bigWigFilename != None):
		bigWigFilename = do_filename_substitutition(bigWigFilename)
		if ("{track}" in bigWigFilename):
			trackTemp = trackName
			if   (trackTemp.endswith(".gz")):  trackTemp = trackTemp[:-3]
			elif (trackTemp.endswith(".dat")): trackTemp = trackTemp[:-4]
			bigWigFilename = bigWigFilename.replace("{track}",trackTemp)
			if (bigWigFilename.endswith(".bw")): infoFilename = bigWigFilename[:-3] + ".info"
			else:                                infoFilename = bigWigFilename      + ".info"

	if (bigWigChromsFilename != None):
		bigWigChromsFilename = do_filename_substitutition(bigWigChromsFilename)

	if (bigWigUrl != None):
		bigWigTemp = bigWigFilename
		slashIx = bigWigTemp.rfind("/")
		if (slashIx >= 0): bigWigTemp = bigWigTemp[slashIx+1:]
		bigWigUrl = bigWigUrl.replace("{bigwig}",bigWigTemp)

	if (bigWigLink != None):
		bigWigSave = bigWigLink
		bigWigTemp = bigWigFilename
		slashIx = bigWigTemp.rfind("/")
		if (slashIx >= 0): bigWigTemp = bigWigTemp[slashIx+1:]
		bigWigLink = bigWigLink.replace("{bigwig}",bigWigTemp)

		infoTemp = infoFilename
		slashIx = infoTemp.rfind("/")
		if (slashIx >= 0): infoTemp = infoTemp[slashIx+1:]
		infoLink = bigWigSave.replace("{bigwig}",infoTemp)

	# catalog name

	if (catalogFilename != None):
		catalogFilename = do_filename_substitutition(catalogFilename)
		if ("{track}" in catalogFilename):
			trackTemp = trackName
			if   (trackTemp.endswith(".gz")):  trackTemp = trackTemp[:-3]
			elif (trackTemp.endswith(".dat")): trackTemp = trackTemp[:-4]
			catalogFilename = catalogFilename.replace("{track}",trackTemp)

	if (catalogLink != None):
		catalogSave = catalogLink
		catalogTemp = catalogFilename
		slashIx = catalogTemp.rfind("/")
		if (slashIx >= 0): catalogTemp = catalogTemp[slashIx+1:]
		catalogLink = catalogLink.replace("{catalog}",catalogTemp)

	if (catalogSpecs != None):
		for (ix,(nick,name)) in enumerate(catalogSpecs):
			name = do_filename_substitutition(name)
			catalogSpecs[ix] = (nick,name)

	##########
	# get values from the control file
	##########

	closureLength = None
	proximity     = None

	f = file(controlFilename,"rt")

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 3), \
		       "not enough fields at control file line %d (%d, expected at least 3)" \
		     % (lineNumber,len(fields))
		assert (fields[1] == "="), \
		       "can't understand control file line %d:\n%s" \
		     % (lineNumber,line)

		(name,_,val) = fields[:3]
		if (name == "call_insertions.closure"):   closureLength = val
		if (name == "call_insertions.proximity"): proximity = val

	f.close()

	if (closureLength == None): assert(False), "control file lacks closure"
	if (proximity     == None): assert(False), "control file lacks proximity"

	if (closureLength == "NA"): closureLength = None

	##########
	# create the job's shell script
	##########

	# write bash intitializers

	if (bashInitializers != None):
		for (ix,bashInitializer) in enumerate(bashInitializers):
			if (bashInitializer != bashShebang): continue
			print bashInitializer
			bashInitializers[ix] = None
		for (ix,bashInitializer) in enumerate(bashInitializers):
			if (bashInitializer != "set -eu"): continue
			print bashInitializer
			bashInitializers[ix] = None
		for bashInitializer in bashInitializers:
			if (bashInitializer == None): continue
			print do_filename_substitutition(bashInitializer)
		print

	print "today=`today {mmm}/{d}/{yyyy}`"

	# write commands describing the files the script will create

	if (tempInputFilenames != None):
		for tempInputFilename in tempInputFilenames:
			print "echo \"will use %s for temporary input files\"" % tempInputFilename

	if (tempFilename != None):
		print "echo \"will write temporary files to %s\"" % tempFilename

	print "echo \"will write track file to      %s\"" % trackFilename

	if (bigWigFilename != None):
		print "echo \"will write bigwig file to     %s\"" % bigWigFilename

	if (catalogFilename != None):
		print "echo \"will write catalog file to    %s\"" % catalogFilename

	# write command(s) to create track file
	#
	# the process is as follows:
	#
	#	previous stage output is track 1
	#
	#	for track N in 1..  (note that this skips the first input track)
	# 		if track N is not zipped
	#			stage input is track N
	#		else
	#			unzip track N to temp2
	#			stage input is temp2
	#		cat previous stage output | proximal stage input > temp3
	#		mv temp3 to temp1
	#		if track N was zipped
	#			destroy temp2
	#		previous stage output is temp1

	print
	print "echo \"=== creating track %s ===\"" % trackId

	tempInputFilename1 = tempInputFilenames[0]
	if (len(tempInputFilenames) > 1): tempInputFilename2 = tempInputFilenames[1]
	else:                             tempInputFilename2 = None

	for (stageIx,inputFilename) in enumerate(pipelineFilenames):
		if (stageIx == 0):
			prevStageOutput = inputFilename
			continue

		isFinalStage = (stageIx == len(pipelineFilenames)-1)

		# do we need to unzip the next input track, or not?

		if   ((not inputFilename.endswith(".gz"))
		  and (not inputFilename.endswith(".gzip"))):
			trackSourceFilename = inputFilename
		else:
			commands =  []
			command  =  ["time gzip -dc %s" % inputFilename]
			commands += [command]
			command  =  ["> %s" % tempInputFilename2]
			commands += [command]

			print
			print commands_to_pipeline(commands)

			trackSourceFilename = tempInputFilename2

		# what's the output at this stage, intermediate or final?

		if (isFinalStage) and (gzipOutput):
			trackDestFilename  = trackFilename
			trackCountFilename = tempFilename
		elif (isFinalStage):
			trackDestFilename  = trackFilename
			trackCountFilename = trackFilename
		else:
			trackDestFilename  = tempFilename
			trackCountFilename = tempFilename

		# do we need to unzip the previous output stage, or not?  (the only
		# time we will unzip here is if the first track input is compressed)

		commands = []
		if   ((not prevStageOutput.endswith(".gz"))
		  and (not prevStageOutput.endswith(".gzip"))):
			command = ["time cat %s" % prevStageOutput]
			commands += [command]
		else:
			command = ["time gzip -dc %s" % prevStageOutput]
			commands += [command]

		command  =  ["proximal_feature_intervals --positive"]
		command  += [trackSourceFilename]
		command  += ["--proximity=%s" % proximity]
		commands += [command]

		command  =  ["awk '{ print $1,$2,$3 }'"]
		commands += [command]

		command  =  ["keep_first"]
		commands += [command]

		if (isFinalStage):
			if (closureLength != None):
				command  =  ["close_intervals %s" % closureLength]
				commands += [command]
			command  =  ["fill_genomic_interval_gaps --chroms=%s" % chromsFilename]
			commands += [command]

		if (isFinalStage) and (gzipOutput):
			if (tempFilename != None):
				command  =  ["tee %s" % tempFilename]
				commands += [command]
			command  =  ["gzip"]
			commands += [command]

		command  =  ["> %s" % trackDestFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		if (not isFinalStage):
			commands =  []
			command  =  ["echo \"stage %d: `wc -l %s | awk '{print $1}'` intervals\"" \
			           % (stageIx+1,trackCountFilename)]
			commands += [command]
			print commands_to_pipeline(commands)
		else:
			commands =  []
			command  =  ["echo \"stage %d: `cat %s | awk '{if ($4>0) n++} END {print n}'` intervals\"" \
			           % (stageIx+1,trackCountFilename)]
			commands += [command]
			print commands_to_pipeline(commands)

		# if this isn't the final stage, move the temporary output file to the
		# temporary input file

		if (not isFinalStage):
			commands =  []
			command  =  ["mv %s %s" % (tempFilename,tempInputFilename1)]
			commands += [command]
			print commands_to_pipeline(commands)
			prevStageOutput = tempInputFilename1

	# if we created a temporary input file, get rid of it

	if (tempInputFilename2 != None):
		commands =  []
		command  =  ["rm %s" % (tempInputFilename2)]
		commands += [command]

		print
		print commands_to_pipeline(commands)

	# write command(s) to convert track file to bigwig

	if (bigWigFilename != None):
		print
		print "echo \"=== converting track %s to bigwig ===\"" % trackId

		if (gzipOutput): trackInput = tempFilename
		else:            trackInput = trackFilename

		commands =  []
		command  =  ["time bedGraphToBigWig"]
		command  += [trackInput]
		command  += [bigWigChromsFilename]
		command  += [bigWigFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		if (gzipOutput):
			commands =  []
			command  =  ["rm %s" % tempFilename]
			commands += [command]
			print
			print commands_to_pipeline(commands)

		commands =  []
		command  =  ["make_bigwig_info"]
		command  += ["--url=%s" % bigWigUrl]
		command  += ["--name=\"%s called insertions\"" % runId]
		command  += ["--desc=\"%s called insertions (${today}) (%s)\"" \
				   % (runId,pipelineText)]
		command  += ["--autoscale=\"on\""]
		command  += ["--alwayszero=\"on\""]
		command  += ["--maxheight=\"10:10:10\""]
		command  += ["--color=250,30,100"]
		if (bigWigPosition != None): command  += ["--pos=\"%s\"" % bigWigPosition]
		command  += ["> %s" % infoFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		if (bigWigLink != None):
			print
			print "rm -f %s" % infoLink
			print "ln -s %s %s" % (infoFilename,infoLink)
			print "rm -f %s" % bigWigLink
			print "ln -s %s %s" % (bigWigFilename,bigWigLink)

			infoUrl = bigWigUrl
			slashIx = infoUrl.rfind("/")
			if (slashIx >= 0): infoUrl = infoUrl[:slashIx]
			infoTemp = infoFilename
			slashIx = infoTemp.rfind("/")
			if (slashIx >= 0): infoTemp = infoTemp[slashIx+1:]
			infoUrl = infoUrl + "/" + infoTemp
			print >>stderr, infoUrl
			print
			print "echo \"track URL is %s\"" % (infoUrl)

	# write command(s) to create a catalog of the track file

	if (catalogFilename != None):
		print
		print "echo \"=== converting track %s to catalog ===\"" % trackId

		if (catalogSpecs == None):
			commands = []
			if (gzipOutput): commands += [["gzip -dc %s" % trackFilename]]
			else:            commands += [["cat %s" % trackFilename]]
			commands += [["awk '{ if ($4>0) print $0 }'"]]

			command  =  ["intervals_to_ucsc_catalog"]
			command  += ["--show:numbers"]
			command  += ["--genome=%s" % catalogGenome]
			command  += ["--split=20,80"]
			command  += ["--center=20K"]
			command  += ["--center=150K"]
			command  += ["--catalogonly"]
			if (catalogTitle != None): command += ["--title=\"%s (${today})\"" % catalogTitle]
			command  += [catalogFilename]
			commands += [command]

			print
			print commands_to_pipeline(commands)

		else: # (catalogSpecs != None):
			print
			specs = ["%s:%s" % (nick,name) for (nick,name) in catalogSpecs]
			print "specs=\"%s\"" % "\n       ".join(specs)

			commands = []
			commands += [["echo ${specs}"]]
			commands += [["tr \" :\" \"\\n \""]]
			commands += [["while read nick trackname ; do"]]
			print commands_to_pipeline(commands)

			commands = []
			if (gzipOutput): commands += [[" gzip -dc %s" % trackFilename]]
			else:            commands += [[" cat %s" % trackFilename]]

			command  =  [" proximal_feature_intervals --positive"]
			command  += ["${trackname}"]
			command  += ["--proximity=%s" % proximity]
			commands += [command]

			commands += [[" awk '{ print $1\"~\"$2\"~\"$3,nick,1 }' nick=${nick}"]]
			commands += [[" keep_first"]]
			print commands_to_pipeline(commands)

			commands =  []
			commands += [[" done"]]
			commands += [["collect_tags --separator=~"]]
			commands += [["awk '{ print $1,\"#\",$2 }'"]]
			commands += [["sed \"s/~/ /g\""]]
			commands += [["encodachrom | env LC_ALL=C sort -k 1,1n -k 2,2n | decodachrom"]]

			command  =  ["intervals_to_ucsc_catalog"]
			command  += ["--show:numbers"]
			command  += ["--show:comments"]
			command  += ["--genome=%s" % catalogGenome]
			command  += ["--split=20,80"]
			command  += ["--center=20K"]
			command  += ["--center=150K"]
			command  += ["--catalogonly"]
			if (catalogTitle != None): command += ["--title=\"%s (${today})\"" % catalogTitle]
			command  += [catalogFilename]
			commands += [command]
			print commands_to_pipeline(commands)

		if (catalogLink != None):
			print
			print "rm -f %s" % catalogLink
			print "ln -s %s %s" % (catalogFilename,catalogLink)


def commands_to_pipeline(commands):
	pipeline = []
	for (cmdNum,cmd) in enumerate(commands):
		if (cmdNum == 0): prefix = ""
		else:             prefix = "  | "
		indent = ""

		if (cmd[0].startswith(" ")):
			indent = "      "
			prefix = indent + prefix
			cmd[0] = cmd[0][1:]

		if (cmd[0].startswith(">")):
			assert (cmdNum != 0)
			assert (len(cmd) == 1)
			prefix = "  "

		pipeline += [prefix + cmd[0]] 
		for line in cmd[1:]:
			pipeline += [indent + "      " + line]

	return " \\\n".join(pipeline)


def do_filename_substitutition(s):
	if ("{base}" in s):
		assert (basePath != None)
		s = s.replace("{base}",basePath)
	if ("{run}" in s):
		assert (runId != None)
		s = s.replace("{run}",runId)
	if ("{perun}" in s) or ("{runpe}" in s):
		assert (runId != None)
		s = s.replace("{perun}",peRunName)
	if ("{mprun}" in s) or ("{runmp}" in s):
		assert (runId != None)
		s = s.replace("{mprun}",mpRunName)
	return s


if __name__ == "__main__": main()

