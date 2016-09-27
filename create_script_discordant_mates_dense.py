#!/usr/bin/env python
"""
Create a cluster job file to create a discordant-mates density discriminator
track.
"""

from sys  import argv,stdin,stdout,stderr,exit


def usage(s=None):
	message = """

usage: create_script_discordant_mates_dense [options] > discordant_mates_dense.sh
  <sub>_<samp>_<type>        (required) run descriptor; for example, CS_NORM_PE
                             means subject "CS", sample "NORM", and type "PE";
                             other filenames can use "{run}" to refer to this
                             string
  --control=<filename>       read control values from a file (see list below)
  --base=<path>              path prefix; other filenames can use "{base}" to
                             refer to this path
  --chromosomes=<filename>   read chromosome names and lengths from a file
                             (default is {base}/data/hg19.chrom_lengths)
  --input=<filename>         (required) track file to process
  --track=<filename>         (required) track file to create
                             (default is {base}/tracks/{run}.discordant_mates.dense)
  --tempinput=<filename>     temporary file to hold input track file, if
                             needed; only needed if the input track was gzipped
                             (default is {base}/tracks/{run}.discordant_mates.dense.scratch)
  --temp=<filename>          temporary file to hold track file, if needed; only
                             needed if --gzip and --bigwig are both used
                             (default is {base}/tracks/{run}.discordant_mates.dense.temp)
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
  --initialize=<text>        (cumulative) shell command to add to job beginning
                             "shebang:bash" is mapped "#!/usr/bin/env bash"
                             other commands are copied "as is"

values read from control file:
  discordant_mates_dense.minLength.{run}
  discordant_mates_dense.density
  discordant_mates_dense.samplingStep"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global basePath,runName
	global debug

	bashShebang = "#!/usr/bin/env bash"

	# parse args

	runName              = None
	controlFilename      = None
	basePath             = None
	inputFilename        = None
	chromsFilename       = None
	trackName            = None
	scratchFilename      = None
	tempFilename         = None
	gzipOutput           = False
	bigWigFilename       = None
	bigWigChromsFilename = None
	bigWigUrl            = None
	bigWigLink           = None
	bigWigPosition       = None
	bashInitializers     = ["set -eu"]
	debug                = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1].strip()

		if (arg.startswith("--control=")):
			controlFilename = argVal
		elif (arg.startswith("--base=")) or (arg.startswith("--basepath=")) or (arg.startswith("--path=")):
			basePath = argVal
		elif (arg.startswith("--input=")):
			inputFilename = argVal
		elif (arg.startswith("--chromosomes=")) or (arg.startswith("--chroms=")):
			chromsFilename = argVal
		elif (arg.startswith("--track=")):
			trackName = argVal
		elif (arg.startswith("--tempinput=")):
			tempInputFilename = argVal
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
		elif (runName == None):
			fields = arg.split(":",2)
			if (len(fields) != 3):
				fields = arg.split("_")
			if (len(fields) < 3) or (fields[-1] not in ["PE","MP"]):
				usage("\"%s\" is not a valid run descriptor" % arg)
			runName = "_".join(fields)
		else:
			usage("unrecognized option: %s" % arg)

	if (runName == None):
		usage("you have to give me a run descriptor")

	if (controlFilename == None):
		usage("you have to give me a control filename")

	if (inputFilename == None):
		usage("you have to give me an input track filename")

	if (chromsFilename == None):
		chromsFilename = "{base}/data/hg19.chrom_lengths"

	if (trackName == None):
		trackName = "{base}/tracks/{run}.discordant_mates.dense"

	if (scratchFilename == None):
		scratchFilename = trackName + ".scratch"

	if (tempFilename == None) and (bigWigFilename != None) and (gzipOutput):
		tempFilename = trackName + ".temp"

	if (bigWigFilename != None):
		if (bigWigChromsFilename == None):
			bigWigChromsFilename = "{base}/temp/ucsc.hg19.chrom_lengths"
		if (bigWigUrl == None):
			usage("you have to give me a url for the bigwig file")

	trackId = "%s.discordant_mates.dense" % runName

	##########
	# perform filename substitution
	##########

	if   (basePath == None):       basePath = "."
	elif (basePath.endswith("/")): basePath = basePath[:-1]

	controlFilename = do_filename_substitutition(controlFilename)
	chromsFilename  = do_filename_substitutition(chromsFilename)

	# input track name

	inputFilename = do_filename_substitutition(inputFilename)

	gzipInput = False
	if   (inputFilename.endswith(".gz")):       gzipInput = True
	elif (inputFilename.endswith(".gzip")):     gzipInput = True
	elif (inputFilename.endswith(".bedgraph")): pass
	elif (not inputFilename.endswith(".dat")):  inputFilename += ".dat"

	# track name

	trackName = do_filename_substitutition(trackName)

	trackFilename = trackName
	if (gzipOutput):
		if (not trackFilename.endswith(".gz")): trackFilename += ".gz"
	else:
		if (not trackFilename.endswith(".dat")): trackFilename += ".dat"

	if (scratchFilename != None):
		scratchFilename = do_filename_substitutition(scratchFilename)

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

	##########
	# get values from the control file
	##########

	minLength        = None
	densityThreshold = None
	samplingStep     = None

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
		if (name == "discordant_mates_dense.minLength." + runName): minLength        = int(val)
		if (name == "discordant_mates_dense.density"):              densityThreshold = val
		if (name == "discordant_mates_dense.samplingStep"):         samplingStep     = int(val)

	f.close()

	if (minLength        == None): assert(False), "control file lacks minLength"
	if (samplingStep     == None): assert(False), "control file lacks samplingStep"
	if (densityThreshold == None): assert(False), "control file lacks density"

	if ("." in densityThreshold):
		while (densityThreshold.endswith("0")):
			densityThreshold = densityThreshold[:-1]
		if (densityThreshold.endswith(".")):
			densityThreshold = densityThreshold[:-1]

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

	if (scratchFilename != None):
		print "echo \"will use %s as a scratch file\"" % scratchFilename

	if (tempFilename != None):
		print "echo \"will write temporary files to %s\"" % tempFilename

	print "echo \"will write track file to     %s\"" % trackFilename

	if (bigWigFilename != None):
		print "echo \"will write bigwig file to    %s\"" % bigWigFilename

	# write command(s) to create track file

	print
	print "echo \"=== creating track %s ===\"" % trackId

	if (gzipInput):
		commands =  []
		command  =  ["time gzip -dc %s" % inputFilename]
		commands += [command]
		command  =  ["> %s" % scratchFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		trackSourceFilename = scratchFilename
	else:
		trackSourceFilename = inputFilename

	commands =  []
	command  =  ["time genodsp"]
	command  += ["--chromosomes=%s" % chromsFilename]
	command  += ["--show:uncovered"]
	command  += ["= input %s" % trackSourceFilename]
	command  += ["= binarize"]
	command  += ["= slidingsum W=%d D=W" % minLength]
	command  += ["= percentile %s W=%d --min=1/inf --quiet" \
			   % (densityThreshold,samplingStep)]
	if (gzipInput): command += ["= input %s --destroy" % trackSourceFilename]
	else:           command += ["= input %s"           % trackSourceFilename]
	command  += ["= binarize"]
	command  += ["= clump --average=percentile%s L=%d" \
			   % (densityThreshold,minLength)]
	commands += [command]

	if (gzipOutput):
		if (tempFilename != None):
			command  =  ["tee %s" % tempFilename]
			commands += [command]
		command  =  ["gzip"]
		commands += [command]

	command  =  ["> %s" % trackFilename]
	commands += [command]

	print
	print commands_to_pipeline(commands)

	# write command(s) to convert track file to bigwig

	if (bigWigFilename != None):
		print
		print "echo \"=== converting track %s to bigwig ===\"" % trackId

		if (tempFilename != None): trackInput = tempFilename
		else:                      trackInput = trackFilename

		commands =  []
		command  =  ["time bedGraphToBigWig"]
		command  += [trackInput]
		command  += [bigWigChromsFilename]
		command  += [bigWigFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		if (tempFilename != None):
			commands =  []
			command  =  ["rm %s" % tempFilename]
			commands += [command]
			print
			print commands_to_pipeline(commands)

		commands =  []
		command  =  ["make_bigwig_info"]
		command  += ["--url=%s" % bigWigUrl]
		command  += ["--name=\"%s discordant mates dense\"" % runName]
		command  += ["--desc=\"%s dense intervals in discordant mates mapped occupancy (${today})\"" \
				   % runName]
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


def commands_to_pipeline(commands):
	pipeline = []
	for (cmdNum,cmd) in enumerate(commands):
		if (cmdNum == 0): prefix = ""
		else:             prefix = "  | "

		if (cmd[0].startswith(">")):
			assert (cmdNum != 0)
			assert (len(cmd) == 1)
			prefix = "  "

		pipeline += [prefix + cmd[0]] 
		for line in cmd[1:]:
			pipeline += ["      " + line]

	return " \\\n".join(pipeline)


def do_filename_substitutition(s):
	if ("{base}" in s):
		assert (basePath != None)
		s = s.replace("{base}",basePath)
	if ("{run}" in s):
		assert (runName != None)
		s = s.replace("{run}",runName)
	return s


if __name__ == "__main__": main()

