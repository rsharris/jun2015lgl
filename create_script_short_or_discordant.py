#!/usr/bin/env python
"""
Create a cluster job file to combine the short insert coverage depth density
and discordant-mates density discriminator tracks.
"""

from sys  import argv,stdin,stdout,stderr,exit
from time import strftime


def usage(s=None):
	message = """

usage: create_script_short_or_discordant [options] > short_or_discordant.sh
  <sub>_<samp>_<type>        (required) run descriptor; for example, CS_NORM_PE
                             means subject "CS", sample "NORM", and type "PE";
                             other filenames can use "{run}" to refer to this
                             string
  --base=<path>              path prefix; other filenames can use "{base}" to
                             refer to this path
  --chromosomes=<filename>   read chromosome names and lengths from a file
                             (default is {base}/data/hg19.chrom_lengths)
  --input=<filename>         (required) first track file to process
  --input2=<filename>        (required) second track file to process
  --track=<filename>         (required) track file to create
                             (default is {base}/tracks/{run}.short_or_discordant)
  --tempinput=<filename>     temporary file to hold an input track file, if
                             needed; only needed if the second input track was
                             gzipped
                             (default is {base}/tracks/{run}.short_or_discordant.scratch)
  --temp=<filename>          temporary file to hold track file, if needed; only
                             needed if --gzip and --bigwig are both used
                             (default is {base}/tracks/{run}.short_or_discordant.temp)
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
                             other commands are copied "as is" """

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global basePath,runName
	global debug

	bashShebang = "#!/usr/bin/env bash"

	# parse args

	runName              = None
	basePath             = None
	chromsFilename       = None
	inputFilename1       = None
	inputFilename2       = None
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

		if (arg.startswith("--base=")) or (arg.startswith("--basepath=")) or (arg.startswith("--path=")):
			basePath = argVal
		elif (arg.startswith("--input=")) or (arg.startswith("--input1=")):
			inputFilename1 = argVal
		elif (arg.startswith("--input2=")):
			inputFilename2 = argVal
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

	if (inputFilename1 == None) or (inputFilename2 == None):
		usage("you have to give me two input track filenames")

	if (chromsFilename == None):
		chromsFilename = "{base}/data/hg19.chrom_lengths"

	if (trackName == None):
		trackName = "{base}/tracks/{run}.short_or_discordant"

	if (scratchFilename == None):
		scratchFilename = trackName + ".scratch"

	if (tempFilename == None) and (bigWigFilename != None) and (gzipOutput):
		tempFilename = trackName + ".temp"

	if (bigWigFilename != None):
		if (bigWigChromsFilename == None):
			bigWigChromsFilename = "{base}/temp/ucsc.hg19.chrom_lengths"
		if (bigWigUrl == None):
			usage("you have to give me a url for the bigwig file")

	trackId = "%s.short_or_discordant" % runName

	##########
	# perform filename substitution
	##########

	if   (basePath == None):       basePath = "."
	elif (basePath.endswith("/")): basePath = basePath[:-1]

	chromsFilename  = do_filename_substitutition(chromsFilename)

	# input track names

	inputFilename1 = do_filename_substitutition(inputFilename1)
	inputFilename2 = do_filename_substitutition(inputFilename2)

	gzipInput1 = False
	if   (inputFilename1.endswith(".gz")):      gzipInput1 = True
	elif (inputFilename1.endswith(".gzip")):    gzipInput1 = True
	elif (not inputFilename1.endswith(".dat")): inputFilename1 += ".dat"

	gzipInput2 = False
	if   (inputFilename2.endswith(".gz")):      gzipInput2 = True
	elif (inputFilename2.endswith(".gzip")):    gzipInput2 = True
	elif (not inputFilename2.endswith(".dat")): inputFilename2 += ".dat"

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

	if (gzipInput2):
		commands =  []
		command  =  ["time gzip -dc %s" % inputFilename2]
		commands += [command]
		command  =  ["> %s" % scratchFilename]
		commands += [command]

		print
		print commands_to_pipeline(commands)

		track2SourceFilename = scratchFilename
	else:
		track2SourceFilename = inputFilename2

	commands =  []
	if (gzipInput1): command = ["time gzip -dc %s" % inputFilename1]
	else:            command = ["time cat %s" % inputFilename1]
	commands += [command]

	command  =  ["genodsp"]
	command  += ["--chromosomes=%s" % chromsFilename]
	command  += ["--show:uncovered"]
	if (gzipInput2): command += ["= add %s --destroy" % track2SourceFilename]
	else:            command += ["= add %s"           % track2SourceFilename]
	command  += ["= binarize"]
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
		command  += ["--name=\"%s short OR discordant\"" % runName]
		command  += ["--desc=\"%s dense intervals in short inserts OR discordant mates (${today})\"" \
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

