#!/usr/bin/env python
"""
Create cluster job files to align a batch of reads to a genome, using bwa.
"""

from sys  import argv,stdin,stdout,stderr,exit
from re   import compile
from math import ceil


def usage(s=None):
	message = """

usage: make_bwa_jobs [options]
  --base=<path>            path prefix;  other filenames can use "{base}" to
                           refer to this path
  --ref=<filename>         (mandatory) name of reference fasta file
                           we also expect a bwa reference index in
                           <filename>.bwt, <filename>.pac, <filename>.ann,
                           <filename>.amb, and <filename>.sa
  --reads=<filename>       (mandatory,cumulative) name of query fastq file(s)
                           these must contain {mate} to indicate "1" or "2",
                           and may optionally contain a short identifier as a
                           prefix ending in ":";  e.g.
                             W1:reads/MQ854_K763WCJIIW_s_1_{mate}.fastq
  --phred64to33            convert fastq files from phred+64 to phred+33 'on
                           the fly'
  --readgroupinfo=<string> read group info string to shove into the sam header,
                           e.g. "ID:SRR748091 SM:SAMN01920490 PL:Illumina LB:lib1 PU:unknown"
                           (this is only supported for sampe)
  --id=<string>            unique identifier for this project
  --job=<path>             (mandatory) where to create job files
  --results=<path>         (mandatory) where to write alignment results
  --subblocks=[<L>:]<filespec>  create several jobs for each reads file, with
                           each job dealing with on "block" of the reads;
                           <filespec> is the name of a file suitable for use
                           with extract_blocks (a tabulated index); <L> is
                           the number of lines (in that index) to process in
                           each job;  if there is more than one reads file,
                           <filespec> must contain "{reads}" which is
                           replaced by the reads file name;  a typical
                           <filespec> is:
                             {reads}.tabulated
  --exe:bwa=<filename>     location of bwa executable
  --aln=<params>           (cumulative) parameters to pass to bwa aln
  --sampe=<params>         (cumulative) parameters to pass to bwa sampe
  --mem[=<params>]         (cumulative) parameters to pass to bwa mem
                           (by default we run "bwa aln" and "bwa sampe"
  --mem:unpaired[=<params>] (cumulative) parameters to pass to bwa mem, and
                           we're aligning reads as unpaired
  --sort=<params>          (cumulative) parameters to pass to samtools sort
  --memory:aln=<bytes>     run-time memory limit for aln (e.g. "5.6G")
                           add ":ulimit" to enforce this in the shell script
  --memory:sampe=<bytes>   run-time memory limit for sampe (e.g. "5.6G")
                           add ":ulimit" to enforce this in the shell script
  --memory:mem=<bytes>     run-time memory limit for mem (e.g. "5.6G")
                           add ":ulimit" to enforce this in the shell script
  --threads:aln=<number>   number of threads for aln
  --threads:sampe=<number> number of threads for sampe
  --threads:mem=<number>   number of threads for mem
  --initialize=<text>      (cumulative) shell command to add to job beginning
                           "shebang:bash" is mapped "#!/usr/bin/env bash"
                           other commands are copied "as is"
  --nojobnames             inhibit job names in submission list
  --outputas=sam           output in sam format (this is the default)
  --outputas=bam           output in bam format
  --outputas=bamsorted     output in sorted bam format

Typical command line:

    make_bwa_jobs \\ 
      --base="   /home/username/projects/orange" \\ 
      --ref="    {base}/reference/MalusDomestica0.fa" \\ 
      --id="     apple" \
      --reads="  {id}1:reads/MQ854_K763WCJIIW_s_1_{mate}.fastq" \\ 
      --job="    {base}/jobs" \\ 
      --results="{base}/results" \\ 
      --aln="    -I -q 15 -l 35 -k 2 -n 0.04 -o 2 -e 6" \\ 
      --sampe="  -P -n 100 -N 100" \\ 
      --memory:aln="  5.6G" \\ 
      --memory:sampe="6.2G" """

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global basePath,refFilename,readsFilespecs,jobDirectory,resultsDirectory
	global bwaProgramName,commandParams,phred64to33
	global bashInitializers,maxMemory,numThreads,readGroupInfo
	global jobId,giveJobsNames,jobNumber
	global outputAs
	global debug

	# parse the command line

	basePath           = None
	refFilename        = None
	readsFilespecs     = None
	phred64to33        = False
	readGroupInfo      = None
	jobId              = None
	jobDirectory       = None
	resultsDirectory   = None
	subsampleFileSpec  = None
	subsampleBlockSize = None
	bwaProgramName     = None
	aligner            = "aln"
	commandParams      = {}
	maxMemory          = {}
	numThreads         = {}
	bashInitializers   = ["set -eu"]
	giveJobsNames      = True
	outputAs           = "sam"
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1].strip()

		if (arg.startswith("--base=")) or (arg.startswith("--basepath=")) or (arg.startswith("--path=")):
			basePath = argVal
		elif (arg.startswith("--reference=")) or (arg.startswith("--ref=")):
			refFilename = argVal
		elif (arg.startswith("--exe:bwa=")) or (arg.startswith("--bwaexe=")):
			bwaProgramName = argVal
		elif (arg.startswith("--reads=")):
			if (readsFilespecs == None): readsFilespecs = []
			readsFilespecs += argVal.split(",")
		elif (arg == "--phred64to33"):
			phred64to33 = True
		elif (arg.startswith("--readgroupinfo=")):
			readGroupInfo = "\\t".join(argVal.split())
		elif (arg.startswith("--id=")):
			jobId = argVal
		elif (arg.startswith("--job=")):
			jobDirectory = argVal
		elif (arg.startswith("--results=")):
			resultsDirectory = argVal
		elif (arg.startswith("--subblocks=")):
			# we allow --subblocks=spec or --subblocks=L:spec
			subsampleFileSpec = subsampleBlockSize = None
			if (":" not in argVal):
				subsampleBlockSize = 1
				subsampleFileSpec  = argVal
			else:
				(subsampleBlockSize,subsampleFileSpec) = argVal.split(":",1)
				subsampleBlockSize = int(subsampleBlockSize)
				assert (subsampleBlockSize >= 1)
		elif (arg.startswith("--aln=")) or (arg.startswith("--params:aln=")):
			params = argVal
			if ("aln" not in commandParams): commandParams["aln"] = []
			commandParams["aln"] += [params]
		elif (arg.startswith("--sampe=")) or (arg.startswith("--params:sampe=")):
			params = argVal
			if ("sampe" not in commandParams): commandParams["sampe"] = []
			commandParams["sampe"] += [params]
		elif (arg.startswith("--mem=")) or (arg.startswith("--params:mem=")):
			params = argVal
			if ("mem" not in commandParams): commandParams["mem"] = []
			commandParams["mem"] += [params]
			aligner = "mem"
		elif (arg == "--mem"):
			aligner = "mem"
		elif (arg.startswith("--mem:unpaired=")) or (arg.startswith("--params:memunpaired=")):
			params = argVal
			if ("mem" not in commandParams): commandParams["mem"] = []
			commandParams["mem"] += [params]
			aligner = "mem unpaired"
		elif (arg == "--mem:unpaired"):
			aligner = "mem unpaired"
		elif (arg.startswith("--sort=")) or (arg.startswith("--params:sort=")):
			params = argVal
			if ("sort" not in commandParams): commandParams["sort"] = []
			commandParams["sort"] += [params]
		elif (arg.startswith("--memory:aln=")):
			if (argVal.endswith(":ulimit")):
				argVal = argVal[:-len(":ulimit")]
				maxMemory["aln"]        = argVal
				maxMemory["aln:ulimit"] = True
			else:
				maxMemory["aln"] = argVal
		elif (arg.startswith("--memory:sampe=")):
			if (argVal.endswith(":ulimit")):
				argVal = argVal[:-len(":ulimit")]
				maxMemory["sampe"]        = argVal
				maxMemory["sampe:ulimit"] = True
			else:
				maxMemory["sampe"] = argVal
		elif (arg.startswith("--memory:mem=")):
			if (argVal.endswith(":ulimit")):
				argVal = argVal[:-len(":ulimit")]
				maxMemory["mem"]        = argVal
				maxMemory["mem:ulimit"] = True
			else:
				maxMemory["mem"] = argVal
		elif (arg.startswith("--threads:aln=")):
			val = int(argVal)
			assert (val > 0)
			if (val > 1): numThreads["aln"] = val
		elif (arg.startswith("--threads:sampe=")):
			val = int(argVal)
			assert (val > 0)
			if (val > 1): numThreads["sampe"] = val
		elif (arg.startswith("--threads:mem=")):
			val = int(argVal)
			assert (val > 0)
			if (val > 1): numThreads["mem"] = val
		elif (arg.startswith("--initialize=")) or (arg.startswith("--init=")):
			if (argVal == "shebang:bash"):
				argVal = "#!/usr/bin/env bash"
			if (argVal == "set -eu"):
				bashInitializers = [x for x in bashInitializers if (x != "set -eu")]
			bashInitializers += [argVal]
		elif (arg == "--nojobnames"):
			giveJobsNames = False
		elif (arg.startswith("--outputas=")) and (argVal == "sam"):
			outputAs = "sam"
		elif (arg.startswith("--outputas=")) and (argVal == "bam"):
			outputAs = "bam"
		elif (arg.startswith("--outputas=")) and (argVal in ["bamsorted","sortedbam"]):
			outputAs = "sorted bam"
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# validate params

	assert (refFilename      != None)
	assert (readsFilespecs   != None)
	assert (jobDirectory     != None)
	assert (resultsDirectory != None)

	if (len(readsFilespecs) > 1) and (subsampleFileSpec != None):
		assert ("{reads}" in subsampleFileSpec)

	# perform filename substitution

	if (basePath != None) and (basePath.endswith("/")):
		basePath = basePath[:-1]

	refFilename = do_filename_substitutition(refFilename)

	for (ix,filename) in enumerate(readsFilespecs):
		readsFilespecs[ix] = do_filename_substitutition(filename)
		if (aligner != "mem unpaired"):
			assert ("{mate}" in readsFilespecs[ix])

	if (subsampleFileSpec != None): 
		subsampleFileSpec = do_filename_substitutition(subsampleFileSpec)

	jobDirectory = do_filename_substitutition(jobDirectory)
	if (jobDirectory.endswith("/")):
		jobDirectory = jobDirectory[:-1]

	resultsDirectory = do_filename_substitutition(resultsDirectory)
	if (resultsDirectory.endswith("/")):
		resultsDirectory = resultsDirectory[:-1]

	if (bwaProgramName == None):
		bwaProgramName = "bwa"
	else:
		bwaProgramName = do_filename_substitutition(bwaProgramName)

	for command in commandParams:
		for (ix,param) in enumerate(commandParams[command]):
			commandParams[command][ix] = do_filename_substitutition(param)

	for (ix,bashInitializer) in enumerate(bashInitializers):
		bashInitializer = do_filename_substitutition(bashInitializer)
		bashInitializers[ix] = bashInitializer

	if (aligner == "aln"):
		assert ("mem" not in commandParams)
		assert ("mem" not in maxMemory)
		assert ("mem" not in numThreads)
	elif (aligner == "mem"):
		assert ("aln"   not in commandParams)
		assert ("aln"   not in maxMemory)
		assert ("aln"   not in numThreads)
		assert ("sampe" not in commandParams)
		assert ("sampe" not in maxMemory)
		assert ("sampe" not in numThreads)
	elif (aligner == "mem unpaired"):
		assert ("aln"   not in commandParams)
		assert ("aln"   not in maxMemory)
		assert ("aln"   not in numThreads)
		assert ("sampe" not in commandParams)
		assert ("sampe" not in maxMemory)
		assert ("sampe" not in numThreads)

	if (phred64to33):
		assert (aligner in ["mem","mem unpaired"])

	# convert memory specifications

	if ("aln:ulimit" in maxMemory):
		multiplier = 1
		val = maxMemory["aln"]
		if (val.startswith("{threads}*")):
			if ("aln" in numThreads): multiplier = numThreads["aln"]
			val = val.replace("{threads}*","")
		elif (val.endswith("*{threads}")):
			if ("aln" in numThreads): multiplier = numThreads["aln"]
			val = val.replace("*{threads}","")
		maxMemory["aln"] = int_with_unit(val) * multiplier
	elif ("aln" in maxMemory):
		val = maxMemory["aln"]
		maxMemory["aln"] = int_with_unit(val)

	if ("sampe:ulimit" in maxMemory):
		multiplier = 1
		val = maxMemory["sampe"]
		if (val.startswith("{threads}*")):
			if ("sampe" in numThreads): multiplier = numThreads["sampe"]
			val = val.replace("{threads}*","")
		elif (val.endswith("*{threads}")):
			if ("sampe" in numThreads): multiplier = numThreads["sampe"]
			val = val.replace("*{threads}","")
		maxMemory["sampe"] = int_with_unit(val) * multiplier
	elif ("sampe" in maxMemory):
		val = maxMemory["sampe"]
		maxMemory["sampe"] = int_with_unit(val)

	if ("mem:ulimit" in maxMemory):
		multiplier = 1
		val = maxMemory["mem"]
		if (val.startswith("{threads}*")):
			if ("mem" in numThreads): multiplier = numThreads["mem"]
			val = val.replace("{threads}*","")
		elif (val.endswith("*{threads}")):
			if ("mem" in numThreads): multiplier = numThreads["mem"]
			val = val.replace("*{threads}","")
		maxMemory["mem"] = int_with_unit(val) * multiplier
	elif ("mem" in maxMemory):
		val = maxMemory["mem"]
		maxMemory["mem"] = int_with_unit(val)

	# read the sub-block tabulation files, to determine their lengths

	if ("subblocks" in debug):
		print >>stderr, "subsampleFileSpec  = %s" % subsampleFileSpec
		print >>stderr, "subsampleBlockSize = %d" % subsampleBlockSize

	readsToTabulation = None
	if (subsampleFileSpec == None):
		subsampleN = None
	else:
		readsToTabulation = {}
		subsampleN = {}
		subsampleW = 1
		for readsFilespec in readsFilespecs:
			# nota bene: we assume the tabulation files for both mates have the
			#            same number of lines
			(_,_,readsMatespec) = interpret_filespec(readsFilespec,replaceMate=False)
			subsampleFilename = subsampleFileSpec.replace("{reads}",readsMatespec)
			subsampleFilename = subsampleFilename.replace("{mate}", "1")
			numLines = number_of_lines_in(subsampleFilename)
			readsToTabulation[readsFilespec] = subsampleFilename
			subLinesList = []
			for startLine in xrange(0,numLines,subsampleBlockSize):
				endLine = min(startLine+subsampleBlockSize,numLines)
				subLinesList += ["%d..%d" % (startLine+1,endLine)]
			subsampleN[readsFilespec] = subLinesList
			subsampleW = max(subsampleW,len(str(len(subLinesList))))
			if ("subblocks" in debug):
				print >>stderr, "subsampleN[%s] = #%d [%s]" \
				              % (readsFilespec,len(subLinesList),",".join(subLinesList))

	# write the jobs

	jobs = []
	jobNumber = -1

	if (aligner == "aln"):
		for jobSpec in job_specs(readsFilespecs,mateSet=[1,2,"P"],subsampleN=subsampleN):
			jobNumber += 1

			(readsFilespec,mate,subK,subLines) = jobSpec
			if (subK != None): subK = "%0*d" % (subsampleW,subK)

			if   (mate == 1):   dependencies =  [jobNumber]
			elif (mate != "P"): dependencies += [jobNumber]

			if (type(mate) == int):
				jobInfo = create_aln_job(readsFilespec,mate,subK,subLines)
			elif (mate == "P"):
				jobInfo = create_sampe_job(readsFilespec,subK,subLines,dependencies)
			else:
				assert (False), "Internal Error: mate = \"%s\"" % mate

			jobs += [jobInfo]

	elif (aligner == "mem"):
		for jobSpec in job_specs(readsFilespecs,mateSet=["P"],subsampleN=subsampleN):
			jobNumber += 1

			(readsFilespec,mate,subK,subLines) = jobSpec
			if (subK != None): subK = "%0*d" % (subsampleW,subK)
			jobInfo = create_mem_job(readsFilespec,subK,subLines)
			jobs += [jobInfo]

	elif (aligner == "mem unpaired"):
		for jobSpec in job_specs(readsFilespecs,mateSet=[None],subsampleN=subsampleN):
			jobNumber += 1

			(readsFilespec,mate,subK,subLines) = jobSpec
			if (subK != None): subK = "%0*d" % (subsampleW,subK)
			jobInfo = create_mem_unpaired_job(readsFilespec,subK,subLines)
			jobs += [jobInfo]

	# print the jobs list file

	fn = [jobDirectory,"/"]
	if (jobId != None): fn += [jobId]
	fn += [".bwa_map"]
	fn += [".job_list"]
	fn += [".txt"]
	fn = "".join(fn)
	print >>stderr, "writing \"%s\"" % fn
	f  = file(fn, "wt")

	print >>f, "\n".join(jobs)

	f.close()


# create a "bwa aln" job

def create_aln_job(readsFilespec,mate,subK,subLines):
	(dataset,readsId,_) = interpret_filespec(readsFilespec)
	readsFilename = do_reads_filename_substitutition(readsFilespec,mate)

	# determine job name

	jobName = []
	if (jobId != None):  jobName += [jobId]
	jobName += [dataset]
	if (subK != None):   jobName += [subK]
	jobName += [str(mate)]
	jobName = ".".join(jobName)

	resultsFilename = resultsDirectory + "/" + jobName + ".sai"

	if (giveJobsNames):
		shortJobName = []
		if (readsId != None): shortJobName += [readsId]
		else:                 shortJobName += [str(jobNumber)]
		if (subK != None):    shortJobName += [subK]
		shortJobName += [str(mate)]
		shortJobName = "_".join(shortJobName)

	# create the job file, and the entry for the jobs list

	fn = jobDirectory + "/" + jobName + ".sh"
	print >>stderr, "writing \"%s\"" % fn
	jobF = file(fn,"wt")

	jobInfo = []
	if (giveJobsNames):       jobInfo += ["%d %s %s" % (jobNumber,shortJobName,fn)]
	else:                     jobInfo += ["%d %s" % (jobNumber,fn)]
	if ("aln" in maxMemory):  jobInfo += ["--memory=%d" % int(ceil(maxMemory["aln"]/1000000000.0))]
	if ("aln" in numThreads): jobInfo += ["--cores=%d" % numThreads["aln"]]

	# write the job file

	if (bashInitializers == None): myBashInitializers = []
	else:                          myBashInitializers = [x for x in bashInitializers]

	if ("aln" in maxMemory):
		memPages = int(ceil(maxMemory["aln"]/1024.0))
		myBashInitializers += ["ulimit -v %d # %s bytes" \
		                    % (memPages,commatize(memPages*1024))]

	if (myBashInitializers != []):
		print >>jobF, "\n".join(myBashInitializers)
		print >>jobF

	pipe = []
	if (subLines != None):
		command =  ["extract_block"]
		command += ["--block=%s.tabulated:%s" % (readsFilename,subLines)]
		command += ["--path=%s" % readsFilename]
		pipe += [" \\\n     ".join(command)]

	command = ["%s aln " % bwaProgramName]
	if ("aln" in commandParams): command += [" ".join(commandParams["aln"])]
	if ("aln" in numThreads):    command += ["-t %d" % numThreads["aln"]]
	command += [refFilename]
	if (pipe == []): command += [readsFilename]
	else:            command += ["/dev/stdin"]
	if (pipe == []): pipe += [" \\\n     ".join(command)]
	else:            pipe += ["  | %s" % "\\\n     ".join(command)]

	pipe += ["  > %s" % resultsFilename]

	print >>jobF, "time " + " \\\n".join(pipe)

	jobF.close()

	return " ".join(jobInfo)


# create a "bwa sampe" job

def create_sampe_job(readsFilespec,subK,subLines,dependencies):
	(dataset,readsId,_) = interpret_filespec(readsFilespec)

	readsFilename1 = do_reads_filename_substitutition(readsFilespec,1)
	readsFilename2 = do_reads_filename_substitutition(readsFilespec,2)

	if (dependencies == []): dependencies = None

	# determine job name

	jobName = []
	if (jobId != None):  jobName += [jobId]
	jobName += [dataset]
	if (subK != None):   jobName += [subK]
	jobName += ["{mate}"]
	jobName = ".".join(jobName)

	if   (outputAs == "bam"):        resultsExt = ".bam"
	elif (outputAs == "sorted bam"): resultsExt = "" # (samtools will add ".bam")
	else:                            resultsExt = ".sam"

	resultsFilename1 = resultsDirectory + "/" + jobName.replace(".{mate}",".1") + ".sai"
	resultsFilename2 = resultsDirectory + "/" + jobName.replace(".{mate}",".2") + ".sai"
	resultsFilename  = resultsDirectory + "/" + jobName.replace(".{mate}","")   + resultsExt
	jobName          =                          jobName.replace(".{mate}",".P")

	if (giveJobsNames):
		shortJobName = []
		if (readsId != None): shortJobName += [readsId]
		else:                 shortJobName += [str(jobNumber)]
		if (subK != None):    shortJobName += [subK]
		shortJobName += ["P"]
		shortJobName = "_".join(shortJobName)

	# create the job file, and the entry for the jobs list

	fn = jobDirectory + "/" + jobName + ".sh"
	print >>stderr, "writing \"%s\"" % fn
	jobF = file(fn,"wt")

	jobInfo = []
	if (giveJobsNames):         jobInfo += ["%d %s %s" % (jobNumber,shortJobName,fn)]
	else:                       jobInfo += ["%d %s" % (jobNumber,fn)]
	if (dependencies != None):  jobInfo += ["--depend=%s" % ",".join([str(x) for x in dependencies])]
	if ("sampe" in maxMemory):  jobInfo += ["--memory=%d" % int(ceil(maxMemory["sampe"]/1000000000.0))]
	if ("sampe" in numThreads): jobInfo += ["--cores=%d" % numThreads["sampe"]]

	# write the job file

	if (bashInitializers == None): myBashInitializers = []
	else:                          myBashInitializers = [x for x in bashInitializers]

	if ("sampe" in maxMemory):
		memPages = int(ceil(maxMemory["sampe"]/1024.0))
		myBashInitializers += ["ulimit -v %d # %s bytes" \
		                    % (memPages,commatize(memPages*1024))]

	if (myBashInitializers != []):
		print >>jobF, "\n".join(myBashInitializers)
		print >>jobF

	pipe = []
	command = ["%s sampe" % bwaProgramName]
	if ("sampe" in commandParams): command += [" ".join(commandParams["sampe"])]
	if (readGroupInfo != None):    command += ["-r \"@RG\\t%s\"" % readGroupInfo]
	if ("sampe" in numThreads):    command += ["-t %d" % numThreads["sampe"]]
	command += [refFilename]
	command += [resultsFilename1]
	command += [resultsFilename2]
	if (subLines == None):
		command += [readsFilename1]
		command += [readsFilename2]
	else:
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename1,subLines)]
		command += ["   --path=%s)" % readsFilename1]
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename2,subLines)]
		command += ["   --path=%s)" % readsFilename2]
	pipe += [" \\\n      ".join(command)]

	if (outputAs == "bam"):
		pipe += ["  | samtools view -Sb /dev/stdin"]
		pipe += ["  > %s" % resultsFilename]
	elif (outputAs == "sorted bam"):
		pipe += ["  | samtools view -Su -"]
		sortParams = ""
		if ("sort" in commandParams):
			jobIdForSort = jobId
			if (subK != None): jobIdForSort += "." + subK
			sortParams = " ".join(commandParams["sort"])
			sortParams = sortParams.replace("{base}",basePath)
			sortParams = sortParams.replace("{id}",jobIdForSort)
		pipe += ["  | samtools sort %s - -o %s.bam" % (sortParams,resultsFilename)]
	else:
		pipe += ["  > %s" % resultsFilename]

	print >>jobF, "time " + " \\\n".join(pipe)

	jobF.close()

	return " ".join(jobInfo)


# create a "bwa mem" job

def create_mem_job(readsFilespec,subK,subLines):
	(dataset,readsId,_) = interpret_filespec(readsFilespec)

	readsFilename1 = do_reads_filename_substitutition(readsFilespec,1)
	readsFilename2 = do_reads_filename_substitutition(readsFilespec,2)

	# determine job name

	jobName = []
	if (jobId != None):  jobName += [jobId]
	jobName += [dataset]
	if (subK != None):   jobName += [subK]
	jobName += ["{mate}"]
	jobName = ".".join(jobName)

	if   (outputAs == "bam"):        resultsExt = ".bam"
	elif (outputAs == "sorted bam"): resultsExt = "" # (samtools will add ".bam")
	else:                            resultsExt = ".sam"

	resultsFilename = resultsDirectory + "/" + jobName.replace(".{mate}","")   + resultsExt
	jobName         =                          jobName.replace(".{mate}",".P")

	if (giveJobsNames):
		shortJobName = []
		if (readsId != None): shortJobName += [readsId]
		else:                 shortJobName += [str(jobNumber)]
		if (subK != None):    shortJobName += [subK]
		shortJobName += ["P"]
		shortJobName = "_".join(shortJobName)

	# create the job file, and the entry for the jobs list

	fn = jobDirectory + "/" + jobName + ".sh"
	print >>stderr, "writing \"%s\"" % fn
	jobF = file(fn,"wt")

	jobInfo = []
	if (giveJobsNames):       jobInfo += ["%d %s %s" % (jobNumber,shortJobName,fn)]
	else:                     jobInfo += ["%d %s" % (jobNumber,fn)]
	if ("mem" in maxMemory):  jobInfo += ["--memory=%d" % int(ceil(maxMemory["mem"]/1000000000.0))]
	if ("mem" in numThreads): jobInfo += ["--cores=%d" % numThreads["mem"]]

	# write the job file

	if (bashInitializers == None): myBashInitializers = []
	else:                          myBashInitializers = [x for x in bashInitializers]

	if ("mem" in maxMemory):
		memPages = int(ceil(maxMemory["mem"]/1024.0))
		myBashInitializers += ["ulimit -v %d # %s bytes" \
		                    % (memPages,commatize(memPages*1024))]

	if (myBashInitializers != []):
		print >>jobF, "\n".join(myBashInitializers)
		print >>jobF

	mustUnzip = (readsFilename1.endswith(".gz") or readsFilename1.endswith(".gzip"))

	pipe = []
	command = ["%s mem" % bwaProgramName]
	if ("mem" in commandParams): command += [" ".join(commandParams["mem"])]
	if ("mem" in numThreads):    command += ["-t %d" % numThreads["mem"]]
	command += [refFilename]
	if (subLines == None) and (not phred64to33) and (not mustUnzip):
		command += [readsFilename1]
		command += [readsFilename2]
	elif (subLines == None) and (phred64to33) and (not mustUnzip):
		command += ["<(cat %s" % readsFilename1]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
		command += ["<(cat %s" % readsFilename2]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	if (subLines == None) and (not phred64to33) and (mustUnzip):
		command += ["<(gzip -dc %s)" % readsFilename1]
		command += ["<(gzip -dc %s)" % readsFilename2]
	elif (subLines == None) and (phred64to33) and (mustUnzip):
		command += ["<(gzip -dc %s" % readsFilename1]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
		command += ["<(gzip -dc %s" % readsFilename2]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	elif (subLines != None) and (not phred64to33) and (not mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename1,subLines)]
		command += ["   --path=%s)" % readsFilename1]
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename2,subLines)]
		command += ["   --path=%s)" % readsFilename2]
	elif (subLines != None) and (phred64to33) and (not mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename1,subLines)]
		command += ["   --path=%s" % readsFilename1]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename2,subLines)]
		command += ["   --path=%s" % readsFilename2]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	elif (subLines != None) and (not phred64to33) and (mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename1,subLines)]
		command += ["   --path=%s" % readsFilename1]
		command += ["   | gzip -dc)"]
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename2,subLines)]
		command += ["   --path=%s" % readsFilename2]
		command += ["   | gzip -dc)"]
	elif (subLines != None) and (phred64to33) and (mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename1,subLines)]
		command += ["   --path=%s" % readsFilename1]
		command += ["   | gzip -dc"]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename2,subLines)]
		command += ["   --path=%s" % readsFilename2]
		command += ["   | gzip -dc"]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	pipe += [" \\\n      ".join(command)]

	if (outputAs == "bam"):
		pipe += ["  | samtools view -Sb /dev/stdin"]
		pipe += ["  > %s" % resultsFilename]
	elif (outputAs == "sorted bam"):
		pipe += ["  | samtools view -Su -"]
		sortParams = ""
		if ("sort" in commandParams):
			jobIdForSort = jobId
			if (subK != None): jobIdForSort += "." + subK
			sortParams = " ".join(commandParams["sort"])
			sortParams = sortParams.replace("{base}",basePath)
			sortParams = sortParams.replace("{id}",jobIdForSort)
		pipe += ["  | samtools sort %s - -o %s.bam" % (sortParams,resultsFilename)]
	else:
		pipe += ["  > %s" % resultsFilename]

	print >>jobF, "time " + " \\\n".join(pipe)

	jobF.close()

	return " ".join(jobInfo)


# create a "bwa mem" job for unpaired reads

def create_mem_unpaired_job(readsFilespec,subK,subLines):
	(dataset,readsId,_) = interpret_filespec(readsFilespec,replaceMate=False)

	readsFilename = do_reads_filename_substitutition(readsFilespec,None)

	# determine job name

	jobName = []
	if (jobId != None):  jobName += [jobId]
	jobName += [dataset]
	if (subK != None):   jobName += [subK]
	jobName = ".".join(jobName)

	if   (outputAs == "bam"):        resultsExt = ".bam"
	elif (outputAs == "sorted bam"): resultsExt = "" # (samtools will add ".bam")
	else:                            resultsExt = ".sam"

	resultsFilename = resultsDirectory + "/" + jobName + resultsExt

	if (giveJobsNames):
		shortJobName = []
		if (readsId != None): shortJobName += [readsId]
		else:                 shortJobName += [str(jobNumber)]
		if (subK != None):    shortJobName += [subK]
		shortJobName = "_".join(shortJobName)

	# create the job file, and the entry for the jobs list

	fn = jobDirectory + "/" + jobName + ".sh"
	print >>stderr, "writing \"%s\"" % fn
	jobF = file(fn,"wt")

	jobInfo = []
	if (giveJobsNames):       jobInfo += ["%d %s %s" % (jobNumber,shortJobName,fn)]
	else:                     jobInfo += ["%d %s" % (jobNumber,fn)]
	if ("mem" in maxMemory):  jobInfo += ["--memory=%d" % int(ceil(maxMemory["mem"]/1000000000.0))]
	if ("mem" in numThreads): jobInfo += ["--cores=%d" % numThreads["mem"]]

	# write the job file

	if (bashInitializers == None): myBashInitializers = []
	else:                          myBashInitializers = [x for x in bashInitializers]

	if ("mem" in maxMemory):
		memPages = int(ceil(maxMemory["mem"]/1024.0))
		myBashInitializers += ["ulimit -v %d # %s bytes" \
		                    % (memPages,commatize(memPages*1024))]

	if (myBashInitializers != []):
		print >>jobF, "\n".join(myBashInitializers)
		print >>jobF

	mustUnzip = (readsFilename.endswith(".gz") or readsFilename.endswith(".gzip"))

	pipe = []
	command = ["%s mem" % bwaProgramName]
	if ("mem" in commandParams): command += [" ".join(commandParams["mem"])]
	if ("mem" in numThreads):    command += ["-t %d" % numThreads["mem"]]
	command += [refFilename]
	if (subLines == None) and (not phred64to33) and (not mustUnzip):
		command += [readsFilename]
	elif (subLines == None) and (phred64to33) and (not mustUnzip):
		command += ["<(cat %s" % readsFilename]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	if (subLines == None) and (not phred64to33) and (mustUnzip):
		command += ["<(gzip -dc %s)" % readsFilename]
	elif (subLines == None) and (phred64to33) and (mustUnzip):
		command += ["<(gzip -dc %s" % readsFilename]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	elif (subLines != None) and (not phred64to33) and (not mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename,subLines)]
		command += ["   --path=%s)" % readsFilename]
	elif (subLines != None) and (phred64to33) and (not mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename,subLines)]
		command += ["   --path=%s" % readsFilename]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	elif (subLines != None) and (not phred64to33) and (mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename,subLines)]
		command += ["   --path=%s" % readsFilename]
		command += ["   | gzip -dc)"]
	elif (subLines != None) and (phred64to33) and (mustUnzip):
		command += ["<(extract_block"]
		command += ["   --block=%s.tabulated:%s" % (readsFilename,subLines)]
		command += ["   --path=%s" % readsFilename]
		command += ["   | gzip -dc"]
		command += ["   | fastq_convert_phred --from=phred+64 --to=phred+33)"]
	pipe += [" \\\n      ".join(command)]

	if (outputAs == "bam"):
		pipe += ["  | samtools view -Sb /dev/stdin"]
		pipe += ["  > %s" % resultsFilename]
	elif (outputAs == "sorted bam"):
		pipe += ["  | samtools view -Su -"]
		sortParams = ""
		if ("sort" in commandParams):
			jobIdForSort = jobId
			if (subK != None): jobIdForSort += "." + subK
			sortParams = " ".join(commandParams["sort"])
			sortParams = sortParams.replace("{base}",basePath)
			sortParams = sortParams.replace("{id}",jobIdForSort)
		pipe += ["  | samtools sort %s - -o %s.bam" % (sortParams,resultsFilename)]
	else:
		pipe += ["  > %s" % resultsFilename]

	print >>jobF, "time " + " \\\n".join(pipe)

	jobF.close()

	return " ".join(jobInfo)


# perform filename substitutions

def do_filename_substitutition(s):
	if ("{base}" in s):
		assert (basePath != None)
		s = s.replace("{base}",basePath)
	return s


def do_reads_filename_substitutition(filespec,mate):
	if (":" in filespec):
		saveFilespec = filespec
		numFields = len(filespec.split(":"))
		assert (numFields == 2)
		filespec = filespec.split(":",1)[1]
		if (filespec == ""): filespec = saveFilespec

	if (mate != None):
		assert ("{mate}" in filespec)
		filespec = filespec.replace("{mate}",str(mate))

	return filespec


# extract info from a reads filespec

fileSpecRe = compile("(?P<mate>_*\{mate\}_*)")

def interpret_filespec(filespec,replaceMate=True):
	readsId = None
	if (":" in filespec):
		saveFilespec = filespec
		numFields = len(filespec.split(":"))
		assert (numFields == 2)
		(readsId,filespec) = filespec.split(":",1)
		if (readsId  == ""): readsId  = None
		if (filespec == ""): (readsId,filespec) = (None,saveFilespec)

	if (readsId != None):
		if (jobId != None): readsId = readsId.replace("{id}",jobId)

	saveFilespec = filespec
	if ("/" in filespec):
		filespec = filespec.split("/")[-1]

	if (filespec.endswith(".fastq")):
		filespec = filespec[:filespec.rfind(".fastq")]
	elif (filespec.endswith(".fq")):
		filespec = filespec[:filespec.rfind(".fq")]

	dataset = filespec
	if (replaceMate):
		m = fileSpecRe.search(filespec)			# find first match
		assert (m != None)
		(sIx,eIx) = (m.start(),m.end())
		m = fileSpecRe.search(filespec[eIx:])	# make sure there's no other match
		assert (m == None)

		replacement = ""
		if (sIx != 0) and (eIx != len(filespec)): replacement = "_"
		dataset = filespec[:sIx] + replacement + filespec[eIx:]

		while (dataset[-1] in [".","_"]):
			dataset = dataset[:-1]

	return (dataset,readsId,saveFilespec)


# generate job specs, tuples of the format (filespec,mate,subsample,subInfo)

def job_specs(readsFilespecs,mateSet=2,subsampleN=None):
	if (type(mateSet) == int):
		assert (mateSet > 0)
	else:
		assert (mateSet != [])

	if (subsampleN != None):
		assert (type(subsampleN) == dict)

	if (type(mateSet) != int):
		for readFilespec in readsFilespecs:
			if (subsampleN == None):
				for mate in mateSet:
					yield (readFilespec,mate,None,None)
			else:
				for (subK,subLines) in enumerate(subsampleN[readFilespec]):
					for mate in mateSet:
						yield (readFilespec,mate,subK+1,subLines)

	elif (mateSet == 1):
		for readFilespec in readsFilespecs:
			if (subsampleN == None):
				yield (readFilespec,None,None,None)
			else:
				for (subK,subLines) in enumerate(subsampleN[readFilespec]):
					yield (readFilespec,None,subK+1,subLines)

	else:
		for readFilespec in readsFilespecs:
			if (subsampleN == None):
				for mate in xrange(1,mateSet+1):
					yield (readFilespec,mate,None,None)
			else:
				for (subK,subLines) in enumerate(subsampleN[readFilespec]):
					for mate in xrange(1,mateSet+1):
						yield (readFilespec,mate,subK+1,subLines)


# number_of_lines_in--
#	Count the number of lines in a file.

def number_of_lines_in(filename):
	f = file(filename,"rt")
	numLines = 0
	for line in f:
		numLines += 1
	f.close()
	return numLines


# int_with_unit--
#	Parse a strings as an integer, allowing unit suffixes

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

	try:               return            int(s) * multiplier
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

