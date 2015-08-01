#!/usr/bin/env python
"""
Create a track info file to reference a bigwig file, suitable for display at
the UCSC browser.
"""

from sys import argv,stdin,stderr,exit


def usage(s=None):
	message = """
usage: make_bigwig_info [options] > bigwig_info_file
  --url=<string>                (required) URL of track's bigwig file
  --name=<string>               track's name field
  --description=<string>        track's description field
  --visibility=<string>         track's visibility field
  --autoscale=<string>          track's autoScale field ("on" or "off")
  --alwayszero=<string>         track's alwaysZero field ("on" or "off")
  --viewlimits=<string>         track's viewLimits field
  --maxheight=<string>          track's maxHeightPixels field
  --color=<r,g,b>               track's color field (e.g. red is 255,0,0)
  --position=<chrom:start-end>  track's initial browser position
  --graphtype=<string>          track's graphType field ("bar" or "points")"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse args

	url         = None
	name        = None
	description = None
	visibility  = "full"
	autoScale   = "off"
	alwaysZero  = None
	viewLimits  = "0:100"
	maxHeight   = None
	color       = None
	position    = None
	graphType   = None
	debug       = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--url=")) or (arg.startswith("--URL=")):
			url = argVal
		elif (arg.startswith("--name=")):
			name = argVal
		elif (arg.startswith("--description=")) or (arg.startswith("--desc=")):
			description = argVal
		elif (arg.startswith("--visibility=")):
			visibility = argVal
		elif (arg.startswith("--autoscale=")) or (arg.startswith("--autoScale=")):
			autoScale = argVal
		elif (arg.startswith("--alwayszero=")) or (arg.startswith("--alwaysZero=")):
			alwaysZero = argVal
		elif (arg.startswith("--viewlimits=")) or (arg.startswith("--viewLimits=")):
			viewLimits = argVal
		elif (arg.startswith("--maxheight=")) or (arg.startswith("--maxHeightPixels=")):
			maxHeight = argVal
		elif (arg.startswith("--color=")):
			color = argVal
			if (color.startswith("(")) and (color.endswith(")")):
				color = color[1:-1].split(",")
				color = [float(v) for v in color]
				color = [max(0.0,min(1.0,v)) for v in color]
				color = [int(round(255*v)) for v in color]
				color = ",".join([str(v) for v in color])
		elif (arg.startswith("--position=")) or (arg.startswith("--pos=")):
			position = argVal
		elif (arg.startswith("--graphtype=")) or (arg.startswith("--graphType=")):
			graphType = argVal
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (url == None):
		usage("You have to give me the URL of the bigwig file")

	# output the track info

	if (position != None):
		print "browser position %s" % position

	line =  ["track"]
	line += ["type=bigWig"]
	if (name        != None): line += ["name=\"%s\""        % name]
	if (True):                line += ["bigDataUrl=%s"      % url]
	if (description != None): line += ["description=\"%s\"" % description]
	if (visibility  != None): line += ["visibility=%s"      % visibility]
	if (autoScale   != None): line += ["autoScale=%s"       % autoScale]
	if (alwaysZero  != None): line += ["alwaysZero=%s"      % alwaysZero]
	if (viewLimits  != None): line += ["viewLimits=%s"      % viewLimits]
	if (maxHeight   != None): line += ["maxHeightPixels=%s" % maxHeight]
	if (color       != None): line += ["color=%s"           % color]
	if (graphType   != None): line += ["graphType=%s"       % graphType]

	print " ".join(line)


if __name__ == "__main__": main()
