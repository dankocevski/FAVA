#!/usr/bin/env python

print "\nFAVA Followup Analysis Tool v1.0"
print "Support Contact: Daniel Kocevski (daniel.kocevski@nasa.gov)\n"
print "Importing modules:"
import os
import subprocess
import pyfits
import fileinput
import time
import random
from GtApp import GtApp
import sys
import glob
import numpy
import fileinput
import pyLikelihood
from UnbinnedAnalysis import *
import subprocess
from functools import wraps
import errno
import os
import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plot
from matplotlib import dates as mdates
import matplotlib.colors as mcolors
import getpass
from make2FGLxml import *
#from make2FGLxml_FAVA import *
#from make3FGLxml import *
import traceback


##########################################################################################

class TimeoutError(Exception):
	pass

class timeout:
	def __init__(self, seconds=1, error_message='Timeout'):
		self.seconds = seconds
		self.error_message = error_message
	def handle_timeout(self, signum, frame):
		raise TimeoutError(self.error_message)
	def __enter__(self):
		signal.signal(signal.SIGALRM, self.handle_timeout)
		signal.alarm(self.seconds)
	def __exit__(self, type, value, traceback):
		signal.alarm(0)


##########################################################################################

def astro_query(options, store='store', verbose=True):
	command = "/afs/slac/u/gl/glast/astroserver/prod/astro"
	for item in options:
		if options[item] is not None:
			command += ' --%s %s' % (item, options[item])
	command += ' ' + store
	if verbose:
		print command
	subprocess.call(command, shell=True)
	for item in glob.glob('*_fix_checksums.sh'):
		os.remove(item)

##########################################################################################

def getFT1(ra, dec, radius, tmin, tmax, emin=100., emax= 3e5,
		   evclassname="Source", outfile='FT1.fits', verbose=True):
#    options = {'event-sample' : 'P7.6_P120_ALL',
#    options = {'event-sample' : 'P7.6_P130_ALL',
	options = {'event-sample' : 'P7_P202_ALL',
			   'output-ft1' : outfile,
			   'minEnergy' : emin,
			   'maxEnergy' : emax,
			   'minTimestamp' : tmin,
			   'maxTimestamp' : tmax,
			   'ra' : ra,
			   'dec' : dec,
			   'radius' : radius}
	astro_query(options, store='store', verbose=verbose)


##########################################################################################

def getFT2(tmin, tmax, outfile='FT2.fits', verbose=True):
	options = {'event-sample' : 'P7_P202_ALL',
			   'output-ft2-30s' : outfile,
			   'minTimestamp' : tmin,
			   'maxTimestamp' : tmax}
	astro_query(options, store='storeft2', verbose=verbose)


##########################################################################################

def make_ds9_reg(ra, dec, outfile='CandidateSource.reg'):
	tpl = """# Region file format: DS9 version 4.0
# Filename: cmap_sum.fits
global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source
fk5
point(%.3f,%.3f) # point=circle
"""
	output = open(outfile, 'w')
	output.write(tpl % (ra, dec))
	output.close()


##########################################################################################

def writeXML(self, ra, dec, galpropModel, isotropicModel, outfile='SourceModel.xml'):
	output = open(outfile, 'w')
	output.write(Source_xml % (ra, dec, galpropModel, isotropicModel))
	output.close()


##########################################################################################

def AddCandidateSource(ra, dec, xmlModel):

	CandidateSource = """
<source name="CandidateSource" type="PointSource">
	<spectrum type="PowerLaw2">
	  <parameter free="1" max=".10000E+07" min=".10000E-09" name="Integral" scale="1" value="1e-4"/>
	  <parameter free="1" max="-1.0" min="-5.0" name="Index" scale="1.0" value="-2.1"/>
	  <parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="100.0"/>
	  <parameter free="0" max="300000.0" min="20.0" name="UpperLimit" scale="1.0" value="1e5"/>
	</spectrum>
	<spatialModel type="SkyDirFunction">
	  <parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="%.3f"/>
	  <parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="%.3f"/>
	</spatialModel>
  </source>
</source_library>
""" % (float(ra), float(dec))

	print 'Adding a candidate source to the xml model...'
	xmlModelContents = open(xmlModel).read()
	xmlModelContents = xmlModelContents.replace('</source_library>',CandidateSource)
	xmlModelContentsAppended = open(xmlModel,'w')
	xmlModelContentsAppended.write(xmlModelContents)
	xmlModelContentsAppended.close()
	print 'Done.'


##########################################################################################

# def ModifySourceModel(xmlModel, RemoveSource=True):

# 	xmlModelModified = xmlModel.replace('.xml','_Modified.xml')

# 	# Open the files
# 	infile = open(xmlModel,'r')
# 	outfile = open(xmlModelModified,'w')

# 	# Read in the lines from the original flare file
# 	lines = infile.readlines()

# 	# Loop through each line 
# 	doFix = True
# 	eraseLine = False

# 	for line in lines:

# 		newline = line

# 		# Fix all non-diffuse model components
# 		if 'type="DiffuseSource">' in newline:
# 			print 'Keeping Component:'
# 			print newline.rstrip()
# 			doFix = False

# 		# Make sure the CandidateSource is not affected
# 		if 'name="CandidateSource"' in newline:
# 			print 'Keeping Component:'
# 			print newline.rstrip()
# 			doFix = False

# 		if '<source name=' in line:
# 			SourceName = line

# 		if 'free="1"' in line and doFix == True:
# 			print 'Fixing Component:'
# 			print SourceName.rstrip()
# 			newline = line.replace('free="1"','free="0"')

# 		if '</source>' in line and doFix == False:
# 			doFix = True

# 		# Remove the candidate source from the model    
# 		if 'CandidateSource' in line and RemoveSource == True:
# 			newline = ''
# 			eraseLine = True

# 		if '</source>' in line and eraseLine == True:
# 			newline = ''
# 			eraseLine = False

# 		if eraseLine == True:
# 			newline = ''

# 		# Save the modified line
# 		outfile.write(newline)

# 	print "Writing modified xml file to: %s" % xmlModelModified
# 	infile.close()
# 	outfile.close()
# 	print "Done."


##########################################################################################

def ModifySourceModel(xmlModel, xmlModelModified, candidateSource='CandidateSource', RemoveSource=True):

	# Open the files
	infile = open(xmlModel,'r')
	outfile = open(xmlModelModified,'w')

	# Read in the lines from the original flare file
	lines = infile.readlines()

	# Loop through each line 
	doFix = True
	eraseLine = False

	for line in lines:

		newline = line

		# Fix all non-diffuse model components
		if 'type="DiffuseSource">' in newline:
			print 'Keeping Component:'
			print newline.rstrip()
			doFix = False

		# Make sure the CandidateSource is not affected
		if 'name="%s"' % candidateSource in newline:
			print 'Keeping Component:'
			print newline.rstrip()
			doFix = False

		if '<source name=' in line:
			SourceName = line

		if 'free="1"' in line and doFix == True:
			print 'Fixing Component:'
			print SourceName.rstrip()
			newline = line.replace('free="1"','free="0"')

		if '</source>' in line and doFix == False:
			doFix = True

		# Remove the candidate source from the model    
		if candidateSource in line and RemoveSource == True:
			print 'Removing Component:'
			print candidateSource
			newline = ''
			eraseLine = True

		if '</source>' in line and eraseLine == True:
			newline = ''
			eraseLine = False

		if eraseLine == True:
			newline = ''

		# Save the modified line
		outfile.write(newline)

	print "Writing modified xml file to: %s" % xmlModelModified
	infile.close()
	outfile.close()
	print "Done."

##########################################################################################

def ExtractSources(xmlModel, sourceNames):

	# Open the files
	infile = open(xmlModel,'r')
	outfile = open(sourceNames,'w')

	# Read in the lines from the original flare file
	lines = infile.readlines()

	# Go through each line and extract the source names
	SourceNames =[]
	Values = []
	for line in lines:
		if '<source name=' in line:

			# Extract the source name
			SourceName = line
			SourceName = SourceName.replace('  <source name="','')
			SourceName = SourceName[0:SourceName.find('"')]
			SourceName = SourceName.replace('_2FGL', '2FGL_')

			# Extract the url value to the passed
			Value = SourceName
			Value = SourceName.replace('+','%2B')

			# Add the source names and values to their respective arrays
			SourceNames.append(SourceName)
			Values.append(Value)



	# Remove the candidate source, and the galactic, and isotropic components
	SourceNames = SourceNames[1:-2]

	# Add a none source
	SourceNames.insert(0,'none')

	# Write the source names to a file
	for Value, Source in zip(Values, SourceNames):
		line = '<option value="' + Value + '">' + Source + '</option>\n'
		outfile.write(line)

	print "Writing source names file to: %s" % sourceNames
	infile.close()
	outfile.close()
	print "Done."



##########################################################################################

def SetPfilesDirectory(pfile_dir):
	"""each thread/job which uses FTOOLS must have its own
	PFILES dir"""
	
	try:
		hea_pdir = os.getenv("HEADAS")+"/syspfiles/"
		gt_pdir  = os.getenv("INST_DIR")+"/syspfiles/"
	except:
		print '\n*** Science tools has not be inialized! ***\n'
		print 'Exiting!'
		sys.exit()
	

	if(os.path.isdir(pfile_dir)==False):
		print "\nCreating custom pfile directory:\n%s" % pfile_dir        
		cmd = "mkdir -p "+pfile_dir
		os.system(cmd)

	# --- now copy all the .par files
	cmd = "cp %s/*.par %s"%(hea_pdir,pfile_dir)
	os.system(cmd)
	cmd = "cp %s/*.par %s"%(gt_pdir,pfile_dir)
	os.system(cmd)
	
	# ---Set the new environmental PFILES    
	#os.putenv("PFILES",pfile_dir)
	os.environ['PFILES'] = pfile_dir

	# --- test
#   print "Testing: temporary PFILES location "
#   print os.listdir(pfile_dir)

	return pfile_dir


##########################################################################################	

def customRALabel(deg):
	if (deg == 360):
		deg = ''
	return deg

def customDecLabel(deg):
	return deg

##########################################################################################	

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

##########################################################################################

# def MakeButterflyPlot(sourceName):
# 	E1 = like2.model[sourceName].funcs['Spectrum'].getParam('LowerLimit').value()
#     E2 = like2.model[sourceName].funcs['Spectrum'].getParam('UpperLimit').value()
# 	gamma = like2.model[sourceName].funcs['Spectrum'].getParam('Index').value()
# 	I = like2.model[sourceName].funcs['Spectrum'].getParam('Integral').value()
# 	cov_gg = like2.covariance[16][16]
# 	cov_II = like2.covariance[15][15]
# 	cov_Ig = like2.covariance[15][16]
# 	print "Index: " + str(gamma) + " +/- " + str(math.sqrt(cov_gg))
# 	print "Integral: " + str(I) + " +/- " + str(math.sqrt(cov_II))
# 	print E1,E2
# 	epsilon = (E2/E1)**(1-gamma)
# 	logE0 = (math.log(E1) - epsilon*math.log(E2))/(1-epsilon) + 1/(gamma-1) + cov_Ig/(I*cov_gg)
# 	E0 = math.exp(logE0)
# 	print E0


	
##########################################################################################

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

##########################################################################################

def RemoveSources(DuplicateSources, Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique):

	for Source in DuplicateSources:
		i = numpy.where(Names_SourcesUnique == Source)
		Names_SourcesUnique = numpy.delete(Names_SourcesUnique,i)
		RA_SourcesUnique = numpy.delete(RA_SourcesUnique,i)
		DEC_SourcesUnique = numpy.delete(DEC_SourcesUnique,i)
		pass

	return Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique

##########################################################################################

def ExtractCoordinates(xmlModel, Source):

	# Define some default values
	SourceFound = False
	SourceRA = 'NA'
	SourceDEC = 'NA'

	# Loop through the xml file and look for the requested source 
	for line in fileinput.input([xmlModel]):
		if Source in line:
			SourceFound = True

		# Extract the RA
		if ('name="RA"' in line) and (SourceFound == True):
			SourceRA = line[line.find('value='):]
			SourceRA = SourceRA.replace('/>','')
			SourceRA = SourceRA.replace('value=','')
			SourceRA = SourceRA.replace('"','')

		# Extract the Dec			
		if ('name="DEC"' in line) and (SourceFound == True):
			SourceDEC = line[line.find('value='):]
			SourceDEC = SourceDEC.replace('/>','')			
			SourceDEC = SourceDEC.replace('value=','')
			SourceDEC = SourceDEC.replace('"','')

		# Stop recording values
		if ('</source>' in line):
			SourceFound = False

	fileinput.close()

	return SourceRA, SourceDEC

##########################################################################################

def SubmitJobs(flarefile,Test=False):

	# Read in the flare information
	nums, ras, decs, lbin, bbin, gall, galb, tmins, tmaxs, sigma, avnev, nev, he_nev, he_avnev, he_sigma, sundist, varindex, favasrc, fglassoc, assoc = numpy.loadtxt(flarefile, unpack=True, comments='#', dtype = str)

	# Find the flares with greater than a positive 3 sigma significance at either low or high energy
	positiveFlares = numpy.where( (sigma.astype(float) > 3) | (he_sigma.astype(float) > 3))[0]
	negativeFlares = numpy.where( (sigma.astype(float) < -3) | (he_sigma.astype(float) < -3))[0]

	# Find the unassociated positive flares
	unassociatedFlares = numpy.where(fglassoc[positiveFlares] == 'none')[0]

	# Find the associated positive flares
	associatedFlares = numpy.where(fglassoc[positiveFlares] != 'none')[0]


	# Keep only the positive 3 sigma flares
	nums = nums[positiveFlares]
	ras = ras[positiveFlares]
	decs = decs[positiveFlares]
	tmins = tmins[positiveFlares]
	tmaxs = tmaxs[positiveFlares]
	fglassoc = fglassoc[positiveFlares]

	# Get the start and stop METs for the run
	tmin = tmins[0]
	tmax = tmaxs[0]	

	print "\nLaunching followup analyses on sources in FAVA run: %s to %s" % (tmin, tmax)
	print "Number of sources with positive flares: %s" % nums[-1]

	# Define some defaults
	FAVAHome = '/nfs/slac/g/ki/ki08/kocevski/FAVA'
	RunDirectory = "%s/favaonline/data/weekly/P7SOURCE_V6MC/maps/%s_%s/" % (FAVAHome, tmin, tmax)
	OutputDirectory = "%s/favaonline/data/weekly/P7SOURCE_V6MC/maps/%s_%s/SourceAnalysis" % (FAVAHome, tmin, tmax)
	LogDirectory = "%s/favaonline/data/weekly/P7SOURCE_V6MC/maps/%s_%s/Logs" % (FAVAHome, tmin, tmax)
	JobsDirectory = "%s/favaonline/data/weekly/P7SOURCE_V6MC/maps/%s_%s/Jobs" % (FAVAHome, tmin, tmax)

	# Creating the necessary directories
	if(os.path.isdir(OutputDirectory)==False):
		print "\n >> Creating Directory: " + OutputDirectory
		cmd = "mkdir " + OutputDirectory
		os.system(cmd)
		
	if(os.path.isdir(LogDirectory)==False):
		print "\n >> Creating Directory: " + LogDirectory
		cmd = "mkdir " + LogDirectory
		os.system(cmd)	
               
  	if(os.path.isdir(JobsDirectory)==False):
  		print "\n >> Creating Directory: " + JobsDirectory
  		cmd = "mkdir " + JobsDirectory
  		os.system(cmd)

	# Move into the output directory 
	os.chdir(OutputDirectory)
  		
  	# Mark the start time
  	startTime = time.time()
	
	# Prepare the data for each flare
	counter = 1
	for num, ra, dec, tmin, tmax, association in zip(nums, ras, decs, tmins, tmaxs, fglassoc):

		# Get the name
		sourceName = 'Source%s' % num
		print "\nSource %s:" % num
		
		# Construct the command
		command = FAVAHome + "/python/FAVAFollowupAnalysis.py sourceName=%s ra=%s dec=%s tmin=%s tmax=%s association=%s" % (sourceName, ra, dec, tmin, tmax, association)
		logfile = LogDirectory + "/%s.log" % sourceName

		# Construct the process call
		process = 'bsub -o ' + logfile + ' -J ' + sourceName + ' -q xlong -R rhel60 -g ' + JobsDirectory + ' "' + command + '"'
	
		if (Test == False):

			# Start the process and wait 600 seconds before submitting another request
			print process
			subprocess.call(process, shell=True)

		else:

			# Print the command, but don't submit the process (for testing purposes only)
			print command

		counter = counter + 1



	time.sleep(30)
	print "\n"

	nJobs = -1
	while nJobs != 0:		

		# Determine the number of jobs still running
		command = "bjobs -g %s | grep 'RUN' | wc" % JobsDirectory	
		process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
		lines = process.stdout.readlines()
	
		# Display the results and the elapsed time since the start of the analysis
		if 'No unfinished job found' in lines[0]:
			nJobs = 0
			print '\nDone.\n'
							
		else:
		
			nJobs = lines[0].split()[0]
			print "Jobs Remaining: %s" % nJobs
			
			# Display the elapsed time
			elapsedTimeSeconds = time.time() - startTime
			elapsedTimeMinutes = elapsedTimeSeconds/60.0
			elapsedTimeHours = elapsedTimeMinutes/60.0

			if elapsedTimeMinutes > 60:				
				print "Elapsed Time: %.2f Hours" % elapsedTimeHours
			else:
				print "Elapsed Time: %.2f Minutes" % elapsedTimeMinutes 

			time.sleep(120)

		
	# # Gather the results
	# command = "grep 'TS value' %s/likelihoodResults_MaxTS*.txt" % OutputDirectory	
	# process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
	# lines = process.stdout.readlines()
	
	# # Read in the results
	# sourceNumbers = []
	# TSs = []
	# for line in lines:
	# 	try:
	# 		sourceNumber = line.split()[0].split('/')[-1].replace(".txt:'TS",'').replace('likelihoodResults_MaxTS_Source','').strip()
	# 		TS = "%.2f" % float(line.split()[2].replace("'",'').replace(',','').strip())
	# 		sourceNumbers.append(sourceNumber)
	# 		TSs.append(TS)
	# 	except:
	# 		TSs.append('NA')

	# Gather the results
	command = "grep 'MaxTS:' %s/likelihoodResults_*.txt" % OutputDirectory	
	process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
	lines = process.stdout.readlines()
	
	# Read in the results
	sourceNumbers = []
	TSs = []
	for line in lines:
		try:
			sourceNumber = line.split()[0].split('/')[-1].replace(".txt:MaxTS:",'').replace('likelihoodResults_Source','').strip()
			TS = "%.2f" % float(line.split()[1].replace("'",'').replace(',','').strip())
			sourceNumbers.append(sourceNumber)
			TSs.append(TS)
		except:
			TSs.append('NA')



	# Put the results in a dictionary for easy retrieval later
	Results = {key:value for key, value in zip(sourceNumbers,TSs)}
	
	# Define a new modified flare list
	flarefileModified = flarefile.replace('.ff','_Modified.ff')
	
	# Open the flare file and the modified flare file
	infile = open(flarefile,'r')
	outfile = open(flarefileModified,'w')
	
	# Read in the lines from the original flare file
	lines = infile.readlines()
	
	# Loop through each line in the flare file and add the TS values by retrieving them from the results dictionary
	print "\nSummary of likelihood results:"
	for line in lines:
	
		if '#' in line:
			outfile.write(line)
			continue
			
		sourceNumber = line.split()[0].strip()
		
		try:
			TS = Results[sourceNumber]
		except:
			TS = 'NA'
			
		print "Source %s  %s" % (sourceNumber, TS) 
		newline = line.replace('\n'," %s\n" % TS)
		outfile.write(newline)
	
	infile.close()
	outfile.close()


	print "\nResults saved to: %s" % flarefileModified

	# Clean up
	print "\nCleaning up..."

	# Delete working files
	# os.system("rm %s/LC_*.*" % OutputDirectory)
	# os.system("rm %s/ft1_*.*" % OutputDirectory)
	# os.system("rm %s/expmap_*.*" % OutputDirectory)
	# os.system("rm %s/bexpmap_*.*" % OutputDirectory)
	# os.system("rm %s/ltcube_*.*" % OutputDirectory)
	# os.system("rm %s/ccube_*.*" % OutputDirectory)
	# os.system("rm %s/cmap_*.*" % OutputDirectory)
	# os.system("rm %s/srcmap_*.*" % OutputDirectory)
	# os.system("rm %s/modelmap_*.fits" % OutputDirectory)
	# os.system("rm %s/ModelSource*_Fit.xml" % OutputDirectory)
	# os.system("rm %s/ModelSource*_Fit_Modified.xml" % OutputDirectory)

	# Change the group ownership of the output directory so that the users can edit the results from the ground-glast machine
	os.system("chgrp -R glground %s" % OutputDirectory)


	print 'Done.\n'

	
	return


##########################################################################################

def UnbinnedLikelihood(sourceName, ra, dec, tmin, tmax, association='none', usegttsmap=True):

	import traceback
	import pyfits

	print "\nPerforming analysis on %s:" % sourceName
	print "tmin = %s, tmax = %s" % (tmin, tmax)
	print "RA = %s, Dec = %s" % (ra, dec)

	# Wait a period of time before starting in order to not crash the asf/nsf disks at SLAC
	try:
		sourceNumber = int(sourceName.replace('Source'))
		waitTime = sourceNumber * 20
	except:
		waitTime = random.random()*600

	print "\nWaiting %i seconds before starting..." % waitTime
	#time.sleep(waitTime)
	print "Proceeding."

	# Setup some defaults directories
	FAVAHome = '/nfs/slac/g/ki/ki08/kocevski/FAVA'
	DataDirectory = FAVAHome + '/favaonline/data/weekly/P7SOURCE_V6MC/ft1'
	OutputDirectory = "%s/favaonline/data/weekly/P7SOURCE_V6MC/maps/%s_%s/SourceAnalysis" % (FAVAHome, tmin, tmax)
	FT1Directory = FAVAHome + '/favaonline/data/weekly/P7SOURCE_V6MC/ft1'
	FT2Directory = FAVAHome + '/favaonline/data/weekly/ft2'
	extendedSourcesDirectory = FAVAHome + '/extendedSources/Templates/'
	dtsmapDirectory = "%s/favaonline/data/weekly/P7SOURCE_V6MC/maps/%s_%s/SourceAnalysis/dtsmap_%s" % (FAVAHome, tmin, tmax, sourceName)

	# Move into the output directory 
	os.chdir(OutputDirectory)

	# Defind the temporary scratch directories
	JobID = os.environ.get('LSB_JOBID')
	Username = getpass.getuser()
	ScratchDirectory = "/scratch/%s/%s/" % (Username, JobID)

	# Define the pfile directory
	if JobID == None:
		PFILESDirectory = "%s/pfiles_%s/" % (OutputDirectory, sourceName)	
	else:
		PFILESDirectory = "%s/pfiles/" % ScratchDirectory

	# Remove any pre-existing pfiles
	if(os.path.isdir(PFILESDirectory)==True):
		cmd = "rm -r %s" % PFILESDirectory
		os.system(cmd)			

	# Set the new pfiles directory
	SetPfilesDirectory(PFILESDirectory)
	print os.environ.get('PFILES')

	# Import the necessary gtapps	
	gtselect = GtApp('gtselect')
	gtmktime = GtApp('gtmktime')
	gtexpmap = GtApp('gtexpmap')
	gtbin = GtApp('gtbin')
	gtltcube = GtApp('gtltcube')
	gtexpcube2 = GtApp('gtexpcube2')
	gtdiffrsp = GtApp('gtdiffrsp')
	gtlike = GtApp('gtlike')
	gttsmap = GtApp('gttsmap')
	gtsrcmaps = GtApp('gtsrcmaps')
	gtmodel = GtApp('gtmodel')
	gtfindsrc = GtApp('gtfindsrc')

	# Setup extraction parameters
	radius = 20					# Degrees	
	emin = 100					# MeV
	emax = 1e5					# MeV
	evclass = 2					# Source Class
	zmax = 100					# MeV
	binsize = 0.25				# Degrees
	irfs = 'P7REP_SOURCE_V15'	# LAT Instrument Response Function
	optimizer = 'MINUIT'		# Minimization Algorithm
	
	# Setup the FT1 and FT2 files
	ft1file = "%s/ft1_%s_%s.fits" % (FT1Directory, tmin, tmax)
	ft2file = "%s/ft2_%s_%s.fits" % (FT2Directory, tmin, tmax)
	
	# Define the working files
	filteredEvents = '%s/ft1_filtered_%s.fits' % (OutputDirectory, sourceName)
	filteredEventsGTI = '%s/ft1_filteredGTI_%s.fits' % (OutputDirectory, sourceName)
	cmap = '%s/cmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	ccube = '%s/ccube_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	xmlModel = '%s/Model%s.xml' % (OutputDirectory,sourceName)
	xmlModelFit = '%s/Model%s_Fit.xml' % (OutputDirectory,sourceName)
	xmlModelFitModified = '%s/Model%s_Fit_Modified.xml' % (OutputDirectory,sourceName)
	Catalog2FGL = FAVAHome + '/catalogs/gll_psc_v08.fit'
	Catalog3FGL = FAVAHome + '/catalogs/gll_psc4yearsource_v12r2_Modified.fit'
	galpropModel = FAVAHome + '/diffuseModels/gll_iem_v05.fit'
	isotropicModel = FAVAHome + '/diffuseModels/iso_source_v05.txt'		
	ltcube = '%s/ltcube_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	expmap = '%s/expmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	bexpmap = '%s/bexpmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	bexpmap_allSky = '%s/bexpmap_allSky_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	srcmap = '%s/srcmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	likelihoodResults = '%s/likelihoodResults_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)
	tsmap = '%s/tsmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	dtsmapFile = '%s/dtsmap_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	srcmap = '%s/srcmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	modelMap = '%s/modelmap_%s_%s_%s.fits' % (OutputDirectory, tmin, tmax, sourceName)
	gtfindsrcResults = '%s/gtfindsrc_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)
	modelMapPlot = '%s/modelmap_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	modelMapPlotAnnotated = '%s/modelmapAnnotated_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	lightCurvePlot = '%s/lightcurve_%s_%s_%s.png' % (OutputDirectory, tmin, tmax, sourceName)
	lightCurveData = '%s/lightcurve_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)	
	sourceNames = '%s/ModelSources_%s_%s_%s.txt' % (OutputDirectory, tmin, tmax, sourceName)
	xmlModelBackground = '%s/Model%s_Background.xml' % (OutputDirectory,sourceName)
	xmlModelFitMaxTS = "%s/Model%s_MaxTS.xml" % (OutputDirectory, sourceName)


	# Select the photons
	print '\nSelecting the photons:'
	gtselect.run(infile=ft1file,
				outfile=filteredEvents,
				ra=ra, dec=dec, rad=radius,
				emin=emin, emax=emax,
				tmin=tmin, tmax=tmax,
				zmax=zmax, evclass=evclass)

			
	# Check to see if any photons survived the cut			 
	ft1 = pyfits.open(filteredEvents)
	if ft1['EVENTS'].header['NAXIS2'] == 0:
		print 'No photons survived the gtselect cut'
		return
		

	# Select the photons
	print '\nSelecting the good time intervals:'		
	gtmktime.run(scfile=ft2file,
				filter='DATA_QUAL>0 && LAT_CONFIG==1',
				roicut='yes',
				evfile=filteredEvents,
				outfile=filteredEventsGTI)
	

	# Check to see if any photons survived the cut			 
	ft1 = pyfits.open(filteredEventsGTI)
	if ft1['EVENTS'].header['NAXIS2'] == 0:
		print 'No photons survived the gtmktime cut'
		return				


	# Create the source model
	print '\nGenerating the xml model:'
	mymodel=srcList(Catalog2FGL,filteredEventsGTI, xmlModel)
	mymodel.makeModel(galpropModel,'gll_iem_v05',isotropicModel,'iso_source_v05',radLim=5,extDir=extendedSourcesDirectory)

	# Determine the object of interest
	if 'none' in association:
		Source = 'CandidateSource'
		AddCandidateSource(float(ra), float(dec), xmlModel)

	else:

		# Generate the source name
		Source = association
		Source = Source.replace('_','')
		Source = '_' + Source

		# Check to see if our associated source is in the xml model.  If not, add a candidate point source.
		SourceRA, SourceDec = ExtractCoordinates(xmlModel,Source)

		# Add a point source to the xml model if there is no association
		if 'NA' in SourceRA:
			association = 'none'
			Source = 'CandidateSource'
			AddCandidateSource(float(ra), float(dec), xmlModel)		
		


	# Compute the diffuse response
	print '\nComputing the diffuse response:'				
	gtdiffrsp.run(evfile=filteredEventsGTI,
				scfile=ft2file,
				srcmdl=xmlModel,
				irfs=irfs,
				evclass=evclass,
				clobber='yes')


	# Make a counts map from the event data
	print '\nCreating the counts map:'
	gtbin.run(evfile=filteredEventsGTI,
				scfile=ft2file,
				outfile=cmap,
				algorithm='CMAP',
				nxpix=160, nypix=160, binsz=binsize, coordsys='CEL',
#				nxpix=200, nypix=200, binsz=0.2, coordsys='CEL',
#				nxpix=300, nypix=300, binsz=0.2, coordsys='CEL',					
				xref=ra, yref=dec, axisrot=0, proj='AIT')


	# Make a counts map from the event data
	print '\nCreating the counts cube:'
	gtbin.run(evfile=filteredEventsGTI,
				scfile=ft2file,
				outfile=ccube,
				algorithm='CCUBE',
				nxpix=160, nypix=160, binsz=binsize, coordsys='CEL',
#				nxpix=200, nypix=200, binsz=0.2, coordsys='CEL',
				xref=ra, yref=dec, axisrot=0, proj='AIT',
				emin=emin, emax=emax, enumbins=30)
				

	# Generate the livetime cube
	print '\nGenerating the livetime cube:'		
	gtltcube.run(evfile=filteredEventsGTI,
				scfile=ft2file,
				outfile=ltcube,
				dcostheta=0.025,
				binsz=1)
			  

	# Generate an exposure map.  This is for unbinned likelihood only!	
	print '\nGenerating the exposure map:'				
	gtexpmap.run(evfile=filteredEventsGTI,
				scfile=ft2file,
				expcube=ltcube,
				outfile=expmap,
				irfs=irfs,
				srcrad=30,
				nlong=120, nlat=120, nenergies=30, coordsys='CEL')		


	# Generate an exposure map.  This is for binned likelihood only, but is also used in the srcmap generation
	print '\nGenerating the binned exposure map:'				
	gtexpcube2.run(infile=ltcube,
				cmap=ccube,
				outfile=bexpmap,
				irfs=irfs,
				nxpix=400, nypix=400, binsz=binsize, coordsys='CEL',
				xref=ra, yref=dec, axisrot=0, proj='AIT',
				emin=emin, emax=emax, nenergies=30)


	# Generate an all sky exposure map.  This is for binned likelihood only, but is also used in the srcmap generation
	print '\nGenerating the binned exposure map:'				
	gtexpcube2.run(infile=ltcube,
				cmap=ccube,
				outfile=bexpmap_allSky,
				irfs=irfs,
				nxpix=int(360/binsize), nypix=int(180/binsize), binsz=binsize, coordsys='CEL',
				xref=ra, yref=dec, axisrot=0, proj='AIT',
				emin=emin, emax=emax, nenergies=30)


	# Run the likelihood analysis
	# gtLikeFailed = False
	# try:
	# 	print '\nPerforming the likelihood fit:'		
	# 	gtlike.run(statistic='UNBINNED',
	# 				scfile=ft2file,
	# 				evfile=filteredEventsGTI,
	# 				expmap=expmap,
	# 				expcube=ltcube,
	# 				srcmdl=xmlModel,
	# 				sfile=xmlModelFit,
	# 				irfs=irfs,
	# 				optimizer=optimizer,
	# 				results=likelihoodResults,
	# 				plot='no',
	# 				save='yes')

	# except Exception, message:
	#     print traceback.format_exc()
	#     gtLikeFailed = True	


	# Setup the unbinned likelihood object
	print '\nPerforming the likelihood fit:'

	try:
	
		obs = UnbinnedObs(filteredEventsGTI,ft2file,expMap=expmap,expCube=ltcube,irfs=irfs)
	
		# Define the likelihood object
		like = UnbinnedAnalysis(obs,xmlModel,optimizer=optimizer)
			
		# Setup the likelihood parameters
		if Source == 'CandidateSource':

			# Set the flux scale
			Integral = like.par_index(Source, 'Integral')
			like[Integral].setScale(1e-3)

			# Set limits on the photon index
			Index = like.par_index(Source, 'Index')
			like[Index].setBounds(-5, -0.5)

		# Set the energy limits of the fit
		#LowerLimit = like.par_index(Source, 'LowerLimit')
		#UpperLimit = like.par_index(Source, 'UpperLimit')
		#like[LowerLimit] = emin
		#like[UpperLimit] = emax		
	
		# Perform the likelihood fit
		optObject = pyLike.NewMinuit(like.logLike)				  
		like.fit(verbosity=0,covar=True,tol=0.01,optObject=optObject)
	
		# Extract likelihood fit results
		try:
			print '\nLikelihood Results:'
			print like.model[Source]
			print "TS = %s" % like.Ts(Source)
			TSValue = like.Ts(Source)
		except:
			TSValue = 'NA'

		# Extract the photon flux
		try:
			PhotonFlux = "%.2e" % like.flux(Source, emin=100, emax=3e5)
			PhotonFluxError = "%.2e" % like.fluxError(Source, emin=100, emax=3e5)
		except:
			PhotonFlux = 'NA'
			PhotonFluxError = 'NA'

		# Extract the photon index
		try:
			Index = like.par_index(Source, 'Index')
			PhotonIndex = "%.2f" % like[Index].value()
			PhotonIndexError = "%.2f" % like[Index].error()
		except:
			PhotonIndex = 'NA'
			PhotonIndexError = 'NA'

		# Save the xml file
		like.writeXml(xmlFile=xmlModelFit)

	except Exception, message:
		print traceback.format_exc()	


	# Produce a TS map, but only if the object is unassociated
	if 'none' in association:

		# Check if the user wants to use gttsmap rather than dtsmap
		if usegttsmap == False:

			# Fix the free parameters in the best fit source model
			print "\nModifying likelihood xml model..."
			try:
				ModifySourceModel(xmlModelFit, xmlModelFitModified, RemoveSource=False)
			except Exception, message:
				print traceback.format_exc()

			# Generate a distributed TS map
			print '\nGenerating a distributed TS map:'
			try:

				from dtsmap import dtsmap
				dtsmapObject = dtsmap()
				dtsmapObject.scfile = ft2file
				dtsmapObject.evfile = filteredEventsGTI
				dtsmapObject.expmap = expmap
				dtsmapObject.expcube = ltcube
				dtsmapObject.cmap = cmap
				dtsmapObject.srcmdl = xmlModelFitModified
				dtsmapObject.irfs = irfs
				dtsmapObject.source = 'CandidateSource'
				dtsmapObject.optimizer = optimizer
				dtsmapObject.ftol = 1e-5
				dtsmapObject.nxdeg = 5
				dtsmapObject.nydeg = 5
				dtsmapObject.binsz = 0.25
				dtsmapObject.coordsys = 'CEL'
				dtsmapObject.xref = ra
				dtsmapObject.yref = dec
				dtsmapObject.proj = 'AIT'
				dtsmapObject.statistic='UNBINNED'
				dtsmapObject.outdir = dtsmapDirectory
				dtsmapObject.outfile = dtsmapFile
				dtsmapObject.maxJobs = 20
				dtsmapObject.batch = True

				tsmapResults = dtsmapObject.run(Test=False)

				print '\nMoving Results...'

				# Moving the max bin source model
				if os.path.isfile("%s/ModelSource_MaxTS.xml" % dtsmapDirectory):
					cmd = "cp %s/ModelSource_MaxTS.xml %s" % (dtsmapDirectory, xmlModelFitMaxTS)
					os.system(cmd)

				# Moving the max bin results file
				if os.path.isfile("%s/likelihoodResults_MaxTS.txt" % dtsmapDirectory):
					cmd = "cp %s/likelihoodResults_MaxTS.txt %s/likelihoodResults_%s.txt" % (dtsmapDirectory, OutputDirectory, sourceName)
					os.system(cmd)

				# Moving the max bin results file
				if os.path.isfile("%s/dtsmap_MaxTS.log" % dtsmapDirectory):
					cmd = "cp %s/dtsmap_MaxTS.log %s/dtsmap_MaxTS_%s.log" % (dtsmapDirectory, OutputDirectory, sourceName)
					os.system(cmd)

				# Move and compress the dtsmap log file
				if os.path.isfile("%s/dtsmap.log" % dtsmapDirectory):
					cmd = "cp %s/dtsmap.log %s/dtsmap_%s.log" % (dtsmapDirectory, OutputDirectory, sourceName)
					os.system(cmd)
					cmd = "gzip %s/dtsmap_%s.log" % (OutputDirectory, sourceName)
					os.system(cmd)
				
				# Remove the dtsmap directory
				cmd = "rm -R %s/" % dtsmapDirectory
				#os.system(cmd)

				print 'Done.\n'

			except Exception, message:
				print traceback.format_exc()

		# Use gttsmap
		else:

			# Create a xml model without the point source for use with gttsmap
			print "\nModifying likelihood xml model..."
			try:
				
				ModifySourceModel(xmlModelFit, xmlModelBackground, candidateSource='CandidateSource', RemoveSource=True)

			except Exception, message:
				print traceback.format_exc()

			# Find a localization using a TS Map
			print '\nGenerating a distributed TS map:'
			gttsmap.run(evfile=filteredEventsGTI,
						scfile=ft2file,
						expmap=expmap,
						expcube=ltcube,
						cmap=cmap,
						bexpmap=bexpmap,
						srcmdl=xmlModelBackground,
						outfile=tsmap,
						irfs=irfs,
						optimizer=optimizer,
						ftol=1e-8,
						nxpix=50, nypix=50, binsz=0.1, coordsys='CEL',
						xref=ra, yref=dec, proj='AIT')

			# Extract the burst localization
			print '\nExtracting TS map localization:'	
			try:

				import PlotTSMap
				MaxTS, MaxRa, MaxDec, errorRadius68, errorRadius90, errorRadius95, errorRadius99 = PlotTSMap.PlotTS(tsmap)

				# copy the background only model		
				cmd = "cp %s %s" % (xmlModelBackground, xmlModelFitMaxTS)
				os.system(cmd)

				# Create a new xml model with a point source at the max TS location
				AddCandidateSource(MaxRa, MaxDec, xmlModelFitMaxTS)

				# Perform an unbinned likelihood analysis on the max ts location to obtain the source index and flux
				print '\nPerforming the likelihood fit:'

				try:
				
					# Setup the observation object
					obs = UnbinnedObs(filteredEventsGTI,ft2file,expMap=expmap,expCube=ltcube,irfs=irfs)
				
					# Define the likelihood object
					like = UnbinnedAnalysis(obs,xmlModelFitMaxTS,optimizer=optimizer)
						
					# Setup the likelihood parameters
					if Source == 'CandidateSource':

						# Set the flux scale
						Integral = like.par_index(Source, 'Integral')
						like[Integral].setScale(1e-3)

						# Set limits on the photon index
						Index = like.par_index(Source, 'Index')
						like[Index].setBounds(-5, -0.5)

					# Setup the optimizer
					optObject = pyLike.NewMinuit(like.logLike)				

					# Perform the likelihood fit  
					like.fit(verbosity=0,covar=True,tol=0.01,optObject=optObject)
				
					# Extract likelihood fit results
					try:
						print '\nLikelihood Results:'
						print like.model[Source]
						print "TS = %s" % like.Ts(Source)
						TSValue = like.Ts(Source)
					except:
						TSValue = 'NA'

					# Extract the photon flux
					try:
						PhotonFlux = "%.2e" % like.flux(Source, emin=100, emax=3e5)
						PhotonFluxError = "%.2e" % like.fluxError(Source, emin=100, emax=3e5)
					except:
						PhotonFlux = 'NA'
						PhotonFluxError = 'NA'

					# Extract the photon index
					try:
						Index = like.par_index(Source, 'Index')
						PhotonIndex = "%.2f" % like[Index].value()
						PhotonIndexError = "%.2f" % like[Index].error()
					except:
						PhotonIndex = 'NA'
						PhotonIndexError = 'NA'

					# Save the xml file
					like.writeXml(xmlFile=xmlModelFitMaxTS)

				except Exception, message:
					print traceback.format_exc()


				# Save the likelihood results
				likelihoodResults = "%s/likelihoodResults_%s.txt" % (OutputDirectory, sourceName)
				output = open(likelihoodResults, 'w')
				output.write("Maximum TS Likelihood Results\n")
				output.write("MaxTS: %s\n" % like.Ts(Source))
				output.write("RaDec: %.3f %.3f +/- %.3f\n" % (float(SourceRA), float(SourceDec), 0))
				output.write("PhotonIndex: %s +/- %s\n" % (PhotonIndex, PhotonIndexError))
				output.write("PhotonFlux: %s +/- %s\n" % (PhotonFlux,PhotonFluxError))
				output.close()

				# Save the results to a dictionary for later use
				tsmapResults = {'MaxTS':TSValue, 'MaxRA':float(SourceRA), 'MaxDec':float(SourceDec), 'Error':0, 'Index':PhotonIndex, 'IndexError':PhotonIndexError, 'Flux':PhotonFlux,'FluxError':PhotonFluxError, 'MaxBin':None}

			
			except Exception, message:
				print traceback.format_exc()	

	else:
		
		# Save the likelihood results
		likelihoodResults = "%s/likelihoodResults_%s.txt" % (OutputDirectory, sourceName)
		output = open(likelihoodResults, 'w')
		output.write("Maximum TS Likelihood Results\n")
		output.write("MaxTS: %s\n" % like.Ts(Source))
		output.write("RaDec: %.3f %.3f +/- %.3f\n" % (float(SourceRA), float(SourceDec), 0))
		output.write("PhotonIndex: %s +/- %s\n" % (PhotonIndex, PhotonIndexError))
		output.write("PhotonFlux: %s +/- %s\n" % (PhotonFlux,PhotonFluxError))
		output.close()

		# Save the results to a dictionary for later use
		tsmapResults = {'MaxTS':TSValue, 'MaxRA':float(SourceRA), 'MaxDec':float(SourceDec), 'Error':0, 'Index':PhotonIndex, 'IndexError':PhotonIndexError, 'Flux':PhotonFlux,'FluxError':PhotonFluxError, 'MaxBin':None}


	# Find a localization using gtfindsrc
	# print '\nRunning gtfindsrc:'
	# Timeout the operation after 10 minutes
	# with timeout(seconds=18000):
	# 	try:
	# 		gtfindsrc.run(evfile=filteredEventsGTI,
	# 							   scfile=ft2file,
	# 							   outfile=gtfindsrcResults,
	# 							   irfs=irfs,
	# 							   expcube=ltcube,
	# 							   expmap=expmap,
	# 							   srcmdl=xmlModelBest,
	# 							   target='CandidateSource',
	# 							   ra=ra,dec=dec,optimizer='MINUIT',ftol=1e-2,reopt='yes',atol=0.03,toltype='ABS',
	# 							   posacc=0.001, chatter=4)
	# 	except Exception, message:
	# 	    print traceback.format_exc()
	# 	    sys.exit()

	# Find a localization using a TS Map
	# print '\nGenerating a distributed TS map:'
	# gttsmap.run(evfile=filteredEventsGTI,
	# 			scfile=ft2file,
	# 			expmap=expmap,
	# 			expcube=ltcube,
	# 			cmap=cmap,
	# 			bexpmap=bexpmap,
	# 			srcmdl=xmlModelFitModified,
	# 			outfile=tsmap,
	# 			irfs=irfs,
	# 			optimizer=optimizer,
	# 			ftol=0.01,
	# 			nxpix=40, nypix=40, binsz=0.1, coordsys='CEL',
	# 			xref=ra, yref=dec, proj='AIT')


	# Set the best fit source model to that obtined at the maximum of the ts map 
	try:
		if os.path.isfile(xmlModelFitMaxTS):
			xmlModelBest = xmlModelFitMaxTS
		else:
			xmlModelBest = xmlModelFit
	except:
			xmlModelBest = xmlModelFit


	# # Resetting the pfiles directory (?)
	# SetPfilesDirectory(PFILESDirectory)		


	# Generate a source map of the best fit likelihood model
	print '\nGenerating a source map:'
	try:
		gtsrcmaps.run(scfile=ft2file,
						cmap=ccube,
						expcube=ltcube,
						srcmdl=xmlModelBest,
						outfile=srcmap,
						bexpmap=bexpmap_allSky,
						irfs=irfs,
						emapbnds='no')
	except Exception, message:
	    print traceback.format_exc()


	# Generate model map of the best fit likelihood model	
	print '\nGenerating a model map:'	
	try:	
		gtmodel.run(srcmaps=srcmap,
						srcmdl=xmlModelBest,
						outfile=modelMap,
						irfs=irfs,
						expcube=ltcube,
						bexpmap=bexpmap_allSky)
	except Exception, message:
		print traceback.format_exc()
	

	# Plot the model map without catalog annotations
	try:

		import matplotlib
		import matplotlib.pylab as plot
		import numpy
		import pyfits
		import traceback
		from matplotlib.ticker import FuncFormatter

		fits = pyfits.open(modelMap)
	
		# Extract the wcs data from the fits file
		xReferencePixel = fits[0].header['CRPIX1']
		xReferencePixelRA = fits[0].header['CRVAL1']
		xIncrementPerPixel = fits[0].header['CDELT1']
		yReferencePixel = fits[0].header['CRPIX2']
		yReferencePixelDEC = fits[0].header['CRVAL2']
		yIncrementPerPixel = fits[0].header['CDELT2']

		# Calculate the extent of the image
		xmin = xReferencePixelRA - (xReferencePixel* xIncrementPerPixel)
		xmax = xReferencePixelRA + (xReferencePixel* xIncrementPerPixel)
		ymin = yReferencePixelDEC - (yReferencePixel* yIncrementPerPixel)
		ymax = yReferencePixelDEC + (yReferencePixel* yIncrementPerPixel)
		xRange = [xmin,xmax]
		yRange = [ymin,ymax]

		# Make sure that we don't have any ra values below zero or greater than 360, they should wrap ra instead.
		for i in range(len(xRange)):
			if xRange[i] < 0:
				xRange[i] = xRange[i] + 360.0
			if xRange[i] > 360:
				xRange[i] = xRange[i] - 360.0

		# Make sure that we don't have any dec values below or above +/- 90, they should instead wrap in both ra and dec.
		for i in range(len(yRange)):
			if yRange[i] < -90:
				yRange[i] = ((yRange[i] + 90) + 90)*-1
				#xRange[i] = xRange[i] + 180.0
			if yRange[i] > 90:
				yRange[i] = 90 - (yRange[i] - 90)
				#xRange[i] = xRange[i] + 180


		# Extract the model map
		array = fits[0].data

		# Explicitly create the figure.  This helps with forcing a specific aspect ratio later.
		fig = plot.figure()
		ax = fig.add_subplot(111)		

		# Plot the model map using a Lambert Azimuthal Equal Area Projection
		sys.path.append("/nfs/slac/g/ki/ki08/kocevski/LATBA/lib/python_rhel6-64/")
		from mpl_toolkits.basemap import Basemap

		# Create a base map on which to plot the results
		m = Basemap(height=4.3e6,width=4.3e6, projection='laea', lon_0 = float(ra)*-1, lat_0 = float(dec), resolution ='l',area_thresh=1000., celestial=True)

		# plot the image map
		m.imshow(numpy.log10(array), origin='lower', cmap=matplotlib.cm.get_cmap('seismic'))

		# Setup the map grid
		m.drawmapboundary(fill_color='#ffffff')
		m.drawparallels(numpy.arange(-90,95,5),labels=[1,0,0,0], fmt=customDecLabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
		m.drawmeridians(numpy.arange(0,365,5),labels=[0,0,0,1],fmt=customRALabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
		m.suppress_ticks = False

		# Force the aspect ratio to be 1:1
		try:
			# Force the aspect ratio to be 1:1
			forceAspect(ax,aspect=1)
		except Exception, message:
			print traceback.format_exc()

		# Annotate the plot with the original FAVA location
		x_FAVA, y_FAVA = m(float(ra), float(dec))
		#m.scatter(x_FAVA, y_FAVA, marker='+',color='skyblue', alpha=1, clip_on=True)
		#m.scatter(x_FAVA, y_FAVA, marker='+', s=75, facecolors='none', edgecolors='w')

		# Add the candidate source annotation
		try:
			# Convert the max ts ra and dec to map coordinates
			(mMaxRA, mMaxDec) = m(float(tsmapResults['MaxRA'])+0.5,float(tsmapResults['MaxDec'])+0.5)
			plot.annotate('Candidate Source', xy=(mMaxRA,mMaxDec), xytext=(-25,25), textcoords='offset points', ha='center', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='black'))
		except:
			print "TSMap results uavailable. Skipping candidate source annotation."

		# Setup the plot
		plot.xlabel('RA')
		plot.ylabel('Dec')
		plot.gca().xaxis.labelpad = 20
		plot.gca().yaxis.labelpad = 20

		print "\nSaving model map plot to: %s" % modelMapPlot
		plot.savefig(modelMapPlot, bbox_inches='tight', dpi=100)
		plot.close()

	except Exception, message:
		print traceback.format_exc()


	# Plot the model map with catalog annotations
	try:

		fig = plot.figure(figsize=(6.2,6.0))
		fig.clip_on = True
		ax = fig.add_subplot(111)	

		# Create a base map on which to plot the results
		m = Basemap(height=4.3e6,width=4.3e6, projection='laea', lon_0 = float(ra)*-1, lat_0 = float(dec), resolution ='l',area_thresh=1000., celestial=True)

		# plot the image map
		m.imshow(numpy.log10(array), origin='lower', cmap=matplotlib.cm.get_cmap('seismic'))

		# Setup the map grid
		m.drawmapboundary(fill_color='#ffffff')
		m.drawparallels(numpy.arange(-90,95,5),labels=[1,0,0,0], fmt=customDecLabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
		m.drawmeridians(numpy.arange(0,365,5),labels=[0,0,0,1],fmt=customRALabel, linewidth=0.25, color='gray', latmax=89, ax=ax)
		m.suppress_ticks = False

		# Force the aspect ratio to be 1:1
		try:
			# Force the aspect ratio to be 1:1
			forceAspect(ax,aspect=1)
		except Exception, message:
			print traceback.format_exc()

		# Setup the plot
		plot.xlabel('RA')
		plot.ylabel('Dec')
		plot.gca().xaxis.labelpad = 20
		plot.gca().yaxis.labelpad = 20

		# Extract 2FGL sources
		try:
			print '\nOpening 2FGL catalog...'
			LAT2FGL = pyfits.open('/nfs/slac/g/ki/ki08/kocevski/LATBA/Catalogs/2FGL.fits')
			print 'Done.'

			print 'Opening monitored source list..'
			MonitoredSourceList = pyfits.open('/nfs/slac/g/ki/ki08/kocevski/LATBA/Catalogs/MonitoredSourceList.fits')
			print 'Done.'
			
			print 'Opening ATel list..'
			Name_ATels = numpy.array([])
			Name_2FGL_ATels = numpy.array([])
			RA_ATels = numpy.array([])
			DEC_ATels = numpy.array([])
			ATelNumber = numpy.array([])
			ATelCatalog = '/nfs/slac/g/ki/ki08/kocevski/LATBA/Catalogs/ascd_atels_feb142013.dat'
			for line in fileinput.input([ATelCatalog]):
				if 'RA (J2000.0)' not in line:
					LineContents = line.split('|')
					Name_ATels = numpy.append(Name_ATels,LineContents[1].strip())
					Name_2FGL_ATels = numpy.append(Name_2FGL_ATels,LineContents[2].strip())
					RA_ATels = numpy.append(RA_ATels,float(LineContents[3].strip()))
					DEC_ATels = numpy.append(DEC_ATels,float(LineContents[4].strip()))
					ATelNumber = numpy.append(ATelNumber,LineContents[5].strip())
			print 'Done.'
			fileinput.close()
			
		except Exception, message:
			print message

		# 2FGL Sources
		RA_2FGL = LAT2FGL[1].data.RAJ2000
		DEC_2FGL = LAT2FGL[1].data.DEJ2000
		Flux_2FGL = LAT2FGL[1].data.Flux100_300
		SourceName_2FGL = LAT2FGL[1].data.Source_Name
		ASSOC1 = LAT2FGL[1].data.ASSOC1

		# Determine the best 2FGL source name
		BestName_2FGL = SourceName_2FGL
		#BestName_2FGL = numpy.array([])
		#for x,y in zip(ASSOC1, SourceName_2FGL):
		#		if len(x) > 0:
		#				BestName_2FGL = numpy.append(BestName_2FGL, x)
		#				pass
		#		else:
		#				BestName_2FGL = numpy.append(BestName_2FGL, y)
		#				pass
		#		pass

		# Extract the monitored sources
		Name_MonitoredSources = MonitoredSourceList[1].data['NAME']
		RA_MonitoredSources = MonitoredSourceList[1].data['RA']
		DEC_MonitoredSources = MonitoredSourceList[1].data['DEC']

		# Keep only the unique entries
		Name_MonitoredSourcesUnique, index = numpy.unique(Name_MonitoredSources,return_index=True)
		RA_MonitoredSourcesUnique = RA_MonitoredSources[index]
		DEC_MonitoredSourcesUnique = DEC_MonitoredSources[index]

		# Combine all sources
		RA_Sources = numpy.concatenate((RA_2FGL,RA_MonitoredSourcesUnique,RA_ATels))
		DEC_Sources = numpy.concatenate((DEC_2FGL,DEC_MonitoredSourcesUnique,DEC_ATels))
		Name_Sources = numpy.concatenate((BestName_2FGL,Name_MonitoredSourcesUnique,Name_ATels))

		# Remove any duplicate sources
		Names_SourcesUnique, index = numpy.unique(Name_Sources, return_index=True)
		RA_SourcesUnique = RA_Sources[index]
		DEC_SourcesUnique = DEC_Sources[index]

		# Remove individual sources
		DuplicateSources1 = ['2FGL J1745.6-2858','BL LAC', 'CGRABS J1848+3219', 'CGRABS J1849+6705','CGRABS J0211+1051','LS I+61 303','Mkn 421','PKS 0727-11','PKS 1424-41']
		DuplicateSources2 = ['PKS 1510-08','S5 1803+78','SUN','0FGLJ1641.4+3939','0235+164','S5 0716+71','GB6 J0742+5444','0827+243','PKS B0906+015','0FGL J0910.2-5044']
		DuplicateSources3 = ['1150+497','TON 599','PKS B1222+216','4C +21.35','J123939+044409','PSRB1259-63','1510-089','FERMI J1532-1321','PKS B 1622-297','4C +38.41']
		DuplicateSources4 = ['1633+382','J1717-5156','1730-130','3EG J2033+4118']
		DuplicateSources = DuplicateSources1 + DuplicateSources2 + DuplicateSources3 + DuplicateSources4
		Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique = RemoveSources(DuplicateSources, Names_SourcesUnique, RA_SourcesUnique, DEC_SourcesUnique)

		# Add the sources
		sort = [numpy.argsort(RA_SourcesUnique)]
		RA_SourcesUnique = RA_SourcesUnique[sort]
		DEC_SourcesUnique = DEC_SourcesUnique[sort] 
		Names_SourcesUnique = Names_SourcesUnique[sort]

		# Annotate the plot with the original FAVA location
		x_FAVA, y_FAVA = m(float(ra), float(dec))
		#m.scatter(x_FAVA, y_FAVA, marker='+',color='skyblue', s=50, alpha=1, clip_on=True)
		m.scatter(x_FAVA, y_FAVA, marker='+', s=75, facecolors='none', edgecolors='w')

		# Annotate the plot with the source names
		x_Sources, y_Sources = m(RA_SourcesUnique, DEC_SourcesUnique)
		m.scatter(x_Sources, y_Sources, marker='x',color='skyblue', alpha=1, clip_on=True)
		for i,j,k in zip(x_Sources, y_Sources, Names_SourcesUnique):
			plot.text(i+5e4,j,k,clip_on=True, alpha=1,size=8,color='skyblue')
			pass

		# Add the candidate source annotation
		try:
			# Convert the max ts ra and dec to map coordinates
			(mMaxRA, mMaxDec) = m(float(tsmapResults['MaxRA'])+0.5,float(tsmapResults['MaxDec'])+0.5)
			plot.annotate('Candidate Source', xy=(mMaxRA,mMaxDec), xytext=(-25,25), textcoords='offset points', ha='center', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='black'))
		except:
			print "TSMap results uavailable. Skipping candidate source annotation."


		print "\nSaving annotated model map plot to: %s" % modelMapPlotAnnotated
		plot.savefig(modelMapPlotAnnotated, bbox_inches=matplotlib.transforms.Bbox([[0.28, 0.13], [5.67, 5.5]]), dpi=100)
		plot.close()

	except Exception, message:
		print traceback.format_exc()


	# Produce a light curve
	try:
		print '\nGenerating a light curve:'
		sys.path.append("/nfs/slac/g/ki/ki08/kocevski/FAVA/python/")

		from LightCurve import LightCurve

		# Setup the light curve object
		lightCurve = LightCurve(time_start=float(tmin),time_end=float(tmax), delta_t = 86400.)
		lightCurve.ra = float(ra)
		lightCurve.dec = float(dec)
		lightCurve.evclass = evclass
		lightCurve.rad = radius
		lightCurve.emin = float(emin)
		lightCurve.emax = float(emax)
		lightCurve.irfs = irfs
		lightCurve.SC = ft2file
		lightCurve.source = Source
		lightCurve.model = xmlModelBest
		lightCurve.infile = filteredEventsGTI
		lightCurve.suffix = sourceName

		# Generate the necessary files
		lightCurve.generate_files()

		# Start the analysis
		lightCurve.runLikelihood()

		# Check for failed bins
		good = lightCurve.retCode

		# Get the range of the x-axis
		times = mdates.epoch2num(978307200+(lightCurve.times[0:-1]+lightCurve.times[1:])/2.)
		tminDate = mdates.epoch2num(978307200+float(tmin)-(86400))
		tmaxDate = mdates.epoch2num(978307200+float(tmax)+(86400))

		# Plot the results
		fig = plot.figure(figsize=(7.0,5.162))
		ax = fig.add_subplot(111)
		plot.errorbar(mdates.num2date(times),lightCurve.IntFlux,yerr=lightCurve.IntFluxErr,fmt='o')
		plot.gcf().autofmt_xdate()
		plot.ylabel('Photon Flux')
		plot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plot.ylabel('photons cm$^{-2}$ s$^{-1}$')
		plot.xlim(mdates.num2date(tminDate),mdates.num2date(tmaxDate))

		# Set the format of the x-axis
		ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b %d %Y'))

		# Save the light curve plot
		print "\nSaving light curve plot to: %s" % lightCurvePlot
		plot.savefig(lightCurvePlot, bbox_inches='tight', dpi=100)
		plot.close()

		# Save the light curve data
		print "\nSaving light curve data to: %s" % lightCurveData
		output = open(lightCurveData, 'w')
		output.write("# Day Month DayNumber Time Year PhotonFlux PhotonFluxError")
		for time, flux, error, in zip(times,lightCurve.IntFlux, lightCurve.IntFluxErr):
			date = mdates.num2date(times)
			output.write("%s %s %s\n" % (mdates.num2date(times).ctime(), flux, error))
		output.close()

	except Exception, message:
		print traceback.format_exc()

	# Save a list of model sources
	print "\nExtracting Model Sources..."	
	ExtractSources(xmlModelBest, sourceNames)

		
	# Copy all of the results from the temporary output directory to the final output directory
	# cmd = "cp -R %s/* %s/" % (OutputDirectory, OutputDirectoryFinal)
	# os.system(cmd)	

	# cmd = "rm -R %s" % OutputDirectory
	# os.system(cmd)	


	# Clean up
	print "\nCleaning up..."

	# Remove any dumped cores
	cores = glob.glob('core.*')
	if len(cores) > 0:
		cmd = "rm %s/core.*" % OutputDirectory
		os.system(cmd)	

	# Delete the pfiles directory
	if(os.path.isdir(PFILESDirectory)==True):
		cmd = "rm -R %s" % PFILESDirectory
		os.system(cmd)

	# Delete working files
	os.system("rm %s/LC_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/ft1_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/expmap_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/bexpmap_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/ltcube_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/ccube_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/cmap_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/srcmap_*%s.*" % (OutputDirectory, sourceName))
	os.system("rm %s/modelmap_*%s.fits" % (OutputDirectory, sourceName))
	os.system("rm %s/Model%s_Fit.xml" % (OutputDirectory, sourceName))
	os.system("rm %s/Model%s_Fit_Modified.xml" % (OutputDirectory, sourceName))

	print 'Done.\n'
	return
		
				

##########################################################################################
				
if __name__ == '__main__':


	if len(sys.argv) > 1:

		# Extact the keywords
		kwargs = {}
		for keyword in sys.argv:
			if '=' in keyword:
				key, value = keyword.split('=', 1)
				kwargs[key] = value

		
		if 'flarefile' in kwargs:

			# Get the user supplied flarefile name
			flarefile = kwargs['flarefile']

			# Parse the flarefile and submit the individual batch jobs for each source
			SubmitJobs(flarefile, Test=True)

		
		elif 'sourceName' in kwargs:

			# Get the user supplied keywords
			sourceName = kwargs['sourceName']
			ra = kwargs['ra']
			dec = kwargs['dec']
			tmin = kwargs['tmin']
			tmax = kwargs['tmax']
			association = kwargs['association']

			# Run the followup analysis
			UnbinnedLikelihood(sourceName, ra, dec, tmin, tmax, association=association)


	else:	

		print "usage: python RunFAVAFollowupAnalysis.py flarefile=flarefile"
	 	sys.exit()


		
