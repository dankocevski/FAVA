

import numpy
import subprocess
import matplotlib.pylab as plot
import fileinput
from math import asin, cos, radians, sin, sqrt
import matplotlib.pyplot as plot
import numpy
import pyfits
import sys

##########################################################################################	

def great_circle_distance(pnt1, pnt2, radius):
			""" Returns distance on sphere between points given as (latitude, longitude) in degrees. """
			lat1 = radians(pnt1[0])
			lat2 = radians(pnt2[0])
			dLat = lat2 - lat1
			dLon = radians(pnt2[1]) - radians(pnt1[1])
			a = sin(dLat / 2.0) ** 2 + cos(lat1) * cos(lat2) * sin(dLon / 2.0) ** 2
			return 2 * asin(min(1, sqrt(a))) * radius * 57.2957795


##########################################################################################	

def customRALabel(deg):
	if (deg == 360):
		deg = ''
	return deg

def customDecLabel(deg):
	return deg

##########################################################################################	

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

 ##########################################################################################	

def Plot(UseGtlike=False):

#	ra = 260.21
#	dec = -77.92
#	binsize = 0.25
#	LogDirectory = '/Users/kocevski/Research/Analysis/FAVA/Results/dtsmap_Source18'

	ra = 99.2
	dec = -3.92
	binsize = 0.25
	LogDirectory = '/Users/kocevski/Research/Analysis/FAVA/Results/dtsmap_Source27'

	nxdeg = 5
	nydeg = 5

	ra_min = ra - nxdeg/2.0
	ra_max = ra + nxdeg/2.0
	dec_min = dec - nydeg/2.0
	dec_max = dec + nydeg/2.0

	# Calcualate the range of ra and dec values
	ra_range = numpy.arange(ra_min,ra_max+binsize,binsize)
	dec_range = numpy.arange(dec_min,dec_max+binsize,binsize)
	xsize = len(ra_range)
	ysize = len(dec_range)

	# Make sure that we don't have any ra values below zero or greater than 360, they should wrap ra instead.
	for i in range(len(ra_range)):
		if ra_range[i] < 0:
			ra_range[i] = ra_range[i] + 360.0
		if ra_range[i] > 360:
			ra_range[i] = ra_range[i] - 360.0

	# Make sure that we don't have any dec values below or above +/- 90, they should instead wrap in both ra and dec.
	for i in range(len(dec_range)):
		if dec_range[i] < -90:
			dec_range[i] = ((dec_range[i] + 90) + 90)*-1
			ra_range[i] = ra_range[i] + 180.0
		if dec_range[i] > 90:
			dec_range[i] = 90 - (dec_range[i] - 90)
			ra_range[i] = ra_range[i] + 180

	RaDecMap = numpy.zeros(shape=(xsize, ysize))
	RaDecMap = RaDecMap.astype(str)
	RaDecPairs = {}
	binNumbers = []
	binNumber = 0
	xStep = 0
	yStep = 0
	for raStep in ra_range:
		yStep = 0
		for decStep in dec_range:
			RaDecPairs[binNumber] = [raStep,decStep]
			RaDecMap[xStep][yStep] = "%s,%s" % (raStep,decStep)
			binNumbers.append(binNumber)
			binNumber = binNumber + 1
			yStep = yStep + 1
			pass
		xStep = xStep + 1
		pass

	# Gather the results (pylikelihood)
	command = "grep 'TS = ' %s/dtsmap_bin*.log" % LogDirectory	
	process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
	lines = process.stdout.readlines()
	
	# Read in the results
	binNumbersReadIn = []
	TSs = []
	for line in lines:
		if ('No match' in line) or ('grep' in line):
			print "Error: *** No available likelihood results ***"
			print "Check your science tools setup and try again..."
			print "Exiting.\n"

			Results = {'MaxTS':None, 'MaxRA':None, 'MaxDec':None, 'Error':None, 'Index':None, 'IndexError':None, 'Flux':None,'FluxError':None, 'MaxBin':None}
			return Results
		else:
			# For use with pylikelihood results
			binNumber = int(line.split()[0].split('/')[-1].replace('.log:TS','').replace('dtsmap_bin','').strip())
			TS = float("%.2f" % float(line.split()[2].strip().replace("'","").replace(',','')))
			binNumbersReadIn.append(binNumber)
			TSs.append(TS)


	# Put the results in a dictionary for easy retrieval later
	LikelihoodResults = {key:value for key, value in zip(binNumbersReadIn,TSs)}

	# Make an matrix to store the values.  Saving any missing values as NaNs
	TSMap = numpy.zeros(shape=(xsize, ysize))
	binNumber = 0
	xStep = 0
	yStep = 0
	for raStep in ra_range:
		yStep = 0
		for decStep in dec_range:
			if binNumber in LikelihoodResults:
				TSMap[xStep][yStep] = LikelihoodResults[binNumber]
				pass
			else:
				print "Bad bin: %s" % binNumber
				TSMap[xStep][yStep] = numpy.nan
				pass
			if TSMap[xStep][yStep] < -1:
				TSMap[xStep][yStep] = numpy.nan
				pass
			binNumber = binNumber + 1
			yStep = yStep + 1
			pass
		xStep = xStep + 1
		pass


	# Loop through the results matrix and fill any missing values with interpolations from nearby pixels
	binNumber = 0
	xStep = 0
	yStep = 0
	for raStep in ra_range:
		yStep = 0
		for decStep in dec_range:
			if numpy.isnan(TSMap[xStep][yStep]) == True:

	#			try:
	#				TSMap[xStep][yStep] = numpy.mean([TSMap[xStep-1,yStep-1],TSMap[xStep-1,yStep],TSMap[xStep-1,yStep+1],TSMap[xStep,yStep-1],TSMap[xStep,yStep+1],TSMap[xStep+1,yStep-1],TSMap[xStep+1,yStep],TSMap[xStep+1,yStep+1]])
				try:					
					if numpy.isnan(TSMap[xStep-1][yStep]) == False and  numpy.isnan(TSMap[xStep+1][yStep]) == False:
						TSMap[xStep][yStep] = (TSMap[xStep-1][yStep] + TSMap[xStep+1][yStep]) / 2.0
					elif numpy.isnan(TSMap[xStep][yStep-1]) == False and  numpy.isnan(TSMap[xStep][yStep+1]) == False:
						TSMap[xStep][yStep] = (TSMap[xStep][yStep-1] + TSMap[xStep][yStep+1]) / 2.0
					else:
						TSMap[xStep][yStep] = numpy.nan
				except:
					TSMap[xStep][yStep] = numpy.nan

			yStep = yStep + 1
			pass
		xStep = xStep + 1

	# Finding the maximum TS
	MaxTS = numpy.nanmax(TSMap)
	MaxBin = numpy.nanargmax(TSMap)

	# Check to see if a maximum couldn't be found.
	if numpy.isnan(MaxBin) == True:
		print "\nAnalysis Complete."
		print "Maximum TS: %s" % 'None'
		print "Coordinates: RA = %s, Dec = %s" % ('NA', 'NA')
		print "*** Unable to locate maximum TS ***"

		print '\nCleaning up...'
		# Cat the gtlike log files
		cmd = "cat %s/likelihoodResults_bin*.txt > dtsmap_LikelihoodResults.txt" % self.outdir
#			print cmd
		os.system(cmd)

		# Cat the log files
		cmd = "cat %s/dtsmap_bin*.log > dtsmap_LikelihoodFits.log" % self.outdir
#			print cmd
		os.system(cmd)

		# Erase the individual xml files
		cmd = "rm %s/ModelSource_bin*.xml" % self.outdir
#			print cmd
		os.system(cmd)

		print 'Done.'

		Results = {'MaxTS':None, 'MaxRA':None, 'MaxDec':None, 'Error':None, 'Index':None, 'IndexError':None, 'Flux':None,'FluxError':None, 'MaxBin':None}
		return Results

	else:

		MaxRa = RaDecPairs[MaxBin][0]
		MaxDec = RaDecPairs[MaxBin][1]

	# Define default spectral parameters
	index = 'NA'
	indexError = 'NA'
	flux = 'NA'
	fluxError = 'NA'

	if UseGtlike == True:

		# Extract the fit parameters for the bin with the maximum TS (gtlike)
		MaxBinFile = "%s/dtsmap_bin%s.log" % (LogDirectory, MaxBin)
		for line in fileinput.input([MaxBinFile]):
			if 'Index:' in line: 
				lineContents = line.split()	
				index = lineContents[1]
				indexError = lineContents[3]
			if 'Flux:' in line: 
				lineContents = line.split()	
				flux = lineContents[1]
				fluxError = lineContents[3]
				break
		fileinput.close()

	else:

		# Extract the spectral fit parameters for the bin with the maximum TS (pyLikelihood)
		MaxBinFile = "%s/dtsmap_bin%s.log" % (LogDirectory, MaxBin)
		for line in fileinput.input([MaxBinFile]):
			if 'Flux =' in line: 
				lineContents = line.split()	
				flux = float(lineContents[2])
				fluxError = float(lineContents[4])
			if 'Index =' in line: 
				lineContents = line.split()	
				index = float(lineContents[2])
				indexError = float(lineContents[4])
		fileinput.close()


	# Rotate and flip the matrix in order to it to match ra and dec ordering conventions
	TSMapRotated = numpy.rot90(TSMap)
	TSMapFlippedUp = numpy.flipud(TSMapRotated)
	TSMapFlippedLeft2Right = numpy.fliplr(TSMapFlippedUp)
	MaxRa = RaDecPairs[numpy.nanargmax(TSMap)][0]
	MaxDec = RaDecPairs[numpy.nanargmax(TSMap)][1]

	# Import the basemap module
	sys.path.append("/nfs/slac/g/ki/ki08/kocevski/LATBA/lib/python_rhel6-64/")
	from mpl_toolkits.basemap import Basemap

	# Create the figure
	fig = plot.figure()
	ax = fig.add_subplot(111)

	# Create a base map on which to plot the results
	m = Basemap(height=5.5e5,width=5.5e5, projection='laea', lon_0 = ra*-1, lat_0 = dec, resolution ='l', area_thresh=1000., celestial=True )

	# Set the plot limits (in map coordinates)
	xMin, yMin = m(ra_range[0], dec_range[0])
	xMax, yMax = m(ra_range[-1], dec_range[-1])
	m.lonmin = xMin
	m.lonmax = xMax
	m.latmin = yMin
	m.latmax = yMax

	# Plot the matrix as an image
	extent=[xMax, xMin, yMin, yMax]
	#m.imshow(TSMapFlippedLeft2Right, origin='lower', extent=extent)
	m.imshow(TSMapFlippedLeft2Right, origin='lower')

	# Setup the map grid
	m.drawmapboundary(fill_color='#ffffff')
	m.drawparallels(numpy.arange(181)-90,labels=[1,0,0,0], fmt=customDecLabel, color='grey', linewidth=0.25)
	m.drawmeridians(numpy.arange(0,360,1),labels=[0,0,0,1], fmt=customRALabel, color='grey', linewidth=0.25)
	m.suppress_ticks = False
	m.fix_aspect = False

	# Force the aspect ratio to be 1:1
	try:
		forceAspect(ax,aspect=1)
		#extent=[xMax, xMin, yMin, yMax]
		#ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2])))
	except Exception, message:
		print traceback.format_exc()

	# Add a color bar
	colorBarObject = plot.colorbar()
	colorBarObject.set_label('TS')

	# Setup the plot
	plot.xlabel('RA')
	plot.ylabel('Dec')
	plot.gca().xaxis.labelpad = 20
	plot.gca().yaxis.labelpad = 20

	# Get the ra and dec of the max TS in map coordinates
	mMaxRa, mMaxDec = m(MaxRa, MaxDec)
	mRa, mDec = m(ra, dec)


	# Define the 68%, 95%, and 99% confidence levels
	MaxTS_68CL = (MaxTS-1)
	MaxTS_95CL = (MaxTS-4)
	MaxTS_99CL = (MaxTS-9)

	# Add contour lines
	try:
		contourObject = plot.contour(TSMapFlippedLeft2Right, levels=[MaxTS_68CL,MaxTS_95CL,MaxTS_99CL], origin='lower', extent=[ra_range[-1], ra_range[0], dec_range[0], dec_range[-1]])

		# Extract the 68% contour interval data
		contourLine68 = numpy.array(contourObject.allsegs[0])[0]
		contourLine68_ra = contourLine68[:,0]
		contourLine68_dec = contourLine68[:,1]

		# Extract the 95% contour interval data	
		contourLine95 = numpy.array(contourObject.allsegs[1])[0]
		contourLine95_ra = contourLine95[:,0]
		contourLine95_dec = contourLine95[:,1]

		# Extract the 99% contour interval data	
		contourLine99 = numpy.array(contourObject.allsegs[2])[0]
		contourLine99_ra = contourLine99[:,0]
		contourLine99_dec = contourLine99[:,1]		

		# Convert the contour line data into map coordinates
		mContourLine68_ra, mContourLine68_dec = m(contourLine68_ra, contourLine68_dec)
		mContourLine95_ra, mContourLine95_dec = m(contourLine95_ra, contourLine95_dec)
		mContourLine99_ra, mContourLine99_dec = m(contourLine99_ra, contourLine99_dec)

		# Plot the contour lines in map coordinates
		plot.plot(mContourLine68_ra, mContourLine68_dec, linestyle='dotted', color='black')	
		plot.plot(mContourLine95_ra, mContourLine95_dec, linestyle='dotted', color='black')	
		plot.plot(mContourLine99_ra, mContourLine99_dec, linestyle='dotted', color='black')	

		# Find the spherical distance from the contour interval data to the starting ra and dec
		sphericalDistances68 = []
		for contour_ra, contour_dec in zip(contourLine68_ra,contourLine68_dec):
			sphericalDistance = great_circle_distance([MaxRa,MaxDec],[contour_ra,contour_dec],1)
			sphericalDistances68.append(sphericalDistance)
			pass

		sphericalDistances95 = []
		for contour_ra, contour_dec in zip(contourLine95_ra,contourLine95_dec):
			sphericalDistance = great_circle_distance([MaxRa,MaxDec],[contour_ra,contour_dec],1)
			sphericalDistances95.append(sphericalDistance)
			pass

		sphericalDistances99 = []
		for contour_ra, contour_dec in zip(contourLine99_ra,contourLine99_dec):
			sphericalDistance = great_circle_distance([MaxRa,MaxDec],[contour_ra,contour_dec],1)
			sphericalDistances99.append(sphericalDistance)
			pass		

		# Find the median of the spherical distances.  This will be considered as the error in the localization
		sphericalDistances68 = numpy.array(sphericalDistances68)
		errorRadius68 = numpy.median(sphericalDistances68)
		sphericalDistances95 = numpy.array(sphericalDistances95)
		errorRadius95 = numpy.median(sphericalDistances95)
		sphericalDistances99 = numpy.array(sphericalDistances99)
		errorRadius99 = numpy.median(sphericalDistances99)

		print ra_range[-1],dec_range[0]

		# Annotate the plot
		m.scatter(mMaxRa, mMaxDec, marker='x', s=75, facecolors='none', edgecolors='w')
		m.scatter(mRa, mDec, marker='+', s=75, facecolors='none', edgecolors='w')

		plot.annotate("Max TS = %s\nRA = %s, Dec = %s\nError = +/-%0.3f (99%% CL)" % (MaxTS, MaxRa, MaxDec, errorRadius99), xy=(0,0), xycoords='axes points',  xytext=(10,10), textcoords='offset points', ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3))#, arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='b'))

		print "\nAnalysis Complete."
		print "Maximum TS: %s" % MaxTS
		print "Coordinates: RA = %s, Dec = %s" % (MaxRa, MaxDec)
	#	print "Error Radius (68%%): %.4f deg" % errorRadius68
		print "Error Radius (95%%): %.4f deg" % errorRadius95
		print "Error Radius (99%%): %.4f deg" % errorRadius99

		print ''
		print 'Best Fit Parameters:'
		print "Index = %.2f +/- %.2f" % (index, indexError)
		print "Flux = %.2e +/- %.2e" % (flux, fluxError)

	except:

		# Set the error radius to None
		errorRadius95 = None

		# Annotate the plot
		m.scatter(mMaxRa, mMaxDec, marker='x', s=75, facecolors='none', edgecolors='w')
		plot.annotate("Max TS = %s\nRA=%s, Dec=%s\nError = NA" % (MaxTS, MaxRa, MaxDec), xy=(0,0), xycoords='axes points', xytext=(10,10), textcoords='offset points', ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3))#, arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='b'))

		print "\nAnalysis Complete."
		print "Maximum TS: %s" % MaxTS
		print "Coordinates: RA = %s, Dec = %s" % (MaxRa, MaxDec)
	#	print "Error Radius (68%%): NA"
		print "Error Radius (95%%): NA"
		print "Error Radius (99%%): NA"

		print ''
		print 'Best Fit Parameters:'
		print "Index = %s +/- %s" % (index, indexError)
		print "Flux = %s +/- %s" % (flux, fluxError)	



	plot.show()
#	plot.savefig(lightCurvePlot, bbox_inches='tight', dpi=100)

	return




