#!/usr/bin/env python

import os
import pyfits
import fileinput
import sys
import numpy
import fileinput
import os
import matplotlib.colors as mcolors
import traceback
from math import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plot
import glob

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

def run():

	# Define some default values
	runDirectory = '/nfs/slac/g/ki/ki06/buehler/nobackup/favaonline/data/weekly/P7SOURCE_V6MC/maps/'
	flarefile = 'flares_allsky_GAL_0.5dpb_AIT_100-300000MeV_P7SOURCE_V6MC.ff'
	diffuseModel = '/nfs/slac/g/ki/ki08/kocevski/FAVA/diffuseModels/gll_iem_v05.fit'

	# Open the signifiance map
	fits = pyfits.open(diffuseModel)

	# Extract the data
	diffuseMap = fits[0].data[0]

	# Flip the map
	diffuseMap = numpy.flipud(diffuseMap)

	# Define a custom color map
	c = mcolors.ColorConverter().to_rgb
	cmap_FermiSkyMap = make_colormap([c('#00004c'), c('#1111ff'), 0.35,c('#1111ff'),c('#ff1a1a'),0.60,c('#ff1a1a'),c('#fffc17')])
	cmap_BlackBlueWhite = make_colormap([c('#000000'), c('#000a78'), 0.20, c('#000a78'), c('#0b5fef'), 0.5, c('#0b5fef'),c('#19b6ff'), 0.70, c('#19b6ff'),c('#ffffff')])

	# Open the figure and set the projection to aitoff
	figure = plot.figure(figsize=(10, 8))
	ax = figure.add_subplot(111, projection='aitoff', adjustable='box', frame_on=True)

	# Plot the signifiance map
	plot.imshow(numpy.log10(diffuseMap), extent=[-pi,pi,-pi/2,pi/2], cmap=cmap_FermiSkyMap)

	# Add a color bar
	#colorBarObject = plot.colorbar(orientation='horizontal')
	#colorBarObject.set_label('Significance')

	# Force the aspect ratio to be 2:1
	ax.set_aspect(0.5)

	# Remove the axes labels
	ax.axes.xaxis.set_ticklabels([])
	ax.axes.yaxis.set_ticklabels([])

	# Add a grid
	ax.grid(True, color='gray')

	# Get the directory names for all FAVA runs
	directories = glob.glob(runDirectory + '*')

	# Record the number of sources plotted
	numberOfSources = 0
	numberOfPositiveSources = 0
	numberOfNegativeSources = 0
	numberOfWeeks = 0
	NumberOfUnassociatedSources = 0
	NumberOfAssociatedSources = 0

	# Loop through the directories and read the flare information
	for directory in directories:

		# Read the flare file
		print 'Reading flare file: %s' % directory + '/' + flarefile

		try:
			nums, ras, decs, lbin, bbin, gall, galb, tmins, tmaxs, sigma, avnev, nev, he_nev, he_avnev, he_sigma, sundist, varindex, favasrc, fglassoc, assoc = numpy.loadtxt(directory + '/' + flarefile, unpack=True, comments='#', dtype = str)

			# Seperate the negative and positive flares
			positiveFlares = numpy.where( (sigma.astype(float) > 5) | (he_sigma.astype(float) > 5))[0]
			negativeFlares = numpy.where( (sigma.astype(float) < -3) | (he_sigma.astype(float) < -3))[0]

			# Keep only the positive 3 sigma flares
			gall = gall[positiveFlares]
			galb  = galb[positiveFlares]
			sigmaPositive = sigma[positiveFlares]
			sigmaPositive = sigmaPositive.astype(float)
			fglassoc = fglassoc[positiveFlares]
			
			# Make a custom color array
			color = numpy.arange(len(gall)).astype('S25')
			associated = numpy.where(fglassoc != 'none')
			unassociated = numpy.where(fglassoc == 'none')
			
			color[associated] = 'skyblue'
			color[unassociated] = '#e14169'


			# Covnert the data to floats
			gall = gall.astype(float)
			galb = galb.astype(float)

			# Convert the galactic latitude to +/- 180 scale
			for i in range(len(gall)):
				if (gall[i] <= 180):
					gall[i] = -1 * gall[i]
				if (gall[i] > 180):
					gall[i] = -1 * (gall[i] - 360)

			# Convert from degrees to radians
			Deg2Rad = 0.0174532925
			StretchFactor = 1.00
			gall = gall * Deg2Rad * StretchFactor
			galb = galb * Deg2Rad * StretchFactor


			# Plot the flares
			#plot.scatter(gall, galb, s=8, marker='o', linewidth=0.25, edgecolor='black', color='white', alpha=0.85)
			#plot.scatter(gall, galb, s=10, marker='+', linewidth=0.75, color='w')
			plot.scatter(gall, galb, c=color, s=8, marker='o', linewidth=0.25, edgecolor='black', alpha=1.0)
			#plot.scatter(gall, galb, c=sigmaPositive, s=8, marker='o', linewidth=0.25, edgecolor='black', alpha=0.85, cmap=matplotlib.cm.get_cmap('Reds'))

			# Save the total number of sources plotted
			numberOfPositiveSources = numberOfPositiveSources + len(positiveFlares)
			numberOfNegativeSources = numberOfNegativeSources + len(negativeFlares)	
			numberOfSources = numberOfSources + len(positiveFlares) + len(negativeFlares)
			NumberOfUnassociatedSources = NumberOfUnassociatedSources + len(numpy.where(fglassoc == 'none')[0])
			NumberOfAssociatedSources = NumberOfAssociatedSources + len(numpy.where(fglassoc != 'none')[0])

			numberOfWeeks = numberOfWeeks + 1

		except Exception, message:
			print traceback.format_exc()
			continue

	# Annotate the plot
	plot.annotate("Weeks Analyzed: %s\nFAVA Detections (>5$\sigma$) : %s" % (numberOfWeeks,numberOfPositiveSources), xy=(0,-290), xycoords='axes points', xytext=(0,0), textcoords='offset points', ha='left', va='bottom', size=10)
	#plot.annotate("Positive Flares: %s\nQuiescent Sources: %s" % (numberOfPositiveSources, numberOfNegativeSources), xy=(550,-285), xycoords='axes points', xytext=(0,0), textcoords='offset points', ha='right', va='bottom', size=10)
	plot.annotate("Associated: %s\nUnassociated: %s" % (NumberOfAssociatedSources,NumberOfUnassociatedSources ), xy=(550,-285), xycoords='axes points', xytext=(0,0), textcoords='offset points', ha='right', va='bottom', size=10)


	# Print a flare summary
	print 'Total number of detections: %s.' % numberOfSources	
	print '\nTotal number of positive detections: %s.' % numberOfPositiveSources
	print 'Total number of negative detections: %s.' % numberOfNegativeSources
	print 'Total number of associated sources: %s.' % NumberOfAssociatedSources
	print 'Total number of unassociated sources: %s.' % NumberOfUnassociatedSources



	# Display the map
	FAVAFlareSummaryPlot = runDirectory + 'FAVAFlareSummary.png'
	print '\nSaving FAVA flares map: %s' % FAVAFlareSummaryPlot
	plot.savefig(FAVAFlareSummaryPlot, bbox_inches='tight', dpi=100)
	plot.close()



##########################################################################################
				
if __name__ == '__main__':


	# if len(sys.argv) > 1:

	# 	# Extact the keywords
	# 	kwargs = {}
	# 	for keyword in sys.argv:
	# 		if '=' in keyword:
	# 			key, value = keyword.split('=', 1)
	# 			kwargs[key] = value

	# 	try:	
	# 		# Get the user supplied keywords
	# 		significanceMap = kwargs['sigmap']
	# 		significanceMapHE = kwargs['sigmapHE']
	# 		flareFile = kwargs['flarefile']


	# 		# Run the followup analysis
	# 		run(significanceMap, significanceMapHE, flareFile)

	# 	except:

	# 		print "usage: python PlotFAVAFlares.py flarefile=flarefile"
	# 	 	sys.exit()			

	# else:	

	# 	print "usage: python PlotFAVAFlares.py flarefile=flarefile"
	#  	sys.exit()

	run()

