#!/usr/bin/env python

import os
import subprocess
import pyfits
import fileinput
import sys
import numpy
import fileinput
import os
from matplotlib import dates as mdates
import matplotlib.colors as mcolors
import traceback
from math import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plot
import traceback


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

def run(significanceMap, significanceMapHE, flareFile):

	# Plot the model map using a Lambert Azimuthal Equal Area Projection
	sys.path.append("/nfs/slac/g/ki/ki08/kocevski/LATBA/lib/python_rhel6-64/")
	from mpl_toolkits.basemap import Basemap

	# Open the signifiance map
	sigmap = pyfits.open(significanceMap)
	sigmapHE = pyfits.open(significanceMapHE)

	# Extract the data
	map = sigmap[0].data
	mapHE = sigmapHE[0].data

	# Remove empty portions of the map
	#mapModified = map[18:-19,36:-37]
	mapModified = map[11:-12,22:-23]
	mapHEModified = mapHE[11:-12,22:-23]

	# Flip the map
	mapModified = numpy.flipud(mapModified)
	mapHEModified = numpy.flipud(mapHEModified)

	# Define a custom color map
	c = mcolors.ColorConverter().to_rgb
	cmap_FermiSkyMap = make_colormap([c('#00004c'), c('#1111ff'), 0.35,c('#1111ff'),c('#ff1a1a'),0.60,c('#ff1a1a'),c('#fffc17')])
	cmap_BlackBlueWhite = make_colormap([c('#000000'), c('#000a78'), 0.20, c('#000a78'), c('#0b5fef'), 0.5, c('#0b5fef'),c('#19b6ff'), 0.70, c('#19b6ff'),c('#ffffff')])

	# Open the figure and set the projection to aitoff
	figure = plot.figure()
	ax = figure.add_subplot(111, projection='aitoff', adjustable='box', frame_on=True)

	# Plot the signifiance map
	plot.imshow(mapModified, extent=[-pi,pi,-pi/2,pi/2], cmap=cmap_FermiSkyMap, vmin=-4)

	# Add a color bar
	colorBarObject = plot.colorbar(orientation='horizontal')
	colorBarObject.set_label('Significance', size=10)
	colorBarObject.ax.tick_params(labelsize=10) 

	# Force the aspect ratio to be 2:1
	ax.set_aspect(0.5)

	# Remove the axes labels
	ax.axes.xaxis.set_ticklabels([])
	ax.axes.yaxis.set_ticklabels([])

	# Add a grid
	ax.grid(True, color='w')

	# Read in the flare information
	nums, ras, decs, lbin, bbin, gall, galb, tmins, tmaxs, sigma, avnev, nev, he_nev, he_avnev, he_sigma, sundist, varindex, favasrc, fglassoc, assoc = numpy.loadtxt(flareFile, unpack=True, comments='#', dtype = str)

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
	plot.scatter(gall, galb, s=20, marker='+', linewidth=0.75, color='w')

	# Label the flares
	for num, x, y in zip(nums,gall,galb):
		plot.annotate(num,xy=(x,y), color='w', xytext=(5,5), textcoords='offset points')

	# Make a copy of the original significance map
	outdir = significanceMap[0:significanceMap.rfind('/')+1]
	significanceMapPlot = outdir + 'flares_allsky_GAL_0.5dpb_AIT_100-300000MeV_P7SOURCE_V6MC.png'

	# Make a copy of the original significance map before overwritting it
	significanceMapPlotOriginal = significanceMapPlot.replace('.png','.original.png')
	os.system("cp %s %s" % (significanceMapPlot, significanceMapPlotOriginal))

	# Save the plot
	print '\nSaving low energy map: %s' % significanceMapPlot
	plot.savefig(significanceMapPlot, bbox_inches='tight', dpi=100)
	plot.close()

	# Open the figure and set the projection to aitoff
	figure = plot.figure()
	ax = figure.add_subplot(111, projection='aitoff', adjustable='box', frame_on=True)

	# Plot the signifiance map
	plot.imshow(mapHEModified, extent=[-pi,pi,-pi/2,pi/2], cmap=cmap_FermiSkyMap, vmin=-4)

	# Add a color bar
	colorBarObject = plot.colorbar(orientation='horizontal')
	colorBarObject.set_label('Significance', size=10)
	colorBarObject.ax.tick_params(labelsize=10) 

	# Force the aspect ratio to be 2:1
	ax.set_aspect(0.5)

	# Remove the axes labels
	ax.axes.xaxis.set_ticklabels([])
	ax.axes.yaxis.set_ticklabels([])

	# Add a grid
	ax.grid(True, color='w')

	# Read in the flare information
	nums, ras, decs, lbin, bbin, gall, galb, tmins, tmaxs, sigma, avnev, nev, he_nev, he_avnev, he_sigma, sundist, varindex, favasrc, fglassoc, assoc = numpy.loadtxt(flareFile, unpack=True, comments='#', dtype = str)

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
	plot.scatter(gall, galb, s=20, marker='+', linewidth=0.75, color='w')

	# Label the flares
	for num, x, y in zip(nums,gall,galb):
		plot.annotate(num,xy=(x,y), color='w', xytext=(5,5), textcoords='offset points')

	# Construct the output filename
	outdir = significanceMapHE[0:significanceMapHE.rfind('/')+1]
	significanceMapHEPlot = outdir + 'flares_allsky_GAL_0.5dpb_AIT_100-300000MeV_P7SOURCE_V6MC_he.png'

	# Make a copy of the original significance map before overwritting it
	significanceMapHEPlotOriginal = significanceMapHEPlot.replace('.png','.original.png')
	os.system("cp %s %s" % (significanceMapPlot, significanceMapPlotOriginal))

	# Save the plot
	print 'Saving high energy map: %s' % significanceMapHEPlot
	plot.savefig(significanceMapHEPlot, bbox_inches='tight', dpi=100)
	plot.close()


##########################################################################################
				
if __name__ == '__main__':


	if len(sys.argv) > 1:

		# Extact the keywords
		kwargs = {}
		for keyword in sys.argv:
			if '=' in keyword:
				key, value = keyword.split('=', 1)
				kwargs[key] = value

		try:	
			# Get the user supplied keywords
			significanceMap = kwargs['sigmap']
			significanceMapHE = kwargs['sigmapHE']
			flareFile = kwargs['flarefile']


			# Run the followup analysis
			run(significanceMap, significanceMapHE, flareFile)

		except Exception, message:
			print traceback.format_exc()

			print "usage: python PlotSignificanceMap.py sigmap=sigmap, sigmapHE=sigmapHE, flarefile=flareFile"
		 	sys.exit()			

	else:	

		print "usage: python PlotSignificanceMap.py sigmap=sigmap, sigmapHE=sigmapHE, flarefile=flareFile"
	 	sys.exit()

