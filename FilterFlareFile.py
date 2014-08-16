#!/usr/bin/env python

import os
import subprocess
import pyfits
import fileinput
import sys
import numpy
import fileinput

##########################################################################################

def run(flarefile, Save=True):

	# Read in the flare information
	nums, ras, decs, lbin, bbin, gall, galb, tmins, tmaxs, sigma, avnev, nev, he_nev, he_avnev, he_sigma, sundist, varindex, favasrc, fglassoc, assoc = numpy.loadtxt(flarefile, unpack=True, comments='#', dtype = str)

	# Find the negative and positive flares
	positiveFlares = numpy.where( (sigma.astype(float) > 3) | (he_sigma.astype(float) > 3))[0]
	negativeFlares = numpy.where( (sigma.astype(float) < -3) | (he_sigma.astype(float) < -3))[0]

	# Find the unassociated positive flares
	unassociatedFlares = numpy.where(fglassoc[positiveFlares] == 'none')[0]

	# Find the associated positive flares
	associatedFlares = numpy.where(fglassoc[positiveFlares] != 'none')[0]

	# Renumber the bursts
	sourceNumber = 1

	# Read in the lines from the original flare file
	infile = open(flarefile,'r')
	lines = infile.readlines()

	if Save == True:

		# Open a new file
		flarefileFiltered = flarefile.replace('.ff','_Filtered.ff')
		outfile = open(flarefileFiltered,'w')

		# Print the header of the original flare file
		numberOfLines = 0
		for line in lines:
			outfile.write(line)
			numberOfLines = numberOfLines + 1
			if numberOfLines == 3:
				break			

		# Write out the new lines
		outfile.write("#Unassociated Flares\n")
		for numsi, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ in zip(nums[positiveFlares][unassociatedFlares], ras[positiveFlares][unassociatedFlares], decs[positiveFlares][unassociatedFlares], lbin[positiveFlares][unassociatedFlares], bbin[positiveFlares][unassociatedFlares], gall[positiveFlares][unassociatedFlares], galb[positiveFlares][unassociatedFlares], tmins[positiveFlares][unassociatedFlares], tmaxs[positiveFlares][unassociatedFlares], sigma[positiveFlares][unassociatedFlares], avnev[positiveFlares][unassociatedFlares], nev[positiveFlares][unassociatedFlares], he_nev[positiveFlares][unassociatedFlares], he_avnev[positiveFlares][unassociatedFlares], he_sigma[positiveFlares][unassociatedFlares], sundist[positiveFlares][unassociatedFlares], varindex[positiveFlares][unassociatedFlares], favasrc[positiveFlares][unassociatedFlares], fglassoc[positiveFlares][unassociatedFlares], assoc[positiveFlares][unassociatedFlares]):
			outfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % (sourceNumber, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ))
			sourceNumber = sourceNumber + 1

		outfile.write("#Associated Flares\n")
		for numsi, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ in zip(nums[positiveFlares][associatedFlares], ras[positiveFlares][associatedFlares], decs[positiveFlares][associatedFlares], lbin[positiveFlares][associatedFlares], bbin[positiveFlares][associatedFlares], gall[positiveFlares][associatedFlares], galb[positiveFlares][associatedFlares], tmins[positiveFlares][associatedFlares], tmaxs[positiveFlares][associatedFlares], sigma[positiveFlares][associatedFlares], avnev[positiveFlares][associatedFlares], nev[positiveFlares][associatedFlares], he_nev[positiveFlares][associatedFlares], he_avnev[positiveFlares][associatedFlares], he_sigma[positiveFlares][associatedFlares], sundist[positiveFlares][associatedFlares], varindex[positiveFlares][associatedFlares], favasrc[positiveFlares][associatedFlares], fglassoc[positiveFlares][associatedFlares], assoc[positiveFlares][associatedFlares]):
			outfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % (sourceNumber, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ))
			sourceNumber = sourceNumber + 1

		outfile.write("#Quiescent Sources\n")
		for numsi, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ in zip(nums[negativeFlares], ras[negativeFlares], decs[negativeFlares], lbin[negativeFlares], bbin[negativeFlares], gall[negativeFlares], galb[negativeFlares], tmins[negativeFlares], tmaxs[negativeFlares], sigma[negativeFlares], avnev[negativeFlares], nev[negativeFlares], he_nev[negativeFlares], he_avnev[negativeFlares], he_sigma[negativeFlares], sundist[negativeFlares], varindex[negativeFlares], favasrc[negativeFlares], fglassoc[negativeFlares], assoc[negativeFlares]):
			outfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % (sourceNumber, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ))
			sourceNumber = sourceNumber + 1

		# Close the file
		print "Saving filtered flare file: %s" % flarefileFiltered
		outfile.close()	

	else:

		# Print the header of the original flare file
		numberOfLines = 0
		for line in lines:
			print line
			numberOfLines = numberOfLines + 1
			if numberOfLines == 3:
				break	

		print "Unassociated Flares"
		for numsi, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ in zip(nums[positiveFlares][unassociatedFlares], ras[positiveFlares][unassociatedFlares], decs[positiveFlares][unassociatedFlares], lbin[positiveFlares][unassociatedFlares], bbin[positiveFlares][unassociatedFlares], gall[positiveFlares][unassociatedFlares], galb[positiveFlares][unassociatedFlares], tmins[positiveFlares][unassociatedFlares], tmaxs[positiveFlares][unassociatedFlares], sigma[positiveFlares][unassociatedFlares], avnev[positiveFlares][unassociatedFlares], nev[positiveFlares][unassociatedFlares], he_nev[positiveFlares][unassociatedFlares], he_avnev[positiveFlares][unassociatedFlares], he_sigma[positiveFlares][unassociatedFlares], sundist[positiveFlares][unassociatedFlares], varindex[positiveFlares][unassociatedFlares], favasrc[positiveFlares][unassociatedFlares], fglassoc[positiveFlares][unassociatedFlares], assoc[positiveFlares][unassociatedFlares]):
			line = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (sourceNumber, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ)
			sourceNumber = sourceNumber + 1

		print "\nAssociated Flares"
		for numsi, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ in zip(nums[positiveFlares][associatedFlares], ras[positiveFlares][associatedFlares], decs[positiveFlares][associatedFlares], lbin[positiveFlares][associatedFlares], bbin[positiveFlares][associatedFlares], gall[positiveFlares][associatedFlares], galb[positiveFlares][associatedFlares], tmins[positiveFlares][associatedFlares], tmaxs[positiveFlares][associatedFlares], sigma[positiveFlares][associatedFlares], avnev[positiveFlares][associatedFlares], nev[positiveFlares][associatedFlares], he_nev[positiveFlares][associatedFlares], he_avnev[positiveFlares][associatedFlares], he_sigma[positiveFlares][associatedFlares], sundist[positiveFlares][associatedFlares], varindex[positiveFlares][associatedFlares], favasrc[positiveFlares][associatedFlares], fglassoc[positiveFlares][associatedFlares], assoc[positiveFlares][associatedFlares]):
			print "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (sourceNumber, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ)
			sourceNumber = sourceNumber + 1

		print "\nQuiescent Sources"
		for numsi, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ in zip(nums[negativeFlares], ras[negativeFlares], decs[negativeFlares], lbin[negativeFlares], bbin[negativeFlares], gall[negativeFlares], galb[negativeFlares], tmins[negativeFlares], tmaxs[negativeFlares], sigma[negativeFlares], avnev[negativeFlares], nev[negativeFlares], he_nev[negativeFlares], he_avnev[negativeFlares], he_sigma[negativeFlares], sundist[negativeFlares], varindex[negativeFlares], favasrc[negativeFlares], fglassoc[negativeFlares], assoc[negativeFlares]):
			print "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (sourceNumber, rasi, decsi, lbini, bbini, galli, galbi, tminsi, tmaxsi, sigmai, avnevi, nevi, he_nevi, he_avnevi, he_sigmai, sundisti, varindexi, favasrci, fglassoci, associ)
			sourceNumber = sourceNumber + 1





##########################################################################################
				
if __name__ == '__main__':


	# Check to see if any arguments were passed
	if(len(sys.argv) > 0):

		# Extract the first arugument
		flarefile = sys.argv[1]
		
		# Run the followup analysis
		run(flarefile)

	# If no arguments were passed, tell the user how to use the script		
	else:
		print "usage: python FilterFlareFile.py flarefile=flarefile"
		sys.exit(0)




