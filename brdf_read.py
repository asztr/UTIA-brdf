"""
**************************************************************************
\file    read_brdf.py
\author  Alejandro Sztrajman
\date    March 2018
\version 1.00

Script for the UTIA BRDF data loading and value interpolation

################################################################
run using:
   python read_brdf.py
   (BRDF bin files are assumed to be in subdir 'data/')
################################################################

Based on the code by Jiri Filip.
******************************************************************************
"""

import imageio
import sys
import os
import math
import numpy as np
import pdb
import struct
import matplotlib.pyplot as plt

class BRDF:
	def __init__(self, step_t=15, step_p=7.5, nti=6, ntv=6, planes=3): #FIXME: make all inputs the same units (radians or degrees)
		"""Initialization of BRDF measurement parameters
		:param step_t (float): step of elevation angles theta_i, theta_v
		:param step_p (float): step of azimuthal angles phi_i, phi_v
		:param nti (int): number of theta_i directions
		:param ntv (int): number of theta_v directions
		"""
		self.step_t = step_t
		self.step_p = step_p
		self.nti = nti
		self.ntv = ntv
		self.npi = int(360.0/step_p) #number of phi_i, phi_v directions
		self.npv = int(360.0/step_p)
		self.planes = planes

	def load_brdf(self, fname):
		"""Load .bin file from UTIA BRDF binary database
		:param fname (string): filename of .bin BRDF measurement file
		:returns: numpy 1D array with BRDF values
		"""
		self.floats = np.fromfile(fname, dtype=np.float64) #FIXME: maybe translate .bin files to .npy which is platform independent
		return self.floats

	def get_brdf_slice(self, theta_i, theta_v, res=None):
		"""
		:param theta_i (float): incident light elevation angle theta_i (in radians)
		:param theta_v (float): outgoing light elevation angle theta_v (in radians)
		:param res ([int, int]): resolution / number of phi_i, phi_v directions
		:returns: numpy 3D array with patch of BRDF values (single theta_i, theta_v directions)
		"""
		if (res is None):
			res = [self.npi, self.npv]

		refl = []
		for p_i in np.linspace(0, 2*np.pi, res[0]):
			for p_v in np.linspace(0, 2*np.pi, res[1]):
				refl += [brdf.lookup_brdf_val(theta_i, p_i, theta_v, p_v)]
		return np.array(refl).reshape(res[0], res[1], 3)

	def get_brdf_mosaic(self, res_theta=None, res_phi=None):
		"""
		:param res_theta ([int, int]): resolution / number of theta_i, theta_v directions
		:param res_phi ([int, int]): resolution / number of phi_i, phi_v directions
		:returns: numpy 3D array with mosaic of BRDF values (each patch is a single theta_i, theta_v direction)
		"""
		if (res_theta is None):
			res_theta = [self.nti, self.ntv]
		if (res_phi is None):
			res_phi = [self.npi, self.npv]

		rows = []
		for t_i in np.linspace(0, np.pi/2, res_theta[0]):
			row = []
			for t_v in np.linspace(0, np.pi/2, res_theta[1]):
				slc = self.get_brdf_slice(t_i, t_v, res_phi)
				if (len(row) == 0):
					row = slc
				else:
					row = np.hstack((row, slc))
			if (len(rows) == 0):
				rows = row
			else:
				rows = np.vstack((rows, row))
		return rows

	def lookup_brdf_val(self, theta_i, phi_i, theta_v, phi_v):
		"""
		:param theta_i (float): incident light elevation angle theta_i (in radians)
		:param theta_v (float): outgoing light elevation angle theta_v (in radians)
		:param res ([int, int]): resolution / number of phi_i, phi_v directions
		:returns: 3-channel BRDF values for single ingoing and outgoing directions, computed by interpolation
		"""
		pi2 = np.pi/2.
		if (theta_i > pi2) or (theta_v > pi2):
			return [0,0,0]
		d2r = 180.0/np.pi
		theta_i *= d2r
		theta_v *= d2r
		phi_i *= d2r
		phi_v *= d2r
		if (phi_i >= 360.0):
			phi_i = 0.0
		if (phi_v >= 360.0):
			phi_v = 0.0

		iti, itv, ipi, ipv = [0,0], [0,0], [0,0], [0,0]
		iti[0] = int(math.floor(theta_i/self.step_t))
		iti[1] = iti[0]+1
		if (iti[0] > self.nti-2):
			iti[0] = self.nti-2
			iti[1] = self.nti-1
		itv[0] = int(math.floor(theta_v/self.step_t))
		itv[1] = itv[0]+1
		if (itv[0] > self.ntv-2):
			itv[0] = self.ntv-2
			itv[1] = self.ntv-1
		ipi[0] = int(math.floor(phi_i/self.step_p))
		ipi[1] = ipi[0]+1
		ipv[0] = int(math.floor(phi_v/self.step_p))
		ipv[1] = ipv[0]+1

		wti, wtv, wpi, wpv = [0,0], [0,0], [0,0], [0,0]
		wti[1] = theta_i - float(self.step_t*iti[0])
		wti[0] = float(self.step_t*iti[1]) - theta_i
		sum = wti[0]+wti[1]
		wti[0] /= sum
		wti[1] /= sum

		wtv[1] = theta_v - float(self.step_t*itv[0])
		wtv[0] = float(self.step_t*itv[1]) - theta_v
		sum = wtv[0]+wtv[1]
		wtv[0] /= sum
		wtv[1] /= sum

		wpi[1] = phi_i - float(self.step_p*ipi[0])
		wpi[0] = float(self.step_p*ipi[1]) - phi_i
		sum = wpi[0]+wpi[1]
		wpi[0] /= sum
		wpi[1] /= sum
		wpv[1] = phi_v - float(self.step_p*ipv[0])
		wpv[0] = float(self.step_p*ipv[1]) - phi_v
		sum = wpv[0]+wpv[1]
		wpv[0] /= sum
		wpv[1] /= sum

		if (ipi[1] == self.npi):
			ipi[1] = 0
		if (ipv[1] == self.npv):
			ipv[1] = 0

		nc = self.npv*self.ntv
		nr = self.npi*self.nti
		RGB = [0,0,0]
	
		for isp in range(self.planes):
			for i in range(2):
				for j in range(2):
					for k in range(2):
						for l in range(2):
							idx = isp*nr*nc + nc*(self.npi*iti[i] + ipi[k]) + self.npv*itv[j] + ipv[l]
							RGB[isp] += self.floats[idx] * wti[i]*wtv[j]*wpi[k]*wpv[l]
		return RGB

	

if __name__ == "__main__":
	binfile = 'data/m003_carpet01.bin'

	#Example 1 - Basic loading of one material and value lookup
	brdf = BRDF()
	brdf.load_brdf(binfile)

	d2r = np.pi/180.
	theta_i = 50*d2r
	phi_i = 100*d2r
	theta_v = 50*d2r
	phi_v = 0*d2r
	rgb = brdf.lookup_brdf_val(theta_i, phi_i, theta_v, phi_v)
	print('RGB:', rgb)
	
	#Example 2 - restoring azimuthal subspace for fixed elevation angles theta_i=60/theta_v=60 and saving to PNG
	theta_i = 60*d2r
	theta_v = 60*d2r

	brdf = BRDF()
	brdf.load_brdf(binfile)
	slc = brdf.get_brdf_slice(theta_i, theta_v, res=[120, 120])
	slc = np.where(slc < 0, 0, slc)
	slc = np.where(slc > 1.0, 1.0, slc)
	imageio.write_image(255*slc, binfile[:-4]+'.png')
	print('wrote file '+binfile[:-4]+'.png')

	#Example 3 - saving mosaic of BRDF slices to an EXR file (linear gamma RGB)
	#For more detail of the output see: http://btf.utia.cas.cz/?data_des and http://btf.utia.cas.cz/img/brdf/database/BRDF150_25.jpg
	brdf = BRDF()
	brdf.load_brdf(binfile)
	mosaic = brdf.get_brdf_mosaic()
	imageio.writeEXR(mosaic.astype(np.float32), binfile[:-4]+'.exr')
	print('wrote file '+binfile[:-4]+'.exr')
