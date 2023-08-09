from __future__ import print_function

# ---- System Libraries ---- #
from builtins import str
from builtins import range
from builtins import object
import sys
import os
import datetime as dt
import time
import glob
import platform
import struct
import re

import subprocess as subproc

import multiprocessing as mp

from collections import defaultdict

from itertools import groupby

# ---- Math Libraries ---- #
import scipy as sp
import numpy as np
import math
from pyfftw.interfaces import scipy_fftpack as fftw

# ---- Data Classes ---- #
class Experiment(object):
	def __init__(self, pulseseq='slaser', b0=297.2):

		# General Experiment Information
		self.name = ''
		self.author = ''
		self.date = ''
		self.description = ''
		self.type = ''

		# Metabolite Information
		self.insysfiles = []
		self.metabolites = []

		# Experimental Parameters Information
		self.b0 = b0				# default b0 is 7T
		self.obs_iso = '1H'			# proton spectroscopy by default
		self.acq_time = 341E-3		# default acq_time is 341 msec
		self.dwell_time = 0.000166 	# default dwell_time is 0.000166
		self.TE = 0					# echo time

		# Experimental Parameters for Simulation
		if pulseseq == 'slaser_bruker':
			self.TE1 = 0
			self.TE2 = 0

			# amplitudes for pulse calibration
			self.A_90s = []
			self.A_180s = []

			# default is 4.7 ppm (center frequency of excitation pulse)
			self.RF_OFFSET = 4.7

			# RF pulse files
			self.inpulse90file  = ''
			self.inpulse180file = ''

			# RF pulse lengths
			self.PULSE_90_LENGTH  = 0
			self.PULSE_180_LENGTH = 0

			# To store calibrated RF amplitudes
			self.A_180 = 1.0
			self.A_90 = 1.0

		elif pulseseq == 'slaser':
			# amplitudes for pulse calibration
			self.A_90s = []
			self.A_180s = []

			# default is 4.7 ppm (center frequency of excitation pulse)
			self.RF_OFFSET = 4.7

			# RF pulse files
			if self.b0 == 123.3:
				self.inpulse90file  = 'pints/pulses/hsinc_400_8750.pta'
				self.inpulse180file = 'pints/pulses/HS1_R20.pta'
			else:
				self.inpulse90file  = 'pints/pulses/P10.P10NORM.pta'
				self.inpulse180file = 'pints/pulses/HS4_R25.HS4_R25.pta'
			
			# RF pulse lengths
			self.PULSE_90_LENGTH = 0
			self.PULSE_180_LENGTH = 0
			
			# To store calibrated RF amplitudes
			self.A_180 = 1.0
			self.A_90 = 1.0
			self.fudge_factor = 0	# SLR fudge factor
		
		elif pulseseq == 'laser':
			# amplitudes for pulse calibration
			self.A_90s = []
			self.A_180s = []

			# default is 4.7 ppm (center frequency of excitation pulse)
			self.RF_OFFSET = 4.7

			# RF pulse files
			self.inpulse90file = 'pints/pulses/at60.n29.RF'
			self.inpulse180file = 'pints/pulses/HS2_R15_512.RF'
			
			# RF pulse lengths
			self.PULSE_90_LENGTH = 0
			self.PULSE_180_LENGTH = 0
			
			# To store calibrated RF amplitudes
			self.A_180 = 1.0
			self.A_90 = 1.0

		# binning parameters
		self.tolppm = 0.0015
		self.tolpha = 50.0
		self.ppmlo = 0.0
		self.ppmhi = 10.0

		# Experimental Data
		self.data = self.tree()

		# Relationships for .CST and .GES Files
		self.shift_linked_metabs = []
		self.amp_linked_metabs = []
		self.cst_groups = []

		self.LW_limit = 15
		self.LW_linkall = False

		self.splitNAA = False

	def tree(self): return defaultdict(self.tree)

	def getMetabolitesNum(self):
		return sp.size(metabolites)

	def getGyratio(self):
		if self.obs_iso == '1H':
			gyratio = 42576.0
		elif self.obs_iso == '13C':
			gyratio = 67262.0
		elif self.obs_iso == '19F':
			gyratio = 251662.0
		elif self.obs_iso == '23Na':
			gyratio = 70761.0
		elif self.obs_iso == '31P':
			gyratio = 108291.0

		return gyratio

	def getTime(self):
		return sp.arange(0, self.acq_time, self.dwell_time)

	def getFs(self):
		return 1/self.dwell_time

class Pulse(object):
	def __init__(self, inpulsefile, pulse_length, scanner='siemens'):

		if scanner == 'siemens':
			self.waveform = []

			# read pulse (magnitude of pulse is in mT)
			siemens_file = open(inpulsefile, 'r')
			for line in siemens_file:
				if "AMPINT" in line:
					self.AMPINT = float(line.split('\t')[1])
				elif "POWERINT" in line:
					self.POWERINT = float(line.split('\t')[1])
				elif "ABSINT" in line:
					self.ABSINT = float(line.split('\t')[1])
				elif ";" in line:
					mag = float(line.split('\t')[0]) * 1E-3
					phase = float(line.split('\t')[1])
					self.waveform.append(mag * np.exp(1j * phase))
			self.waveform.pop(-1)	# delete last element
			self.waveform = np.array(self.waveform)
			siemens_file.close()

			# define paramaters
			self.pulsestep = pulse_length / len(self.waveform)

		elif scanner == 'varian':
			self.waveform = []
			phases = []
			mags   = []

			# read pulse
			varian_file = open(inpulsefile, 'r')
			for line in varian_file:
				if "#" in line:
					if "TYPE" in line:
						self.TYPE = line.split(' ')[-1]
					elif "MODULATION" in line:
						self.MODULATION = line.split(' ')[-1]
					elif "EXCITEWIDTH" in line:
						self.EXCITEWIDTH = float(line.split(' ')[-1])
					elif "INVERTWIDTH" in line:
						self.INVERTWIDTH = float(line.split(' ')[-1])
					elif "INTEGRAL" in line:
						self.INTEGRAL = float(line.split(' ')[-1])
				else:
					phases.append(np.deg2rad(float([_f for _f in line.replace('\t', ' ').split(' ') if _f][0])))
					mags.append(float([_f for _f in line.replace('\t', ' ').split(' ') if _f][1]))

			# create waveform
			for (i, phase) in enumerate(phases):
				mag   = (mags[i] / np.max(mags))*1E-3
				self.waveform.append(mag * np.exp(1j * phase))
			self.waveform = np.array(self.waveform)
			varian_file.close()

			#define parameters
			self.pulsestep = pulse_length / len(self.waveform)

		elif scanner == 'bruker':
			self.header, data_mag, data_pha = self.__parseBrukerPulseFile__(inpulsefile)
			data_mag = [d for d in data_mag]
			self.waveform = np.array(data_mag)*np.exp(1j*np.deg2rad(data_pha))
			self.header_params = list(self.header.keys())
			self.pulsestep = pulse_length / len(self.waveform)

	def __parseBrukerPulseFile__(self, pulse_file):
		f = open(pulse_file, 'r')
		lines = []
		for line in f:
			lines.append(line.rstrip('\n'))
		f.close()

		header = {}
		data_mag = []
		data_pha = []
		for line in lines:
			if '##' in line and not('##END=' in line):
				key = line.replace('##','').replace('$','').replace(' ','').split('=')[0]
				val = line.replace('##','').replace('$','').replace(' ','').split('=')[1]
				header[key] = val
			elif '##END=' in line:
				
				return header, data_mag, data_pha
			else:
				mag = float(line.replace(' ', '').split(',')[0])
				pha = float(line.replace(' ', '').split(',')[1])
				data_mag.append(mag)
				data_pha.append(pha)
		
		return header, data_mag, data_pha

class RefSignal(object):
	def __init__(self, signal, n, fs, t, b0):
		self.n = n
		self.signal = signal[0:n]
		self.fs = fs
		self.t = t
		self.b0 = b0

class Metabolite(object):
	def __init__(self):

		# General Metabolite Information
		self.name = ''

		# Metabolite Data
		self.peak = []
		self.ppm = []		# i.e. shift
		self.width_L = []
		self.area = []		# i.e. amplitude
		self.phase = []
		self.delay = []
		self.width_G = []
		
		self.crlb = []
		
		# Metabolite Properties
		self.A_m = 1.0	# scaling factor for FID
		self.T2 = 0.0	# T2 relaxation in tissue

		# Properties Specifically for Fitting
		self.T1_GM = 0.0 # T1 relaxation in gray matter
		self.T1_WM = 0.0 # T1 relaxation in white matter
		self.T2_GM = 0.0 # T2 relaxation in gray matter
		self.T2_WM = 0.0 # T2 relaxation in white matter
		self.protons = 0.0 # number of protons

		# A calibration variable
		self.var = 0.0

	def name_short(self):
		return self.name.split('_')[0]

	def num_peaks(self):
		return sp.size(self.area)

	def getFID(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb):

		# b0 should be float [MHz]
		# t should be a numpy array

		# FUDGE FACTOR DUE TO FT1=0.000498 ... something in FITMAN doesn't quite add up?!
		# if pfactor != 0:
		# 	pfactor = pfactor + 4.98

		ppmtoHz = b0    # conversion between ppm and Hz is the main field

		# generate exponentially decaying sinusoids for each peak
		sinusoids = sp.empty([self.num_peaks(), sp.size(t)], dtype=complex)
		for i in range (0, self.num_peaks()):

			# ---- MODEL BASED ON FITMAN ---- #
			c_k = self.area[i] * afactor
			w_k = -(-self.ppm[i] + sfactor) * ppmtoHz
			a_k = lb if sp.size(self.width_L) == 0 else self.width_L[i] + lb
			b_k = 0 if sp.size(self.width_G) == 0 else self.width_G[i]
			phi_k = sp.deg2rad(self.phase[i]) + pfactor	# radians
			t_0 = dfactor if sp.size(self.delay) == 0 else self.delay[i]

			sinusoids[i,] = c_k * \
							sp.exp(1j*(2*sp.pi*w_k*(t+t_0)+phi_k)) * \
							sp.exp(-sp.pi*a_k*np.abs(t+t_0)) * \
							sp.exp(-np.power(sp.pi,2)/(4*np.log(2))*np.power(b_k,2)*np.power(t+t_0,2))


		# ---- OLD MODEL BASED ON VESPA ---- #
		# 	A_k = self.area[i] * afactor
		# 	w_k = 2 * sp.pi * ppmtoHz * (self.ppm[i] - sfactor)
		# 	phi_k = sp.deg2rad(self.phase[i]) + pfactor
		# 	sinusoids[i,] = A_k * sp.exp(-1j*(w_k*t) + 1j*(phi_k))

		# ---- NOTE: BOTH MODELS ARE CONSISTENT ... CAN USE EITHER ONE ---- #

		# generate FID
		FID = sp.sum(sinusoids, axis=0)

		# scale FID by A_m
		FID = self.A_m * FID

		# scale FID by T2 effects (only if a T2 value is specified)
		FID = sp.exp(-(t+TE*10**(-3))/self.T2) * FID if self.T2 > 0 else FID

		# # add line broadening
		# FID = sp.exp(-sp.pi*lb*t) * FID

		return FID

	def getSpec(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb, fs):
		n = sp.size(self.getFID(TE, b0, t, sfactor, afactor, pfactor, dfactor, lb))
		f = sp.arange(+n//2,-n//2,-1)*(fs/n)*(1/b0)

		return (-f, fftw.fftshift(fftw.fft(self.getFID(TE, b0, t, sfactor, afactor, pfactor, dfactor, lb))))

	def energy_FID(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb):
		return np.sum(np.power(np.absolute(self.getFID(TE, b0, t, sfactor, afactor, pfactor, dfactor, lb)), 2))

	def energy_spec(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb):
		return np.sum(np.power(np.absolute(fftw.fftshift(fftw.fft(self.getFID(TE, b0, t, sfactor, afactor, pfactor, dfactor, lb)))), 2))

	def sumAmp(self):
		return np.sum(self.area)

	def calcRatio(self, ref_value):
		return np.sumAmp()/ref_value

class Macromolecule(object):
	def __init__(self, name, shift, line_type, lw, area, phase):

		self.name = name
		self.line_type = line_type
		self.lw = lw
		self.A_m = 1.0
		self.T2 = 0.0

		self.ppm = [shift]
		self.area = [area]
		self.phase = [phase]

	def name_short(self):
		return self.name

	def getFID(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb):
		# b0, lb should be floats (MHz, Hz)
		# t should be a numpy array

		# FUDGE FACTOR DUE TO FT1=0.000498 ... something in FITMAN doesn't quite add up?!
		# if pfactor != 0:
		# 	pfactor = pfactor + 4.98

		ppmtoHz = b0    # conversion between ppm and Hz is the main field

		lw  = self.lw
		A_m = self.A_m
		A_k = self.area[0] * afactor
		w_k = 2 * sp.pi * ppmtoHz * (self.ppm[0] - sfactor)
		phi_k = self.phase[0] + pfactor

		if self.line_type == 'L':
			FID = A_m * A_k * sp.exp(1j*(w_k*t) + 1j*(phi_k)) * sp.exp(-sp.pi*lw*t)
		elif self.line_type == 'G':
			FID = A_m * A_k * sp.exp(1j*(w_k*t) + 1j*(phi_k)) * sp.exp(-np.power(sp.pi,2)/(4*np.log(2)) * np.power(lw,2) * np.power(t,2))

		return FID * sp.exp(-sp.pi*lb*t)

	def getSpec(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb, fs):
		n = sp.size(self.getFID(TE, b0, t, sfactor, afactor, pfactor, lb))
		f = sp.arange(+n//2,-n//2,-1)*(fs/n)*(1/b0)

		return (-f, fftw.fftshift(fftw.fft(self.getFID(TE, b0, t, sfactor, afactor, pfactor, dfactor, lb))))

	def energy_FID(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb):
		return np.sum(np.power(np.absolute(self.getFID(TE, b0, t, sfactor, afactor, pfactor, lb)), 2))

	def energy_spec(self, TE, b0, t, sfactor, afactor, pfactor, dfactor, lb):
		return np.sum(np.power(np.absolute(fftw.fftshift(fftw.fft(self.getFID(TE, b0, t, sfactor, afactor, pfactor, dfactor, lb)))), 2))

class CSTGroup(object):
	def __init__(self, typeCST, name, members, minCST, maxCST):
		self.typeCST = typeCST
		self.name = name

		self.members = members

		self.maxCST = maxCST
		self.minCST = minCST

		self.ref_value = 0
		self.ref_T2 = 0

class OutputFile(object):
	def __init__(self, filename):
		self.filename = filename
		self.metabolites = self.tree()
		self.metabolites_list = []
		self.loadOutputFile()

	def loadOutputFile(self):
		print('=======================================')
		print('Reading data from: ' + self.filename + ' ...')
		
		with open(self.filename, 'r') as f_in:
			lines = (line.rstrip() for line in f_in) 
			lines = list(line for line in lines if line) # Non-blank lines in a list

		output = self.tree()
		crlb   = self.tree()
		constraints_file_flag = False
		for (i, line) in enumerate(lines):
			if "NUMBER_PEAKS" in line:
				self.num_peaks = [_f for _f in line.replace('\r','').replace('\n','').split('\t') if _f][1]
				print(self.num_peaks)
			elif "noise_STDEV_real" in line:
				self.noise_STDEV_real = float([_f for _f in line.split('\t') if _f][1].replace('\r','').replace('\n',''))
				print(self.noise_STDEV_real)
			elif "noise_STDEV_imag" in line:
				self.noise_STDEV_imag = float([_f for _f in line.split('\t') if _f][1].replace('\r','').replace('\n',''))
				print(self.noise_STDEV_imag)
			elif i > 10 and i < int(self.num_peaks) + 11:
				peak = int([_f for _f in line.replace('\r','').replace('\n','').split(' ') if _f][0])
				output[peak] = [_f for _f in line.replace('\r','').replace('\n','').split(' ') if _f]
				print(str(peak) + ' ' + str(output[peak]))
			elif i > 10 and i > int(self.num_peaks) + 11 and i < (2*int(self.num_peaks) + 12):
				peak = int([_f for _f in line.replace('\r','').replace('\n','').replace(';','').split(' ') if _f][0])
				crlb[peak] = [_f for _f in line.replace('\r','').replace('\n','').replace(';','').split(' ') if _f]
				print(str(peak) + ' ' + str(crlb[peak]))
			elif ";PEAK#" in line or ";PEAK #" in line or ";Peak#" in line or ";Peak #" in line:
				constraints_file_flag = True
			elif "Constraints_File_Ends" in line:
				constraints_file_flag = False
			elif constraints_file_flag == True:
				if (";" in line) and not("\t" in line):
					current_metabolite = line.replace(';', '').replace('\r','').replace('\n','')
					print(current_metabolite)
					self.metabolites_list.append(current_metabolite)
					self.metabolites[current_metabolite] = Metabolite()
					self.metabolites[current_metabolite].name = current_metabolite
				else:
					if not(line[0] == ';') and len(line.replace('\r','').replace('\n','')) > 0:
						self.metabolites[current_metabolite].peak.append(int(line.split('\t')[0]))

		print(self.metabolites_list)

		for metabolite in self.metabolites_list:
			for peak in self.metabolites[metabolite].peak:
				self.metabolites[metabolite].ppm.append(float(output[peak][1]))
				self.metabolites[metabolite].width_L.append(float(output[peak][2]))
				self.metabolites[metabolite].area.append(float(output[peak][3]))
				self.metabolites[metabolite].phase.append(sp.rad2deg(float(output[peak][4])))
				self.metabolites[metabolite].delay.append(float(output[peak][5]))
				self.metabolites[metabolite].width_G.append(float(output[peak][6]))
				if float(output[peak][3]) == 0:
					self.metabolites[metabolite].crlb.append(np.nan)
				else:
					self.metabolites[metabolite].crlb.append((float(crlb[peak][3])/float(output[peak][3])) * 100)

				# print metabolite, peak, self.metabolites[metabolite].ppm[-1], self.metabolites[metabolite].width_L[-1], \
				# self.metabolites[metabolite].area[-1], self.metabolites[metabolite].phase[-1], self.metabolites[metabolite].delay[-1], \
				# self.metabolites[metabolite].width_G[-1], self.metabolites[metabolite].crlb[-1]

		print('=======================================')

		# this attribute allows user to reference output data structure directly via peak number 
		# (instead of indexing by metabolite)
		self.output = output
		self.crlbs  = crlb

	def tree(self): return defaultdict(self.tree)

class DatFile(object):
	def __init__(self, filename):
		self.filename = filename
		self.loadDatFile()

	def loadDatFile(self):
		print('=======================================')
		print('Reading dat from ', self.filename, ' ...')		
	
		in_file = open(self.filename,'r')

		self.n = 0
		self.fs = 0
		self.b0 = 0
		dat = []

		for (i, line) in enumerate(in_file):
			if i == 0:
				self.n = int(line.replace(' ', '').replace('\n', ''))//2
			elif i == 2:
				self.fs = 1/float(line.replace(' ', '').replace('\n', ''))
			elif i == 3:
				self.b0 = float(line.replace(' ', '').replace('\n', ''))
			elif i == 7:
				self.ConvS = float(line.replace('\n','').split(' ')[1].split('=')[1])
			elif i == 8:
				self.gain = float(line.replace('\n','').split(' ')[-1].split('=')[1])
				self.TE   = float(line.replace('\n','').split(' ')[0].split('=')[1])
			elif i == 11:
				pass
			elif i > 11:
				dat.append(float(line.replace(' ', '').replace('\n', '')))

		real = sp.array(dat)[0::2]
		imag = sp.array(dat)[1::2]
		print('n', self.n, 'fs', self.fs, 'b0', self.b0)
		print('ConvS', self.ConvS, 'gain', self.gain)
		print('Real', real, sp.size(real))
		print('Imag', imag, sp.size(imag))
		self.signal = real+1j*imag

		self.t = sp.arange(0, self.n, 1) * (1/self.fs)
		print('=======================================')

	def energy_FID(self):
		return np.sum(np.power(np.absolute(self.signal),2))

	def getSpec(self):
		n = sp.size(self.signal)
		f = sp.arange(-n//2,+n//2, 1)*(self.fs/n)*(1/self.b0)
		return (f, fftw.fftshift(fftw.fft(self.signal)))

class RDAFile(object):
	def __init__(self, filename, **kwargs):
		self.filename = filename
		self.kwargs = kwargs
		self.load_rda()

	def load_rda(self):

		filename = self.filename
		kwargs = self.kwargs

		try:
			f = open(filename, "rb")
		except IOError as e:
			print("Could not load {} -- exiting".format(filename))
			return

		flags = dict()
		defaults= {'scale_fid': True}
		for name in list(defaults.keys()):
			flags[name] = kwargs.get(name, defaults[name])

		# Read header of RDA file
		rawfid = f.read()
		hstart = rawfid.find(b">>> Begin of header <<<") + len(b">>> Begin of header <<<\n") + 1
		hend = rawfid.find(b">>> End of header <<<")
		hstring = rawfid[hstart:hend]; hstring = hstring.decode()
		hlist = hstring.split("\r\n")
		hlist2 = [h.split(": ") for h in hlist]
		fid_hdr = dict(hlist2[h] for h in range(0, len(hlist2)-1, 1))
		fid_data = rawfid[hend + len(">>> End of header <<<\n") + 1:]
		f.close()

		# Convert raw spectral data and scale
		element_size = 16  # complex double
		element_count = len(fid_data) // element_size
		format = "dd"

		try:
			fid_data = struct.unpack(format * element_count, fid_data)
		except struct.error:
			print("Unexpected input encountered while reading raw data.")
			return

		fid_data = [complex(fid_data[w], fid_data[w+1]) for w in range(0, len(fid_data), 2)] 

		fid_data_type = "complex128"
		fid_data = np.fromiter(fid_data, fid_data_type)

		# Some acquisition parameters from header
		larmor = float(fid_hdr['MRFrequency'])
		fid_date = dt.datetime.strptime(fid_hdr['StudyDate'], '%Y%m%d').strftime('%d%b%y')
		fid_nt = int(float(fid_hdr["NumberOfAverages"]))	#.rstrip('.0')) <-- if number of averages end with "0", this will be a problem
		fid_dt = float(fid_hdr["DwellTime"]) / 1e6
		# fid_time = np.arange(fid_dt, (len(fid_data) + 1) * fid_dt, fid_dt)

		# Reorder data by size of CSI grid (SVS has grid size of 1 x 1)
		vect_size = int(fid_hdr['VectorSize'])
		fid_time = np.arange(fid_dt, (vect_size + 1) * fid_dt, fid_dt)
		fid_data = fid_data.reshape([len(fid_data)//vect_size, vect_size])
		spec_data = np.zeros(fid_data.shape, dtype=complex)

		final_scale = np.ones(fid_data.shape[0])
		for i_csi in range(fid_data.shape[0]):
			if flags['scale_fid']:
				scale = np.ones(fid_data.shape[1])
				for p in range(48):
					data_pt = fid_data[i_csi, p]
					if np.abs(data_pt) > 10:
						while np.abs(data_pt) > 10:
							scale[p] = scale[p] / 10
							data_pt = data_pt / 10
						final_scale[i_csi] = min(scale)
					elif np.abs(data_pt) < 1:
						while np.abs(data_pt) < 1:
							scale[p] = scale[p] * 10
							data_pt = data_pt * 10
						final_scale[i_csi] = max(scale)
			fid_data[i_csi,:] = fid_data[i_csi,:] * final_scale[i_csi]
			spec_data[i_csi,:], freq = self.fid_to_spec(fid_data[i_csi,:], fid_time)

		csi_size = np.array([fid_hdr['CSIMatrixSize[0]'], fid_hdr['CSIMatrixSize[1]'], fid_hdr['CSIMatrixSize[2]'], vect_size], dtype=int)
		csi_size = csi_size[csi_size > 1]  # Remove singleton dimensions
		# fid_data = fid_data.reshape(csi_size)
		# spec_data = spec_data.reshape(csi_size)
		freq = -freq
		ppm = freq / larmor + 4.7


		vox_size = np.array([fid_hdr["PixelSpacingRow"], fid_hdr["PixelSpacingCol"],
							  fid_hdr["PixelSpacing3D"]],'f8')
		fid_pvec = np.array([fid_hdr["VOIPositionSag"], fid_hdr["VOIPositionCor"],
							  fid_hdr["VOIPositionTra"]],'f8')
		fid_pvec = [-fid_pvec[0], -fid_pvec[1], fid_pvec[2]]
		fid_rvec = np.array([fid_hdr["RowVector[0]"], fid_hdr["RowVector[1]"],
							  fid_hdr["RowVector[2]"]],'f8')
		fid_cvec = np.array([fid_hdr["ColumnVector[0]"], fid_hdr["ColumnVector[1]"],
							  fid_hdr["ColumnVector[2]"]],'f8')
		fid_nvec = np.cross(fid_rvec, fid_cvec)

		R_fid = np.c_[fid_rvec, fid_cvec, fid_nvec]  # concatenate column vectors
		R_fid[:2,2] = -R_fid[:2,2]
		R_fid[2,:2] = -R_fid[2,:2]
		
		angle = np.pi/2
		R_z = np.c_[[np.cos(angle), np.sin(angle), 0],[-np.sin(angle), np.cos(angle), 0], [0, 0, 1]]
		R_fid = np.dot(R_z, R_fid)

		self.patient = fid_hdr["PatientName"]
		self.institute = fid_hdr['InstitutionName']
		self.date = fid_date
		self.larmor = larmor
		self.fid = fid_data
		self.TE = float(fid_hdr['TE'])
		self.TR = float(fid_hdr['TR'])
		self.time = fid_time
		self.spec = spec_data
		self.freq = freq
		self.ppm = ppm
		self.vect_size = vect_size
		self.ConvS = final_scale
		self.vox_center = fid_pvec
		self.n_averages = fid_nt
		self.vox_size = vox_size
		self.vox_affine = R_fid

	def fid_to_spec(self, fid_data, time):
		spec = sp.fftpack.fftshift(sp.fftpack.fft(fid_data)) 
		freq = sp.fftpack.fftshift(sp.fftpack.fftfreq(fid_data.size, time[1] - time[0])) 
		return spec, freq

class Procpar(object):
	def __init__(self, filename):
		self.filename = filename
		self.load_procpar()

	def load_procpar(self):

		f = open(self.filename, 'rb')
		rawproc = f.read().decode()
		rawproc1 = rawproc.split("\n")
		rawproc2 = [x for x in rawproc1 if x != '0 ']
		f.close()

		# Remove exccess characters and data
		y = 0
		rawproc3 = rawproc2
		while (y < len(rawproc2)):
			ele = rawproc2[y]
			if ele[0:2] == '1 ':
				rawproc3[y] = ele[2:]
				y = y+1
			else:
				name = ele.split(' ')
				rawproc3[y] = name[0]
				y = y+1

		self.gain      = float(rawproc3[rawproc3.index("gain")+1])
		self.acqcycles = float(rawproc3[rawproc3.index("acqcycles")+1]) 

		self.lro     = float(rawproc3[rawproc3.index("lro")+1])
		self.lpe     = float(rawproc3[rawproc3.index("lpe")+1])
		self.lpe2    = float(rawproc3[rawproc3.index("lpe2")+1])
		# self.fovunit = float(rawproc3[rawproc3.index("fovunit")+1]) 

		self.pro  = float(rawproc3[rawproc3.index("pro")+1])
		self.ppe  = float(rawproc3[rawproc3.index("ppe")+1])
		self.ppe2 = float(rawproc3[rawproc3.index("ppe2")+1])
		self.pss  = float(rawproc3[rawproc3.index("pss")+1])
		self.pss0 = float(rawproc3[rawproc3.index("pss0")+1])

		self.sw  = float(rawproc3[rawproc3.index("sw")+1])
		self.sw1 = float(rawproc3[rawproc3.index("sw1")+1])
		self.sw2 = float(rawproc3[rawproc3.index("sw2")+1])

		self.psi   = float(rawproc3[rawproc3.index("psi")+1])
		self.phi   = float(rawproc3[rawproc3.index("phi")+1])
		self.theta = float(rawproc3[rawproc3.index("theta")+1])

		self.vpsi   = float(rawproc3[rawproc3.index("vpsi")+1])
		self.vphi   = float(rawproc3[rawproc3.index("vphi")+1])
		self.vtheta = float(rawproc3[rawproc3.index("vtheta")+1])

		self.pos1 = float(rawproc3[rawproc3.index("pos1")+1])
		self.pos2 = float(rawproc3[rawproc3.index("pos2")+1])
		self.pos3 = float(rawproc3[rawproc3.index("pos3")+1])

		self.vox1    = float(rawproc3[rawproc3.index("vox1")+1])
		self.vox2    = float(rawproc3[rawproc3.index("vox2")+1])
		self.vox3    = float(rawproc3[rawproc3.index("vox3")+1])

class FDF2D(object):
	def __init__(self, fdfdir, size):
		self.size    = size
		self.fdfdir  = fdfdir
		print('Loading ... ' + fdfdir)
		self.load_fdf2D()
		self.load_imginfo()
		print('')

	def load_fdf2D(self):
		self.fseimg = np.empty(self.size)
		self.header = self.tree()

		for sl in range(1, self.size[2]):
			f = open(self.fdfdir + '/slice%03dimage001echo001.fdf' % sl, 'rb')
			line = f.readline().decode()

			while not(line == '') and len(line) > 1 and not('checksum' in line):
				line = f.readline().decode()
				var = [_f for _f in line.replace('\n','').replace('=','').replace(';','').replace(',','').replace('*','').replace('[]','').replace('"','').replace('{', '').replace('}', '').split(' ') if _f]
				varval = []
				for i, el in enumerate(var):
					if i == 0:
						vartype = var[i]
					elif i == 1:
						varname = var[i]
					else:
						if vartype == 'float':
							varval.append(float(var[i]))
						elif vartype == 'char':
							varval.append(str(var[i]))
						elif vartype == 'int':
							varval.append(int(var[i]))
				self.header[sl][varname] = np.squeeze(varval)

				if sl == self.size[2]//2:
					print(vartype, varname, varval, f.tell())

			# compute data size
			dataSize = int(np.prod(self.header[sl]['matrix']) * self.header[sl]['bits'] / 8)
			# print ''
			# print dataSize

			# see how much data is left in file
			currentPos = f.tell();
			f.seek(0,2);
			bytesInFile = f.tell()-currentPos;
			# print bytesInFile

			# seek back from end the number of bytes needed
			f.seek(-dataSize, 2)

			# now read the data and reshape
			# if header['storage'] == 'float' and header['bits'] == 32:
			self.fseimg[:,:,sl] = np.fromfile(f, np.float32).reshape(np.int_(self.header[sl]['matrix']), order="F")

			f.close()

		# middle slice of volume
		mid = self.size[2]//2; self.mid = mid

		# signed resolution
		res = np.array([10*(self.header[mid]['span'][0]) / (self.header[mid]['matrix'][0]) , 
						10*(self.header[mid]['span'][1]) / (self.header[mid]['matrix'][1]) ,
						10*(self.header[mid]['roi'][2])])
		self.res = res

	def load_imginfo(self):
		print('')
		
		mid = self.mid
		res = self.res
		fseimg = self.fseimg
		header = self.header

		fse2d_procpar = Procpar(self.fdfdir + '/procpar'); self.procpar = fse2d_procpar

		X_VARIAN = np.zeros(3); X_VARIAN[0] = 1; self.X_VARIAN = X_VARIAN
		Y_VARIAN = np.zeros(3); Y_VARIAN[1] = 1; self.Y_VARIAN = Y_VARIAN
		Z_VARIAN = np.zeros(3); Z_VARIAN[2]	= 1; self.Z_VARIAN = Z_VARIAN

		a = np.deg2rad(90 - fse2d_procpar.psi)
		b = np.deg2rad(fse2d_procpar.theta)
		v = np.deg2rad(fse2d_procpar.phi)
		print('a, b, v:', a,b,v)

		R_img = np.eye(3)
		R_img = np.dot(self.R_z(-a), R_img)
		R_img = np.dot(self.R_y(-b), R_img)
		R_img = np.dot(self.R_z(v),  R_img)
		R_img = np.dot(self.R_y(b),  R_img)
		R_img = np.dot(self.R_z(a),  R_img)
		self.R_img = R_img
		print('R_img:', R_img)

		R_img_scaled = np.array([[res[0],0,0],[0,res[1],0],[0,0,res[2]]])
		R_img_scaled = np.dot(self.R_z(-a), R_img_scaled)
		R_img_scaled = np.dot(self.R_y(-b), R_img_scaled)
		R_img_scaled = np.dot(self.R_z(v),  R_img_scaled)
		R_img_scaled = np.dot(self.R_y(b),  R_img_scaled)
		R_img_scaled = np.dot(self.R_z(a),  R_img_scaled)
		self.R_img_scaled = R_img_scaled

		x_img = np.dot(R_img, X_VARIAN); self.x_img = x_img
		y_img = np.dot(R_img, Y_VARIAN); self.y_img = y_img
		z_img = np.dot(R_img, Z_VARIAN); self.z_img = z_img

		x_img_scaled = np.dot(R_img_scaled, X_VARIAN); self.x_img_scaled = x_img_scaled
		y_img_scaled = np.dot(R_img_scaled, Y_VARIAN); self.y_img_scaled = y_img_scaled
		z_img_scaled = np.dot(R_img_scaled, Z_VARIAN); self.z_img_scaled = z_img_scaled

		p_img_varian    = np.zeros(3)
		p_img_varian[0] = 10*(header[mid]['location'][0])
		p_img_varian[1] = 10*(header[mid]['location'][1])
		p_img_varian[2] = 10*(header[mid]['location'][2])
		self.p_img_varian = p_img_varian

		o_img_varian    = np.zeros(3)
		o_img_varian[0] = 10*(header[1]['origin'][0])
		o_img_varian[1] = 10*(header[1]['origin'][1])
		o_img_varian[2] = 10*(header[1]['location'][2])
		self.o_img_varian = o_img_varian

		fseimg_ijk = []
		fseimg_xyz = []

		for i in range(0,np.size(fseimg, 0)):
			for j in range(0,np.size(fseimg, 1)):
				for k in range(0,np.size(fseimg, 2)):
					xyz_img = i*x_img_scaled + \
								 j*y_img_scaled + \
								 k*z_img_scaled + \
								 o_img_varian
					fseimg_ijk.append([i,j,k])
					fseimg_xyz.append(xyz_img)

		self.fseimg_ijk = fseimg_ijk
		self.fseimg_xyz = fseimg_xyz
		self.fseimg_xyz_kdt = sp.spatial.KDTree(fseimg_xyz)

		R_img_nifti = self.varian_to_nifti(R_img_scaled); self.R_img_nifti = R_img_nifti

		o_img_nifti = o_img_varian - p_img_varian
		o_img_nifti = np.dot(R_img, o_img_nifti)
		o_img_nifti = self.varian_to_nifti(o_img_nifti)
		self.o_img_nifti = o_img_nifti

		p_img_nifti    = self.varian_to_nifti(p_img_varian)
		p_img_nifti[0] = p_img_nifti[0] + o_img_nifti[0]
		p_img_nifti[1] = p_img_nifti[1] + o_img_nifti[1]
		p_img_nifti[2] = p_img_nifti[2] + o_img_nifti[2]
		self.p_img_nifti = p_img_nifti

		affine = np.eye(4)

		affine[0][0] = R_img_nifti[0][0]
		affine[0][1] = R_img_nifti[0][1]
		affine[0][2] = R_img_nifti[0][2]
		affine[0][3] = p_img_nifti[0]

		affine[1][0] = R_img_nifti[1][0]
		affine[1][1] = R_img_nifti[1][1]
		affine[1][2] = R_img_nifti[1][2]
		affine[1][3] = p_img_nifti[1]

		affine[2][0] = R_img_nifti[2][0]
		affine[2][1] = R_img_nifti[2][1]
		affine[2][2] = R_img_nifti[2][2]
		affine[2][3] = p_img_nifti[2]

		self.affine = affine

	def tree(self): return defaultdict(self.tree)

	def R_x(self, angle, handedness='r'):
		if handedness == 'r':
			return np.c_[[1, 0, 0], [0, np.cos(angle), np.sin(angle)],[0, -np.sin(angle), np.cos(angle)]]
		else:
			return np.c_[[1, 0, 0], [0, np.cos(angle), -np.sin(angle)],[0, np.sin(angle), np.cos(angle)]]

	def R_y(self, angle, handedness='r'):
		if handedness == 'r':
			return np.c_[[np.cos(angle), 0, np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]]
		else:
			return np.c_[[np.cos(angle), 0, -np.sin(angle)],[0,1,0],[np.sin(angle),0,np.cos(angle)]]

	def R_z(self, angle, handedness='r'):
		if handedness == 'r':
			return np.c_[[np.cos(angle), np.sin(angle), 0],[-np.sin(angle),np.cos(angle),0],[0,0,1]]
		else:
			return np.c_[[np.cos(angle), -np.sin(angle), 0],[np.sin(angle),np.cos(angle),0],[0,0,1]]

	def R_r(self, angle, r, handedness='r'):
		if handedness == 'r':
			angle = angle
		else:
			angle = -angle

		rx = r[0]; ry = r[1]; rz = r[2]
		cv = np.cos(angle); sv = np.sin(angle)

		R = np.empty([3,3])
		R[0][0] = rx*rx*(1-cv)+cv;		R[0][1] = rx*ry*(1-cv)-rz*sv;	R[0][2] = rx*rz*(1-cv)+ry*sv
		R[1][0] = rx*ry*(1-cv)+rz*sv;	R[1][1] = ry*ry*(1-cv)+cv;		R[1][2] = ry*rz*(1-cv)-rx*sv
		R[2][0] = rx*rz*(1-cv)-ry*sv;	R[2][1] = ry*rz*(1-cv)+rx*sv;	R[2][2] = rz*rz*(1-cv)+cv

		return R

	def varian_to_nifti(self, m_varian):
		m_nifti = m_varian
		m_nifti = np.dot(self.R_x(np.pi/2, 'l'), m_nifti)
		m_nifti = np.dot(self.R_y(np.pi/2, 'l'), m_nifti)
		return m_nifti

class VarianVoxel(object):
	def __init__(self, fiddir, size, X_VARIAN, Y_VARIAN, Z_VARIAN, fseimg_ijk, fseimg_xyz_kdt):
		spec_procpar = Procpar(fiddir + '/procpar'); self.procpar = spec_procpar

		R_vox = np.eye(3)
		a = np.deg2rad(90 - spec_procpar.vpsi)
		b = np.deg2rad(spec_procpar.vtheta)
		v = np.deg2rad(spec_procpar.vphi)

		R_vox = np.eye(3)
		R_vox = np.dot(self.R_z(-a), R_vox)
		R_vox = np.dot(self.R_y(-b), R_vox)
		R_vox = np.dot(self.R_z(v),  R_vox)
		R_vox = np.dot(self.R_y(b),  R_vox)
		R_vox = np.dot(self.R_z(a),  R_vox)
		self.R_vox = R_vox

		vox_size = np.array([spec_procpar.vox1, spec_procpar.vox2, spec_procpar.vox3])
		self.vox_size = vox_size

		x_vox = np.dot(R_vox, X_VARIAN); self.x_vox = x_vox
		y_vox = -np.dot(R_vox, Y_VARIAN); self.y_vox = y_vox
		z_vox = np.dot(R_vox, Z_VARIAN); self.z_vox = z_vox

		p_vox_varian    = np.zeros(3)
		p_vox_varian[0] = np.dot(10*spec_procpar.pos1*x_vox, X_VARIAN)
		p_vox_varian[1] = np.dot(10*spec_procpar.pos2*y_vox, Y_VARIAN)
		p_vox_varian[2] = np.dot(10*spec_procpar.pos3*z_vox, Z_VARIAN)
		self.p_vox_varian = p_vox_varian

		print(10*spec_procpar.pos1, 10*spec_procpar.pos2, 10*spec_procpar.pos3)
		print('p_vox_varian', p_vox_varian)

		voximg = np.zeros(size)
		vox_res = np.array([0.1, 0.1, 0.1]) # arbitrary
		lims = (np.array(vox_size / vox_res / 2)).astype(int)

		for i in range(-lims[0]+1, lims[0]):
			for j in range(-lims[1]+1,lims[1]):
				for k in range(-lims[2]+1,lims[2]):
					xyz_vox = i*x_vox*vox_res[0] + \
								j*y_vox*vox_res[1] +\
								k*z_vox*vox_res[2] + \
								p_vox_varian

					# Find closest image coordinate
					dist, idx = fseimg_xyz_kdt.query(xyz_vox)

					vox_i = fseimg_ijk[idx][0]
					vox_j = fseimg_ijk[idx][1]
					vox_k = fseimg_ijk[idx][2]

					# Set image intensity to 1
					voximg[vox_i][vox_j][vox_k] = 1

		self.voximg = voximg
		print('')

	def R_x(self, angle, handedness='r'):
		if handedness == 'r':
			return np.c_[[1, 0, 0], [0, np.cos(angle), np.sin(angle)],[0, -np.sin(angle), np.cos(angle)]]
		else:
			return np.c_[[1, 0, 0], [0, np.cos(angle), -np.sin(angle)],[0, np.sin(angle), np.cos(angle)]]

	def R_y(self, angle, handedness='r'):
		if handedness == 'r':
			return np.c_[[np.cos(angle), 0, np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]]
		else:
			return np.c_[[np.cos(angle), 0, -np.sin(angle)],[0,1,0],[np.sin(angle),0,np.cos(angle)]]

	def R_z(self, angle, handedness='r'):
		if handedness == 'r':
			return np.c_[[np.cos(angle), np.sin(angle), 0],[-np.sin(angle),np.cos(angle),0],[0,0,1]]
		else:
			return np.c_[[np.cos(angle), -np.sin(angle), 0],[np.sin(angle),np.cos(angle),0],[0,0,1]]

	def R_r(self, angle, r, handedness='r'):
		if handedness == 'r':
			angle = angle
		else:
			angle = -angle

		rx = r[0]; ry = r[1]; rz = r[2]
		cv = np.cos(angle); sv = np.sin(angle)

		R = np.empty([3,3])
		R[0][0] = rx*rx*(1-cv)+cv;		R[0][1] = rx*ry*(1-cv)-rz*sv;	R[0][2] = rx*rz*(1-cv)+ry*sv
		R[1][0] = rx*ry*(1-cv)+rz*sv;	R[1][1] = ry*ry*(1-cv)+cv;		R[1][2] = ry*rz*(1-cv)-rx*sv
		R[2][0] = rx*rz*(1-cv)-ry*sv;	R[2][1] = ry*rz*(1-cv)+rx*sv;	R[2][2] = rz*rz*(1-cv)+cv

		return R

	def varian_to_nifti(self, m_varian):
		m_nifti = m_varian
		m_nifti = np.dot(self.R_x(np.pi/2, 'l'), m_nifti)
		m_nifti = np.dot(self.R_y(np.pi/2, 'l'), m_nifti)
		return m_nifti

class BrukerFID(object):
	def __init__(self, file_dir):
		self.file_dir = file_dir

		# READ RAW DATA
		data_real, data_imag = self.__parseFID__(file_dir)

		# READ ACQUISITION PARAMS
		self.header = self.__parseMethodFile__(file_dir)
		self.header_params = list(self.header.keys())

		print(self.header)

		DigShift = int(68)
		DigDw    = float(0.1248)

		# CHOP OFF ADC DELAY
		self.signal = np.array(data_real) + 1j*np.array(data_imag)
		self.signal = self.signal[DigShift:]
		self.n = np.size(self.signal, 0)

		# SCALE SIGNAL SUCH THAT MAGNITUDE OF FID IS BETWEEN 1 AND 10
		self.ConvS = 1
		# | Find the maximum value of the first 100 points of the time domain signal
		scaled_point = np.max(np.abs(self.signal[0:100]))
		# | Scale the maximum value such that it is between 1 and 10
		if scaled_point < 1:
			while (scaled_point < 1):
				scaled_point = scaled_point * 10.
				self.ConvS = self.ConvS * 10.
		elif scaled_point > 10:
			while scaled_point > 10:
				scaled_point = scaled_point / 10.
				self.ConvS = self.ConvS / 10.
		# | Apply scaling factor to signal
		self.signal = np.real(self.signal) * self.ConvS + 1j*np.imag(self.signal) * self.ConvS

		self.fs = 1/(DigDw/1000)
		self.t = sp.arange(0, self.n, 1) * (1/self.fs)

	def __parseFID__(self, file_dir):
		f = open(file_dir + '/fid', 'r')
		data = np.fromfile(f, np.double)
		data_real = []
		data_imag = []
		for (i, num) in enumerate(data):
			if i % 2 == 0:
				data_real.append(num)
			else:
				data_imag.append(num)
		f.close()
		return data_real, data_imag

	def __parseMethodFile__(self, file_dir):
		f = open(file_dir + '/method', 'r')
		
		lines = []
		for line in f:
			if not('$$ ' in line):
				lines.append(line.rstrip('\n'))
		f.close()

		info = ''
		for line in lines:
			info += line
		param_list = np.array(info.split('##'))

		param_dict = {}
		for param in param_list:
			if param == '':
				continue

			key = param.split('=')[0].replace('$','')
			val = param.split('=')[1]
			
			shape = re.findall(r'\( \d{0,4}\,{0,1}\ {0,1}\d{0,4}\,{0,1}\ {0,1}\d{0,4} \)', val)
			if len(shape) == 0:
				shape = ''
			else:
				shape = shape[0]
			val = val.replace(shape, '')
			
			shape = [int(s) for s in shape.replace('(', '').replace(')', '').replace('$','').replace(' ','').split(',') if len(s) > 0]

			if len(shape) > 0:
				try:
					val = np.reshape(np.array(val.split(' ')), shape)
					val = val.squeeze()
				except ValueError:
					pass
			param_dict[key] = {'value': val, 'shape': shape}

		return param_dict

	def writeDAT(self, out_name, suffix=''):
		if not(suffix is ''):
			out_name = out_name + '_' + suffix + '.dat'
		else:
			out_name = out_name + '.dat'

		if self.Gain is None or len(self.Gain) == 0:
			self.Gain = [1.0]

		now = dt.datetime.now()
		
		data_real = np.real(self.signal)
		data_imag = -np.imag(self.signal)

		DigDw = float(0.1248)
		FrqRef = float(self.header['PVM_FrqRef']['value'][0])
		NAvgs = float(self.header['PVM_NAverages']['value'][0])

		# ... assumes single-voxel spectroscopy here ...
		VoxArrSize = [float(v) for v in self.header['PVM_VoxArrSize']['value']]
		VoxArrPosition = [float(v) for v in self.header['PVM_VoxArrPosition']['value']]

		if np.size(self.header['PVM_EncChanScaling']['value']) > 1:
			EncChanScaling = [float(v) for v in self.header['PVM_EncChanScaling']['value']]
		else:
			EncChanScaling = [float(self.header['PVM_EncChanScaling']['value'])]
		self.Gain = EncChanScaling

		EchoTime = float(self.header['PVM_EchoTime']['value'])
		RepetitionTime = float(self.header['PVM_RepetitionTime']['value'])

		o = open(out_name, 'w')
		o.write(str(np.size(data_real) + np.size(data_imag)) + '\n')
		o.write('1\n')
		o.write(str(DigDw/1000.) + '\n')
		o.write(str(FrqRef) + '\n')
		o.write(str(NAvgs) + '\n')
		o.write(self.file_dir + '/fid\n')
		o.write(now.strftime("%Y %m %d") + '\n')
		o.write('MachS=0 ConvS=' + str(self.ConvS) + ' ')
		o.write('V1=' + str(VoxArrSize[0]) + ' ' + 'V2=' + str(VoxArrSize[1]) + ' ' + 'V3=' + str(VoxArrSize[2]) + '\n')
		o.write('TE=' + str(EchoTime / 1000.) + ' s ')
		o.write('TR=' + str(RepetitionTime / 1000.) + ' s ')
		o.write('P1=' + str(VoxArrPosition[0]) + ' P2=' + str(VoxArrPosition[1]) + ' P3=' + str(VoxArrPosition[2]) + ' Gain=' + str(EncChanScaling[0]) + '\n')
		o.write('SIMULTANEOUS\n0.0\n')
		o.write('EMPTY\n')

		for i, p in enumerate(self.signal):
			o.write(str(data_real[i]) + '\n')
			o.write(str(data_imag[i]) + '\n')
		o.close()

	def getSpec(self):
		n = sp.size(self.signal)
		f = sp.arange(-n//2,+n//2, 1)*(self.fs/n)*(1/float(self.header['PVM_FrqRef']['value'][0]))
		return (-f, fftw.fftshift(fftw.fft(self.signal)))

	def print_params(self, console):
		# PRINT SOME USEFUL PARAMS TO CONSOLE

		EchoTime = self.header['PVM_EchoTime']['value']
		RepetitionTime = self.header['PVM_RepetitionTime']['value']
		NAverages = self.header['PVM_NAverages']['value']
		FrqRef = self.header['PVM_FrqRef']['value'][0]
		DigDw = 0.1248
		DigShift = 68

		# ... assumes single-voxel spectroscopy here ...
		VoxArrSize = [float(v) for v in self.header['PVM_VoxArrSize']['value']]
		VoxArrPosition = [float(v) for v in self.header['PVM_VoxArrPosition']['value']]
		VoxArrPositionRPS = [float(v) for v in self.header['PVM_VoxArrPositionRPS']['value']]

		if np.size(self.header['PVM_EncChanScaling']['value']) > 1:
			EncChanScaling = [float(v) for v in self.header['PVM_EncChanScaling']['value']]
		else:
			EncChanScaling = [float(self.header['PVM_EncChanScaling']['value'])]
		self.Gain = EncChanScaling

		console.append(' | file_dir ' + str(self.file_dir))
		console.append(' | EchoTime ' + str(EchoTime))
		console.append(' | RepetitionTime ' + str(RepetitionTime))
		console.append(' | NAverages ' + str(NAverages))
		console.append(' | FrqRef ' + str(FrqRef))
		console.append(' | DigDw ' + str(DigDw))
		console.append(' | DigShift ' + str(DigShift))
		console.append(' | VoxArrSize ' + str(VoxArrSize))
		console.append(' | VoxArrPosition ' + str(VoxArrPosition))
		console.append(' | VoxArrPositionRPS ' + str(VoxArrPositionRPS))
		console.append(' | EncChanScaling ' + str(EncChanScaling))
		console.append(' | n ' + str(self.n))
		console.append(' | ConvS ' + str(self.ConvS))
		console.append(' | fs ' + str(self.fs))
		console.append(' | t' + str(self.t))