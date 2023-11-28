from __future__ import print_function

# ---- System Libraries ---- #
from builtins import zip
from builtins import str
from builtins import range
import sys
import os
import datetime
import time
import platform

# Check for PyQt6 (native mac M1 compatibility), otherwise continue using PyQt5
import importlib.util
PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore, QtGui, QtWidgets, uic
else:
	from PyQt5 import QtCore, QtGui, QtWidgets, uic

from collections import defaultdict

# ---- Math Libraries ---- #
import scipy as sp
import scipy.signal as spsg

import numpy as np
import math
from pyfftw.interfaces import scipy_fftpack as fftw

# ---- Plotting Libraries ---- #
import matplotlib as mpl;
mpl.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (
	FigureCanvasQTAgg as FigureCanvas,
	NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt; plt.style.use('ggplot')

# ---- Data Classes ---- #
from magiqdataclasses import *

# ---- Simulation Libraries ---- #
import pygamma as pg
from simclasses import *

qtCreatorFile = "pints/ui/PINTS.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

# ---- Main Application Class ---- #
class MyApp(QtWidgets.QWidget, Ui_MainWindow):
	def __init__(self):
		QtWidgets.QWidget.__init__(self)
		Ui_MainWindow.__init__(self)
		self.setupUi(self)

		# Bind buttons to methods in each tab
		self.setBindings('Simulate')
		self.setBindings('Load Experiment')
		self.setBindings('Edit Metabolites')
		self.setBindings('Plot')
		self.setBindings('Export')
		self.setBindings('EditCST')
		self.setBindings('Guess')

		# Setup for plotting in the "Plot Metabolites" and "Generate Guess File" Tabs
		self.canvas = [None]*3
		self.toolbar = [None]*3
		self.setPlot('Plot')
		self.setPlot('Guess')

	def setBindings(self, tab):
		if tab == 'Simulate':

			self.confirmSimParamsButton.clicked.connect(self.confirmSimParams)
			
			self.calibrationMetaboliteComboBox.setCurrentIndex(35)
			self.calibrationMetaboliteComboBox_bruker.setCurrentIndex(35)
			self.calibrationMetaboliteComboBox_laser.setCurrentIndex(35)

			self.loadBrukerDataBrowseButton.clicked.connect(self.loadBrukerData)
			self.workingDirectory_bruker = os.path.expanduser('~')
			self.loadBrukerDataLoadButton.clicked.connect(self.loadBrukerParams)
			self.loadBrukerDataLoadButton.setEnabled(False)

			self.runSimulationButton.setEnabled(False)
			self.runSimulationButton.clicked.connect(self.runSimulation)

		elif tab == 'Load Experiment':

			self.filenameBrowseButton.clicked.connect(self.chooseExperimentFile)

			self.filenameConfirmButton.clicked.connect(self.loadExperimentFile)
			self.filenameConfirmButton.setEnabled(False)
			
			self.specExpButton.clicked.connect(self.specExp)
			self.specExpButton.setEnabled(False)

		elif tab == 'Edit Metabolites':

			self.editMetabOKButton.clicked.connect(self.saveMetabInfo)
			self.editMetabOKButton.setEnabled(False)
			
			self.loadMetabConfirmButton.clicked.connect(self.loadMetabInfoFile)
			self.loadMetabConfirmButton.setEnabled(False)

			self.loadMetabBrowseButton.clicked.connect(self.chooseMetabInfoFile)
			self.loadMetabBrowseButton.setEnabled(False)

		elif tab == 'Plot':

			self.plotButton.clicked.connect(self.plot)
			self.plotButton.setEnabled(False)

		elif tab == 'Export':

			self.plotButton.clicked.connect(self.plot)

			self.outputFilenameConfirmButton_rawTime.clicked.connect(self.setOutputFile_rawTime)
			self.outputFilenameConfirmButton_rawTime.setEnabled(False)

			self.outputFilenameConfirmButton_rawFreq.clicked.connect(self.setOutputFile_rawFreq)
			self.outputFilenameConfirmButton_rawFreq.setEnabled(False)

			self.outputFilenameConfirmButton_dat.clicked.connect(self.setOutputFile_dat)
			self.outputFilenameConfirmButton_dat.setEnabled(False)

			self.exportButton_rawTime.clicked.connect(self.export_rawTime)
			self.exportButton_rawTime.setEnabled(False)
			
			self.exportButton_rawFreq.clicked.connect(self.export_rawFreq)
			self.exportButton_rawFreq.setEnabled(False)
			
			self.exportButton_dat.clicked.connect(self.export_dat)
			self.exportButton_dat.setEnabled(False)

		elif tab == 'EditCST':

			self.cstOKButton.clicked.connect(self.saveCSTInfo)
			self.cstOKButton.setEnabled(False)
			
			self.cstLoadButton.clicked.connect(self.loadCSTInfo)
			self.cstLoadButton.setEnabled(False)
			
			self.cstFilenameConfirmButton.clicked.connect(self.setCSTFile)
			self.cstFilenameConfirmButton.setEnabled(False)
			
			self.genCSTButton.clicked.connect(self.genCST)
			self.genCSTButton.setEnabled(False)

		elif tab == 'Guess':

			self.filenameBrowseButton_dat.clicked.connect(self.chooseDATFile)
			self.filenameBrowseButton_dat.setEnabled(False)

			self.filenameConfirmButton_dat.clicked.connect(self.loadDATFile)
			self.filenameConfirmButton_dat.setEnabled(False)

			self.plotGESButton.clicked.connect(self.plotGES)
			self.plotGESButton.setEnabled(False)

			self.gesFilenameConfirmButton.clicked.connect(self.setGESFile)
			self.gesFilenameConfirmButton.setEnabled(False)

			self.genGESButton.clicked.connect(self.genGES)
			self.genGESButton.setEnabled(False)

	# ---- Methods for Simulate Tab ---- #
	def loadBrukerData(self):
		prev = str(self.loadBrukerDataLineEdit.text())
		self.loadBrukerDataLineEdit.setText(str(QtWidgets.QFileDialog.getExistingDirectory(self, 'Open Bruker Data Directory', self.workingDirectory_bruker)))
		if str(self.loadBrukerDataLineEdit.text()) == '':
			if str(prev) == '':
				self.loadBrukerDataLineEdit.setText(str(prev))
				self.loadBrukerDataLoadButton.setEnabled(False)
				
				self.T3Button.setEnabled(True)
				self.T7Button.setEnabled(True)
				self.T9Button.setEnabled(True)

				self.dwellTimeInput_sim.setEnabled(True)
				self.acqLengthInput_sim.setEnabled(True)
				self.echoTimeInput_sim.setEnabled(True)

				self.sLASERradioButton.setEnabled(True)
				self.sLASERradioButton_bruker.setEnabled(True)
				self.LASERradioButton.setEnabled(False)
		else:
			self.loadBrukerDataLoadButton.setEnabled(True)

			self.T3Button.setEnabled(False)
			self.T7Button.setEnabled(False)
			self.T9Button.setEnabled(False)

			self.dwellTimeInput_sim.setEnabled(False)
			self.acqLengthInput_sim.setEnabled(False)
			self.echoTimeInput_sim.setEnabled(False)

			self.sLASERradioButton.setEnabled(False)
			self.sLASERradioButton_bruker.setEnabled(False)
			self.LASERradioButton.setEnabled(False)

			self.workingDirectory_bruker = os.path.abspath(os.path.join(os.path.expanduser(str(prev)), os.pardir))

	def loadBrukerParams(self):
		file_dir = str(self.loadBrukerDataLineEdit.text())
		data = BrukerFID(file_dir)

		self.sLASERradioButton.setChecked(False)
		self.sLASERradioButton_bruker.setChecked(True)
		self.LASERradioButton.setChecked(False)

		# Experiment Info
		b0 = float(data.header['PVM_FrqRef']['value'][0])
		# This if-statement is just a sanity check.
		if b0 > 350:
			self.T3Button.setChecked(False)
			self.T7Button.setChecked(False)
			self.T9Button.setChecked(True)
		elif b0 > 200:
			self.T3Button.setChecked(False)
			self.T7Button.setChecked(True)
			self.T9Button.setChecked(False)
		else:
			self.T3Button.setChecked(True)
			self.T7Button.setChecked(False)
			self.T9Button.setChecked(False)

		dt   = float(data.header['PVM_SpecDwellTime']['value'].item())*2/1e6  # in seconds
		self.dwellTimeInput_sim.setText(str(dt))

		acqt = dt * int(data.header['PVM_DigNp']['value'])
		self.acqLengthInput_sim.setText(str(acqt))

		digshift = int(data.digShift)
		self.digShiftInput_bruker.setText(str(digshift))

		te   = float(data.header['PVM_EchoTime']['value'])
		self.echoTimeInput_sim.setText(str(te))

		te1  = int(float(data.header['TE1']['value']))
		te2  = int(float(data.header['TE2']['value']))
		self.TE1_bruker.setText(str(te1))
		self.TE2_bruker.setText(str(te2))

		# sLASER Pulse Info
		# | Pulse Structure (pulse length, pulse bandwith, flip angle, excitation, ~, ~, ~, ~, ~, amplitude, shape)
		
		# Excitation
		ep  = data.header['VoxPul1']['value'].replace('(','').replace(' ','').replace(')','').split(',')
		epn = data.header['VoxPul1Enum']['value'].replace('>','').replace('<','')
		self.excPulse_bruker = {'name': epn, 'params': ep}

		plen = float(ep[0]) * 1000
		self.excPulseLength_bruker.setText(str(plen))

		pamp = float(ep[-2]) / 10
		self.excAmpMin_bruker.setText(str(0))
		self.excAmpMax_bruker.setText(str(pamp))

		# Refocussing
		rp  = data.header['VoxPul2']['value'].replace('(','').replace(' ','').replace(')','').split(',')
		rpn = data.header['VoxPul2Enum']['value'].replace('>','').replace('<','')
		self.rfcPulse_bruker = {'name': rpn, 'params': rp}

		plen = float(rp[0]) * 1000
		self.afpPulseLengthInput_bruker.setText(str(plen))

		pamp = float(rp[-2]) / 10
		self.afpAmpMin_bruker.setText(str(0))
		self.afpAmpMax_bruker.setText(str(pamp))

		# Editing Pulse Info
		# --- TBD ---

	def fit_sin(self, tt, yy):
		'''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
		tt = np.array(tt)
		yy = np.array(yy)
		ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
		Fyy = abs(np.fft.fft(yy))
		guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
		guess_amp = np.std(yy) * 2.**0.5
		guess_offset = np.mean(yy)
		guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

		popt, pcov = sp.optimize.curve_fit(self.sinfunc, tt, yy, p0=guess)
		A, w, p, c = popt
		f = w/(2.*np.pi)
		fitfunc = lambda t: A * np.sin(w*t + p) + c
		return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}

	def sinfunc(self, t, A, w, p, c):
		return A * np.sin(w*t + p) + c

	def logs_func(self, x, A, B, K, s):
		return A + (K - A)/(1 + np.exp(-B*(x-s)))

	def confirmSimParams(self):
		if self.sLASERradioButton_bruker.isChecked():

			if self.T7Button.isChecked():
				self.sim_experiment = Experiment('slaser_bruker', 297.2)
			elif self.T3Button.isChecked():
				self.sim_experiment = Experiment('slaser_bruker', 123.3)
			elif self.T9Button.isChecked():
				self.sim_experiment = Experiment('slaser_bruker', 400.2)
			self.sim_experiment.name = 'semi-LASER (Bruker)'

			if 'Calculated' in self.excPulse_bruker['name']:
				self.sim_experiment.inpulse90file = 'pints/pulses/rect.exc'
			else:
				self.sim_experiment.inpulse90file = 'pints/pulses/'+self.excPulse_bruker['name']+'.exc'
			
			self.sim_experiment.inpulse180file = 'pints/pulses/'+self.rfcPulse_bruker['name']+'.inv'

		elif self.sLASERradioButton.isChecked():

			if self.T7Button.isChecked():
				self.sim_experiment = Experiment('slaser', 297.2)
			elif self.T3Button.isChecked():
				self.sim_experiment = Experiment('slaser', 123.3)
			elif self.T9Button.isChecked():
				self.sim_experiment = Experiment('slaser', 400.2)
			self.sim_experiment.name = 'semi-LASER'

		elif self.LASERradioButton.isChecked():

			if self.T7Button.isChecked():
				self.sim_experiment = Experiment('laser', 297.2)
			elif self.T3Button.isChecked():
				self.sim_experiment = Experiment('laser', 123.3)
			elif self.T9Button.isChecked():
				self.sim_experiment = Experiment('laser', 400.2)
			self.sim_experiment.name = 'LASER'
		
		self.sim_experiment.author = platform.node()
		self.sim_experiment.date = datetime.datetime.now().strftime("%Y_%m_%d")
		self.sim_experiment.description = 'Density-matrix simulations of metabolites using PyGAMMA and metabolite parameters from V. Govindaraju et. al, NMR Biomed. 2000;13:129-153'
		self.sim_experiment.type = 'PINTS'

		self.insysfiles = [     
					'alanine.sys',
					'aspartate.sys',
					'choline_1-CH2_2-CH2.sys',
					'choline_N(CH3)3_a.sys',
					'choline_N(CH3)3_b.sys',
					'creatine_N(CH3).sys',
					'creatine_X.sys',
					'd-glucose-alpha.sys',
					'd-glucose-beta.sys',
					'eth.sys',
					'gaba.sys',
					'glutamate.sys',
					'glutamine.sys',
					'glutathione_cysteine.sys',
					'glutathione_glutamate.sys',
					'glutathione_glycine.sys',
					'glycine.sys',
					'gpc_7-CH2_8-CH2.sys',
					'gpc_glycerol.sys',
					'gpc_N(CH3)3_a.sys',
					'gpc_N(CH3)3_b.sys',
					'lactate.sys',
					'myoinositol.sys',
					'naa_acetyl.sys',
					'naa_aspartate.sys',
					'naag_acetyl.sys',
					'naag_aspartyl.sys',
					'naag_glutamate.sys',
					'pcho_N(CH3)3_a.sys',
					'pcho_N(CH3)3_b.sys',
					'pcho_X.sys',
					'pcr_N(CH3).sys',
					'pcr_X.sys',
					'peth.sys',
					'scyllo-inositol.sys',
					'taurine.sys']

		if self.sim_experiment.b0 == 123.3:
			self.macromolecules = [ 
						'MM09', 
						'MM12', 
						'MM14', 
						'MM16', 
						'MM20', 
						'MM21', 
						'MM23', 
						'MM26', 
						'MM30', 
						'MM31', 
						'MM37', 
						'MM38', 
						'MM40' ]

			# From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5215417/
			self.macromolecules_data = self.tree()
			self.macromolecules_data['MM09'] = Macromolecule('MM09',  0.90, 'L', 21.20, 0.72, 0)
			self.macromolecules_data['MM12'] = Macromolecule('MM12',  1.21, 'L', 19.16, 0.28, 0)
			self.macromolecules_data['MM14'] = Macromolecule('MM14',  1.38, 'L', 15.90, 0.38, 0)
			self.macromolecules_data['MM16'] = Macromolecule('MM16',  1.63, 'L', 07.50, 0.05, 0)
			self.macromolecules_data['MM20'] = Macromolecule('MM20',  2.01, 'L', 29.03, 0.45, 0)
			self.macromolecules_data['MM21'] = Macromolecule('MM21',  2.09, 'L', 20.53, 0.36, 0)
			self.macromolecules_data['MM23'] = Macromolecule('MM23',  2.25, 'L', 17.89, 0.36, 0)
			self.macromolecules_data['MM26'] = Macromolecule('MM26',  2.61, 'L', 05.30, 0.04, 0)
			self.macromolecules_data['MM30'] = Macromolecule('MM30',  2.96, 'L', 14.02, 0.20, 0)
			self.macromolecules_data['MM31'] = Macromolecule('MM31',  3.11, 'L', 17.89, 0.11, 0)
			self.macromolecules_data['MM37'] = Macromolecule('MM37',  3.67, 'L', 33.52, 0.64, 0)
			self.macromolecules_data['MM38'] = Macromolecule('MM38',  3.80, 'L', 11.85, 0.07, 0)
			self.macromolecules_data['MM40'] = Macromolecule('MM40',  3.96, 'L', 37.48, 1.00, 0)
		else:
			self.macromolecules = [ 'lm1', 
						'lm2', 
						'lm3', 
						'lm4', 
						'lm5', 
						'lm6', 
						'lm7', 
						'lm8', 
						'lm9', 
						'lm10', 
						'lm11', 
						'lm12', 
						'lm13', 
						'lm14' ]

			self.macromolecules_data = self.tree()
			self.macromolecules_data['lm1']  = Macromolecule('lm1',  0.9457810, 'L', 50.72009, 3.5501260, -1.13820100)
			self.macromolecules_data['lm2']  = Macromolecule('lm2',  1.4757810, 'L', 54.12614, 2.8557010, 0.054570260)
			self.macromolecules_data['lm3']  = Macromolecule('lm3',  1.7057810, 'L', 68.32043, 3.9072800, -0.06015777)
			self.macromolecules_data['lm4']  = Macromolecule('lm4',  2.1157810, 'L', 96.43580, 4.7191390, -0.09620739)
			self.macromolecules_data['lm5']  = Macromolecule('lm5',  2.3157810, 'L', 69.70896, 2.7845620, -0.57759460)
			self.macromolecules_data['lm6']  = Macromolecule('lm6',  3.0157810, 'L', 42.72932, 1.0176800, -0.65827850)
			self.macromolecules_data['lm7']  = Macromolecule('lm7',  3.9854041, 'L', 29.64625, 1.2094070, 0.288971500)
			self.macromolecules_data['lm8']  = Macromolecule('lm8',  7.1257810, 'L', 274.4760, 4.0466410, 0.302700400)
			self.macromolecules_data['lm9']  = Macromolecule('lm9',  7.8557810, 'L', 37.37627, 1.3135740, 0.020221810)
			self.macromolecules_data['lm10'] = Macromolecule('lm10', 2.6557810, 'L', 121.4557, 2.3499680, -0.58261020)
			self.macromolecules_data['lm11'] = Macromolecule('lm11', 1.2357810, 'L', 30.68139, 0.7444791, -1.13820100)
			self.macromolecules_data['lm12'] = Macromolecule('lm12', -0.628240, 'L', 1789.202, 21.518780, -3.65296700)
			self.macromolecules_data['lm13'] = Macromolecule('lm13', 3.2957810, 'L', 109.4958, 1.6655390, 0.163575600)
			self.macromolecules_data['lm14'] = Macromolecule('lm14', 3.7057810, 'L', 169.6304, 4.8735330, -1.35431900)

		self.sim_experiment.obs_iso = '1H'

		try:
			self.sim_experiment.acq_time = float(self.acqLengthInput_sim.text())
			self.sim_experiment.dwell_time = float(self.dwellTimeInput_sim.text())
			self.sim_experiment.TE = float(self.echoTimeInput_sim.text())

			if self.sLASERradioButton_bruker.isChecked():

				self.sim_experiment.TE1 = float(self.TE1_bruker.text())
				self.sim_experiment.TE2 = float(self.TE2_bruker.text())

				self.sim_experiment.PULSE_90_LENGTH = float(self.excPulseLength_bruker.text())*10**(-6)
				self.sim_experiment.PULSE_180_LENGTH = float(self.afpPulseLengthInput_bruker.text())*10**(-6)
				self.sim_experiment.RF_OFFSET = float(self.rfOffsetInput_bruker.text())
				
				amin = int(np.floor(float(self.excAmpMin_bruker.text()))); amax = int(np.ceil(float(self.excAmpMax_bruker.text())))
				self.sim_experiment.A_90s  = np.linspace(amin, amax, (amax-amin)*10 + 1)
				amin = int(np.floor(float(self.afpAmpMin_bruker.text()))); amax = int(np.ceil(float(self.afpAmpMax_bruker.text())))
				self.sim_experiment.A_180s = np.linspace(amin, amax, (amax-amin)*10 + 1)

				self.sim_experiment.DigShift = int(self.digShiftInput_bruker.text())*self.sim_experiment.dwell_time

				# Editing Pulse Info
				# --- TBD ---

			elif self.sLASERradioButton.isChecked():
				
				self.sim_experiment.PULSE_90_LENGTH = float(self.slrPulseLengthInput.text())*10**(-6)
				self.sim_experiment.fudge_factor = int(self.fudgeFactorInput.text())
				self.sim_experiment.PULSE_180_LENGTH = float(self.afpPulseLengthInput.text())*10**(-6)
				self.sim_experiment.RF_OFFSET = float(self.rfOffsetInput.text())
				
				self.sim_experiment.A_90s  = np.linspace(int(self.slrAmpMinInput.text()), int(self.slrAmpMaxInput.text()), (int(self.slrAmpMaxInput.text())-int(self.slrAmpMinInput.text()))*10 + 1)
				self.sim_experiment.A_180s = np.linspace(int(self.afpAmpMinInput.text()), int(self.afpAmpMaxInput.text()), (int(self.afpAmpMaxInput.text())-int(self.afpAmpMinInput.text()))*10 + 1)

			elif self.LASERradioButton.isChecked():
				
				self.sim_experiment.PULSE_90_LENGTH = float(self.ahpPulseLengthInput_laser.text())*10**(-6)
				self.sim_experiment.PULSE_180_LENGTH = float(self.afpPulseLengthInput_laser.text())*10**(-6)
				self.sim_experiment.RF_OFFSET = float(self.rfOffsetInput_laser.text())

				self.sim_experiment.A_90s = np.linspace(int(self.ahpAmpMinInput_laser.text()), int(self.ahpAmpMaxInput_laser.text()), (int(self.ahpAmpMaxInput_laser.text())-int(self.ahpAmpMinInput_laser.text()))*10 + 1)
				self.sim_experiment.A_180s = np.linspace(int(self.afpAmpMinInput_laser.text()), int(self.afpAmpMaxInput_laser.text()), (int(self.afpAmpMaxInput_laser.text())-int(self.afpAmpMinInput_laser.text()))*10 + 1)

			self.sim_experiment.tolppm = float(self.tolppmInput.text())
			self.sim_experiment.tolpha = float(self.tolphaInput.text())
			self.sim_experiment.ppmlo = float(self.ppmloInput.text()) - self.sim_experiment.RF_OFFSET
			self.sim_experiment.ppmhi = float(self.ppmhiInput.text()) - self.sim_experiment.RF_OFFSET

			self.simConsole.clear()
			self.simConsole.append('Confirming simulation experiment parameters ...')
			self.simConsole.append(' | b0:               ' + str(self.sim_experiment.b0) + ' MHz')
			self.simConsole.append(' | observed species: ' + str(self.sim_experiment.obs_iso))
			self.simConsole.append(' | gyratio:          ' + str(self.sim_experiment.getGyratio()) + ' MHz/mT')
			self.simConsole.append(' | dwell time:       ' + str(self.sim_experiment.dwell_time) + ' sec')
			self.simConsole.append(' | acqusition time:  ' + str(self.sim_experiment.acq_time) + ' sec')
			self.simConsole.append(' | TE:               ' + str(self.sim_experiment.TE) + ' msec')
			self.simConsole.append(' | A_90s:           ' + str(np.array(self.sim_experiment.A_90s)))
			self.simConsole.append(' | A_180s:           ' + str(np.array(self.sim_experiment.A_180s)))
			self.simConsole.append(' | rf_offset:        ' + str(self.sim_experiment.RF_OFFSET) + ' ppm')

			self.confirmSimParamsButton.setEnabled(False)
			self.runSimulationButton.setEnabled(True)
		except Exception as e:
			self.simConsole.append('\nERROR: ' + str(e) + '\n>> Please check your simulation parameters.\n')

	def runSimulation(self):

		self.runSimulationButton.setEnabled(False)

		self.simProgressBar.setMaximum(int(np.size(self.insysfiles)) + 2)
		self.simProgressBar.setValue(0)

		if self.sLASERradioButton_bruker.isChecked():

			# Create a new simulations results file
			self.save_dir_sim = 'pints/experiments/sLASER_sim_bruker_' \
								+ datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S") \
								+ '_TE' + str(self.sim_experiment.TE)
			
			if self.T7Button.isChecked():
				self.save_dir_sim += '_7T/'
			elif self.T3Button.isChecked():
				self.save_dir_sim += '_3T/'
			elif self.T9Button.isChecked():
				self.save_dir_sim += '_9T/'

			if not os.path.exists(self.save_dir_sim):
				os.makedirs(self.save_dir_sim)
			self.sim_results = open(self.save_dir_sim + 'sLASER_sim_results.txt', 'w')

			# Write parameters in file
			self.sim_results.write(';PINTS for FITMAN Simulation Output\n')
			self.sim_results.write(';Experiment Information\n')
			self.sim_results.write(';---------------------------------------------------------------------------\n')
			self.sim_results.write(';Name: ' + self.sim_experiment.name + '\n')
			self.sim_results.write(';Created: ' + self.sim_experiment.date + '\n')
			self.sim_results.write(';Comment: ' + self.sim_experiment.description + '\n')
			self.sim_results.write(';PI: ' + self.sim_experiment.author + '\n')
			self.sim_results.write(';b0: ' + str(self.sim_experiment.b0) + '\n')
			self.sim_results.write(';' + str(int(np.size(self.insysfiles)+int(np.size(self.macromolecules)))) + ' Metabolites: ' + str(self.insysfiles).replace('[','').replace(']','').replace('\'','').replace('.sys',''))
			if self.macroIncludeButton.isChecked(): self.sim_results.write(', ' + str(self.macromolecules).replace('[','').replace(']','').replace('\'',''))
			self.sim_results.write('\n')

			self.sim_results.write(';Simulation Results\n')
			self.sim_results.write(';---------------------------------------------------------------------------\n')

			# Set up for calibration experiments
			insysfile = str(self.calibrationMetaboliteComboBox.currentText()) + '.sys'
			if self.sim_experiment.b0 == 123.3:
				insysfile = 'pints/metabolites/3T_' + insysfile
			elif self.sim_experiment.b0 == 297.2:
				insysfile = 'pints/metabolites/7T_' + insysfile
			elif self.sim_experiment.b0 == 400.2:
				insysfile = 'pints/metabolites/9.4T_' + insysfile
			print('')
			print(insysfile)

			spin_system = pg.spin_system()
			spin_system.read(insysfile)
			for i in range(spin_system.spins()):
				spin_system.PPM(i, spin_system.PPM(i) - self.sim_experiment.RF_OFFSET)

			# Run the 180 calibration
			self.simConsole.append('\n1. 180-degree calibration sLASER (w/ ideal 90, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			A180_calibration_data = []

			for A_180 in self.sim_experiment.A_180s:

				TE = self.sim_experiment.TE * 1E-3
				TE1 = self.sim_experiment.TE1 * 1E-3
				TE2 = self.sim_experiment.TE2 * 1E-3

				# build 90 degree pulse
				pulse_dur_90 = 0

				# build 180 degree pulse
				A_180, pulse180, pulse_dur_180, Ureal180 = self.slaser_build180(self.sim_experiment.inpulse180file, A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, False, 'bruker')

				H = pg.Hcs(spin_system) + pg.HJ(spin_system)
				D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
				ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
				ACQ = ac

				delay1 = TE1/2.0 - pulse_dur_90/2.0 - pulse_dur_180/2.0
				delay2 = TE1/2.0 + TE2/2.0 - pulse_dur_180
				delay3 = TE2 - pulse_dur_180
				delay4 = delay2
				delay5 = TE1/2.0 - pulse_dur_180 + self.sim_experiment.DigShift

				# TE_fill = TE - 2.*TE1 - 2.*TE2
				# delay1 = TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay2 = TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay3 = TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay4 = TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay5 = TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + self.sim_experiment.DigShift

				Udelay1 = pg.prop(H, delay1)
				Udelay2 = pg.prop(H, delay2)
				Udelay3 = pg.prop(H, delay3)
				Udelay4 = pg.prop(H, delay4)
				Udelay5 = pg.prop(H, delay5)

				sigma0 = pg.sigma_eq(spin_system)	# init
				sigma1 = pg.Ixpuls(spin_system, sigma0, self.sim_experiment.obs_iso, 90.0)		# apply ideal 90-degree pulse
				sigma0 = pg.evolve(sigma1, Udelay1)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP1
				sigma0 = pg.evolve(sigma1, Udelay2)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP2
				sigma0 = pg.evolve(sigma1, Udelay3)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP3
				sigma0 = pg.evolve(sigma1, Udelay4)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP4
				sigma0 = pg.evolve(sigma1, Udelay5)

				# acquire
				mx = pg.TTable1D(ACQ.table(sigma0))

				# binning to remove degenerate peaks
				outf, outa, outp = self.binning_code(mx, self.sim_experiment.b0, spin_system, self.sim_experiment.obs_iso, self.sim_experiment.tolppm, self.sim_experiment.tolpha, self.sim_experiment.ppmlo, self.sim_experiment.ppmhi, self.sim_experiment.RF_OFFSET)

				metab_calib180 = self.apply_metab_properties(str(self.calibrationMetaboliteComboBox.currentText()), A_180, outf, outa, outp, insysfile)
				lb = 15 if (self.sim_experiment.b0 == 297.2 or self.sim_experiment.b0 == 400.2) else 6
				f, spectra = metab_calib180.getSpec(TE, self.sim_experiment.b0, self.sim_experiment.getTime(), 0, 1, 0, 0, lb, self.sim_experiment.getFs())

				A180_calibration_data.append(np.real(spectra)[np.argmax(np.abs(spectra))])

			# Fit a logistic function
			initial_guess_A = np.amin(A180_calibration_data)
			initial_guess_B = np.abs(np.amax(A180_calibration_data)-np.amin(A180_calibration_data))/np.abs(A180_calibration_data[np.argmax(A180_calibration_data)]-A180_calibration_data[np.argmin(A180_calibration_data)])
			initial_guess_K = np.amax(A180_calibration_data)
			initial_guess_s = self.sim_experiment.A_180s[np.argmin(A180_calibration_data)]
			fit_x = self.sim_experiment.A_180s
			fit_y = np.pad(A180_calibration_data[np.argmin(A180_calibration_data):], (np.size(A180_calibration_data[:np.argmin(A180_calibration_data)]), 0), 'constant', constant_values=(np.amin(A180_calibration_data),0))
			params, params_covariance = sp.optimize.curve_fit(self.logs_func, fit_x, fit_y, p0=[initial_guess_A,initial_guess_B,initial_guess_K, initial_guess_s])#, bounds=([-np.inf, 0, 0, 0], [0, np.inf, np.inf, self.sim_experiment.A_180s[-1]]))
			self.simConsole.append('       | A + (K - A)/(1 + np.exp(-B*x)): ')
			self.simConsole.append('       | A0 = ' + str(initial_guess_A) + ', B0 = ' + str(initial_guess_B) + ', K0 = ' + str(initial_guess_K) + ', s0 = ' + str(initial_guess_s))
			self.simConsole.append('       | A = ' + str(params[0]) + ', B = ' + str(params[1]) + ', K = ' + str(params[2]) + ', s = ' + str(params[3]))
			A180_calibration_data_init   = self.logs_func(self.sim_experiment.A_180s, initial_guess_A, initial_guess_B, initial_guess_K, initial_guess_s)
			A180_calibration_data_fitted = self.logs_func(self.sim_experiment.A_180s, params[0], params[1], params[2], params[3])

			# save results
			if self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_180s[-1]/6) < self.sim_experiment.A_180s[-1]:
				self.sim_experiment.A_180 = self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_180s[-1]/6) # pad to be sure of adiabicity
			else:
				self.sim_experiment.A_180 = self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()]

			plt.figure()
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data, '.', color='blue')
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data_init, color='green')
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data_fitted, color='red')
			plt.plot(self.sim_experiment.A_180, A180_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_180s)-self.sim_experiment.A_180)).argmin()], 'x', color='black')
			plt.text(self.sim_experiment.A_180, A180_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_180s)-self.sim_experiment.A_180)).argmin()]*1.2, str(self.sim_experiment.A_180), rotation='vertical', color='black')
			plt.title('180-degree calibration sLASER (w/ ideal 90, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			plt.xlabel('A_180 [mT]')
			plt.ylabel('Signal Intensity'); plt.ylim([np.amin(A180_calibration_data) + 0.25*np.amin(A180_calibration_data),np.amax(A180_calibration_data) + 0.25*np.amax(A180_calibration_data)])
			plt.savefig(self.save_dir_sim + '180-DEGREE-CALIBRATION' + '_' + str(self.sim_experiment.PULSE_180_LENGTH) + 'secAFP.pdf')
			plt.close()

			# update progress bar
			self.simProgressBar.setValue(1)

			self.simConsole.append('       | CALIBRATED 180 AFP AMPLITUDE: ' + str(self.sim_experiment.A_180))

			self.slaser_build180(self.sim_experiment.inpulse180file, self.sim_experiment.A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, True, 'bruker')

			# Run the 90 calibration
			self.simConsole.append('\n2. 90-degree calibration sLASER (w/ AFP 180, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			A90_calibration_data = []

			for A_90 in self.sim_experiment.A_90s:

				TE = self.sim_experiment.TE * 1E-3
				TE1 = self.sim_experiment.TE1 * 1E-3
				TE2 = self.sim_experiment.TE2 * 1E-3

				# build 90 degree pulse
				A_90, pulse90, pulse_dur_90, peak_to_end_90, Ureal90 = self.slaser_build90(self.sim_experiment.inpulse90file, A_90, self.sim_experiment.PULSE_90_LENGTH, self.sim_experiment.getGyratio(), spin_system, False, 'bruker')

				# build 180 degree pulse
				A_180, pulse180, pulse_dur_180, Ureal180 = self.slaser_build180(self.sim_experiment.inpulse180file, A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, False, 'bruker')

				H = pg.Hcs(spin_system) + pg.HJ(spin_system)
				D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
				ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
				ACQ = ac

				delay1 = TE1/2.0 - pulse_dur_90/2.0 - pulse_dur_180/2.0
				delay2 = TE1/2.0 + TE2/2.0 - pulse_dur_180
				delay3 = TE2 - pulse_dur_180
				delay4 = delay2
				delay5 = TE1/2.0 - pulse_dur_180 + self.sim_experiment.DigShift

				# TE_fill = TE - 2.*TE1 - 2.*TE2
				# delay1 = TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 - pulse_dur_90/2.0
				# delay2 = TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay3 = TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay4 = TE2/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0
				# delay5 = TE1/2.0 - pulse_dur_180/2.0 + TE_fill/8.0 + self.sim_experiment.DigShift

				Udelay1 = pg.prop(H, delay1)
				Udelay2 = pg.prop(H, delay2)
				Udelay3 = pg.prop(H, delay3)
				Udelay4 = pg.prop(H, delay4)
				Udelay5 = pg.prop(H, delay5)

				sigma0 = pg.sigma_eq(spin_system)	# init
				sigma1 = Ureal90.evolve(sigma0)		# apply 90-degree pulse
				sigma0 = pg.evolve(sigma1, Udelay1)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP1
				sigma0 = pg.evolve(sigma1, Udelay2)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP2
				sigma0 = pg.evolve(sigma1, Udelay3)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP3
				sigma0 = pg.evolve(sigma1, Udelay4)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP4
				sigma0 = pg.evolve(sigma1, Udelay5)

				# acquire
				mx = pg.TTable1D(ACQ.table(sigma0))

				# binning to remove degenerate peaks
				outf, outa, outp = self.binning_code(mx, self.sim_experiment.b0, spin_system, self.sim_experiment.obs_iso, self.sim_experiment.tolppm, self.sim_experiment.tolpha, self.sim_experiment.ppmlo, self.sim_experiment.ppmhi, self.sim_experiment.RF_OFFSET)

				metab_calib90 = self.apply_metab_properties(str(self.calibrationMetaboliteComboBox.currentText()), A_90, outf, outa, outp, insysfile)
				lb = 15 if (self.sim_experiment.b0 == 297.2 or self.sim_experiment.b0 == 400.2) else 6
				f, spectra = metab_calib90.getSpec(TE, self.sim_experiment.b0, self.sim_experiment.getTime(), 0, 1, 0, 0, lb, self.sim_experiment.getFs())

				A90_calibration_data.append(np.real(spectra)[np.argmax(np.abs(spectra))])

			# Fit a sine function
			fitted_sine     = self.fit_sin(self.sim_experiment.A_90s, A90_calibration_data)
			params          = [fitted_sine[key] for key in fitted_sine.keys()]
			self.simConsole.append('       | A * sin(w*x + p) + c:')
			self.simConsole.append('       | A = ' + str(params[0]) + ', w = ' + str(params[1]) + ', p = ' + str(params[2]) + ', c = ' + str(params[3]))
			A90_calibration_data_fitted = self.sinfunc(self.sim_experiment.A_90s, params[0], params[1], params[2], params[3])

			# save results
			self.sim_experiment.A_90 = self.sim_experiment.A_90s[np.argmax(A90_calibration_data_fitted)]
			plt.figure()
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data, '.', color='blue')
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data_fitted, color='red')
			plt.plot(self.sim_experiment.A_90, A90_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_90s)-self.sim_experiment.A_90)).argmin()], 'x', color='black')
			plt.text(self.sim_experiment.A_90, A90_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_90s)-self.sim_experiment.A_90)).argmin()]*1.2, str(self.sim_experiment.A_90), rotation='vertical', color='black')
			plt.title('90-degree calibration sLASER (w/ AFP 180, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			plt.xlabel('A_90 [mT]')
			plt.ylabel('Signal Intensity'); plt.ylim([np.amin(A90_calibration_data) + 0.25*np.amin(A90_calibration_data),np.amax(A90_calibration_data) + 0.25*np.amax(A90_calibration_data)])
			plt.savefig(self.save_dir_sim + '90-DEGREE-CALIBRATION' + '_' + str(self.sim_experiment.PULSE_180_LENGTH) + 'secAFP.pdf')
			plt.close()

			self.simConsole.append('       | CALIBRATED 90 SLR AMPLITUDE: ' + str(self.sim_experiment.A_90))

			self.slaser_build90(self.sim_experiment.inpulse90file, self.sim_experiment.A_90, self.sim_experiment.PULSE_90_LENGTH, self.sim_experiment.getGyratio(), spin_system, True, 'bruker')

		elif self.sLASERradioButton.isChecked():
			
			# Create a new simulations results file
			self.save_dir_sim = 'pints/experiments/sLASER_sim_siemens_' \
								+ datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S") \
								+ '_TE' + str(self.sim_experiment.TE)
			
			if self.T7Button.isChecked():
				self.save_dir_sim += '_7T/'
			elif self.T3Button.isChecked():
				self.save_dir_sim += '_3T/'
			elif self.T9Button.isChecked():
				self.save_dir_sim += '_9T/'

			if not os.path.exists(self.save_dir_sim):
				os.makedirs(self.save_dir_sim)
			self.sim_results = open(self.save_dir_sim + 'sLASER_sim_results.txt', 'w')

			# Write parameters in file
			self.sim_results.write(';PINTS for FITMAN Simulation Output\n')
			self.sim_results.write(';Experiment Information\n')
			self.sim_results.write(';---------------------------------------------------------------------------\n')
			self.sim_results.write(';Name: ' + self.sim_experiment.name + '\n')
			self.sim_results.write(';Created: ' + self.sim_experiment.date + '\n')
			self.sim_results.write(';Comment: ' + self.sim_experiment.description + '\n')
			self.sim_results.write(';PI: ' + self.sim_experiment.author + '\n')
			self.sim_results.write(';b0: ' + str(self.sim_experiment.b0) + '\n')
			self.sim_results.write(';' + str(int(np.size(self.insysfiles)+int(np.size(self.macromolecules)))) + ' Metabolites: ' + str(self.insysfiles).replace('[','').replace(']','').replace('\'','').replace('.sys',''))
			if self.macroIncludeButton.isChecked(): self.sim_results.write(', ' + str(self.macromolecules).replace('[','').replace(']','').replace('\'',''))
			self.sim_results.write('\n')

			self.sim_results.write(';Simulation Results\n')
			self.sim_results.write(';---------------------------------------------------------------------------\n')

			# Set up for calibration experiments
			insysfile = str(self.calibrationMetaboliteComboBox.currentText()) + '.sys'
			if self.sim_experiment.b0 == 123.3:
				insysfile = 'pints/metabolites/3T_' + insysfile
			elif self.sim_experiment.b0 == 297.2:
				insysfile = 'pints/metabolites/7T_' + insysfile
			elif self.sim_experiment.b0 == 400.2:
				insysfile = 'pints/metabolites/9.4T_' + insysfile
			print('')
			print(insysfile)

			spin_system = pg.spin_system()
			spin_system.read(insysfile)
			for i in range(spin_system.spins()):
				spin_system.PPM(i, spin_system.PPM(i) - self.sim_experiment.RF_OFFSET)

			# Run the 180 calibration
			self.simConsole.append('\n1. 180-degree calibration sLASER (w/ ideal 90, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			A180_calibration_data = []

			for A_180 in self.sim_experiment.A_180s:

				TE = self.sim_experiment.TE
				TE1 = float((TE * 0.31) / 1000.0)
				TE3 = float((TE * 0.31) / 1000.0)
				TE2 = float(TE/1000.0 - TE1 - TE3)
				TE_fill = TE/1000.0 - TE1 - TE2 - TE3

				# build 90 degree pulse
				peak_to_end_90 = 0

				# build 180 degree pulse
				A_180, pulse180, pulse_dur_180, Ureal180 = self.slaser_build180(self.sim_experiment.inpulse180file, A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, False, 'siemens')

				H = pg.Hcs(spin_system) + pg.HJ(spin_system)
				D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
				ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
				ACQ = ac

				delay1 = TE1/2.0 + TE_fill/8.0                         - pulse_dur_180/2.0 - peak_to_end_90
				delay2 = TE1/2.0 + TE_fill/8.0 + TE2/4.0 + TE_fill/8.0 - pulse_dur_180
				delay3 = TE2/4.0 + TE_fill/8.0 + TE2/4.0 + TE_fill/8.0 - pulse_dur_180
				delay4 = TE2/4.0 + TE_fill/8.0 + TE3/2.0 + TE_fill/8.0 - pulse_dur_180
				delay5 = TE3/2.0 + TE_fill/8.0                         - pulse_dur_180/2.0

				Udelay1 = pg.prop(H, delay1)
				Udelay2 = pg.prop(H, delay2)
				Udelay3 = pg.prop(H, delay3)
				Udelay4 = pg.prop(H, delay4)
				Udelay5 = pg.prop(H, delay5)

				sigma0 = pg.sigma_eq(spin_system)	# init
				sigma1 = pg.Ixpuls(spin_system, sigma0, self.sim_experiment.obs_iso, 90.0)		# apply ideal 90-degree pulse
				sigma0 = pg.evolve(sigma1, Udelay1)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP1
				sigma0 = pg.evolve(sigma1, Udelay2)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP2
				sigma0 = pg.evolve(sigma1, Udelay3)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP3
				sigma0 = pg.evolve(sigma1, Udelay4)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP4
				sigma0 = pg.evolve(sigma1, Udelay5)

				# acquire
				mx = pg.TTable1D(ACQ.table(sigma0))

				# binning to remove degenerate peaks
				outf, outa, outp = self.binning_code(mx, self.sim_experiment.b0, spin_system, self.sim_experiment.obs_iso, self.sim_experiment.tolppm, self.sim_experiment.tolpha, self.sim_experiment.ppmlo, self.sim_experiment.ppmhi, self.sim_experiment.RF_OFFSET)

				metab_calib180 = self.apply_metab_properties(str(self.calibrationMetaboliteComboBox.currentText()), A_180, outf, outa, outp, insysfile)
				lb = 15 if (self.sim_experiment.b0 == 297.2 or self.sim_experiment.b0 == 400.2) else 6
				f, spectra = metab_calib180.getSpec(TE, self.sim_experiment.b0, self.sim_experiment.getTime(), 0, 1, 0, 0, lb, self.sim_experiment.getFs())

				A180_calibration_data.append(np.real(spectra)[np.argmax(np.abs(spectra))])

			# Fit a logistic function
			initial_guess_A = np.amin(A180_calibration_data)
			initial_guess_B = np.abs(np.amax(A180_calibration_data)-np.amin(A180_calibration_data))/np.abs(A180_calibration_data[np.argmax(A180_calibration_data)]-A180_calibration_data[np.argmin(A180_calibration_data)])
			initial_guess_K = np.amax(A180_calibration_data)
			initial_guess_s = self.sim_experiment.A_180s[np.argmin(A180_calibration_data)]
			fit_x = self.sim_experiment.A_180s
			fit_y = np.pad(A180_calibration_data[np.argmin(A180_calibration_data):], (np.size(A180_calibration_data[:np.argmin(A180_calibration_data)]), 0), 'constant', constant_values=(np.amin(A180_calibration_data),0))
			params, params_covariance = sp.optimize.curve_fit(self.logs_func, fit_x, fit_y, p0=[initial_guess_A,initial_guess_B,initial_guess_K, initial_guess_s])#, bounds=([-np.inf, 0, 0, 0], [0, np.inf, np.inf, self.sim_experiment.A_180s[-1]]))
			self.simConsole.append('       | A + (K - A)/(1 + np.exp(-B*x)): ')
			self.simConsole.append('       | A0 = ' + str(initial_guess_A) + ', B0 = ' + str(initial_guess_B) + ', K0 = ' + str(initial_guess_K) + ', s0 = ' + str(initial_guess_s))
			self.simConsole.append('       | A = ' + str(params[0]) + ', B = ' + str(params[1]) + ', K = ' + str(params[2]) + ', s = ' + str(params[3]))
			A180_calibration_data_init   = self.logs_func(self.sim_experiment.A_180s, initial_guess_A, initial_guess_B, initial_guess_K, initial_guess_s)
			A180_calibration_data_fitted = self.logs_func(self.sim_experiment.A_180s, params[0], params[1], params[2], params[3])

			# save results
			if self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_180s[-1]/6) < self.sim_experiment.A_180s[-1]:
				self.sim_experiment.A_180 = self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_180s[-1]/6) # pad to be sure of adiabicity
			else:
				self.sim_experiment.A_180 = self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()]
			plt.figure()
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data, '.', color='blue')
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data_init, color='green')
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data_fitted, color='red')
			plt.plot(self.sim_experiment.A_180, A180_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_180s)-self.sim_experiment.A_180)).argmin()], 'x', color='black')
			plt.text(self.sim_experiment.A_180, A180_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_180s)-self.sim_experiment.A_180)).argmin()]*1.2, str(self.sim_experiment.A_180), rotation='vertical', color='black')
			plt.title('180-degree calibration sLASER (w/ ideal 90, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			plt.xlabel('A_180 [mT]')
			plt.ylabel('Signal Intensity'); plt.ylim([np.amin(A180_calibration_data) + 0.25*np.amin(A180_calibration_data),np.amax(A180_calibration_data) + 0.25*np.amax(A180_calibration_data)])
			plt.savefig(self.save_dir_sim + '180-DEGREE-CALIBRATION' + '_' + str(self.sim_experiment.PULSE_180_LENGTH) + 'secAFP.pdf')
			plt.close()

			# update progress bar
			self.simProgressBar.setValue(1)

			self.simConsole.append('       | CALIBRATED 180 AFP AMPLITUDE: ' + str(self.sim_experiment.A_180))

			self.slaser_build180(self.sim_experiment.inpulse180file, self.sim_experiment.A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, True, 'siemens')

			# Run the 90 calibration
			self.simConsole.append('\n2. 90-degree calibration sLASER (w/ AFP 180, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			A90_calibration_data = []

			for A_90 in self.sim_experiment.A_90s:
				TE = self.sim_experiment.TE
				TE1 = float((TE * 0.31) / 1000.0)
				TE3 = float((TE * 0.31) / 1000.0)
				TE2 = float(TE/1000.0 - TE1 - TE3)
				TE_fill = TE/1000.0 - TE1 - TE2 - TE3

				# build 90 degree pulse
				A_90, pulse90, pulse_dur_90, peak_to_end_90, Ureal90 = self.slaser_build90(self.sim_experiment.inpulse90file, A_90, self.sim_experiment.PULSE_90_LENGTH, self.sim_experiment.getGyratio(), spin_system, False, 'siemens')

				# build 180 degree pulse
				A_180, pulse180, pulse_dur_180, Ureal180 = self.slaser_build180(self.sim_experiment.inpulse180file, self.sim_experiment.A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, False, 'siemens')

				H = pg.Hcs(spin_system) + pg.HJ(spin_system)
				D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
				ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
				ACQ = ac

				delay1 = TE1/2.0 + TE_fill/8.0                         - pulse_dur_180/2.0 - peak_to_end_90
				delay2 = TE1/2.0 + TE_fill/8.0 + TE2/4.0 + TE_fill/8.0 - pulse_dur_180
				delay3 = TE2/4.0 + TE_fill/8.0 + TE2/4.0 + TE_fill/8.0 - pulse_dur_180
				delay4 = TE2/4.0 + TE_fill/8.0 + TE3/2.0 + TE_fill/8.0 - pulse_dur_180
				delay5 = TE3/2.0 + TE_fill/8.0                         - pulse_dur_180/2.0

				Udelay1 = pg.prop(H, delay1)
				Udelay2 = pg.prop(H, delay2)
				Udelay3 = pg.prop(H, delay3)
				Udelay4 = pg.prop(H, delay4)
				Udelay5 = pg.prop(H, delay5)

				sigma0 = pg.sigma_eq(spin_system)							# init
				sigma1 = Ureal90.evolve(sigma0)								# apply 90-degree pulse
				sigma0 = pg.evolve(sigma1, Udelay1)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP1
				sigma0 = pg.evolve(sigma1, Udelay2)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP2
				sigma0 = pg.evolve(sigma1, Udelay3)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP3
				sigma0 = pg.evolve(sigma1, Udelay4)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP4
				sigma0 = pg.evolve(sigma1, Udelay5)

				# acquire
				mx = pg.TTable1D(ACQ.table(sigma0))

				# binning to remove degenerate peaks
				outf, outa, outp = self.binning_code(mx, self.sim_experiment.b0, spin_system, self.sim_experiment.obs_iso, self.sim_experiment.tolppm, self.sim_experiment.tolpha, self.sim_experiment.ppmlo, self.sim_experiment.ppmhi, self.sim_experiment.RF_OFFSET)

				metab_calib90 = self.apply_metab_properties(str(self.calibrationMetaboliteComboBox.currentText()), A_90, outf, outa, outp, insysfile)
				lb = 15 if (self.sim_experiment.b0 == 297.2 or self.sim_experiment.b0 == 400.2) else 6
				f, spectra = metab_calib90.getSpec(TE, self.sim_experiment.b0, self.sim_experiment.getTime(), 0, 1, 0, 0, lb, self.sim_experiment.getFs())

				A90_calibration_data.append(np.real(spectra)[np.argmax(np.abs(spectra))])

			# Fit a sine function
			fitted_sine     = self.fit_sin(self.sim_experiment.A_90s, A90_calibration_data)
			params          = [fitted_sine[key] for key in fitted_sine.keys()]
			self.simConsole.append('       | A * sin(w*x + p) + c:')
			self.simConsole.append('       | A = ' + str(params[0]) + ', w = ' + str(params[1]) + ', p = ' + str(params[2]) + ', c = ' + str(params[3]))
			A90_calibration_data_fitted = self.sinfunc(self.sim_experiment.A_90s, params[0], params[1], params[2], params[3])

			# save results
			self.sim_experiment.A_90 = self.sim_experiment.A_90s[np.argmax(A90_calibration_data_fitted)]
			plt.figure()
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data, '.', color='blue')
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data_fitted, color='red')
			plt.plot(self.sim_experiment.A_90, A90_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_90s)-self.sim_experiment.A_90)).argmin()], 'x', color='black')
			plt.text(self.sim_experiment.A_90, A90_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_90s)-self.sim_experiment.A_90)).argmin()]*1.2, str(self.sim_experiment.A_90), rotation='vertical', color='black')
			plt.title('90-degree calibration sLASER (w/ AFP 180, ' + str(self.calibrationMetaboliteComboBox.currentText()) + ') experiment')
			plt.xlabel('A_90 [mT]')
			plt.ylabel('Signal Intensity'); plt.ylim([np.amin(A90_calibration_data) + 0.25*np.amin(A90_calibration_data),np.amax(A90_calibration_data) + 0.25*np.amax(A90_calibration_data)])
			plt.savefig(self.save_dir_sim + '90-DEGREE-CALIBRATION' + '_' + str(self.sim_experiment.PULSE_180_LENGTH) + 'secAFP.pdf')
			plt.close()

			self.simConsole.append('       | CALIBRATED 90 SLR AMPLITUDE: ' + str(self.sim_experiment.A_90))

			self.slaser_build90(self.sim_experiment.inpulse90file, self.sim_experiment.A_90, self.sim_experiment.PULSE_90_LENGTH, self.sim_experiment.getGyratio(), spin_system, True, 'siemens')

		elif self.LASERradioButton.isChecked():
			
			# Create a new simulations results file
			self.save_dir_sim = 'pints/experiments/LASER_sim_' \
								+ datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S") \
								+ '_TE' + str(self.sim_experiment.TE)
			
			if self.T7Button.isChecked():
				self.save_dir_sim += '_7T/'
			elif self.T3Button.isChecked():
				self.save_dir_sim += '_3T/'
			elif self.T9Button.isChecked():
				self.save_dir_sim += '_9T/'

			if not os.path.exists(self.save_dir_sim):
				os.makedirs(self.save_dir_sim)
			self.sim_results = open(self.save_dir_sim + 'LASER_sim_results.txt', 'w')

			# Write parameters in file
			self.sim_results.write(';PINTS for FITMAN Simulation Output\n')
			self.sim_results.write(';Experiment Information\n')
			self.sim_results.write(';---------------------------------------------------------------------------\n')
			self.sim_results.write(';Name: ' + self.sim_experiment.name + '\n')
			self.sim_results.write(';Created: ' + self.sim_experiment.date + '\n')
			self.sim_results.write(';Comment: ' + self.sim_experiment.description + '\n')
			self.sim_results.write(';PI: ' + self.sim_experiment.author + '\n')
			self.sim_results.write(';b0: ' + str(self.sim_experiment.b0) + '\n')
			self.sim_results.write(';' + str(int(np.size(self.insysfiles)+int(np.size(self.macromolecules)))) + ' Metabolites: ' + str(self.insysfiles).replace('[','').replace(']','').replace('\'','').replace('.sys',''))
			if self.macroIncludeButton.isChecked(): self.sim_results.write(', ' + str(self.macromolecules).replace('[','').replace(']','').replace('\'',''))
			self.sim_results.write('\n')

			self.sim_results.write(';Simulation Results\n')
			self.sim_results.write(';---------------------------------------------------------------------------\n')

			# Set up for calibration experiments
			insysfile = str(self.calibrationMetaboliteComboBox_laser.currentText()) + '.sys'
			if self.sim_experiment.b0 == 123.3:
				insysfile = 'pints/metabolites/3T_' + insysfile
			elif self.sim_experiment.b0 == 297.2:
				insysfile = 'pints/metabolites/7T_' + insysfile
			elif self.sim_experiment.b0 == 400.2:
				insysfile = 'pints/metabolites/9.4T_' + insysfile

			spin_system = pg.spin_system()
			spin_system.read(insysfile)
			for i in range(spin_system.spins()):
				spin_system.PPM(i, spin_system.PPM(i) - self.sim_experiment.RF_OFFSET)

			# Run the 180 calibration
			self.simConsole.append('\n1. 180-degree calibration LASER (w/ ideal 90, ' + str(self.calibrationMetaboliteComboBox_laser.currentText())  + ') experiment')
			A180_calibration_data = []

			for A_180 in self.sim_experiment.A_180s:

				# build 90 degree pulse
				pulse_dur_90 = 0

				# build 180 degree pulse
				A_180, pulse180, pulse_dur_180, Ureal180 = self.laser_buildafp(self.sim_experiment.inpulse180file, A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, False)

				# calculate pulse timings
				ROF1 = 100E-6   #sec
				ROF2 = 10E-6    #sec
				TCRUSH1 = 0.0008 #sec
				TCRUSH2 = 0.0008 #sec

				ss_grad_rfDelayFront = TCRUSH1 - ROF1
				ss_grad_rfDelayBack  = TCRUSH2 - ROF2
				ro_grad_atDelayFront = 0
				ro_grad_atDelayBack  = 0

				TE  = self.sim_experiment.TE / 1000.
				ipd = (TE - pulse_dur_90 \
						  - 6*(ss_grad_rfDelayFront + pulse_dur_180 + ss_grad_rfDelayBack) \
						  - ro_grad_atDelayFront) / 12

				# initialize acquisition
				H = pg.Hcs(spin_system) + pg.HJ(spin_system)
				D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
				ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
				ACQ = ac

				delay1 = ipd + ss_grad_rfDelayFront
				delay2 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay3 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay4 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay5 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay6 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay7 = ss_grad_rfDelayBack + ipd + ro_grad_atDelayFront

				Udelay1 = pg.prop(H, delay1)
				Udelay2 = pg.prop(H, delay2)
				Udelay3 = pg.prop(H, delay3)
				Udelay4 = pg.prop(H, delay4)
				Udelay5 = pg.prop(H, delay5)
				Udelay6 = pg.prop(H, delay6)
				Udelay7 = pg.prop(H, delay7)

				sigma0 = pg.sigma_eq(spin_system)	# init
				sigma1 = pg.Ixpuls(spin_system, sigma0, self.sim_experiment.obs_iso, 90.0)		# apply ideal 90-degree pulse
				sigma0 = pg.evolve(sigma1, Udelay1)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP1
				sigma0 = pg.evolve(sigma1, Udelay2)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP2
				sigma0 = pg.evolve(sigma1, Udelay3)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP3
				sigma0 = pg.evolve(sigma1, Udelay4)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP4
				sigma0 = pg.evolve(sigma1, Udelay5)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP5
				sigma0 = pg.evolve(sigma1, Udelay6)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP6
				sigma0 = pg.evolve(sigma1, Udelay7)

				# acquire
				mx = pg.TTable1D(ACQ.table(sigma0))

				# binning to remove degenerate peaks
				outf, outa, outp = self.binning_code(mx, self.sim_experiment.b0, spin_system, self.sim_experiment.obs_iso, self.sim_experiment.tolppm, self.sim_experiment.tolpha, self.sim_experiment.ppmlo, self.sim_experiment.ppmhi, self.sim_experiment.RF_OFFSET)

				metab_calib180 = self.apply_metab_properties(str(self.calibrationMetaboliteComboBox_laser.currentText()), A_180, outf, outa, outp, insysfile)
				lb = 15 if (self.sim_experiment.b0 == 297.2 or self.sim_experiment.b0 == 400.2) else 6
				f, spectra = metab_calib180.getSpec(TE, self.sim_experiment.b0, self.sim_experiment.getTime(), 0, 1, 0, 0, lb, self.sim_experiment.getFs())

				A180_calibration_data.append(np.real(spectra)[np.argmax(np.abs(spectra))])

			# Fit a logistic function
			initial_guess_A = np.amin(A180_calibration_data)
			initial_guess_B = np.abs(np.amax(A180_calibration_data)-np.amin(A180_calibration_data))/np.abs(A180_calibration_data[np.argmax(A180_calibration_data)]-A180_calibration_data[np.argmin(A180_calibration_data)])
			initial_guess_K = np.amax(A180_calibration_data)
			initial_guess_s = self.sim_experiment.A_180s[np.argmin(A180_calibration_data)]
			fit_x = self.sim_experiment.A_180s
			fit_y = np.pad(A180_calibration_data[np.argmin(A180_calibration_data):], (np.size(A180_calibration_data[:np.argmin(A180_calibration_data)]), 0), 'constant', constant_values=(np.amin(A180_calibration_data),0))
			params, params_covariance = sp.optimize.curve_fit(self.logs_func, fit_x, fit_y, p0=[initial_guess_A,initial_guess_B,initial_guess_K, initial_guess_s])#, bounds=([-np.inf, 0, 0, 0], [0, np.inf, np.inf, self.sim_experiment.A_180s[-1]]))
			self.simConsole.append('       | A + (K - A)/(1 + np.exp(-B*x)): ')
			self.simConsole.append('       | A0 = ' + str(initial_guess_A) + ', B0 = ' + str(initial_guess_B) + ', K0 = ' + str(initial_guess_K) + ', s0 = ' + str(initial_guess_s))
			self.simConsole.append('       | A = ' + str(params[0]) + ', B = ' + str(params[1]) + ', K = ' + str(params[2]) + ', s = ' + str(params[3]))
			A180_calibration_data_init   = self.logs_func(self.sim_experiment.A_180s, initial_guess_A, initial_guess_B, initial_guess_K, initial_guess_s)
			A180_calibration_data_fitted = self.logs_func(self.sim_experiment.A_180s, params[0], params[1], params[2], params[3])

			# save results
			if self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_180s[-1]/6) < self.sim_experiment.A_180s[-1]:
				self.sim_experiment.A_180 = self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_180s[-1]/6) # pad to be sure of adiabicity
			else:
				self.sim_experiment.A_180 = self.sim_experiment.A_180s[(np.abs(np.asarray(fit_y)-params[2])).argmin()]
			plt.figure()
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data, '.', color='blue')
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data_init, color='green')
			plt.plot(self.sim_experiment.A_180s, A180_calibration_data_fitted, color='red')
			plt.plot(self.sim_experiment.A_180, A180_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_180s)-self.sim_experiment.A_180)).argmin()], 'x', color='black')
			plt.text(self.sim_experiment.A_180, A180_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_180s)-self.sim_experiment.A_180)).argmin()]*1.2, str(self.sim_experiment.A_180), rotation='vertical', color='black')
			plt.title('180-degree calibration sLASER (w/ ideal 90, ' + str(self.calibrationMetaboliteComboBox_laser.currentText()) + ') experiment')
			plt.xlabel('A_180 [mT]')
			plt.ylabel('Signal Intensity'); plt.ylim([np.amin(A180_calibration_data) + 0.25*np.amin(A180_calibration_data),np.amax(A180_calibration_data) + 0.25*np.amax(A180_calibration_data)])
			plt.savefig(self.save_dir_sim + '180-DEGREE-CALIBRATION' + '_' + str(self.sim_experiment.PULSE_180_LENGTH) + 'secAFP.pdf')
			plt.close()

			# update progress bar
			self.simProgressBar.setValue(1)

			self.simConsole.append('       | CALIBRATED 180 AFP AMPLITUDE: ' + str(self.sim_experiment.A_180))

			self.laser_buildafp(self.sim_experiment.inpulse180file, self.sim_experiment.A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, True)

			# Run the 90 calibration
			self.simConsole.append('\n2. 90-degree calibration LASER (w/ AFP 180, ' + str(self.calibrationMetaboliteComboBox_laser.currentText()) + ') experiment')
			A90_calibration_data = []

			for A_90 in self.sim_experiment.A_90s:

				# build 90 degree pulse
				A_90, pulse90, pulse_dur_90, Ureal90 = self.laser_buildahp(self.sim_experiment.inpulse90file, A_90, self.sim_experiment.PULSE_90_LENGTH, self.sim_experiment.getGyratio(), spin_system, False)

				# build 180 degree pulse
				A_180, pulse180, pulse_dur_180, Ureal180 = self.laser_buildafp(self.sim_experiment.inpulse180file, self.sim_experiment.A_180, self.sim_experiment.PULSE_180_LENGTH, self.sim_experiment.getGyratio(), spin_system, False)

				# calculate pulse timings
				ROF1 = 100E-6   #sec
				ROF2 = 10E-6    #sec
				TCRUSH1 = 0.0008 #sec
				TCRUSH2 = 0.0008 #sec

				ss_grad_rfDelayFront = TCRUSH1 - ROF1
				ss_grad_rfDelayBack  = TCRUSH2 - ROF2
				ro_grad_atDelayFront = 0
				ro_grad_atDelayBack  = 0

				TE  = self.sim_experiment.TE / 1000.
				ipd = (TE - pulse_dur_90 \
						  - 6*(ss_grad_rfDelayFront + pulse_dur_180 + ss_grad_rfDelayBack) \
						  - ro_grad_atDelayFront) / 12

				# initialize acquisition
				H = pg.Hcs(spin_system) + pg.HJ(spin_system)
				D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
				ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
				ACQ = ac

				delay1 = ipd + ss_grad_rfDelayFront
				delay2 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay3 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay4 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay5 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay6 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
				delay7 = ss_grad_rfDelayBack + ipd + ro_grad_atDelayFront

				Udelay1 = pg.prop(H, delay1)
				Udelay2 = pg.prop(H, delay2)
				Udelay3 = pg.prop(H, delay3)
				Udelay4 = pg.prop(H, delay4)
				Udelay5 = pg.prop(H, delay5)
				Udelay6 = pg.prop(H, delay6)
				Udelay7 = pg.prop(H, delay7)

				sigma0 = pg.sigma_eq(spin_system)	# init
				sigma1 = Ureal90.evolve(sigma0)		# apply 90-degree pulse
				sigma0 = pg.evolve(sigma1, Udelay1)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP1
				sigma0 = pg.evolve(sigma1, Udelay2)
				sigma1 = Ureal180.evolve(sigma0)	# apply AFP2
				sigma0 = pg.evolve(sigma1, Udelay3)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP3
				sigma0 = pg.evolve(sigma1, Udelay4)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP4
				sigma0 = pg.evolve(sigma1, Udelay5)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP5
				sigma0 = pg.evolve(sigma1, Udelay6)
				sigma1 = Ureal180.evolve(sigma0) 	# apply AFP6
				sigma0 = pg.evolve(sigma1, Udelay7)

				# acquire
				mx = pg.TTable1D(ACQ.table(sigma0))

				# binning to remove degenerate peaks
				outf, outa, outp = self.binning_code(mx, self.sim_experiment.b0, spin_system, self.sim_experiment.obs_iso, self.sim_experiment.tolppm, self.sim_experiment.tolpha, self.sim_experiment.ppmlo, self.sim_experiment.ppmhi, self.sim_experiment.RF_OFFSET)

				metab_calib90 = self.apply_metab_properties(str(self.calibrationMetaboliteComboBox_laser.currentText()), A_90, outf, outa, outp, insysfile)
				lb = 15 if (self.sim_experiment.b0 == 297.2 or self.sim_experiment.b0 == 400.2) else 6
				f, spectra = metab_calib90.getSpec(TE, self.sim_experiment.b0, self.sim_experiment.getTime(), 0, 1, 0, 0, lb, self.sim_experiment.getFs())

				A90_calibration_data.append(np.real(spectra)[np.argmax(np.abs(spectra))])

			# Fit a logistic function
			initial_guess_A = 0 # np.amin(A90_calibration_data)
			initial_guess_B = np.abs(np.amax(A90_calibration_data)-np.amin(A90_calibration_data))/np.abs(A90_calibration_data[np.argmax(A90_calibration_data)]-A90_calibration_data[np.argmin(A90_calibration_data)])
			initial_guess_K = np.amax(A90_calibration_data)
			initial_guess_s = self.sim_experiment.A_90s[np.argmin(A90_calibration_data)]
			fit_x = self.sim_experiment.A_90s
			fit_y = np.pad(A90_calibration_data[np.argmin(A90_calibration_data):], (np.size(A90_calibration_data[:np.argmin(A90_calibration_data)]), 0), 'constant', constant_values=(np.amin(A90_calibration_data),0))
			params, params_covariance = sp.optimize.curve_fit(self.logs_func, fit_x, fit_y, p0=[initial_guess_A,initial_guess_B,initial_guess_K, initial_guess_s])#, bounds=([-np.inf, 0, 0, 0], [0, np.inf, np.inf, self.sim_experiment.A_90s[-1]]))
			self.simConsole.append('       | A + (K - A)/(1 + np.exp(-B*x)): ')
			self.simConsole.append('       | A0 = ' + str(initial_guess_A) + ', B0 = ' + str(initial_guess_B) + ', K0 = ' + str(initial_guess_K) + ', s0 = ' + str(initial_guess_s))
			self.simConsole.append('       | A = ' + str(params[0]) + ', B = ' + str(params[1]) + ', K = ' + str(params[2]) + ', s = ' + str(params[3]))
			A90_calibration_data_init   = self.logs_func(self.sim_experiment.A_90s, initial_guess_A, initial_guess_B, initial_guess_K, initial_guess_s)
			A90_calibration_data_fitted = self.logs_func(self.sim_experiment.A_90s, params[0], params[1], params[2], params[3])

			# save results
			if self.sim_experiment.A_90s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_90s[-1]/6) < self.sim_experiment.A_90s[-1]:
				self.sim_experiment.A_90 = self.sim_experiment.A_90s[(np.abs(np.asarray(fit_y)-params[2])).argmin()] + int(self.sim_experiment.A_90s[-1]/6) # pad to be sure of adiabicity
			else:
				self.sim_experiment.A_90 = self.sim_experiment.A_90s[(np.abs(np.asarray(fit_y)-params[2])).argmin()]
			plt.figure()
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data, '.', color='blue')
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data_init, color='green')
			plt.plot(self.sim_experiment.A_90s, A90_calibration_data_fitted, color='red')
			plt.plot(self.sim_experiment.A_90, A90_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_90s)-self.sim_experiment.A_90)).argmin()], 'x', color='black')
			plt.text(self.sim_experiment.A_90, A90_calibration_data[(np.abs(np.asarray(self.sim_experiment.A_90s)-self.sim_experiment.A_90)).argmin()]*1.2, str(self.sim_experiment.A_90), rotation='vertical', color='black')
			plt.title('90-degree calibration LASER (w/ AFP 180, ' + str(self.calibrationMetaboliteComboBox_laser.currentText()) + ') experiment')
			plt.xlabel('A_90 [mT]')
			plt.ylabel('Signal Intensity'); plt.ylim([np.amin(A90_calibration_data) + 0.25*np.amin(A90_calibration_data),np.amax(A90_calibration_data) + 0.25*np.amax(A90_calibration_data)])
			plt.savefig(self.save_dir_sim + '90-DEGREE-CALIBRATION' + '_' + str(self.sim_experiment.PULSE_180_LENGTH) + 'secAFP.pdf')
			plt.close()

			self.simConsole.append('       | CALIBRATED 90 AHP AMPLITUDE: ' + str(self.sim_experiment.A_90))

			self.laser_buildahp(self.sim_experiment.inpulse90file, self.sim_experiment.A_90, self.sim_experiment.PULSE_90_LENGTH, self.sim_experiment.getGyratio(), spin_system, True)

		self.simProgressBar.setValue(2)

		self.simConsole.append('\n3. Density Matrix Simulations of Metabolites')

		self.thread_pool = self.tree()
		self.sim_objects = self.tree()

		self.youngest_thread = 0
		self.threads_processed = 0

		for (i, insysfile) in enumerate(self.insysfiles):

			# create a new thread to simulate metabolite
			simulation_thread = QtCore.QThread()
			self.thread_pool[i] = simulation_thread

			# create a new simulation
			simulation = MetaboliteSimulation(i, insysfile, self.sim_experiment)
			self.sim_objects[i] = simulation

			# events to output to console
			self.sim_objects[i].postToConsole.connect(self.postToConsole)

			# events to output results
			self.sim_objects[i].outputResults.connect(self.outputResults)

			# events to signal that simulation is complete
			self.sim_objects[i].finished.connect(self.simFinished)

			self.sim_objects[i].moveToThread(self.thread_pool[i])
			self.thread_pool[i].started.connect(self.sim_objects[i].simulate)

		self.startThreads()

	# ---- Methods for Simulation ---- #
	def tree(self): return defaultdict(self.tree)

	def startThreads(self):
		if self.youngest_thread + 3 > len(self.insysfiles)-1:
			for i in range(self.youngest_thread, len(self.insysfiles)):
				self.thread_pool[i].start()
		else:
			self.thread_pool[self.youngest_thread + 3].start()
			self.thread_pool[self.youngest_thread + 2].start()
			self.thread_pool[self.youngest_thread + 1].start()
			self.thread_pool[self.youngest_thread + 0].start()

	def simFinished(self, thread_num):
		self.thread_pool[thread_num].quit()
		self.threads_processed = self.threads_processed + 1
		self.simProgressBar.setValue(self.simProgressBar.value() + 1)

		if self.simProgressBar.value() == self.simProgressBar.maximum():
			
			if self.macroIncludeButton.isChecked():
				# add lipids
				for lipid in self.macromolecules:
					macromolecule = self.macromolecules_data[lipid]
					# name, TE, A_m, shift, line_type, lw, area, phase
					self.sim_results.write(macromolecule.name + '\t' + str(float(self.sim_experiment.TE)) \
						+ '\t' + str(float(macromolecule.A_m)) + '\t' + str(float(macromolecule.ppm[0])) \
						+ '\t' + str(macromolecule.line_type) + '\t' + str(macromolecule.lw) \
						+ '\t' + str(macromolecule.area[0]) + '\t' + str(macromolecule.phase[0]) + '\n')

			self.simConsole.append('  | Saving simulation data ...')
			time.sleep(2)
			self.sim_results.close()
			self.simConsole.append('  | Simulation finished.')

			# save console output to text file
			with open(self.save_dir_sim + 'console.txt', 'w') as console_output_file:
				console_output_file.write(str(self.simConsole.toPlainText()))

			self.confirmSimParamsButton.setEnabled(True)
			self.runSimulationButton.setEnabled(False)
		else:
			if self.threads_processed == 4:
				self.youngest_thread = self.youngest_thread + 4
				self.threads_processed = 0
				self.startThreads()

	def postToConsole(self, string):
		self.simConsole.append(string)

	def outputResults(self, metabolite):
		for (i, area) in enumerate(metabolite.area):
			self.sim_results.write(metabolite.name + '\t' + str(float(self.sim_experiment.TE)) + '\t' + str(float(metabolite.A_m)) + '\t' + str(float(metabolite.T2)) + '\t' + str(i) + '\t' + str(metabolite.ppm[i]) + '\t' + str(metabolite.area[i]) + '\t' + str(metabolite.phase[i]) + '\n')

	def binning_code(self, mx, b0, spin_system, obs_iso, tolppm, tolpha, ppmlo, ppmhi, rf_off):

		# BINNING
		# Note: Metabolite Peak Normalization and Blending

		# The transition tables calculated by the GAMMA density matrix simulations frequently contain a
		# large number of transitions caused by degenerate splittings and other processes. At the
		# conclusion of each simulation run a routine is called to extract lines from the transition table.
		# These lines are then normalized using a closed form calculation based on the number of spins.
		# To reduce the number of lines required for display, multiple lines are blended by binning them
		# together based on their PPM locations and phases. The following parameters are used to
		# customize these procedures:

		# Peak Search Range -- Low/High (PPM): the range in PPM that is searched for lines from the
		# metabolite simulation.

		# Peak Blending Tolerance (PPM and Degrees): the width of the bins (+/- in PPM and +/- in
		# PhaseDegrees) that are used to blend the lines in the simulation. Lines that are included in the
		# same bin are summed using complex addition based on Amplitude and Phase.

		field  = b0
		nspins = spin_system.spins()

		nlines = mx.size()

		tmp = pg.Isotope(obs_iso)
		obs_qn = tmp.qn()

		qnscale = 1.0
		for i in range(nspins):
			qnscale *= 2*spin_system.qn(i)+1
		qnscale = qnscale / (2.0 * (2.0*obs_qn+1))

		freqs = []
		outf = []
		outa = []
		outp = []
		nbin = 0
		found = False

		PI = 3.14159265358979323846
		RAD2DEG = 180.0/PI

		indx = mx.Sort(0,-1,0)

		for i in range(nlines):
			freqs.append(-1 * mx.Fr(indx[i])/(2.0*PI*field))

		for i in range(nlines):
			freq = freqs[i]
			if (freq > ppmlo) and (freq < ppmhi):
				val = mx.I(indx[i])
				tmpa = np.sqrt(val.real()**2 + val.imag()**2) / qnscale
				tmpp = -RAD2DEG * np.angle(val.real()+1j*val.imag())

			if nbin == 0:
				outf.append(freq)
				outa.append(tmpa)
				outp.append(tmpp)
				nbin += 1
			else:
				for k in range(nbin):
					if (freq >= outf[k]-tolppm) and (freq <= outf[k]+tolppm):
						if (tmpp >= outp[k]-tolpha) and (tmpp <= outp[k]+tolpha):
							ampsum   =  outa[k]+tmpa
							outf[k]  = (outa[k]*outf[k] + tmpa*freq)/ampsum
							outp[k]  = (outa[k]*outp[k] + tmpa*tmpp)/ampsum
							outa[k] +=  tmpa;
							found = True 
				if not found:
					outf.append(freq)
					outa.append(tmpa)
					outp.append(tmpp)
					nbin += 1
				found = False

		for i, item in enumerate(outf):
			outf[i] = item + rf_off
			outp[i] = outp[i] - 90.0

		return outf, outa, outp

	def apply_metab_properties(self, metab_name, var, outf, outa, outp, insysfile):

		metab = Metabolite()
		metab.name = metab_name
		metab.var = var

		for i in range(np.size(outf)):
			if outf[i] <= 5:
				metab.ppm.append(outf[i])
				metab.area.append(outa[i])
				metab.phase.append(-1.0*outp[i])

		insysfile = insysfile.replace('pints/metabolites/3T_', '')
		insysfile = insysfile.replace('pints/metabolites/7T_', '')
		insysfile = insysfile.replace('pints/metabolites/9.4T_', '')

		if insysfile == 'alanine.sys': #
			metab.A_m = 0.078
			metab.T2 = (87E-3)
		elif insysfile == 'aspartate.sys':
			metab.A_m = 0.117
			metab.T2 = (87E-3)
		elif insysfile == 'choline_1-CH2_2-CH2.sys': #
			metab.A_m = 0.165
			metab.T2 = (87E-3)
		elif insysfile == 'choline_N(CH3)3_a.sys' or insysfile == 'choline_N(CH3)3_b.sys': #
			metab.A_m = 0.165
			metab.T2 = (121E-3)
		elif insysfile == 'creatine_N(CH3).sys':
			metab.A_m = 0.296
			metab.T2 = (90E-3)
		elif insysfile == 'creatine_X.sys':
			metab.A_m = 0.296
			metab.T2 = (81E-3)
		elif insysfile == 'd-glucose-alpha.sys': #
			metab.A_m = 0.049
			metab.T2 = (87E-3)
		elif insysfile == 'd-glucose-beta.sys': #
			metab.A_m = 0.049
			metab.T2 = (87E-3)
		elif insysfile == 'eth.sys': #
			metab.A_m = 0.320
			metab.T2 = (87E-3)
		elif insysfile == 'gaba.sys': #
			metab.A_m = 0.155
			metab.T2 = (82E-3)
		elif insysfile == 'glutamate.sys':
			metab.A_m = 0.898
			metab.T2 = (88E-3)
		elif insysfile == 'glutamine.sys':
			metab.A_m = 0.427
			metab.T2 = (87E-3)
		elif insysfile == 'glutathione_cysteine.sys':
			metab.A_m = 0.194
			metab.T2 = (87E-3)
		elif insysfile == 'glutathione_glutamate.sys':
			metab.A_m = 0.194
			metab.T2 = (87E-3)
		elif insysfile == 'glutathione_glycine.sys':
			metab.A_m = 0.194
			metab.T2 = (87E-3)
		elif insysfile == 'glycine.sys':
			metab.A_m = 0.068
			metab.T2 = (87E-3)
		elif insysfile == 'gpc_7-CH2_8-CH2.sys': #
			metab.A_m = 0.097
			metab.T2 = (87E-3)
		elif insysfile == 'gpc_glycerol.sys': #
			metab.A_m = 0.097
			metab.T2 = (87E-3)
		elif insysfile == 'gpc_N(CH3)3_a.sys': #
			metab.A_m = 0.097
			metab.T2 = (121E-3)
		elif insysfile == 'gpc_N(CH3)3_b.sys': #
			metab.A_m = 0.097
			metab.T2 = (121E-3)
		elif insysfile == 'lactate.sys': #
			metab.A_m = 0.039
			metab.T2 = (87E-3)
		elif insysfile == 'myoinositol.sys':
			metab.A_m = 0.578
			metab.T2 = (87E-3)
		elif insysfile == 'naa_acetyl.sys':
			metab.A_m = 1.000
			metab.T2 = (130E-3)
		elif insysfile == 'naa_aspartate.sys':
			metab.A_m = 1.000
			metab.T2 = (69E-3)
		elif insysfile == 'naag_acetyl.sys':
			metab.A_m = 0.160
			metab.T2 = (130E-3)
		elif insysfile == 'naag_aspartyl.sys':
			metab.A_m = 0.160
			metab.T2 = (87E-3)
		elif insysfile == 'naag_glutamate.sys':
			metab.A_m = 0.160
			metab.T2 = (87E-3)
		elif insysfile == 'pcho_N(CH3)3_a.sys': #
			metab.A_m = 0.058
			metab.T2 = (121E-3)
		elif insysfile == 'pcho_N(CH3)3_b.sys': #
			metab.A_m = 0.058
			metab.T2 = (121E-3)
		elif insysfile == 'pcho_X.sys': #
			metab.A_m = 0.058
			metab.T2 = (87E-3)
		elif insysfile == 'pcr_N(CH3).sys':
			metab.A_m = 0.422
			metab.T2 = (90E-3)
		elif insysfile == 'pcr_X.sys':
			metab.A_m = 0.422
			metab.T2 = (81E-3)
		elif insysfile == 'peth.sys':
			metab.A_m = 0.126
			metab.T2 = (87E-3)
		elif insysfile == 'scyllo-inositol.sys':
			metab.A_m = 0.044
			metab.T2 = (87E-3)
		elif insysfile == 'taurine.sys':
			metab.A_m = 0.117
			metab.T2 = (85E-3)
		elif insysfile == 'water.sys':
			metab.A_m = 1.000
			metab.T2 = (43.60E-3)

		return metab

	def slaser_build180(self, inpulse180file, A_180, PULSE_180_LENGTH, gyratio, spin_system, plot_flag, scanner):
		# A: amplitude in milliTesla

		pulse180 = Pulse(inpulse180file, PULSE_180_LENGTH, scanner)

		n_old = np.linspace(0, PULSE_180_LENGTH, np.size(pulse180.waveform))
		n_new = np.linspace(0, PULSE_180_LENGTH, np.size(pulse180.waveform)+1)

		waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse180.waveform)*A_180)(n_new)
		waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse180.waveform)*A_180)(n_new)
		pulse180.waveform = waveform_real + 1j*(waveform_imag)

		if scanner == 'bruker':
			ampl_arr = np.abs(pulse180.waveform)
		else:
			ampl_arr = np.abs(pulse180.waveform) *gyratio

		phas_arr = np.unwrap(np.angle(pulse180.waveform))*180.0/math.pi
		freq_arr = np.gradient(phas_arr)

		if plot_flag:
			plt.figure()
			plt.subplot(3,1,1)
			plt.plot(n_new, waveform_real)
			plt.plot(n_new, waveform_imag)
			plt.subplot(3,1,2)
			plt.plot(n_new, ampl_arr)
			plt.subplot(3,1,3)
			plt.plot(n_new, freq_arr)
			plt.savefig(self.save_dir_sim + '180AFP' + '_' + str(PULSE_180_LENGTH) + 'sec.png')
			plt.close()

		pulse = pg.row_vector(len(pulse180.waveform))
		ptime = pg.row_vector(len(pulse180.waveform))
		for j, val in enumerate(zip(ampl_arr, phas_arr)):
			pulse.put(pg.complex(val[0],val[1]), j)
			ptime.put(pg.complex(n_new[1],0), j)

		pulse_dur_180 = pulse.size() * pulse180.pulsestep
		pwf_180 = pg.PulWaveform(pulse, ptime, "180afp")
		pulc_180 = pg.PulComposite(pwf_180, spin_system, self.sim_experiment.obs_iso)

		Ureal180 = pulc_180.GetUsum(-1)

		return A_180, pulse180, pulse_dur_180, Ureal180

	def laser_buildafp(self, inpulse180file, A_180, PULSE_180_LENGTH, gyratio, spin_system, plot_flag):
		# A: amplitude in milliTesla

		pulse180 = Pulse(inpulse180file, PULSE_180_LENGTH, 'varian')

		n_new = np.linspace(0, PULSE_180_LENGTH, 512)

		waveform_real = np.real(pulse180.waveform)*A_180
		waveform_imag = np.imag(pulse180.waveform)*A_180
		pulse180.waveform = waveform_real + 1j*(waveform_imag)

		ampl_arr = np.abs(pulse180.waveform)*gyratio
		phas_arr = np.unwrap(np.angle(pulse180.waveform))*180.0/math.pi
		freq_arr = np.gradient(phas_arr)

		if plot_flag:
			plt.figure()
			plt.subplot(3,1,1)
			plt.plot(n_new, waveform_real)
			plt.plot(n_new, waveform_imag)
			plt.subplot(3,1,2)
			plt.plot(n_new, ampl_arr)
			plt.subplot(3,1,3)
			plt.plot(n_new, freq_arr)
			plt.savefig(self.save_dir_sim + '180AFP' + '_' + str(PULSE_180_LENGTH) + 'sec.png')
			plt.close()

		pulse = pg.row_vector(len(pulse180.waveform))
		ptime = pg.row_vector(len(pulse180.waveform))
		for j, val in enumerate(zip(ampl_arr, phas_arr)):
			pulse.put(pg.complex(val[0],val[1]), j)
			ptime.put(pg.complex(n_new[1],0), j)

		pulse_dur_180 = pulse.size() * pulse180.pulsestep
		pwf_180 = pg.PulWaveform(pulse, ptime, "180afp")
		pulc_180 = pg.PulComposite(pwf_180, spin_system, self.sim_experiment.obs_iso)

		Ureal180 = pulc_180.GetUsum(-1)

		return A_180, pulse180, pulse_dur_180, Ureal180

	def slaser_build90(self, inpulse90file, A_90, PULSE_90_LENGTH, gyratio, spin_system, plot_flag, scanner):
		pulse90 = Pulse(inpulse90file, PULSE_90_LENGTH, scanner)

		n_old = np.linspace(0, PULSE_90_LENGTH, np.size(pulse90.waveform))
		n_new = np.linspace(0, PULSE_90_LENGTH, np.size(pulse90.waveform)+1)

		waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse90.waveform)*A_90)(n_new)
		waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse90.waveform)*A_90)(n_new)
		pulse90.waveform = waveform_real + 1j*(waveform_imag)

		if scanner == 'bruker':
			ampl_arr = np.abs(pulse90.waveform)
		else:
			ampl_arr = np.abs(pulse90.waveform)*gyratio

		phas_arr = np.unwrap(np.angle(pulse90.waveform))*180.0/math.pi
		
		if plot_flag:
			plt.figure()
			plt.subplot(3,1,1)
			plt.plot(n_new, waveform_real)
			plt.plot(n_new, waveform_imag)
			plt.subplot(3,1,2)
			plt.plot(n_new, ampl_arr)
			plt.subplot(3,1,3)
			plt.plot(n_new, phas_arr)
			plt.savefig(self.save_dir_sim + '90EXCITE' + '_' + str(PULSE_90_LENGTH) + 'sec.png')
			plt.close()

		pulse = pg.row_vector(len(pulse90.waveform))
		ptime = pg.row_vector(len(pulse90.waveform))
		for j, val in enumerate(zip(ampl_arr, phas_arr)):
			pulse.put(pg.complex(val[0],val[1]), j)
			ptime.put(pg.complex(pulse90.pulsestep,0), j)

		pulse_dur_90 = pulse.size() * pulse90.pulsestep
		if scanner == 'bruker':
			peak_to_end_90 = 0
		else:
			peak_to_end_90 = pulse_dur_90 - (209 + self.sim_experiment.fudge_factor) * pulse90.pulsestep
		pwf_90 = pg.PulWaveform(pulse, ptime, "90excite")
		pulc_90 = pg.PulComposite(pwf_90, spin_system, self.sim_experiment.obs_iso)

		Ureal90 = pulc_90.GetUsum(-1)

		return A_90, pulse90, pulse_dur_90, peak_to_end_90, Ureal90

	def laser_buildahp(self, inpulse90file, A_90, PULSE_90_LENGTH, gyratio, spin_system, plot_flag):
		pulse90 = Pulse(inpulse90file, PULSE_90_LENGTH, 'varian')

		n_old = np.linspace(0, PULSE_90_LENGTH, np.size(pulse90.waveform))
		n_new = np.linspace(0, PULSE_90_LENGTH, np.size(pulse90.waveform)+1)

		waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse90.waveform)*A_90)(n_new)
		waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse90.waveform)*A_90)(n_new)
		pulse90.waveform = waveform_real + 1j*(waveform_imag)

		ampl_arr = np.abs(pulse90.waveform)*gyratio
		phas_arr = np.unwrap(np.angle(pulse90.waveform))*180.0/math.pi
		
		if plot_flag:
			plt.figure()
			plt.subplot(3,1,1)
			plt.plot(n_new, waveform_real)
			plt.plot(n_new, waveform_imag)
			plt.subplot(3,1,2)
			plt.plot(n_new, ampl_arr)
			plt.subplot(3,1,3)
			plt.plot(n_new, phas_arr)
			plt.savefig(self.save_dir_sim + '90AHP' + '_' + str(PULSE_90_LENGTH) + 'sec.png')
			plt.close()

		pulse = pg.row_vector(len(pulse90.waveform))
		ptime = pg.row_vector(len(pulse90.waveform))
		for j, val in enumerate(zip(ampl_arr, phas_arr)):
			pulse.put(pg.complex(val[0],val[1]), j)
			ptime.put(pg.complex(pulse90.pulsestep,0), j)

		pulse_dur_90 = pulse.size() * pulse90.pulsestep
		pwf_90 = pg.PulWaveform(pulse, ptime, "90excite")
		pulc_90 = pg.PulComposite(pwf_90, spin_system, self.sim_experiment.obs_iso)

		Ureal90 = pulc_90.GetUsum(-1)

		return A_90, pulse90, pulse_dur_90, Ureal90

	# ---- Methods for Plotting in the GUI ---- #
	def setPlot(self, tab):
		if tab == 'Plot':
			fig = plt.figure(1)
			self.addmpl(0, fig, self.plotTime_mplvl)

			fig = plt.figure(2)
			self.addmpl(1, fig, self.plotFreq_mplvl)

		elif tab == 'Guess':
			fig = plt.figure(3)
			self.addmpl(2, fig, self.plotGES_mplvl)

	def addmpl(self, canvas_index, fig, vertical_layout):
		self.canvas[canvas_index] = FigureCanvas(fig)
		vertical_layout.addWidget(self.canvas[canvas_index])
		self.canvas[canvas_index].draw()

		self.toolbar[canvas_index] = NavigationToolbar(self.canvas[canvas_index], self, coordinates=True)
		vertical_layout.addWidget(self.toolbar[canvas_index])

	# ---- Methods to Load and Specify Experiments ---- #
	def chooseExperimentFile(self):
		self.filenameConfirmButton.setEnabled(False)
		prev = str(self.filenameInput.text())
		self.filenameInput.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Experiment File', os.getcwd() + '/pints/experiments', 'Text files (*.txt)')[0]))
		self.filenameConfirmButton.setEnabled(True)
		if str(self.filenameInput.text()) == '':
			self.filenameInput.setText(str(prev))
			if str(prev) == '':
				self.filenameConfirmButton.setEnabled(False)

	def loadExperimentFile(self):
		try:
			self.specExpButton.setEnabled(False)

			self.editMetabOKButton.setEnabled(False)

			self.plotButton.setEnabled(False)
			
			self.exportButton_rawTime.setEnabled(False)
			self.exportButton_rawFreq.setEnabled(False)
			self.exportButton_dat.setEnabled(False)

			self.cstOKButton.setEnabled(False)
			self.cstLoadButton.setEnabled(False)

			self.genCSTButton.setEnabled(False)
			self.cstFilenameConfirmButton.setEnabled(False)

			self.filenameBrowseButton_dat.setEnabled(False)
			self.filenameConfirmButton_dat.setEnabled(False)
			self.plotGESButton.setEnabled(False)

			self.gesFilenameConfirmButton.setEnabled(False)
			self.genGESButton.setEnabled(False)
			
			filename = str(self.filenameInput.text())
			in_file = open(filename, 'r')
			self.experiment = self.loadExperiment(filename, in_file)
			
			self.metabListWidget.clear()
			self.metabListWidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
			self.metabListWidget.addItems(self.experiment.metabolites)

			self.populateMetabTable()
			self.populateMetabList()

			self.loop1Input.setText(str(self.experiment.TE))
			self.loop1Input.setEnabled(False)

			self.outputTEInput.setText(str(self.experiment.TE))
			self.outputTEInput.setEnabled(False)
			
			self.filenameInfoLabel.setText("Experiment loaded! \n >> " + self.experiment.name)

			self.outputFilenameInfoLabel_dat.append(">> Output .dat Filename: " + str(self.filenameInput.text()).replace(".txt", ".dat"))
			self.outputFilenameInfoLabel_rawTime.append(">> Output .dat Filename: " + str(self.filenameInput.text()).rsplit('/',1)[0] + "/bf_raw_t/" + str(self.filenameInput.text()).rsplit('/',1)[-1].replace(".txt", "") + "_METAB_t.raw")
			self.outputFilenameInfoLabel_rawFreq.append(">> Output .dat Filename: " + str(self.filenameInput.text()).rsplit('/',1)[0] + "/bf_raw_f/" + str(self.filenameInput.text()).rsplit('/',1)[-1].replace(".txt", "") + "_METAB_f.raw")

			self.cstFilenameInfoLabel.append(">> Output .cst Filename:" + str(self.filenameInput.text()).replace(".txt", ".cst"))
			self.gesFilenameInfoLabel.append(">> Output .ges Filename:" + str(self.filenameInput.text()).replace(".txt", ".ges"))
			
			self.outputFilenameInput_dat.setText(str(self.filenameInput.text()).replace(".txt", ".dat"))
			self.outputFilenameInput_rawTime.setText(str(self.filenameInput.text()).rsplit('/',1)[0] + "/bf_raw_t/" + str(self.filenameInput.text()).rsplit('/',1)[-1].replace(".txt", "") + "_METAB_t.raw")
			self.outputFilenameInput_rawFreq.setText(str(self.filenameInput.text()).rsplit('/',1)[0] + "/bf_raw_f/" + str(self.filenameInput.text()).rsplit('/',1)[-1].replace(".txt", "") + "_METAB_f.raw")

			self.cstFilenameInput.setText(str(self.filenameInput.text()).replace(".txt", ".cst"))
			self.gesFilenameInput.setText(str(self.filenameInput.text()).replace(".txt", ".ges"))

			self.specExpInfoLabel.setText("Experiment not fully specified!")    

			self.specExpButton.setEnabled(True)  
			
			in_file.close()
		except Exception as e:
			self.filenameInfoLabel.setText("ERROR: " + str(e) + "\n>> Error reading file.\n>> Please select another file.")

	def loadExperiment(self, filename, in_file):
		# Create a new blank experiment.
		experiment = Experiment()

		# Read data from file.
		metab_info = []
		print('Reading data from: ' + filename + ' ...')
		for line in in_file:
			if ";Name:" in line:
				experiment.name = line.replace('\n', '').replace(': ', ':').split(':')[-1]
				print('  | Name: ' + experiment.name)
			elif ";Created:" in line:
				experiment.date = line.replace('\n', '').split(' ')[-1]
				print('  | Date: ' + experiment.date)
			elif ";Comment:" in line:
				experiment.description = line.replace('\n', '').replace(': ', ':').split(':')[-1]
				print('  | Description: ' + experiment.description)    
			elif ";PI:" in line:
				experiment.author = line.replace('\n', '').replace(': ', ':').split(':')[-1]
				print('  | Author: ' + experiment.author)
			elif ";b0:" in line:
				experiment.b0 = float(line.replace('\n', '').split(' ')[-1])
				print('  | b0 = ' + str(experiment.b0))
			elif " Metabolites:" in line:
				experiment.metabolites_num = int(float(line.replace('\n', '').replace(';', '').split(' ')[0]))
				print('  | ' + str(experiment.metabolites_num) + ' metabolites simulated in this experiment.')
			elif ";PINTS for FITMAN" in line:
				experiment.type = 'PINTS'
				print(experiment.type)
			elif ";Vespa-Simulation Mixed Metabolite Output" in line:
				experiment.type = 'VESPA'
				print(experiment.type)
			elif (";" in line) or ("http:" in line):
				pass
			else:
				metab_info.append(line.replace('\n', '').split('\t'))
		print('  | Reading complete.\n')
		print('Importing data ...')

		# load information from each data line
		for currline in metab_info:
			if experiment.type == 'VESPA':
				print('  | Importing ' + currline[0] + ' at (' + currline[1] + ', ' + currline[2] + ', ' + currline[3] + ', ' + currline[4] + ')')
				name = currline[0].replace('vgov_','').replace('kbehar_','').replace('seeger_','').split('_')[0]
				loop1 = float(currline[1])
				# loop2 = float(currline[2])
				# loop3 = float(currline[3])

				experiment.TE = loop1
				if experiment.TE >= 135:
					experiment.splitNAA = True

				# index into data structure and check object type
				if type(experiment.data[name]) is defaultdict:
					# no metabolite object at this node ... so create it
					experiment.data[name] = Metabolite()
					experiment.data[name].name = name
					# if metabolite name is not already recorded, append to list of available metabolties
					if not(name in experiment.metabolites):
						experiment.metabolites.append(name)
						print('    | ' + name + ' added to list of available metabolites.')

				# append data
				experiment.data[name].area.append(float(currline[6]))
				experiment.data[name].ppm.append(float(currline[5]))
				experiment.data[name].phase.append(float(currline[7]))

			elif experiment.type =='PINTS':
				name = currline[0]
				experiment.TE = float(currline[1])

				if name in ['lm1', 'lm2', 'lm3', 'lm4', 'lm5', 'lm6', 'lm7', 'lm8', 'lm9', 'lm10', 'lm11', 'lm12', 'lm13', 'lm14', 'MM09', 'MM12', 'MM14', 'MM16', 'MM20', 'MM21', 'MM23', 'MM26', 'MM30', 'MM31', 'MM37', 'MM38', 'MM40']:
					# name, TE, A_m, shift, line_type, lw, area, phase
					experiment.data[name] = Macromolecule(name, float(currline[3]), str(currline[4]), float(currline[5]), float(currline[6]), float(currline[7]))
					# if metabolite name is not already recorded, append to list of available metabolties
					if not(name in experiment.metabolites):
						experiment.metabolites.append(name)
						print('    | ' + name + ' added to list of available metabolites.')
				else:
					# index into data structure and check object type
					if type(experiment.data[name]) is defaultdict:
						# no metabolite object at this node ... so create it
						experiment.data[name] = Metabolite()
						experiment.data[name].name = name
						experiment.data[name].A_m = float(currline[2])
						experiment.data[name].T2  = float(currline[3])
						# if metabolite name is not already recorded, append to list of available metabolties
						if not(name in experiment.metabolites):
							experiment.metabolites.append(name)
							print('    | ' + name + ' added to list of available metabolites.')

					# append data
					experiment.data[name].area.append(float(currline[6]))
					experiment.data[name].ppm.append(float(currline[5]))
					experiment.data[name].phase.append(float(currline[7]))

		print('  | Total metabolites in file: ' + str(experiment.metabolites_num))
		print('  | Total metabolites after combination: ' + str(np.size(experiment.metabolites)))
		print('       | ', experiment.metabolites.sort())
		print('  | Import complete.\n')

		return experiment

	def specExp(self):		
		try:
			self.experiment.dwell_time = float(self.dwellTimeInput.text())
			self.experiment.acq_length = float(self.acqLengthInput.text())
			self.experiment.t = np.arange(0,self.experiment.acq_length,self.experiment.dwell_time)
			self.experiment.fs = 1/self.experiment.dwell_time
			self.specExpInfoLabel.setText("Experiment fully specified! \n >> Dwell-time = " + str(self.experiment.dwell_time) + " sec \n >> Acq. Length = " + str(self.experiment.acq_length) + " sec")

			self.loadMetabBrowseButton.setEnabled(True)

			self.editMetabOKButton.setEnabled(True)
			self.plotButton.setEnabled(True)

			self.exportButton_rawTime.setEnabled(True)
			self.outputFilenameConfirmButton_rawTime.setEnabled(True)

			self.exportButton_rawFreq.setEnabled(True)
			self.outputFilenameConfirmButton_rawFreq.setEnabled(True)

			self.exportButton_dat.setEnabled(True)
			self.outputFilenameConfirmButton_dat.setEnabled(True)

			self.cstOKButton.setEnabled(False)
			self.cstLoadButton.setEnabled(True)

			self.linkAllCheckBox.setChecked(self.experiment.LW_linkall)
			self.splitNAACheckBox.setChecked(self.experiment.splitNAA)

			self.genCSTButton.setEnabled(False)
			self.cstFilenameConfirmButton.setEnabled(True)

			self.filenameBrowseButton_dat.setEnabled(True)
			self.filenameConfirmButton_dat.setEnabled(False)
		except Exception as e:
			self.specExpInfoLabel.setText("ERROR: " + str(e) + "\n>> Experiment not specified!\n>> Please try again.")

	# ---- Methods for 'Edit Metabolites' Tab ---- #
	def populateMetabTable(self):
		metabolites = self.experiment.metabolites
		
		self.metabTableWidget.setRowCount(np.size(metabolites))
		self.metabTableWidget.setColumnCount(3)

		for i in range(0, np.size(metabolites)):
			metab = self.experiment.data[metabolites[i]]
			self.metabTableWidget.setItem(i, 0, QtWidgets.QTableWidgetItem(str(metab.name)))
			self.metabTableWidget.setItem(i, 1, QtWidgets.QTableWidgetItem(str(metab.A_m)))
			self.metabTableWidget.setItem(i, 2, QtWidgets.QTableWidgetItem(str(metab.T2)))

	def saveMetabInfo(self):
		metabolites = self.experiment.metabolites
		for i in range(0, np.size(metabolites)):
			metab = self.experiment.data[metabolites[i]]
			metab.A_m = float(self.metabTableWidget.item(i, 1).text())
			metab.T2  = float(self.metabTableWidget.item(i, 2).text())
			print(str(metabolites[i]) + ',' + str(metab.A_m) + ',' + str(metab.T2))
		self.loadMetabInfoLabel.setText("Current metabolite information applied to experiment.")

	def loadMetabInfoFile(self):
		try:
			in_file = open(str(self.loadMetabFileInput.text()), 'r')

			A_m_params = []
			T2_params = []
			for line in in_file:
				params = line.replace('\n', '').split('\t')
				if np.size(params) == 3:
					print(params)
					A_m_params.append(float(params[1]))
					T2_params.append(float(params[2]))

			metabolites = self.experiment.metabolites
			if np.size(metabolites) == np.size(A_m_params) and np.size(metabolites) == np.size(T2_params):
				for i in range(0, np.size(metabolites)):
					self.metabTableWidget.setItem(i, 1, QtWidgets.QTableWidgetItem(str(A_m_params[i])))
					self.metabTableWidget.setItem(i, 2, QtWidgets.QTableWidgetItem(str(T2_params[i])))
				self.loadMetabInfoLabel.setText("Parameter file loaded successfully!")
			else:
				self.loadMetabInfoLabel.setText("Error in parameter file! >> num of parameters " + str((np.size(A_m_params), np.size(T2_params))) + " != num of metabolites " + " " + str(np.size(metabolites)))
		except Exception as e:
			self.loadMetabInfoLabel.setText("ERROR: " + str(e) + " >> Please try another parameter file.")

	def chooseMetabInfoFile(self):
		self.loadMetabConfirmButton.setEnabled(False)
		prev = str(self.loadMetabFileInput.text())
		self.loadMetabFileInput.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Metabolite File', os.getcwd() + '/pints/metabinfo', 'Metabolite files (*.metab)')[0]))
		self.loadMetabConfirmButton.setEnabled(True)
		if str(self.loadMetabFileInput.text()) == '':
			self.loadMetabFileInput.setText(str(prev))
			if str(prev) == '':
				self.loadMetabConfirmButton.setEnabled(False)

	# ---- Methods to Calculate and Plot Spectra ---- #
	def plot(self):
		try:
			metabolites = []
			print('----- Plot a Sum of Metabolites -----')            
			print('Metabolites chosen:')

			for item in self.metabListWidget.selectedItems():
				metabolites.append(str(item.text()))
			print(metabolites)

			w0 =   0
			phi0 = np.deg2rad(float(self.phaseInput.text()))
			lb =   float(self.lbInput.text())        
			shift = float(self.shiftInput.text())

			total_signal = self.calcSignal(metabolites, self.experiment, w0, phi0, lb, shift) 
			
			self.plotFIDs(self.experiment, total_signal)
			self.plotSpectra(self.experiment, total_signal)

			self.plotErrorLabel.setText("Selected metabolites successfully plotted.")
		except Exception as e:
			self.plotErrorLabel.setText("ERROR: " + str(e) + " >> Please try again.")

	def calcSignal(self, metabolites, experiment, w0, phi0, lb, shift):
		
		signals = np.empty([np.size(metabolites), np.size(experiment.t)], dtype=complex)
		for i in range(0, np.size(metabolites)):
			signals[i,] = np.exp(-1j * (w0 * experiment.t) + 1j * (phi0)) \
						  * experiment.data[metabolites[i]].getFID(experiment.TE, experiment.b0, experiment.t, 0, 1, 0, 0, lb) \
						  * np.exp(-1j * 2 * sp.pi * experiment.b0 * shift * experiment.t)
		total_signal = np.sum(signals, axis=0)
		
		print('')
		
		return total_signal

	def plotFIDs(self, experiment, total_signal):
		plt.figure(1)
		plt.clf()
		plt.subplot(2,1,1)
		plt.plot(experiment.t, np.real(total_signal))
		plt.xlim(0,experiment.acq_length)
		plt.xlabel('t [s]')
		plt.title('FID (Real)')
		plt.subplot(2,1,2)
		plt.plot(experiment.t, np.imag(total_signal))
		plt.xlim(0,experiment.acq_length)
		plt.xlabel('t [s]')
		plt.title('FID (Imaginary)')
		self.canvas[0].draw()

	def plotSpectra(self, experiment, total_signal):
		n = np.size(total_signal)
		print(n)
		f = np.arange(+n//2,-n//2,-1)*(self.experiment.fs/n)*(1/self.experiment.b0)
		spectra = fftw.fftshift(fftw.fft(total_signal))
		plt.figure(2)
		plt.clf()
		plt.subplot(2,1,1)
		plt.plot(-f, np.real(spectra))
		plt.title('Spectra (Real)')
		plt.xlim(5,0)
		plt.xlabel('ppm')    
		plt.subplot(2,1,2)
		plt.plot(-f, np.imag(spectra))
		plt.title('Spectra (Imaginary)')    
		plt.xlim(5,0)
		plt.xlabel('ppm')
		self.canvas[1].draw()

	# ---- Methods to Export Data ---- #
	def setOutputFile_dat(self):
		self.outputFilenameInfoLabel_dat.append(">> Output file confirmed:\n" + str(self.outputFilenameInput_dat.text()))

	def setOutputFile_rawTime(self):
		self.outputFilenameInfoLabel_rawTime.append(">> Output file confirmed:\n" + str(self.outputFilenameInput_rawTime.text()))
		self.outputFilenameInfoLabel_rawTime.append(">> Note: 'METAB' in the filename will be replaced with the full metabolite name.")

	def setOutputFile_rawFreq(self):
		self.outputFilenameInfoLabel_rawFreq.append(">> Output file confirmed:\n" + str(self.outputFilenameInput_rawFreq.text()))
		self.outputFilenameInfoLabel_rawFreq.append(">> Note: 'METAB' in the filename will be replaced with the full metabolite name.")

	def export_dat(self):
		try:
			print('----- Export current experiment to .dat file -----')        
			metabolites = []
			for item in self.metabListWidget.selectedItems():
				metabolites.append(str(item.text()))
			print(metabolites)

			w0 =   0
			phi0 = np.deg2rad(float(self.phaseInput.text()))
			lb =   float(self.lbInput.text())
			shift = float(self.shiftInput.text())

			total_signal = self.calcSignal(metabolites, self.experiment, w0, phi0, lb, shift)
			
			self.saveToDAT(self.experiment, total_signal)
		except Exception as e:
			self.outputFilenameInfoLabel_dat.append("\n>> ERROR " + str(e) + "\nPLEASE TRY AGAIN.\n")

	def export_rawTime(self):
		try:
			print('----- Export current experiment to .raw files -----')
			print('MODE: Time domain data')

			metabolites = []
			for item in self.metabListWidget.selectedItems():
				metabolites.append(str(item.text()))
			print(metabolites)

			w0 =   0
			phi0 = np.deg2rad(float(self.phaseInput.text()))
			lb =   float(self.lbInput.text()) 
			shift = float(self.shiftInput.text())

			directory = str(self.outputFilenameInput_rawTime.text()).rsplit('/', 1)[0]
			print(directory)

			if not os.path.exists(directory):
				os.makedirs(directory)

			for metabolite in metabolites:
				
				filename_out = str(self.outputFilenameInput_rawTime.text().replace('METAB', str(metabolite)))
				raw_file = open(filename_out, 'w')

				raw_file.write('Robarts Research Institute' + '\n')

				# output experiment information
				raw_file.write('Simulation Name = ' + str(self.experiment.name) + '\n')
				raw_file.write('Simulation Date = ' + str(self.experiment.name) + '\n')
				raw_file.write('Simulation Description = ' + str(self.experiment.description) + '\n')
				raw_file.write('B0 = ' + str(self.experiment.b0) + 'MHz \n')
				raw_file.write('TE = ' + str(self.experiment.TE) + 'ms \n')
				raw_file.write('Dwell Time = ' + str(self.experiment.dwell_time) + 'sec \n')
				raw_file.write('Acquisition Length = ' + str(self.experiment.acq_length) + 'sec \n')
				raw_file.write('Sampling Frequnecy = ' + str(self.experiment.fs) + 'Hz \n')
				raw_file.write('Lorentzian Line Broadening = ' + str(self.lbInput.text()) + 'Hz \n')
				raw_file.write('\n')
				raw_file.write('ID = ' + metabolite + '\n')

				# begin data output
				raw_file.write('REAL \t IMAGINARY + \n')
				signal = self.calcSignal([metabolite], self.experiment, w0, phi0, lb, shift)
				for element in signal:
					raw_file.write("{0:.6E}".format(float(np.real(element))) + '\t')
					raw_file.write("{0:.6E}".format(float(np.imag(element))) + '\t')
					raw_file.write('\n')

				self.outputFilenameInfoLabel_rawTime.append(">> Successfully output " + metabolite + " to " + filename_out)

				raw_file.close()
		except Exception as e:
			self.outputFilenameInput_rawTime.append("\n>> ERROR " + str(e) + "\nPLEASE TRY AGAIN.\n")

	def export_rawFreq(self):
		try:
			print('----- Export current experiment to .raw files -----')
			print('MODE: Frequency domain data')

			metabolites = []
			for item in self.metabListWidget.selectedItems():
				metabolites.append(str(item.text()))
			print(metabolites)

			w0 =   0
			phi0 = np.deg2rad(float(self.phaseInput.text()))
			lb =   float(self.lbInput.text()) 
			shift = float(self.shiftInput.text())

			directory = str(self.outputFilenameInput_rawFreq.text()).rsplit('/', 1)[0]
			print(directory)

			if not os.path.exists(directory):
				os.makedirs(directory)

			for metabolite in metabolites:
				
				filename_out = str(self.outputFilenameInput_rawFreq.text().replace('METAB', str(metabolite)))
				raw_file = open(filename_out, 'w')

				raw_file.write('Robarts Research Institute' + '\n')

				# output experiment information
				raw_file.write('Simulation Name = ' + str(self.experiment.name) + '\n')
				raw_file.write('Simulation Date = ' + str(self.experiment.name) + '\n')
				raw_file.write('Simulation Description = ' + str(self.experiment.description) + '\n')
				raw_file.write('B0 = ' + str(self.experiment.b0) + 'MHz \n')
				raw_file.write('TE = ' + str(self.experiment.TE) + 'ms \n')
				raw_file.write('Dwell Time = ' + str(self.experiment.dwell_time) + 'sec \n')
				raw_file.write('Acquisition Length = ' + str(self.experiment.acq_length) + 'sec \n')
				raw_file.write('Sampling Frequnecy = ' + str(self.experiment.fs) + 'Hz \n')
				raw_file.write('Lorentzian Line Broadening = ' + str(self.lbInput.text()) + 'Hz \n')
				raw_file.write('\n')
				raw_file.write('ID = ' + metabolite + '\n')

				# begin data output
				raw_file.write('REAL \t IMAGINARY + \n')
				total_signal = self.calcSignal([metabolite], self.experiment, w0, phi0, lb, shift)
				n = np.size(total_signal)
				f = np.arange(+n//2,-n//2,-1)*(self.experiment.fs/n)*(1/self.experiment.b0)
				spectra = fftw.fftshift(fftw.fft(total_signal))
				for i in range(0, n):
					element = spectra[i]
					shift = f[i]

					if f[i] <= 5 and f[i] >= 0:
						raw_file.write("{0:.6E}".format(float(np.real(element))) + '\t')
						raw_file.write("{0:.6E}".format(float(np.imag(element))) + '\t')
						raw_file.write('\n')

				self.outputFilenameInfoLabel_rawFreq.append(">> Successfully output " + metabolite + " to " + filename_out)

				raw_file.close()
		except Exception as e:
			self.outputFilenameInfoLabel_rawFreq.append("\n>> ERROR " + str(e) + "\nPLEASE TRY AGAIN.\n")

	def saveToDAT(self, experiment, signal):
		filename_out = str(self.outputFilenameInput_dat.text())
		dat_file = open(filename_out, 'w')
		
		dat_file.write(str(np.size(signal)*2) + '\n') # num of data lines
		dat_file.write(str(1) + '\n') # num of averages
		dat_file.write(str(experiment.dwell_time) + '\n') # ADC dwell time
		dat_file.write(str(experiment.b0) + '\n') # main field strength
		
		dat_file.write(str(1) + '\n') # (?) NOT SURE WHAT THIS LINE IS (?)
		
		dat_file.write("Robarts Research Institute, " + experiment.name + '\n')
		dat_file.write(' \"' + experiment.date + '\"' + '\n')
		
		dat_file.write("MachS=0 ConvS=1.00e-03 V1=20.000 V2=20.000 V3=20.000 vtheta=0.0" + '\n')
		# ^arbitrary line, since not a real experiment off the scanner
	
		dat_file.write("TE=" + str(self.outputTEInput.text()) + " s TR=7.500 s P1=9.02056 P2=37.26239 P3=6.46689 Gain=1.00" + '\n')
		
		# (?) NOT SURE WHAT THESE LINES ARE FOR (?)
		dat_file.write("SIMULTANEOUS" + '\n')
		dat_file.write('0.0' + '\n')
		dat_file.write('EMPTY' + '\n')
		
		# begin to write data
		for element in signal:
			dat_file.write("{0:.6f}".format(float(np.real(element))) + '\n')
			dat_file.write("{0:.6f}".format(float(np.imag(element))) + '\n')
		
		self.outputFilenameInfoLabel_dat.append(">> Successfully output to " + filename_out)

		dat_file.close()

	# ---- Methods for 'Edit Constraint Metabolites' Tab ---- #
	def populateMetabList(self):
		metabolites = self.experiment.metabolites
		for metabolite in metabolites:
			self.metabTextBrowser.append(str(metabolite))

	def loadCSTInfo(self):
		try:
			prev_filename = str(self.cstFilenameLineEdit.text())
			cst_filename = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Constraint Information File', os.getcwd() + '/pints/cstinfo', 'Constraint information files (*.cstinfo)')[0])
			if cst_filename == '':
				cst_filename = prev_filename

			if cst_filename != '':
				self.cstFilenameInfoLabel.append('Loading constraints from: ' + cst_filename)
				cst_file = open(cst_filename, 'r')

				self.shiftsTextEdit.clear()
				self.ampsTextEdit.clear()
				self.phasesTextEdit.clear()

				for line in cst_file:
					if '[Shifts]' in line:
						dest = self.shiftsTextEdit
					elif '[Amplitudes]' in line:
						dest = self.ampsTextEdit
					elif '[Phases]' in line:
						dest = self.phasesTextEdit
					else:
						dest.appendPlainText(line.replace('\n',''))

				self.cst_filename = cst_filename
				self.cstFilenameLineEdit.setText(str(cst_filename))

				cst_file.close()
				self.cstInfoErrorLabel.setText('File loaded!')
				self.cstOKButton.setEnabled(True)
			else:
				self.cstInfoErrorLabel.setText('Please specify a .cstinfo file!')
		except Exception as e:
			self.cstInfoErrorLabel.setText('ERROR: ' + str(e))

	def saveCSTInfo(self):
		try:
			self.experiment.cst_groups = []
			
			prev_filename = str(self.cstFilenameLineEdit.text())		
			cst_filename = str(QtWidgets.QFileDialog.getSaveFileName(self, 'Save Constraint Information File', os.getcwd() + '/pints/cstinfo', 'Constraint information files (*.cstinfo)')[0])
			if cst_filename == '':
				cst_filename = prev_filename

			if cst_filename != '':
				cst_file = open(cst_filename, 'w')

				self.cstFilenameInfoLabel.append('Reading shift groups ...')
				self.readGroups('shift')
				cst_file.write('[Shifts]\n')
				cst_file.write(self.shiftsTextEdit.toPlainText())

				self.cstFilenameInfoLabel.append('Reading amp groups ...')
				self.readGroups('amp')
				cst_file.write('\n[Amplitudes]\n')
				cst_file.write(self.ampsTextEdit.toPlainText())

				self.cstFilenameInfoLabel.append('Reading phase groups ...')
				self.readGroups('phase')
				cst_file.write('\n[Phases]\n')
				cst_file.write(self.phasesTextEdit.toPlainText())

				if self.linkLipidsCheckBox.isChecked():
					self.cstFilenameInfoLabel.append('Linking lipid linewidths ...')
					self.readGroups('lipid_LWs')
					# self.readGroups('lipid_WLs') # (new macromolecule model has gaussian lineshapes)

				self.experiment.LW_limit = float(self.maxLWInput.text())
				self.cstFilenameInfoLabel.append('LW_limit = ' + str(self.experiment.LW_limit))

				self.experiment.LW_linkall = self.linkAllCheckBox.isChecked()
				self.cstFilenameInfoLabel.append('LW_linkall = ' + str(self.experiment.LW_linkall))

				self.experiment.splitNAA = self.splitNAACheckBox.isChecked()
				self.cstFilenameInfoLabel.append('splitNAA  = ' + str(self.experiment.splitNAA))

				self.cst_filename = cst_filename

				self.cstFilenameInfoLabel.append(str(len(self.experiment.cst_groups)) + ' group(s) saved!')

				self.cst_filename = cst_filename
				self.cstFilenameLineEdit.setText(str(cst_filename))
				self.cstFilenameInfoLabel.append('Writing to: ' + str(self.cst_filename))

				cst_file.close()
				self.cstInfoErrorLabel.setText('File saved!')
				self.genCSTButton.setEnabled(True)
			else:
				self.cstInfoErrorLabel.setText('Please specify a .cstinfo file!')
		except Exception as e:
			self.cstInfoErrorLabel.setText('ERROR: ' + str(e) + '\n' + str(sys.exc_info()[2].tb_lineno))

	def readGroups(self, group_type):
		if group_type == 'shift':
			source = self.shiftsTextEdit.toPlainText()
		elif group_type == 'amp':
			source = self.ampsTextEdit.toPlainText()
		elif group_type == 'phase':
			source = self.phasesTextEdit.toPlainText()
		elif group_type == 'lipid_LWs':
			if self.experiment.b0 == 123.3:
				source = 'LWlipMM, -inf, +inf, MM09, MM12, MM14, MM16, MM20, MM21, MM23, MM26, MM30, MM31, MM37, MM38, MM40;'
			else:
				source = 'LWlmX, -inf, +inf, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11, lm12, lm13, lm14;'

		# parse groups in a column
		defined_groups = str(source).replace('\n', '').replace(' ', '').split(';')
		defined_groups = defined_groups[0:-1]

		if group_type == 'amp':
			# Reading amplitude groups (this must always come after reading shift groups!!!)
			
			# Get metabolites
			metabolites = []
			for group in self.experiment.cst_groups:
				if group.typeCST == 'shift':
					for member in group.members:
						if not(member in metabolites):
							metabolites.append(member)
			self.cstFilenameInfoLabel.append('Grabbing metabolites ...' + str(metabolites))

			amp_groups = []

			for i in range(0, np.size(metabolites)):
				metab = self.experiment.data[metabolites[i]]

				# Get submetabolites
				sub_metabs = [m for m in metabolites if str(metab.name_short()) == m.split('_')[0]]

				# If sub-metabolites have different T2s (e.g naa_acetyl and naa_aspartate),
				# split into different constraint groups so their amplitudes can be fit separately
				sub_metabs_grouped = {}
				for sub_metab in sub_metabs:
					key = self.experiment.data[sub_metab].T2
					sub_metabs_grouped.setdefault(key,[]).append(sub_metab)
				sub_metabs_groups = list(sub_metabs_grouped.values())

				# Define sub_metabolite group
				for (group_num, sub_metab_group) in enumerate(sub_metabs_groups):
					new_defined_group = str(metab.name_short()) + '_' + str(group_num) + ',-inf,+inf,'
					for group_member in sub_metab_group:
						USER_DEFINED_FLAG = False
						for defined_group in defined_groups:
							# If the group member is already defined by an amplitude group created by the user,
							# don't add to the automatically created group
							
							# print '    |', group_member, defined_group, group_member in defined_group.split(',')
							if group_member in defined_group.split(','):
								# print '        |', group_member, ' is already part of a user-defined group.'
								USER_DEFINED_FLAG = True

						if not(USER_DEFINED_FLAG):
							new_defined_group = new_defined_group + str(group_member) + ','
							# print '        | Adding', group_member, 'to auto-defined group.'
					
					self.cstFilenameInfoLabel.append('  | Auto-defined group: ' + str(new_defined_group) + ' ... ' + str(np.size(new_defined_group[0:-1].split(','))))
					if not(new_defined_group[0:-1] in amp_groups) and np.size(new_defined_group[0:-1].split(',')) > 3:
						amp_groups.append(new_defined_group[0:-1])

			# Add the amplitude groups to the list of defined_groups
			for group in amp_groups:
				defined_groups.append(group)		

		self.cstFilenameInfoLabel.append('  | Defined groups:')
		for group in defined_groups:
			# parse each group's info
			self.cstFilenameInfoLabel.append('     | ' + str(group))
			group_info = group.split(',')
			group_name = group_info[0]
			group_members = group_info[3:]
			group_min = -sp.inf if group_info[1] == '-inf' else group_info[1]
			group_max = +sp.inf if group_info[2] == '+inf' else group_info[2]

			# FITMAN SEEMS TO SET THE FIRST PEAK IT SEES AS THE REFERENCE VALUE FOR LINKED METABOLITES IN A GROUP
			# CODE BELOW MATCHES THIS
			ref_values = [0,0,0]
			member = group_members[0]
			metab = self.experiment.data[member]
			ref_values[0] = metab.ppm[0]
			ref_values[1] = metab.area[0] * metab.A_m
			ref_values[2] = metab.phase[0]

			new_cstgroup = CSTGroup(group_type, group_name, group_members, float(group_min), float(group_max))
			new_cstgroup.ref_T2 = metab.T2
			if group_type == 'shift':
				new_cstgroup.ref_value = ref_values[0]
			elif group_type == 'amp':
				new_cstgroup.ref_value = ref_values[1]
			elif group_type == 'phase':
				new_cstgroup.ref_value = ref_values[2]
			elif group_type == 'lipid_LWs' or group_type == 'lipid_WLs':
				new_cstgroup.ref_value = metab.lw
			
			self.experiment.cst_groups.append(new_cstgroup)

	def setCSTFile(self):
		self.cstFilenameInfoLabel.append(">> Output .cst File Confirmed:\n" + str(self.cstFilenameInput.text()))

	def genCST(self):
		try:
			print('----- Generating .cst file from current experiment -----')
			
			experiment = self.experiment

			self.savetoCST(self.experiment)
		except Exception as e:
			self.cstFilenameInfoLabel.append('\n>> ERROR ' + str (e) + '\nCould not generate constraint file.\n')

	def savetoCST(self, experiment):
		filename_out = str(self.cstFilenameInput.text())
		cst_file = open(filename_out, 'w')

		pk_count = 0
		output = []

		output.append("****_Constraints_File_Begins_****  <--- Do not remove this line" + '\n')
		output.append("[Parameters]" + '\n')
		output.append("number_peaks\t\t" + '<PLACEHOLDER>' + '\n')
		output.append("tolerence\t\t0.001" + '\n')
		output.append("minimum_iterations\t\t5" + '\n')
		output.append("maximum_iterations\t\t500" + '\n')
		output.append("range\t\t0\t" + str(int(experiment.acq_length / experiment.dwell_time * 0.3) * 2) + '\n')
		output.append("fwhm_exp_weighting\t\t0" + '\n')
		output.append("shift_units\t\tppm" + '\n')
		output.append("output_shift_unit\t\tppm" + '\n')
		output.append("alambda_inc\t\t1000" + '\n')
		output.append("alambda_dec\t\t0.5" + '\n')
		output.append("noise_equal" + '\n')

		output.append('\n')

		output.append("[Variables]" + '\n')

		output.append('\n')

		output.append("[Peaks]" + '\n')
		output.append(";PEAK#\tSHIFT\tLORENTZIAN WIDTH\tAMPLITUDE\tPHASE\tDELAY\tGAUSSIAN WIDTH" + '\n')

		# ---- Get from selected items from 'Plot Metabolites' tab ----
		# metabolites = []
		# for item in self.metabListWidget.selectedItems():
		# 	metabolites.append(str(item.text()))
		# print metabolites

		# ---- Get from shift constraint groups ---- #
		metabolites = []
		for group in experiment.cst_groups:
			if group.typeCST == 'shift':
				for member in group.members:
					if not(member in metabolites):
						metabolites.append(member)
		print(metabolites)

		for i in range(0, np.size(metabolites)):
			
			# get metabolite
			metab = experiment.data[metabolites[i]]
			print('')
			print(metab.name + ',', end=' ')

			output.append(";" + str(metab.name) + '\n')

			# set default reference values
			ref_shift = float(metab.ppm[0])
			ref_amp = float(metab.area[0] * metab.A_m)
			ref_T2 = float(metab.T2)
			ref_phase = float(sp.deg2rad(metab.phase[0]))

			# set default max/min values
			min_shift = -sp.inf
			max_shift = sp.inf

			min_amp = -sp.inf
			max_amp = sp.inf

			min_phase = -sp.inf
			max_phase = sp.inf

			# set default constraint group names
			group_name_shift = "S" + str(metab.name_short())
			group_name_amp = str(metab.name_short())

			group_name_phase = "P" + str(metab.name_short())
			group_name_lw = "LW" + str(metab.name_short())
			group_name_wl = "WL" + str(metab.name_short())

			# check if metabolite is part of a constraint group
			# and override defaults with group parameters
			for group in experiment.cst_groups:
				if metab.name in group.members:
					print(group.name + ',', end=' ')
					if group.typeCST == 'shift':
						group_name_shift = group.name
						ref_shift = float(group.ref_value)
						max_shift = group.maxCST
						min_shift = group.minCST
					elif group.typeCST == 'amp':
						group_name_amp = group.name
						ref_amp = float(group.ref_value)
						ref_T2 = group.ref_T2
						max_amp = group.maxCST
						min_amp = group.minCST
					elif group.typeCST == 'phase':
						group_name_phase = group.name
						ref_phase = float(sp.deg2rad(group.ref_value))
						max_phase = group.maxCST
						min_phase = group.minCST
					elif group.typeCST == 'lipid_LWs':
						group_name_lw = group.name
						ref_lw = float(group.ref_value)
						max_lw = group.maxCST
						min_lw = group.minCST
					elif group.typeCST == 'lipid_WLs':
						group_name_wl = group.name
						ref_lw = float(group.ref_value)
						max_lw = group.maxCST
						min_lw = group.minCST

			print(str(ref_shift) + ',', end=' ')
			print(str(ref_amp)+ ',', end=' ')
			print(str(ref_phase))				

			# set flag to group metabolite peaks > minWaterPPMInput with a separate amplitude param
			minWaterPPMInput_FLAG = False

			# output metabolite into constraint file
			for j in range(0, np.size(metab.area)):
				pk_count = pk_count + 1

				# don't use any peaks > 5 ppm (or whatever is specified in the input field)
				if metab.ppm[j] > float(self.maxPPMInput.text()):
					output.append(';')

				output.append(str(pk_count) + '\t')

				rel_shift = float(metab.ppm[j]) - ref_shift

				if metab.T2 == 0 or ref_T2 == 0:
					rel_amp = (float(metab.area[j] * metab.A_m)/ref_amp)
				else:
					rel_amp = (float(metab.area[j] * metab.A_m)/ref_amp) * (np.exp(-experiment.TE/(metab.T2*1000)) / np.exp(-experiment.TE/(ref_T2*1000)))
				rel_phase = float(sp.deg2rad(metab.phase[j]))-ref_phase

				# NOTE: with shifts, FITMAN describes shifts in the positive ppm direction to be negative! So we must reverse the signs!
				if rel_shift > 0:
					output.append("{" + str(group_name_shift) + "}-" + "{0:.6f}".format(float(rel_shift)))
				elif rel_shift < 0:
					output.append("{" + str(group_name_shift) + "}+" + "{0:.6f}".format(float(rel_shift)).replace('-',''))
				else:
					output.append("{" + str(group_name_shift) + "}")
					if min_shift != -sp.inf:
						output.append(" >" + "{0:.6f}".format(float(min_shift)))
					if max_shift != sp.inf:
						output.append(" <" + "{0:.6f}".format(float(max_shift)))
				output.append('\t')

				# ACCOUNTING FOR MACROMOLECULE BASLINE MODEL
				if ('LWlmX' in group_name_lw):
					print(group_name_lw, float(metab.lw))
					rel_lw = float(metab.lw)-ref_lw
					
					if rel_lw == 1.0:
						output.append("{" + 'LWlmX' + "}")
					else:
						output.append("{" + 'LWlmX' + "}+" + "{0:.6f}".format(rel_lw))

					if min_lw != -sp.inf:
						output.append(" >" + "{0:.6f}".format(min_lw))
					if max_lw != +sp.inf:
						output.append(" <" + "{0:.6f}".format(max_lw))
				else:
					if float(experiment.LW_limit) == 0:
						output.append("{" + "LW" + "} >0.000000" if experiment.LW_linkall else "{" + str(group_name_lw) + "} >0.000000")
					else:
						output.append("{" + "LW" + "} >0.000000 <" + "{0:.6f}".format(float(experiment.LW_limit)) if experiment.LW_linkall else "{" + str(group_name_lw) + "} >0.000000 <" + "{0:.6f}".format(float(experiment.LW_limit)))
				output.append('\t')

				if rel_amp == 1:
					output.append("{" + str(group_name_amp) + "}" + " >" + "{0:.6f}".format(float(0)))
					if min_amp != -sp.inf:
						output.append(" > {0:.6f}".format(float(min_amp)) if float(min_amp) > 0 else " ")
					if max_amp != sp.inf:
						output.append(" < {0:.6f}".format(float(max_amp)))
				else:
					if metab.ppm[j] > float(self.minWaterPPMInput.text()):
						if self.linkWaterRemovalAmpsCheckBox.isChecked():
							# Let all peaks greater than 4.1 ppm have their own amplitude. This is because HSVD water removal will affect these peaks.
							if minWaterPPMInput_FLAG == False:
								ref_amp_new = metab.area[j]
								group_name_amp_new = group_name_amp + '_' + str(pk_count)
								output.append("{" + str(group_name_amp_new) + "}" + " >" + "{0:.6f}".format(float(0)))
								if min_amp != -sp.inf:
									output.append(" > {0:.6f}".format(float(min_amp)) if float(min_amp) > 0 else " ")
								if max_amp != sp.inf:
									output.append(" < {0:.6f}".format(float(max_amp)))
								minWaterPPMInput_FLAG = True
							else:
								output.append("{" + str(group_name_amp_new) + "}*" + "{0:.6f}".format(float(metab.area[j]/ref_amp_new)) + " >" + "{0:.6f}".format(float(0)))
						else:
							group_name_amp_new = group_name_amp + '_' + str(pk_count)
							output.append("{" + str(group_name_amp_new) + "}" + " >" + "{0:.6f}".format(float(0)))
							if min_amp != -sp.inf:
								output.append(" > {0:.6f}".format(float(min_amp)) if float(min_amp) > 0 else " ")
							if max_amp != sp.inf:
								output.append(" < {0:.6f}".format(float(max_amp)))
					elif metab.name == 'NAA' and experiment.splitNAA == True:
						if float(metab.ppm[j]) > 2.1 and float(metab.ppm[j]) < 2.6:
							ref_amp_new = metab.area[[index for index,value in enumerate(metab.ppm) if value > 2.1 and value < 2.6][0]]*metab.A_m
							rel_amp = float(metab.area[j] * metab.A_m)/ref_amp_new
							output.append("{" + str(group_name_amp) + "1}*" + "{0:.6f}".format(float(rel_amp)) + " >" + "{0:.6f}".format(float(0)))
						elif float(metab.ppm[j]) > 2.6 and float(metab.ppm[j]) < 3.5:
							ref_amp_new = metab.area[[index for index,value in enumerate(metab.ppm) if value > 2.6 and value < 3.5][0]]*metab.A_m
							rel_amp = float(metab.area[j] * metab.A_m)/ref_amp_new
							output.append("{" + str(group_name_amp) + "2}*" + "{0:.6f}".format(float(rel_amp)) + " >" + "{0:.6f}".format(float(0)))
						elif float(metab.ppm[j]) > 3.5:
							ref_amp_new = metab.area[[index for index,value in enumerate(metab.ppm) if value > 3.5][0]]*metab.A_m
							rel_amp = float(metab.area[j] * metab.A_m)/ref_amp_new
							output.append("{" + str(group_name_amp) + "3}*" + "{0:.6f}".format(float(rel_amp)) + " >" + "{0:.6f}".format(float(0)))
					else:
						output.append("{" + str(group_name_amp) + "}*" + "{0:.6f}".format(float(rel_amp)) + " >" + "{0:.6f}".format(float(0)))
				output.append('\t')

				if rel_phase > 0:
					output.append("{" + str(group_name_phase) + "}+" + "{0:.6f}".format(float(rel_phase)))
				elif rel_phase < 0:
					output.append("{" + str(group_name_phase) + "}" + "{0:.6f}".format(float(rel_phase)))
				else:
					output.append("{" + str(group_name_phase) + "}")
					if min_phase != -sp.inf:
						output.append(" >" + "{0:.6f}".format(float(min_phase)))
					if max_phase != sp.inf:
						output.append(" <" + "{0:.6f}".format(float(max_phase)))
				output.append('\t')

				# output.append("@{DT}" + '\t')
				output.append("{DT}" + '\t')

				# ACCOUNTING FOR MM/LIPIDS WITH GAUSSIAN LINE SHAPES (ADD THIS IN AN IMPORT TAB LATER)
				if ('lip3' in group_name_wl) or ('lip4' in group_name_wl) or ('lip5' in group_name_wl) or ('mm2' in group_name_wl) or ('mm3' in group_name_wl) or ('mm4' in group_name_wl):
					output.append("{" + str(group_name_wl) + "} >10.000000 <50.000000")
				elif ('WLlipX' in group_name_wl):
					print(group_name_wl, float(metab.lw))
					rel_lw = float(metab.lw)-ref_lw
					if rel_lw == 1.0:
						output.append("{" + 'WLlipX' + "} >" + "{0:.6f}".format(min_lw) + " <" + "{0:.6f}".format(max_lw))
					else:
						output.append("{" + 'WLlipX' + "}+" + "{0:.6f}".format(rel_lw) + " >" + "{0:.6f}".format(min_lw) + " <" + "{0:.6f}".format(max_lw))
				elif ('lipmm1' in group_name_wl) or ('lip2' in group_name_wl):
					output.append("@{" + str(group_name_wl) + "}")
				else:
					output.append("@{WL}")
				output.append('\t')

				# output.append("@{WL}" + '\t')
				
				output.append('\n')
		
		output.append("****_Constraints_File_Ends_****")

		for line in output:
			if 'number_peaks' in line:
				cst_file.write(line.replace('<PLACEHOLDER>', str(pk_count)))
				print(str(pk_count))
			else:
				cst_file.write(line)
		cst_file.close()

		self.cstFilenameInfoLabel.append(">> Constraint file successfully output to: \n" + str(self.cstFilenameInput.text()))

	# ----- METHODS FOR 'GENERATE GUESS FILE' TAB ----- #
	def chooseDATFile(self):
		self.filenameConfirmButton_dat.setEnabled(False)
		prev = str(self.filenameInput_dat.text())
		self.filenameInput_dat.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open *.DAT File', '/home', 'DAT files (*.dat)')[0]))
		self.filenameConfirmButton_dat.setEnabled(True)
		if str(self.filenameInput_dat.text()) == '':
			self.filenameInput_dat.setText(str(prev))
			if str(prev) == '':
				self.filenameConfirmButton_dat.setEnabled(False)

	def loadDATFile(self):
		try:
			filename = str(self.filenameInput_dat.text())
			in_file = open(filename, 'r')

			n = 0
			fs = 0
			b0 = 0
			dat = []

			for i, line in enumerate(in_file):
				if i == 0:
					n = int(line.replace(' ', '').replace('\n', ''))//2
				elif i == 2:
					fs = 1/float(line.replace(' ', '').replace('\n', ''))
				elif i == 3:
					b0 = float(line.replace(' ', '').replace('\n', ''))
				elif i == 11:
					print(line)
				elif i > 11:
					dat.append(float(line.replace(' ', '').replace('\n', '')))

			real = sp.array(dat)[0::2]
			imag = sp.array(dat)[1::2]

			t = np.arange(0, n, 1) * (1/fs)
			
			self.ref_signal = RefSignal(real+1j*imag, n, fs, t, b0)
			print('n: ' + str(n) + ', fs: ' + str(fs) + ', dt: ' + str(1/fs) + ', b0: ' + str(b0))

			self.filenameInfoLabel_dat.setText("Acquired data loaded! \n >> " + str(self.filenameInput_dat.text()).rsplit('/', 1)[-1])
			self.plotGESButton.setEnabled(True)
		except Exception as e:
			self.filenameInfoLabel_dat.setText("ERROR: " + str(s) + "\n >> Data could not be loaded.\n >> Please try another dat file.")

	def plotGES(self):
		try:
			# REFERENCE
			n_ref = self.ref_signal.n
			
			print('ref_n: ' + str(n_ref))
			print('ref_signal: ' + str(self.ref_signal.signal))
			
			f_ref = np.arange(+n_ref//2,-n_ref//2,-1)*(self.ref_signal.fs/n_ref)*(1/self.ref_signal.b0)
			spectra_ref = fftw.fftshift(fftw.fft(self.ref_signal.signal))

			
			# GUESS

			sfactor = float(self.sFactorLineEdit.text())
			afactor = float(self.aFactorLineEdit.text())
			dfactor = float(self.dFactorLineEdit.text())	# sec
			wfactor = float(self.wFactorLineEdit.text())	# Hz
			pfactor = float(self.pFactorLineEdit.text())	# rad

			w0 = 0
			phi0 = 0

			signals = np.empty([np.size(self.metabListWidget.selectedItems()), np.size(self.experiment.t)], dtype=complex)
			i = 0

			# ---- Get from selected items from 'Plot Metabolites' tab ----
			# metabolites = []
			# for item in self.metabListWidget.selectedItems():
			# 	metabolites.append(str(item.text()))
			# print metabolites

			# ---- Get from shift constraint groups ---- #
			metabolites = []
			for group in self.experiment.cst_groups:
				if group.typeCST == 'shift':
					for member in group.members:
						if not(member in metabolites):
							metabolites.append(member)
			print(metabolites)
		
			w0 =   0
			phi0 = 0
			lb =   wfactor      
			
			signals = np.empty([np.size(metabolites), np.size(self.experiment.t)], dtype=complex)
			for i in range(0, np.size(metabolites)):
				signals[i,] = np.exp(-1j * (w0 * self.experiment.t) + 1j * (phi0)) \
							  * self.experiment.data[metabolites[i]].getFID(self.experiment.TE, self.experiment.b0, self.experiment.t, sfactor, afactor, pfactor, dfactor, lb)
			total_signal = np.sum(signals, axis=0)

			n = np.size(total_signal)
			
			print(n)
			print(total_signal)
			
			f = np.arange(+n//2,-n//2,-1)*(self.experiment.fs/n)*(1/self.experiment.b0)
			spectra = fftw.fftshift(fftw.fft(total_signal))

			print(np.size(f))
			print(np.size(np.real(spectra)))

			plt.figure(3)
			plt.clf()
			
			print(np.size(f_ref))
			print(np.size(np.real(spectra_ref)))

			plt.plot(f_ref, np.real(spectra_ref), label='Acquired Data') # 7.5 is a fudge factor to get the amplitudes in line with fitMAN
			plt.plot(-f, np.real(spectra), label='Guess') 
			plt.title('Spectra (Real)')
			plt.xlim(5-sfactor,0-sfactor)
			plt.xlabel('ppm')
			plt.legend()

			self.canvas[2].draw()

			self.gesFilenameConfirmButton.setEnabled(True)
			self.genGESButton.setEnabled(True)
		except Exception as e:
			self.gesFilenameInfoLabel.append('\n>> ERROR: ' + str(e) + '\n>> Could not plot GES. Please try again.\n')

	def setGESFile(self):
		self.gesFilenameInfoLabel.append(">> Output .ges File \n" + str(self.gesFilenameInput.text()))

	def genGES(self):
		try:
			print('----- Generating .ges file from current experiment -----')
			
			experiment = self.experiment

			self.savetoGES(self.experiment)
		except Exception as e:
			self.gesFilenameInfoLabel.append("\n>> ERROR: " + str (e) + '\n>> Could not output guess file. Please try again.\n')

	def savetoGES(self, experiment):
		filename_out = str(self.gesFilenameInput.text())
		ges_file = open(filename_out, 'w')

		pk_count = 0
		output = []

		output.append("****_Guess_File_Begins_**** <--- Do not remove this line" + '\n')
		output.append("[Parameters]" + '\n')
		output.append("number_peaks\t\t" + '<PLACEHOLDER>' + '\n')
		output.append("shift_units\t\tppm" + '\n')

		output.append('\n')

		output.append("[Variables]" + '\n')
		output.append('sfactor ' + "{0:.6f}".format(float(self.sFactorLineEdit.text())) + '\n')
		output.append('afactor ' + "{0:.6f}".format(float(self.aFactorLineEdit.text())) + '\n')
		output.append('dfactor ' + "{0:.6f}".format(float(self.dFactorLineEdit.text())) + '\n')
		# output.append('dfactor -' + str(experiment.dwell_time*3) + '\n')
		output.append('wfactor ' + "{0:.6f}".format(float(self.wFactorLineEdit.text())) + '\n')
		output.append('pfactor ' + "{0:.6f}".format(float(self.pFactorLineEdit.text())) + '\n')

		output.append('\n')

		output.append("[Peaks]" + '\n')
		output.append(";PEAK#\tSHIFT\tLORENTZIAN WIDTH\tAMPLITUDE\tPHASE\tDELAY\tGAUSSIAN WIDTH" + '\n')

		# ---- Get from selected items from 'Plot Metabolites' tab ----
		# metabolites = []
		# for item in self.metabListWidget.selectedItems():
		# 	metabolites.append(str(item.text()))
		# print metabolites

		# ---- Get from shift constraint groups ---- #
		metabolites = []
		for group in experiment.cst_groups:
			if group.typeCST == 'shift':
				for member in group.members:
					if not(member in metabolites):
						metabolites.append(member)
		print(metabolites)

		for i in range(0, np.size(metabolites)):
			
			# get metabolite
			metab = experiment.data[metabolites[i]]
			print('')
			print(metab.name + ',', end=' ')

			output.append(";" + str(metab.name) + '\n')

			# set default reference values
			ref_shift = float(metab.ppm[0])
			ref_amp = float(metab.area[0] * metab.A_m)
			ref_T2 = float(metab.T2)
			ref_phase = float(sp.deg2rad(metab.phase[0]))

			# set default max/min values
			min_shift = -sp.inf
			max_shift = sp.inf

			min_amp = -sp.inf
			max_amp = sp.inf

			min_phase = -sp.inf
			max_phase = sp.inf

			group_name_lw = "LW" + str(metab.name_short())
			group_name_wl = "WL" + str(metab.name_short())			

			# check if metabolite is part of a constraint group
			# and override defaults with group parameters
			for group in experiment.cst_groups:
				if metab.name in group.members:
					print(group.name + ',', end=' ')					
					if group.typeCST == 'shift':
						group_name_shift = group.name
						ref_shift = float(group.ref_value)
						max_shift = group.maxCST
						min_shift = group.minCST
					elif group.typeCST == 'amp':
						group_name_amp = group.name
						ref_amp = float(group.ref_value)
						ref_T2 = group.ref_T2
						max_amp = group.maxCST
						min_amp = group.minCST
					elif group.typeCST == 'phase':
						group_name_phase = group.name
						ref_phase = float(sp.deg2rad(group.ref_value))
						max_phase = group.maxCST
						min_phase = group.minCST

			print(str(ref_shift) + ',', end=' ')
			print(str(ref_amp)+ ',', end=' ')
			print(str(ref_phase)+ ',', end=' ')

			# set flag to group metabolite peaks > minWaterPPMInput with a separate amplitude param
			minWaterPPMInput_FLAG = False

			# output metabolite into guess file
			for j in range(0, np.size(metab.area)):
				pk_count = pk_count + 1

				# don't use any peaks > 5 ppm (or whatever is specified in the input field)
				if metab.ppm[j] > float(self.maxPPMInput.text()):
					output.append(';')

				output.append(str(pk_count) + '\t')

				rel_shift = float(metab.ppm[j]) - ref_shift

				if metab.T2 == 0 or ref_T2 == 0:
					rel_amp = (float(metab.area[j] * metab.A_m)/ref_amp)
				else:
					rel_amp = (float(metab.area[j] * metab.A_m)/ref_amp) * (np.exp(-experiment.TE/(metab.T2*1000)) / np.exp(-experiment.TE/(ref_T2*1000)))
				rel_phase = float(sp.deg2rad(metab.phase[j]))-ref_phase

				rel_phase = float(sp.deg2rad(metab.phase[j]))-ref_phase

				# NOTE: with shifts, FITMAN describes shifts in the positive ppm direction to be negative! So we must reverse the signs!
				if rel_shift > 0:
					output.append("0")
				elif rel_shift < 0:
					output.append("0")
				else:
					output.append("-{0:.6f}".format(float(metab.ppm[j])) + '+sfactor')
				output.append('\t')

				# --- Legacy code that still works (lipmm1 ... mm4 defined the old macromolecule model) --- #
				# ACCOUNTING FOR MM/LIPIDS WITH LORENTZIAN LINE SHAPES (ADD THIS IN AS AN IMPORT TAB LATER)
				if ('lipmm1' in group_name_lw) or ('lip2' in group_name_lw):
					output.append(str(metab.lw))
				elif ('lip3' in group_name_lw) or ('lip4' in group_name_lw) or ('lip5' in group_name_lw) or ('mm2' in group_name_lw) or ('mm3' in group_name_lw) or ('mm4' in group_name_lw):
					output.append("0")
				else:
					output.append("0+wfactor")
				output.append('\t')

				if metab.name == 'NAA' and experiment.splitNAA == True:
					if float(metab.ppm[j]) > 2.1 and float(metab.ppm[j]) < 2.6:
						ref_amp_new = metab.area[[index for index,value in enumerate(metab.ppm) if value > 2.1 and value < 2.6][0]]*metab.A_m
						rel_amp = float(metab.area[j] * metab.A_m)/ref_amp_new
					elif float(metab.ppm[j]) > 2.6 and float(metab.ppm[j]) < 3.5:
						ref_amp_new = metab.area[[index for index,value in enumerate(metab.ppm) if value > 2.6 and value < 3.5][0]]*metab.A_m
						rel_amp = float(metab.area[j] * metab.A_m)/ref_amp_new
					elif float(metab.ppm[j]) > 3.5:
						ref_amp_new = metab.area[[index for index,value in enumerate(metab.ppm) if value > 3.5][0]]*metab.A_m
						rel_amp = float(metab.area[j] * metab.A_m)/ref_amp_new

				if rel_amp == 1:
					output.append("{0:.6f}".format(float(metab.area[j] * metab.A_m)) + '*afactor')
				else:
					if metab.ppm[j] > float(self.minWaterPPMInput.text()):
						if self.linkWaterRemovalAmpsCheckBox.isChecked():
							# Let all peaks greater than 4.1 ppm have their own amplitude. This is because HSVD water removal will affect these peaks.
							if minWaterPPMInput_FLAG == False:
								output.append("{0:.6f}".format(float(metab.area[j] * metab.A_m)) + '*afactor')
								minWaterPPMInput_FLAG = True
							else:
								output.append("0")
						else:
							output.append("{0:.6f}".format(float(metab.area[j] * metab.A_m)) + '*afactor')
					else:
						output.append("0")
				output.append('\t')

				if rel_phase > 0:
					output.append("0 \t 0")
				elif rel_phase < 0:
					output.append("0 \t 0")
				else:
					output.append("{0:.6f}".format(float(sp.deg2rad(metab.phase[j]))) + '+pfactor \t 0+dfactor')
					
					# - IGNORE ABOSLUTE PHASE DATA FROM SIMULATION
					# - RELATIVE PHASE DATA CAN BE KEPT BY LINKING THE PHASE OF ALL METABOLITES TOGETHER (RECOMMENDED)
					# - PFACTOR IS SET IN GUI ... PFACTOR SHOULD EQUAL THE ZERO ORDER PHASE CORRECTION PERFORMED IN FITMAN
					# output.append('0+pfactor \t 0+dfactor')
				output.append('\t')

				# ACCOUNTING FOR MM/LIPIDS WITH GAUSSIAN LINE SHAPES (ADD THIS IN AN IMPORT TAB LATER)
				if ('lip3' in group_name_wl) or \
				   ('lip4' in group_name_wl) or \
				   ('lip5' in group_name_wl) or \
				   ('mm2' in group_name_wl)  or \
				   ('mm3' in group_name_wl)  or \
				   ('mm4' in group_name_wl):
					output.append(str(metab.lw))
				else:
					output.append("0")
				output.append('\t')

				output.append('\n')
		
		output.append("****_Guess_File_Ends_****  <--- Do not remove this line")

		for line in output:
			if 'number_peaks' in line:
				ges_file.write(line.replace('<PLACEHOLDER>', str(pk_count)))
				print(str(pk_count))
			else:
				ges_file.write(line)
		ges_file.close()
		self.gesFilenameInfoLabel.append(">> Successfully output to .ges file \n" + str(self.gesFilenameInput.text()))


# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())