from __future__ import print_function

# ---- System Libraries ---- #
from builtins import zip
from builtins import str
from builtins import range
import sys
import errno
import os
import copy
import traceback

# Check for PyQt6 (for native mac M1 compatibility), otherwise continue using PyQt5
import importlib
PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore, QtGui, QtWidgets, uic
else:
	from PyQt5 import QtCore, QtGui, QtWidgets, uic

# ---- Math Libraries ---- #
import numpy as np

# ---- Plotting Libraries ---- #
import matplotlib as mpl;
mpl.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (
	FigureCanvasQTAgg as FigureCanvas,
	NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt

# ---- Data Classes ---- #
from magiqdataclasses import *

# ---- Pre-Processing Functions ---- #
from preproc import *

qtCreatorFile = "apps/ui/apps.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

# ---- Main Application Class ---- #
class MyApp(QtWidgets.QWidget, Ui_MainWindow):

	# ---- Methods to Set Up UI ---- #
	def __init__(self):
		'''
			Method initializes the application 
			and binds methods to UI buttons.
		'''
		QtWidgets.QWidget.__init__(self)
		Ui_MainWindow.__init__(self)
		self.setupUi(self)

		self.workingDirectory_bruker = os.path.expanduser('~')
		self.workingDirectory_wr  = os.path.expanduser('~')
		self.workingDirectory_mmr = os.path.expanduser('~')
		
		# Setup for plotting
		self.canvas = [None]*3
		self.toolbar = [None]*3

		# Bind buttons to methods in each tab
		self.setBindings('Water Removal')
		self.setBindings('Macromolecule Removal')
		self.setBindings('Bruker File Conversion')

	def setBindings(self, tab):
		'''
			Method binds methods to UI buttons in each UI tab.
		'''
		if tab == 'Water Removal':

			self.filenameBrowseButton_dat.clicked.connect(self.chooseDatFile_wr)
			self.filenameConfirmButton_dat.clicked.connect(self.loadDatFile_wr)
			self.filenameConfirmButton_dat.setEnabled(False)

			self.runWaterRemovalButton.clicked.connect(self.runWaterRemoval)
			self.runWaterRemovalButton.setEnabled(False)

			self.saveWaterRemovalButton.clicked.connect(self.saveWaterRemoval)
			self.saveWaterRemovalButton.setEnabled(False)

		elif tab == 'Macromolecule Removal':

			self.filenameBrowseButton_fulldat_mmr.clicked.connect(self.chooseFullDatFile_mmr)
			self.filenameConfirmButton_fulldat_mmr.clicked.connect(self.loadFullDatFile_mmr)
			self.filenameConfirmButton_fulldat_mmr.setEnabled(False)

			self.filenameBrowseButton_mmdat_mmr.clicked.connect(self.chooseMMDatFile_mmr)
			self.filenameConfirmButton_mmdat_mmr.clicked.connect(self.loadMMDatFile_mmr)
			self.filenameConfirmButton_mmdat_mmr.setEnabled(False)

			self.runMMRemovalButton.clicked.connect(self.runMMRemoval)
			self.runMMRemovalButton.setEnabled(False)

			self.saveMMRemovalButton.clicked.connect(self.saveMMRemoval)
			self.saveMMRemovalButton.setEnabled(False)

		elif tab == 'Bruker File Conversion':

			self.inputFilenameButton_bruker.clicked.connect(self.chooseInputFile_bruker)
			self.referenceFilenameButton_bruker.clicked.connect(self.chooseRefFile_bruker)

			self.runConversionButton_bruker.clicked.connect(self.runConversion_bruker)
			self.runConversionButton_bruker.setEnabled(False)

		self.setPlot(tab)
	
	# ---- Bruker File Conversion ---- #
	def chooseInputFile_bruker(self):
		'''
			Method allows you to choose a Bruker dataset and 
			creates a 'converted' folder within the dataset directory.
		'''
		prev = str(self.inputFilenameInput_bruker.text())
		
		self.inputFilenameInput_bruker.setText(str(QtWidgets.QFileDialog.getExistingDirectory(self, 'Open Bruker Data Directory', self.workingDirectory_bruker)))
		
		if str(self.inputFilenameInput_bruker.text()) == '':
			if str(prev) == '':
				self.inputFilenameInput_bruker.setText(str(prev))
				self.runConversionButton_bruker.setEnabled(False)
		else:
			if str(self.outputFilenameInput_bruker.text()) == '' or  str(self.referenceFilenameInput_bruker.text()) == '':
				self.runConversionButton_bruker.setEnabled(False)
			else:
				self.runConversionButton_bruker.setEnabled(True)
			self.outputFilenameInput_bruker.setText(self.inputFilenameInput_bruker.text() + '/converted/')

			try:
				os.mkdir(str(self.outputFilenameInput_bruker.text()))
			except OSError as exc:
				msg = QtWidgets.QMessageBox()
				if exc.errno == errno.EEXIST:
					msg.setIcon(QtWidgets.QMessageBox.Warning)
					msg.setText('Warning')
					msg.setInformativeText('"' + self.inputFilenameInput_bruker.text() + '/converted/" already exists.')
					msg.setWindowTitle('Warning')
				elif exc.errno == errno.EROFS:
					msg.setIcon(QtWidgets.QMessageBox.Critical)
					msg.setText('Error')
					msg.setInformativeText('Read-only filesystem error occurred while making "' + self.inputFilenameInput_bruker.text() + '/converted/"')
					msg.setWindowTitle('Error')
				else:
					msg.setIcon(QtWidgets.QMessageBox.Critical)
					msg.setText('Error ' + str(exc) + 'has occurred!')
					msg.setInformativeText(traceback.format_exc())
					msg.setWindowTitle('Error')
					msg.exec_()
					raise exc
				msg.exec_()
			self.workingDirectory_bruker = os.path.abspath(os.path.join(os.path.expanduser(str(prev)), os.pardir))

	def chooseRefFile_bruker(self):
		'''
			Method allows you to choose a Bruker dataset directory.
		'''
		prev = str(self.referenceFilenameInput_bruker.text())
		self.referenceFilenameInput_bruker.setText(str(QtWidgets.QFileDialog.getExistingDirectory(self, 'Open Bruker Data Directory', self.workingDirectory_bruker)))
		if str(self.referenceFilenameInput_bruker.text()) == '':
			if str(prev) == '':
				self.referenceFilenameInput_bruker.setText(str(prev))
				self.runConversionButton_bruker.setEnabled(False)
		else:
			if str(self.outputFilenameInput_bruker.text()) == '' or  str(self.inputFilenameInput_bruker.text()) == '':
				self.runConversionButton_bruker.setEnabled(False)
			else:
				self.runConversionButton_bruker.setEnabled(True)
			
	def runConversion_bruker(self):
		'''
			This method runs the conversion process.
		'''
		self.conversionConsole_bruker.clear()

		out_name = str(self.outputFilenameInput_bruker.text())
		out_name_sup = out_name + 'corr_sup'
		out_name_uns = out_name + 'corr_uns'

		# Read suppressed file
		sup_file   = BrukerFID(str(self.inputFilenameInput_bruker.text()))
		self.conversionConsole_bruker.append('sup_file: ' + sup_file.file_dir)
		sup_file.print_params(self.conversionConsole_bruker)
		self.conversionConsole_bruker.append('')
		self.scaleFactorInput_bruker.setText(str(sup_file.ConvS))
		self.timeDelayInput_bruker.setText(str(float(sup_file.digShift) * 1/sup_file.fs))

		# Read unsuppressed file
		uns_file = BrukerFID(str(self.referenceFilenameInput_bruker.text()))
		self.conversionConsole_bruker.append('uns_file: ' + uns_file.file_dir)
		uns_file.print_params(self.conversionConsole_bruker)
		self.conversionConsole_bruker.append('')

		# Write files as fitMAN dat files.
		sup_file.writeDAT(out_name + 'raw_sup', '')
		uns_file.writeDAT(out_name + 'raw_uns', '')

		# Plot Suppressed File
		plt.figure(3)

		ax = plt.subplot(411)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		plt.title('Raw Signal')
		plt.plot(sup_file.t, np.real(sup_file.signal))

		ax = plt.subplot(412)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		f_sup, spec_sup = sup_file.getSpec()
		plt.plot(f_sup[0:sup_file.n], np.real(spec_sup[0:sup_file.n]))


		# Run baseline correction.
		if self.baselineCorrectionButton1_bruker.isChecked():
			sup_file.signal = baseline_corr(sup_file.signal)
			uns_file.signal = baseline_corr(uns_file.signal)
			out_name_sup = out_name_sup + '_bc'
			out_name_uns = out_name_uns + '_bc'

		# Run post-processing.
		out_file_sup = copy.deepcopy(sup_file)
		out_file_uns = copy.deepcopy(uns_file)
		if self.queccRadioButton_bruker.isChecked():
			# QUECC
			quecc_points = int(self.qualityPointsInput_bruker.text())
			out_file_sup.signal = np.array(quecc(sup_file.signal, uns_file.signal, quecc_points, sup_file.t))
			out_file_uns.signal = np.array(quecc(uns_file.signal, uns_file.signal, quecc_points, uns_file.t))
			out_file_sup.writeDAT(out_name_sup, 'quecc' + str(quecc_points))
			out_file_uns.writeDAT(out_name_uns, 'quecc' + str(quecc_points))
		elif self.qualityRadioButton_bruker.isChecked():
			# QUALITY
			out_file_sup.signal = np.array(quality(sup_file.signal, uns_file.signal))
			out_file_uns.signal = np.array(quality(uns_file.signal, uns_file.signal))
			out_file_sup.writeDAT(out_name_sup, 'quality')
			out_file_uns.writeDAT(out_name_uns, 'quality')
		elif self.eccRadioButton_bruker.isChecked():
			# ECC
			out_file_sup.signal = np.array(ecc(sup_file.signal, uns_file.signal))
			out_file_uns.signal = np.array(ecc(uns_file.signal, uns_file.signal))
			out_file_sup.writeDAT(out_name_sup, 'ecc')
			out_file_uns.writeDAT(out_name_uns, 'ecc')

		ax = plt.subplot(413)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		plt.title('Corrected Signal')
		plt.plot(out_file_sup.t, np.real(out_file_sup.signal))

		ax = plt.subplot(414)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)
		
		f_out, spec_out = out_file_sup.getSpec()
		plt.plot(f_out[0:out_file_sup.n], np.real(spec_out[0:out_file_sup.n]))

		plt.tight_layout()
		self.canvas[2].draw()

	# ---- Water Removal ---- #
	def runWaterRemoval(self):
		'''
			This method performs removal of the residual water signal 
			using Hankel Singular Value decomposition.
		'''
		
		#1. Fit specturm with HSVD.
		peak, width_L, ppm, area, phase = self.hsvd(self.dat,int(self.hsvdPointsLineEdit_wr.text()),float(self.hsvdRatioLineEdit_wr.text()),int(self.hsvdComponentsLineEdit_wr.text()),'water')
		hsvd_fit = Metabolite()

		# take only the peaks within the specified frequency range
		for (i, comp) in enumerate(peak):
			if (ppm[i]/self.dat.b0) >= float(self.XminLineEdit.text()) and (ppm[i]/self.dat.b0) <= float(self.XmaxLineEdit.text()):
				hsvd_fit.peak.append(peak[i])
				hsvd_fit.width_L.append(width_L[i])
				hsvd_fit.ppm.append(ppm[i]/self.dat.b0)
				hsvd_fit.area.append(area[i])
				hsvd_fit.phase.append(phase[i])

		hsvd_fid = hsvd_fit.getFID(self.dat.TE, self.dat.b0, self.dat.t, 0, 1, 0, 0, 0)

		self.dat_hsvd = copy.deepcopy(self.dat)
		self.dat_hsvd.signal = hsvd_fid

		ax = plt.subplot(211)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		# plot raw spectrum
		f_dat, spec_dat = self.dat.getSpec()
		plt.plot(f_dat[0:self.dat.n], np.real(spec_dat[0:self.dat.n]), label='Data')

		# plot fitted spectrum
		f_hsvd, spec_hsvd = self.dat_hsvd.getSpec()
		plt.plot(f_hsvd, np.real(spec_hsvd), label='Fit')

		# plot residual
		VSHIFT = 0.25
		plt.plot(f_hsvd, np.real(spec_dat[0:self.dat.n]) - np.real(spec_hsvd) + VSHIFT*np.amax(np.real(spec_hsvd)), label='Residual')

		# legend and title
		ax.legend(loc="upper left")
		plt.title('Original Spectrum')

		# 2. Subtract from original spectrum.
		scale = 1.0
		self.dat_wr = copy.deepcopy(self.dat)
		self.dat_wr.n = np.min([np.size(self.dat.signal,0), np.size(self.dat_hsvd.signal, 0)])
		self.dat_wr.signal = self.dat.signal[0:self.dat_wr.n] - scale*self.dat_hsvd.signal[0:self.dat_wr.n]

		ax = plt.subplot(212)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		f_wr, spec_wr = self.dat_wr.getSpec()
		plt.plot(f_wr, np.real(spec_wr))
		plt.title('Water Removed Spectrum')

		plt.tight_layout()
		self.canvas[0].draw()

		self.saveWaterRemovalButton.setEnabled(True)

	def saveWaterRemoval(self):
		'''
			This method saves the water removed signal as a *.dat file.
		'''

		self.dat_wr.filename = self.dat.filename.replace('.dat', '_wr.dat')

		out_file = open(self.dat_wr.filename, 'w')
		in_file  = open(self.dat.filename, 'r')

		for (i, line) in enumerate(in_file):
			if i > 11:
				for element in self.dat_wr.signal:
					out_file.write("{0:.6f}".format(float(np.real(element))) + '\n')
					out_file.write("{0:.6f}".format(float(np.imag(element))) + '\n')
				break
			else:
				out_file.write(line)

		out_file.close()
		in_file.close()

	# ---- Macromolecule Removal ---- #
	def runMMRemoval(self):
		'''
			This method performs macromolecule removing using Hankel
			singular value decomposition.
		'''

		# 1. Fit macromolecule spectrum with HSVD.
		peak, width_L, ppm, area, phase = self.hsvd(self.MMDat_mmr,int(self.hsvdPointsLineEdit_mmr.text()),float(self.hsvdRatioLineEdit_mmr.text()),int(self.hsvdComponentsLineEdit_mmr.text()),'mm')
		hsvd_fit = Metabolite()
		hsvd_fit.peak = peak
		hsvd_fit.width_L = width_L
		hsvd_fit.ppm = ppm / self.MMDat_mmr.b0
		hsvd_fit.area = area
		hsvd_fit.phase = phase
		hsvd_fid = hsvd_fit.getFID(self.MMDat_mmr.TE, self.MMDat_mmr.b0, self.MMDat_mmr.t, 0, 1, 0, 0, 0)

		self.MMDatHSVD_mmr = copy.deepcopy(self.MMDat_mmr)
		self.MMDatHSVD_mmr.signal = hsvd_fid

		ax = plt.subplot(312)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		# plot raw MM signal
		f_dat, spec_dat = self.MMDat_mmr.getSpec()
		plt.plot(f_dat[0:self.MMDat_mmr.n], np.real(spec_dat[0:self.MMDat_mmr.n]), label='Data')

		# plot fitted MM signal
		f_hsvd, spec_hsvd = self.MMDatHSVD_mmr.getSpec()
		plt.plot(f_hsvd, np.real(spec_hsvd), label='Fit')

		# plot residual
		VSHIFT = 1.0
		plt.plot(f_hsvd, np.real(spec_dat[0:self.MMDat_mmr.n]) - np.real(spec_hsvd) - VSHIFT*np.amax(np.real(spec_hsvd)), label='Residual')

		# legend and title
		ax.legend(loc="upper left")
		plt.title('MM Spectrum')

		# 2. Subtract fitted macromolecule spectrum from full spectrum.
		scale = 1.3 * (float(self.fullDat_mmr.ConvS) / float(self.MMDat_mmr.ConvS))
		print(self.fullDat_mmr.ConvS, self.MMDat_mmr.ConvS, (float(self.fullDat_mmr.ConvS) / float(self.MMDat_mmr.ConvS)), scale)
		
		self.metabDat_mmr = copy.deepcopy(self.fullDat_mmr)
		self.metabDat_mmr.n = np.min([np.size(self.fullDat_mmr.signal,0), np.size(self.MMDatHSVD_mmr.signal, 0)])
		self.metabDat_mmr.signal = self.fullDat_mmr.signal[0:self.metabDat_mmr.n] - scale*self.MMDatHSVD_mmr.signal[0:self.MMDatHSVD_mmr.n]

		ax = plt.subplot(313)
		ax.clear()
		ax.spines['top'].set_visible(False)    
		ax.spines['bottom'].set_visible(True)    
		ax.spines['right'].set_visible(False)    
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().set_visible(False)

		f_metab, spec_metab = self.metabDat_mmr.getSpec()
		plt.plot(f_metab, np.real(spec_metab))
		plt.title('Metabolite Spectrum')

		plt.tight_layout()
		self.canvas[1].draw()

		self.saveMMRemovalButton.setEnabled(True)

	def saveMMRemoval(self):
		'''
			This method saves the macromolecule removed result
			as a *.dat file. Additionally, it also saves the HSVD fit as a
			*.dat file.
		'''

		# 1. Save metabolite spectrum.
		self.metabDat_mmr.filename = self.fullDat_mmr.filename.replace(self.fullDat_mmr.filename.split('/')[-1], 'sup.dat')

		out_file = open(self.metabDat_mmr.filename, 'w')
		in_file  = open(self.fullDat_mmr.filename, 'r')

		for (i, line) in enumerate(in_file):
			if i > 11:
				for element in self.metabDat_mmr.signal:
					out_file.write("{0:.6f}".format(float(np.real(element))) + '\n')
					out_file.write("{0:.6f}".format(float(np.imag(element))) + '\n')
				break
			else:
				out_file.write(line)

		out_file.close()
		in_file.close()

		# 2. Save HSVD fit.
		self.MMDatHSVD_mmr.filename = self.MMDat_mmr.filename.replace('.dat', '_hsvd.dat')
		out_file = open(self.MMDatHSVD_mmr.filename, 'w')
		in_file  = open(self.MMDat_mmr.filename, 'r')

		for (i, line) in enumerate(in_file):
			if i > 11:
				for element in self.MMDatHSVD_mmr.signal:
					out_file.write("{0:.6f}".format(float(np.real(element))) + '\n')
					out_file.write("{0:.6f}".format(float(np.imag(element))) + '\n')
				break
			else:
				out_file.write(line)

		out_file.close()
		in_file.close()

	# ---- Methods to Load Files ---- #
	def chooseDatFile_wr(self):
		'''
			This method launches a dialog window allowing you to
			choose a *.dat file.
		'''
		self.filenameConfirmButton_dat.setEnabled(False)
		prev = str(self.filenameInput_dat.text())
		self.filenameInput_dat.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Dat File', self.workingDirectory_wr, 'Dat files (*.dat)')[0]))
		self.filenameConfirmButton_dat.setEnabled(True)
		if str(self.filenameInput_dat.text()) == '':
			self.filenameInput_dat.setText(str(prev))
			if str(prev) == '':
				self.filenameConfirmButton_dat.setEnabled(False)

	def loadDatFile_wr(self):
		'''
			This method loads in the specified *.dat file
			and plots the signal.
		'''
		try:
			# load dat file
			datFile  = str(self.filenameInput_dat.text())
			self.dat = DatFile(datFile)
			self.filenameInfoLabel_dat.setText(datFile.split('/')[-1] + " successfully loaded.")

			# plot
			plt.figure(1)

			ax = plt.subplot(211)
			ax.clear()
			ax.spines['top'].set_visible(False)    
			ax.spines['bottom'].set_visible(True)    
			ax.spines['right'].set_visible(False)    
			ax.spines['left'].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().set_visible(False)

			plt.tick_params(axis="both", which="both", bottom="on", top="off",    
				labelbottom="on", left="off", right="off", labelleft="off")

			f_dat, spec_dat = self.dat.getSpec()
			plt.plot(f_dat[0:self.dat.n], np.real(spec_dat[0:self.dat.n]))
			plt.title('Original Spectrum')
			
			plt.tight_layout()
			self.canvas[0].draw()

			# reset buttons
			self.filenameBrowseButton_dat.setEnabled(True)
			self.filenameConfirmButton_dat.setEnabled(False)

			self.runWaterRemovalButton.setEnabled(True)

		except Exception as e:
			traceback.print_exc()
			self.filenameInfoLabel_dat.setText("ERROR: " + str(e) + " >> Please try again.")
			
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()

	def chooseFullDatFile_mmr(self):
		self.filenameConfirmButton_fulldat_mmr.setEnabled(False)
		prev = str(self.filenameInput_fulldat_mmr.text())
		self.filenameInput_fulldat_mmr.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Dat File', self.workingDirectory_mmr, 'Dat files (*.dat)')[0]))
		self.filenameConfirmButton_fulldat_mmr.setEnabled(True)
		if str(self.filenameInput_fulldat_mmr.text()) == '':
			self.filenameInput_fulldat_mmr.setText(str(prev))
			if str(prev) == '':
				self.filenameConfirmButton_fulldat_mmr.setEnabled(False)

		fullDatFile_mmr     = str(self.filenameInput_fulldat_mmr.text())
		self.workingDirectory_mmr = fullDatFile_mmr.replace(fullDatFile_mmr.split('/')[-1], '')

	def loadFullDatFile_mmr(self):
		'''
			This method loads in the specified *.dat file
			and plots the signal.
		'''
		try:
			# load dat file
			fullDatFile_mmr     = str(self.filenameInput_fulldat_mmr.text())
			self.fullDat_mmr = DatFile(fullDatFile_mmr)
			self.filenameInfoLabel_fulldat_mmr.setText(fullDatFile_mmr.split('/')[-1] + " successfully loaded.")

			# plot
			plt.figure(2)

			ax = plt.subplot(311)
			ax.clear()
			ax.spines['top'].set_visible(False)    
			ax.spines['bottom'].set_visible(True)    
			ax.spines['right'].set_visible(False)    
			ax.spines['left'].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().set_visible(False)

			plt.tick_params(axis="both", which="both", bottom="on", top="off",    
				labelbottom="on", left="off", right="off", labelleft="off")

			f_dat, spec_dat = self.fullDat_mmr.getSpec()
			plt.plot(f_dat[0:self.fullDat_mmr.n], np.real(spec_dat[0:self.fullDat_mmr.n]))
			plt.title('Full Spectrum')
			
			plt.tight_layout()
			self.canvas[1].draw()
			
			# reset buttons
			self.filenameBrowseButton_fulldat_mmr.setEnabled(True)
			self.filenameConfirmButton_fulldat_mmr.setEnabled(False)

		except Exception as e:
			traceback.print_exc()
			self.filenameInfoLabel_fulldat_mmr.setText("ERROR: " + str(e) + " >> Please try again.")

			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()

	def chooseMMDatFile_mmr(self):
		'''
			This method launches a dialog window allowing you to
			choose a *.dat file.
		'''
		self.filenameConfirmButton_mmdat_mmr.setEnabled(False)
		prev = str(self.filenameInput_mmdat_mmr.text())
		self.filenameInput_mmdat_mmr.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Dat File', self.workingDirectory_mmr, 'Dat files (*.dat)')[0]))
		self.filenameConfirmButton_mmdat_mmr.setEnabled(True)
		if str(self.filenameInput_mmdat_mmr.text()) == '':
			self.filenameInput_mmdat_mmr.setText(str(prev))
			if str(prev) == '':
				self.filenameConfirmButton_mmdat_mmr.setEnabled(False)

		MMDatFile_mmr     = str(self.filenameInput_mmdat_mmr.text())
		self.workingDirectory_mmr = MMDatFile_mmr.replace(MMDatFile_mmr.split('/')[-1], '')

	def loadMMDatFile_mmr(self):
		'''
			This method loads in the specified *.dat file
			and plots the signal.
		'''
		try:
			# load dat file
			MMDatFile_mmr     = str(self.filenameInput_mmdat_mmr.text())
			self.MMDat_mmr = DatFile(MMDatFile_mmr)
			self.filenameInfoLabel_mmdat_mmr.setText(MMDatFile_mmr.split('/')[-1] + " successfully loaded.")

			# plot
			plt.figure(2)

			ax = plt.subplot(312)
			ax.clear()
			ax.spines['top'].set_visible(False)    
			ax.spines['bottom'].set_visible(True)    
			ax.spines['right'].set_visible(False)    
			ax.spines['left'].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().set_visible(False)

			plt.tick_params(axis="both", which="both", bottom="on", top="off",    
				labelbottom="on", left="off", right="off", labelleft="off")

			f_dat, spec_dat = self.MMDat_mmr.getSpec()
			plt.plot(f_dat[0:self.MMDat_mmr.n], np.real(spec_dat[0:self.MMDat_mmr.n]))
			plt.title('MM Spectrum')

			plt.tight_layout()
			self.canvas[1].draw()

			# reset buttons
			self.filenameBrowseButton_mmdat_mmr.setEnabled(True)
			self.filenameConfirmButton_mmdat_mmr.setEnabled(False)

			self.runMMRemovalButton.setEnabled(True)

		except Exception as e:
			traceback.print_exc()
			self.filenameInfoLabel_mmdat_mmr.setText("ERROR: " + str(e) + " >> Please try again.")

			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()

	# ---- Methods for Plotting ---- #
	def setPlot(self, tab):
		'''
			This method sets references to the matplotlib figures.
		'''
		if tab == 'Water Removal':
			fig = plt.figure(1)
			self.addmpl(0, fig, self.plotResult_mplvl_wr)
		elif tab == 'Macromolecule Removal':
			fig = plt.figure(2)
			self.addmpl(1, fig, self.plotResult_mplvl_mmr)
		elif tab == 'Bruker File Conversion':
			fig = plt.figure(3)
			self.addmpl(2, fig, self.plotResult_mplvl_bruker)

	def addmpl(self, canvas_index, fig, vertical_layout):
		'''
			This method adds to matplotlib Figure Canvas and Toolbar
			to the UI.
		'''
		self.canvas[canvas_index] = FigureCanvas(fig)
		vertical_layout.addWidget(self.canvas[canvas_index])
		self.canvas[canvas_index].draw()

		self.toolbar[canvas_index] = NavigationToolbar(self.canvas[canvas_index], self, coordinates=True)
		vertical_layout.addWidget(self.toolbar[canvas_index])

	# ---- Methods for HSVD Fitting ---- #
	def lorentzian(self, time_axis, frequency, phase, fwhm):
		'''
			This method returns a time-domain lorentzian function.
				Inputs:	
					time_axis	discrete time array
					frequency	frequency of the Lorentzian function
					phase		phase of the Lorentzian function
					fwhm		Lorentzian linewidth
				
				Outputs:
					time-domian lorentzian function
		'''
		oscillatory_term = np.exp(1j * (2 * np.pi * frequency * time_axis + phase))
		damping = np.exp(-time_axis * np.pi * fwhm)
		fid = oscillatory_term * damping
		return fid / len(time_axis)

	def hsvd(self, dat, n, ratio, comp, mode):
		'''
			This method's code is based (heavily) on: https://github.com/openmrslab/suspect/blob/master/suspect/processing/water_suppression.py
			It performs an HSVD decomposition and returns a list of all parameters for each Lorentzian component.
			
				Inputs:
					dat		DatFile object containing the signal on which to perform the decomposition 
					n		number of points
					ratio	Hankel matrix row/col ratio
					comp	number of signal related singular values

				Outputs:
					An array containing the following elements:
					[0] array of component numbers
					[1] array of the damping coefficient of each component
					[2] array of the frequency coeffienct of each component
					[3] array of the amplitude of each component
					[4] array of the phase of each component
		'''
		
		if mode == 'mm':
			console = self.hsvdFitConsole_mmr
		elif mode == 'water':
			console = self.hsvdFitConsole_wr

		cols = n/(ratio+1)
		rows = n-cols
		L = int(np.ceil(rows))
		
		# build the Hankel matrix
		hankel_matrix = np.zeros((L, n-L+1), "complex")
		console.clear()
		console.append('Fitting:\t' + str(dat.filename))
		console.append('Creating Hankel matrix:\t[' + str(np.size(hankel_matrix, 0)) + 'x' + str(np.size(hankel_matrix, 1)) + ']')
		console.append('Points:\t' + str(n))
		console.append('Ratio:\t' + str(ratio))
		console.append('Components:\t' + str(comp))
		console.append('')
		for i in range(int(n-L)):
			hankel_matrix[:, i] = dat.signal[i:(i + L)]
		
		# perform the singular value decomposition
		U, s, V = np.linalg.svd(np.matrix(hankel_matrix))
		V = V.H # numpy returns the Hermitian conjugate of V
		
		# truncate the matrixes to the given rank (number of components)
		U_K = U[:, :comp]
		V_K = V[:, :comp]
		s_K = np.matrix(np.diag(s[:comp]))
		
		# because of the structure of the Hankel matrix, each row of U_K is the
		# result of multiplying the previous row by the delta t propagator matrix
		# Z' (a similar result holds for V as well). This gives us U_Kb * Z' = U_Kt
		# where U_Kb is U_K without the bottom row and U_Kt is U_K without the top
		# row.
		U_Kt = U_K[1:, :]
		U_Kb = U_K[:-1, :]
		
		# this gives us a set of linear equations which can be solved to find Z'.
		# Because of the noise in the system we solve with least-squares
		Zp = np.linalg.inv(U_Kb.H * U_Kb) * U_Kb.H * U_Kt
		
		# in the right basis, Zp is just the diagonal matrix describing the
		# evolution of each frequency component, by diagonalising the matrix we can
		# find that basis and get the z = exp((-damping + j*2pi * f) * dt) terms
		
		# alternatively we can just get the eigenvalues instead
		val, vec = np.linalg.eig(Zp)

		# the magnitude gives the damping and the angle gives the frequency
		damping_coeffs = np.zeros(comp)    # corresponds to width_L
		frequency_coeffs = np.zeros(comp)  # corresponds to ppm
		for i in range(comp):
			damping_coeffs[i]   = -np.log(abs(val[i])) / (1/dat.fs) / np.pi
			frequency_coeffs[i] = (np.angle(val[i]) / ((1/dat.fs) * 2 * np.pi))
		
		# we can calculate the magnitude of each signal from the
		# RHS decomposition, linalg.inv(vec) * (S_K * V_K.H)[:, 0] but
		# a simpler but more expensive way is to construct a basis set from the
		# known damping and frequency components and fit to the original data to
		# get the amplitudes and phase data
		X = np.zeros((dat.n, comp), "complex")
		# TODO this should use the singlet fitting module to make the basis
		for i in range(comp):
			X[:, i] = self.lorentzian(dat.t,
									  frequency_coeffs[i],
									  0,
									  damping_coeffs[i]) * dat.n

		# we use the linear non-iterative least squares again
		U2, s2, V2 = np.linalg.svd(np.matrix(X), full_matrices=False)
		s2_inv = np.diag(1 / s2)
		beta = V2.H * s2_inv * U2.H * np.matrix(np.reshape(dat.signal[0:dat.n], (dat.n, 1)))
		amplitudes = np.squeeze(np.array(np.abs(beta)))
		phases     = np.squeeze(np.rad2deg(np.array(np.angle(beta))))
		
		amplitudes = [x for _,x in sorted(zip(frequency_coeffs,amplitudes))]
		damping_coeffs = [x for _,x in sorted(zip(frequency_coeffs,damping_coeffs))]
		phases = [x for _,x in sorted(zip(frequency_coeffs,phases))]
		frequency_coeffs = sorted(frequency_coeffs)
		
		console.append('peak \t feq \t ampl \t damp \t phase')
		for i in range(0,comp):
			console.append(str(i) + '\t' + str(frequency_coeffs[i]) + '\t' + str(amplitudes[i]) + '\t' + str(damping_coeffs[i]) + '\t' + str(phases[i]))
		
		return np.array(list(range(0, comp))), np.array(damping_coeffs), np.array(frequency_coeffs), np.array(amplitudes), np.array(phases)

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())