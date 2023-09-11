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

from apps_processing import AppsProcessing

qtCreatorFile = "apps/ui/apps.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

# ---- Main Application Class ---- #
class MyApp(QtWidgets.QWidget, Ui_MainWindow, AppsProcessing):

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

		post_processing = "quecc" if self.queccRadioButton_bruker.isChecked() \
		else "quality" if self.qualityRadioButton_bruker.isChecked() \
		else "ecc" if self.eccRadioButton_bruker.isChecked() \
		else None

		sup_file, uns_file = MyApp.run_bruker_conversion(out_name,
										sup_file_path=str(self.inputFilenameInput_bruker.text()),
										ref_file_path=str(self.referenceFilenameInput_bruker.text()),
										baseline_correction=self.baselineCorrectionButton1_bruker.isChecked(),
										post_processing=post_processing,
										quality_points=self.qualityPointsInput_bruker.text()
										)
		self.conversionConsole_bruker.append('sup_file: ' + sup_file.file_dir)
		sup_file.print_params(self.conversionConsole_bruker)
		self.conversionConsole_bruker.append('')
		self.scaleFactorInput_bruker.setText(str(sup_file.ConvS))
		self.timeDelayInput_bruker.setText(str(float(sup_file.digShift) * 1/sup_file.fs))

		self.conversionConsole_bruker.append('uns_file: ' + uns_file.file_dir)
		uns_file.print_params(self.conversionConsole_bruker)
		self.conversionConsole_bruker.append('')
		self.canvas[2].draw()

	# ---- Water Removal ---- #
	def runWaterRemoval(self):
		'''
			This method performs removal of the residual water signal 
			using Hankel Singular Value decomposition.
		'''

		self.dat_hsvd, self.dat_wr = MyApp.run_water_removal(
			input_dat=self.dat,
			hsvd_points=int(self.hsvdPointsLineEdit_wr.text()),
			hsvd_ratio=float(self.hsvdRatioLineEdit_wr.text()),
			hsvd_components=int(self.hsvdComponentsLineEdit_wr.text()),
			frequency_range_xmin=float(self.XminLineEdit.text()),
			frequency_range_xmax=float(self.XmaxLineEdit.text()),
			console=self.hsvdFitConsole_wr
		)
		self.canvas[0].draw()

		self.saveWaterRemovalButton.setEnabled(True)

	def saveWaterRemoval(self):
		'''
			This method saves the water removed signal as a *.dat file.
		'''
		return MyApp.save_water_removal(self.dat, self.dat_wr)

	# ---- Macromolecule Removal ---- #
	def runMMRemoval(self):
		'''
			This method performs macromolecule removing using Hankel
			singular value decomposition.
		'''

		# 1. Fit macromolecule spectrum with HSVD.
		peak, width_L, ppm, area, phase = MyApp.hsvd(
			self.MMDat_mmr,
			int(self.hsvdPointsLineEdit_mmr.text()),
			float(self.hsvdRatioLineEdit_mmr.text()),
			int(self.hsvdComponentsLineEdit_mmr.text()),
			console=self.hsvdFitConsole_mmr
		)
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

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())