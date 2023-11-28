from __future__ import print_function

# ---- System Libraries ---- #
from builtins import str
from builtins import range
import sys
import os
import copy

# Check for PyQt6 (for native mac M1 compatibility), otherwise continue using PyQt5
import importlib.util
PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore, QtGui, QtWidgets, uic
else:
	from PyQt5 import QtCore, QtGui, QtWidgets, uic

# ---- Math Libraries ---- #
import scipy as sp
import numpy as np
from pyfftw.interfaces import scipy_fftpack as fftw

# ---- Plotting Libraries ---- #
import matplotlib as mpl;
mpl.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (
	FigureCanvasQTAgg as FigureCanvas,
	NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt

# ---- Data Classes ---- #
from magiqdataclasses import *

qtCreatorFile = "spices/ui/SPICeS.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

# ---- Main Application Class ---- #
class MyApp(QtWidgets.QWidget, Ui_MainWindow):

	# ---- Methods to Set Up UI ---- #
	def __init__(self):
		QtWidgets.QWidget.__init__(self)
		Ui_MainWindow.__init__(self)
		self.setupUi(self)

		# Bind Buttons
		self.setDirectoryButton.clicked.connect(self.setWorkingDirectory)
		self.setDirectoryButton.setEnabled(True)

		self.outFileBrowseButton.clicked.connect(self.chooseOutputFile)
		self.outFileBrowseButton.setEnabled(False)

		self.outFileConfirmButton.clicked.connect(self.loadOutputFile)
		self.outFileConfirmButton.setEnabled(False)

		self.datFileBrowseButton.clicked.connect(self.chooseDatFile)
		self.datFileBrowseButton.setEnabled(False)

		self.datFileConfirmButton.clicked.connect(self.loadDatFile)
		self.datFileConfirmButton.setEnabled(False)

		self.loadGroupsButton.clicked.connect(self.loadVisGroupFile)
		self.loadGroupsButton.setEnabled(False)

		self.saveGroupsButton.clicked.connect(self.saveVisGroupFile)
		self.saveGroupsButton.setEnabled(False)

		self.plotButton.clicked.connect(self.plot)
		self.plotButton.setEnabled(False)

		# Setup for plotting in the "Plot Metabolites" and "Generate Guess File" Tabs
		self.canvas = [None]*1
		self.toolbar = [None]*1
		self.setPlot('Plot')

	def setPlot(self, tab):
		if tab == 'Plot':
			fig = plt.figure(1)
			self.addmpl(0, fig, self.plotFreq_mplvl)

	# ---- Methods to Load Files ---- #
	def setWorkingDirectory(self):
		self.workingDirectory = str(QtWidgets.QFileDialog.getExistingDirectory(self, 'Set Working Directory', os.path.expanduser('~')))
		if self.workingDirectory == '':
			self.workingDirectory = os.path.expanduser('~')

		self.outFileBrowseButton.setEnabled(True)
		self.datFileBrowseButton.setEnabled(True)

	def chooseOutputFile(self):
		self.outFileConfirmButton.setEnabled(False)
		prev = str(self.outFileInput.text())
		self.outFileInput.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Output File', self.workingDirectory, 'Output files (*.out)')[0]))
		self.outFileConfirmButton.setEnabled(True)
		if str(self.outFileInput.text()) == '':
			self.outFileInput.setText(str(prev))
			if str(prev) == '':
				self.outFileConfirmButton.setEnabled(False)

	def chooseDatFile(self):
		self.datFileConfirmButton.setEnabled(False)
		prev = str(self.datFileInput.text())
		self.datFileInput.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Dat File', self.workingDirectory, 'Dat files (*.dat)')[0]))
		self.datFileConfirmButton.setEnabled(True)
		if str(self.datFileInput.text()) == '':
			self.datFileInput.setText(str(prev))
			if str(prev) == '':
				self.datFileConfirmButton.setEnabled(False)

	def loadOutputFile(self):
		try:
			# load output file
			fit_file     = str(self.outFileInput.text())
			self.fit_out = OutputFile(fit_file)
			self.outFileInfoLabel.setText(fit_file.split('/')[-1] + "\nsuccessfully loaded.")

			# populate available metabolite list
			self.metabTextBrowser.clear()
			for metabolite in self.fit_out.metabolites_list:
				self.metabTextBrowser.append(str(metabolite) + ',')

			# reset buttons
			self.outFileBrowseButton.setEnabled(True)
			self.outFileConfirmButton.setEnabled(False)

		except Exception as e:
			self.outFileInfoLabel.setText("ERROR: " + str(e) + "\n>> Please try again.")

	def loadDatFile(self):
		try:
			# load dat file
			invivo_file     = str(self.datFileInput.text())
			self.invivo_dat = DatFile(invivo_file)
			self.datFileInfoLabel.setText(invivo_file.split('/')[-1] + "\nsuccessfully loaded.")

			# populate plotting parameters
			self.b0Input.setText(str(self.invivo_dat.b0))
			self.dwellTimeInput.setText(str(1/self.invivo_dat.fs))
			self.nInput.setText(str(self.invivo_dat.n))

			# reset buttons
			self.datFileBrowseButton.setEnabled(True)
			self.datFileConfirmButton.setEnabled(False)

			self.loadGroupsButton.setEnabled(True)
			self.saveGroupsButton.setEnabled(True)
			self.plotButton.setEnabled(True)

		except Exception as e:
			self.datFileInfoLabel.setText("ERROR: " + str(e) + "\n>> Please try again.")

	def loadVisGroupFile(self):
		try:
			vis_filename = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Visualization Information File', os.getcwd() + '/spices/vis', 'Visualization information files (*.vis)')[0])

			if vis_filename != '':
				print('Loading vis groups from: ' + vis_filename)
				vis_file = open(vis_filename, 'r')

				self.groupsTextEdit.clear()

				for line in vis_file:
					self.groupsTextEdit.appendPlainText(line.replace('\n',''))

				vis_file.close()

		except Exception as e:
			print(e)

	def saveVisGroupFile(self):
		try:
			vis_filename = str(QtWidgets.QFileDialog.getSaveFileName(self, 'Save Visualization Information File', os.getcwd() + '/spices/vis', 'Visualization information files (*.vis)')[0])

			if vis_filename != '':
				print('Writing vis groups to: ' + vis_filename)
				vis_file = open(vis_filename, 'w')
				vis_file.write(self.groupsTextEdit.self.shiftsTextEdit.toPlainText())
				vis_file.close()

		except Exception as e:
			print(e)

	# ---- Methods to Calculate and Plot Spectra ---- #
	def addmpl(self, canvas_index, fig, vertical_layout):
		self.canvas[canvas_index] = FigureCanvas(fig)
		vertical_layout.addWidget(self.canvas[canvas_index])
		self.canvas[canvas_index].draw()

		self.toolbar[canvas_index] = NavigationToolbar(self.canvas[canvas_index], self, coordinates=True)
		vertical_layout.addWidget(self.toolbar[canvas_index])

	def plot(self):
		try:
			self.plotErrorLabel.setText('')

			plot_groups = []

			# Read Groups
			defined_groups = str(self.groupsTextEdit.toPlainText()).replace('\n','').replace(' ','').split(';')
			defined_groups = defined_groups[0:-1]

			# Parse each group as if defining a constraint group (recylcing code from PINTS)
			for group in defined_groups:
				group_info = group.split(',')
				group_name = group_info[0]
				group_members = group_info[1:]

				new_group = CSTGroup('', group_name, group_members, 0, 0)

				plot_groups.append(new_group)


			# Gather plotting variables
			sfactor = float(self.sfactorInput.text())
			pfactor = np.deg2rad(float(self.pfactorInput.text()))
			VSHIFT  = float(self.VSHIFTInput.text())
			HSHIFT  = float(self.HSHIFTInput.text())
			FT1     = int(self.ft1Input.text())
			lb      = float(self.lbInput.text())

			b0 = float(self.b0Input.text())
			dwell_time = float(self.dwellTimeInput.text())
			acq_time = int(self.nInput.text())*dwell_time
			fs = 1/dwell_time
			t = np.arange(0, acq_time, dwell_time)

			tableau10 = [(31, 119, 180), (255, 127, 14),    
						 (44, 160, 44), (214, 39, 40),    
						 (148, 103, 189), (140, 86, 75),    
						 (227, 119, 194), (127, 127, 127),    
						 (188, 189, 34), (23, 190, 207)]

			for i in range(len(tableau10)):    
				r, g, b = tableau10[i]    
				tableau10[i] = (r / 255., g / 255., b / 255.)

			# Calculate Summed Spectra
			fit_spec_sum = []; fit_fid_sum = []
			for metabolite in self.fit_out.metabolites_list:
				fid  = self.fit_out.metabolites[metabolite].getFID(0, b0, t, 0, 1, pfactor, 0, lb)
				if not(self.extrap0CheckBox.isChecked()):
					fid  = fid[FT1:]
				spec = fftw.fftshift(fftw.fft(fid))

				n = np.size(fid)
				f = np.arange(+n//2,-n//2,-1)*(fs/n)*(1/b0)
				fit_f = -f

				# fit_f, spec = self.fit_out.metabolites[metabolite].getSpec(0, b0, t, 0, 1, pfactor, 0, lb, fs)
				fit_fid_sum.append(fid)
				fit_spec_sum.append(spec)
			fit_fid_sum  =         np.sum(np.array(fit_fid_sum), axis=0)
			fit_spec_sum = np.real(np.sum(np.array(fit_spec_sum), axis=0))

			# Calculate In-Vivo Spectra
			invivo_dat_temp        = copy.copy(self.invivo_dat);
			invivo_dat_temp.signal = invivo_dat_temp.signal[0:np.size(t)] * np.exp(1j*pfactor) * np.exp(-sp.pi*lb*t)
			if not(self.extrap0CheckBox.isChecked()):
				invivo_dat_temp.signal = invivo_dat_temp.signal[FT1:]
			else:
				invivo_dat_temp.signal = np.hstack((fit_fid_sum[0:FT1], invivo_dat_temp.signal[FT1:]))
			invivo_f, invivo_spec  = invivo_dat_temp.getSpec()
			invivo_spec = np.real(invivo_spec)

			# Calculate Fitted Spectra
			fit_spec = []; fit_spec_names = []
			for group in plot_groups:
				group_spec = []
				for member in group.members:
					fid  = self.fit_out.metabolites[member].getFID(0, b0, t, 0, 1, pfactor, 0, lb)
					if not(self.extrap0CheckBox.isChecked()):
						fid  = fid[FT1:]
					spec = fftw.fftshift(fftw.fft(fid))

					n = np.size(fid)
					f = np.arange(+n//2,-n//2,-1)*(fs/n)*(1/b0)
					fit_f = -f
					# fit_f, spec = self.fit_out.metabolites[member].getSpec(0, b0, t, 0, 1, pfactor, 0, lb, fs)
					group_spec.append(spec)
				group_spec = np.real(np.sum(np.array(group_spec), axis=0))
				fit_spec.append(group_spec)
				fit_spec_names.append(group.name)

			plt.figure(1)
			plt.clf()
			
			ax = plt.subplot(111)
			ax.spines['top'].set_visible(False)    
			ax.spines['bottom'].set_visible(True)    
			ax.spines['right'].set_visible(False)    
			ax.spines['left'].set_visible(False)
			ax.get_xaxis().tick_bottom()
			ax.get_yaxis().set_visible(False)
			plt.xlim(5,0)

			for i, spec in enumerate(fit_spec):
				plt.plot(np.array(-fit_f)+sfactor, np.real(fit_spec[-(i+1)])-i*VSHIFT*np.amax(np.real(fit_spec_sum)), color=tableau10[-((i)%np.size(tableau10, axis=0)+1)], lw=1.5, alpha=0.8)
				ax.text(5, -i*VSHIFT*np.amax(np.real(fit_spec_sum)), fit_spec_names[-(i+1)])

			print(np.size(fit_f), np.size(fit_spec_sum), np.size(invivo_spec))

			plt.plot(np.array(-fit_f)+sfactor, (np.real(fit_spec_sum)-np.roll(np.real(invivo_spec)[0:np.size(fit_f)], int(HSHIFT)))+np.amax(np.real(invivo_spec))+3.0*VSHIFT*np.amax(np.real(invivo_spec)), color=tableau10[2], lw=1.5, alpha=0.8)
			plt.plot(np.array(-fit_f)+sfactor, np.roll(np.real(invivo_spec)[0:np.size(fit_f)], int(HSHIFT))+1.5*VSHIFT*np.amax(np.real(invivo_spec)), color=tableau10[0], lw=1.5)
			plt.plot(np.array(-fit_f)+sfactor, np.real(fit_spec_sum)+1.5*VSHIFT*np.amax(np.real(invivo_spec)), color=tableau10[1], lw=1.5, alpha=0.8)

			plt.xlabel('ppm')

			self.canvas[0].draw()

		except Exception as e:
			self.plotErrorLabel.setText('>> ERROR: ' + str(e) + '\n>> Could not plot. Please try again.')

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())