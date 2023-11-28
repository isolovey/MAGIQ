from __future__ import print_function

# ---- System Libraries ---- #
from builtins import str
from builtins import range
import sys
import os
import subprocess
import traceback

# Check for PyQt6 (for native mac M1 compatibility), otherwise continue using PyQt5
import importlib.util
PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore, QtGui, QtWidgets, uic
else:
	from PyQt5 import QtCore, QtGui, QtWidgets, uic

from collections import defaultdict
from itertools import groupby

# ---- Math Libraries ---- #
import scipy as sp
import numpy as np

# ---- Plotting Libraries ---- #
import matplotlib as mpl;
mpl.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (
	FigureCanvasQTAgg as FigureCanvas,
	NavigationToolbar2QT as NavigationToolbar)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec

# ---- Image Libraries ---- #
import nibabel as nib
from autozoom import *

# ---- Data Classes ---- #
from magiqdataclasses import *

# ---- Spectroscopy Calculations ---- #
import metabcalcs as mc

qtCreatorFile = "barstool/ui/BARSTOOL.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QtWidgets.QWidget, Ui_MainWindow):
	def __init__(self):
		'''
			This method initializes the UI and binds methods to UI buttons.
		'''
		QtWidgets.QWidget.__init__(self)
		Ui_MainWindow.__init__(self)
		self.setupUi(self)

		# Bind buttons to methods in each tab
		self.setBindings('Sum Amplitudes')
		self.setBindings('Brain Extraction')
		self.setBindings('Brain Segmentation')
		self.setBindings('Set Parameters')
		self.setBindings('Quantify Metabolites')

		# If running on Windows, set the WSL environment up so FSL can be used
		if os.name == 'nt':
			print('Running on Windows ... setting up WSL environment.')
			proc = subprocess.Popen(["wsl", "bash", "-c", "grep -oE 'gcc version ([0-9]+)' /proc/version"], stdout=subprocess.PIPE, shell=True)
			(out, err) = proc.communicate()
						
			if int(out.decode().split(' ')[-1]) > 5:
				print('WSL2 detected.')
				proc = subprocess.Popen(["wsl", "echo", "$(cat /etc/resolv.conf | grep nameserver)"], stdout=subprocess.PIPE, shell=True)
				(out, err) = proc.communicate()
				os.environ["DISPLAY"] = out.decode().split(' ')[-1].rstrip() + ":0"
			else:
				print('WSL1 detected.')
				os.environ["DISPLAY"] = ":0"
			os.environ["FSLDIR"] = "/usr/local/fsl"
			os.environ["FSLOUTPUTTYPE"] = "NIFTI_GZ"
			os.environ["WSLENV"] = "FSLDIR/u:FSLOUTPUTTYPE/u:DISPLAY/u"

			print("DISPLAY", os.environ["DISPLAY"])
			print("FSLDIR", os.environ["FSLDIR"])
			print("FSLOUTPUTTYPE", os.environ["FSLOUTPUTTYPE"])
			print("WSLENV", os.environ["WSLENV"])

	def tree(self): return defaultdict(self.tree)

	def setBindings(self, tab):
		'''
			This method binds methods to the UI buttons.
		'''
		if tab == 'Sum Amplitudes':
			
			self.setWorkingDirectoryButton.clicked.connect(self.setWorkingDirectory)
			self.setWorkingDirectoryButton.setEnabled(True)

			self.loadOutputsButton.clicked.connect(self.loadOutputs)
			self.loadOutputsButton.setEnabled(False)

			self.confirmIDsButton.clicked.connect(self.confirmIDs)
			self.confirmIDsButton.setEnabled(False)

			self.confirmSaveFileButton.clicked.connect(self.confirmSaveFileName)
			self.confirmSaveFileButton.setEnabled(False)

			self.calculateButton.clicked.connect(self.calculate)
			self.calculateButton.setEnabled(False)

		elif tab == 'Brain Extraction':

			self.BETIMAGELIST_INDEX = 0

			self.selectBetImagesButton.clicked.connect(self.loadBetImages)
			self.selectBetImagesButton.setEnabled(False)

			self.reorientBetImagesButton.clicked.connect(self.reorientBetImages)
			self.reorientBetImagesButton.setEnabled(False)

			self.identifyIsoButton.clicked.connect(self.identifyIso)
			self.identifyIsoButton.setEnabled(False)

			self.confirmIsoButton.clicked.connect(self.confirmIso)
			self.confirmIsoButton.setEnabled(False)

			self.confirmFthreshListButton.clicked.connect(self.confirmFthreshList)
			self.confirmFthreshListButton.setEnabled(False)

			self.runBetTestButton.clicked.connect(self.runBetTest)
			self.runBetTestButton.setEnabled(False)

			self.checkFthreshButton.clicked.connect(self.checkFthresh)
			self.checkFthreshButton.setEnabled(False)

			self.selectFthreshButton.clicked.connect(self.selectFthresh)
			self.selectFthreshButton.setEnabled(False)

			self.runBetButton.clicked.connect(self.runBet)
			self.runBetButton.setEnabled(False)

		elif tab == 'Brain Segmentation':

			# Default values for the FAST segmentation.
			self.FAST_T = 1
			self.FAST_N = 3
			self.FAST_H = 0.1
			self.FAST_I = 4
			self.FAST_L = 20.0

			self.FAST_COMMAND = []

			self.selectFastImagesButton.clicked.connect(self.loadFastImages)
			self.selectFastImagesButton.setEnabled(False)

			self.confirmImageTypeButton.clicked.connect(self.confirmImageType)
			self.confirmImageTypeButton.setEnabled(False)

			self.confirmFastSettingsButton.clicked.connect(self.confirmFastSettings)
			self.confirmFastSettingsButton.setEnabled(False)

			self.runFastButton.clicked.connect(self.runFast)
			self.runFastButton.setEnabled(False)

		elif tab == 'Set Parameters':

			self.loadMetabParamsButton.clicked.connect(self.loadMetabParams)
			
			self.saveMetabParamsButton.clicked.connect(self.saveMetabParams)
			self.saveMetabParamsButton.setEnabled(False)

			self.confirmParamsButton.clicked.connect(self.verifyParams)

		elif tab == 'Quantify Metabolites':

			self.loadOutputsButton_quant.clicked.connect(self.loadQuantFiles)
			self.loadOutputsButton_quant.setEnabled(False)

			self.confirmSaveFileButton_quant.clicked.connect(self.setQuantSaveFile)
			self.confirmSaveFileButton_quant.setEnabled(False)

			self.runQuantButton.clicked.connect(self.runQuant)
			self.runQuantButton.setEnabled(False)

			# Set Up Plotting Region
			fig = plt.figure(1)
			self.canvas = FigureCanvas(fig)
			self.plotQuant_mplvl.addWidget(self.canvas)
			self.canvas.draw()

	# ---- Methods for 'Sum Amplitudes' Tab ---- #
	def setWorkingDirectory(self):
		'''
			This method sets the default working directory for this instance of the application.
		'''
		self.workingDirectory = str(QtWidgets.QFileDialog.getExistingDirectory(self, 'Set Working Directory', os.path.expanduser('~')))
		if self.workingDirectory == '':
			self.workingDirectory = os.path.expanduser('~')

		self.consoleOutputText.append('Working directory set to:')
		self.consoleOutputText.append(' >> ' + str(self.workingDirectory))
		self.consoleOutputText.append('')

		self.setWorkingDirectoryButton.setEnabled(False)
		self.loadOutputsButton.setEnabled(True)
		self.selectBetImagesButton.setEnabled(True)
		self.selectFastImagesButton.setEnabled(True)
		self.loadOutputsButton_quant.setEnabled(True)

	def loadOutputs(self):
		'''
			This method presents a dialog allowing the user to select which *.out files to analyze.
		'''
		self.consoleOutputText.append('===== CALCULATE AMPLITUDES AND CRLBS =====')

		self.fileList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open Output File', self.workingDirectory, 'fitMAN Suppressed Output Files (*_sup.out)')[0]
		self.outputs = []

		if len(self.fileList) > 0:
			
			self.consoleOutputText.append('The following output files were loaded:')
			for (i, file) in enumerate(self.fileList):
				self.fileList[i] = str(file)
				self.consoleOutputText.append(' >> ' + str(file))
				self.studyIDsTextEdit.appendPlainText(str(file).replace('.out','').split('/')[-1])
			self.consoleOutputText.append('')

			self.loadOutputsButton.setEnabled(False)
			self.confirmIDsButton.setEnabled(True)

		else:
			
			self.consoleOutputText.append('No output files selected ... try again.')
			self.consoleOutputText.append('')
			self.loadOutputsButton.setEnabled(True)

	def confirmIDs(self):
		'''
			This method checks the number of IDs inputted by the user
			against the number of *.out files loaded into the application.
		'''

		self.IDsList = self.studyIDsTextEdit.toPlainText().split('\n')
		for (i, ID) in enumerate(self.IDsList):
			self.IDsList[i] = str(ID)

		if len(self.IDsList) != len(self.fileList):
			self.consoleOutputText.append('Number of study IDs do not equal number of files ... specify the IDs again.')
			self.consoleOutputText.append('')
		else:
			self.consoleOutputText.append('The following IDs were selected:')
			self.consoleOutputText.append(' >> ' + str(self.IDsList))
			self.consoleOutputText.append('')
			self.confirmIDsButton.setEnabled(False)
			self.confirmSaveFileButton.setEnabled(True)
			self.saveFileLineEdit.setText(self.workingDirectory + '/' + '_____.csv')

	def confirmSaveFileName(self):
		'''
			This method sets the *.csv filename to which the amplitudes/crlb information will be saved.
		'''
		self.saveFileName = self.saveFileLineEdit.text()
		self.calculateButton.setEnabled(True)
		self.consoleOutputText.append('Calculations will be saved to: ' + str(self.saveFileName))
		self.consoleOutputText.append('')

	def calculate(self):
		'''
			This method calculates the amplitudes and CRLBs from the the specified *.out files and outputs
			them to the specified *.csv file.
		'''
		saveFile = open(self.saveFileName, 'w')
		for (i, file) in enumerate(self.fileList):

			# Load output file
			self.outputs.append(OutputFile(file))

			# Get metabolite list from output file
			ID = self.IDsList[i]
			metabs = self.outputs[i].metabolites_list

			# For each metabolite calculate sum amps
			amps = []
			for metab in metabs: amps.append(self.outputs[i].metabolites[metab].sumAmp())

			# For each metabolite calculate crlbs
			crlbs = []
			for metab in metabs:
				# nanmean -- fix bug where a single nan crlb in a metabolite 
				# due to 0 amplitude in water suppression region results in nan output
				# -> except the real "fix" needs to be ignoring the peaks 
				# -> OR outputing CRLBs as part of quantification (where you define which peaks to include)
				crlb_result = np.nanmean(self.outputs[i].metabolites[metab].crlb) 
				print(metab, np.nanmean(self.outputs[i].metabolites[metab].crlb))
				crlbs.append(crlb_result)

			# Write heading to file if we're at the first line
			if i == 0:
				saveFile.write('ID,')
				print('ID', end=' ')
				for metab in metabs: 
					saveFile.write(str(self.outputs[i].metabolites[metab].name)+',')
					print(self.outputs[i].metabolites[metab].name, end=' ')
				for metab in metabs: 
					saveFile.write(str(self.outputs[i].metabolites[metab].name)+',')
					print(self.outputs[i].metabolites[metab].name, end=' ')
				saveFile.write('\n')
				print('')

			# Write ID and amplitudes to file
			saveFile.write(ID+',')
			saveFile.write(str(amps).replace(' ', '').replace('[', '').replace(']',''))
			saveFile.write(',')
			saveFile.write(str(crlbs).replace(' ', '').replace('[', '').replace(']',''))
			print(ID, str(amps).replace(' ', '').replace('[', '').replace(']',''))
			print(ID, str(crlbs).replace(' ', '').replace('[', '').replace(']',''))
			saveFile.write('\n')
			print('')

		saveFile.close()

		# Reset interface
		self.setWorkingDirectoryButton.setEnabled(True)
		self.loadOutputsButton.setEnabled(True)
		self.confirmIDsButton.setEnabled(False)
		self.confirmSaveFileButton.setEnabled(False)
		self.calculateButton.setEnabled(False)

		self.studyIDsTextEdit.clear()
		self.saveFileLineEdit.clear()

	# ---- Methods for 'Brain Extraction' Tab ---- #
	def loadBetImages(self):
		'''
			This method presents a dialog that allows the user to select *.nii.gz files
			which will be processed by BET for brain extraction.
		'''
		self.consoleOutputText.append('===== BRAIN EXTRACTION =====')

		self.betImageList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open Image File', self.workingDirectory, 'Image Files (*.nii.gz)')[0]

		if len(self.betImageList) > 0:
			
			self.consoleOutputText.append('The following image files were loaded:')
			for (i, image) in enumerate(self.betImageList):
				self.betImageList[i] = str(image)
				self.consoleOutputText.append(' >> ' + str(image))
			self.consoleOutputText.append('')

			self.selectBetImagesButton.setEnabled(False)
			self.reorientBetImagesButton.setEnabled(True)

		else:
			
			self.consoleOutputText.append('No image files selected ... try again.')
			self.consoleOutputText.append('')

			self.selectBetImagesButton.setEnabled(True)

	def reorientBetImages(self):
		'''
			This method applies fslreorient2std for all specfied images.
			fslreorient2std requries that file paths don't have a space in it.
		'''
		self.consoleOutputText.append('Reorienting ...')
		for (i, image) in enumerate(self.betImageList):
			if os.name == 'nt':
				image_win = os.popen('wsl wslpath "' + str(image) + '"').read().replace('\n','')
				command = 'wsl ' + '$FSLDIR/bin/fslreorient2std "' + str(image_win) + '" "' + str(image_win).replace('.nii.gz', '_std.nii.gz') + '"'
			else:
				command = ['fslreorient2std', str(image), str(image).replace('.nii.gz', '_std.nii.gz')]
			self.consoleOutputText.append(' >> ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
			print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
			print(subprocess.check_output(command, stderr=subprocess.STDOUT).decode())
		self.consoleOutputText.append('')

		self.identifyIsoButton.setEnabled(True)
		self.reorientBetImagesButton.setEnabled(False)

	def identifyIso(self):
		'''
			This method launches FSL eyes allowing the user to find the isocentre of the volume.
		'''
		self.consoleOutputText.append('---- Processing Image ' + str(self.BETIMAGELIST_INDEX+1) + ' of ' + str(np.size(self.betImageList)) + ' ---- ')

		image = self.betImageList[self.BETIMAGELIST_INDEX].replace('.nii.gz', '_std.nii.gz')
		
		if os.name == 'nt':
			image_win = os.popen('wsl wslpath "' + str(image) + '"').read().replace('\n','')
			command = 'wsl ' + '$FSLDIR/bin/fsleyes "' + str(image_win) + '"'
		else:
			command = ['fsleyes', str(image)]

		self.consoleOutputText.append('Finding isocenter of image ...')
		self.consoleOutputText.append(' >> ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
		print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
		subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

		self.consoleOutputText.append('')
		self.confirmIsoButton.setEnabled(True)

	def confirmIso(self):
		'''
			This method loads the specified isocentre coordinates into the application.
		'''
		try:
			self.isoX = int(self.isoXLineEdit.text())
			self.isoY = int(self.isoYLineEdit.text())
			self.isoZ = int(self.isoZLineEdit.text())

			self.consoleOutputText.append('Brain Isocenter Selected: ' + str(self.isoX) + ', ' + str(self.isoY) + ', ' + str(self.isoZ))
			self.consoleOutputText.append('')

			self.confirmFthreshListButton.setEnabled(True)
			self.identifyIsoButton.setEnabled(False)
			self.confirmIsoButton.setEnabled(False)
		
		except Exception as e:
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()

	def confirmFthreshList(self):
		'''
			This method adds the thresholds specified from the user into the combo box.
		'''
		# Get threshold list from user
		try:
			fthreshList = self.fthreshListLineEdit.text().replace(' ', '').split(',')
			self.fthreshList = []
			for fthresh in fthreshList:
				self.fthreshList.append(float(fthresh)) # add to global list
				self.selectFthreshComboBox.addItem(str(fthresh)) # add to combo box

			self.consoleOutputText.append('The following thresholds for brain extraction will be tested ...')
			self.consoleOutputText.append(' >> ' + str(self.fthreshList))
			self.consoleOutputText.append('')

			self.runBetTestButton.setEnabled(True)
			self.confirmFthreshListButton.setEnabled(False)
		except Exception as e:
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc() + '\n\n' + 'Check your intensity thresholds.')
			msg.setWindowTitle('Error')
			msg.exec_()

	def runBetTest(self):
		'''
			This method runs BET using the thresholds specified by the user.
			BET requires that filenames have no spaces in them.
		'''
		try:
			image = self.betImageList[self.BETIMAGELIST_INDEX].replace('.nii.gz', '')
			brainctr = str(self.isoX) + ' ' + str(self.isoY) + ' ' + str(self.isoZ)
			luts = ["Red-Yellow", "Blue-Lightblue", "Red", "Blue", "Green", "Yellow", "Pink", "Hot", "Cool", "Copper"]

			self.consoleOutputText.append('Running BET with at most ' + str(np.size(luts)) + ' test thresholds on the image ...')
			for fthresh in self.fthreshList[0:np.size(luts)]:
				if os.name == 'nt':
					image_win = os.popen('wsl wslpath "' + str(image) + '"').read().replace('\n','')
					command = ['wsl', '$FSLDIR/bin/bet', '"'+str(image_win)+'_std"']
					command.append('"'+str(image_win)+'_std_brain_f'+str(fthresh)+'"')
				else:
					command = ['bet', '"'+str(image)+'_std"']
					command.append('"'+str(image)+'_std_brain_f'+str(fthresh)+'"')
				command.append('-f'); command.append(str(fthresh))
				command.append('-g'); command.append('0')
				command.append('-c'); command.append(str(brainctr))
				self.consoleOutputText.append(' >> ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
				print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))

				if os.name == 'nt':
					print(subprocess.check_output(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''), stderr=subprocess.STDOUT).decode())
				else:
					print(subprocess.check_output(command, stderr=subprocess.STDOUT).decode())
			self.consoleOutputText.append('')

			self.runBetTestButton.setEnabled(False)
			self.checkFthreshButton.setEnabled(True)
		except Exception as e:
			traceback.print_exc()
				
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()

	def checkFthresh(self):
		'''
			This method launches fsleyes in a separate window, with the results of BET using varying thresholds preloaded.
			This allows the user to decide which threshold is best to use.
		'''
		image = self.betImageList[self.BETIMAGELIST_INDEX].replace('.nii.gz', '')
		brainctr = str(self.isoX) + ' ' + str(self.isoY) + ' ' + str(self.isoZ)
		luts = ["Red-Yellow", "Blue-Lightblue", "Red", "Blue", "Green", "Yellow", "Pink", "Hot", "Cool", "Copper"]

		self.consoleOutputText.append('Select threshold that produced the best brain extraction ...')
		if os.name == 'nt':
			image_win = os.popen('wsl wslpath "' + str(image) + '"').read().replace('\n','')
			command = ['wsl', '$FSLDIR/bin/fsleyes', '"'+str(image_win)+'_std"']
			for (i, fthresh) in enumerate(self.fthreshList[0:np.size(luts)]):
				command.append('"'+str(image_win)+'_std_brain_f'+str(fthresh)+'"')
				command.append('-cm'); command.append('"' + str(luts[i]) +'"')
				command.append('-a'); command.append('0.0')
		else:
			command = ['fsleyes', '"'+str(image)+'_std"']
			for (i, fthresh) in enumerate(self.fthreshList[0:np.size(luts)]):
				command.append('"'+str(image)+'_std_brain_f'+str(fthresh)+'"')
				command.append('-cm'); command.append('"' + str(luts[i]) +'"')
				command.append('-a'); command.append('0.0')
		
		self.consoleOutputText.append(' >> ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
		self.consoleOutputText.append('')

		print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))

		if os.name == 'nt':
			subprocess.Popen(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''), shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		else:
			subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

		self.checkFthreshButton.setEnabled(False)
		self.selectFthreshButton.setEnabled(True)

	def selectFthresh(self):
		'''
			This method loads the user selected threshold into the application..
		'''
		self.fthresh_best = self.selectFthreshComboBox.currentText()
		self.consoleOutputText.append('Selected threshold: ' + str(self.fthresh_best))
		self.consoleOutputText.append('')

		self.selectFthreshButton.setEnabled(False)
		self.runBetButton.setEnabled(True)

	def runBet(self):
		'''
			This method runs BET using the user specified isocentre and threshold.
		'''
		image = self.betImageList[self.BETIMAGELIST_INDEX].replace('.nii.gz', '')
		brainctr = str(self.isoX) + ' ' + str(self.isoY) + ' ' + str(self.isoZ)
		luts = ["Red-Yellow", "Blue-Lightblue", "Red", "Blue", "Green", "Yellow", "Pink", "Hot", "Cool", "Copper"]

		# Remove test BET files
		self.consoleOutputText.append('Removing test BET files ...')
		for (i, fthresh) in enumerate(self.fthreshList[0:np.size(luts)]):
			try:
				print('os.remove('+str(image)+'_std_brain_f'+str(fthresh)+'.nii.gz'+')')
				os.remove(str(image)+'_std_brain_f'+str(fthresh)+'.nii.gz')
			except Exception as e:
				traceback.print_exc()
				
				msg = QtWidgets.QMessageBox()
				msg.setIcon(QtWidgets.QMessageBox.Critical)
				msg.setText('Error ' + str(e) + 'has occurred!')
				msg.setInformativeText(traceback.format_exc())
				msg.setWindowTitle('Error')
				msg.exec_()
				
				print('Continuing ...')

		self.consoleOutputText.append('')

		self.consoleOutputText.append('Running BET with selected threshold (' + str(self.fthresh_best) + ')')

		image = self.betImageList[self.BETIMAGELIST_INDEX].replace('.nii.gz', '')

		# Run BET
		if os.name == 'nt':
			image_win = os.popen('wsl wslpath "' + str(image) + '"').read().replace('\n','')
			command = ['wsl', '$FSLDIR/bin/bet', '"'+str(image_win)+'_std"']
			command.append('"'+str(image_win)+'_std_brain_f'+str(self.fthresh_best)+'"')
		else:
			command = ['bet', '"'+str(image)+'_std"']
			command.append('"'+str(image)+'_std_brain_f'+str(self.fthresh_best)+'"')
		command.append('-f'); command.append(str(self.fthresh_best))
		command.append('-g'); command.append('0')
		command.append('-c'); command.append(str(brainctr))
		command.append('-m')
		self.consoleOutputText.append(' >> ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
		print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))

		if os.name == 'nt':
			print(subprocess.check_output(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''), stderr=subprocess.STDOUT).decode())
		else:
			print(subprocess.check_output(command, stderr=subprocess.STDOUT).decode())

		# Rename files
		try:
			print('os.rename('+str(image)+'_std_brain_f'+str(self.fthresh_best)+'.nii.gz' +','+ str(image)+'_std_brain.nii.gz' +')')
			os.rename(str(image)+'_std_brain_f'+str(self.fthresh_best)+'.nii.gz', str(image)+'_std_brain.nii.gz')
		except Exception as e:
			traceback.print_exc()
				
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()
				
			print('Continuing ...')

		try:
			print('os.rename('+str(image)+'_std_brain_f'+str(self.fthresh_best)+'_mask.nii.gz' +','+ str(image)+'_std_brain_mask.nii.gz'+')')
			os.rename(str(image)+'_std_brain_f'+str(self.fthresh_best)+'_mask.nii.gz', str(image)+'_std_brain_mask.nii.gz')
		except Exception as e:
			traceback.print_exc()
				
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc())
			msg.setWindowTitle('Error')
			msg.exec_()
				
			print('Continuing ...')

		self.consoleOutputText.append('---- Processed Image ' + str(self.BETIMAGELIST_INDEX+1) + ' of ' + str(np.size(self.betImageList)) + ' ---- ') 
		self.consoleOutputText.append('')

		if self.BETIMAGELIST_INDEX < np.size(self.betImageList) - 1:
			# Move on to next image
			self.BETIMAGELIST_INDEX = self.BETIMAGELIST_INDEX + 1
			
			# Reset buttons
			self.selectBetImagesButton.setEnabled(False)
			self.reorientBetImagesButton.setEnabled(False)
			self.identifyIsoButton.setEnabled(True)
			self.confirmIsoButton.setEnabled(False)
			self.confirmFthreshListButton.setEnabled(False)
			self.checkFthreshButton.setEnabled(False)
			self.selectFthreshButton.setEnabled(False)
			self.runBetButton.setEnabled(False)
			
			# Reset fthresh choices
			self.selectFthreshComboBox.clear()
		else:
			# Reset buttons
			self.selectBetImagesButton.setEnabled(True)
			self.reorientBetImagesButton.setEnabled(False)
			self.identifyIsoButton.setEnabled(False)
			self.confirmIsoButton.setEnabled(False)
			self.confirmFthreshListButton.setEnabled(False)
			self.checkFthreshButton.setEnabled(False)
			self.selectFthreshButton.setEnabled(False)
			self.runBetButton.setEnabled(False)
			
			# Reset fthresh choices
			self.selectFthreshComboBox.clear()

	# ---- Methods for Brain Segmentation Tab ---- #
	def loadFastImages(self):
		'''
			This method shows a dialog in which the user can load multiple images to be segmented
			into GM, WM, CSF using FSL FAST.
		'''
		self.consoleOutputText.append('===== BRAIN SEGMENTATION =====')

		self.fastImageList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open Image File', self.workingDirectory, 'Image Files (*.nii.gz)')[0]

		if len(self.fastImageList) > 0:
			
			self.consoleOutputText.append('The following image files were loaded:')
			for (i, image) in enumerate(self.fastImageList):
				self.fastImageList[i] = str(image)
				self.consoleOutputText.append(' >> ' + str(image))
			self.consoleOutputText.append('')

			self.selectFastImagesButton.setEnabled(False)
			self.confirmImageTypeButton.setEnabled(True)

		else:
			
			self.consoleOutputText.append('No image files selected ... try again.')
			self.consoleOutputText.append('')

			self.selectFastImagesButton.setEnabled(True)

	def confirmImageType(self):
		'''
			This method sets the FAST parameter for the 
			corresponding image type (t1, t2, pd-weighted).
		'''
		if self.T1RadioButton.isChecked():
			self.FAST_T = 1
			self.consoleOutputText.append('Images are T1 Weighted.')
		elif self.T2RadioButton.isChecked():
			self.FAST_T = 2
			self.consoleOutputText.append('Images are T2 Weighted.')
		elif self.PdRadioButton.isChecked():
			self.FAST_T = 3
			self.consoleOutputText.append('Images are Proton Density Weighted.')

		self.consoleOutputText.append('')

		self.confirmImageTypeButton.setEnabled(False)
		self.confirmFastSettingsButton.setEnabled(True)

	def confirmFastSettings(self):
		'''
			This method loads the FAST parameters specified by the user in the
			application's UI.
		'''
		self.FAST_H = float(self.fast_hLineEdit.text())
		self.FAST_I = float(self.fast_iLineEdit.text())
		self.FAST_L = float(self.fast_lLineEdit.text())

		self.consoleOutputText.append('The following FAST command will be run on the images ...')
		self.consoleOutputText.append('>> fast -t ' + str(self.FAST_T) + ' -n ' + str(self.FAST_N) + ' -H ' + str(self.FAST_H) + ' -I ' + str(self.FAST_I) + ' -l ' + str(self.FAST_L) + ' -o <filename> <filename>')
		self.consoleOutputText.append('')

		self.runFastButton.setEnabled(True)
		self.confirmFastSettingsButton.setEnabled(False)

	def runFast(self):
		'''
			This method runs FSL FAST on the user specified images using the 
			user specified FAST parameters.
		'''
		try:
			self.consoleOutputText.append('Running FAST on images ...')

			for image in self.fastImageList:
				image = image.replace('.nii.gz', '')
				image_win = os.popen('wsl wslpath "' + str(image) + '"').read().replace('\n','')
			
				if os.name == 'nt':
					command = ['wsl $FSLDIR/bin/fast']
				else:
					command = ['fast']
			
				command.append('-t'); command.append(str(int(self.FAST_T)))
				command.append('-n'); command.append(str(int(self.FAST_N)))
				command.append('-H'); command.append(str(self.FAST_H))
				command.append('-I'); command.append(str(int(self.FAST_I)))
				command.append('-l'); command.append(str(self.FAST_L))
				command.append('-v')
			
				if os.name == 'nt':
					command.append('-o "' + str(image_win) + '" "' + str(image_win) + '"')
				else:
					command.append('-o "' + str(image) + '" "' + str(image) + '"')
			
				self.consoleOutputText.append('>> ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
				print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))

				if os.name == 'nt':
					print(subprocess.check_output(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''), stderr=subprocess.STDOUT).decode())
				else:
					print(subprocess.check_output(command, stderr=subprocess.STDOUT).decode())

			self.consoleOutputText.append('Done.')

			self.runFastButton.setEnabled(False)
			self.selectFastImagesButton.setEnabled(True)
		except Exception as e:
			traceback.print_exc()
				
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc() + '\n\n\n' + 'Check your FAST settings and try running it again.')
			msg.setWindowTitle('Error')
			msg.exec_()

	# ---- Methods for Setting Quantification Parameters ---- #
	def loadMetabParams(self):
		'''
			This method allows the user to select the parameter files that specify the parameters
			used in the quantification of metabolites.
		'''
		prev = str(self.metabParamsFileLineEdit.text())
		self.metabParamsFileLineEdit.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Quantification Information File', os.getcwd() + '/barstool/qinfo', 'Quantification Info Files (*.qinfo)')[0]))
		self.saveMetabParamsButton.setEnabled(True)
		if str(self.metabParamsFileLineEdit.text()) == '':
			self.metabParamsFileLineEdit.setText(str(prev))
			if str(prev) == '':
				self.saveMetabParamsButton.setEnabled(False)
			else:
				self.populateMetabTable()
				self.consoleOutputText.append('===== SETTING QUANTIFICATION PARAMETERS ====')
				self.consoleOutputText.append('Quantification information loaded from: ')
				self.consoleOutputText.append('>> ' + str(self.metabParamsFileLineEdit.text()))
				self.consoleOutputText.append('')
		else:
			self.populateMetabTable()
			self.consoleOutputText.append('===== SETTING QUANTIFICATION PARAMETERS ====')
			self.consoleOutputText.append('Quantification information loaded from: ')
			self.consoleOutputText.append('>> ' + str(self.metabParamsFileLineEdit.text()))
			self.consoleOutputText.append('')

	def populateMetabTable(self):
		'''
			This method parses the parameter files and populates
			the table of parameters in the UI.
		'''
		try:
			in_file = open(str(self.metabParamsFileLineEdit.text()), 'r')

			rows = []
			for line in in_file:
				if not('#' in line):
					params = line.replace('\n','').split('\t')
					if params[0] == 'water':
						print(params)
						self.protonsLineEdit_water.setText(str(params[1]))
						self.T1GMLineEdit_water.setText(str(params[2]))
						self.T2GMLineEdit_water.setText(str(params[3]))
						self.T1WMLineEdit_water.setText(str(params[4]))
						self.T2WMLineEdit_water.setText(str(params[5]))
						self.T1CSFLineEdit_water.setText(str(params[6]))
						self.T2CSFLineEdit_water.setText(str(params[7]))
					elif params[0] == 'exp':
						print(params)
						self.TRLineEdit.setText(str(params[1]))
						self.TELineEdit.setText(str(params[2]))
						self.waterConcLineEdit.setText(str(params[3]))
						self.waterConcGMLineEdit.setText(str(params[4]))
						self.waterConcWMLineEdit.setText(str(params[5]))
						if bool(params[6]):
							self.tissueConcButton.checked = True
							self.voxelConcButton.checked = False
						else:
							self.tissueConcButton.checked = False
							self.voxelConcButton.checked = True
					else:
						print(params)
						rows.append(params)

			self.metabParamsTableWidget.setRowCount(np.size(rows,0))
			self.metabParamsTableWidget.setColumnCount(8)

			for i in range(0, np.size(rows, 0)):
				self.metabParamsTableWidget.setItem(i, 0, QtWidgets.QTableWidgetItem(rows[i][0])) # metabolite
				self.metabParamsTableWidget.setItem(i, 1, QtWidgets.QTableWidgetItem(rows[i][1])) # protons
				self.metabParamsTableWidget.setItem(i, 2, QtWidgets.QTableWidgetItem(rows[i][2])) # T1 GM
				self.metabParamsTableWidget.setItem(i, 3, QtWidgets.QTableWidgetItem(rows[i][3])) # T2 GM
				self.metabParamsTableWidget.setItem(i, 4, QtWidgets.QTableWidgetItem(rows[i][4])) # T1 WM
				self.metabParamsTableWidget.setItem(i, 5, QtWidgets.QTableWidgetItem(rows[i][5])) # T2 WM
				self.metabParamsTableWidget.setItem(i, 6, QtWidgets.QTableWidgetItem(rows[i][6])) # first peak
				self.metabParamsTableWidget.setItem(i, 7, QtWidgets.QTableWidgetItem(rows[i][7])) # last peak

			self.consoleOutputText.append('Quantification information saved to: ')
			self.consoleOutputText.append('>> ' + str(self.metabParamsFileLineEdit.text()))
			self.consoleOutputText.append('')

			in_file.close()
		except Exception as e:
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc() + '\n\n\n Consider checking your *.qinfo file for errors.')
			msg.setWindowTitle('Error')
			msg.exec_()

	def saveMetabParams(self):
		'''
			This method saves the parameter values currently extant in the UI to a 
			user specified file that can be loaded and parsed in the future.
		'''
		prev = str(self.metabParamsFileLineEdit.text())
		out_filename = str(QtWidgets.QFileDialog.getSaveFileName(self, 'Save Quantification Information File', os.getcwd() + '/barstool/qinfo', 'Quantification Info Files (*.qinfo)')[0])

		self.consoleOutputText.append('Quantification information saved to: ')
		self.consoleOutputText.append('>> ' + str(self.metabParamsFileLineEdit.text()))
		self.consoleOutputText.append('')

		if out_filename == '':
			out_filename = prev

		if out_filename != '':
			out_file = open(out_filename, 'w')
			out_file.write('#\n')
			out_file.write('# Columns:\n')
			out_file.write('#     1.  Metabolite\n')
			out_file.write('#     2.  Number of protons for quantifiable singlet or whole signal sum\n')
			out_file.write('#     3.  T1 values in GM (in sec)\n')
			out_file.write('#     4.  T2 values in GM (in msec)\n')
			out_file.write('#     5.  T1 values in WM (in sec)\n')
			out_file.write('#     6.  T2 values in WM (in msec)\n')
			out_file.write('#     7.  First Peak\n')
			out_file.write('#     8.  Last Peak\n')
			out_file.write('#\n')
			for i in range(0, self.metabParamsTableWidget.rowCount()):
				out_file.write(self.metabParamsTableWidget.item(i,0).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,1).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,2).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,3).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,4).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,5).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,6).text() + '\t')
				out_file.write(self.metabParamsTableWidget.item(i,7).text())
				out_file.write('\n')
			out_file.write('#\n')
			out_file.write('# Water:\n')
			out_file.write('#\tprotons\tT1_GM\tT2_GM\tT1_WM\tT2_WM\tT1_CSF\tT2_CSF\n')
			out_file.write('water\t' \
			+ self.protonsLineEdit_water.text() + '\t' \
			+ self.T1GMLineEdit_water.text() + '\t' \
			+ self.T2GMLineEdit_water.text() + '\t' \
			+ self.T1WMLineEdit_water.text() + '\t' \
			+ self.T2WMLineEdit_water.text() + '\t' \
			+ self.T1CSFLineEdit_water.text() + '\t' \
			+ self.T2CSFLineEdit_water.text() + '\n')
			out_file.write('#\n')
			out_file.write('# Experiment:\n')
			out_file.write('#\tTR\tTE\tConc\tConcGM\tConcWM\tConcVox\n')
			out_file.write('exp\t' \
			+ self.TRLineEdit.text() + '\t' \
			+ self.TELineEdit.text() + '\t' \
			+ self.waterConcLineEdit.text() + '\t' \
			+ self.waterConcGMLineEdit.text() + '\t' \
			+ self.waterConcWMLineEdit.text() + '\t' \
			+ str(int(self.voxelConcButton.isChecked())) + '\n')
			out_file.close()

	def verifyParams(self):
		'''
			This method prints the entered parameters to the UI console,
			allowing the user to check over the parameters carefully.
		'''
		self.consoleOutputText.append('The following parameters were entered. Please check them carefully:')
		self.consoleOutputText.append('')
		self.consoleOutputText.append('\t\tProtons\tT1 (GM) [sec]\tT2 (GM) [ms]\tT1 (WM) [ms]\tT2 (WM) [ms]\tFirst Peak\tLast Peak')
		for i in range(0, self.metabParamsTableWidget.rowCount()):
			self.consoleOutputText.append( 
				self.metabParamsTableWidget.item(i,0).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,1).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,2).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,3).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,4).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,5).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,6).text() + '\t' \
				+ self.metabParamsTableWidget.item(i,7).text())
		self.consoleOutputText.append('')
		self.consoleOutputText.append('water\t' \
			+ self.protonsLineEdit_water.text() + '\t' \
			+ self.T1GMLineEdit_water.text() + '\t' \
			+ self.T2GMLineEdit_water.text() + '\t' \
			+ self.T1WMLineEdit_water.text() + '\t' \
			+ self.T2WMLineEdit_water.text() + '\t' \
			+ self.T1CSFLineEdit_water.text() + '\t' \
			+ self.T2CSFLineEdit_water.text())
		self.consoleOutputText.append('')
		self.consoleOutputText.append('             TR [ms]: ' + self.TRLineEdit.text())
		self.consoleOutputText.append('             TE [ms]: ' + self.TELineEdit.text())
		self.consoleOutputText.append('             [water]: ' + self.waterConcLineEdit.text() + ' M')
		self.consoleOutputText.append('[water] scaling (GM): ' + self.waterConcGMLineEdit.text())
		self.consoleOutputText.append('[water] scaling (WM): ' + self.waterConcWMLineEdit.text())
		self.consoleOutputText.append('     voxel conc flag: ' + str(int(self.voxelConcButton.isChecked())))
		self.consoleOutputText.append('')

	# ---- Methods for Actual Quantification ---- #
	def loadQuantFiles(self):
		'''
			This method presents a dialog to the user allowing them to 
			select the *.out files for quantification. Note that the method automatically
			looks for associated files required for quantification based on an 
			assumed naming convention:

			NAMING CONVENTION FOR SPECTROSCOPY FILES
				- allows for multiple regions, but assumes _sup and _uns suffixes
				- ID_VISIT_<region description>_sup.rda / .out
				- ID_VISIT_<region description>_uns.rda / .out

			NAMING CONVENTION FOR IMAGES
				- assumes a whole brain T1 is acquired for each spectroscopy session
				- ID_VISIT_<file suffixes>.nii.gz
					- ID_VISIT.nii.gz = original T1
			 		- ID_VISIT_std.nii.gz = result after fslreorient2std
			 		- ID_VISIT_std_brain.nii.gz = result after fsl BET
					- ID_VISIT_std_brain_seg.nii.gz = result after fsl FAST
		'''
		self.consoleOutputText.append('===== METABOLITE QUANTIFICATION =====')

		self.quantFileList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open Output File', self.workingDirectory, 'fitMAN Suppressed Output Files (*_sup.out)')[0]

		self.quantIDs = []
		self.quantSupRDAs = []
		self.quantUnsupRDAs = []
		self.quantSupFiles = []
		self.quantUnsupFiles = []
		self.quantBrainFiles = []
		self.quantBrainSegFiles = []

		if len(self.quantFileList) > 0:

			self.consoleOutputText.append('The following output files were loaded:')
			for (i, file) in enumerate(self.quantFileList):

				ID        = str(file).replace('.out','').split('/')[-1].replace('_sup',''); self.quantIDs.append(ID)
				sup       = str(file); self.quantSupFiles.append(sup)
				unsup     = str(file).replace('_sup.out', '_uns.out'); self.quantUnsupFiles.append(unsup)
				sup_rda   = str(file).replace('.out', '.rda'); self.quantSupRDAs.append(sup_rda)
				unsup_rda = str(file).replace('_sup.out', '_uns.rda'); self.quantUnsupRDAs.append(unsup_rda)

				brain     = str(str(file).split('/')[0:-1]).replace('[','').replace(']','').replace(',','/').replace("'",'').replace('/ ','/') \
				            + '/' + ID.split('_')[0] + '_' + ID.split('_')[1] + '_std_brain.nii.gz'; self.quantBrainFiles.append(brain)
				brainseg  = str(str(file).split('/')[0:-1]).replace('[','').replace(']','').replace(',','/').replace("'",'').replace('/ ','/') \
						    + '/' + ID.split('_')[0] + '_' + ID.split('_')[1] + '_std_brain_seg.nii.gz'; self.quantBrainSegFiles.append(brainseg)

				self.consoleOutputText.append(str(i+1) + '.' + self.quantIDs[-1])
				self.consoleOutputText.append('  |' + self.quantSupRDAs[-1])
				self.consoleOutputText.append('  |' + self.quantUnsupRDAs[-1])
				self.consoleOutputText.append('  |' + self.quantSupFiles[-1])
				self.consoleOutputText.append('  |' + self.quantUnsupFiles[-1])
				self.consoleOutputText.append('  |' + self.quantBrainFiles[-1])
				self.consoleOutputText.append('  |' + self.quantBrainSegFiles[-1])

			self.loadOutputsButton_quant.setEnabled(False)
			self.saveFileLineEdit_quant.setText(self.workingDirectory + '/' + '_____.csv')
			self.confirmSaveFileButton_quant.setEnabled(True)

			self.consoleOutputText.append('')

		else:

			self.consoleOutputText.append('No output files selected ... try again.')
			self.consoleOutputText.append('')
			self.loadOutputsButton_quant.setEnabled(True)

	def setQuantSaveFile(self):
		'''
			This method sets the filename of the *.csv file where 
			the quantification results will be saved.
		'''
		self.quantSaveFileName = self.saveFileLineEdit_quant.text()
		self.runQuantButton.setEnabled(True)
		self.consoleOutputText.append('Quantification results will be saved to: ' + str(self.quantSaveFileName))
		self.consoleOutputText.append('')

	def runQuant(self):
		'''
			This method runs the quantification calculation 
			and saves the results to file.
		'''
		try:
			self.confirmSaveFileButton_quant.setEnabled(False)

			out_file = open(self.quantSaveFileName, 'w')

			# Write header
			out_file.write('ID,GM,WM,CSF,N_AVG_SUP,N_AVG_UNS,SCALE_SUP,SCALE_UNS,SCANNER,')
			for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
				out_file.write(str(self.metabParamsTableWidget.item(metab_index,0).text())+',')
			for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
				out_file.write(str(self.metabParamsTableWidget.item(metab_index,0).text())+'_CRLB,')
			out_file.write('\n')

			if self.siemensScanner.isChecked():
			
				# Loop through each ID and analyse.
				for (ID_index, ID) in enumerate(self.quantIDs):
				
					out_file.write(str(ID)+',')

					# Load sup and unsup RDAs
					sup_rda   = RDAFile(self.quantSupRDAs[ID_index], scale_fid=True)
					unsup_rda = RDAFile(self.quantUnsupRDAs[ID_index], scale_fid=True)

					# Load output files
					sup_out   = OutputFile(self.quantSupFiles[ID_index])
					unsup_out = OutputFile(self.quantUnsupFiles[ID_index])

					# Load brain image
					imgloadname = self.quantBrainFiles[ID_index]
					img = nib.nifti1.load(imgloadname)
					img_studyIDs = str(img.header['aux_file']).split(', ')  # grabs name and studydate header output from mri_convert_with_id.convert
					img_data = np.asanyarray(img.dataobj).squeeze()
				
					img_qform = img.get_qform()
					img_pvec = img_qform[:3,3]  # this is the isocenter position for the 3D nifti volume
					img_rvec = img_qform[:3,0]
					img_cvec = img_qform[:3,1]
					img_nvec = img_qform[:3,2]
					R_img = img_qform[:3,:3]

					i = np.arange(img_data.shape[0])
					j = np.arange(img_data.shape[1])
					k = np.arange(img_data.shape[2])
					i2 = np.transpose(np.tile(i,(img_data.shape[2], img_data.shape[1], 1)), (2, 1, 0)).reshape(len(i)*len(j)*len(k))
					j2 = np.transpose(np.tile(j,(img_data.shape[0], img_data.shape[2], 1)), (0, 2, 1)).reshape(len(i)*len(j)*len(k))
					k2 = np.transpose(np.tile(k,(img_data.shape[0], img_data.shape[1], 1)), (0, 1, 2)).reshape(len(i)*len(j)*len(k))
					img_ijk_lin = np.zeros((3,len(i2)),dtype=float)
					img_ijk_lin[0,:] = i2
					img_ijk_lin[1,:] = j2
					img_ijk_lin[2,:] = k2
					img_xyz_lin = np.dot(R_img, img_ijk_lin) + np.transpose(np.tile(img_pvec,(len(i2), 1)), (1, 0))
					img_xyz = img_xyz_lin.reshape(3,len(i),len(j),len(k))

					# Load segmentation image
					segloadname = self.quantBrainSegFiles[ID_index]
					seg = nib.nifti1.load(segloadname)
					seg_data = np.asanyarray(seg.dataobj).squeeze()
					seg_csf = np.array(seg_data == 1, dtype=int)  # check these tissue type assignments
					seg_gm = np.array(seg_data == 2, dtype=int)
					seg_wm = np.array(seg_data == 3, dtype=int)

					# File sanity check
					print('')
					print("anat:\t{}".format(imgloadname))
					print(" seg:\t{}".format(segloadname))

					# Calculate voxel tissue fractions
					spec_xyz_lin = np.dot(sup_rda.vox_affine, img_xyz_lin - np.transpose(np.tile(sup_rda.vox_center,(len(i2), 1)), (1, 0))) + np.transpose(np.tile(sup_rda.vox_center,(len(i2), 1)), (1, 0))
					spec_xyz = spec_xyz_lin.reshape(3,len(i),len(j),len(k))

					xgrab = np.array(abs(spec_xyz_lin[0,:] - sup_rda.vox_center[0]) < sup_rda.vox_size[0]/2,dtype=int)
					ygrab = np.array(abs(spec_xyz_lin[1,:] - sup_rda.vox_center[1]) < sup_rda.vox_size[1]/2,dtype=int)
					zgrab = np.array(abs(spec_xyz_lin[2,:] - sup_rda.vox_center[2]) < sup_rda.vox_size[2]/2,dtype=int)
					vox = np.array((xgrab * ygrab * zgrab).reshape(len(i),len(j),len(k)),dtype=np.int16)

					voxfile = nib.Nifti1Image(vox, img_qform)
					voxfile.to_filename(self.quantFileList[ID_index].replace('_sup.out', "_voxel_overlay.nii.gz"))

					vox_frac = [0,0,0]
					vox_frac[0] = (seg_gm * vox).sum() / float(vox.sum())
					vox_frac[1] = (seg_wm * vox).sum() / float(vox.sum())
					vox_frac[2] = (seg_csf * vox).sum() / float(vox.sum())

					print('')
					print("voxfrac:\t", vox_frac)
					print('')
					out_file.write(str(vox_frac[0]) + ',' + str(vox_frac[1]) + ',' + str(vox_frac[2]) + ',')

					# Get number of averages
					n_avg_sup = sup_rda.n_averages
					n_avg_uns = unsup_rda.n_averages

					# Get scaling factors
					scale_sup = sup_rda.ConvS[0]
					scale_uns = unsup_rda.ConvS[0]

					print('n_avg_sup:\t', n_avg_sup)
					print('n_avg_uns:\t', n_avg_uns)
					print('scale_sup:\t', scale_sup)
					print('scale_uns:\t', scale_uns)
					print('')
					out_file.write(str(n_avg_sup) + ',' + str(n_avg_uns) + ',' + str(scale_sup) + ',' + str(scale_uns) + ',')

					# Get metabolite parameters
					metab_params = self.tree()
					num_params   = self.metabParamsTableWidget.rowCount(); print('num_params:\t', num_params)
					for param_index in range(0, num_params):
						metab_params[param_index][0] = str(self.metabParamsTableWidget.item(param_index,0).text())
						metab_params[param_index][1] = float(self.metabParamsTableWidget.item(param_index,1).text())
						metab_params[param_index][2] = float(self.metabParamsTableWidget.item(param_index,2).text())
						metab_params[param_index][3] = float(self.metabParamsTableWidget.item(param_index,3).text())
						metab_params[param_index][4] = float(self.metabParamsTableWidget.item(param_index,4).text())
						metab_params[param_index][5] = float(self.metabParamsTableWidget.item(param_index,5).text())
						metab_params[param_index][6] = int(self.metabParamsTableWidget.item(param_index,6).text())
						metab_params[param_index][7] = int(self.metabParamsTableWidget.item(param_index,7).text())

					# Get water parameters
					water_params = ('water\t' \
						+ self.protonsLineEdit_water.text() + '\t' \
						+ self.T1GMLineEdit_water.text() + '\t' \
						+ self.T2GMLineEdit_water.text() + '\t' \
						+ self.T1WMLineEdit_water.text() + '\t' \
						+ self.T2WMLineEdit_water.text() + '\t' \
						+ self.T1CSFLineEdit_water.text() + '\t' \
						+ self.T2CSFLineEdit_water.text()).split('\t')
					water_params[0] = str(water_params[0])
					water_params[1] = float(water_params[1])
					water_params[2] = float(water_params[2])
					water_params[3] = float(water_params[3])
					water_params[4] = float(water_params[4])
					water_params[5] = float(water_params[5])
					water_params[6] = float(water_params[6])
					water_params[7] = float(water_params[7])

					# Get experimental parameters
					exp_params = ('exp\t' \
						+ self.TRLineEdit.text() + '\t' \
						+ self.TELineEdit.text() + '\t' \
						+ self.waterConcLineEdit.text() + '\t' \
						+ self.waterConcGMLineEdit.text() + '\t' \
						+ self.waterConcWMLineEdit.text() + '\t' \
						+ str(int(self.voxelConcButton.isChecked()))).split('\t')
					exp_params[0] = str(exp_params[0])
					exp_params[1] = float(exp_params[1])
					exp_params[2] = float(exp_params[2])
					exp_params[3] = float(exp_params[3])
					exp_params[4] = float(exp_params[4])
					exp_params[5] = float(exp_params[5])
					exp_params[6] = int(exp_params[6])
					exp_params.append(float(self.waterConcCSFLineEdit.text()))

					# Get scanner type
					scanner_type = 'siemens' if self.siemensScanner.isChecked() else 'varian'
					out_file.write(str(scanner_type) + ',')

					# Calculate absolute metabolite levels (in mM)
					f_conc, f_crlb = mc.calc(sup_out, unsup_out, \
								 vox_frac, n_avg_sup, n_avg_uns, scale_sup, scale_uns, \
								 1, 1, \
								 metab_params, num_params, water_params, exp_params, \
								 scanner_type)
					#   | note: gain_sup and gain_uns are set to 1 for Siemens scanners

					# Find image voxel closest to MRS voxel isocenter
					dist_to_isocenter = np.sqrt(np.sum((img_xyz_lin - np.transpose(np.tile(sup_rda.vox_center,(len(i2), 1)), (1, 0)))**2, 0))
					k_vox_isocenter = np.unravel_index(dist_to_isocenter.argmin(), img_data.shape)

					# Create masked array versions of vox to overlay on top of image
					voxmask = np.ma.masked_where(vox == 0, vox * img_data.max() + 2)

					# Display images of MRS voxel overlay for cor, sag, tra through isocenter point
					transpose_indices = [0, 1]
					n_rot90 = [1, 1, 1]
					panes = [2, 0, 1]

					fig = plt.figure(1)
					fig.patch.set_facecolor('black')
					gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1])
					gs.update(wspace=0.1)
					ax1 = plt.subplot(gs[panes[0]])
					ax2 = plt.subplot(gs[panes[1]])
					ax3 = plt.subplot(gs[panes[2]])
					ax4 = plt.subplot(gs[3])
					fig.axes[panes.index(1)].set_title(ID + "\n", size=12, color='w', family='monospace')
					fig.axes[panes.index(1)].set_xlabel("\nGM: {:.3f}, WM: {:.3f}, CSF: {:.3f}".format(vox_frac[0], vox_frac[1], vox_frac[2]), size=11, color='w', family='monospace')

					palette = cm.Greys_r
					palette.set_over('g', 0.6)

					anat_zoom1, k1 = autozoom(np.transpose(img_data[:,k_vox_isocenter[1],:].squeeze(),transpose_indices))
					anat_zoom2, k2 = autozoom(np.transpose(img_data[:,:,k_vox_isocenter[2]].squeeze(),transpose_indices))
					anat_zoom3, k3 = autozoom(np.transpose(img_data[k_vox_isocenter[0],:,:].squeeze(),transpose_indices))

					ax1.imshow(np.rot90(anat_zoom1, n_rot90[0]), 
								 cmap = palette, norm = colors.Normalize(vmin = img_data.min() - 1, 
								 vmax = img_data.max() + 1, clip = False), aspect='equal')  #, extent=[i_min, i_max, k_min, k_max])
					ax2.imshow(np.rot90(anat_zoom2, n_rot90[1]), 
								 cmap = palette, norm = colors.Normalize(vmin = img_data.min() - 1, 
								 vmax = img_data.max() + 1, clip = False), aspect='equal')
					ax3.imshow(np.rot90(anat_zoom3, n_rot90[2]), 
								 cmap = palette, norm = colors.Normalize(vmin = img_data.min() - 1, 
								 vmax = img_data.max() + 1, clip = False), aspect='equal')
					ax4.imshow(np.zeros((10,10)), cmap = palette)

					ax4.text(0, 0, "Metabolite Levels\n_________________\n\n" + 
						''.join(["{}: {:6.3f} mM\n".format(key, f_conc[key]) for key in f_conc]), 
						ha="left", va="top", size=8, color='w', family='monospace')


					ax1.imshow(np.rot90(np.transpose(voxmask[min(k1[0]):max(k1[0])+1,k_vox_isocenter[1],min(k1[1]):max(k1[1])+1].squeeze(),transpose_indices), n_rot90[0]), 
								 cmap = palette, norm = colors.Normalize(vmin = img_data.min() - 1, 
								 vmax = img_data.max() + 1, clip = False), aspect='equal')
					ax2.imshow(np.rot90(np.transpose(voxmask[min(k2[0]):max(k2[0])+1,min(k2[1]):max(k2[1])+1,k_vox_isocenter[2]].squeeze(),transpose_indices), n_rot90[1]), 
								 cmap = palette, norm = colors.Normalize(vmin = img_data.min() - 1, 
								 vmax = img_data.max() + 1, clip = False), aspect='equal')
					ax3.imshow(np.rot90(np.transpose(voxmask[k_vox_isocenter[0],min(k3[0]):max(k3[0])+1,min(k3[1]):max(k3[1])+1].squeeze(),transpose_indices), n_rot90[2]),
								 cmap = palette, norm = colors.Normalize(vmin = img_data.min() - 1, 
								 vmax = img_data.max() + 1, clip = False), aspect='equal')

					plt.setp([a.set_xticks([]) for a in fig.axes])
					plt.setp([a.set_yticks([]) for a in fig.axes])
					plt.savefig(self.quantSupFiles[ID_index].replace('_sup.out','') + "_barstool_output.png", dpi=300, facecolor='k', bbox_inches='tight', pad_inches = 0.2)

					for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
						out_file.write("{:6.6f},".format(f_conc[str(self.metabParamsTableWidget.item(metab_index,0).text())]))
					for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
						out_file.write("{:6.6f},".format(f_crlb[str(self.metabParamsTableWidget.item(metab_index,0).text())]))
					out_file.write('\n')

			else:
				self.consoleOutputText.append('Analysis of data from varian scanners are not implemented yet.')
				self.consoleOutputText.append('Please try again.')
				self.runQuantButton.setEnabled(False)
				self.loadOutputsButton_quant.setEnabled(True)
		except Exception as e:
			msg = QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Critical)
			msg.setText('Error ' + str(e) + 'has occurred!')
			msg.setInformativeText(traceback.format_exc() + '\n\n' \
				+ '1. Please check that all necessary files for quantification are present in the working directory.\n' \
				+ '2. Please check that proper parameters have been loaded in the "Set Parameters" tab.')
			msg.setWindowTitle('Error')
			msg.exec_()

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())