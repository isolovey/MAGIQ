from __future__ import print_function

# ---- System Libraries ---- #
from builtins import str
from builtins import range
import sys
import os
import subprocess

# Check for PyQt6 (for native mac M1 compatibility), otherwise continue using PyQt5
import importlib

from nibabel.filebasedimages import ImageFileError

PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore, QtGui, QtWidgets, uic
else:
	from PyQt5 import QtCore, QtGui, QtWidgets, uic

from collections import defaultdict

# ---- Math Libraries ---- #
import scipy as sp
from scipy import stats
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

# ---- Matlab ---- #
import matlab.engine

qtCreatorFile = "barstoolrv/ui/BARSTOOLRV.ui"
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
		self.setBindings('Brain Extraction and Segmentation')
		self.setBindings('Set Parameters')
		self.setBindings('Quantify Metabolites')

		#If running on Windows, set the WSL environment up so FSL can be used
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

		elif tab == 'Brain Extraction and Segmentation':

			# --- For legacy VARIAN datasets --- #
			self.FDFIMAGELIST_INDEX = 0

			self.selectFDFImagesButton.clicked.connect(self.loadMouseDirs)
			self.selectFDFImagesButton.setEnabled(False)

			self.fdf2niftiButton.clicked.connect(self.fdf2nifti)
			self.fdf2niftiButton.setEnabled(False)

			self.runVoxAlignButton.clicked.connect(self.runVoxAlign)
			self.runVoxAlignButton.setEnabled(False)

			self.runPCNNButton.clicked.connect(self.runPCNN)
			self.runPCNNButton.setEnabled(False)

			self.runSegButton.clicked.connect(self.runSeg)
			self.runSegButton.setEnabled(False)

			# --- For new BRUKER datasets --- #
			self.selectBrukerDatasetsButton.clicked.connect(self.loadBrukerDatasets)
			self.selectBrukerDatasetsButton.setEnabled(False)

			self.runBrainExtButton.clicked.connect(self.runBrainExtBruker)
			self.runBrainExtButton.setEnabled(False)

			self.runVoxAlignButton_bruker.clicked.connect(self.runVoxAlignBruker)
			self.runVoxAlignButton_bruker.setEnabled(False)

			self.runSegButton_bruker.clicked.connect(self.runSeg_bruker)
			self.runSegButton_bruker.setEnabled(False)

		elif tab == 'Set Parameters':

			self.loadMetabParamsButton.clicked.connect(self.loadMetabParams)
			
			self.saveMetabParamsButton.clicked.connect(self.saveMetabParams)
			self.saveMetabParamsButton.setEnabled(False)

			self.confirmParamsButton.clicked.connect(self.verifyParams)

		elif tab == 'Quantify Metabolites':

			# --- For legacy VARIAN datasets --- #
			self.loadOutputsButton_quant.clicked.connect(self.loadMouseDirs)
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

			# --- For new BRUKER datasets --- #
			self.loadOutputsButton_quant_bruker.clicked.connect(self.loadBrukerDatasets)
			self.loadOutputsButton_quant_bruker.setEnabled(False)

			self.confirmSaveFileButton_quant_bruker.clicked.connect(self.setQuantSaveFile_bruker)
			self.confirmSaveFileButton_quant_bruker.setEnabled(False)

			self.runQuantButton_bruker.clicked.connect(self.runQuant_bruker)
			self.runQuantButton_bruker.setEnabled(False)

			# Set Up Plotting Region
			fig_bruker = plt.figure(2)
			self.canvas_bruker = FigureCanvas(fig)
			self.plotQuant_mplvl_bruker.addWidget(self.canvas_bruker)
			self.canvas_bruker.draw()

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

		# --- For legacy VARIAN datasets --- #
		self.selectFDFImagesButton.setEnabled(True)
		self.loadOutputsButton_quant.setEnabled(True)

		# --- For new BRUKER datasets --- #
		self.selectBrukerDatasetsButton.setEnabled(True)
		self.loadOutputsButton_quant_bruker.setEnabled(True)

	def loadOutputs(self):
		'''
			This method presents a dialog allowing the user to select which *.out files to analyze.
		'''
		self.consoleOutputText.append('===== CALCULATE AMPLITUDES AND CRLBS =====')

		self.outputs = []

		file_dialog = QtWidgets.QFileDialog(directory=self.workingDirectory)
		file_dialog.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
		file_dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
		file_view = file_dialog.findChild(QtWidgets.QListView, 'listView')

		# to make it possible to select multiple directories:
		if file_view:
			file_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
		f_tree_view = file_dialog.findChild(QtWidgets.QTreeView)
		
		if f_tree_view:
			f_tree_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

		if file_dialog.exec_():
			paths = file_dialog.selectedFiles()

			self.fileList = paths

			self.consoleOutputText.append('The following output files were loaded:')
			for (i, file) in enumerate(self.fileList):
				self.fileList[i] = str(file) + '/sup.out'
				self.consoleOutputText.append(' >> ' + str(self.fileList[i]))
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

	# ---- Methods for 'Brain Extraction' Tab (BRUKER) ---- #
	def loadBrukerDatasets(self):

		if self.mainTabWidget.currentIndex() == 1:
			self.consoleOutputText.append('==== BRAIN EXTRACTION (BRUKER) ====')
		elif self.mainTabWidget.currentIndex() == 3:
			self.consoleOutputText.append('==== QUANTIFICATION (BRUKER) ====')

		file_dialog = QtWidgets.QFileDialog(directory=self.workingDirectory)
		file_dialog.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
		file_dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
		file_view = file_dialog.findChild(QtWidgets.QListView, 'listView')

		# to make it possible to select multiple directories:
		if file_view:
			file_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
		f_tree_view = file_dialog.findChild(QtWidgets.QTreeView)

		if file_dialog.exec():
			paths = file_dialog.selectedFiles()
			self.mouseDirsBruker = paths
			
			self.consoleOutputText.append('The following files were loaded:')
			for (i, mouse) in enumerate(self.mouseDirsBruker):
				self.mouseDirsBruker[i] = str(mouse)
				self.consoleOutputText.append(' >> ' + str(mouse))
			self.mouseIDsBruker = [mouse.split('/')[-1] for mouse in self.mouseDirsBruker]
			self.consoleOutputText.append('')

			if self.mainTabWidget.currentIndex() == 1:
				self.selectBrukerDatasetsButton.setEnabled(False)
				self.runBrainExtButton.setEnabled(True)
			elif self.mainTabWidget.currentIndex() == 3:
				self.loadOutputsButton_quant_bruker.setEnabled(True)
				self.saveFileLineEdit_quant_bruker.setText(self.workingDirectory + '/' + '_____.csv')
				self.confirmSaveFileButton_quant_bruker.setEnabled(True)
		else:
			if self.mainTabWidget.currentIndex() == 1:
				self.consoleOutputText.append('No image files selected ... try again.')
				self.consoleOutputText.append('')
				self.selectBrukerDatasetsButton.setEnabled(True)
			elif self.mainTabWidget.currentIndex() == 3:
				self.consoleOutputText.append('No output files selected ... try again.')
				self.consoleOutputText.append('')
				self.loadOutputsButton_quant_bruker.setEnabled(False)

	def runBrainExtBruker(self):
		for d, directory in enumerate(self.mouseDirsBruker):
			img_file = directory + '/' + self.mouseIDsBruker[d] + '.nii.gz'
			mas_file = directory + '/' + self.mouseIDsBruker[d] + '_mask.nii.gz'
			bra_file = directory + '/' + self.mouseIDsBruker[d] + '_brain.nii.gz'

			self.consoleOutputText.append('==== RATS_MM ====')
			self.consoleOutputText.append('Loading ' + str(img_file))

			try:
				img = nib.nifti1.load(img_file)
			except FileNotFoundError:
				file_dialog = QtWidgets.QFileDialog(directory=self.workingDirectory)
				file_dialog.setWindowTitle("Select structural image in Nifti or DICOM format")
				file_dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
				file_view = file_dialog.findChild(QtWidgets.QListView, 'listView')

				if file_dialog.exec():
					paths = file_dialog.selectedFiles()
				else:
					raise RuntimeError("No valid structural Nifti file selected for brain extraction")

				try:
					img = nib.nifti1.load(paths[0])
					os.link(paths[0], img_file)
				except ImageFileError:
					cmd = [
						'dcm2niix', '-z', 'y', '-o', os.path.dirname(img_file),
						'-f', os.path.basename(img_file)[:-len('.nii.gz')],
						paths[0]
					]
					self.consoleOutputText.append(' Converting to DICOM to Nifti: ' + ' '.join(cmd))
					self.consoleOutputText.append(subprocess.check_output(cmd).decode() + '\n')
					img = nib.nifti1.load(img_file)

			img_data = np.asanyarray(img.dataobj)
	
			k = int(self.ratsMM_k.text())
			t = int(np.mean(img_data)); self.ratsMM_t.setText(str(int(np.mean(img_data))))
			v = int(self.ratsMM_v.text())

			cmd = [f'{os.getcwd()}/barstoolrv/RATS_MM', '-k', str(k), '-t', str(t), '-v', str(v)]
			if os.name == 'nt':
				img_file_path = subprocess.check_output('wsl wslpath -u "' + img_file + '"').decode().rstrip()
				mas_file_path = subprocess.check_output('wsl wslpath -u "' + mas_file + '"').decode().rstrip()
				self.consoleOutputText.append(' | ' + str(img_file_path))
				cmd = ['wsl'] + cmd + [img_file_path, mas_file_path]
			else:
				cmd += [img_file, mas_file]
			self.consoleOutputText.append(' >> ' + ' '.join(cmd))
			self.consoleOutputText.append(subprocess.check_output(cmd).decode() + '\n')

			mas = nib.nifti1.load(mas_file)
			mas_data = np.asanyarray(mas.dataobj)
			bra_data = img_data * mas_data
			bra = nib.Nifti1Image(bra_data, img.affine, img.header)
			bra.to_filename(bra_file)

		self.runBrainExtButton.setEnabled(False)
		self.runVoxAlignButton_bruker.setEnabled(True)

	def runVoxAlignBruker(self):
		self.consoleOutputText.append('===== alignvoxel_bruker =====')
		for d, directory in enumerate(self.mouseDirsBruker):
			file_dir = directory + '/sup'
			print('Loading ', file_dir)
			
			data = BrukerFID(file_dir)

			VoxArrSize = np.array([float(v) for v in data.header['PVM_VoxArrSize']['value']])
			VoxArrPosition = np.array([float(v) for v in data.header['PVM_VoxArrPosition']['value']])
			VoxArrPositionRPS = np.array([float(v) for v in data.header['PVM_VoxArrPositionRPS']['value']])
			VoxArrCSDisplacement = np.array([float(v) for v in data.header['PVM_VoxArrCSDisplacement']['value']])
			VoxArrGradOrient = np.array([list(map(float, sublist)) for sublist in data.header['PVM_VoxArrGradOrient']['value']])

			img_file = directory + '/' + self.mouseIDsBruker[d] + '.nii.gz'
			bra_file = directory + '/' + self.mouseIDsBruker[d] + '_brain.nii.gz'
			vox_file = directory + '/' + self.mouseIDsBruker[d] + '_voxel_overlay.nii.gz'
			
			print(' | Getting ', img_file, bra_file)
			img = nib.nifti1.load(img_file)
			bra = nib.nifti1.load(bra_file)

			bra_ijk = []
			bra_xyz = []
			bra_size = np.array(bra.header.get_data_shape())
			bra_affine_inv = np.linalg.inv(bra.affine)

			M = bra.affine[:3,:3]; abc = bra.affine[:3,3]
			M_inv = bra_affine_inv[:3,:3]; abc_inv = bra_affine_inv[:3,3]

			print('   | M', M)
			print('   | M_inv', M_inv)
			print('   | abc', abc)
			print('   | abc_inv', abc)

			R_vox = VoxArrGradOrient
			p_vox_water = R_vox.dot(VoxArrPosition)
			# p_vox_fat   = R_vox.dot(VoxArrPosition + VoxArrCSDisplacement)

			print('   | p_vox_water', p_vox_water)
			# print('   | p_vox_fat', p_vox_fat)
			print('   | VoxArrSize', VoxArrSize)

			vox_size = VoxArrSize
			vox_img_water  = np.zeros(np.array(img.header.get_data_shape())) 
			# vox_img_fat    = np.zeros(np.array(img.header.get_data_shape()))
			vox_res = np.array(img.header.get_zooms())

			vlminw = p_vox_water - VoxArrSize/2; # vlminf = p_vox_fat - VoxArrSize/2
			vlmaxw = p_vox_water + VoxArrSize/2; # vlmaxf = p_vox_fat + VoxArrSize/2

			print('   | vlminw', vlminw)
			print('   | vlmaxw', vlmaxw)
			# print('   | vlminf', vlminf)
			# print('   | vlmaxf', vlmaxf)

			# Build Water Voxel
			vlminw_ijk = np.round(M_inv.dot(vlminw) + abc_inv).astype(int)
			vlmaxw_ijk = np.round(M_inv.dot(vlmaxw) + abc_inv).astype(int)

			print('   | vlminw_ijk', vlminw_ijk)
			print('   | vlmaxw_ijk', vlmaxw_ijk)

			i1 = np.min([vlminw_ijk[0], vlmaxw_ijk[0]]); i2 = np.max([vlminw_ijk[0], vlmaxw_ijk[0]])
			j1 = np.min([vlminw_ijk[1], vlmaxw_ijk[1]]); j2 = np.max([vlminw_ijk[1], vlmaxw_ijk[1]])
			k1 = np.min([vlminw_ijk[2], vlmaxw_ijk[2]]); k2 = np.max([vlminw_ijk[2], vlmaxw_ijk[2]])

			print(i1, j1, k1)
			print(i2, j2, k2)

			vox_img_water[i1:i2+1, j1:j2+1, k1:k2+1] = int(1)

			print(' | Saving ', vox_file)
			nifti_vox_water = nib.Nifti1Image(vox_img_water, bra.affine, bra.header)
			nifti_vox_water.to_filename(vox_file)

			self.consoleOutputText.append(' >> ' + str(directory))
		print('')
		self.consoleOutputText.append('')

		self.runVoxAlignButton_bruker.setEnabled(False)
		self.runSegButton_bruker.setEnabled(True)

	def runSeg_bruker(self):
		self.consoleOutputText.append('==== csf_thresh (BRUKER) ====')
		for d, directory in enumerate(self.mouseDirsBruker):

			bra_file = directory + '/' + self.mouseIDsBruker[d] + '_brain.nii.gz'
			mas_file = directory + '/' + self.mouseIDsBruker[d] + '_mask.nii.gz'
			vox_file = directory + '/' + self.mouseIDsBruker[d] + '_voxel_overlay.nii.gz'
			csf_file = directory + '/' + self.mouseIDsBruker[d] + '_csf_mask.nii.gz'

			self.consoleOutputText.append(' >> ' + str(directory))
			print('Processing ', directory, '...')

			brain = nib.load(bra_file)
			mask  = nib.load(mas_file)
			vox   = nib.load(vox_file)

			brain_img = np.asanyarray(brain.dataobj).squeeze()
			mask_img  = np.asanyarray( mask.dataobj).squeeze()
			vox_img   = np.asanyarray(  vox.dataobj).squeeze()

			brain_img_vec = np.reshape(brain_img, np.size(brain_img)).astype(int)
			mask_img_vec  = np.reshape( mask_img, np.size(mask_img )).astype(int)
			vox_img_vec   = np.reshape(  vox_img, np.size(vox_img  ))
			vox_img_vec   = np.array([np.round(v) for v in vox_img_vec]).astype(int)

			brain_img_vec_masked = brain_img_vec[mask_img_vec.astype(bool)]
			brain_kde = stats.gaussian_kde(brain_img_vec_masked)

			brain_img_vec_vox_masked = brain_img_vec[vox_img_vec.astype(bool)]
			vox_kde = stats.gaussian_kde(brain_img_vec_vox_masked)

			if self.csfThreshMode1_bruker.isChecked():
				for i, elem in enumerate(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000)):
					if brain_kde.integrate_box_1d(np.min(brain_img_vec), elem) > float(self.csfThreshLineEdit_bruker.text()):
						csf_thresh = elem
						break
			elif self.csfThreshMode2_bruker.isChecked():
				for i, elem in enumerate(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000)):
					if vox_kde.integrate_box_1d(np.min(brain_img_vec), elem) > float(self.csfThreshLineEdit_bruker.text()):
						csf_thresh = elem
						break


			plt.figure()
			plt.plot(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000), brain_kde(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000)))
			plt.plot(csf_thresh, brain_kde(csf_thresh), '.')
			plt.title('Gaussian Kernel Density Estimate of PDF')
			plt.xlim(0, np.max(brain_img_vec))
			plt.xlabel('Image Intensity')
			plt.ylabel('Probability')
			plt.legend(['PDF', 'Threshold: '+str(int(csf_thresh))])

			if self.csfThreshMode1_bruker.isChecked():
				plt.savefig(directory + '/' + self.mouseIDsBruker[d] + '_brain_gkde.pdf')
			elif self.csfThreshMode2_bruker.isChecked():
				plt.savefig(directory + '/' + self.mouseIDsBruker[d] + '_vox_gkde.pdf')

			brain_csf_mask = brain_img >= int(csf_thresh)
			brain_csf_mask = brain_csf_mask.astype(int)
			brain_csf_mask_file = nib.Nifti1Image(brain_csf_mask, brain.affine, brain.header)
			brain_csf_mask_file.to_filename(csf_file)

		print('')
		self.runSegButton_bruker.setEnabled(False)
		self.selectBrukerDatasetsButton.setEnabled(True)

	# ---- Methods for 'Brain Extraction' Tab (VARIAN) ---- #
	def loadMouseDirs(self):
		
		if self.mainTabWidget.currentIndex() == 1:
			self.consoleOutputText.append('===== BRAIN EXTRACTION =====')
		elif self.mainTabWidget.currentIndex() == 3:
			self.consoleOutputText.append('===== QUANTIFICATION =====')

		file_dialog = QtWidgets.QFileDialog(directory=self.workingDirectory, filter='Image Files (*.img)')
		file_dialog.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
		file_dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
		file_view = file_dialog.findChild(QtWidgets.QListView, 'listView')

		# to make it possible to select multiple directories:
		if file_view:
			file_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
		f_tree_view = file_dialog.findChild(QtWidgets.QTreeView)
		
		if f_tree_view:
			f_tree_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

		if file_dialog.exec_():
			paths = file_dialog.selectedFiles()

			self.mouseDirs = paths #QtWidgets.QFileDialog.getOpenFileNames(self, 'Open FDF Image File', self.workingDirectory, 'Image Files (*.img)')[0]
		
			self.consoleOutputText.append('The following files were loaded:')
			for (i, mouse) in enumerate(self.mouseDirs):
				self.mouseDirs[i] = str(mouse)
				self.consoleOutputText.append(' >> ' + str(mouse))
			self.consoleOutputText.append('')

			if self.mainTabWidget.currentIndex() == 1:
				self.selectFDFImagesButton.setEnabled(False)
				self.fdf2niftiButton.setEnabled(True)
			elif self.mainTabWidget.currentIndex() == 3:
				self.loadOutputsButton_quant.setEnabled(False)
				self.saveFileLineEdit_quant.setText(self.workingDirectory + '/' + '_____.csv')
				self.confirmSaveFileButton_quant.setEnabled(True)
		else:
			
			if self.mainTabWidget.currentIndex() == 1:
				self.consoleOutputText.append('No image files selected ... try again.')
				self.consoleOutputText.append('')
				self.selectFDFImagesButton.setEnabled(True)
			elif self.mainTabWidget.currentIndex() == 3:
				self.consoleOutputText.append('No output files selected ... try again.')
				self.consoleOutputText.append('')
				self.loadOutputsButton_quant.setEnabled(False)

	def fdf2nifti(self):

		self.consoleOutputText.append('===== fdf2nifti =====')
		self.fdf_imgs = []
		for mouse in self.mouseDirs:
			fdf_img   = FDF2D(mouse + '/fse2D.img', (int(self.matXLineEdit.text()), int(self.matYLineEdit.text()), int(self.matZLineEdit.text())))
			self.fdf_imgs.append(fdf_img)
			nifti_img = nib.Nifti1Image(fdf_img.fseimg, fdf_img.affine)
			nifti_img.to_filename(mouse + '/fse2d.nii.gz')
			self.consoleOutputText.append(' >> ' + str(mouse))
		self.consoleOutputText.append('')

		self.fdf2niftiButton.setEnabled(False)
		self.runVoxAlignButton.setEnabled(True)

	def runVoxAlign(self):
		self.consoleOutputText.append('===== alignvoxel_varian =====')
		for i, mouse in enumerate(self.mouseDirs):
			fdf_img = self.fdf_imgs[i]
			try:
				print('Loading ', mouse + '/metab.fid ...')
				voxel = VarianVoxel(mouse + '/metab.fid', fdf_img.size, fdf_img.X_VARIAN, fdf_img.Y_VARIAN, fdf_img.Z_VARIAN, fdf_img.fseimg_ijk, fdf_img.fseimg_xyz_kdt)
			except Exception as e:
				print('Error: ', e)
				print('Loading ', mouse + '/water.fid ...')
				voxel = VarianVoxel(mouse + '/water.fid', fdf_img.size, fdf_img.X_VARIAN, fdf_img.Y_VARIAN, fdf_img.Z_VARIAN, fdf_img.fseimg_ijk, fdf_img.fseimg_xyz_kdt)
			nifti_img = nib.Nifti1Image(voxel.voximg, fdf_img.affine)
			nifti_img.to_filename(mouse + '/mrsvoxel.nii.gz')
			self.consoleOutputText.append(' >> ' + str(mouse))
		print('')
		self.consoleOutputText.append('')

		self.runVoxAlignButton.setEnabled(False)
		self.runPCNNButton.setEnabled(True)

	def runPCNN(self):
		self.consoleOutputText.append('==== PCNN 3D ====')
		
		# Start Matlab Engine
		print('Starting Matlab Engine ...')
		cwd = os.getcwd()
		eng = matlab.engine.start_matlab()
		eng.cd(cwd + '/barstoolrv')
		matlab_wd = eng.pwd(); print(matlab_wd)

		failed_mice = []
		successful_mice = []
		
		for i, mouse in enumerate(self.mouseDirs):
			self.consoleOutputText.append(' >> ' + str(mouse))

			print('Processing ', mouse, '...')
			fse2d_nifti = mouse + '/fse2d.nii.gz'
			fse2d_mask  = mouse + '/fse2d_mask.nii.gz'
			fse2d_brain = mouse + '/fse2d_brain.nii.gz'

			try:
				eng.runPCNN3D(fse2d_nifti, nargout=0)
				self.consoleOutputText.append('    runPCNN3D(' + str(mouse) + ')')

				eng.clear('all', nargout=0)
				eng.close('all', nargout=0)

				command = ['fslmaths', fse2d_nifti]
				command.append('-mas ' + fse2d_mask)
				command.append(fse2d_brain)
				print(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
				os.system(str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))
				self.consoleOutputText.append('    ' + str(command).replace('[','').replace(']','').replace(',','').replace("'", ''))

				successful_mice.append(mouse)
			except Exception as e:
				print(e)
				failed_mice.append(mouse)

		print('')
		self.consoleOutputText.append('')
		self.consoleOutputText.append('The following files could not be processed successfully:')
		for mouse in failed_mice:
			self.consoleOutputText.append(' >> ' + str(mouse))
		self.consoleOutputText.append('')

		# Stop Matlab Engine
		print('Stopping Matlab Engine ...')
		eng.quit()

		self.runPCNNButton.setEnabled(False)
		self.runSegButton.setEnabled(True)
		self.mouseDirs = successful_mice # remove files that could not be processed successfully

	def runSeg(self):
		self.consoleOutputText.append('==== CSF EXTRACT ====')
		for i, mouse in enumerate(self.mouseDirs):
			self.consoleOutputText.append(' >> ' + str(mouse))
			print('Processing ', mouse, '...')

			brain = nib.load(mouse + '/fse2d_brain.nii.gz')
			mask  = nib.load(mouse + '/fse2d_mask.nii.gz')

			brain_img = np.asanyarray(brain.dataobj).squeeze()
			mask_img  = np.asanyarray( mask.dataobj).squeeze()

			brain_img_vec = np.reshape(brain_img, np.size(brain_img)).astype(int)
			mask_img_vec  = np.reshape(mask_img,  np.size(mask_img )).astype(int)

			brain_img_vec_masked = brain_img_vec[mask_img_vec.astype(bool)]
			brain_kde = stats.gaussian_kde(brain_img_vec_masked)

			for i, elem in enumerate(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000)):
				if brain_kde.integrate_box_1d(np.min(brain_img_vec), elem) > float(self.csfThreshLineEdit.text()):
					print('  | ', i, elem)
					csf_thresh = elem
					break

			plt.figure()
			plt.plot(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000), brain_kde(np.linspace(np.min(brain_img_vec), np.max(brain_img_vec), 1000)))
			plt.plot(csf_thresh, brain_kde(csf_thresh), '.')
			plt.title('Gaussian Kernel Density Estimate of PDF')
			plt.xlim(0, np.max(brain_img_vec))
			plt.xlabel('Image Intensity')
			plt.ylabel('Probability')
			plt.legend(['PDF', 'Threshold: '+str(int(csf_thresh))])
			plt.savefig(mouse + '/fse2d_brain_gkde.pdf')

			brain_csf_mask = brain_img >= int(csf_thresh)
			brain_csf_mask = brain_csf_mask.astype(int)
			brain_csf_mask_file = nib.Nifti1Image(brain_csf_mask, brain.get_qform())
			brain_csf_mask_file.to_filename(mouse + '/fse2d_csf_mask.nii.gz')

		print('')
		self.runSegButton.setEnabled(False)
		self.selectFDFImagesButton.setEnabled(True)

		self.consoleOutputText.append('')

	# ---- Methods for Setting Quantification Parameters ---- #
	def loadMetabParams(self):
		prev = str(self.metabParamsFileLineEdit.text())
		self.metabParamsFileLineEdit.setText(str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open Quantification Information File', os.getcwd() + '/barstoolrv/qinfo', 'Quantification Info Files (*.qinfo)')[0]))
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

		self.metabParamsTableWidget.setRowCount(sp.size(rows,0))
		self.metabParamsTableWidget.setColumnCount(8)

		for i in range(0, sp.size(rows, 0)):
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

	def saveMetabParams(self):
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
		self.consoleOutputText.append('The following parameters were entered. Please check them carefully:')
		self.consoleOutputText.append('')
		self.consoleOutputText.append('\tProtons\tT1 (GM) [sec]\tT2 (GM) [ms]\tT1 (WM) [ms]\tT2 (WM) [ms]\tFirst Peak\tLast Peak')
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

	# ---- Methods for Actual Quantification (VARIAN) ---- #
	def setQuantSaveFile(self):
		self.quantSaveFileName = self.saveFileLineEdit_quant.text()
		self.runQuantButton.setEnabled(True)
		self.consoleOutputText.append('Quantification results will be saved to: ' + str(self.quantSaveFileName))
		self.consoleOutputText.append('')

	def runQuant(self):
		self.confirmSaveFileButton_quant.setEnabled(False)

		out_file = open(self.quantSaveFileName, 'w')

		# Write header
		out_file.write('ID,TISSUE,CSF,N_AVG_SUP,N_AVG_UNS,SCALE_SUP,SCALE_UNS,SCANNER,')
		for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
			out_file.write(str(self.metabParamsTableWidget.item(metab_index,0).text())+',')
		for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
			out_file.write(str(self.metabParamsTableWidget.item(metab_index,0).text())+'_CRLB,')
		out_file.write('\n')

		self.consoleOutputText.append('==== METAB QUANT ====')
		failed_mice = []
		for i, mouse in enumerate(self.mouseDirs):

			try:
				ID = mouse.split('/')[-1]
				out_file.write(str(ID)+',')

				self.consoleOutputText.append(' >> ' + str(mouse))
				print('Processing ', mouse, '...')

				brain = nib.load(mouse + '/fse2d_brain.nii.gz')
				csf   = nib.load(mouse + '/fse2d_csf_mask.nii.gz')
				vox   = nib.load(mouse + '/mrsvoxel.nii.gz')
				sup_out     = OutputFile(mouse + '/sup.out')
				unsup_out   = OutputFile(mouse + '/uns.out')
				sup_dat   = DatFile(mouse + '/sup.dat')
				unsup_dat = DatFile(mouse + '/metab_uns.dat')

				brain_img = brain.get_fdata()
				csf_img   = csf.get_fdata()
				vox_img   = vox.get_fdata()

				vox_img_vec = np.reshape(vox_img, np.size(vox_img)).astype(int)
				csf_img_vec = np.reshape(csf_img, np.size(csf_img)).astype(int)

				vox_n = np.sum(vox_img_vec)
				csf_n = np.sum(vox_img_vec[csf_img_vec.astype(bool)])

				tissue_frac = 1-float(csf_n)/float(vox_n)
				vox_frac = [tissue_frac/2, tissue_frac/2, 1-tissue_frac]

				print("voxfrac:\t", tissue_frac, vox_frac)
				out_file.write(str(tissue_frac) + ',' + str(vox_frac[2]) + ',')

				# Get number of averages
				procpar_sup = Procpar(mouse + '/metab.fid/procpar')
				n_avg_sup = int(procpar_sup.acqcycles)//2

				procpar_uns = Procpar(mouse + '/water.fid/procpar')
				n_avg_uns = int(procpar_uns.acqcycles)//2

				# procpar_sup = open(mouse + '/metab.fid/procpar', 'r')
				# for line in procpar_sup:
				# 	if 'acqcycles' in line:
				# 		n_avg_sup = int(procpar_sup.next().split(' ')[1])/2

				# procpar_uns = open(mouse + '/water.fid/procpar', 'r')
				# for line in procpar_uns:
				# 	if 'acqcycles' in line:
				# 		n_avg_uns = int(procpar_uns.next().split(' ')[1])/2

				# Get scaling factors -- CHECK WITH BARTHA!
				scale_sup = sup_dat.ConvS
				scale_uns = unsup_dat.ConvS

				gain_sup = procpar_sup.gain
				gain_uns = procpar_uns.gain

				print('n_avg_sup:\t', n_avg_sup, '\t', end=' ')
				print('n_avg_uns:\t', n_avg_uns)
				print('gain_sup:\t', gain_sup, '\t', end=' ')
				print('gain_uns:\t', gain_uns)
				print('sup_ConvS:\t', sup_dat.ConvS, '\t', end=' ')
				print('uns_ConvS:\t', unsup_dat.ConvS)
				print('scale_sup:\t', scale_sup, '\t', end=' ')
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
				scanner_type = 'varian'
				out_file.write(str(scanner_type) + ',')

				# Calculate absolute metabolite levels (in mM)
				f_conc, f_crlb = mc.calc(sup_out, unsup_out, \
							 vox_frac, n_avg_sup, n_avg_uns, scale_sup, scale_uns, \
							 gain_sup, gain_uns, \
							 metab_params, num_params, water_params, exp_params, \
							 scanner_type)

				# Figure out centroid of voxel to slice image appropriately
				vox_indices = np.where(vox_img == 1)
				vox_centroid = np.round([np.mean(vox_indices[0]), np.mean(vox_indices[1]), np.mean(vox_indices[2])]).astype(int)

				# Create masked array version of vox_img to overal on top of brain_img
				vox_mas = np.ma.masked_where(vox_img == 0, vox_img * brain_img.max() + 2)

				# Display images of MRS voxel overlay
				transpose_indices = [0, 1]
				panes = [0, 1, 2]
				n_rot90 = [0, 0, 0]

				fig = plt.figure(1)
				fig.patch.set_facecolor('black')

				ax1 = plt.subplot(1,3,1)
				ax2 = plt.subplot(1,3,2)
				ax3 = plt.subplot(1,3,3)
				fig.axes[panes.index(1)].set_title(ID + "\n", size=12, color='w', family='monospace')
				fig.axes[panes.index(1)].set_xlabel("\nTISSUE: {:.3f}, CSF: {:.3f}".format(tissue_frac, vox_frac[2]), size=11, color='w', family='monospace')

				palette = cm.Greys_r
				palette.set_over('g', 0.6)

				asp = brain.header['pixdim'][1]/brain.header['pixdim'][3]
				print('asp', asp)
				# raw_input()

				print('brain_img', brain_img[vox_centroid[0],:,:].squeeze())
				print('vox_centroid', vox_centroid)
				# raw_input()

				anat_zoom1, k1 = autozoom(np.transpose(brain_img[vox_centroid[0],:,:].squeeze(), transpose_indices))
				# raw_input()

				anat_zoom2, k2 = autozoom(np.transpose(brain_img[:,vox_centroid[1],:].squeeze(), transpose_indices))
				# raw_input()

				anat_zoom3, k3 = autozoom(np.transpose(brain_img[:,:,vox_centroid[2]].squeeze(), transpose_indices))
				# raw_input()
				
				ax1.imshow(np.rot90(anat_zoom1, n_rot90[0]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax2.imshow(np.rot90(anat_zoom2, n_rot90[1]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax3.imshow(np.rot90(anat_zoom3, n_rot90[2]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect='equal')

				ax1.imshow(np.rot90(np.transpose(vox_mas[vox_centroid[0],min(k1[0]):max(k1[0])+1,min(k1[1]):max(k1[1])+1].squeeze(),transpose_indices), n_rot90[0]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax2.imshow(np.rot90(np.transpose(vox_mas[min(k2[0]):max(k2[0])+1,vox_centroid[1],min(k2[1]):max(k2[1])+1].squeeze(),transpose_indices), n_rot90[1]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax3.imshow(np.rot90(np.transpose(vox_mas[min(k3[0]):max(k3[0])+1,min(k3[1]):max(k3[1])+1,vox_centroid[2]].squeeze(),transpose_indices), n_rot90[2]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect='equal')

				plt.setp([a.set_xticks([]) for a in fig.axes])
				plt.setp([a.set_yticks([]) for a in fig.axes])
				plt.savefig(mouse + "/barstool_output.png", dpi=300, facecolor='k', bbox_inches='tight', pad_inches = 0.2)

				# raw_input()

				for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
					out_file.write("{:6.6f},".format(f_conc[str(self.metabParamsTableWidget.item(metab_index,0).text())]))
				for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
					out_file.write("{:6.6f},".format(f_crlb[str(self.metabParamsTableWidget.item(metab_index,0).text())]))
			except Exception as e:
				print(e)
				failed_mice.append(mouse)			
			out_file.write('\n')
		
		out_file.close()

		self.consoleOutputText.append('')
		self.consoleOutputText.append('The following files could not be processed successfully:')
		for mouse in failed_mice:
			self.consoleOutputText.append(' >> ' + str(mouse))

		self.loadOutputsButton_quant.setEnabled(True)
		self.confirmSaveFileButton_quant.setEnabled(False)
		self.runQuantButton.setEnabled(False)

	# ---- Methods for Actual Quantification (BRUKER) ---- #
	def setQuantSaveFile_bruker(self):
		self.quantSaveFileName_bruker = self.saveFileLineEdit_quant_bruker.text()
		self.runQuantButton_bruker.setEnabled(True)
		self.consoleOutputText.append('Quantification results will be saved to: ' + str(self.quantSaveFileName_bruker))
		self.consoleOutputText.append('')

	def runQuant_bruker(self):
		self.confirmSaveFileButton_quant_bruker.setEnabled(False)

		out_file = open(self.quantSaveFileName_bruker, 'w')

		# Write header
		out_file.write('ID,TISSUE,CSF,N_AVG_SUP,N_AVG_UNS,SCALE_SUP,SCALE_UNS,SCANNER,')
		for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
			out_file.write(str(self.metabParamsTableWidget.item(metab_index,0).text())+',')
		for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
			out_file.write(str(self.metabParamsTableWidget.item(metab_index,0).text())+'_CRLB,')
		out_file.write('\n')

		self.consoleOutputText.append('==== METAB QUANT (Bruker) ====')
		failed_mice = []
		for i, mouse in enumerate(self.mouseDirsBruker):

			try:
				ID = mouse.split('/')[-1]
				out_file.write(str(ID)+',')

				self.consoleOutputText.append(' >> ' + str(mouse))
				print('Processing ', mouse, '...')

				brain = nib.load(mouse + '/' + mouse + '_brain.nii.gz')
				csf   = nib.load(mouse + '/' + mouse + '_csf_mask.nii.gz')
				vox   = nib.load(mouse + '/' + mouse + '_voxel_overlay.nii.gz')
				sup_out     = OutputFile(mouse + '/' + mouse + '_sup.out')
				unsup_out   = OutputFile(mouse + '/' + mouse + '_uns.out')
				sup_dat   = DatFile(mouse + '/' + mouse + '_sup.dat')
				unsup_dat = DatFile(mouse + '/' + mouse + '_uns.dat')

				brain_img = brain.get_fdata()
				csf_img   = csf.get_fdata()
				vox_img   = vox.get_fdata()

				vox_img_vec = np.reshape(vox_img, np.size(vox_img)).astype(int)
				csf_img_vec = np.reshape(csf_img, np.size(csf_img)).astype(int)

				vox_n = np.sum(vox_img_vec)
				csf_n = np.sum(vox_img_vec[csf_img_vec.astype(bool)])

				tissue_frac = 1-float(csf_n)/float(vox_n)
				vox_frac = [tissue_frac/2, tissue_frac/2, 1-tissue_frac]

				print("voxfrac:\t", tissue_frac, vox_frac)
				out_file.write(str(tissue_frac) + ',' + str(vox_frac[2]) + ',')

				# Get number of averages
				fid_sup = BrukerFID(mouse + '/sup')
				n_avg_sup = int(fid_sup.header['PVM_NAverages']['value'])//2

				fid_uns = BrukerFID(mouse + '/uns')
				n_avg_uns = int(fid_uns.header['PVM_NAverages']['value'])//2

				# Get scaling factors -- CHECK WITH BARTHA!
				scale_sup = fid_sup.ConvS
				scale_uns = fid_uns.ConvS

				gain_sup = fid_sup.gain()[0]
				gain_uns = fid_uns.gain()[0]

				print('n_avg_sup:\t', n_avg_sup, '\t', end=' ')
				print('n_avg_uns:\t', n_avg_uns)
				print('gain_sup:\t', gain_sup, '\t', end=' ')
				print('gain_uns:\t', gain_uns)
				print('sup_ConvS:\t', sup_dat.ConvS, '\t', end=' ')
				print('uns_ConvS:\t', unsup_dat.ConvS)
				print('scale_sup:\t', scale_sup, '\t', end=' ')
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
				scanner_type = 'bruker'
				out_file.write(str(scanner_type) + ',')

				# Calculate absolute metabolite levels (in mM)
				f_conc, f_crlb = mc.calc(sup_out, unsup_out, \
							 vox_frac, n_avg_sup, n_avg_uns, scale_sup, scale_uns, \
							 gain_sup, gain_uns, \
							 metab_params, num_params, water_params, exp_params, \
							 scanner_type)

				# Figure out centroid of voxel to slice image appropriately
				vox_indices = np.where(vox_img == 1)
				vox_centroid = np.round([np.mean(vox_indices[0]), np.mean(vox_indices[1]), np.mean(vox_indices[2])]).astype(int)

				# Create masked array version of vox_img to overal on top of brain_img
				vox_mas = np.ma.masked_where(vox_img == 0, vox_img * brain_img.max() + 2)

				# Display images of MRS voxel overlay
				transpose_indices = [0, 1]
				panes = [0, 1, 2]
				n_rot90 = [0, 0, 0]

				fig = plt.figure(1)
				fig.patch.set_facecolor('black')

				ax1 = plt.subplot(1,3,1)
				ax2 = plt.subplot(1,3,2)
				ax3 = plt.subplot(1,3,3)
				fig.axes[panes.index(1)].set_title(ID + "\n", size=12, color='w', family='monospace')
				fig.axes[panes.index(1)].set_xlabel("\nTISSUE: {:.3f}, CSF: {:.3f}".format(tissue_frac, vox_frac[2]), size=11, color='w', family='monospace')

				palette = cm.Greys_r
				palette.set_over('g', 0.6)

				asp = brain.header['pixdim'][1]/brain.header['pixdim'][3]
				print('asp', asp)
				# raw_input()

				print('brain_img', brain_img[vox_centroid[0],:,:].squeeze())
				print('vox_centroid', vox_centroid)
				# raw_input()

				anat_zoom1, k1 = autozoom(np.transpose(brain_img[vox_centroid[0],:,:].squeeze(), transpose_indices))
				# raw_input()

				anat_zoom2, k2 = autozoom(np.transpose(brain_img[:,vox_centroid[1],:].squeeze(), transpose_indices))
				# raw_input()

				anat_zoom3, k3 = autozoom(np.transpose(brain_img[:,:,vox_centroid[2]].squeeze(), transpose_indices))
				# raw_input()
				
				ax1.imshow(np.rot90(anat_zoom1, n_rot90[0]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax2.imshow(np.rot90(anat_zoom2, n_rot90[1]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax3.imshow(np.rot90(anat_zoom3, n_rot90[2]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect='equal')

				ax1.imshow(np.rot90(np.transpose(vox_mas[vox_centroid[0],min(k1[0]):max(k1[0])+1,min(k1[1]):max(k1[1])+1].squeeze(),transpose_indices), n_rot90[0]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax2.imshow(np.rot90(np.transpose(vox_mas[min(k2[0]):max(k2[0])+1,vox_centroid[1],min(k2[1]):max(k2[1])+1].squeeze(),transpose_indices), n_rot90[1]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect=asp)
				ax3.imshow(np.rot90(np.transpose(vox_mas[min(k3[0]):max(k3[0])+1,min(k3[1]):max(k3[1])+1,vox_centroid[2]].squeeze(),transpose_indices), n_rot90[2]), cmap = palette, norm = colors.Normalize(vmin = brain_img.min() - 1, vmax = brain_img.max() + 1, clip = False), aspect='equal')

				plt.setp([a.set_xticks([]) for a in fig.axes])
				plt.setp([a.set_yticks([]) for a in fig.axes])
				plt.savefig(mouse + "/barstool_output.png", dpi=300, facecolor='k', bbox_inches='tight', pad_inches = 0.2)

				# raw_input()

				for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
					out_file.write("{:6.6f},".format(f_conc[str(self.metabParamsTableWidget.item(metab_index,0).text())]))
				for metab_index in range(0,self.metabParamsTableWidget.rowCount()):
					out_file.write("{:6.6f},".format(f_crlb[str(self.metabParamsTableWidget.item(metab_index,0).text())]))
			except Exception as e:
				print(e)
				failed_mice.append(mouse)			
			out_file.write('\n')
		
		out_file.close()

		self.consoleOutputText.append('')
		self.consoleOutputText.append('The following files could not be processed successfully:')
		for mouse in failed_mice:
			self.consoleOutputText.append(' >> ' + str(mouse))

		self.loadOutputsButton_quant_bruker.setEnabled(True)
		self.confirmSaveFileButton_quant_bruker.setEnabled(False)
		self.runQuantButton_bruker.setEnabled(False)

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())
