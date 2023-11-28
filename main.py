from __future__ import print_function

# Check for PyQt6 (for native mac M1 compatibility), otherwise continue using PyQt5
import importlib.util
PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore, QtGui, QtWidgets, uic
else:
	from PyQt5 import QtCore, QtGui, QtWidgets, uic

import os
import sys
from threading import Thread

qtCreatorFile = "MAGIQ.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

# ---- Main Application Object ---- #
class MyApp(QtWidgets.QWidget, Ui_MainWindow):
	def __init__(self):
		QtWidgets.QWidget.__init__(self)
		Ui_MainWindow.__init__(self)
		self.setupUi(self)

		# Bind buttons
		self.pintsButton.clicked.connect(self.launchPINTS)
		self.fitmanButton.clicked.connect(self.launchFITMAN)
		self.spicesButton.clicked.connect(self.launchSPICeS)
		self.barstoolButton.clicked.connect(self.launchBARSTOOL)
		self.barstoolrvButton.clicked.connect(self.launchBASRTOOLRV)
		self.appsButton.clicked.connect(self.launchAPPS)

	def launchPINTS(self):
		print('Launching PINTS')
		command = 'python pints.py'
		print(command)
		t = Thread(target = lambda: os.system(command))
		t.start()
		

	def launchFITMAN(self):
		print('Launching FITMAN')
		command = 'python fitman.py'
		print(command)
		t = Thread(target = lambda: os.system(command))
		t.start()

	def launchSPICeS(self):
		print('Launching SPICeS')
		command = 'python spices.py'
		print(command)
		t = Thread(target = lambda: os.system(command))
		t.start()

	def launchBARSTOOL(self):
		print('Launching BARSTOOL')
		command = 'python barstool.py'
		print(command)
		t = Thread(target = lambda: os.system(command))
		t.start()

	def launchBASRTOOLRV(self):
		print('Launching BARSTOOL (Rodent Version)')
		command = 'python barstoolrv.py'
		print(command)
		t = Thread(target = lambda: os.system(command))
		t.start()

	def launchAPPS(self):
		print('Launch APPS')
		command = 'python apps.py'
		print(command)
		t = Thread(target = lambda: os.system(command))
		t.start()

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec())