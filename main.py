from PyQt5 import QtCore, QtGui, QtWidgets, uic
import subprocess
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

	def launchPINTS(self):
		print 'Launching PINTS'
		command = 'python pints.py'
		print command
		t = Thread(target = lambda: os.system(command))
		t.start()
		

	def launchFITMAN(self):
		print 'Launching FITMAN'
		command = 'python fitman.py'
		print command
		t = Thread(target = lambda: os.system(command))
		t.start()

	def launchSPICeS(self):
		print 'Launching SPICeS'
		command = 'python spices.py'
		print command
		t = Thread(target = lambda: os.system(command))
		t.start()

	def launchBARSTOOL(self):
		print 'Launching BARSTOOL'
		command = 'python barstool.py'
		print command
		t = Thread(target = lambda: os.system(command))
		t.start()

# ---- Launch Application ---- #
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec_())