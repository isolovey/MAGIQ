# ---- System Libraries ---- #
from builtins import zip
from builtins import range

# Check for PyQt6 (for native mac M1 compatibility), otherwise continue using PyQt5
import importlib.util
PyQt6_spec = importlib.util.find_spec("PyQt6")
if PyQt6_spec != None:
	from PyQt6 import QtCore
else:
	from PyQt5 import QtCore

# ---- Math Libraries ---- #
import scipy as sp
import numpy as np
import math

# ---- Simulation Libraries ---- #
import pygamma as pg

# ---- Data Classes ---- #
from magiqdataclasses import *

# ---- Simulation Class ---- #
class MetaboliteSimulation(QtCore.QObject):

	postToConsole = QtCore.pyqtSignal(str)
	outputResults = QtCore.pyqtSignal(object)
	finished = QtCore.pyqtSignal(int)

	def __init__(self, thread_num, insysfile, sim_experiment):
		
		QtCore.QObject.__init__(self)

		self.thread_num = thread_num
		self.insysfile = insysfile
		self.sim_experiment = sim_experiment

	def simulate(self):
		self.postToConsole.emit('   | Simulating ... ' + self.insysfile)
		print('    | Simulating ...' + self.insysfile)

		metab_name = self.insysfile.replace('.sys','')

		if self.sim_experiment.b0 == 123.3:
			self.insysfile = 'pints/metabolites/3T_' + self.insysfile
		elif self.sim_experiment.b0 == 297.2:
			self.insysfile = 'pints/metabolites/7T_' + self.insysfile
		elif self.sim_experiment.b0 == 400.2:
			self.insysfile = 'pints/metabolites/9.4T_' + self.insysfile

		if self.sim_experiment.name == "semi-LASER (Bruker)":
			spin_system = pg.spin_system()
			spin_system.read(self.insysfile)
			for i in range(spin_system.spins()):
				spin_system.PPM(i, spin_system.PPM(i) - self.sim_experiment.RF_OFFSET)

			TE = self.sim_experiment.TE * 1E-3
			TE1 = self.sim_experiment.TE1 * 1E-3
			TE2 = self.sim_experiment.TE2 * 1E-3

			# build 90 degree pulse
			inpulse90file = self.sim_experiment.inpulse90file
			A_90 = self.sim_experiment.A_90
			PULSE_90_LENGTH = self.sim_experiment.PULSE_90_LENGTH
			gyratio = self.sim_experiment.getGyratio()

			pulse90 = Pulse(inpulse90file, PULSE_90_LENGTH, 'bruker')

			n_old = np.linspace(0, PULSE_90_LENGTH, sp.size(pulse90.waveform))
			n_new = np.linspace(0, PULSE_90_LENGTH, sp.size(pulse90.waveform)+1)

			waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse90.waveform)*A_90)(n_new)
			waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse90.waveform)*A_90)(n_new)
			pulse90.waveform = waveform_real + 1j*(waveform_imag)

			ampl_arr = np.abs(pulse90.waveform)
			phas_arr = np.unwrap(np.angle(pulse90.waveform))*180.0/math.pi
			
			pulse = pg.row_vector(len(pulse90.waveform))
			ptime = pg.row_vector(len(pulse90.waveform))
			for j, val in enumerate(zip(ampl_arr, phas_arr)):
				pulse.put(pg.complex(val[0],val[1]), j)
				ptime.put(pg.complex(pulse90.pulsestep,0), j)

			pulse_dur_90 = pulse.size() * pulse90.pulsestep
			pwf_90 = pg.PulWaveform(pulse, ptime, "90excite")
			pulc_90 = pg.PulComposite(pwf_90, spin_system, self.sim_experiment.obs_iso)

			Ureal90 = pulc_90.GetUsum(-1)

			# build 180 degree pulse
			inpulse180file = self.sim_experiment.inpulse180file
			A_180 = self.sim_experiment.A_180
			PULSE_180_LENGTH = self.sim_experiment.PULSE_180_LENGTH
			gyratio = self.sim_experiment.getGyratio()

			pulse180 = Pulse(inpulse180file, PULSE_180_LENGTH, 'bruker')

			n_old = np.linspace(0, PULSE_180_LENGTH, sp.size(pulse180.waveform))
			n_new = np.linspace(0, PULSE_180_LENGTH, sp.size(pulse180.waveform)+1)

			waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse180.waveform)*A_180)(n_new)
			waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse180.waveform)*A_180)(n_new)
			pulse180.waveform = waveform_real + 1j*(waveform_imag)

			ampl_arr = np.abs(pulse180.waveform)
			phas_arr = np.unwrap(np.angle(pulse180.waveform))*180.0/math.pi
			freq_arr = np.gradient(phas_arr)

			pulse = pg.row_vector(len(pulse180.waveform))
			ptime = pg.row_vector(len(pulse180.waveform))
			for j, val in enumerate(zip(ampl_arr, phas_arr)):
				pulse.put(pg.complex(val[0],val[1]), j)
				ptime.put(pg.complex(n_new[1],0), j)

			pulse_dur_180 = pulse.size() * pulse180.pulsestep
			pwf_180 = pg.PulWaveform(pulse, ptime, "180afp")
			pulc_180 = pg.PulComposite(pwf_180, spin_system, self.sim_experiment.obs_iso)

			Ureal180 = pulc_180.GetUsum(-1)

			H = pg.Hcs(spin_system) + pg.HJ(spin_system)
			D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
			ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
			ACQ = ac

			delay1 = TE1/2.0 - pulse_dur_90/2.0 - pulse_dur_180/2.0
			delay2 = TE1/2.0 + TE2/2.0 - pulse_dur_180
			delay3 = TE2 - pulse_dur_180
			delay4 = delay2
			delay5 = TE1/2.0 - pulse_dur_180 + self.sim_experiment.DigShift

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

		elif self.sim_experiment.name == "semi-LASER":
			spin_system = pg.spin_system()
			spin_system.read(self.insysfile)
			for i in range(spin_system.spins()):
				spin_system.PPM(i, spin_system.PPM(i) - self.sim_experiment.RF_OFFSET)

			TE = self.sim_experiment.TE
			TE1 = float((TE * 0.31) / 1000.0)
			TE3 = float((TE * 0.31) / 1000.0)
			TE2 = float(TE/1000.0 - TE1 - TE3)
			TE_fill = TE/1000.0 - TE1 - TE2 - TE3

			# build 90 degree pulse
			inpulse90file = self.sim_experiment.inpulse90file
			A_90 = self.sim_experiment.A_90
			PULSE_90_LENGTH = self.sim_experiment.PULSE_90_LENGTH
			gyratio = self.sim_experiment.getGyratio()

			pulse90 = Pulse(inpulse90file, PULSE_90_LENGTH)

			n_old = np.linspace(0, PULSE_90_LENGTH, sp.size(pulse90.waveform))
			n_new = np.linspace(0, PULSE_90_LENGTH, sp.size(pulse90.waveform)+1)

			waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse90.waveform)*A_90)(n_new)
			waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse90.waveform)*A_90)(n_new)
			pulse90.waveform = waveform_real + 1j*(waveform_imag)

			ampl_arr = np.abs(pulse90.waveform)*gyratio
			phas_arr = np.unwrap(np.angle(pulse90.waveform))*180.0/math.pi
			
			pulse = pg.row_vector(len(pulse90.waveform))
			ptime = pg.row_vector(len(pulse90.waveform))
			for j, val in enumerate(zip(ampl_arr, phas_arr)):
				pulse.put(pg.complex(val[0],val[1]), j)
				ptime.put(pg.complex(pulse90.pulsestep,0), j)

			pulse_dur_90 = pulse.size() * pulse90.pulsestep
			peak_to_end_90 = pulse_dur_90 - (209 + self.sim_experiment.fudge_factor) * pulse90.pulsestep
			pwf_90 = pg.PulWaveform(pulse, ptime, "90excite")
			pulc_90 = pg.PulComposite(pwf_90, spin_system, self.sim_experiment.obs_iso)

			Ureal90 = pulc_90.GetUsum(-1)

			# build 180 degree pulse
			inpulse180file = self.sim_experiment.inpulse180file
			A_180 = self.sim_experiment.A_180
			PULSE_180_LENGTH = self.sim_experiment.PULSE_180_LENGTH
			gyratio = self.sim_experiment.getGyratio()

			pulse180 = Pulse(inpulse180file, PULSE_180_LENGTH)

			n_old = np.linspace(0, PULSE_180_LENGTH, sp.size(pulse180.waveform))
			n_new = np.linspace(0, PULSE_180_LENGTH, sp.size(pulse180.waveform)+1)

			waveform_real = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.real(pulse180.waveform)*A_180)(n_new)
			waveform_imag = sp.interpolate.InterpolatedUnivariateSpline(n_old, np.imag(pulse180.waveform)*A_180)(n_new)
			pulse180.waveform = waveform_real + 1j*(waveform_imag)

			ampl_arr = np.abs(pulse180.waveform)*gyratio
			phas_arr = np.unwrap(np.angle(pulse180.waveform))*180.0/math.pi
			freq_arr = np.gradient(phas_arr)

			pulse = pg.row_vector(len(pulse180.waveform))
			ptime = pg.row_vector(len(pulse180.waveform))
			for j, val in enumerate(zip(ampl_arr, phas_arr)):
				pulse.put(pg.complex(val[0],val[1]), j)
				ptime.put(pg.complex(n_new[1],0), j)

			pulse_dur_180 = pulse.size() * pulse180.pulsestep
			pwf_180 = pg.PulWaveform(pulse, ptime, "180afp")
			pulc_180 = pg.PulComposite(pwf_180, spin_system, self.sim_experiment.obs_iso)

			Ureal180 = pulc_180.GetUsum(-1)

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

		elif self.sim_experiment.name == "LASER":
			spin_system = pg.spin_system()
			spin_system.read(self.insysfile)
			for i in range(spin_system.spins()):
				spin_system.PPM(i, spin_system.PPM(i) - self.sim_experiment.RF_OFFSET)

			# build 90 degree AHP pulse
			inpulse90file = self.sim_experiment.inpulse90file
			A_90 = self.sim_experiment.A_90
			PULSE_90_LENGTH = self.sim_experiment.PULSE_90_LENGTH
			gyratio = self.sim_experiment.getGyratio()

			pulse90 = Pulse(inpulse90file, PULSE_90_LENGTH, 'varian')

			n_new = np.linspace(0, PULSE_90_LENGTH, 256)

			waveform_real = np.real(pulse90.waveform)*A_90
			waveform_imag = np.imag(pulse90.waveform)*A_90
			pulse90.waveform = waveform_real + 1j*(waveform_imag)

			ampl_arr = np.abs(pulse90.waveform)*gyratio
			phas_arr = np.unwrap(np.angle(pulse90.waveform))*180.0/math.pi
			
			pulse = pg.row_vector(len(pulse90.waveform))
			ptime = pg.row_vector(len(pulse90.waveform))
			for j, val in enumerate(zip(ampl_arr, phas_arr)):
				pulse.put(pg.complex(val[0],val[1]), j)
				ptime.put(pg.complex(pulse90.pulsestep,0), j)

			pulse_dur_90 = pulse.size() * pulse90.pulsestep
			pwf_90 = pg.PulWaveform(pulse, ptime, "90excite")
			pulc_90 = pg.PulComposite(pwf_90, spin_system, self.sim_experiment.obs_iso)

			Ureal90 = pulc_90.GetUsum(-1)

			# build 180 degree pulse
			inpulse180file = self.sim_experiment.inpulse180file
			A_180 = self.sim_experiment.A_180
			PULSE_180_LENGTH = self.sim_experiment.PULSE_180_LENGTH
			gyratio = self.sim_experiment.getGyratio()

			pulse180 = Pulse(inpulse180file, PULSE_180_LENGTH, 'varian')

			n_new = np.linspace(0, PULSE_180_LENGTH, 512)

			waveform_real = np.real(pulse180.waveform)*A_180
			waveform_imag = np.imag(pulse180.waveform)*A_180
			pulse180.waveform = waveform_real + 1j*(waveform_imag)

			ampl_arr = np.abs(pulse180.waveform)*gyratio
			phas_arr = np.unwrap(np.angle(pulse180.waveform))*180.0/math.pi
			freq_arr = np.gradient(phas_arr)

			pulse = pg.row_vector(len(pulse180.waveform))
			ptime = pg.row_vector(len(pulse180.waveform))
			for j, val in enumerate(zip(ampl_arr, phas_arr)):
				pulse.put(pg.complex(val[0],val[1]), j)
				ptime.put(pg.complex(n_new[1],0), j)

			pulse_dur_180 = pulse.size() * pulse180.pulsestep
			pwf_180 = pg.PulWaveform(pulse, ptime, "180afp")
			pulc_180 = pg.PulComposite(pwf_180, spin_system, self.sim_experiment.obs_iso)

			Ureal180 = pulc_180.GetUsum(-1)

			# calculate pulse timings
			ROF1 = 100E-6    #sec
			ROF2 = 10E-6     #sec
			TCRUSH1 = 0.0008 #sec
			TCRUSH2 = 0.0008 #sec

			ss_grad_rfDelayFront = 0#TCRUSH1 - ROF1
			ss_grad_rfDelayBack  = 0#TCRUSH2 - ROF2
			ro_grad_atDelayFront = 0
			ro_grad_atDelayBack  = 0

			TE  = self.sim_experiment.TE / 1000.
			ipd = (TE - pulse_dur_90 \
					  - 6*(ss_grad_rfDelayFront + pulse_dur_180 + ss_grad_rfDelayBack) \
					  - ro_grad_atDelayFront) / 12

			delay1 = ipd + ss_grad_rfDelayFront
			delay2 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
			delay3 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
			delay4 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
			delay5 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
			delay6 = ss_grad_rfDelayBack + 2*ipd + ss_grad_rfDelayFront
			delay7 = ss_grad_rfDelayBack + ipd + ro_grad_atDelayFront

			# print A_90, A_180, pulse_dur_90, pulse_dur_180
			# print TE, ipd, pulse_dur_90+6*pulse_dur_180, delay1+delay2+delay3+delay4+delay5+delay6+delay7, pulse_dur_90+6*pulse_dur_180+delay1+delay2+delay3+delay4+delay5+delay6+delay7
			# print ''

			# initialize acquisition
			H = pg.Hcs(spin_system) + pg.HJ(spin_system)
			D = pg.Fm(spin_system, self.sim_experiment.obs_iso)
			ac = pg.acquire1D(pg.gen_op(D), H, self.sim_experiment.dwell_time)
			ACQ = ac

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

		b0 = self.sim_experiment.b0
		obs_iso = self.sim_experiment.obs_iso
		tolppm = self.sim_experiment.tolppm
		tolpha = self.sim_experiment.tolpha
		ppmlo = self.sim_experiment.ppmlo
		ppmhi = self.sim_experiment.ppmhi
		rf_off = self.sim_experiment.RF_OFFSET

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

		metab = Metabolite()
		metab.name = metab_name
		metab.var = 0.0

		for i in range(sp.size(outf)):
			if outf[i] <= 5:
				metab.ppm.append(outf[i])
				metab.area.append(outa[i])
				metab.phase.append(-1.0*outp[i])

		insysfile = self.insysfile.replace('pints/metabolites/3T_', '')
		insysfile = self.insysfile.replace('pints/metabolites/7T_', '')
		insysfile = self.insysfile.replace('pints/metabolites/9.4T_','')

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

		# Send save data signal
		self.outputResults.emit(metab)
		self.postToConsole.emit('        | Simulation completed for ... ' + self.insysfile)
		self.finished.emit(self.thread_num)
