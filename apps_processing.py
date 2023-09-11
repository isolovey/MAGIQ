import matplotlib.pyplot as plt
import numpy as np

import copy
from typing import Literal, Optional, Tuple, Union

from PyQt5.QtWidgets import QTextBrowser

from magiqdataclasses import BrukerFID, DatFile, Metabolite
from preproc import baseline_corr, quecc, quality, ecc

class AppsProcessing(object):

    @staticmethod
    def run_bruker_conversion(out_name: str,
                              sup_file_path: str,
                              ref_file_path: str,
                              baseline_correction: bool = True,
                              post_processing: Literal["quecc", "quality", "ecc", None] = "quecc",
                              quality_points: int = 200,
                              save_plot: Optional[str] = None):
        '''
            This method runs the conversion process.
        '''
        out_name_sup = out_name + 'corr_sup'
        out_name_uns = out_name + 'corr_uns'

        # Read suppressed file
        sup_file = BrukerFID(sup_file_path)
        # Read unsuppressed file
        uns_file = BrukerFID(ref_file_path)

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
        if baseline_correction:
            sup_file.signal = baseline_corr(sup_file.signal)
            uns_file.signal = baseline_corr(uns_file.signal)
            out_name_sup = out_name_sup + '_bc'
            out_name_uns = out_name_uns + '_bc'

        # Run post-processing.
        out_file_sup = copy.deepcopy(sup_file)
        out_file_uns = copy.deepcopy(uns_file)
        if post_processing == "quecc":
            # QUECC
            quecc_points = int(quality_points)
            out_file_sup.signal = np.array(quecc(sup_file.signal, uns_file.signal, quecc_points, sup_file.t))
            out_file_uns.signal = np.array(quecc(uns_file.signal, uns_file.signal, quecc_points, uns_file.t))
            out_file_sup.writeDAT(out_name_sup, 'quecc' + str(quecc_points))
            out_file_uns.writeDAT(out_name_uns, 'quecc' + str(quecc_points))
        elif post_processing == "quality":
            # QUALITY
            out_file_sup.signal = np.array(quality(sup_file.signal, uns_file.signal))
            out_file_uns.signal = np.array(quality(uns_file.signal, uns_file.signal))
            out_file_sup.writeDAT(out_name_sup, 'quality')
            out_file_uns.writeDAT(out_name_uns, 'quality')
        elif post_processing == "ecc":
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

        if save_plot:
            plt.savefig(save_plot, dpi=300, bbox_inches='tight')
            plt.close()

        return sup_file, uns_file



    # ---- Methods for HSVD Fitting ---- #
    @staticmethod
    def lorentzian(time_axis, frequency, phase, fwhm):
        '''
            This method returns a time-domain lorentzian function.
                Inputs:
                    time_axis    discrete time array
                    frequency    frequency of the Lorentzian function
                    phase        phase of the Lorentzian function
                    fwhm        Lorentzian linewidth

                Outputs:
                    time-domian lorentzian function
        '''
        oscillatory_term = np.exp(1j * (2 * np.pi * frequency * time_axis + phase))
        damping = np.exp(-time_axis * np.pi * fwhm)
        fid = oscillatory_term * damping
        return fid / len(time_axis)


    @staticmethod
    def hsvd(dat, n, ratio, comp, console):
        '''
            This method's code is based (heavily) on: https://github.com/openmrslab/suspect/blob/master/suspect/processing/water_suppression.py
            It performs an HSVD decomposition and returns a list of all parameters for each Lorentzian component.

                Inputs:
                    dat        DatFile object containing the signal on which to perform the decomposition
                    n        number of points
                    ratio    Hankel matrix row/col ratio
                    comp    number of signal related singular values

                Outputs:
                    An array containing the following elements:
                    [0] array of component numbers
                    [1] array of the damping coefficient of each component
                    [2] array of the frequency coeffienct of each component
                    [3] array of the amplitude of each component
                    [4] array of the phase of each component
        '''


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
            X[:, i] = AppsProcessing.lorentzian(dat.t,
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

    @staticmethod
    def run_water_removal(input_dat: DatFile,
                          hsvd_points: int = 512,
                          hsvd_ratio: float = 1.25,
                          hsvd_components: int = 35,
                          frequency_range_xmin: float = -1.5,
                          frequency_range_xmax: float = 0.5,
                          console: Optional[Union[QTextBrowser, list]] = None,
                          save_plot: Optional[str] = None,
                          ) -> Tuple[DatFile, DatFile]:
        '''
            This method performs removal of the residual water signal
            using Hankel Singular Value decomposition.
        '''
        if not console:
            console = list()
        # 1. Fit specturm with HSVD.
        peak, width_L, ppm, area, phase = AppsProcessing.hsvd(
            input_dat,
            hsvd_points,
            hsvd_ratio,
            hsvd_components,
            console
        )
        hsvd_fit = Metabolite()

        # take only the peaks within the specified frequency range
        for (i, comp) in enumerate(peak):
            if (ppm[i] / input_dat.b0) >= frequency_range_xmin and (ppm[i] / input_dat.b0) <= frequency_range_xmax:
                hsvd_fit.peak.append(peak[i])
                hsvd_fit.width_L.append(width_L[i])
                hsvd_fit.ppm.append(ppm[i] / input_dat.b0)
                hsvd_fit.area.append(area[i])
                hsvd_fit.phase.append(phase[i])

        hsvd_fid = hsvd_fit.getFID(input_dat.TE, input_dat.b0, input_dat.t, 0, 1, 0, 0, 0)

        dat_hsvd = copy.deepcopy(input_dat)
        dat_hsvd.signal = hsvd_fid

        ax = plt.subplot(211)
        ax.clear()
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().set_visible(False)

        # plot raw spectrum
        f_dat, spec_dat = input_dat.getSpec()
        plt.plot(f_dat[0:input_dat.n], np.real(spec_dat[0:input_dat.n]), label='Data')

        # plot fitted spectrum
        f_hsvd, spec_hsvd = dat_hsvd.getSpec()
        plt.plot(f_hsvd, np.real(spec_hsvd), label='Fit')

        # plot residual
        VSHIFT = 0.25
        plt.plot(f_hsvd, np.real(spec_dat[0:input_dat.n]) - np.real(spec_hsvd) + VSHIFT * np.amax(np.real(spec_hsvd)),
                 label='Residual')

        # legend and title
        ax.legend(loc="upper left")
        plt.title('Original Spectrum')

        # 2. Subtract from original spectrum.
        scale = 1.0
        dat_wr = copy.deepcopy(input_dat)
        dat_wr.n = np.min([np.size(input_dat.signal, 0), np.size(dat_hsvd.signal, 0)])
        dat_wr.signal = input_dat.signal[0:dat_wr.n] - scale * dat_hsvd.signal[0:dat_wr.n]

        ax = plt.subplot(212)
        ax.clear()
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().set_visible(False)

        f_wr, spec_wr = dat_wr.getSpec()
        plt.plot(f_wr, np.real(spec_wr))
        plt.title('Water Removed Spectrum')

        plt.tight_layout()

        if save_plot:
            plt.savefig(save_plot, dpi=300, bbox_inches='tight')
            plt.close()

        return dat_hsvd, dat_wr

    @staticmethod
    def save_water_removal(dat, dat_wr):
        '''
            This method saves the water removed signal as a *.dat file.
        '''

        dat_wr.filename = dat.filename.replace('.dat', '_wr.dat')

        out_file = open(dat_wr.filename, 'w')
        in_file  = open(dat.filename, 'r')

        for (i, line) in enumerate(in_file):
            if i > 11:
                for element in dat_wr.signal:
                    out_file.write("{0:.6f}".format(float(np.real(element))) + '\n')
                    out_file.write("{0:.6f}".format(float(np.imag(element))) + '\n')
                break
            else:
                out_file.write(line)

        out_file.close()
        in_file.close()