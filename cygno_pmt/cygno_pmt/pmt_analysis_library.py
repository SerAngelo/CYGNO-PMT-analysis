"""
@author: SerAngelo
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.signal import resample, butter, filtfilt, find_peaks


class PMT:
    def __init__(self, name, model, RCconstant, channel):
        self.name = name
        self.model = model
        self.RCconstant = RCconstant
        self.channel = channel
        self.events = []
        

class PMT_event:
    
    def __init__(self, PMT, serial_number, waveformf, digitizer_type, event_type, threshold):
        self.PMT = PMT
        self.serial_number = serial_number
        
        if digitizer_type == 'fast':
            conv_factor = 1000.0/(2.0**12)
            self.waveform = np.asarray(waveformf)*conv_factor
        elif digitizer_type == 'slow':
            conv_factor = 2000.0/(2.0**12)
            self.waveform = np.asarray(waveformf)*conv_factor
        
        self.digitizer_type = digitizer_type
        self.event_type = event_type
        self.threshold = threshold
        
        
    def pmt_max_ampl(self):
        return np.argmin(self.waveform)
                
        
    def pmt_peaks(self, t1, t2, dist):
                
        mean, _, _ = self.background()
        
        peaks, _= find_peaks(abs(self.waveform[t1:t2]-mean), height=self.threshold, distance=dist)
        n_peaks = len(peaks)
        
        return [peaks, n_peaks]  
    
    
    def background(self):
        
        t_peak = self.pmt_max_ampl()
          
        fondo = np.concatenate((self.waveform[0:t_peak-100], self.waveform[t_peak+100:]))
            
        mean = np.mean(fondo)
        #print('mean=',mean)
        std = np.std(fondo)
        #print('std=',std)
        return [mean,std, fondo]
    
    
    def __static_find_boundaries(self):

        try:
            peaks, n_peaks = self.pmt_peaks(0, len(self.waveform), 5)
            tleft = peaks[0]
            tright = peaks[-1]
        except:
            t_peak = self.pmt_max_ampl()
            tleft = t_peak
            tright = t_peak

        mean, std, _ = self.background()
        
        eps=0.1
        
        tmin, tmax = 0, len(self.waveform)-1
        
            
            
        for t in reversed(range(0,tleft)):
            if mean-self.waveform[t]-std < eps:
                tmin=t
                break
                

        for t in range(tright,len(self.waveform)):
            if mean-self.waveform[t]-std < eps:
                tmax=t
                break
                
        return [tmin,tmax]
    
    
    
    def __dynamic_find_boundaries(self):
        
        t_peak = self.pmt_max_ampl() 
          
        mean, std, _ = self.background()
        
        eps=0.1
        
        tmin, tmax = 0, len(self.waveform)
        
        while(True):
            
            for t in reversed(range(0,t_peak)):
                if mean-self.waveform[t]-std < eps:
                    tmin=t
                    break
                    
        
            for t in range(t_peak,len(self.waveform)):
                if mean-self.waveform[t]-std < eps:
                    tmax=t
                    break
        
            _, n_peaks = self.pmt_peaks(tmin, tmax, 20)
        
            eps+=1
        
            if n_peaks == 1 or eps>50:
                break
            
        return [tmin,tmax, mean]
    
    

    def gaussian_noise(self, axi):
                
        #determino info sul fondo
        mean, std, fondo = self.background()
        
        n, bins, _ = axi.hist(fondo-mean, bins=8, orientation='horizontal', edgecolor='black', color='gray', alpha=0.6)
        
        
        y = np.linspace(axi.get_ylim()[0], axi.get_ylim()[1], 300)
        
        p = norm.pdf(y, 0, std)
        #print(axi.get_ylim()[0], axi.get_ylim()[1], max(n), max(p))
        p = p * max(n)/max(p)
        
    
        
        axi.plot(p, y, 'k-', linewidth=2)
        
        axi.get_xaxis().set_visible(False)
        axi.get_yaxis().set_visible(True)
        axi.set_title('Statistical noise', fontsize=20)
        
        axi.text(0.9, 0.8, f'$\sigma$={std:.2f} mV', transform=axi.transAxes, fontsize=20,
            verticalalignment='top', horizontalalignment='right', color='black',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', pad=0.5))

            
    
    def PMT_waveforms_integral_info(self, resampling=False, fit=False, verbose=False, axi=None):
        
        if self.digitizer_type == "fast":
            sample_points = 1024
            frequence = 1333e6
            DT = sample_points/frequence
            
            #definisco l'asse temporale in base al tipo di digitizer
            t = np.linspace(0, DT, sample_points)
            
            if resampling:
                sample_points = 512
        
            
        elif self.digitizer_type == "slow":
            sample_points = 4000
            frequence = 250e6 #S/s
            DT = sample_points/frequence
            
            #definisco l'asse temporale in base al tipo di digitizer
            t = np.linspace(0, DT, sample_points)
            
            if resampling:
                sample_points = 2000
                
    
        if resampling:
            #Eseguo il resampling del segnale
            nyquist_freq = 0.5 * sample_points / (t[-1]*frequence - t[0]*frequence)
            b, a = butter(N=5, Wn=nyquist_freq, btype='low', analog=False, output='ba')            
            self.waveform, t = resample(filtfilt(b, a, self.waveform), sample_points, t)

    
            
        #determino la posizione del picco più grande
        amplitude = np.abs(np.min(self.waveform))
    
    
        #determino info sul fondo
        mean, std, fondo = self.background()
        
        #determino l'intervallo di integrazione
        l, r = self.__static_find_boundaries()
        
        
        #determino il numero di picchi nell'intervallo
        peaks, npeaks = self.pmt_peaks(l, r, 10)
        #print("#############", peaks+l, t_peak )
        
        
        
          
        #determino gli array di fondo (necessari per la visualizzazione dell'area)
        fondo = np.array([mean]*sample_points)
        
        intervallo = (t >= t[l]) & (t <= t[r])
        area = np.trapz(np.array(self.waveform)[intervallo]-mean, t[intervallo])
        integrale = (abs(area)*1e12)/(1000*50) #50 Ohm impedenza
        
        intervallo_noise = (t < t[l]) | (t > t[r])
        area_noise = np.trapz(np.array(self.waveform)[intervallo_noise]-mean, t[intervallo_noise])
        integrale_noise = (abs(area_noise)*1e12)/(1000*50) #50 Ohm impedenza
        
   
        
        if fit:
            
            if l-50 > 0:
                xmin_fit = l-50
            else:
                xmin_fit = 0
                
            if l+50 < np.size(t):
                xmax_fit = l+50
            else:
                xmax_fit = -1
               
            try:
                popt, pcov = curve_fit(signal_model, t[xmin_fit:xmax_fit], self.waveform[xmin_fit:xmax_fit]-mean, p0 = [np.min(self.waveform)-mean, 0, t[peaks[0]+l], t[r]-t[l]])
                std_fit = np.sqrt(np.diag(pcov))
                FWHM, half, points = fwhm(t[xmin_fit:xmax_fit], signal_model(t[xmin_fit:xmax_fit], *popt))
                
                
                for p, std in zip(popt, std_fit):
                    if std > abs(p) * 100:  
                        print("Fit non convergente!")
                        fit = False
                        break
                    
            except:
                print("Fit non convergente!")
                fit=False
        
            
                            
    
        
        if verbose:
          
            print("="*50)
            print(f"{self.PMT.name} - event: {self.serial_number} - threshold: {self.threshold}mV")
            print("="*50)
            print(f"Baseline: ({mean:.2f} ± {std:.2f}) mV")
            print(f"Number of peaks: {npeaks}")
            print(f"Released charge: {integrale:.2f} pC")    
            print(f"Peak width: {t[r]-t[l]:.2e} s")
            print(f"Noise charge: {integrale_noise:.2f} pC")    
            
            if fit:
                param_names = ["a", "b", "c", "tau"]
                # Estrai le deviazioni standard dai valori diagonali della matrice di covarianza
                std_fit = np.sqrt(np.diag(pcov))
                
                # Arrotonda i parametri e le deviazioni standard a 2 cifre decimali
                parameters = [f"{param:.2e}" for param in popt]
                std_fit_rounded = [f"{std:.2e}" for std in std_fit]
                dim = ["mV", "mV", "s", "s"]
                # Associa i parametri e le loro deviazioni standard
                params_with_std = [f"{name} = ({param} ± {std}) {d}" for name, param, std, d in zip(param_names, parameters, std_fit_rounded, dim)]
                
                print("="*50)
                print("Fit parameters and FWHM:")
                print("-" * 50)
                for param in params_with_std:
                    print(f"  - {param}")
                print(f"  - FWHM: {FWHM:.2e} s")
                print("-" * 50)
                
            print("\n\n")
        
        
        
        if axi != None:
            
            try:
                xmin = t[l-50]
            except:
                xmin=t[0]
            
            try:
                xmax = t[r+50]
            except:
                xmax = t[-1]
            
            axi.set_title(f"{self.PMT.name} - event: {self.serial_number} - threshold: {self.threshold}mV", fontsize=20)
            axi.plot(t, self.waveform-mean, marker='o', label='Raw signal', alpha=0.5)
    
            axi.plot(t[peaks+l], np.array(self.waveform)[peaks+l]-mean, 'x', markersize = 20, color='orange', linewidth=20, label=f'{npeaks} peaks above {self.threshold}mV')
    
            #linee rosse
            axi.axvline(x=t[l], color='red', linestyle='--') 
            axi.axvline(x=t[r], color='red', linestyle='--')   
            
            axi.fill_between(t, fondo-mean, self.waveform-mean, where=intervallo, color='grey', alpha=0.2, label=f'Charge: {integrale:.2f} $pC$')
            
            
            if fit:
                axi.plot(t[l-10:r+10], signal_model(t[l-10:r+10], *popt), label='Fit of the signal', color='red', linewidth=2)    
                
                if FWHM != 0:
                    axi.plot([t[points[0]+xmin_fit], t[points[-1]+xmin_fit]], [half, half], 'g-', lw=2, label=f'FWHM={FWHM:.2e}')
                
               
                    
            axi.set_xlim(xmin,xmax)
            axi.grid(True, linestyle='--', linewidth=0.5)
            axi.set_ylabel('Signal Value [mV]', fontsize=20)
            axi.set_xlabel('Time (s)', fontsize=20)
            axi.legend(fontsize=15, loc='lower right')
        
        
        return [integrale, integrale/(t[r]-t[l]), t[r]-t[l], integrale_noise, l, r, peaks+l, npeaks, amplitude]
        
        
        # if self.event_type == 'source':
        #     if npeaks == 1:
        #         return [integrale, integrale/(t[r]-t[l]), t[r]-t[l], integrale_noise]
        #     else:
        #         return [-1, -1, t[r]-t[l], integrale_noise]
        # else:
        #     return [integrale, integrale/(t[r]-t[l]), t[r]-t[l], integrale_noise]
        
        
def signal_model(t, a, b, c, taus):
    return np.where(t < c, b, a * (np.exp(-(t - c) / taus) - np.exp(-(t - c) / 7.8e-9)))

    #7.8e-9 costante RC PMT

        
def fwhm(x, y):
    
    height = np.min(y)  # Valore massimo
    half_height = height / 2
    
    # Trova i punti dove y incrocia la mezza altezza
    intersections = np.where(np.diff(np.sign(y - half_height)))[0]
    
    # Assicura che ci siano almeno due punti per calcolare la FWHM
    if len(intersections) >= 2:
        # Calcola la FWHM come la differenza tra i punti estremi di incrocio
        fwhm_val = x[intersections[-1]] - x[intersections[0]]
    else:
        fwhm_val = 0  # Non è possibile calcolare la FWHM se ci sono meno di due incroci
        
    return [fwhm_val, half_height, intersections]    

