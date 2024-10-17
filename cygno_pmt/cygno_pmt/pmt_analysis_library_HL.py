"""
@author: SerAngelo
"""

import cygno as cy 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pmt_analysis_library as pal  
import pickle
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from scipy.stats import chi2

def find_common_elements_within_range(arrays, tolerance=2):
    common_elements = set(arrays[0])
    for array in arrays[1:]:
        current_set = set()
        for element in array:
            for el in common_elements:
                if abs(element-el) <= tolerance:
                    current_set.add(el)
        common_elements = current_set
    return common_elements

def filter_event(serial_number, c, no_tracks = False):
    c = list(zip(*c))
    c = [list(i) for i in c] 
    
    
    # DT = np.asarray(c[2])
    peaks = c[6]
    n_peaks = c[7]
    
    if no_tracks:
        if (np.array(n_peaks) > 1).any():
            print(f"Traccia rivelata nell'evento {serial_number}")
            return False
    
    #controllo eventi con 0 picchi
    if 0 in n_peaks:
        print(f"Il numero di picchi dell'evento {serial_number} è uguale a 0")
        return False
    
    common_peaks = find_common_elements_within_range(peaks)
    
    if len(common_peaks) == 0:
        print(f"Nessun picco coincidente nell'evento {serial_number}")
        return False
    
        
    
    
    
    
    return True





def charge_distribution(mfile, corrected, channels_offsets, PMTs, digitizer_type, n_max, event_type, mode=2, track_filter=True, signal_fit = False, resampling=False, show_plots = False, show_distribution1D=False, show_distribution2D=False, verbose=False, save=False):
    
    carica = []
    DEDx = []
    Dt = []
    carica_noise = []
    amplitude = []
    
    TOT_EVENT=0
    FILTERED_EVENT=0
    
    
    NBLOCKS = []
    NBLOCKS_FILTERED = []
    
    for event in mfile:
        
        
        
        serial_number = event.header.serial_number
        
        if serial_number >= n_max:
            if serial_number < 405:
                break
            continue
            
        bank_names = ", ".join(b.name for b in event.banks.values())
        if ('DGH0' not in bank_names):
            continue
        
        print("Event # %s of type ID %s contains banks %s" % (serial_number, event.header.event_id, bank_names))
        
        full_header= cy.daq_dgz_full2header(event.banks['DGH0'], verbose=False)
        w_fast, w_slow = cy.daq_dgz_full2array(event.banks['DIG0'], full_header, verbose=False, corrected=corrected, ch_offset=channels_offsets)
        
        
        
        offset = 8
        nblocks = len(w_slow)//offset
        FILTERED_EVENT_BLOCK = 0
        NBLOCKS.append(nblocks)
        TOT_EVENT+=nblocks
        #print('################',len(w_slow))

        for b in range(nblocks):
            
            o = b*offset
            print(f"Analisi blocco {b}")
            
            w_sum = np.asarray(w_slow[PMTs[0].channel+o])
            for i in range(1,len(PMTs)):
                w_sum += np.asarray(w_slow[PMTs[i].channel+o])
            w_mean = w_sum/len(PMTs)
              
            if mode == 2:
                N_plots = len(PMTs)
                
            elif mode ==1:
                N_plots = 1
            
            if digitizer_type == 'slow':
                pmt_events = [pal.PMT_event(pmt, serial_number, w_slow[pmt.channel+o], 'slow', event_type, 16) for pmt in PMTs]
                pmt_event_mean = pal.PMT_event(PMTs[0], serial_number, w_mean, 'slow', event_type, 16)
            
            else:
                pmt_events = [pal.PMT_event(pmt, serial_number, w_fast[pmt.channel+o], 'fast', event_type, 16) for pmt in PMTs]
                pmt_event_mean = pal.PMT_event(PMTs[0], serial_number, w_mean, 'fast', event_type, 16)
            
            
            if show_plots:
                fig, ax = plt.subplots(1, N_plots, figsize=(20, 10))
                if mode == 2:
                    c = [pmt_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=ax[i]) for i,pmt_event in enumerate(pmt_events)]    
                    check = filter_event(serial_number, c, no_tracks=track_filter)
                    if check==False:
                        continue
                    #C.append(c)
                elif mode == 1:
                    c = [pmt_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=None) for i,pmt_event in enumerate(pmt_events)]    
                    c_mean = pmt_event_mean.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=ax) 
                    check = filter_event(serial_number, c, no_tracks=track_filter)
                    if check==False:
                        continue
                    #CMEAN.append(c_mean)
                    
                plt.tight_layout()
                plt.show()
                
            else:
                if mode == 2:
                    c = [pmt_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=None) for pmt_event in pmt_events]
                    check = filter_event(serial_number, c, no_tracks=track_filter)
                    if check==False:
                        FILTERED_EVENT+=1
                        FILTERED_EVENT_BLOCK+=1
                        continue
                    #C.append(c)
                elif mode == 1:
                    c = [pmt_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=None) for i,pmt_event in enumerate(pmt_events)]    
                    check = filter_event(serial_number, c, no_tracks=track_filter)
                    if check==True:
                        c_mean = pmt_event_mean.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=None)                        
                    else:
                        FILTERED_EVENT+=1
                        FILTERED_EVENT_BLOCK+=1
                        continue
                              
        # c = list(zip(*c))
        # c = [list(i) for i in c]    
            
        # if event_type != 'background':
        #     if -1 in c[0]:
        #         print(f"Il numero di picchi dell'evento {serial_number} è diverso da 1")
        #         continue
        
        #     if any(abs(pmt_events[0].pmt_max_ampl() - pmt_event.pmt_max_ampl()) > 25 for pmt_event in pmt_events[1:]):
        #         print(f"I {len(PMTs)} picchi dell'evento {serial_number} non sono molto coincidenti")
        #         continue
           
         
                
            if mode == 2:
                c = list(zip(*c))
                c = [list(i) for i in c]   
                carica.append(c[0])
                DEDx.append(c[1])
                Dt.append(c[2])
                carica_noise.append(c[3])
                amplitude.append(c[8])
            elif mode == 1:
                carica.append(c_mean[0])
                DEDx.append(c_mean[1])
                Dt.append(c_mean[2])
                carica_noise.append(c_mean[3])
                amplitude.append(c_mean[8])
    
        NBLOCKS_FILTERED.append(FILTERED_EVENT_BLOCK)
    
    
    if mode == 2:  
        carica = list(zip(*carica))
        DEDx = list(zip(*DEDx))
        Dt = list(zip(*Dt))
        carica_noise = list(zip(*carica_noise))
        amplitude = list(zip(*amplitude))
    
     
        
    print("TOT_EVENTS=", TOT_EVENT)
    print("FILTERED_EVENTS=", FILTERED_EVENT)
    print("NBLOCKS=", NBLOCKS)
    print("NBLOCKS_FILTERED=", NBLOCKS_FILTERED)
    
    if show_distribution1D:        
   
        fig, ax = plt.subplots(1, N_plots, figsize=(20, 10))
        
        if mode == 2:
        
            for i, pmt in enumerate(PMTs):
                ax[i].hist(carica[i], bins=int(len(carica[i])), edgecolor='black', color='orange', alpha=0.2, label=f'{pmt.name} signal charge')
                ax[i].hist(carica_noise[i], bins=int(len(carica_noise[i])), edgecolor='black', color='blue', alpha=0.2, label=f'{pmt.name} noise charge')
                             
                
                ax[i].set_xlim(0,100)
                ax[i].set_xlabel('$\Delta Q$ [pC]')
                ax[i].set_ylabel('Frequency')
                ax[i].set_title(f'PMT{i} $\Delta Q$ distribution')
                
        elif mode == 1:
            
            counts1, bins1, _ = ax.hist(carica, bins=int(len(carica)), edgecolor='black', color='orange', alpha=0.2, label='signal charge', density=True)
            counts2, bins2, _ = ax.hist(carica_noise, bins=int(len(carica_noise)), edgecolor='black', color='blue', alpha=0.2, label='noise charge', density=True)

            ax.set_xlim(0,100)
            ax.set_xlabel('$\Delta Q$ [pC]')
            ax.set_ylabel('Frequency')
            ax.set_title('$\Delta Q$ distribution')
                
        
        plt.tight_layout()
        plt.legend()
        plt.show()
    
    
    if show_distribution2D:        
        
        fig, ax = plt.subplots(1, N_plots, figsize=(20, 10))
        if mode == 2:
            for i, pmt in enumerate(PMTs):
                ax[i].hist2d(Dt[i], carica[i], bins=int(len(carica[i])), cmap='viridis')
                ax[i].set_xlabel('$\Delta T$ [s]')
                ax[i].set_ylabel('$\Delta Q$ [pC]')
                ax[i].set_title(f'PMT{i} $\Delta Q vs \Delta T$')
        elif mode == 1:
            ax.hist2d(Dt, carica, bins=int(len(carica)), cmap='viridis')
            ax.set_xlabel('$\Delta T$ [s]')
            ax.set_ylabel('$\Delta Q$ [pC]')
            ax.set_title('$\Delta Q vs \Delta T$')
            
        plt.tight_layout()
        plt.legend()
        plt.show()
    
    if save:
    
        output_file = event_type+".pkl"
        with open(output_file, 'wb') as file:
            pickle.dump([carica, DEDx, Dt, carica_noise, amplitude, TOT_EVENT, FILTERED_EVENT, NBLOCKS, NBLOCKS_FILTERED] , file)
            
        print(f"I risultati sono stati correttamente salvati nel file binario: {output_file}")
        

    
    
    
    
    
def single_waveform_info(mfile, corrected, channels_offsets, PMTs, n, digitizer_type, event_type, signal_fit=False, resampling=False, verbose=False):
    
    
    for event in mfile:
            
        serial_number = event.header.serial_number
        
        
        if serial_number == n:
        
            bank_names = ", ".join(b.name for b in event.banks.values())
            
            if ('DGH0' not in bank_names):
                continue
            
            print("Event # %s of type ID %s contains banks %s" % (serial_number, event.header.event_id, bank_names))
            
            full_header= cy.daq_dgz_full2header(event.banks['DGH0'], verbose=False)
            w_fast, w_slow = cy.daq_dgz_full2array(event.banks['DIG0'], full_header, verbose=False, corrected=corrected, ch_offset=channels_offsets)
            
            if digitizer_type == 'slow':
                pmt1_event = pal.PMT_event(PMTs[0], serial_number, w_slow[2], digitizer_type, event_type, 16)
                pmt2_event = pal.PMT_event(PMTs[1], serial_number, w_slow[3], digitizer_type, event_type, 16)
                
            elif digitizer_type == 'fast':
                pmt1_event = pal.PMT_event(PMTs[0], serial_number, w_fast[2], digitizer_type, event_type, 16)
                pmt2_event = pal.PMT_event(PMTs[1], serial_number, w_fast[3], digitizer_type, event_type, 16)
                

            fig = plt.figure(figsize=(20, 10))
            gs = gridspec.GridSpec(1, 4, width_ratios=[5, 1, 5, 1])
            
            ax0 = fig.add_subplot(gs[0])
            ax1 = fig.add_subplot(gs[1])
            ax2 = fig.add_subplot(gs[2])
            ax3 = fig.add_subplot(gs[3])
            
            
            
            c1 = pmt1_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=ax0)
            c2 = pmt2_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=ax2)
            
            
            if (digitizer_type == 'fast' or event_type == 'background') and resampling==False:
                print('Attenzione, resampling consigliato')
             
    
            pmt1_event.gaussian_noise(ax1)
            ax1.set_ylim(ax0.get_ylim())
            
            pmt2_event.gaussian_noise(ax3)
            ax3.set_ylim(ax2.get_ylim())
            
            plt.tight_layout()
            plt.show()
            
            fig.savefig(f'{event_type}_{n}.png', dpi=100)
            
            
            break
        
    return [c1,c2]

from scipy.stats import poisson

# Definizione della funzione di Poisson continua per il fit
def poisson_pmf(k, lamb,A):
    return A*poisson.pmf(k, lamb)

def exp(x, l, A):
    return A*np.exp(-x/l)

from landaupy import landau

def landau_pdf(x, mpv, xi):
    return landau.pdf(x,mpv,xi)


def distribution1D_mean_BG(file_name_bg1, file_name_bg2):
    
    plt.rcParams.update({'font.size': 30})

    try:
        with open(file_name_bg1, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data3 = pickle.load(file)
    
        carica_bg1, DEDx_bg1, Dt_bg1, carica_noise_bg1, amplitude_bg1, total_bg1, filtered_bg1, nblocks_bg1, nblocks_filtered_bg1 = data3
            
        with open(file_name_bg2, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data4 = pickle.load(file)
    
        carica_bg2, DEDx_bg2, Dt_bg2, carica_noise_bg2, amplitude_bg2, total_bg2, filtered_bg2, nblocks_bg2, nblocks_filtered_bg2 = data4
        
        carica_bg = carica_bg1+carica_bg2
        # DEDx_bg = np.asarray(DEDx_bg1+DEDx_bg2)*10**(-9)
        # Dt_bg = Dt_bg1+Dt_bg2
        # amplitude_bg = amplitude_bg1+amplitude_bg2
        #carica_noise_bg = carica_noise_bg1+carica_noise_bg2
    
    except FileNotFoundError:
        print(f"'file_name_bg1' o '{file_name_bg2}' non trovati! Effettuare un'analisi su due run di background!")
        return [0,0]
   
    # fig, ax = plt.subplots(1, 1, figsize=(20, 10))

    # ax.set_title('$\Delta Q$ background distribution')
    # ax.set_xlabel('$\Delta Q$ [pC]')
    # ax.set_ylabel('Frequency')
    #ax.set_xlim(0,200)
    
    for n in range(3,10):
        counts1, bins1 = np.histogram(carica_bg, bins=len(carica_bg)//n, density=True)
    
    
      #  counts1, bins1, _ = ax.hist(carica_bg, bins=len(carica_bg)//n, edgecolor='none', color='orange', alpha=0.5, label='background charge', density=True)
    
        x_S = np.linspace(0, 600, 10000)
        
      
        p0 = [0.1, 1]
       
        
        
        bin_centers_S = (bins1[:-1] + bins1[1:]) / 2
        
        popt_S, pcov_S = curve_fit(landau_pdf, bin_centers_S, counts1, p0) 
       
   #     ax.plot(x_S, landau_pdf(x_S, *popt_S), 'k-', label='Landau Fit', lw=3, zorder=5)
    #    plt.legend()
        
        expected_counts = landau_pdf(bin_centers_S, *popt_S) * (bins1[1] - bins1[0]) * len(carica_bg)
    
        observed_counts = counts1 * len(carica_bg)
        chi_squared = np.sum((observed_counts - expected_counts) ** 2 / expected_counts)
        degrees_of_freedom = len(bin_centers_S) - len(popt_S)
        reduced_chi_squared = chi_squared / degrees_of_freedom
        p_value = chi2.sf(chi_squared, df=degrees_of_freedom)
        if p_value>0.05:
            print(f"Chi-squared: {chi_squared}")
            print(f"Degrees of freedom: {degrees_of_freedom}")
            print(f"Reduced Chi-squared: {reduced_chi_squared}")
            print(f"P-value: {p_value}")
            break
        

    
    # plt.savefig("landau.png", dpi=100)
    return popt_S

def distribution1D_mean(file_name1, file_name2, file_name_bg1, file_name_bg2, obs, v):
    
    plt.rcParams.update({'font.size': 30})
  
    # file_name_bg1 = "background1_AUGMENTED.pkl"
    # file_name_bg2 = "background2_AUGMENTED.pkl"
    
    try:
        with open(file_name_bg1, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data3 = pickle.load(file)
    
        carica_bg1, DEDx_bg1, Dt_bg1, carica_noise_bg1, amplitude_bg1, total_bg1, filtered_bg1, nblocks_bg1, nblocks_filtered_bg1 = data3
            
        with open(file_name_bg2, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data4 = pickle.load(file)
    
        carica_bg2, DEDx_bg2, Dt_bg2, carica_noise_bg2, amplitude_bg2, total_bg2, filtered_bg2, nblocks_bg2, nblocks_filtered_bg2 = data4
        
        carica_bg = carica_bg1+carica_bg2
        DEDx_bg = np.asarray(DEDx_bg1+DEDx_bg2)*10**(-9)
        Dt_bg = Dt_bg1+Dt_bg2
        amplitude_bg = amplitude_bg1+amplitude_bg2
        #carica_noise_bg = carica_noise_bg1+carica_noise_bg2
        total_bg = total_bg1 + total_bg2
        filtered_bg = filtered_bg1+filtered_bg2
        perc = filtered_bg/total_bg
        
    
    except FileNotFoundError:
        print(f"'{file_name_bg1}' o '{file_name_bg2}' non trovati! Effettuare un'analisi su due run di background!")
        return [0,0]
    #print("###########",len(carica_bg))
    #print("############", np.mean(carica_S),'+/-',np.std(carica_S)/np.sqrt(len(carica_S)))
        
    
    with open(file_name1, 'rb') as file:
        # Carica gli oggetti dal file pickle
        data1 = pickle.load(file)
    
    carica_S1, DEDx_S1, Dt_S1, carica_noise_S1, amplitude_S1, total_S1, filtered_S1, nblocks_S1, nblocks_filtered_S1 = data1

    
    with open(file_name1, 'rb') as file:
        # Carica gli oggetti dal file pickle
        data2 = pickle.load(file)

    carica_S2, DEDx_S2, Dt_S2, carica_noise_S2, amplitude_S2, total_S2, filtered_S2, nblocks_S2, nblocks_filtered_S2 = data2

    carica_S = carica_S1+carica_S2
    DEDx_S = np.asarray(DEDx_S1+DEDx_S2)*10**(-9)
    Dt_S = Dt_S1+Dt_S2
    carica_noise_S = carica_noise_S1+carica_noise_S2
    amplitude_S = amplitude_S1+amplitude_S2
    filtered_S = filtered_S1+filtered_S2
    nblocks = nblocks_S1 + nblocks_S2
    nblocks_filtered = nblocks_filtered_S1 + nblocks_filtered_S2
    survived_events = int(filtered_S*(1-perc)/perc)
    
    
    
    #print('check1')
    fig, ax = plt.subplots(1, 1, figsize=(20, 15))

    if obs == 'charge':
        popts = distribution1D_mean_BG(file_name_bg1, file_name_bg2)
        f = lambda x: landau_pdf(x, *popts)
        # Definiamo i limiti di integrazione
        

        
        
        
        ax.set_title('$\Delta Q$ distribution')
        ax.set_xlabel('$\Delta Q$ [pC]')
        ax.set_xlim(0,50)
        #ax.set_ylim(0,700)
        counts1, bins1, _ = ax.hist(carica_S, bins=len(carica_S)//40, edgecolor='black', color='red', alpha=1, label='signal charge', zorder=2)
        counts2, bins2, _ = ax.hist(carica_noise_S, bins=len(carica_noise_S)//40, edgecolor='none', color='blue', alpha=1, label='noise charge', zorder=3)
      
        x = bins1
        y = f(x)

        # # Calcolo dell'area sotto la curva usando il metodo dei rettangoli
        # area_originale = np.sum(y * (bins1[1]-bins1[0]))

        # Determinazione del fattore di normalizzazione
        N = survived_events

        # Definizione della funzione normalizzata
        f_norm = lambda x: N * f(x) * (bins1[1]-bins1[0])

        # Generazione dei dati per la curva normalizzata
        y_norm = f_norm(x)
        #print(y_norm)

        # Calcolo dell'area sotto la curva normalizzata per verifica
        area_normalizzata = np.sum(y_norm)
        #print('area_normalizata',area_normalizzata) 
        
        counts3, bins3, _ = ax.hist(x, bins=bins1, weights=y_norm, alpha=0.5, label='Predicted background charge', edgecolor='black', color='black', zorder=5)

        #counts3, bins3, _ = ax.hist(carica_bg, bins=bins1, edgecolor='none', color='black', alpha=1, label='background charge', zorder=4)
        x_S = np.linspace(0, 100, 10000)
        
        p0=[np.mean(carica_S), np.std(carica_S),0.1]
        #p0 = [np.mean(carica_S), np.max(counts1)]
        #print('check2')
        
    elif obs == 'Dt':
        
        ax.set_title('$\Delta T$ distribution')
        ax.set_xlabel('$\Delta T$ [s]')
        #ax.set_ylim(0,2800)
        if v == 440:
            counts1, bins1, _ = ax.hist(Dt_S, bins=150, edgecolor='black', color='red', alpha=1, label='signal charge', zorder=1)
        elif v == 430:
            counts1, bins1, _ = ax.hist(Dt_S, bins=100, edgecolor='black', color='red', alpha=1, label='signal charge', zorder=1)
        elif v == 420:
            counts1, bins1, _ = ax.hist(Dt_S, bins=100, edgecolor='black', color='red', alpha=1, label='signal charge', zorder=1)
        
        #counts3, bins3, _ = ax.hist(Dt_bg, bins=bins1, edgecolor='black', color='black', alpha=1, label='background charge', zorder=4)
        
        counts, bins = np.histogram(Dt_bg, bins=bins1, density=True)
        
        # Determinazione del fattore di normalizzazione
        N = survived_events

        # Definizione della funzione normalizzata
        y_norm = N * counts * (bins1[1]-bins1[0])
        #print(y_norm)
        # Calcolo dell'area sotto la curva normalizzata per verifica
        area_normalizzata = np.sum(y_norm)
        #print('area_normalizata',area_normalizzata) 
        
        counts3, bins3, _ = ax.hist(bins1[:-1], bins=bins1, weights=y_norm, alpha=0.5, label='Predicted background charge', edgecolor='black', color='black', zorder=5)

        
        
        
        ax.set_xlim(0,2e-7)
        x_S = np.linspace(0, 2e-7, 100000)
    
        p0=[np.mean(Dt_S), np.std(Dt_S),50]
        
    elif obs == 'DEDx':
        
        ax.set_title('$\Delta Q/\Delta T$ distribution')
        ax.set_xlabel('$\Delta Q/\Delta T$ [pC/ns]')
        ax.set_xlim(0,0.6)
        ax.set_ylim(0,700)
        
        counts1, bins1, _ = ax.hist(DEDx_S, bins=len(DEDx_S)//100, edgecolor='black', color='red', alpha=1, label='signal charge', zorder=1)
        #counts3, bins3, _ = ax.hist(DEDx_bg, bins=bins1, edgecolor='black', color='black', alpha=1, label='background charge', zorder=2)
        
        counts, bins = np.histogram(DEDx_bg, bins=bins1, density=True)
        
        # Determinazione del fattore di normalizzazione
        N = survived_events

        # Definizione della funzione normalizzata
        y_norm = N * counts * (bins1[1]-bins1[0])
        #print(y_norm)
        # Calcolo dell'area sotto la curva normalizzata per verifica
        area_normalizzata = np.sum(y_norm)
        #print('area_normalizata',area_normalizzata) 
        
        counts3, bins3, _ = ax.hist(bins1[:-1], bins=bins1, weights=y_norm, alpha=0.5, label='Predicted background charge', edgecolor='black', color='black', zorder=5)

        
        
        x_S = np.linspace(0, 1, 1000)
    
        p0=[np.mean(DEDx_S), np.std(DEDx_S),25]
    
    elif obs == 'amplitude':
        
        ax.set_title('Amplitude distribution')
        ax.set_xlabel('Amplitude [mV]')
        ax.set_xlim(900,1050)
        #ax.set_ylim(0,1000)
        counts1, bins1, _ = ax.hist(amplitude_S, bins=len(amplitude_S)//50, edgecolor='black', color='red', alpha=1, label='signal amplitude', zorder=2)
        #counts3, bins3, _ = ax.hist(amplitude_bg, bins=bins1, edgecolor='none', color='black', alpha=1, label='background amplitude', zorder=4)
        
        counts, bins = np.histogram(amplitude_bg, bins=bins1, density=True)
        
        # Determinazione del fattore di normalizzazione
        N = survived_events

        # Definizione della funzione normalizzata
        y_norm = N * counts * (bins1[1]-bins1[0])
        #print(y_norm)
        # Calcolo dell'area sotto la curva normalizzata per verifica
        area_normalizzata = np.sum(y_norm)
        #print('area_normalizata',area_normalizzata) 
        
        counts3, bins3, _ = ax.hist(bins1[:-1], bins=bins1, weights=y_norm, alpha=0.5, label='Predicted background charge', edgecolor='black', color='black', zorder=5)

        
        
        
        x_S = np.linspace(900, 1100, 10000)
        
        p0=[np.mean(amplitude_S), np.std(amplitude_S),np.max(counts1)]
        #p0 = [np.mean(carica_S), np.max(counts1)]
        #print('check2')
        
    
    bin_centers_S = (bins1[:-1] + bins1[1:]) / 2
    bin_width = bins1[1]-bins1[0]
    
    corrected_counts = np.maximum(counts1-counts3,0)
    #corrected_counts=counts1
    sigma = np.sqrt(counts1 + counts3)
    ax.bar(bin_centers_S, corrected_counts, width=np.diff(bins1), edgecolor='black', color='yellow', align='center', alpha=1, label='Corrected charge distribution', zorder=3)
    
    for center, count, error in zip(bin_centers_S, corrected_counts, sigma):
        plt.gca().add_patch(plt.Rectangle((center-bin_width/2, count-error), bin_width, 2*error, edgecolor='black', facecolor='none', hatch='//', alpha=0.5, zorder=5))
    
    popt_S, pcov_S = curve_fit(gaussian, bin_centers_S, corrected_counts, p0) 
   
    ax.plot(x_S, gaussian(x_S, *popt_S), 'k-', label='Fit signal charge', lw=3, zorder=5)
    #print(popt_S, pcov_S)
  
    if obs=='charge':
        bin_centers_noise = (bins2[:-1] + bins2[1:]) / 2
        popt_noise, pcov_noise = curve_fit(gaussian, bin_centers_noise, counts2, p0=[np.mean(carica_noise_S), np.std(carica_noise_S),0.1]) 
        x_noise = np.linspace(0, 10)
        ax.plot(x_noise, gaussian(x_noise, *popt_noise), color='black', linestyle='--', lw=4, label='Fit noise charge')
        #print(popt_noise, pcov_noise)

    error_patch = Patch(facecolor='none', edgecolor='black', hatch='//', label='Error range')
    handles, labels = plt.gca().get_legend_handles_labels()    
    handles.append(error_patch)
    labels.append('Error range')
    
    
    
    ax.set_ylabel('Frequency')
    
     
    
    plt.tight_layout()
    plt.legend(handles=handles, labels=labels)
    
    plt.savefig(f"{obs}_distribution.png", dpi=100)
    
    plt.show()
    
    if obs == 'charge':
        return [popt_S[0], np.sqrt(pcov_S[0][0]+(popt_noise[1]/np.sqrt(len(bins2)))**2)]#,[popt_S[1], pcov_S[1][1]**0.5]
    elif obs == 'Dt':
        return [popt_S[0], pcov_S[0][0]**0.5]
    elif obs == 'DEDx':
        return [popt_S[0], pcov_S[0][0]**0.5]
    elif obs == 'amplitude':
        return [popt_S[0], pcov_S[0][0]**0.5]
    
    

def distribution1D(PMTs, file_name1, file_name2, obs, voltage, pos):
    
    #plt.rcParams.update({'font.size': 30})
  
    with open(file_name1, 'rb') as file:
        # Carica gli oggetti dal file pickle
        data1 = pickle.load(file)
    
    carica_S1, DEDx_S1, Dt_S1, carica_noise_S1, amplitude_S1 = data1

    
    with open(file_name1, 'rb') as file:
        # Carica gli oggetti dal file pickle
        data2 = pickle.load(file)

    carica_S2, DEDx_S2, Dt_S2, carica_noise_S2, amplitude_S2 = data2

    carica_S = carica_S1+carica_S2
    DEDx_S = np.asarray(DEDx_S1+DEDx_S2)
    Dt_S = Dt_S1+Dt_S2
    carica_noise_S = carica_noise_S1+carica_noise_S2
    amplitude_S = amplitude_S1+amplitude_S2
    

    file_name_bg1 = "background1_AUGMENTED.pkl"
    file_name_bg2 = "background2_AUGMENTED.pkl"
    
    try:
        with open(file_name_bg1, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data3 = pickle.load(file)
    
        carica_bg1, DEDx_bg1, Dt_bg1, carica_noise_bg1, amplitude_bg1 = data3
            
        with open(file_name_bg2, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data4 = pickle.load(file)
    
        carica_bg2, DEDx_bg2, Dt_bg2, carica_noise_bg2, amplitude_bg2 = data4
        
        carica_bg = carica_bg1+carica_bg2
        DEDx_bg = np.asarray(DEDx_bg1+DEDx_bg2)
        Dt_bg = Dt_bg1+Dt_bg2
        amplitude_bg = amplitude_bg1+amplitude_bg2
        #carica_noise_bg = carica_noise_bg1+carica_noise_bg2
    
    except FileNotFoundError:
        print("'background1_AUGMENTED.pkl' e/o 'background2_AUGMENTED.pkl' non trovati! Effettuare un'analisi su due run di background!")
        return [0,0]
    #print("###########",len(carica_bg))
    #print("############", np.mean(carica_S),'+/-',np.std(carica_S)/np.sqrt(len(carica_S)))
        
    #print('check1')
    
    Smean = []
    Sstd = []

    if obs == 'charge':
        
        plt.rcParams.update({'font.size': 20})
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        for i, pmt in enumerate(PMTs):
            ax[i].set_title('$\Delta Q$ distribution')
            ax[i].set_xlabel('$\Delta Q$ [pC]')
            ax[i].set_ylabel('Frequency')
            #ax[i].set_xlim(0, 100)
        
            # Determina il numero di bins e il range in base ai dati combinati
            all_data = np.concatenate((carica_S[i], carica_bg[i]))
            bins = np.linspace(0, 200, 50)  # 50 bins tra 0 e 100, puoi modificare a seconda delle tue necessità
        
            # Istogramma del segnale con lo stesso bin size del background
            ax[i].hist(carica_S[i], bins=bins, edgecolor='black', color='orange', alpha=0.5, label=f'{pmt.name} signal')
            
            # Istogramma del rumore di fondo con lo stesso bin size del segnale
            #ax[i].hist(carica_bg[i], bins=bins, edgecolor='black', color='black', alpha=0.5, label=f'{pmt.name} background charge')
            
            # Calcola la media e la deviazione standard del segnale
            Smean.append(np.mean(carica_S[i]))
            Sstd.append(np.std(carica_S[i]) / np.sqrt(len(carica_S[i])))
    
            ax[i].legend()  # Aggiungi la leggenda
        
    elif obs == 'Dt':
        
        plt.rcParams.update({'font.size': 15})
        fig, ax = plt.subplots()

        DDt = np.asarray(Dt_S[0])-np.asarray(Dt_S[1])
        Smean.append(np.mean(DDt))
        Sstd.append(np.std(DDt)/np.sqrt(len(Dt_S[0])))
        
        ax.set_title(f'$\Delta$Amplitude distribution - {voltage}V - P{pos}')
        ax.set_xlabel('Amplitude [mV]')
        ax.set_ylabel('Frequency')
        ax.grid()
        #ax.set_xlim(900,1100)
        ax.hist(DDt, bins=len(Dt_S[0])//100, edgecolor='none', color='orange', alpha=0.5)
        
        textstr = '\n'.join((
            r'$\mu=%.2e\,\mathrm{mV}$' % (0, ),
            r'$\sigma=%.2e\,\mathrm{mV}$' % (Sstd[0]*np.sqrt(len(amplitude_S[0])), )))
        
        # Adjusted properties for the text box
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        
        # Place the text box below the legend
        plt.gca().text(0.70, 0.80, textstr, transform=plt.gca().transAxes, fontsize=15,
                       verticalalignment='top', bbox=props, color='red')
        
        plt.savefig(f"DeltaDt{voltage}.png", dpi=1000)
        
        
    elif obs == 'DEDx':
        
        plt.rcParams.update({'font.size': 20})
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        for i, pmt in enumerate(PMTs):
            ax[i].set_title('$\Delta Q/\Delta T$ distribution')
            ax[i].set_xlabel('$\Delta Q/\Delta T$ [pC/s]')
            ax[i].set_ylabel('Frequency')
            #ax[i].set_xlim(0, 0.5e-6)
        
            # Determina il numero di bins e il range in base ai dati combinati
            #bins = np.linspace(0, 0.2e-7, 50)  # 50 bins tra 0 e 100, puoi modificare a seconda delle tue necessità
        
            # Istogramma del segnale con lo stesso bin size del background
            ax[i].hist(DEDx_S[i], bins=100, edgecolor='black', color='orange', alpha=0.5, label=f'{pmt.name} signal')
            
            # Istogramma del rumore di fondo con lo stesso bin size del segnale
            #ax[i].hist(carica_bg[i], bins=bins, edgecolor='black', color='black', alpha=0.5, label=f'{pmt.name} background charge')
            
            # Calcola la media e la deviazione standard del segnale
            Smean.append(np.mean(DEDx_S[i]))
            Sstd.append(np.std(DEDx_S[i])/np.sqrt(len(DEDx_S[i])))
     
            ax[i].legend()  # Aggiungi la leggenda
        
      
    
    elif obs == 'amplitude':
        
        plt.rcParams.update({'font.size': 15})
        fig, ax = plt.subplots()

        Damplitude = np.asarray(amplitude_S[0])-np.asarray(amplitude_S[1])
        Smean.append(np.mean(Damplitude))
        Sstd.append(np.std(Damplitude)/np.sqrt(len(amplitude_S[0])))
        # for i, pmt in enumerate(PMTs):
        ax.set_title(f'$\Delta$Amplitude distribution - {voltage}V - P{pos}')
        ax.set_xlabel('Amplitude [mV]')
        ax.set_ylabel('Frequency')
        ax.grid()
        #ax.set_xlim(900,1100)
        ax.hist(Damplitude, bins=len(amplitude_S[0])//60, edgecolor='none', color='orange', alpha=0.5)
        
        # ax.axvline(x=0, color='black', linestyle='--', linewidth=1, label='Zero')
        ax.axvline(x=Smean[0], color='red', linestyle='--', linewidth=1)

        #ax.set_xlim(-200,200)
        
        #ax[i].hist(amplitude_bg[i], bins=len(amplitude_bg[i])//20, edgecolor='none', color='blue', alpha=0.2, label=f'{pmt.name} noise charge')
        # Adding mean and standard deviation box
        textstr = '\n'.join((
            r'$\mu=%.2f\,\mathrm{mV}$' % (Smean[0], ),
            r'$\sigma=%.2f\,\mathrm{mV}$' % (Sstd[0]*np.sqrt(len(amplitude_S[0])), )))
        
        # Adjusted properties for the text box
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        
        # Place the text box below the legend
        plt.gca().text(0.70, 0.80, textstr, transform=plt.gca().transAxes, fontsize=15,
                       verticalalignment='top', bbox=props, color='red')
        
        plt.savefig(f"DeltaAmplitude{voltage}.png", dpi=100)

    elif obs == 'Dt&amplitude':
        plt.rcParams.update({'font.size': 20})
       
        # Calcola le differenze
        DDt = np.asarray(Dt_S[0]) - np.asarray(Dt_S[1])
        Damplitude = np.asarray(amplitude_S[0]) - np.asarray(amplitude_S[1])
        
        # Calcola le differenze
        DDt = np.asarray(Dt_S[0]) - np.asarray(Dt_S[1])
        Damplitude = np.asarray(amplitude_S[0]) - np.asarray(amplitude_S[1])
        
        # Calcola le medie e le deviazioni standard
        mean_DDt = np.mean(DDt)
        std_DDt = np.std(DDt)
        mean_Damplitude = np.mean(Damplitude)
        std_Damplitude = np.std(Damplitude)
        
        # Crea il grafico
        gs = GridSpec(4, 4)
        fig = plt.figure(figsize=(10, 8))
        ax_main = fig.add_subplot(gs[1:4, 0:3])
        ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
        ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
        
        # Imposta titolo e etichette
        ax_main.set_title(f'$\Delta A$ vs $\Delta \Delta T$ - {voltage}V')
        ax_main.set_xlabel('$\Delta \Delta T$ [s]')
        ax_main.set_ylabel('$\Delta A$ [mV]')
        
        # Primo hist2d
        h = ax_main.hist2d(DDt, Damplitude, bins=[80, 80], cmap='plasma', alpha=0.8)
        
        # Istogrammi marginali
        ax_xhist.hist(DDt, bins=80, edgecolor='black', color='blue', alpha=0.5)
        ax_yhist.hist(Damplitude, bins=80, orientation='horizontal', edgecolor='black', color='red', alpha=0.5)
        
        # Imposta limiti degli assi
        ax_main.set_xlim(min(DDt), max(DDt))
        ax_main.set_ylim(min(Damplitude), max(Damplitude))
        
        # Nascondi etichette degli assi degli istogrammi marginali
        plt.setp(ax_xhist.get_xticklabels(), visible=False)
        plt.setp(ax_yhist.get_yticklabels(), visible=False)
        
        # Aggiungi una barra colore per il hist2d
        #plt.colorbar(h[3], ax=ax_main)
        
        # Aggiungi linee tratteggiate per le medie
        ax_main.axvline(x=mean_DDt, color='red', linestyle='--', linewidth=2)
        ax_main.axhline(y=mean_Damplitude, color='red', linestyle='--', linewidth=2)
        
        # Aggiungi testo con le medie e le deviazioni standard sui grafici marginali
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        
        # Testo sulla parte superiore (asse x)
        ax_xhist.text(0.95, 0.85, f'Mean: {mean_DDt:.2e}\nStd: {std_DDt:.2e}', transform=ax_xhist.transAxes, fontsize=15,
                      verticalalignment='top', horizontalalignment='right', bbox=props, color='red')
        
        # Testo sulla parte destra (asse y)
        ax_yhist.text(0.05, 0.95, f'Mean: {mean_Damplitude:.2f}\nStd: {std_Damplitude:.2f}', transform=ax_yhist.transAxes, fontsize=15,
                      verticalalignment='top', horizontalalignment='left', bbox=props, color='red')
        
     
        
    
    plt.tight_layout()
    # plt.legend(handles=handles, labels=labels)
    plt.show()

        
        
        
    
       
  
    
    
    return Smean, Sstd
    
    # if obs == 'charge':
    #     return [popt_S[0], np.sqrt(pcov_S[0][0]+popt_noise[1]**2)]#,[popt_S[1], pcov_S[1][1]**0.5]
    # elif obs == 'Dt':
    #     return [np.mean(Dt_S), np.std(Dt_S)/len(Dt_S)]
    # elif obs == 'DEDx':
    #     return [popt_S[0], pcov_S[0][0]**0.5]
    # elif obs == 'amplitude':
    #     return [popt_S[0], pcov_S[0][0]**0.5]
    


def distribution2D_mean(run):
    
    carica_S = []
    DEDx_S = []
    Dt_S = []
    carica_noise_S = []
    amplitude_S = []
    
    for r in run:
        with open(f"source{r}_AUGMENTED.pkl", 'rb') as file:
            # Carica gli oggetti dal file pickle
            data1 = pickle.load(file)
        
        carica_S1, DEDx_S1, Dt_S1, carica_noise_S1, amplitude_S1 = data1
    
        
        with open(f"source{r+1}_AUGMENTED.pkl", 'rb') as file:
            # Carica gli oggetti dal file pickle
            data2 = pickle.load(file)
    
        carica_S2, DEDx_S2, Dt_S2, carica_noise_S2, amplitude_S2 = data2
    
        carica_S = carica_S + carica_S1+carica_S2
        DEDx_S = DEDx_S + DEDx_S1+DEDx_S2
        Dt_S = Dt_S + Dt_S1+Dt_S2
        carica_noise_S = carica_noise_S + carica_noise_S1+carica_noise_S2
        amplitude_S = amplitude_S + amplitude_S1 + amplitude_S2
        

    file_name_bg1 = "background1_AUGMENTED.pkl"
    file_name_bg2 = "background2_AUGMENTED.pkl"
    
    try:
        with open(file_name_bg1, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data3 = pickle.load(file)
    
        carica_bg1, DEDx_bg1, Dt_bg1, carica_noise_bg1, amplitude_bg1 = data3
            
        with open(file_name_bg2, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data4 = pickle.load(file)
    
        carica_bg2, DEDx_bg2, Dt_bg2, carica_noise_bg2, amplitude_bg2 = data4
        
        carica_bg = carica_bg1+carica_bg2
        DEDx_bg = DEDx_bg1+DEDx_bg2
        Dt_bg = Dt_bg1+Dt_bg2
        carica_noise_bg = carica_noise_bg1+carica_noise_bg2
        amplitude_bg = amplitude_bg1 + amplitude_bg2
        
    except FileNotFoundError:
        print("'background1.pkl' o 'background2.pkl' non trovati! Effettuare un'analisi su due run di background!")
       


        
    # fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    
    # ax.hist2d(Dt_S, carica_S, bins=int(len(carica_S)), cmap='viridis')
    # #ax.hist2d(Dt_S, carica_noise_S, bins=int(len(carica_noise_S)), cmap='viridis')
    # ax.hist2d(Dt_bg, carica_bg, bins=int(len(carica_bg)), cmap='viridis')
    
    
    # ax.set_xlabel('$\Delta T$ [s]')
    # ax.set_ylabel('$\Delta Q$ [pC]')
    # ax.set_title('$\Delta Q vs \Delta T$')
    from scipy.stats import gaussian_kde

    # Calcola la densità di probabilità usando la stima della densità del nucleo (KDE)
    kde = gaussian_kde(np.vstack([DEDx_S, amplitude_S]))
    density_bg = kde(np.vstack([DEDx_bg, amplitude_bg]))
    density_S = kde(np.vstack([DEDx_S, amplitude_S]))
    
    # Definisce una griglia su cui valutare la densità KDE
    x = np.linspace(np.asarray(DEDx_S).min(), np.asarray(DEDx_S).max(), 100)
    y = np.linspace(np.asarray(amplitude_S).min(), np.asarray(amplitude_S).max(), 100)
    X, Y = np.meshgrid(x, y)
    Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    
    # Trova il contorno al 95% del livello di confidenza
    sorted_Z = np.sort(Z.ravel())[::-1]
    cumulative_sum = np.cumsum(sorted_Z)
    cumulative_sum /= cumulative_sum[-1]
    threshold = sorted_Z[np.searchsorted(cumulative_sum, 0.95)]
    
    # Traccia il grafico con i contorni e scatter plot con cmap
    plt.figure(figsize=(10, 8))
    plt.contour(X, Y, Z, levels=[threshold], colors='red', linewidths=3)
    plt.scatter(DEDx_S, amplitude_S, c=density_S, s=5, cmap='Blues')
    plt.scatter(DEDx_bg, amplitude_bg, c=density_bg, s=5, cmap='Reds')

    
    plt.title('Contour at 95% CL with Density Colormap')
    plt.xlabel(r'$\Delta T$ [s]')
    plt.ylabel('Amplitude [mV]')
    plt.show()


    
    plt.tight_layout()
    # plt.legend(handles=handles, labels=labels)
    plt.savefig("distribution_2d.png", dpi=100)
    plt.show()
    
    return kde, threshold


def trigger_efficiency(file_name1, file_name2, file_name_bg1, file_name_bg2):
    
    try:
        with open(file_name_bg1, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data3 = pickle.load(file)
    
        carica_bg1, DEDx_bg1, Dt_bg1, carica_noise_bg1, amplitude_bg1, total_bg1, filtered_bg1, nblocks_bg1, nblocks_filtered_bg1 = data3
            
        with open(file_name_bg2, 'rb') as file:
            # Carica gli oggetti dal file pickle
            data4 = pickle.load(file)
    
        carica_bg2, DEDx_bg2, Dt_bg2, carica_noise_bg2, amplitude_bg2, total_bg2, filtered_bg2, nblocks_bg2, nblocks_filtered_bg2 = data4
        
        # carica_bg = carica_bg1+carica_bg2
        # DEDx_bg = np.asarray(DEDx_bg1+DEDx_bg2)*10**(-9)
        # Dt_bg = Dt_bg1+Dt_bg2
        # amplitude_bg = amplitude_bg1+amplitude_bg2
        #carica_noise_bg = carica_noise_bg1+carica_noise_bg2
        total_bg = total_bg1 + total_bg2
        filtered_bg = filtered_bg1+filtered_bg2
        perc = filtered_bg/total_bg
    
    except FileNotFoundError:
        print(f"'{file_name_bg1}' o '{file_name_bg2}' non trovati! Effettuare un'analisi su due run di background!")
        return [0,0]
    
    
    with open(file_name1, 'rb') as file:
        # Carica gli oggetti dal file pickle
        data1 = pickle.load(file)
    
    carica_S1, DEDx_S1, Dt_S1, carica_noise_S1, amplitude_S1, total_S1, filtered_S1, nblocks_S1, nblocks_filtered_S1 = data1

    
    with open(file_name1, 'rb') as file:
        # Carica gli oggetti dal file pickle
        data2 = pickle.load(file)

    carica_S2, DEDx_S2, Dt_S2, carica_noise_S2, amplitude_S2, total_S2, filtered_S2, nblocks_S2, nblocks_filtered_S2 = data2

    # carica_S = carica_S1+carica_S2
    # DEDx_S = np.asarray(DEDx_S1+DEDx_S2)*10**(-9)
    # Dt_S = Dt_S1+Dt_S2
    # carica_noise_S = carica_noise_S1+carica_noise_S2
    # amplitude_S = amplitude_S1+amplitude_S2
    # filtered_S = filtered_S1+filtered_S2
    nblocks = nblocks_S1 + nblocks_S2
    nblocks_filtered = nblocks_filtered_S1 + nblocks_filtered_S2
    nblocks_survived = np.asarray(nblocks_filtered)*(1-perc)/perc
    print(nblocks_survived)
    
    fig, ax = plt.subplots(1, 1, figsize=(20, 15))

    ax.set_title('Trigger events distribution')
    ax.set_xlabel('#Trigger events')
    #ax.set_xlim(0,50)
    #ax.set_ylim(0,700)
    #counts1, bins1, _ = ax.hist(np.asarray(nblocks)-np.asarray(nblocks_filtered), bins=len(nblocks)//20, edgecolor='black', color='red', alpha=0.5, label='signal charge', zorder=1)
    counts2, bins2, _ = ax.hist(np.asarray(nblocks)-np.asarray(nblocks_filtered)-nblocks_survived,  bins=len(nblocks)//20, edgecolor='black', color='red', alpha=0.5, label='signal charge', zorder=2)
    
    return [np.mean(np.asarray(nblocks)-np.asarray(nblocks_filtered)), np.std(np.asarray(nblocks)+np.asarray(nblocks_filtered))/np.sqrt(len(nblocks)), np.mean(np.asarray(nblocks)-np.asarray(nblocks_filtered)-nblocks_survived), np.std(np.asarray(nblocks)+np.asarray(nblocks_filtered)+nblocks_survived)/np.sqrt(len(nblocks))]
















    
'''
from scipy.stats import gaussian_kde

 # Calcola la densità di probabilità usando la stima della densità del nucleo (KDE)
kde = gaussian_kde(np.vstack([Dt_S, amplitude_S]))
density = kde(np.vstack([Dt_S, amplitude_S]))

# Definisce una griglia su cui valutare la densità KDE
x = np.linspace(np.asarray(Dt_S).min(), np.asarray(Dt_S).max(), 100)
y = np.linspace(np.asarray(amplitude_S).min(), np.asarray(amplitude_S).max(), 100)
X, Y = np.meshgrid(x, y)
Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)

# Trova il contorno al 95% del livello di confidenza
sorted_Z = np.sort(Z.ravel())[::-1]
cumulative_sum = np.cumsum(sorted_Z)
cumulative_sum /= cumulative_sum[-1]
threshold = sorted_Z[np.searchsorted(cumulative_sum, 0.95)]

# Traccia il grafico con i contorni e scatter plot con cmap
plt.figure(figsize=(10, 8))
plt.contour(X, Y, Z, levels=[threshold], colors='blue', linewidths=2)
sc = plt.scatter(Dt_S, amplitude_S, c=density, s=5, cmap='viridis')
plt.colorbar(sc, label='Density')
plt.title('Contour at 95% CL with Density Colormap')
plt.xlabel(r'$\Delta T$ [s]')
plt.ylabel(r'$\Delta Q$ [pC]')
plt.show()



plt.tight_layout()
# plt.legend(handles=handles, labels=labels)
plt.savefig("distribution_2d.png", dpi=100)
plt.show()
'''








def gaussian(x, mu, std, A):
    return A*np.exp(-(x-mu)**2/(2*std**2))
    


# if mode == 2:
    
#     for i, pmt in enumerate(PMTs):
#         counts1, bins1, _ = ax[i].hist(carica_S[i], bins=len(carica_S[i]), edgecolor='black', color='orange', alpha=0.2, label='signal charge')
#         counts2, bins2, _ = ax[i].hist(carica_noise_S[i], bins=len(carica_noise_S[i]), edgecolor='black', color='blue', alpha=0.2, label='noise charge')
        
#         counts3, bins3, _ = ax[i].hist(carica_bg[i], bins=len(carica_bg[i]), edgecolor='black', color='orange', alpha=0.2, label='background charge')
        
         
#         bin_centers = (bins1[:-1] + bins1[1:]) / 2
#         mu, std = norm.fit(carica_S[i])
#         pdf = norm.pdf(bin_centers, mu, std)
        
#         plt.plot(bin_centers, pdf, 'r-', linewidth=2, label='Fit gaussiano')

        
#         ax[i].set_xlim(0,100)
#         ax[i].set_xlabel('$\Delta Q$ [pC]')
#         ax[i].set_ylabel('Frequency')
#         ax[i].set_title(f'PMT{i} $\Delta Q$ distribution')