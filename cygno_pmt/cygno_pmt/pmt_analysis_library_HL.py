"""
@author: SerAngelo
"""

import cygno as cy 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pmt_analysis_library as pal  
import pickle



def charge_distribution(mfile, corrected, channels_offsets, PMTs, digitizer_type, n_max, event_type, signal_fit = False, resampling=False, show_plots = False, show_distribution1D=False, show_distribution2D=False, verbose=False, save=False):
    
    carica = []
    DEDx = []
    Dt = []
    carica_noise = []
    
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
        
        
        
        if digitizer_type == 'slow':
            pmt_events = [pal.PMT_event(pmt, serial_number, w_slow[pmt.channel], 'slow', event_type, 16) for pmt in PMTs]
        else:
            pmt_events = [pal.PMT_event(pmt, serial_number, w_fast[pmt.channel], 'fast', event_type, 16) for pmt in PMTs]
        
        
        if show_plots:
            fig, ax = plt.subplots(1, 2, figsize=(20, 10))
            c = [pmt_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=ax[i]) for i,pmt_event in enumerate(pmt_events)]    
            plt.tight_layout()
            plt.show()
            
        else:
            c = [pmt_event.PMT_waveforms_integral_info(resampling=resampling, fit=signal_fit, verbose=verbose, axi=None) for pmt_event in pmt_events]
        
        c = list(zip(*c))
        c = [list(i) for i in c]    
            
        if event_type != 'background':
            if -1 in c[0]:
                print(f"Il numero di picchi dell'evento {serial_number} è diverso da 1")
                continue
        
            if any(abs(pmt_events[0].pmt_max_ampl() - pmt_event.pmt_max_ampl()) > 25 for pmt_event in pmt_events[1:]):
                print(f"I {len(PMTs)} picchi dell'evento {serial_number} non sono molto coincidenti")
                continue
                
    
        carica.append(c[0])
        DEDx.append(c[1])
        Dt.append(c[2])
        carica_noise.append(c[3])
      
    carica = list(zip(*carica))
    DEDx = list(zip(*DEDx))
    Dt = list(zip(*Dt))
    carica_noise = list(zip(*carica_noise))
    
    
    if show_distribution1D:        
   
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        
        for i, pmt in enumerate(PMTs):
            ax[i].hist(carica[i], bins=int(len(carica[i])), edgecolor='black', color='orange', alpha=0.2, label=f'{pmt.name} signal charge')
            ax[i].hist(carica_noise[i], bins=int(len(carica_noise[i])), edgecolor='black', color='blue', alpha=0.2, label=f'{pmt.name} noise charge')
            
            ax[i].set_xlim(0,100)
            ax[i].set_xlabel('$\Delta Q$ [pC]')
            ax[i].set_ylabel('Frequency')
            ax[i].set_title(f'PMT{i} $\Delta Q$ distribution')
            
        plt.tight_layout()
        plt.legend()
        plt.show()
    
    
    if show_distribution2D:        
        #colors=['orange', 'green']
        plt.title('Carica rilasciata')
        
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        
        for i, pmt in enumerate(PMTs):
            ax[i].hist2d(Dt[i], carica[i], bins=int(len(carica[i])), cmap='viridis')
            ax[i].set_xlabel('$\Delta T$ [s]')
            ax[i].set_ylabel('$\Delta Q$ [pC]')
            ax[i].set_title(f'PMT{i} $\Delta Q vs \Delta T$')
            
        plt.tight_layout()
        plt.legend()
        plt.show()
    
    if save:
    
        output_file = event_type+".pkl"
        with open(output_file, 'wb') as file:
            pickle.dump([carica, DEDx, Dt, carica_noise] , file)
            
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
            
           # fig.savefig(f'{event_type}.png', dpi=600)
            
            
            break
        
    return [c1,c2]
