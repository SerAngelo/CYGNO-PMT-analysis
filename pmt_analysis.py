"""
@author: SerAngelo
"""

import cygno as cy 
import pmt_analysis_library as pal  
import pmt_analysis_library_HL as palhl

# =============================================================================
# CARATTERIZZAZIONE PMT ESPERIMENTO CYGNO
# =============================================================================

# =============================================================================
# Modello: PMT Hamamatsu R1894
# =============================================================================


# =============================================================================
# TODO-LIST:
#     > DETERMINARE I BOUNDARY -------> ok
#     > MISURA DEL NOISE STATISTICO ------>ok
#     > PULSE SHAPE DISCRIMINATION (MUONE VS 55Fe)
#        -fit dello spettro -------> ok
#     > MISURA DELL'ENERGIA RILASCIATA (integrale) ------> ok
# =============================================================================


        

PMT1 = pal.PMT('PMT1', 'Hamamatsu-R1894', 7.8e-9, 2)
PMT2 = pal.PMT('PMT2', 'Hamamatsu-R1894', 7.8e-9, 3) 

PMTs = [PMT1, PMT2]


mfile = cy.open_mid(run=12400, path='/home/angelo/Scrivania/Università/MAGISTRALE/LAB HEP 2/dati/', cloud=False, tag='', verbose=False)
odb = cy.get_bor_odb(mfile)
corrected  = odb.data['Configurations']['DRS4Correction']
channels_offsets  = odb.data['Configurations']['DigitizerOffset'] 

palhl.single_waveform_info(mfile, corrected, channels_offsets, PMTs, 9, 'slow', 'source', signal_fit=True, resampling=False, verbose=True)

palhl.charge_distribution(mfile, corrected, channels_offsets, PMTs, 'slow', 10, 'source', signal_fit = False, resampling=False, show_plots=True, show_distribution1D=True, save=False)

