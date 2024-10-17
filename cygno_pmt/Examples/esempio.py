import cygno as cy
import pmt_analysis_library as pal
import pmt_analysis_library_HL as palhl


PMT1 = pal.PMT('PMT1', 'Hamamatsu-R1894', 7.8e-9, 2)#nome,modello,costanteRC,canale
PMT2 = pal.PMT('PMT2', 'Hamamatsu-R1894', 7.8e-9, 3)
PMTs = [PMT1, PMT2]

mfile = cy.open_mid(run=12398, path='percorso', cloud=True, tag='', verbose=False)
odb = cy.get_bor_odb(mfile)
corrected  = odb.data['Configurations']['DRS4Correction']
channels_offsets  = odb.data['Configurations']['DigitizerOffset']

palhl.single_waveform_info(mfile, corrected, channels_offsets, PMTs, 10, 'slow', 'nome_immagine', signal_fit=False, resampling=False, verbose=True)
