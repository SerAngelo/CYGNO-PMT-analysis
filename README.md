# CYGNO-PMT-analysis
Benvenuto! 
Tutta l'analisi fa riferimento a due librerie, una più a basso livello chiamata "pmt_analysis_library.py" ed una più ad alto livello "pmt_analysis_HL.py". Per visualizzare le waveforms bisogna fare uso della funzione "single_waveform_info" definita in "pmt_analysis_HL.py". 
Un esempio del suo utilizzo è stato inserito in "esempio.py". Tuttavia, in seguito metto una breve guida all'utilizzo.
```
pmt_analysis_HL.single_waveform_info(mfile, corrected, channels_offsets, PMTs, n, digitizer_type, event_type, signal_fit=False, resampling=False, verbose=False)
```
PARAMETRI:
- mfile: è l'oggetto generato attraverso la libreria cygno
- corrected e channels_offsets: sono estraibili da mfile
- PMTs: lista di oggetti PMT (definiti in pmt_analysis_library.py)
- n: numero dell'evento da visualizzare
- digitizer_type: 'fast' o 'slow'
- event_type: nome con cui chiamare l'immagine (che verrà automaticamente salvata)
- signal_fit: esegue il fit e lo fa visualizzare
- resampling: attiva il resampling
- verbose: stampa varie informazioni

RETURN:
-Restituisce ~[c1,c2]~ dove ~c_i = [carica, carica/width, width, carica_rumore, estremo_sinistro , estremo_destro, posizione_picchi, numero_picchi, ampiezza]~


