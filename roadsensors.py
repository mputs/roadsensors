import numpy as np
import logging
import sys
import copy

#logging.basicConfig(filename="roadsensors.log", filemode = "a", level = logging.INFO);

def _blockiness(vals):
	##add start and end 1        
	values = copy.copy(vals);
	values.insert(0,1);
	values.append(1);        
	##get indices for zero and positive values        
	ser1 = [ i for i, x in enumerate(values) if x > -1 ];
	##get series with differences above 1
	ser2a = [(x - 1) for x in np.diff(ser1) if x > 1];
	ser2b = [x for x in np.diff(ser1) if x > 1];
	##Calculate B score        
	B = 0.5*sum([i*j for i,j in zip(ser2a,ser2b)]);             
	return(B);

class roadsensorinput:
	"""
		class roadsensorinput.
		klasse voor datapreparatie voor de filter.
		voorbeeld:	
			input = roadsensorinput(data);
			if not input.error:
				indicators = input.Quality();
	"""
	def __init__(self, measurements):
		self.error = False;
		try:
			if len(measurements)<1440:
				self.error = True;
				logging.error("roadsensensorinput: verwacht lijst met 1440 waarnemingen!");
			else:
				self.measurements = copy.copy(measurements);
		except:
			self.error = True;
			logging.error("roadsensorinput: Bij initialisatie hoort lijst met metingen te worden meegegeven!");
	def initialized(self):
		try:
			if len(self.measurements) == 1440:
				return True
			else:
				return False
		except:
			return False;

	def Quality(self):
		"""
			Berekent kwaliteitsindicatoren
			output: dict object met: 
				* indO: aantal nulmetingen in de dagreeks
				* indL: aantal non-missing waarden in de dagreeks
				* indM: gemiddelde waarde
				* indD: gemiddelde absolute verschil met voorgaande waarden
				* indB: blockiness indicator 
			voorbeeld:
				q = input.Quality();
				print "blockiness = ", q.["indB"];
		"""
		self.error = False;
		try:
			vals  = [int(x) for x in self.measurements];
		except:
			logging.error("roadsensorinput: object roadsensorinput niet juist geinitialiseerd!")
			self.error = True;
			return {}
		try:
			ind = {}
			indO = len([x for x in vals if x == 0]);
			indL = len([x for x in vals if x > -1]);
			indM = -1;
			indD = -1;
			if indL > 0:                             
				indM = round(np.mean([x for x in vals if x > -1]), 4);
				indD = round(np.mean(abs(np.diff([x for x in vals if x > -1]))), 4);        
			indB = _blockiness(vals)
			ind["indO"] = indO;
			ind["indL"] = indL;
			ind["indM"] = indM;
			ind["indD"] = indD;
			ind["indB"] = indB;
			return ind;	
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens Quality: %s %s %s", sys.exc_info()[0],sys.exc_info()[1], sys.exc_info()[2])


class roadsensoroutput:
	"""
		class roadsensoroutput:
		klasse voor het bewaren en comprimeren van het resultaat van roadsensorfilter.
		het gefilterde signaal staat opgeslagen in roadsensoroutput.result

		methoden:
			fftcompress: comprimeer de result met "cutoff" frequenties
			uncompress: decomprimeer het resultaat
			getcompressed: haal de gecomprimeerde parameters op
			setcompressed: zet de gecomprimeerde data voor decompressie (zie uncompress)

		voorbeeld:
			input = roadsensorinput([1,1,2,3,2,...,0]);
			f = roadsensorfilter(200, 50);
			output = f.filterdata(input);
			output.compress(12);
			y1 = nieuw.uncompressed();
			p = output.getcompressed();
			nieuw = roadsensoroutput();
			nieuw.setcompressed(p);
			y2 = nieuw.uncompressed();

			y1 en y2 moeten gelijk zijn.
	"""
	def __init__(self):
		self.result = np.zeros((1440));

	def fftcompress(self, cutoff):
		"""
			fftcompress: comprimeer de data door alleen naar de laagste frequenties in het fourierdomein te behouden
			parameter cutoff: het aantal frequentiecomponenten dat moet worden bewaard.
			let op: dit is een lossy compression!
		"""
		self.error = False;
		try:
			ft = np.fft.fft(self.result);
			self.compressed=np.copy(ft[0:cutoff+1]);
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens fftcompress: %s", sys.exc_info()[0])
	def uncompress(self):
		"""
			uncompress: decomprimeer de resultaten van fftcompress.a
			geeft de gedecomprimeerde reeks terug.
		"""
		self.error = False;
		try:
			cutoff = len(self.compressed)-1;
			ft = np.zeros(1440, dtype=np.cfloat);
			ft[0:cutoff+1] = self.compressed
			ft[(1440-cutoff):1440] = np.conj(self.compressed[1:cutoff+1])[::-1]	
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens uncompress: %s", sys.exc_info()[0])
			return [];
		return np.real(np.fft.ifft(ft));
	def getcompressed(self):
		"""
			getcompressed: geeft een lijst met parameters terug (2*cutoff) voor opslag.
		"""
		self.error = False;
		try:
			return [x for x in np.real(self.compressed)]+ [x for x in np.imag(self.compressed)];
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens uncompress: %s", sys.exc_info()[0])
			return[];
	def setcompressed(self, compressed):
		"""
			setcompressed: biedt een lijst met parameters aan voor decompressie.
		"""
		self.error = False;
		try:
			self.compressed = np.copy(compressed[0:len(compressed)/2])+1j*np.copy(compressed[len(compressed)/2:])
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens uncompress: %s", sys.exc_info()[0])

		
		
class roadsensorfilter:
	"""
		class roadsensorfilter: microgaafmaker voor verkeerslusdata
		gebruik: roadsensorfilter(NSamples, maxmeas);
			NSamples: aantal gediscretiseerde punten van de kansdichtheidsfunctie (meer is beter)
			maxmeas: maximaal aantal voertuigen die door een sensor kan worden gemeten
		methoden:
			filterdata(input): voer de filter uit op "input"
		voorbeeld:
			f = roadsensorfilter(200,50);
			output = f.filterdata(input)
			print output.result;
	"""
	def meas2idx(self,meas):
		if meas > self.maxmeas:
			return self.maxmeas+1
		elif meas < -1:
			return 0
		else:
			return meas + 1


	def __init__(self, NSamples, maxmeas):
		self.error = False;
		try:
			self.NSamples = NSamples;
			self.maxmeas = maxmeas;
			if NSamples < 1:
				logging.error("NSamples (eerste parameter) bij initialisatie filterdata moet groter zijn dan nul!");
			if maxmeas <= 0:
				logging.error("maxmeas (tweede parameter) bij initialisatie filterdata moet groter zijn dan nul!");
			if NSamples < 1 or maxmeas <= 0:
				return
			self.rates = np.linspace(0,maxmeas,NSamples);
			self.memberfuns = np.zeros((maxmeas+2,NSamples));
			self.pdf = np.zeros((1440,self.NSamples));
			self.result = np.zeros(1440);
			for meas in range(0, maxmeas+1):
				idx = self.meas2idx(meas);
				self.memberfuns[idx][1:NSamples] = meas * np.log(self.rates[1:NSamples])-self.rates[1:NSamples];
			for meas in range(1, maxmeas+1):
				self.memberfuns[self.meas2idx(meas)][0] = -10000.0;
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens initialisatie filterdata: %s", sys.exc_info()[0])


	def filterdata(self, rsinput):
		"""
			filterdata: methode voor het filteren van de data
			input:
				rsinput: roadsensorinput-object
			return:
				roadsensoroutput-object
		"""
		self.error = False;
		if not rsinput.initialized():
			self.error = True;
			logging.error("functie filterdata heeft een niet geinitialiseerd input-object gekregen!");
			return [];
		
		#estimate
		try:
			for k in range(0,1440):
				idx = self.meas2idx(rsinput.measurements[k]);
				self.pdf[k] = np.copy(self.memberfuns[idx]);
		except IndexError:
			self.error = True;
			logging.error("functie filterdata verwacht een lijst met 1440 metingen i.p.v. %d metingen!", len(measurements))
			return [];
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens filterdata#estimate: %s", sys.exc_info()[0])
			return [];
		#predict
		try:
			for k in range(1,1440):
				self.pdf[k] += self.pdf[k-1]*.6;
		
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens filterdata#predict: %s", sys.exc_info()[0])
			return [];
		#smooth
		try:
			for k in range(1438,-1,-1):
				self.pdf[k] += self.pdf[k+1]*.8;
		except:
			self.error = True;
			logging.error("onverwachte fout tijdens filterdata#smooth: %s", sys.exc_info()[0])
			return [];
		#crisp
		try:
			result = roadsensoroutput();
			for k in range(0,1440):
				result.result[k] = self.rates[np.argmax(self.pdf[k])];
		except:
			self.error = True;
			logging.error("overwachte fout tijdens filterdata#cruso: %s", sys.esc_info()[0])
			return [];
		return result;

