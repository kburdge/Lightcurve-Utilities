import cuvarbase.bls as bls
import cuvarbase.lombscargle as gls
import numpy as np

def BLS(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.12,remove=True):

	t=t-np.mean(t)
	y=y
	dy=dy

# generate data with a transit


# set up search parameters
	search_params = dict(qmin=qmin, qmax=qmax,
                     # The logarithmic spacing of q
         	            dlogq=0.1,

                     # Number of overlapping phase bins
                     # to use for finding the best phi0
     	                noverlap=3)
	
# derive baseline from the data for consistency
	baseline = max(t) - min(t)

# df ~ qmin / baseline
	df = search_params['qmin'] / baseline
	if pmax:
		fmin = 4/baseline
	else:
		fmin=1/pmax

	fmax = 1440.0/pmin

	nf = int(np.ceil((fmax - fmin) / df))
	freqs = fmin + df * np.arange(nf)
	if remove:
		low=142
		high=146
		freqs_to_remove = [[low,high], [low/2,high/2],[low/3,high/3], [low/4,high/4], [low/5,high/5], [low/6,high/6], [low/7,high/7],[low/8,high/8], [low/9,high/9]]
		#freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
		for pair in freqs_to_remove:
			idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
			freqs = freqs[idx]

	bls_power = bls.eebls_gpu_fast(t, y, dy, freqs,
                     	           **search_params)


	best_index = np.argmax(bls_power)
	best_freq=freqs[best_index]
	freqs_to_remove = [[best_freq-0.1*best_freq, best_freq+0.1*best_freq]]
	for pair in freqs_to_remove:
		idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
		ls_power2 = bls_power[idx]
		freqs2 = freqs[idx]
		
	bls_power_best=(np.max(bls_power)-np.mean(bls_power2))/np.std(bls_power2)
	f_best = freqs[np.argmax(bls_power)]
	#bls_power_best=(np.max(bls_power)-np.median(bls_power))/(np.std(bls_power))
	period=1.0/f_best


	return t, y, dy, period, bls_power_best

def BLS_Full(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.12,remove=True,trim=True):

	t=t-np.mean(t)
	y=y
	dy=dy

# generate data with a transit


# set up search parameters
	search_params = dict(qmin=qmin, qmax=qmax,
                     # The logarithmic spacing of q
         	            dlogq=0.1,

                     # Number of overlapping phase bins
                     # to use for finding the best phi0
     	                noverlap=3)
	
# derive baseline from the data for consistency
	baseline = max(t) - min(t)

# df ~ qmin / baseline
	df = search_params['qmin'] / baseline
	if pmax:
		fmin = 4/baseline
	else:
		fmin=1/pmax

	fmax = 1440.0/pmin

	
	nf = int(np.ceil((fmax - fmin) / df))
	freqs = fmin + df * np.arange(nf)
	if remove:
		low=142
		high=146
		freqs_to_remove = [[low,high], [low/2,high/2],[low/3,high/3], [low/4,high/4], [low/5,high/5], [low/6,high/6], [low/7,high/7],[low/8,high/8], [low/9,high/9]]
		#freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
		for pair in freqs_to_remove:
			idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
			freqs = freqs[idx]

	bls_power = bls.eebls_gpu_fast(t, y, dy, freqs,
                     	           **search_params)
	f_best = freqs[np.argmax(bls_power)]             	           
	if trim:       

		best_freq=freqs[np.argmax(bls_power)]
		freqs_to_remove = [[best_freq-0.1*best_freq, best_freq+0.1*best_freq]]
		for pair in freqs_to_remove:
			idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
			bls_power2 = bls_power[idx]
			freqs2 = freqs[idx]
			


		
					
		bls_power_best=(np.max(bls_power)-np.median(bls_power2))/np.std(bls_power2)
	else:


		bls_power_best=(np.max(bls_power)-np.median(bls_power))/(np.std(bls_power))
	period=1.0/f_best


	return t, y, dy, period, bls_power_best, freqs, bls_power

def LS(t,y,dy,pmin=2):
	t=t-np.mean(t)
	lightcurves=[(t,y,dy)]
	baseline=np.max(t)-np.min(t)

	df = 1.0 / (baseline*3.0)
	fmin = 4.0/baseline
	fmax = 1440/pmin

	nf = int(np.ceil((fmax - fmin) / df))
	freqs=np.linspace(fmin,fmax,nf)

	proc = gls.LombScargleAsyncProcess(use_double=True)
	result = proc.run([(t, y, dy)],freqs=freqs,use_fft=True)
	proc.finish()
	freqs, ls_power = result[0]
	freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
	for pair in freqs_to_remove:
		idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
		ls_power = ls_power[idx]
		freqs = freqs[idx]
				
	significance=(np.max(ls_power)-np.mean(ls_power))/np.std(ls_power)

	period=1.0/freqs[np.argmax(ls_power)]

	return t, y, dy, period, significance
	
def LS_Full(t,y,dy,pmin=2,oversample_factor=3.0,trim=True):
	t=t-np.mean(t)
	lightcurves=[(t,y,dy)]
	baseline=np.max(t)-np.min(t)

	df = 1.0 / (baseline*oversample_factor)
	fmin = 4.0/baseline
	fmax = 1440/pmin

	nf = int(np.ceil((fmax - fmin) / df))
	freqs=np.linspace(fmin,fmax,nf)

	proc = gls.LombScargleAsyncProcess(use_double=True)
	result = proc.run([(t, y, dy)],freqs=freqs,use_fft=True)
	proc.finish()
	freqs, ls_power = result[0]
	#freqs_to_remove = [[3e-2,4e-2], [49.99,50.01], [48.99,49.01], [47.99,48.01], [46.99,47.01], [45.99,46.01], [3.95,4.05], [2.95,3.05], [1.95,2.05], [0.95,1.05], [0.48, 0.52], [0.32, 0.34], [0.24, 0.26], [0.19, 0.21]]
	#for pair in freqs_to_remove:
	#	idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
	#	ls_power = ls_power[idx]
	#	freqs = freqs[idx]
	
	best_freq=freqs[np.argmax(ls_power)]
	if trim:
		freqs_to_remove = [[best_freq-0.1*best_freq, best_freq+0.1*best_freq]]
		for pair in freqs_to_remove:
			idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
			ls_power2 = ls_power[idx]
			freqs2 = freqs[idx]
			

		f_best = freqs[np.argmax(ls_power)]
		
					
		significance=(np.max(ls_power)-np.mean(ls_power2))/np.std(ls_power2)
	else:
	
		f_best = freqs[np.argmax(ls_power)]
		
					
		significance=(np.max(ls_power)-np.mean(ls_power))/np.std(ls_power)	

	period=1.0/freqs[np.argmax(ls_power)]

	return t, y, dy, period, significance, freqs, ls_power

