#Stolen from https://github.com/rciurlea/yin-python/blob/master/yin.py

#import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy.fftpack import fft
from scipy.fftpack import ifft
from scipy.io import wavfile
import numpy as np

import matplotlib.pyplot as plt

def yin_raw(fs, data):
    print(fs)
       
    tau_max = min(3000, len(data))
    w_size = 6000
    r = np.zeros(tau_max)

    s = fft(data)
    r = ifft(s*np.conj(s)).real[:tau_max]
    correction = r[0]
    
    # Difference formula, from the paper, has 
    #
    #        Sum[  (y_(t+tau)  - y_t)^2  ]
    #
    # Here, I expand that into
    #
    #        Sum[  (y_(t+tau)  - y_t)^2  ]  = 
    #        Sum[  y_t^2  ] + Sum[  y_(t+tau)^2 ] - 2*Sum[  y_t*y_(t+tau)  ] = 
    #        2*Sum[  y^t^2  ] - 2*Sum[  y_t*y_(t+tau)  ]
    #
    # Notice the last formula is the classic autocorrelation, which can be computed quickly.
    # It's 0th coordinate is the correction term we see above.
    # This gives us a quick-and-dirty way to compute the difference formula.
    #
    # I ignore the unnecessary factor of 2.
    #
    # -- Artem Bolshakov
        
    r = correction - r
    print(r)

    d = np.zeros(tau_max)
    s = r[0]
    d[0] = 1
    for i in range(1,tau_max):
        s += r[i]
        d[i] = r[i] / ((1.0 / i) * s) 
    print(d) 
   
    # find frequency. use 0.5 as threshold
    
    for i in range(1, tau_max):
        if d[i] > 0.3:
            continue
        if d[i-1] > d[i] < d[i+1]:
            return fs/i
            break
    return 999999 

def yin(freqs, powerSpec, threshold = 0.5):
    res = np.zeros(powerSpec.shape[1])
    print(len(freqs))
    print(powerSpec.shape)
    
    for t in range(powerSpec.shape[1]):
    
        r = ifft(powerSpec[:, t]).real
        correction = r[0]
        
        # Difference formula, from the paper, has 
        #
        #        Sum[  (y_(t+tau)  - y_t)^2  ]
        #
        # Here, I expand that into
        #
        #        Sum[  (y_(t+tau)  - y_t)^2  ]  = 
        #        Sum[  y_t^2  ] + Sum[  y_(t+tau)^2 ] - 2*Sum[  y_t*y_(t+tau)  ] = 
        #        2*Sum[  y^t^2  ] - 2*Sum[  y_t*y_(t+tau)  ]
        #
        # Notice the last formula is the classic autocorrelation, which can be computed quickly.
        # It's 0th coordinate is the correction term we see above.
        # This gives us a quick-and-dirty way to compute the difference formula.
        #
        # I ignore the unnecessary factor of 2.
        #
        # -- Artem Bolshakov
        
        r = correction - r
#        print(r)
    
        d = np.zeros(powerSpec.shape[0])
        s = r[0]
        d[0] = 1
        for i in range(1,len(freqs)):
            s += r[i]
            d[i] = r[i] / ((1.0 / i) * s) 
#        print(d) 
       
        # find frequency. use 0.5 as threshold
        
        for i in range(1, len(freqs) - 1):
            if d[i] > threshold:
                continue
            if d[i-1] > d[i] < d[i+1]:
                res[t] = freqs[-1]/float(i)
                break
    return res 

if __name__ == "__main__":
    #dataDir = 'data_parrot1'
    dataDir = '.'
    fn = '100000syll1673.wav'
    
    fs, data = wavfile.read(dataDir + '/' + fn)
    f, t, s = spectrogram(data, fs)
    print(yin_raw(fs, data[20000:20512])) 
    y = yin(f, s)
    print(y)
    plt.plot(t, y)
    plt.show()
    plt.imshow(np.log(s[::-1, :]).clip(min=0))
    plt.show()
    print(f[-1])
    
