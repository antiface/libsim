# -*- coding: utf-8 -*-

### libsim - My Simulation Library
### ==============================
###
### Copyright © 2009, James North
### 
### Permission is hereby granted, free of charge, to any person
### obtaining a copy of this software and associated documentation
### files (the "Software"), to deal in the Software without
### restriction, including without limitation the rights to use,
### copy, modify, merge, publish, distribute, sublicense, and/or sell
### copies of the Software, and to permit persons to whom the
### Software is furnished to do so, subject to the following
### conditions:
### 
### The above copyright notice and this permission notice shall be
### included in all copies or substantial portions of the Software.
### 
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
### EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
### OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
### NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
### HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
### WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
### FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
### OTHER DEALINGS IN THE SOFTWARE.
### 


# I want the new division feature of python
from __future__ import division

# Import lots of libraries
import numpy as np
from numpy import r_
from matplotlib import pyplot
from scipy import integrate, fft
from scipy.fftpack import fftshift, fftfreq
from matplotlib.mlab import find
import math
import operator

def generateSineWave(amplitude, frequency, phase, samplerate, samples):
    """ Generates a sine wave """
    def generateSample(index):
        theta = 2.0*np.pi*(frequency /samplerate)
        sample = math.sin((theta * index) + phase)
        return sample
    # Create the x array
    x = r_[0:samples:1]
    y = np.zeros(samples)
    # Generate the samples
    for n in x:
        y[n] = amplitude * generateSample(n)
    return y

def generateTrianleWave(amplitude, frequency, phase, samplerate, samples):
    """ Generates a triangle wave """
    x = r_[0:samples:1]
    theta = 2.0*np.pi*(frequency / samplerate)
    y = amplitude * np.arcsin(np.sin(theta * x))
    return y

def generateSawtoothWave(amplitude, frequency, phase, samplerate, samples):
    x = r_[0:samples:1]
    y = x
    theta = 2.0*np.pi*(frequency / samplerate)
    for n in x:
        tmp = (n / theta) + phase
        y[n] = amplitude * (tmp - math.floor(tmp))
    return y

def generateSquareWave(amplitude, frequency, phase, samplerate, samples):
    """ Generates a square wave """
    x = r_[0:samples:1]
    theta = 2.0*np.pi*(frequency / samplerate)
    y = np.sin(theta * x)
    for n in range(0, samples, 1):
        y[n] = amplitude * sgn(y[n])
    
    return y

def generateGausianNoise(level, samples):
    pass

def computeSpectrum(signal):
    """ Computes the spectrum of a given signal """
    temp = []
    window = np.hamming(len(signal))
    for n in range(0, len(signal), 1):
        temp.append(signal[n] * window[n])
    
    spectrum = abs(fft(temp))
    
    return spectrum

def detectZeroCrossings(samples, interval=(0,None), detect='all', amplitude=1):
    """
    Detects zero crossings in the array passed to it.
    samples -  The samples you wish to detect
    samplesize - The number of samples
    amplitude - Value to use to indicate a crossing
    """
    zeroCrossings = []
    last = 0
    
    start = interval[0]
    if interval[1] == None:
        end = len(samples)
    else:
        end = interval[1]
    
    for n in range(0, len(samples), 1):
        zeroCrossings.append(0)
    
    for n in range(start, end, 1):
        curr = samples[n]
        if ((last <= 0 and curr > 0)) and \
           detect in ('all', 'posedge', 'pos'):
            zeroCrossings[n] = amplitude
        elif ((last > 0 and curr <= 0)) and \
             detect in ('all', 'negedge', 'neg'):
            zeroCrossings[n] = amplitude
        else:
            zeroCrossings[n] = 0
        last = curr
    
    return zeroCrossings

def measureFrequency(zerocrossings, samplerate):
    """ Calculates the frequency of a waveform from it's zero crossings """
    crossingsCounter = 0
    widthAccumulator = 0
    widthCounter = 0
    crossings = []
    found = False
    
    # Run through the crossings buffer looking for
    # zero crossing events and measure the width
    # between them.
    for sample in zerocrossings:
        if sample >= 1:
            found = True
            if widthCounter > 1:
                crossingsCounter = crossingsCounter + 1
                
                crossings.append(widthCounter)
                widthAccumulator += widthCounter
                widthCounter = 0
            else:
                widthCounter = 0
        elif found == True:
            widthCounter = widthCounter + 1
        else:
            pass
    
    # Average the crossing widths
    if crossingsCounter > 0:
        widthAverage = float(widthAccumulator) / float(crossingsCounter)
    else:
        widthAverage = widthAccumulator / 1
    
    # Calculate the frequency
    if widthAverage > 0:
        frequency = 1.0 / ((1.0/samplerate) * (widthAverage))
    else:
        frequency = 0
    
    return frequency

def DegreeToRadian( degress ):
    """ Converts degress to radians """
    return degress * (np.pi/180.0)

def RadianToDegree( radians ):
    """ Converts radians to degrees """
    return radians * (180.0/np.pi)

def sgn(num):
    if num > 0:
        return 1
    elif num == 0:
        return 0
    elif num < 0:
        return -1

def computeRMS(samples):
    """ Computes the RMS value of the input signal. """
    # Detect the zero crossings of the signal
    crossings = detectZeroCrossings(samples, interval=(0, None), detect='posedge', amplitude=max(samples))
    # Make sure we get several cycles of the sine wave
    f = find(crossings==max(samples))
    a = len(crossings[:f[1]])
    b = len(crossings[:f[len(f)-1:][0]])
    # Create a smaller buffer with only the cycles defined above
    buff = samples[a:b]
    # Calculate the RMS
    rms = math.sqrt((1 / len(buff)) * integrate.trapz([y**2 for y in buff], dx=1))
    return rms

def lineIntersection(lineA, lineB):
    pass

if __name__ == "__main__":
    print "libsim - My Simulation Library"
    print "Copyright © 2009, James North\n\n"
    
    # Create a range to work from, this is mostly for the plotting
    x = r_[0:3000:1]
    
    print "Generating a triangle wave...."
    triwave = generateTrianleWave(1, 50, 0, 44100, 3000)
    
    print "Generating a square wave....."
    squwave = generateSquareWave(1, 50, 0, 44100, 3000)
    
    print "Generating a sine wave...."
    sinwave = generateSineWave(1, 50, 0, 44100, 3000)
    
    print "Computing the spectrum of thr square wave"
    squspec = computeSpectrum(squwave)
    
    print "Computing the spectrum of the trianle"
    trispec = computeSpectrum(triwave)
    
    print "Computing the spectrum of the sine wave"
    sinspec = computeSpectrum(sinwave)
    
    print "Generating Sawtooth Wave......."
    sawwave = generateSawtoothWave(2, 50, 0, 44100, 3000)
    
    print "Plotting the results....."
    
    pyplot.subplot(4,1,1)
    pyplot.plot(x, triwave, label="Trianlge Wave")
    pyplot.plot(x, squwave, label="Square Wave")
    pyplot.plot(x, sinwave, label="Sine Wave")
    #pyplot.plot(x, sawwave, label="Sawtooth Wave")
    pyplot.grid(True)
    pyplot.legend()
    
    pyplot.subplot(4,1,2)
    pyplot.plot(squspec[0:512], label="Square Wave FFT")
    pyplot.grid(True)
    pyplot.legend()
    
    pyplot.subplot(4,1,3)
    pyplot.plot(trispec[0:512], label="Triangle Wave FFT")
    pyplot.grid(True)
    pyplot.legend()
    
    pyplot.subplot(4,1,4)
    pyplot.plot(sinspec[0:512], label="Sine Wave FFT")
    pyplot.grid(True)
    pyplot.legend()
    
    
    pyplot.show()
