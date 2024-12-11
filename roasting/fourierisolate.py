import numpy as np
import librosa
import matplotlib.pyplot as plt

# Load the audio file
cracks = '2 cracks reference.m4a'
blower = '/Users/quinnwilliams/Documents/Code/Fall Dep Study/CA-Math-Dep-Study/roasting/No crack reference.m4a'
y, sr = librosa.load(cracks, sr=None)
x, sr = librosa.load(blower, sr=None)

cracks_fft = np.fft.fft(y)
blower_fft = np.fft.fft(x)

cracks_freq = np.fft.fftfreq(len(cracks_fft), 1/sr)
blower_freq = np.fft.fftfreq(len(blower_fft), 1/sr)

difference = []

for i in range(len(cracks_fft)):
    val = cracks_fft[i] - blower_fft[i]
    difference.append(val)

print(difference)


# Plot the frequency magnitude between 6kHz and 15kHz
plt.figure(figsize=(14, 5))
plt.plot(blower_freq, np.abs(difference))
plt.title('Frequency Spectrum of Crack')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.grid()
plt.show()
