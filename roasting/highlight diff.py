import numpy as np
import librosa
import matplotlib.pyplot as plt

# Load the two audio files
audio_path1 = '/Users/quinnwilliams/Documents/Code/Fall Dep Study/CA-Math-Dep-Study/roasting/2 cracks reference.m4a'  # Replace with your actual file paths
audio_path2 = '/Users/quinnwilliams/Documents/Code/Fall Dep Study/CA-Math-Dep-Study/roasting/No crack reference.m4a'

y1, sr1 = librosa.load(audio_path1, sr=None)  # Load first audio file
y2, sr2 = librosa.load(audio_path2, sr=None)  # Load second audio file

# Make sure both files have the same length, trim or pad if necessary
min_len = min(len(y1), len(y2))
y1 = y1[:min_len]  # Trim or pad both signals to the same length
y2 = y2[:min_len]

# Compute the FFT of both audio signals
fft1 = np.fft.fft(y1)
fft2 = np.fft.fft(y2)

# Compute the magnitude of the FFT (we take the absolute value)
magnitude1 = np.abs(fft1)
magnitude2 = np.abs(fft2)

# Compute the difference between the two FFT magnitudes
fft_difference = magnitude1 - magnitude2
fft_difference = np.abs(fft_difference)

# Compute the corresponding frequency values for plotting
frequencies = np.fft.fftfreq(min_len, 1/sr1)

# Plot the difference in FFT magnitudes
plt.figure(figsize=(14, 5))
plt.plot(frequencies[:min_len // 2], fft_difference[:min_len // 2])  # Plot only positive frequencies
plt.title('Difference in FFT Magnitudes of Two Audio Files')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Difference in Magnitude')
plt.grid()
plt.show()
