import librosa
import matplotlib.pyplot as plt

# Load audio file
audio_path = '/Users/quinnwilliams/Documents/Code/Fall Dep Study/CA-Math-Dep-Study/roasting/2 cracks reference.m4a'
y, sr = librosa.load(audio_path, sr=None)

# Plot the waveform
plt.figure(figsize=(14, 5))
plt.plot(y)
plt.title('Waveform of the Audio')
plt.xlabel('Time (samples)')
plt.ylabel('Amplitude')
plt.show()
