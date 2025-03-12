## Digit Identification Through Audio Signals
*in MATLAB*

&nbsp;

### 📋 Project Description
---

*Data Description*

The dataset for this project consists of speech signals recorded from 60 participants (each subdirectory in the "data" folder corresponds to one participant). Each participant repeated each digit 50 times, meaning that each of the 60 folders contains 500 audio signals in .wav format. Each audio signal was recorded at a sampling rate of 48,000 Hz in mono-channel mode. Source: *https://www.kaggle.com/datasets/sripaadsrinivasan/audio-mnist*

&nbsp;

*Project Guidelines*

1. Download the dataset.
   
2. Choose the participant (folder "40").
   
3. Develop code to import the signals.
   
4. Reproduce/graphically represent an example of the imported signals, identifying the digit each corresponds to.
   
   4.1. Visually identify possible temporal characteristics that allow differentiating pairs of digits.
   
   4.2. Identify temporal characteristics, such as energy, max amplitude, .., that allow differentiating all digits.
   
   4.3. Formulate possible decision rules (e.g., if-then-else) that best separate the different digits in a feature space of up to three dimensions. Represent the decision surfaces.
   
6. Based on the assigned signals, compute, for each digit, the median and normalized amplitude spectrum by the number of samples (i.e., equivalent to the modulus of the complex Fourier series coefficients) only for positive frequencies. Also, calculate the first quartile (25%) and third quartile (75%).
   
   5.1. Compare four different types of windows.
   
   5.2. Identify possible spectral characteristics that allow differentiating pairs of digits.
   
   5.3. Identify spectral characteristics that allow differentiating all digits.
   
   5.4. Formulate possible decision rules that best separate the different digits in a feature space of up to three dimensions. Represent the decision surfaces.
   
7. Repeat the previous steps using the STFT. Use different parameterizations, such as the number of points for FFT computation, overlap points, etc., and identify the one(s) that seem best suited for the objective of the work.
   
8. Consider the best features from time domain analysis, DFT, and STFT.
   
   7.1. Formulate decision rules that allow separating digit pairs in a feature space of up to three dimensions.
    
   7.2. Formulate possible decision rules that allow separating all digits in a feature space of up to three dimensions.

&nbsp;

---

**NOTE:** This project provided insights into how speech signals can be analyzed and classified based on their time and frequency domain characteristics. However, the development remains incomplete, as points 6 and 7 have not been implemented. This project was developed as part of a Computer Science course by Cláudia Torres, André Louro and Maria João Rosa.
