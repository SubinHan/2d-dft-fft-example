# 2d-dft-fft-example

Here is the C code example of 1-d, 2-d Discrete Fourier Transform and Fast Fourier Transform.
The FFT implemented with Cooley-Tukey Algorithm.
For more detail, consider followings: 
[Discrete Fourier Transform Wiki](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
[Fast Fourier Transform Wiki](https://en.wikipedia.org/wiki/Fast_Fourier_transform)

In 2d FFT, It uses the technique as following:
1. do 1D FFT on each row (real to complex)
2. tranpose matrix
3. do 1D FFT on each row resulting from (2) (complex to complex) (That is, 1D FFT on each column)
