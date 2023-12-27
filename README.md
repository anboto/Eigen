# Eigen

[Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a C++ template library for linear algebra, matrices, vectors, numerical solvers and related algorithms.

Eigen is very easy to use. For example, to multiply two matrices:

```
Matrix2d a, b;

a << 1, 2,            // Fill the data
     3, 4;
b << 5, 6,
     7, 8;
Matrix2d res = a*b;    // Just multiply them using *
```

[Eigen U++](https://anboto.github.io/srcdoc$Eigen$Eigen$en-us.html) package is a wrapper of Eigen library, updated to commit fee5d60b (26/12/2023).

examples/Eigen package includes a sample package (Eigen_demo) to ease its use for U++ users. It has many samples from Eigen library and nonlinear equation solving and optimization like [Eckerle4](http://www.itl.nist.gov/div898/strd/nls/data/eckerle4.shtml) and [Thurber](http://www.itl.nist.gov/div898/strd/nls/data/thurber.shtml). To simplify access to these features, simple functions have been added.

It also includes a simple FFT (Fast Fourier Transform) sample that:
- Generates a data series composed by three sinusoidal series of amplitude 2, 5 and 30, and frequencies 1/50, 1/30 and 1/10 Hz:
```
    f(t) = 2*sin(2*PI*t/50 - PI/3) + 5*sin(2*PI*t/30 - PI/2) + 30*sin(2*PI*t/10 - PI/5)
```
- Gets the FFT
- Filters the frequencies between 1/25 and 1/35 Hz
- Gets the filtered data series
- Saves both FFT and series

Eigen packages have been prepared by dolik.rce and koldo.
