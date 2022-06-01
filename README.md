# Non-Subsampled_Shearlet_Transform

[<img src="https://img.shields.io/badge/FFTW-3.3.5-76B900?style=for-the-badge" style="vertical-align:top margin:6px 4px">](http://www.fftw.org/install/windows.html)


## What is NSST ?
NSST is a discrete form of Shearlet transform, and it differs from other multi-scale transformations by avoiding up-down samplers. NSST consists of two main stages, multi-scale and multi-directional separations.

### Multi-scale
<ul>
  <li>
     The Laplacian Pyramid (NSLP) without subsampling produces low and high frequency images whose size is the same as the size of the source image.
  </li>
</ul>

### Multi-directional
<ul>
  <li>
     Versatility is achieved by using different "combinations of Shear Filters" in the so-called polar (pseude-polar) coordinate
  </li>
</ul>

## NSST Steps
The process steps performed to obtain the NSST coefficients of an image of NxN size at a fixed resolution scale are as follows:

<ol>
  <li>
    Laplacian pyramid is applied to the image. Low and High pass sub-images are obtained.
  </li>
  <li>
    The fourier transformations of the high pass sub-images are calculated and transformed into the Polar coordinate system.
  </li>
  <li>
    Bandpass filter is applied to Polar coordinate system transformations and Fourier transforms (FFT) of Shearlet coefficients are obtained.
  </li>
  <li>
    The Inverse Fourier Transform (IFT) is applied to obtain the Shearlet coefficients and the transformation is performed to the Cartesian coordinate system.
  </li>
</ol>

## Pipeline

![image 1](https://github.com/fbasatemur/Non-Subsampled_Shearlet_Transform/blob/main/screenshots/NSST_Design.jpg)

<ol>
  <li>
    Firstly, the input image must be converted to the Intensity channel. At this stage, only Y channel will be used for NSST. You can use <b>ConvertBMPToIntensity()</b> for this process. General formula :
     <h3> Y = (0.11 * Red + 0.59 * Green + 0.3 * Blue) </h3>
  </li>
  <li>
    Intensity image given to NSST function. Low (AFK) and High (YFK) Frequency coefficients are obtained at its output. You can use <b>NsstDec1e()</b> for this process.
  </li>
</ol>

You can get the input Intensity image using TNSST (Inverse NSST). You can use **NsstRec1()** for this process.


## Run

### ***Windows:***
```shell
mkdir build && cd build
copy ..\external_lib\fftw-3.3.5-dll64\libfftw3f-3.dll .
cmake ..
make
```

***Or run with Clion***  
**Note**: A copy of ***external_lib/fftw-3.3.5-dll64/libfftw3f-3.dll*** must be in the executable directory

### ***Ubuntu:***
```shell
mkdir build && cd build
cmake ..
make
```

**Note**: If you want to see the Eigen array values during debug, you can copy the **.gbinit** file to ***/home/<user-name>/***

## Dependencies
***Don't worry, no setup needed !***
- Eigen 3.4.0
- FFTW 3.3.5

## Thanks !
<ul>
  <li>
    <h2><a href="https://scholar.google.com/citations?hl=tr&user=Mq8UBzQAAAAJ" target="_blank">Assoc. Prof. Hülya Doğan</a></h2>
  </li>
</ul>
