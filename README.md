# Non-Subsampled_Shearlet_Transform

:white_check_mark: Faster NSST is ready for use. Necessary explanations will be made soon...

[<img src="https://img.shields.io/badge/FFTW-3.3.5-76B900?style=for-the-badge" style="vertical-align:top margin:6px 4px">](http://www.fftw.org/install/windows.html)

## Run

The FFTW3 library in external_lib is compiled for Windows. For Ubuntu, FFTW 3.3.5 should be installed.  
Run steps for Windows: 
```shell
mkdir build && cd build
copy ..\external_lib\fftw-3.3.5-dll64\libfftw3f-3.dll .
cmake ..
make
```

***Or run with Clion***  
**Note**: A copy of ***external_lib/fftw-3.3.5-dll64/libfftw3f-3.dll*** must be in the executable directory

![alt text](https://github.com/fbasatemur/Non-Subsampled_Shearlet_Transform/blob/develop/screenshots/NSST_Design.jpeg)
