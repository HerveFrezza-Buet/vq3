# vq3

The vq3 library is a C++ implementation of vector quantization, including GNG-T, using C++17 standard. It enables multi-thread acceleration.

Two packages are provided. vq3 is the core library. It is generic. vq3demo provides demo examples using openCV.


# Unix Installation

First, get the files.

``` 
git clone https://github.com/HerveFrezza-Buet/vq3
``` 

Then, you can install the two packages as follows. The commands below concern the installation of the vq3 package, the installation of vq3demo is similar. 

```
mkdir -p vq3/build
cd vq3/build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr
sudo make install
cd ../..
```

The documentation is in /usr/share/doc/vq3

