# vq3

The vq3 library is a C++ implementation of vector quantization, including GNG-T, using C++17 standard. It enables multi-thread acceleration.

Two packages are provided. vq3 is the core library. It is generic. vq3demo provides demo examples using openCV.

# Support

vq3 is a result of the <a href="http://interreg-grone.eu">GRONE project</a>, supported by the Interreg "Grande Région" program of the European Union's European Regional Development Fund.


# Unix Installation

First, get the files.

``` 
git clone https://github.com/HerveFrezza-Buet/vq3
``` 

Then, you can install the two packages as follows. The commands below concern the installation of the vq3 package, the installation of vq3demo is similar. You need to have opencv and demo2d installed to build vq3demo. The demo2d package is also available on my github.

```
mkdir -p vq3/build
cd vq3/build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr
sudo make install
cd ../..
```

The documentation is in /usr/share/doc/vq3

# Demos

You can run the examples of the vq3demo package. Some demos realized by vq3demo are also visible on <a href="https://frezza.pages.centralesupelec.fr/hfb-com/vq/vq3.html">my web pages</a>.


