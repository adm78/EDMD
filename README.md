# EDMD
## Event Driven Molecular Dynamics with Python

This is a Pythonic version of the JavaScript EDMD code that I wrote for the [Visual Chemical Engineering (VCE) project](http://visualchemeng.com/). 

The JavaScript version is hosted on GitHub [here](https://github.com/adm78/visualchemeng_js/blob/master/modules/md/md.js). 

![edmd example image](https://user-images.githubusercontent.com/17439476/33520943-7bb99ff8-d7bc-11e7-9cd6-9e676f3a297d.png)

## Dependencies
- Python 2.*
- Pygame
- Numpy

## Usage
Default conditions
```shell
./md.py
```
Simulate 60 particles within a 500 x 500px simulation box and a maximum time step of size 10.0
```shell
./md.py -n 60 -x 500 -y 500 --dt 10.0
```
Display help and exit
```shell
./md.py --help
```
