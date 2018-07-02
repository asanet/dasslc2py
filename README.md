# Python wrapper for DASSLC

This is a simple wrapper for using the "Differential-Algebraic System Solver in C" by Argimiro R. Secchi (PEQ/COPPE/UFRJ) in python.

## Getting Started

All the information about the usage of this package can be found in example file
```
Dasslc2py/dasslc_example.py
```

### Prerequisites

* C compiler:
```
archlinux: $ sudo pacman -S gcc gcc-libs
```

* Python 2 with numpy and (optionally) matlplotlib:
```
archlinux: $ sudo pacman -S python2 python2-numpy python2-matplotlib
```

* **OR** Python 3 also with numpy and (optionally) matplotlib
```
archlinux: $ sudo pacman -S python3 python3-numpy python3-matplotlib
```

    OBS: matplotlib is needed by the **dasslc_example.py** file

* If you want sparse algebra support, you need to first compile Sparse. On linux, open a terminal in **sparse/src** directory and run:
```
$ make
```


### Installing

* In order to locally build this module, open a terminal at the dir **Dasslc2py** and run the following command
```
$ python setup.py build_ext --inplace
```

* Instead, if you want to install it along with your python distribution, run
```
$ sudo python setup.py install
```

    OBS: if you have multiple versions of python installed, use python2 or python3 as you wish


## Author

* **Ataíde Neto** - ataide@peq.coppe.ufrj.br

