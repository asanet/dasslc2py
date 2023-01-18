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

Python 3 with numpy, pip and (optionally) matplotlib
```
archlinux: $ sudo pacman -S python python-numpy python-matplotlib
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
$ sudo pip install .
```

## Author

* **Ataíde Neto** - ataide@peq.coppe.ufrj.br

