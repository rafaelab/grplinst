# grplinst
Module for the CRPropa code to calculate the effect of plasma instabilities on the development of electromagnetic cascades.
It is described in:

R. Alves Batista, A. Saveliev, E. M. de Gouveia Dal Pino. [MNRAS 489 (2019) 3836](https://academic.oup.com/mnras/article/489/3/3836/5558261). [arXiv:1904.13345](https://arxiv.org/abs/1904.13345)

If you make use of grplinst, please consider citing the above paper.


## Science

This code was designed to study the effects of plasma instabilities on the development of blazar-induced gamma-ray cascades in the intergalactic medium. 

It provides *approximate* descriptions for some plasma instability models found in the literature. 
They are modelled as a cooling term for electrons/positrons. Note that this is a very rough approximation and a detailed calculation would require particle-in-cell (PIC) simulations.



## Installation procedure

To install *grplinst*, you will need to have CRPropa 3 installed. 
Go to https://github.com/CRPropa/CRPropa3/ and follow the instructions.

Now proceed with the installation of this module.

1. Download the latest version of the code.
```
git clone https://github.com/rafaelab/grplinst.git
```

2. Navigate into the downloaded folder and create a folder called "build/".

3. Install the code with CMake and then make:

```
cmake ..
make
```

4. If the code was properly compiled, you are ready to go!
Make sure to add the path where grplinst.py is created to your PYTHONPATH.
Alternatively, this can be hard-coded in your python script.

For consistency with CRPropa, the flag `-DCMAKE_CXX_FLAGS="-std=c++11"` may be required to force compilers to adopt C++11.


## Disclaimer
This program is provided 'as is', without warranties of any kind. 
Please use your discernement to interpret the results obtained with it.

## Acknowledgements
We thank Michael Kachelries and Jacob Benestad for providing feedback that lead us to finding a few bugs. 



