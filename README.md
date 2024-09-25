# digitizationpp

** work in progress **

Digitization code in C++

## Dependencies

* ROOT [compiled with the C++17 or C++20 standard]
* ROOTANA [https://bitbucket.org/tmidas/rootana/src/master/]
* OPENCV [https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html]

Before compiling, set the variables ROOTANASYS and OPENCVSYS in your environment:

`export ROOTANASYS="/path/to/rootana/installation/"`

`export OPENCVSYS="/path/to/opencv/installation/"`

In the file you use to call the source of thisroot.sh (your setup file or .bashrc), add after the source of the thisroot.sh the line:

`export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$OPENCVSYS`

## Installation

`git clone https://github.com/CYGNUS-RD/digitizationpp.git`

`cd digitizationpp`

`export CXX="/path/to/your/c++17or20/compiler"`

`mkdir build-dir && cd build-dir`

`cmake ..`

`cmake --build .`



Generate documentation inside the `doc/html` folder:

`doxygen doc/doxygen.cfg`

Documentation should be then available at `doc/html/index.html`.

## Suggested usage

Put all the MC `.root` files you want to digitize in an input folder (<input_folder>)
Move to a desired folder where to launch the code

`./<path_to_build-dir>/digitizationpp <path_to_digitizationpp-dir>/config/ConfigFile_new.txt -I <path_to_input_folder> -O <path_to_output_folder>`

If not existing, the Outdir will be created. The -I and -O options can be droppped and the code will search inputfiles in
`<path_to_digitizationpp-dir>/input/`
and the outdir will be created in 
`<path_to_digitizationpp-dir>/OutDir/`

# Example
From `<path_to_digitizationpp-dir>`

`cd build-dir`

`./digitizationpp ../config/ConfigFile_new.txt`
