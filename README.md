# The integrated 3PCF of projected cosmic density fields

This repository hosts code for computing theoretical predictions for the integrated 3-point correlation function of projected cosmic density fields. It also computes predictions for the 2-point correlation function of projected fields. 

## Running the code

After you have setup the external libraries (see External dependencies section below!), do the following in order to run the code:

1. According to your needs edit and save the code (e.g. main.cpp etc.) in the include and src directories. If you add your own .hpp (inside include) and .cpp (inside src) files, do not forget to update the src/CMakeLists.txt file accordingly!

2. Go to the build directory and execute the following commands:

cd /path/to/Integrated_3PCF_theory/build/

make clean

cmake ..

make -j

./Integrated_3PCF_theory

## External dependencies

Firstly, you need to have the **gsl** (https://www.gnu.org/software/gsl/) and **cmake** (https://cmake.org/install/) packages installed on your machine. Please do so if they are not already installed.

Furthermore, you will also need the **CLASS** and **cubature** packages. Please follow the instructions below on how to compile these two packages locally on your machine and include them into the **Integrated_3PCF_theory** directory.

### CLASS (Boltzmann solver)

NOTE: Currently the **Integrated_3PCF_theory** code is written based on the **CLASS** v2.9.4. There were major changes in syntax and file naming which **CLASS** introduced in its v3.0 onwards. We plan to update our code to be compatible with these later **CLASS** versions in the future.

1. Somewhere on your machine (e.g. inside a software folder in your home directory) download the v2.9.4 of **CLASS** from

https://github.com/lesgourg/class_public/releases/tag/v2.9.4

unzip the folder

cd /path/to/class_public-2.9.4/

2. Compile and build this repository (please follow the instructions for building and compiling **CLASS** on the class_public github page). One usually simply has to do:

make class -j

Note that you do NOT need to install the **CLASSY** python package for this project.

3. Once compilation is complete without errors, go to the build directory inside class_public-2.9.4. Compile the *.o files present inside and put them together into a single static library e.g. libclass.a

cd build/

ar rcs libclass.a *.o 

4. Copy this libclass.a library into the following directory of **Integrated_3PCF_theory**:

cp /path/to/class_public-2.9.4/build/libclass.a /path/to/Integrated_3PCF_theory/external/class_files/class_lib/.

5. Also copy the header files from the include directory of the **CLASS** folder to the following folder of **Integrated_3PCF_theory**:

cp /path/to/class_public-2.9.4/include/*.h /path/to/Integrated_3PCF_theory/external/class_files/class_include/.


### cubature (for multi-dimensional integration)

1. Somewhere on your machine (e.g. inside a software folder in your home directory) clone the **cubature** repository

git clone https://github.com/stevengj/cubature.git

cd /path/to/cubature/

2. Compile the hcubature.c and pcubature.c files using:

gcc -fPIC -c -O3 hcubature.c

gcc -fPIC -c -O3 pcubature.c

3. Create a build directory inside the **cubature** folder

mkdir build

4. Move the hcubature.o and pcubature.o files into the build directory

mv hcubature.o build/.

mv pcubature.o build/.

5. Go to the build directory and compile the *.o files that you moved in the previous step and put them together into a single static library e.g. libcubature.a

cd build/

ar rcs libcubature.a *.o 

6. Copy this libcubature.a library into the following directory of **Integrated_3PCF_theory**:

cp /path/to/cubature/build/libcubature.a /path/to/Integrated_3PCF_theory/external/cubature_files/cubature_lib/.

7. Also copy the cubature.h header file from the **cubature** directory to the following folder of **Integrated_3PCF_theory**:

cp /path/to/cubature/cubature.h /path/to/Integrated_3PCF_theory/external/cubature_files/cubature_include/.