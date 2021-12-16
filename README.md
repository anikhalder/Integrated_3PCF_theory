# The integrated 3PCF of projected fields

This repository hosts code for computing theoretical predictions for the integrated 3-point correlation function (i3PCF) of projected fields. It also computes predictions for the 2-point correlation function (2PCF). 

## Running the code

After you have setup the external libraries (see EXTERNAL FILES section below!), do the following in order to run the i3PCF code:

1. According to your needs change and save the code files present inside the src (e.g. main.cpp etc.) and include (e.g. constants.h) directories.
2. Go to the build directory and execute the following commands from inside the build directory:

cd /path/to/Integrated_3PCF_theory/build/

make clean

cmake ..

make -j

./Integrated_3PCF_theory

## External Files

Instructions to include the following external libraries (required for our i3PCF code to work). Perform the following steps:

### CLASS (Boltzmannn solver)

NOTE: Currently the i3PCF project is written based on the CLASS v2.9.4. There were major changes in syntax and variable naming which CLASS introduced in its v3.0 onwards. We plan to update our i3PCF code to these later versions of CLASS in the future.

1. Somewhere on your machine (e.g. inside a software folder in your home directory) clone the v2.9 branch of class_public

git clone -b 2.9 https://github.com/lesgourg/class_public.git

cd /path/to/class_public/

2. Compile and build this class_public repository (please follow the instructions for building and compiling CLASS on the class_public github page)

3. Once CLASS has been built, go to the build directory inside class_public. Compile the *.o files present inside that build directory and put them together into a single static library e.g. libclass.a

cd build/

ar rcs libclass.a *.o 

4. Copy this libclass.a library into the following folder of our Integrated_3PCF_theory folder i.e.

cp /path/to/class_public/build/libclass.a /path/to/Integrated_3PCF_theory/external/class_files/class_lib/.

5. Also copy the header files present in the include directory of the class_public folder into the following folder of our Integrated_3PCF_theory folder i.e.

cp /path/to/class_public/include/*.h /path/to/Integrated_3PCF_theory/external/class_files/class_include/.


### cubature (for multi-dimensional integration)

1. Somewhere on your machine (e.g. inside a software folder in your home directory) clone the cubature repository

git clone https://github.com/stevengj/cubature.git

cd /path/to/cubature/

2. Compile the hcubature.c and pcubature.c files using:

gcc -fPIC -c -O3 hcubature.c
gcc -fPIC -c -O3 pcubature.c

3. Create a build directory inside the cubature folder

mkdir build

4. Move the hcubature.o and pcubature.o files into the build directory

mv hcubature.o build/.
mv pcubature.o build/.

5. Go to the build directory and compile the *.o files that you copied inside that build directory and put them together into a single static library e.g. libcubature.a

cd build/
ar rcs libcubature.a *.o 

6. Copy this libcubature.a library into the following folder of our Integrated_3PCF_theory folder i.e.

cp /path/to/cubature/build/libcubature.a /path/to/Integrated_3PCF_theory/external/cubature_files/cubature_lib/.

7. Also copy the cubature.h header file present in the cubature directory into the following folder of our Integrated_3PCF_theory folder i.e.

cp /path/to/cubature/cubature.h /path/to/Integrated_3PCF_theory/external/cubature_files/cubature_include/.