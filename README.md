## Gudrun

### Data reduction for neutron and x-ray scattering experiments

#### Original code &copy; A. Soper 2021

#### CMake'd / tidied by T. Youngs 2021

Refer to the [Gudrun Manual](https://www.isis.stfc.ac.uk/OtherFiles/Disordered%20Materials/Gudrun-Manual-2017-10.pdf) for more details.

## Installation Instructions

Pre-packaged binaries are avaliable in the [GitHub Releases Page](https://github.com/disorderedmaterials/Gudrun/releases), however if you would like to build the binaries yourself you can follow these steps:

#### Platform:
- [Linux - Ubuntu](#linux---ubuntu)
- [Linux - RHEL 7 / CentOS 7](#linux---rhel-7-and-centos-7)
- [MacOS](#macos)

### Linux - Ubuntu

Install prerequisites:
```
sudo apt-get update -q
sudo apt-get install ninja-build gfortran
```

Download HDF5:
```
wget https://github.com/disorderedmaterials/HDF5/releases/download/hdf5-1_10_7/hdf5-1_10_7-linux.zip
unzip hdf5-1_10_7-linux.zip
```

Build:
```
HDF5_DIR="$(pwd)/hdf5-1_10_7"
mkdir build
cd build
cmake ../ -G Ninja -DLOCAL_STATIC_HDF5:bool=true -DHDF5_DIR:path=${HDF5_DIR} -DCMAKE_Fortran_FLAGS:string="-cpp"
ninja
ninja install
```

Binaries should now be avaliable in the `Gudrun/bin` directory

### Linux - RHEL 7 and CentOS 7

Install prerequisites:
```
yum install -y hdf5-devel cmake3
```
Optional: copy needed hdf5 library if not automatically found

```
cp /usr/lib64/gfortran/modules/hdf5.mod src/neutron/libGudrunN/
```

Create compilation directory
```
mkdir Gudrun/build
cd Gudrun/build
```

Build:
```
cmake3 ../ -DBUILD_SZIP:bool=true -DCMAKE_INSTALL_PREFIX:path=`pwd`/..
make install
```
Binaries should now be avaliable in the `Gudrun/bin` directory


### MacOS
Install prerequisites:

```
brew update-reset
brew install ninja hdf5 zlib
```

Build:
If you are building on Apple Silicon, add an additional flag to the `cmake` command: `-DBUILD_HDF=false`.

```
mkdir build
cd build
cmake ../ -G Ninja -DCMAKE_Fortran_COMPILER:string="gfortran-11" -DBUILD_SZIP:bool=true -DCMAKE_Fortran_FLAGS:string="-cpp"
ninja
ninja install
```

After installing the OSX binaries, the following commands must be run to use
Gudrun:

```
cd Gudrun/bin
chmod +x gudrun_dcs purge_det tophatsub calc_corrsx_in_out
xattr -rd com.apple.quarantine gudrun_dcs purge_det tophatsub calc_corrsx_in_out
```
