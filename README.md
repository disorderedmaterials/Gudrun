## Gudrun

### Data reduction for neutron and x-ray scattering experiments

#### Original code &copy; A. Soper 2021

#### CMake'd / tidied by T. Youngs 2021

#### Note For Apple Users

After installing the OSX binaries, the following commands must be run to use
Gudrun:

```
cd /folder/containing/binaries/
chmod +x gudrun_dcs purge_det tophatsub calc_corrsx_in_out
xattr -rd com.apple.quarantine gudrun_dcs purge_det tophatsub calc_corrsx_in_out
```

#### Installation instructions for RHEL 7 / CentOS 7 distributions:

> #install hdf5 library and cmake3 for compilation \
> yum install -y hdf5-devel cmake3 \
> #optional: copy needed hdf5 library if not automatically found \
> cp /usr/lib64/gfortran/modules/hdf5.mod src/neutron/libGudrunN/ \
> #create compilation directory \
> mkdir Gudrun/build \
> cd Gudrun/build \
> #build it \
> cmake3 ../ -DBUILD_SZIP:bool=true -DCMAKE_INSTALL_PREFIX:path=`pwd`/.. \
> make install \
> #clean up \
> cd ../ \
> rm -rf build \
> #binaries are now location in the 'bin' directory
