!     
! File:   datasetlist.f90
! Author: alansoper
!
! Created on April 25, 2012, 11:11 AM
!

MODULE datasetlist
    
    use HDF5
    use inputfilestrings
    
    implicit none
    
    INTEGER(HID_T) :: pl_id         ! Property list identifier
    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: dspace_id     ! Data space identifier
    INTEGER(HID_T) :: dtype_id      ! Data type identifier
    INTEGER(HID_T) :: mspace_id     ! Memory space identifier
    INTEGER :: class                ! Datatype class, possible values are:
                                         !    H5T_NO_CLASS_F
                                         !    H5T_INTEGER_F
                                         !    H5T_FLOAT_F
                                         !    H5T_STRING_F
                                         !    H5T_BITFIELD_F
                                         !    H5T_OPAQUE_F
                                         !    H5T_COMPOUND_F
                                         !    H5T_REFERENCE_F
                                         !    H5T_ENUM_F
    INTEGER(SIZE_T) :: precision    ! Datatype precision
    INTEGER :: ndims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims    ! Array to store dimension sizes
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: maxdims ! Array to store max dimension sizes
    INTEGER(HSIZE_T)                            :: buffersize     ! Conversion buffer size
    INTEGER     ::   error ! Error flag

    integer, parameter :: ncharline=256
    integer :: nmembers=0
    character(len=ncharline), dimension(:), allocatable :: dsetline,dsetlinesave
    integer, dimension(:), allocatable :: linetype,linetypesave
    CHARACTER(LEN=ncharline) :: rootname
    
    CONTAINS
    
    subroutine initialise_list(nlinesnew)

        integer, intent(in) :: nlinesnew
        ! Ensure list has no members
        call reinitialise_list()
        allocate(dsetline(nlinesnew),linetype(nlinesnew))
        nmembers=nlinesnew
    
    end subroutine initialise_list
    
    subroutine reinitialise_list()
    
        !sets the number of members to zero and deallocates the arrays
        nmembers=0
        ndims=0
        if(allocated(dsetline)) deallocate(dsetline)
        if(allocated(linetype)) deallocate(linetype)
        if(allocated(dims)) deallocate(dims)
        if(allocated(maxdims)) deallocate(maxdims)

    end subroutine reinitialise_list
    
    subroutine retrieve_list(filename)

        !Retrieves the list of datasets in a specified hdf5 file
        character(len=*), intent(in) :: filename
!        integer, intent(out) :: ioerror
        CHARACTER(LEN=ncharline) :: dsetname,dsetnameused,dsetnameusedsave
        CHARACTER(LEN=ncharline) :: root_name,group_name !To hold object's name
        INTEGER :: nmembersroot,nmembersextra ! Number of members in group

        INTEGER     ::  i,j
        
        ! Set list to empty
        call reinitialise_list()
        
        ! Initialize FORTRAN interface.

        call h5open_f(error)

        ! Open an existing file.

        call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

        if(error.ne.0) then 
            return
        end if


        ! Get number of members in root group
        call h5gn_members_f(file_id,"/", nmembersroot, error)
        if(error.ne.0) stop

        if(nmembersroot.gt.0) then


            !Search for all the datasets in the file and store the corresponding path

            call h5gget_obj_info_idx_f(file_id,"/",0,root_name, class, error)
            rootname="/"//root_name(1:len_trim(root_name))
            call h5gn_members_f(file_id,rootname, nmembersroot, error)

            call initialise_list(nmembersroot)
            write(6,'(a,1x,a,1x,i5)') 'Root group and no. of members: ', rootname(1:len_trim(rootname)),nmembers
            dsetnameused=""

            i=0
            do while(i.lt.nmembers)
                call h5gget_obj_info_idx_f(file_id,rootname(1:len_trim(rootname))//dsetnameused,i, group_name, class, error) !Determine if this is a group or dataset
                i=i+1
                linetype(i)=class
                dsetline(i)=dsetnameused(1:len_trim(dsetnameused))//"/"//group_name(1:len_trim(group_name))
!            write(6,'(a,1x,i4,1x,i4,1x,a)') 'Group number, type and name = ',i,types(i),name_buffer(i)(1:len_trim(name_buffer(i)))
            end do

            !OK now step through each of these members and check if any are actually groups, in which case
            !we need to stop and get their names

            i=1
            do while(i.le.nmembers)
                if(linetype(i).eq.0) then
                    !Get the number of members in this group
                    dsetnameused=dsetline(i)
                    call h5gn_members_f(file_id,rootname(1:len_trim(rootname))//dsetnameused, nmembersextra, error)
                    !reallocate dsetlist with the new size
                    call extend_list(nmembersextra-1,i)
                    !Insert the new names
                    j=0
                    do while(j < nmembersextra)
                        call h5gget_obj_info_idx_f(file_id,rootname(1:len_trim(rootname))&
                        //dsetnameused,j, group_name, class, error) !Determine if this is a group or dataset
                        linetype(i+j)=class
                        dsetline(i+j)=dsetnameused(1:len_trim(dsetnameused))//"/"//group_name(1:len_trim(group_name))
                        j=j+1
                    end do
                else
                    i=i+1
                end if
            end do
            !Output results
            do i=1,nmembers
                dsetname=dsetline(i)
                !create name of group to be used
                dsetnameused="/"//root_name(1:len_trim(root_name))//dsetname

                !Open dataset

                call h5dopen_f(file_id, dsetnameused, dset_id, error)

                ! Get the associated data space

                call h5dget_space_f(dset_id,dspace_id,error)

                ! Get the datatype

                call h5dget_type_f(dset_id,dtype_id,error)

                ! Get the class of data stored and its precision

                call h5tget_class_f(dtype_id, class, error)

                if(class.ne.H5T_STRING_F) then
                    call h5tget_precision_f(dtype_id, precision, error)
                else
                    call h5tget_size_f(dtype_id, precision, error)
                endif

                call h5sget_simple_extent_ndims_f(dspace_id,ndims,error)

! Allocate the dimension arrays

                allocate(dims(ndims),maxdims(ndims))

! Get the dimensions

                call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)
                write(dsetline(i),'(a,10(1x,i4))') dsetline(i)(1:len_trim(dsetline(i))),class,precision&
                ,ndims,(dims(j),j=1,ndims)

                deallocate(dims,maxdims)
                error=0

                !Close the data type

                call h5tclose_f(dtype_id,error)

                ! Close the dataspace for the dataset.

                CALL h5sclose_f(dspace_id, error)

                ! Close the dataset.

                CALL h5dclose_f(dset_id, error)
              
            end do
        
        end if

        call h5fclose_f(file_id, error)
        call h5close_f(error)
        
    end subroutine retrieve_list
    
    subroutine extend_list(nlinesadd, iin)
        !extends the size of datasetline by nlinesnew.
        !if istart is present then it inserts the new lines starting from istart
        integer, intent(in) :: nlinesadd
        integer, intent(in), optional :: iin
        integer :: i,istart
        allocate(dsetlinesave(nmembers),linetypesave(nmembers))
        if (present(iin)) then
            istart=iin
        else
            istart=nmembers
        end if
        if (istart.lt.1.or.istart.gt.nmembers) istart=nmembers
        do i=1,nmembers
            dsetlinesave(i)=dsetline(i)
            linetypesave(i)=linetype(i)
        end do
        deallocate(dsetline,linetype)
        allocate(dsetline(nmembers+nlinesadd),linetype(nmembers+nlinesadd))
        i=0
        do while(i<istart)
            i=i+1
            dsetline(i)=dsetlinesave(i)
            linetype(i)=linetypesave(i)
        end do
        do while (i<nmembers)
            i=i+1
            dsetline(i+nlinesadd)=dsetlinesave(i)
            linetype(i+nlinesadd)=linetypesave(i)
        end do
        nmembers=nmembers+nlinesadd
        deallocate(dsetlinesave,linetypesave)
    end subroutine extend_list
    
    function retrieve_dims(filename,dsetname)
        !Retrieves the dimensions and precision of the specified dataset.
        !Returns nothing if dataset not found.
        integer                                 :: retrieve_dims
        character(len=ncharline), intent(in)    :: filename
        character(len=ncharline), intent(in)    :: dsetname
        
        integer                                 :: i,j,readj
        logical                                 :: found
        
        ! Deallocate the dims arrays if allocated
        ndims=0
        if (allocated(dims)) deallocate(dims)
        if (allocated(maxdims)) deallocate(maxdims)
        !First check that the list of datasets has been set up
        if(nmembers.eq.0) call retrieve_list(filename)
        !Step through the list and find the selected dataset
        i=0
        found=.false.
        do while (i.lt.nmembers.and..not.found)
            i=i+1
            call parse(dsetline(i))
            found=dsetname(1:len_trim(dsetname)).eq.dsetline(i)(ncf(1):ncl(1))
        end do
        if(found) then
            class=-1
            if(nwords.gt.1) read(dsetline(i)(ncf(2):ncl(2)),*) class
            precision=0
            if(nwords.gt.2) read(dsetline(i)(ncf(3):ncl(3)),*) precision
            ndims=0
            if(nwords.gt.3) read(dsetline(i)(ncf(4):ncl(4)),*) ndims
!            write(6,*) nwords,ndims
            if(ndims.gt.0) allocate(dims(ndims))
            j=0
            do while (j<ndims)
                j=j+1
                readj=j+4
                read(dsetline(i)(ncf(readj):ncl(readj)),*) dims(j)
            end do
            !retrieve the data
!            write(6,'(a,a,1x,i4)') 'Found: ',dsetname(1:len_trim(dsetname)),class
            retrieve_dims=i
        else
            retrieve_dims=0
        end if

    end function retrieve_dims
        
    function retrieve_data_c(filename,dsetname)
        !retrieves a character dataset
        character(len=ncharline) :: retrieve_data_c
        character(len=ncharline), intent(in) :: filename
        character(len=ncharline), intent(in) :: dsetname
        character(len=ncharline) :: dsetnameused
        character(len=ncharline) :: cdataset0d
        
        integer :: i,j,readj
        logical :: found
        
        cdataset0d=' '
        ! Deallocate the dims arrays if allocated
        ndims=0
        if (allocated(dims)) deallocate(dims)
        if (allocated(maxdims)) deallocate(maxdims)
        !First check that the list of datasets has been set up
        if(nmembers.eq.0) call retrieve_list(filename)
        !Step through the list and find the selected dataset
        i=retrieve_dims(filename,dsetname)
        found=i.gt.0
        if(found.and.class.eq.H5T_STRING_F) then
            !retrieve the data
            dsetnameused=rootname(1:len_trim(rootname))//dsetname
            !Open the dataset and read the data
!      Initialize FORTRAN interface.
!
            call h5open_f(error)
!
!     Open the existing file.
!
            call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

            ! Open dataset

            call h5dopen_f(file_id, dsetnameused, dset_id, error)

            ! Get the associated data space
            
            call h5dget_space_f(dset_id,dspace_id,error)

            ! Get the datatype

            call h5dget_type_f(dset_id,dtype_id,error)
            
            if(ndims.eq.1.and.dims(1).eq.1) then

                call h5dread_f(dset_id, dtype_id, cdataset0d, dims, error)

            endif
            
            call h5tclose_f(dtype_id,error)
            call h5sclose_f(dspace_id, error)
            call h5dclose_f(dset_id, error)
            call h5fclose_f(file_id, error)
            call h5close_f(error)
            
        end if

        retrieve_data_c=cdataset0d
    
    end function retrieve_data_c
        
    subroutine retrieve_data_r(filename,dsetname,rdataset1d,rrdataset1d)
        !retrieves a real 1D dataset
        character(len=ncharline), intent(in) :: filename
        character(len=ncharline), intent(in) :: dsetname
        real(kind=4), intent(inout), dimension (:), allocatable :: rdataset1d
        real(kind=8), intent(inout), dimension (:), allocatable :: rrdataset1d
        character(len=ncharline) :: dsetnameused
        
        integer :: i,j,readj
        logical :: found
        
        ! Deallocate the dims arrays if allocated
        ndims=0
        if (allocated(dims)) deallocate(dims)
        if (allocated(maxdims)) deallocate(maxdims)
        !First check that the list of datasets has been set up
        if(nmembers.eq.0) call retrieve_list(filename)
        !Step through the list and find the selected dataset
        i=retrieve_dims(filename,dsetname)
        found=i.gt.0
        if(found.and.class.eq.H5T_FLOAT_F) then
            !retrieve the data
            dsetnameused=rootname(1:len_trim(rootname))//dsetname
!            write(6,*) filename(1:len_trim(filename))
!            write(6,*) dsetnameused(1:len_trim(dsetnameused))
            !Open the dataset and read the data
!      Initialize FORTRAN interface.
!
            call h5open_f(error)
!            write(6,*) 'h5open_f error ',error
!
!     Open the existing file.
!
            call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
!           write(6,*) 'h5fopen_f error ',error

            ! Open dataset

            call h5dopen_f(file_id, dsetnameused, dset_id, error)
!            write(6,*) 'h5dopen_f error ',error

            ! Get the associated data space
            
            call h5dget_space_f(dset_id,dspace_id,error)
!            write(6,*) 'h5dget_space_f error ',error

            ! Get the datatype

            call h5dget_type_f(dset_id,dtype_id,error)
!            write(6,*) 'h5dget_type_f error ',error

            if(precision.eq.32.and.ndims.eq.1) then

!                write(6,*) dims(1)
                if(allocated(rdataset1d)) deallocate(rdataset1d)
                allocate(rdataset1d(dims(1)))
                call h5dread_f(dset_id, dtype_id, rdataset1d, dims, error)
!                write(6,*) 'h5dread_f error ',error

            else if(precision.eq.64.and.ndims.eq.1) then

!                write(6,*) dims(1)
                if(allocated(rrdataset1d)) deallocate(rrdataset1d)
                allocate(rrdataset1d(dims(1)))
                call h5dread_f(dset_id, dtype_id, rrdataset1d, dims, error)
!                write(6,*) 'h5dread_f error ',error

            endif
            
            call h5tclose_f(dtype_id,error)
            call h5sclose_f(dspace_id, error)
            call h5dclose_f(dset_id, error)
            call h5fclose_f(file_id, error)
            call h5close_f(error)
            
        end if

    end subroutine retrieve_data_r
        
    subroutine retrieve_data_i(filename,dsetname&
        ,idataset0d,idataset1d,idataset2d,idataset3d)
        !retrieves a character dataset
        character(len=ncharline), intent(in) :: filename
        character(len=ncharline), intent(in) :: dsetname
        INTEGER, intent(out), TARGET                                :: idataset0d
        INTEGER, intent(out), DIMENSION(:), ALLOCATABLE, TARGET     :: idataset1d
        INTEGER, intent(out), DIMENSION(:,:), ALLOCATABLE, TARGET   :: idataset2d
        INTEGER, intent(out), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: idataset3d
        
        character(len=ncharline) :: dsetnameused
        
        integer :: i,j,readj
        logical :: found
        
        ! Deallocate the dims arrays if allocated
        ndims=0
        if (allocated(dims)) deallocate(dims)
        if (allocated(maxdims)) deallocate(maxdims)
        !First check that the list of datasets has been set up
        if(nmembers.eq.0) call retrieve_list(filename)
        !Step through the list and find the selected dataset
        i=retrieve_dims(filename,dsetname)
        found=i.gt.0
        if(found.and.class.eq.H5T_INTEGER_F) then
            !retrieve the data
            dsetnameused=rootname(1:len_trim(rootname))//dsetname
            ! Open the dataset and read the data
           ! Initialize FORTRAN interface.

            call h5open_f(error)

            ! Open the existing file.

            call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

            ! Open dataset

            call h5dopen_f(file_id, dsetnameused, dset_id, error)

            ! Get the associated data space
            
            call h5dget_space_f(dset_id,dspace_id,error)
            
            ! Get the datatype

            call h5dget_type_f(dset_id,dtype_id,error)
            
            if(precision.eq.32) then

                if(ndims.eq.1) then
                    if(dims(1).eq.1) then
                        call h5dread_f(dset_id, dtype_id, idataset0d, dims, error)                        
                    else if(dims(1).gt.1) then
                        if(allocated(idataset1d)) deallocate(idataset1d)
                        allocate(idataset1d(dims(1)))
                        call h5dread_f(dset_id, dtype_id, idataset1d, dims, error)
                    end if
                else if(ndims.eq.2) then
                    if(allocated(idataset2d)) deallocate(idataset2d)
                    allocate(idataset2d(dims(1),dims(2)))
                    call h5dread_f(dset_id, dtype_id, idataset2d, dims, error)
                else if(ndims.eq.3) then
                    if(allocated(idataset3d)) deallocate(idataset3d)
                    allocate(idataset3d(dims(1),dims(2),dims(3)))
                    call h5dread_f(dset_id, dtype_id, idataset3d, dims, error)
                end if
            endif
            
            
            call h5tclose_f(dtype_id,error)
            call h5sclose_f(dspace_id, error)
            call h5dclose_f(dset_id, error)
            call h5fclose_f(file_id, error)
            call h5close_f(error)
            
        end if

    end subroutine retrieve_data_i
        
    subroutine retrieve_dataslab_i(filename,dsetname,&
        firstcolumn,maxcolumns,dimen3,idataset3d,inputdims,inputoffsets)
        !retrieves a slab of 2d integer data from a 3d dataset.
        !The column offset is given by offset.
        !On entry maxcolumns gives the maximum number of columns to be output.
        !If zero on entry it will be set to the column dimension of the dataset array.
        !On exit maxcolumns gives the column dimension of the output array, which may be 
        !less than input value if the slab would extend off the range of the dataset
        !dimen3 gives the third dimension of the dataset to be output.
        !If this value is specified less than 1, it is set to 1. If greater than dims(3) it is set to dims(3)
        !Output in idataset2d
        character(len=ncharline), intent(in)        :: filename
        character(len=ncharline), intent(in)        :: dsetname
        integer, intent(inout)                      :: dimen3
        integer, intent(inout)                      :: firstcolumn,maxcolumns
!        integer, intent(inout), dimension(:,:)      :: idataset2d
        integer, intent(inout), dimension(:,:,:)                   :: idataset3d
        integer, intent(in), dimension(:)           :: inputdims,inputoffsets

        integer(HSIZE_T), dimension(3) :: dsetoffset,dsetslabsize
        integer :: mspacerank = 3
        integer(HSIZE_T), dimension(3) :: mspacedims,mspaceoffset,mspacesize
        
        character(len=ncharline) :: dsetnameused
        
        integer :: i,j,readj,offset,ms1,ms2,iref,jref
        logical :: found
        
        ! Deallocate the dims arrays if allocated
        ndims=0
        if (allocated(dims)) deallocate(dims)
        if (allocated(maxdims)) deallocate(maxdims)
        !First check that the list of datasets has been set up
        if(nmembers.eq.0) call retrieve_list(filename)
        !Step through the list and find the selected dataset
        i=retrieve_dims(filename,dsetname)
        found=i.gt.0
        if(found.and.class.eq.H5T_INTEGER_F) then
            if (ndims.eq.3.and.firstcolumn.le.dims(2)) then
                !Set up the dimensions of the slab to be extracted
                if(firstcolumn.lt.1) firstcolumn=1
                if(firstcolumn.gt.dims(2)) return !Cannot get non-existent data!
                offset=firstcolumn-1
                if(maxcolumns.le.0) maxcolumns=dims(2) !Get the whole lot if maxcolumns is zero or less
                maxcolumns=min(maxcolumns,dims(2)-offset)
                if(dimen3.lt.1) dimen3=1
                if(dimen3.gt.dims(3)) dimen3=dims(3)
                dsetoffset(1)=0
                dsetoffset(2)=offset
                dsetoffset(3)=dimen3-1
                dsetslabsize(1)=dims(1)
                dsetslabsize(2)=maxcolumns
                dsetslabsize(3)=1
                mspacerank=3
                mspacedims(1)=inputdims(1) !Responsibility for checking these are correct lies with the user
                mspacedims(2)=inputdims(2)
                mspacedims(3)=inputdims(3)
                mspaceoffset(1)=inputoffsets(1)
                mspaceoffset(2)=inputoffsets(2)
                mspaceoffset(3)=inputoffsets(3)
                mspacesize(1)=dims(1)
                mspacesize(2)=maxcolumns
                mspacesize(3)=1
                ! Make the buffer size appropriate to the size of data being recovered
                buffersize=dsetslabsize(1)*dsetslabsize(2)*dsetslabsize(3)*4*2
                !retrieve the data
                dsetnameused=rootname(1:len_trim(rootname))//dsetname
!                call system_clock(ms1)
!                write(6,*) ms1
                ! Open the dataset and read the data
                ! Initialize FORTRAN interface.

                call h5open_f(error)
                
                ! Open the existing file.

                call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

                ! Open dataset

                call h5dopen_f(file_id, dsetnameused, dset_id, error)

                ! Get the associated data space
            
                call h5dget_space_f(dset_id,dspace_id,error)

                ! Get the datatype

                call h5dget_type_f(dset_id,dtype_id,error)
            
                ! Select hyperslab in the dataset.

                call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
                dsetoffset, dsetslabsize, error)
            
                ! Create memory dataspace.

                CALL h5screate_simple_f(mspacerank, mspacedims, mspace_id, error)

                ! Select hyperslab in memory.

                CALL h5sselect_hyperslab_f(mspace_id, H5S_SELECT_SET_F, &
                mspaceoffset, mspacesize, error)

                !Create the data transfer property list identifier
                
                call h5pcreate_f(H5P_DATASET_XFER_F, pl_id, error)
                call h5pset_buffer_f(pl_id, buffersize, error)
                
                if(precision.eq.32) then

                    ! Read data from hyperslab in the file into the hyperslab in
                    ! memory and display.
                    
                    CALL H5dread_f(dset_id, dtype_id, idataset3d, mspacedims, error, &
                    mspace_id, dspace_id, pl_id)

                endif
                call h5pclose_f(pl_id, error)
                call h5sclose_f(mspace_id, error)
                call h5tclose_f(dtype_id,error)
                call h5sclose_f(dspace_id, error)
                call h5dclose_f(dset_id, error)
                call h5fclose_f(file_id, error)
                call h5close_f(error)
            end if
        end if

    end subroutine retrieve_dataslab_i
        
END MODULE datasetlist
