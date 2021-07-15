!     
! File:   reallocation_routines.f90
! Author: aks45
!
! Created on 01 May 2012, 16:33
!

MODULE reallocation_routines

!Use for reallocating memory. 
    
CONTAINS
    
    subroutine reallocate1d_i(array,newsize)
        !reallocates an array, new values are left uninitialised
        integer, dimension(:), intent(inout), allocatable :: array
        integer             :: oldsize
        integer, intent(in) :: newsize
        integer, dimension(:), allocatable :: arraysave
        integer :: i
        
        if(allocated(array)) then
            oldsize=size(array)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size! 
        if(oldsize.ge.newsize) return
        if(oldsize.gt.0) then
            allocate(arraysave(oldsize))
            i=0
            do while (i.lt.oldsize)
                i=i+1
                arraysave(i)=array(i)
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize.gt.0) then
            allocate(array(newsize))
            if(oldsize.gt.0) then
                i=0
                do while (i.lt.min(oldsize,newsize))
                    i=i+1
                    array(i)=arraysave(i)
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate1d_i
    
    subroutine reallocate2d_i(array,newsize1,newsize2)
        !reallocates an array, new values are left uninitialised
        integer, dimension(:,:), intent(inout), allocatable :: array
        integer, intent(in)                                 :: newsize1,newsize2
        integer, dimension(2)                               :: oldsize
        integer, dimension(2)                               :: newsize
        integer, dimension(:,:), allocatable                :: arraysave
        integer :: i,j
        
        newsize(1)=newsize1
        newsize(2)=newsize2
        if(allocated(array)) then
            oldsize(1)=size(array,1)
            oldsize(2)=size(array,2)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize(1).ge.newsize(1).and.oldsize(2).ge.newsize(2)) return
        if(oldsize(1).gt.0.and.oldsize(2).gt.0) then
            allocate(arraysave(oldsize(1),oldsize(2)))
            i=0
            do while (i.lt.oldsize(2))
                i=i+1
                j=0
                do while (j.lt.oldsize(1))
                    j=j+1
                    arraysave(j,i)=array(j,i)
                end do
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize(1).gt.0.and.newsize(2).gt.0) then
            allocate(array(newsize(1),newsize(2)))
            if(oldsize(1).gt.0.and.oldsize(2).gt.0) then
                i=0
                do while(i.lt.min(oldsize(2),newsize(2)))
                    i=i+1
                    j=0
                    do while (j.lt.min(oldsize(1),newsize(1)))
                        j=j+1
                        array(j,i)=arraysave(j,i)
                    end do
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate2d_i
    
    subroutine reallocate3d_i(array,newsize1,newsize2,newsize3)
        !reallocates an array, new values are left uninitialised
        integer, dimension(:,:,:), intent(inout), allocatable :: array
        integer, intent(in)                     :: newsize1,newsize2,newsize3
        integer, dimension(3)                   :: oldsize
        integer, dimension(3)                   :: newsize
        integer, dimension(:,:,:), allocatable  :: arraysave
        integer :: i,j,k
        
        newsize(1)=newsize1
        newsize(2)=newsize2
        newsize(3)=newsize3
        if(allocated(array)) then
            oldsize(1)=size(array,1)
            oldsize(2)=size(array,2)
            oldsize(3)=size(array,3)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize(1).ge.newsize(1).and.oldsize(2).ge.newsize(2).and.oldsize(3).ge.newsize(3)) return
        if(oldsize(1).gt.0.and.oldsize(2).gt.0.and.oldsize(3).gt.0) then
            allocate(arraysave(oldsize(1),oldsize(2),oldsize(3)))
            i=0
            do while (i.lt.oldsize(3))
                i=i+1
                j=0
                do while (j.lt.oldsize(2))
                    j=j+1
                    k=0
                    do while (k.lt.oldsize(1))
                        k=k+1
                        arraysave(k,j,i)=array(k,j,i)
                    end do
                end do
            end do
        end if
        if(allocated(array)) deallocate(array)
        allocate(array(newsize(1),newsize(2),newsize(3)))
        if(oldsize(1).gt.0.and.oldsize(2).gt.0.and.oldsize(3).gt.0) then
            i=0
            do while(i.lt.min(oldsize(3),newsize(3)))
                i=i+1
                j=0
                do while (j.lt.min(oldsize(2),newsize(2)))
                    j=j+1
                    k=0
                    do while (k.lt.min(oldsize(1),newsize(1)))
                        k=k+1
                        array(k,j,i)=arraysave(k,j,i)
                    end do
                end do
            end do
        end if
        if(allocated(arraysave)) deallocate(arraysave)
!        write(6,*) 'reallocate3d_i> ',size(array,1),size(array,2),size(array,3)
        return
    end subroutine reallocate3d_i
    
    subroutine reallocate1d_r(array,newsize)
        !reallocates an array, new values are left uninitialised
        real, dimension(:), intent(inout), allocatable :: array
        integer :: oldsize
        integer, intent(in) :: newsize
        real, dimension(:), allocatable :: arraysave
        integer :: i
        
        if(allocated(array)) then
            oldsize=size(array)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize.ge.newsize) return
        if(oldsize.gt.0) then
            allocate(arraysave(oldsize))
            i=0
            do while (i.lt.oldsize)
                i=i+1
                arraysave(i)=array(i)
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize.gt.0) then
            allocate(array(newsize))
            if(oldsize.gt.0) then
                i=0
                do while(i.lt.min(oldsize,newsize))
                    i=i+1
                    array(i)=arraysave(i)
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate1d_r
    
    subroutine reallocate2d_r(array,newsize1,newsize2)
        !reallocates an array, new values are left uninitialised
        real, dimension(:,:), intent(inout), allocatable :: array
        integer, intent(in)                 :: newsize1,newsize2
        integer, dimension(2)               :: oldsize
        integer, dimension(2)               :: newsize
        real, dimension(:,:), allocatable   :: arraysave
        integer :: i,j
        
        newsize(1)=newsize1
        newsize(2)=newsize2
        if(allocated(array)) then
            oldsize(1)=size(array,1)
            oldsize(2)=size(array,2)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize(1).ge.newsize(1).and.oldsize(2).ge.newsize(2)) return
        if(oldsize(1).gt.0.and.oldsize(2).gt.0) then
            allocate(arraysave(oldsize(1),oldsize(2)))
!            write(6,*) 'reallocate2d_r> 3 ',oldsize(1),oldsize(2),size(arraysave,1),size(arraysave,2)
            
            i=0
            do while (i.lt.oldsize(2))
                i=i+1
                j=0
                do while (j.lt.oldsize(1))
                    j=j+1
                    arraysave(j,i)=array(j,i)
                end do
            end do
        end if
        if(allocated(array)) then
!        write(6,*) 'reallocate2d_r> 5 ',newsize(1),newsize(2)
            deallocate(array)
!        write(6,*) 'reallocate2d_r> 6 ',newsize(1),newsize(2)
        end if
        if(newsize(1).gt.0.and.newsize(2).gt.0) then
            allocate(array(newsize(1),newsize(2)))
            if(oldsize(1).gt.0.and.oldsize(2).gt.0) then
                i=0
                do while(i.lt.min(oldsize(2),newsize(2)))
                    i=i+1
                    j=0
                    do while (j.lt.min(oldsize(1),newsize(1)))
                        j=j+1
                        array(j,i)=arraysave(j,i)
                    end do
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate2d_r
    
    subroutine reallocate3d_r(array,newsize1,newsize2,newsize3)
        !reallocates an array, new values are left uninitialised
        real, dimension(:,:,:), intent(inout), allocatable :: array
        integer, intent(in)                 :: newsize1,newsize2,newsize3
        integer, dimension(3)               :: oldsize
        integer, dimension(3)               :: newsize
        real, dimension(:,:,:), allocatable :: arraysave
        integer :: i,j,k
        
        newsize(1)=newsize1
        newsize(2)=newsize2
        newsize(3)=newsize3
        if(allocated(array)) then
            oldsize(1)=size(array,1)
            oldsize(2)=size(array,2)
            oldsize(3)=size(array,3)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize(1).ge.newsize(1).and.oldsize(2).ge.newsize(2).and.oldsize(3).ge.newsize(3)) return
        if(oldsize(1).gt.0.and.oldsize(2).gt.0.and.oldsize(3).gt.0) then
            allocate(arraysave(oldsize(1),oldsize(2),oldsize(3)))
            i=0
            do while (i.lt.oldsize(3))
                i=i+1
                j=0
                do while (j.lt.oldsize(2))
                    j=j+1
                    k=0
                    do while (k.lt.oldsize(1))
                        k=k+1
                        arraysave(k,j,i)=array(k,j,i)
                    end do                       
                end do
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize(1).gt.0.and.newsize(2).gt.0.and.newsize(3).gt.0) then
            allocate(array(newsize(1),newsize(2),newsize(3)))
            if(oldsize(1).gt.0.and.oldsize(2).gt.0.and.oldsize(3).gt.0) then
                i=0
                do while(i.lt.min(oldsize(3),newsize(3)))
                    i=i+1
                    j=0
                    do while (j.lt.min(oldsize(2),newsize(2)))
                        j=j+1
                        k=0
                        do while (k.lt.min(oldsize(1),newsize(1)))
                            k=k+1
                            array(k,j,i)=arraysave(k,j,i)
                        end do                       
                    end do
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate3d_r
    
    subroutine reallocate4d_r(array,newsize1,newsize2,newsize3,newsize4)
        !reallocates an array, new values are left uninitialised
        real, dimension(:,:,:,:), intent(inout), allocatable :: array
        integer, intent(in)                 :: newsize1,newsize2,newsize3,newsize4
        integer, dimension(4)               :: oldsize
        integer, dimension(4)               :: newsize
        real, dimension(:,:,:,:), allocatable :: arraysave
        integer :: i,j,k,l
        
        newsize(1)=newsize1
        newsize(2)=newsize2
        newsize(3)=newsize3
        newsize(4)=newsize4
        if(allocated(array)) then
            oldsize(1)=size(array,1)
            oldsize(2)=size(array,2)
            oldsize(3)=size(array,3)
            oldsize(4)=size(array,4)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize(1).ge.newsize(1).and.oldsize(2).ge.newsize(2).and.oldsize(3).ge.newsize(3).and.oldsize(4).ge.newsize(4)) return
        if(oldsize(1).gt.0.and.oldsize(2).gt.0.and.oldsize(3).gt.0.and.oldsize(4).gt.0) then
            allocate(arraysave(oldsize(1),oldsize(2),oldsize(3),oldsize(4)))
            i=0
            do while (i.lt.oldsize(4))
                i=i+1
                j=0
                do while (j.lt.oldsize(3))
                    j=j+1
                    k=0
                    do while (k.lt.oldsize(2))
                        k=k+1
                        l=0
                        do while (l.lt.oldsize(1))
                            l=l+1
                            arraysave(l,k,j,i)=array(l,k,j,i)
                        end do
                    end do                       
                end do
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize(1).gt.0.and.newsize(2).gt.0.and.newsize(3).gt.0.and.newsize(4).gt.0) then
            allocate(array(newsize(1),newsize(2),newsize(3),newsize(4)))
            if(oldsize(1).gt.0.and.oldsize(2).gt.0.and.oldsize(3).gt.0.and.oldsize(4).gt.0) then
                i=0
                do while(i.lt.min(oldsize(4),newsize(4)))
                    i=i+1
                    j=0
                    do while (j.lt.min(oldsize(3),newsize(3)))
                        j=j+1
                        k=0
                        do while (k.lt.min(oldsize(2),newsize(2)))
                            k=k+1
                            l=0
                            do while (l.lt.min(oldsize(1),newsize(1)))
                                l=l+1
                                array(l,k,j,i)=arraysave(l,k,j,i)
                            end do
                        end do                       
                    end do
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate4d_r
    
    subroutine reallocate1d_c(array,ncharsize,newsize)
        !reallocates an array, new values are left uninitialised
        character(len=ncharsize), dimension(:), intent(inout), allocatable :: array
        integer, intent(in)                                                 :: ncharsize
        integer                                                             :: oldsize
        integer, intent(in)                                                 :: newsize
        character(len=ncharsize), dimension(:), allocatable                 :: arraysave
        integer                                                             :: i
        
        if(allocated(array)) then
            oldsize=size(array)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize.ge.newsize) return
        if(oldsize.gt.0) then
            allocate(arraysave(oldsize))
            i=0
            do while (i.lt.oldsize)
                i=i+1
                arraysave(i)=array(i)
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize.gt.0) then
            allocate(array(newsize))
            if(oldsize.gt.0) then
                i=0
                do while(i.lt.min(oldsize,newsize))
                    i=i+1
                    array(i)=arraysave(i)
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate1d_c
    
    subroutine reallocate2d_c(array,ncharsize,newsize1,newsize2)
        !reallocates an array, new values are left uninitialised
        character(len=ncharsize), dimension(:,:), intent(inout), allocatable :: array
        integer, intent(in)                 :: ncharsize
        integer, intent(in)                 :: newsize1,newsize2
        integer, dimension(2)               :: oldsize
        integer, dimension(2)               :: newsize
        character(len=ncharsize), dimension(:,:), allocatable   :: arraysave
        integer :: i,j
        
        newsize(1)=newsize1
        newsize(2)=newsize2
        if(allocated(array)) then
            oldsize(1)=size(array,1)
            oldsize(2)=size(array,2)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size!
        if(oldsize(1).ge.newsize(1).and.oldsize(2).ge.newsize(2)) return
        if(oldsize(1).gt.0.and.oldsize(2).gt.0) then
            allocate(arraysave(oldsize(1),oldsize(2)))
            i=0
            do while (i.lt.oldsize(2))
                i=i+1
                j=0
                do while (j.lt.oldsize(1))
                    j=j+1
                    arraysave(j,i)=array(j,i)
                end do
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize(1).gt.0.and.newsize(2).gt.0) then
            allocate(array(newsize(1),newsize(2)))
            if(oldsize(1).gt.0.and.oldsize(2).gt.0) then
                i=0
                do while(i.lt.min(oldsize(2),newsize(2)))
                    i=i+1
                    j=0
                    do while (j.lt.min(oldsize(1),newsize(1)))
                        j=j+1
                        array(j,i)=arraysave(j,i)
                    end do
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate2d_c
    
END MODULE reallocation_routines
