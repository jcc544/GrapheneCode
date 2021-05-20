program Graphene2
  implicit none
  external :: DGEEV, PRINT_EIGENVALUES
  integer, parameter :: dp = selected_real_kind(15,300)
  integer :: nx, ny, n, i, c, j, temp, info, outunit = 40, istat
  
  integer, dimension(:,:), allocatable :: unitcell, nn
  
  real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(kind=dp), dimension(:), allocatable :: Rwork
  real(kind=dp) :: k, t = 1.0_dp
  
  complex*16 :: im = (0,1), te
  complex*16, dimension(:,:), allocatable :: H, VL, VR
  complex*16, dimension(:), allocatable :: W
  complex*16, dimension(66) :: work
  
  
  character(len=100) :: dummy, errormsg
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section reads the compiler flags for nx and ny values

  nx = 1
  ny = 1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section allocates the values to be what they need to be at the start of the code
  
  n = nx * ny                          !The total number of atoms we are dealing with
  allocate(unitcell(ny,nx))            !Unitcell gives a visual representation of the unit cell we are dealing with
  allocate(H(n,n))                     !The Hamiltonian matrix will be the one we need to solve at the end
  allocate(nn(n,4))                    !Each atom can have 4 nearest neighbours
  allocate(VL(n,n))
  allocate(VR(n,n))
  allocate(W(n))     
  allocate(Rwork(2**n))
  c = 1
  do i = 1, ny
    do j = 1, nx
      unitcell(i,j) = c                !Give each element in unit cell a distinct number to help classify each atom
      c = c + 1
    end do
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section will focus on correctly filling the nearest neighbour array
  !The array will be filled as follows : nearestneighbours = (upneighbourvalue, downneighbourvalue, leftneighbourvalue, rightneighbourvalue)
  c = 1
  do i = 1, ny
    do j = 1, nx

      if (i == 1) then
        nn(c,1) = (-1) * unitcell(i,j)                      !If we are on the first row there can be no neighbour above
      else
        nn(c,1) = unitcell(i-1,j)        !If we are not on the first row, the nearest neighbour is simply directly above
      end if
      
      if (i == ny) then                !If we are on the last row, there can be no neighbour below
        nn(c,2) = (-1) * unitcell(i,j)
      else
        nn(c,2) = unitcell(i+1,j)
      end if
      
      if (j == 1) then
        nn(c,3) = (-1) * unitcell(i,nx)  !If we are in the first column there can be a nearest neighbour to the left, but we must make it negative because it is outside the cell
      else
        nn(c,3) = unitcell(i, j-1)
      end if
      
      if (j == nx) then
        nn(c,4) = (-1) * unitcell(i,1)    !If we are in the last column there can be a nearest neighbour to the right, but again it is outside the cell so must make it negative.
      else
        nn(c,4) = unitcell(i,j+1)
      end if
      c = c + 1
    end do   
  end do
  print *, nn
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !The k loop starts here
  open(unit = outunit, file = "Band.dat", status = "replace", action = "write", iostat = istat, iomsg = errormsg)
  k = (-1 * pi)
  do while (k <= pi)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section is for filling the relevant hamiltonian matrix
  
  H = 0                                         !Make sure H is filled with 0s to start with
  do c = 1, n                                    !c goes from 1 to 6, one line for each atom in the unit cell
    do i = 1, 4                                    !i goes from 1 to 4, for each of the possible nearest neighbours
      if (nn(c,i) == 0) then                       !if nn(c,i) is 0 there was no valid nearest neighbour for this column, so we skip past it
        cycle                          
      else if (nn(c,i) > 0) then                   !If a nearest neighbour was present, we put the hopping energy term in the relevant column in the hamiltonian
        H(c,nn(c,i)) = H(c,nn(c,i)) - t            
      else if (nn(c,i) < 0) then                   !If the nearest neighbour was negative, it was in a different cell, so we have to include the exponential term
        temp = int(abs(nn(c,i)))
        te = t * exp(im * k * (-1)**c)
        H(c,temp) = H(c,temp) - te
      end if
    end do
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section uses LAPACK to solve the eigenproblem
  call zgeev('N', 'N', n, H, n, W, VL, n, VR, n, work, 66, Rwork, info)
  if (info /= 0) then
    print *, "Something has gone wrong"
    stop
  end if
  do i = 1, n
    print *, k, W(i)
  end do
  write(outunit,*) k, real(W)
  k = k + 0.01
  
  end do
  close(outunit)
  contains
  function getoption(flag,getval,cvalue)       !Dummy needs to be character
    implicit none
    character(*) :: flag, cvalue
    character(160) :: arg
    logical :: getoption, getval
    integer :: i,j
    
    getoption = .false.
    i = 0
    do j = 1, iargc()
      call getarg(j,arg)
      if (arg .eq. flag) i = j
    end do
    if (i .gt. 0) then
      getoption = .true.
    end if
    if (getval) call getarg(i+1,cvalue)
  end function getoption

end program Graphene2
