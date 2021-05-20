program Graphene2
  implicit none
  external :: DGEEV, PRINT_EIGENVALUES
  integer, parameter :: dp = selected_real_kind(15,300)
  integer :: nx, n, i, c, temp, info, outunit = 40, istat
  
  integer, dimension(:,:), allocatable :: nn
  integer, dimension(:), allocatable :: unitcell, parity
  
  real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(kind=dp), dimension(:), allocatable :: Rwork
  real(kind=dp) :: k, t = 2.7_dp
  
  complex*16 :: im = (0,1), te
  complex*16, dimension(:,:), allocatable :: H, VL, VR
  complex*16, dimension(:), allocatable :: W
  complex*16, dimension(66) :: work
  
  
  character(len=100) :: dummy, errormsg
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section reads the compiler flags for nx and ny values

!   if (getoption('-nx',.true.,dummy)) then   !These if constructs check if the sizes were defined in the compiling. Helps with bash scripting
!     read(dummy,*) nx
!   else
!     print *, "Please enter a value for nx"  !If the values necessary aren't provided, the code makes sure to ask for them to prevent breakage
!     read *, nx                          
!   end if

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section allocates the values to be what they need to be at the start of the code
  do nx = 1, 15
  n = nx * 2                          !The total number of atoms we are dealing with
  allocate(unitcell(n))            !Unitcell gives a visual representation of the unit cell we are dealing with
  allocate(H(n,n))                     !The Hamiltonian matrix will be the one we need to solve at the end
  allocate(nn(n,3))                    !Each atom can have 4 nearest neighbours
  allocate(VL(n,n))
  allocate(VR(n,n))
  allocate(W(n))     
  allocate(Rwork(2**n))
  allocate(parity(n))
  do i = 1, n
      unitcell(i) = i                !Give each element in unit cell a distinct number to help classify each atom
  end do
  print *, unitcell

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section will focus on correctly filling the nearest neighbour array
  !The array will be filled as follows : nearestneighbours = (leftneighbour, rightneighbour, outsideneighbour)
  c = 1
  parity(1) = -1                  !The parity determines whether the outsideneighbour is in the above cell (parity = 1) or below cell (parity = -1)
  do i = 2, n-1, 2
    parity(i) = c
    parity(i+1) = c
    parity = parity * (-1)
  end do
  parity(n) = c
  
  do i = 1, n                       !The left neighbour = i-1, unless i = 1, in which case there is no left neighbour
    if (i == 1) then
      nn(i,1) = 0
    else
      nn(i,1) = unitcell(i-1)
    end if
    
    if (i == n) then                !The right neighbour = i+1, unless i = n, in which case there is no right neighbour
      nn(i,2) = 0
    else
      nn(i,2) = unitcell(i+1)
    end if
    
    if (mod(i,2) == 0) then           !The even atoms outside neighbour has index i-1, and the odd atomss have neighbours in index i+1
      nn(i,3) = (-1) * unitcell(i-1)
    else if (mod(i,2) == 1) then
      nn(i,3) = (-1) * unitcell(i+1)
    end if
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !The k loop starts here
  open(unit = outunit, file = "Band.dat", status = "replace", action = "write", iostat = istat, iomsg = errormsg)
  k = 0
  !do while (k <= pi)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !This section is for filling the relevant hamiltonian matrix
  
  H = 0                                         !Make sure H is filled with 0s to start with
  do c = 1, n                                    !c goes from 1 to 6, one line for each atom in the unit cell
    do i = 1, 3                                    !i goes from 1 to 4, for each of the possible nearest neighbours
      if (nn(c,i) == 0) then                       !if nn(c,i) is 0 there was no valid nearest neighbour for this column, so we skip past it
        cycle                          
      else if (nn(c,i) > 0) then                   !If a nearest neighbour was present, we put the hopping energy term in the relevant column in the hamiltonian
        H(c,nn(c,i)) = H(c,nn(c,i)) - t            
      else if (nn(c,i) < 0) then                   !If the nearest neighbour was negative, it was in a different cell, so we have to include the exponential term
        temp = int(abs(nn(c,i)))
        te = t * exp(im * k * parity(c))
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
  
  write(outunit,*) k, nx, 2 * real(W(nx))
!   k = k + 0.01
  
  deallocate(unitcell)            !Unitcell gives a visual representation of the unit cell we are dealing with
  deallocate(H)                     !The Hamiltonian matrix will be the one we need to solve at the end
  deallocate(nn)                    !Each atom can have 4 nearest neighbours
  deallocate(VL)
  deallocate(VR)
  deallocate(W)     
  deallocate(Rwork)
  deallocate(parity )
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
