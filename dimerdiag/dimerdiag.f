      Program diagonalizer 
         
      implicit none
         
!     Indexes
      integer                                     :: i,j,k,l,m,n 
       
!     Species.                     
      character(2), allocatable, dimension(:)     :: spec_id   ! Name of species (eg: Ru...)
      integer,      allocatable, dimension(:)     :: spec_no   ! Number of ions per species
      logical                                     :: spec_flag 
      integer                                     :: spec_ts   ! Total species
      integer                                     :: spec_ti   ! Total ions
        
!     Coordinates and degrees of freedom        
      double precision,allocatable,dimension(:,:) :: coor  
      character(1),    allocatable,dimension(:,:) :: free 
!     logical,         allocatable,dimension(:,:) :: free 
        
!     Diagonalization variables         
      integer                                     :: info,lwork,lwmax
      double precision,allocatable,dimension(:)   :: W,WK
      double precision,allocatable,dimension(:,:) :: H,P,dd
      
        
! # # 1.1 - Count number of species in poscar
      open(11,file="diagonalizer_poscar.tmp",err=999)
         
!     No important data
      read(11,*,err=999,end=999) ! title
      read(11,*,err=999,end=999) ! scaling
      read(11,*,err=999,end=999) ! b1
      read(11,*,err=999,end=999) ! b2
      read(11,*,err=999,end=999) ! b3
      read(11,*,err=999,end=999) ! name of species
        
      allocate(spec_no(200))
      spec_no = 0
      read(11,*,err=500) spec_no
  500 i=0 ; j=0 ; spec_flag = .FALSE.  ! Initializing flag and counters
      do while (spec_flag.eqv..FALSE.)
       if (spec_no(i+1)>0) then
        i=i+1          ! Counts number of species
        j=j+spec_no(i) ! Counts number of atoms
       else
        spec_flag = .TRUE.
       end if
      end do
      deallocate(spec_no)
      allocate(spec_no(i),spec_id(i))
      spec_ts = i      ! Total species
      spec_ti = j      ! Total ions
!     write(*,*) spec_ts, spec_ti    
       
      close(11)        
         
!     1.2 - Get and count degrees of freedom      
      open(11,file="diagonalizer_poscar.tmp",err=999)
         
!     No important data
      read(11,*,err=999,end=999) ! title
      read(11,*,err=999,end=999) ! scaling
      read(11,*,err=999,end=999) ! b1
      read(11,*,err=999,end=999) ! b2
      read(11,*,err=999,end=999) ! b3
      read(11,*,err=999,end=999) ! name of species
      read(11,*,err=999,end=999) ! number of each species
      read(11,*,err=999,end=999) ! Selective
      read(11,*,err=999,end=999) ! Direct/Cartesiannumber of species species
        
      allocate(coor(spec_ti,1:3),free(spec_ti,1:3),dd(spec_ti,1:3))
        
      do i=1,spec_ti
       read(11,*,err=999,end=999) coor(i,1:3),free(i,1:3)
!      write(*,*) coor(i,:) , free(i,:)
      end do
        
      close(11) 
         
!     1.3 - Calculate total number of degrees of freedom 
      m=0                      
      do i=1,spec_ti ; do j=1,3  
       if (free(i,j).eq."T") then
        m=m+1                   
       end if                    
      end do         ; end do    
         
!     1.4 - Compare the number of degrees of freedom with the Hessian matrix
           
      open(12,file="diagonalizer_matrix.tmp",err=998)
      read(12,*,err=998,end=998) n
               
      if (n/=m) then
       write(*,*) "The POSCAR has ",m,"degrees of freedom."
       write(*,*) "While the Hessian has ",n,"degrees of freedom."
       go to 997
      end if
                  
!     1.5 - Allocate and read Hessian matrix
      
      lwmax=3*n-1
      allocate(H(n,n),P(n,n),W(n),WK(lwmax))
          
      do i=1,n    
       read(12,*,err=998,end=998) H(i,:)
      end do  
         
      close(12) 
      
!     2.1 - Make symmetric the Hessian matrix and store it in P
      do i=1,n ; do j=1,n   
      P(i,j)=0.5D0*(H(i,j) + H(j,i)) 
      end do   ; end do
         
!     2.1 - Determining optimun lwork previous to diagonalization  
      lwork=-1
      call dsyev('V','U',n,P,n,W,WK,lwork,info)
      lwork=min(lwmax,int(WK(1)))
          
!     2.2 - Diagonalize the matrix
      call dsyev('V','U',n,P,n,W,WK,lwork,info)
           
!     2.3 - Obtain the leading eigenvalue for the reaction (the most positive)
!     l will contain the position of the leading eigenvalue
      l=1
      do i=2,n    
       if (W(i)>W(l)) then
        l=i
       end if   
      end do  
      
!     2.4 - Assign dimer direction 
      k=0
      do i=1,spec_ti ; do j=1,3
       if (free(i,j).eq."T") then 
        k=k+1
        dd(i,j)=P(k,l)
       else
        dd(i,j)=0D0
       end if
      end do         ; end do
          
!     3.1 - Write eigenvalues and eigenvectors
      open(unit=21,file="diagonalizer_output.tmp",action="write",
     &status="unknown")
      write(21,*) 'Info: ', info
      write(21,*) 'Eigenvalues'
      write(21,*) W(:)
      write(21,*)
      write(21,*) 'Eigenvectors, as columns (P-matrix)'
      do i=1,n
       write(21,*) P(i,:)
      end do
      close(21)
          
!     3.2 - Write dimer input
      open(unit=22,file="diagonalizer_taildm.tmp",action="write",
     &status="unknown")
      do i=1,spec_ti
       write(22,"(3f23.16,3a3)") coor(i,:) , free(i,1:3)
      end do
      write(22,*)   
      do i=1,spec_ti      
       write(22,"(3f20.15)") dd(i,1:3)
      end do        
      close(22)     
      stop
              
! # # DEBUG 
999   print *, "Error reading diagonalizer_poscar.tmp"
      stop
998   print *, "Error reading diagonalizer_matrix.tmp"
      stop
997   print *, "Error: The degrees of freedom and Hessian matrix are not
     & consistent."
      stop
996   print *, "Error xxx"
      stop
995   print *, "Error xxx"
      stop
       
       
      end program diagonalizer
        
