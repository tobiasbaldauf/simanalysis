
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! module for reading in the GADGET-2 data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module snapshot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! numerical recipes
  !USE nrtype
	USE eval
! defined type snapINFO

  type :: snapFILES
     integer       :: SNAPNFILES    ! Number of files per snapshot
     character*200 :: SNAPDATA      ! Directory of snapshot data
     character*200 :: SNAPBASE      ! BASE of the file names
     character*200 :: SNAPEXT       ! Extension for this snapshot
  end type snapFILES

! --- GADGET VARIABLES FOR SNAPSHOT THAT WILL BE GLOBAL --- 

! --- GADGET VARIABLES FOR SNAPSHOT THAT ARE LOCAL --- 

  integer(4) ::  NN, Ntot  

  integer  :: Nall(0:5),Npart(0:5),FlagSfr,FlagFeedback,FlagCooling
  integer  :: FlagStellarAge,FlagMetals,NallHW(0:5),flagentrics
  integer  :: NpF,NumFiles
  real(4) :: unused(15)
  real(8) :: Massarr(0:5)
  real(8) :: pmass,afactor,redshift
  real(8) :: Omega0,OmegaL0,hlittle,BoxSize

!--------------------------  
! Global data arrays
  
  real(4), allocatable    :: PartPos(:,:),PartVel(:,:)
  integer(4), allocatable :: PartIDs(:) 

!--------------------------

! Memory buffering

  integer :: iBuffer

  integer*8 :: Nalloc

  real :: BufferFac=0.35

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

  !subroutine read_snap(sF,idat,iPos,iVel,iID,iRed,NCells,ddm,ddvx,ddvy,ddvz,ddv2x,ddv2y,ddv2z,ddvxy,ddvyz,ddvzx)
  subroutine read_snap(sF,idat,iPos,iVel,iID,iRed,NCells,ddm)

!-------------------------

    ! sF    ==  This is the structure of the gadget snapshot file names
    ! idat  ==  Read data == 1 or Read Header Only == 0 and return  
    ! iPos  ==  Read postions == 1 or not ==0
    ! iVel  ==  Read velocities == 1 or not ==0
    ! iID   ==  Read particle IDs ==1 or not ==0
    ! iRed  ==  Redshift space: =0 Real Space; (1--3) distort in x, y, z, respectively
    ! iRan  ==  Generate a Random Catalogue == 0 NO; ==1 YES; ==2 
    ! rsamp ==  sampling rate 

!-------------------------

! MODULES
  
! numerical recipes
      
    !USE nrtype
		USE eval
!------------------------

! DECLARE SNAPSHOT VARAIBLES

    implicit none

! snap file info

    type(snapFILES) :: sf

! local file variables

    integer :: idat  ! if =1 read data, if =0 return.
    integer :: iPos  
    integer :: iVel
    integer :: iID
    integer :: iRed
    integer :: iRan

    integer       :: SNAPNFILES    ! Number of files per snapshot
    character*200 :: SNAPDATA      ! Directory of snapshot data
    character*200 :: SNAPBASE      ! BASE of the file names
    character*200 :: SNAPEXT       ! Extension for this snapshot

    integer(4), parameter :: FILES = 256          ! number of files per snapshot
    character*260 :: filename(FILES), file1, file2

    integer(4) ::  i, ii, kk, noffset

!--------------------------  
! local particle data arrays
  
    real(4), allocatable    :: pos(:,:),vel(:,:)
    integer, allocatable    :: IDs(:)
		
		integer :: NCells
		real(8), dimension(:,:,:) :: ddm
		
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! diagnostics

    real(4) :: posmax1, posmax2, posmax3
    real(4) :: posmin1, posmin2, posmin3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! SNAPSHOT OF INTEREST

    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(*,*) '              >>> GET SNAPSHOT <<<'
    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! make some local variables

    SNAPNFILES=sF%SNAPNFILES
    SNAPDATA=sF%SNAPDATA
    SNAPBASE=sF%SNAPBASE
    SNAPEXT=sF%SNAPEXT

!print*, SNAPEXT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! find total number of particles of each type summing over snapsots

    do i=1,SNAPNFILES
         
       file1=trim(SNAPBASE)//trim(SNAPEXT)//'.'
       
       if(i.gt.1) then
          ii=int(log10(1.0*(i-1)))+1
       else
          ii=1
       endif
       
       SELECT CASE(ii)
       CASE(1) 
          write(file2,'(a,i1)') trim(file1),i-1
       CASE(2)
          write(file2,'(a,i2)') trim(file1),i-1
       CASE(3)
          write(file2,'(a,i3)') trim(file1),i-1
       END SELECT
       
       filename(i)=trim(SNAPDATA)//trim(file2)
 !      write(*,'(a)') filename(i)
         
    enddo
  
    write(*,*) 'opening file:='
    write(*,'(a)') trim(filename(1))
      
! now, read in the header
  
    open (1, file=filename(1), form='unformatted', status='old')
    read (1) Npart, Massarr ,afactor, redshift, FlagSfr, FlagFeedback, Nall, &
         FlagCooling,NumFiles,BoxSize,Omega0,OmegaL0,hlittle
    close (1)
    
    if(NumFiles.ne.SNAPNFILES) then
       write(*,*) 'FUNNY! Number of files per snapshot does not match input!'
       stop
    endif
      
    write(*,*) 
    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(*,*) 'Total number of particles of each type:='
    write(*,'(6(2x,I14))') Nall(0:5)
    write(*,*) 'Snapshot Number of particles of each type'
    write(*,'(6(2x,I14))') Npart(0:5)
    write(*,*) 'Masses of particles of each type'
    write(*,'(6(2x,e14.7))') Massarr(0:5) 
    write(*,'(a,2x,f14.7,2x,a,2x,f14.7)') 'aexp:=',afactor,'redshift:=',redshift 
    write(*,*) 'BOXSIZE:=',BoxSize
    write(*,*) 'om_m0, om_DE0, h',Omega0,OmegaL0,hlittle
    write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(*,*) 

!-----------------------------------------
! if we only need buffer info then return
      
    if(idat.eq.0) return

! Allocate the memory to collect all particle info from the snapshots

    Ntot=Nall(1)
    
    if(iPos==1) then
       allocate(PartPos(1:3,Ntot))
    endif
    if(ivel==1) then
       allocate(PartVel(1:3,Ntot))
    endif
    if(iID==1) then
       allocate(PartIds(1:Ntot))
    endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

! Now we just read in the coordinates, and only keep those of type=1
! The array PartPos(1:3,1:Ntot) will contain all these particles

    noffset=0
  
    do i = 1, SNAPNFILES
         
       write(*,*) 'reading snapshot...  ',i
       write(*,'(a)') trim(filename(i))
         
! --- now, read in the file ---

       open (1, file=filename(i), form='unformatted', status='old')
       read (1) Npart, Massarr ,afactor, redshift, FlagSfr, FlagFeedback, Nall, &
            FlagCooling,NumFiles,BoxSize,Omega0,OmegaL0,hlittle
                
       NN=sum(Npart)
         
       write(*,*) 'reading N particles',NN
        !write(*,*) 'reading N particles',Npart(1)
       allocate(pos(1:3,1:NN))
			allocate(vel(1:3,1:NN))
     
       read (1) pos
       read (1) vel
     
       close (1)

!--- rescale velocities so that they are comoving peculiar velocities [km/s]

       vel = vel*(sqrt(afactor))
			if (doRed==1) then
			!!	call rshift4(pos,vel)
			  pos(3,:)=pos(3,:)+vel(3,:)/100.0
				write(*,*) 'Redshift Space'
          WHERE( pos(3,:) >= BoxSize ) pos(3,:)=pos(3,:)-BoxSize
          WHERE( pos(3,:) <  BoxSize ) pos(3,:)=pos(3,:)+BoxSize
			endif
			call cicmass(NN,pos(1:3,1:NN),ddm)
			!call cicvel2(NN,dble(pos(1:3,1:NN)),dble(vel(1:3,1:NN)),ddm,ddvx,ddvy,ddvz,ddv2x,ddv2y,ddv2z,ddvxy,ddvyz,ddvzx)
			!call cicvel(NN,dble(pos(1:3,1:NN)),dble(vel(1:3,1:NN)),ddm,ddvx,ddvy,ddvz)
!--- move snapshot data into arrays ---

!-- Real or redshift space data required ---
   
!   if ired = 0 => real space positions        
!   if ired = 1 => x-axis distortion
!   if ired = 2 => y-axis distortion 
!   if ired = 3 => z-axis distortion

!   pos_s = pos_r + vel/H0

       SELECT CASE(ired)
       CASE(1)
          
          pos(1,:)=pos(1,:)+vel(1,:)/100.0

          WHERE( pos(1,:) >= BoxSize ) pos(1,:)=pos(1,:)-BoxSize
          WHERE( pos(1,:) <  BoxSize ) pos(1,:)=pos(1,:)+BoxSize
          
       CASE(2)

          pos(2,:)=pos(2,:)+vel(2,:)/100.0

          WHERE( pos(2,:) >= BoxSize ) pos(2,:)=pos(2,:)-BoxSize
          WHERE( pos(2,:) <  BoxSize ) pos(2,:)=pos(2,:)+BoxSize

       CASE(3)

          pos(3,:)=pos(3,:)+vel(3,:)/100.0

          WHERE( pos(3,:) >= BoxSize ) pos(3,:)=pos(3,:)-BoxSize
          WHERE( pos(3,:) <  BoxSize ) pos(3,:)=pos(3,:)+BoxSize

       END SELECT

!---Move local data to particle data arrays------------

       if(iPos==1) then
          PartPos(1:3,1+noffset:noffset+Npart(1))=pos(1:3, 1 + Npart(0):sum(Npart(0:1)))
       endif
       if(iVel==1) then
          PartVel(1:3,1+noffset:noffset+Npart(1))=vel(1:3, 1 + Npart(0):sum(Npart(0:1)))
       endif
       if(iID==1) then
          PartIds(1+noffset:noffset+Npart(1))=IDs(1 + Npart(0):sum(Npart(0:1)))
       endif
     
       noffset=noffset+Npart(1)
     
       write(*,*) 'freeing pos, vel, IDs'

       deallocate(pos,vel)
     
    end do

    print *,'read snapshot data'

    pmass=maxval(Massarr(:))

    write(*,*) 'mass of particles:=',pmass

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

! At this point, the coordinates of all particles of type 1 will
! be stored in the array PartPos(,)

    write(*,*) 'Total number of particles found in sub snapshots:=',noffset
    write(*,*) 'Total number of particles expected:=',Ntot

!    posmax1=maxval(PartPos(1,:)) ;  posmin1=minval(PartPos(1,:))
 !   posmax2=maxval(PartPos(2,:)) ;  posmin2=minval(PartPos(2,:))
  !  posmax3=maxval(PartPos(3,:)) ;  posmin3=minval(PartPos(3,:))
      
  !  write(*,*) 'max positions',posmax1,posmax2,posmax3
  !  write(*,*) 'min positions',posmin1,posmin2,posmin3
    write(*,*) 

  end subroutine read_snap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module snapshot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!loglog(2*pi/1500*HaloMass_NODE1_CICC0_SN0_DMP_512dpar_RES1_013(:,1),HaloMass_NODE1_CICC0_SN0_DMP_512dpar_RES1_013(:,2),'x')
