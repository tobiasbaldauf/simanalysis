program power

USE eval           ! Frequently used Routines for Correlation and Power Spectrum
USE snapshot       ! Reads in Dark Matter Particles 
implicit none


! All that in and out Files
character*200 :: Folder
character*200 :: InFileBase
character*200 :: PsFileBase, PsOutFile,Extension
character*10 :: nodestr,snapstr,ncstr

integer :: idat                ! =1 read data; =0 get header
integer :: iPos                ! =1 read positions; =0 do not.
integer :: iVel                ! =1 read velocities; =0 do not.
integer :: iID                 ! =1 read the ID's; =0 do not.
integer :: iRed                ! =1 no redshift space distortions
                                 ! = (2,3,4) apply distortion in (x-,y-,z-) 
type(snapFILES) :: sF
	
! Switches
integer :: doCorrect						! do CIC correction or not
integer :: NodeNumber
integer :: FileNumber
	
! Box dimensions
real(8), dimension(:,:,:), allocatable :: deltadm
real(8), dimension(:,:,:), allocatable :: deltar

integer, dimension(NkBins) :: kBinCnt
real(8), dimension(NkBins) :: kTrueBinC
real(8), dimension(NkBins) :: powerdfdf
real :: time1,time2
integer :: OMP_GET_NUM_THREADS,TID,OMP_GET_THREAD_NUM

integer :: i




	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
call cpu_time(time1)
! Size of FFT Grid	
NCell=256
box=1500.
! Redshift Output
FileNumber=9
! Realization
NodeNumber=1

call genbink(0.003d0,1.0d0,'lin')


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
!				> > > FILENAMES < < <
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (FileNumber<10) then	
    write(snapstr,'(a2,i1)') '00',FileNumber
else
    write(snapstr,'(a1,i2)') '0',FileNumber
endif

!simulation
if (NodeNumber<10) then
    write(nodestr,'(i1)') NodeNumber
else
    write(nodestr,'(i2)') NodeNumber
endif

call getarg(1,nodestr)

if (Ncell<100) then
    write(ncstr,'(i2)') Ncell
elseif (Ncell<1000) then
    write(ncstr,'(i3)') Ncell
elseif (Ncell<10000) then
    write(ncstr,'(i4)') Ncell
endif

doCorrect=1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Output File Names
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
Extension='NODE'//trim(nodestr)//'_'//assign//'C1_'//trim(ncstr)//'_'
Folder='/home/cosmos/users/dc-bald1/Projects/AnalysisCode/Spectra/'
		
PsFileBase='power_'
PsOutFile=trim(Folder)//trim(PsFileBase)//trim(Extension)//trim(snapstr)//'.dat'
write(*,*) 'Writing Power-Spectrum to ',PsOutFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
! First FFT - Determines FFT strategy
allocate(deltar(1:NCell+2,1:NCell,1:NCell))
call fft3d(deltar,deltar,'f')	
deallocate(deltar)

	
!///////////////////////////////////////////////////////////////////////////////
!__________________________________DARK MATTER__________________________________
!///////////////////////////////////////////////////////////////////////////////		

write(*,'(a)') '\n\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
write(*,'(a)') '          > > > Loading dark matter < < <'
write(*,'(a)') ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'
		


allocate(deltadm(1:NCell+2,1:NCell,1:NCell))

        
deltadm=0.

iDat=1; iPos=0; iVel=0; iID=0; iRed=0


! Fiducial Cosmology Final
sf%SNAPDATA = '/slow/space/cosmos/lss-nongauss/baldauf/SimSuite/wmap7_fid_run'//trim(nodestr)//'_PM/DATA/'
sf%SNAPBASE = 'wmap7_fid_run'//trim(nodestr)//'_PM_'
sf%SNAPEXT = snapstr
sf%SNAPNFILES = 24


call read_snap(sF,iDat,iPos,iVel,iID,iRed,NCell,deltadm)
call normalize(deltadm)


call fft3d(deltadm,deltadm,'f')


if (doCorrect==1) then
    call ciccorrect(deltadm,1.0d0)
endif



    allocate(deltar(1:NCell+2,1:NCell,1:NCell))
    deltar=0.
    call comprod(deltadm,deltadm,deltar)
    call powerspectrum(deltar,kBinCnt,kTrueBinC,powerdfdf)	
    deallocate(deltar)

! Write Output
open(20,file=PsOutFile,form='formatted',status='replace')
    do i=1,NkBins
	write(20,'(2ES30.10,1I30,1ES30.10)') kbinc(i),ktruebinc(i),kbincnt(i),powerdfdf(i)
    enddo
close(20)
					
deallocate(deltadm)



call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') 'All done in', time2


end program power

