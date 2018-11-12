module HaloTYPE

  type Halo

     integer  :: HN                       ! number of particles in the halo                  
     real(8) :: Hmass                    ! Halo Mass fof=0.2 [10**10 Mo/h]

     real(8), dimension(3) :: HPos       ! centre of mass position [Mpc/h]
     real(8), dimension(3) :: HVel       ! centre of mass velocity [km/s]

     real(8) :: Hrad                     ! radius of halo = mean of last 4 part.
     real(8) :: Hradhalf                 ! Half-mass radius
     real(8) :: HDisp1D                  ! 1D velocity dispersion [km/s]

     real(8), dimension(6) :: HInert     ! Inertia matrix: 
                                          ! I11, I22, I33, I12, I13, I23
  end type Halo

end module HaloTYPE

program power
USE HaloTYPE
USE eval           ! Frequently used Routines for Correlation and Power Spectrum
USE snapshot       ! Reads in Dark Matter Particles 
implicit none



! All that in and out Files
character*200 :: Folder
character*200 :: InFileBase
character*200 :: PsFileBase, PsOutFile,Extension
character*200 :: datadir,HalFileBase,HalInFile,HaloFileBase
character*10 :: nodestr,snapstr,ncstr

integer :: idat                ! =1 read data; =0 get header
integer :: iPos                ! =1 read positions; =0 do not.
integer :: iVel                ! =1 read velocities; =0 do not.
integer :: iID                 ! =1 read the ID's; =0 do not.
integer :: iRed                ! =1 no redshift space distortions
                                 ! = (2,3,4) apply distortion in (x-,y-,z-) 
type(snapFILES) :: sF
type(halo) :: HaloFINAL	

! Switches
integer :: doCorrect						! do CIC correction or not
integer :: NodeNumber
integer :: FileNumber
	
! Box dimensions
real(8), dimension(:,:,:), allocatable :: deltadm,deltag
real(8), dimension(:,:,:), allocatable :: deltar

integer, dimension(NkBins) :: kBinCnt,biequicnt
real(8), dimension(NkBins) :: kTrueBinC
real(8), dimension(NkBins) :: powerdfdf
real(8), dimension(5,NkBins) :: biequihhh,biequisn

integer, dimension(NkBins,NkBins,NmuBins) :: bispectcnt
real(8), dimension(NkBins,NkBins,NmuBins) :: bispectdmdmdh
real(8), dimension(NkBins,NkBins,NmuBins,3) :: bispectk

real :: time1,time2
integer :: OMP_GET_NUM_THREADS,TID,OMP_GET_THREAD_NUM

integer :: i

integer :: NHalBin,iMBin,cnt,Nhaltot
integer, parameter :: NMassBins=5
real, dimension(NMassBins+1) :: MBinB
integer, dimension(NMassBins+1) :: NBinB
real, dimension(NMassBins) :: MBinC
integer, dimension(NMassBins) :: NHalosBin
integer, dimension(:), allocatable :: IndM
real(4), dimension(:), allocatable :: HalM
real(8), dimension(:,:), allocatable :: HalPos,HalVel
real(8), dimension(5) :: b1vec
real(8), parameter :: gr=0.013134836482296224d0


NbinB=(/20,60,180,540,1620,4860/) ! Old bins	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
call cpu_time(time1)
! Size of FFT Grid	
NCell=64

b1vec=(/6.484850908844908002d-02,3.118385029929335861d-01,8.176368939457424601d-01,1.687716449102554028d+00,3.171918504181470144d+00/)
box=1500.
! Redshift Output
FileNumber=9
! Realization
NodeNumber=1

call genbink(0.003d0,0.5d0,'log')
call genbinmu()

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


datadir='/slow/space/cosmos/lss-nongauss/baldauf/SimSuite/wmap7_fid_run'//trim(nodestr)//'/DATA/reducedGROUPS/'
!HalInFile=trim(datadir)//'GRPS_TrBack002_CMS_wmap7_fid_run'//trim(nodestr)//'_'//trim(snapstr)
HalInFile=trim(datadir)//'GRPS_CMS_wmap7_fid_run'//trim(nodestr)//'_'//trim(snapstr)
write(*,'(a)') HalInFile
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




open(25,file='/slow/space/cosmos/lss-nongauss/baldauf/SimSuite/wmap7_fid_run'//trim(nodestr)//'_ZELD/ICs/delta_64.dat',form='unformatted',status='old')
read(25) deltadm
close(25)

deltadm=deltadm/gr

allocate(deltar(1:NCell+2,1:NCell,1:NCell))
deltar=0.
call comprod(deltadm,deltadm,deltar)
call powerspectrum(deltar,kBinCnt,kTrueBinC,powerdfdf)	
deallocate(deltar)

! Write Power Spectrum Output
open(20,file=PsOutFile,form='formatted',status='replace')
    do i=1,NkBins
	write(20,'(2ES30.10,1I30,1ES30.10)') kbinc(i),ktruebinc(i),kbincnt(i),powerdfdf(i)
    enddo
close(20)


	do iMBin=1,5
		write(*,'(a)') '\n\n ***********************************************'
		write(*,'(a)') '          > > > Halos Mass-Binned < < <'
		write(*,'(a)') ' ***********************************************\n\n'
		
		open(11,file=trim(HalInFile),status='old',form='unformatted')
		read(11) NHalTot
		
		print*,'NHalTot',NHalTot
		cnt=0
		do i=1,NHalTot
			read(11) HaloFINAL
			if (HaloFINAL%HN>=NBinB(iMBin) .and. HaloFINAL%HN<NBinB(iMBin+1)) then
					cnt=cnt+1
			endif
		end do
		close(11)
		NHalBin=cnt
		print*,'NHalBin',NHalBin

		allocate(HalM(NHalBin))
		allocate(HalPos(3,NHalBin))
		allocate(HalVel(3,NHalBin))
		
		
		
		open(11,file=trim(HalInFile),status='old',form='unformatted')
		read(11) NHalTot
		cnt=0
		do i=1,NHalTot
			read(11) HaloFINAL
			if (HaloFINAL%HN>=NBinB(iMBin) .and. HaloFINAL%HN<NBinB(iMBin+1)) then
			    cnt=cnt+1
			    HalM(cnt)=HaloFINAL%HMass
			    HalPos(:,cnt)=HaloFINAL%HPos(:)
			    HalVel(:,cnt)=HaloFINAL%HVel(:)/sqrt(0.0d0+1.0d0)
			endif
		end do
		close(11)
		
		print*,'Haloes Loaded'
		
			allocate(deltag(1:NCell+2,1:NCell,1:NCell))
			!allocate(thetag(1:NCell+2,1:NCell,1:NCell))
			
			deltag=0.d0
			call cicmass(NHalBin,real(HalPos),deltag)
			MBinC(iMBin)=sum(HalM)/real(NHalBin)

			call normalize(deltag)
			

			call fft3d(deltag,deltag,'f')
			
		
			if (doCorrect==1) then
 				call ciccorrect(deltag,1.0d0)
			endif
		
			print*,'Calling Bispectrum routine'
			!call bispectrum_mu_threefield(deltadm,deltadm,deltag,.False.,20,bispectcnt,bispectdmdmdh,bispectk)

			!call bispectmuhdf5export('Spectra/bispect_hr_NODE'//trim(nodestr)//'.h5','bispect',bispectcnt,bispectdmdmdh,bispectk)


			!do i=1,nkbins
				!write(*,'(1i10,3f12.2,1es20.4)'), bispectcnt(i,i,i),bispectk(i,i,i,:)/((2.0*pi)/box),bispectdmdmdh(i,i,i)
			!enddo
			call bispectrum_fft_equi(deltag-0*b1vec(imbin)*deltadm,biequicnt,biequihhh(imbin,:))
			call bispectrum_fft_equi(deltag-b1vec(imbin)*deltadm,biequicnt,biequisn(imbin,:))
					
			deallocate(deltag)
			
			deallocate(HalM)
			deallocate(HalVel)
			deallocate(HalPos)
	
		enddo
	




				
deallocate(deltadm)
open(21,file='SNbispect/bsn_run'//trim(nodestr)//'.dat',status='replace',form='formatted')
do i=1,nkbins
	write(21,'(20es20.5)') kbinc(i),biequihhh(:,i),biequisn(:,i)
enddo
close(21)

call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') 'All done in', time2


end program power

