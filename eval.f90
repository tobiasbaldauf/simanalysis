module eval
implicit none
	
integer :: NCell
integer, parameter :: NkBins=60
integer, parameter :: NrBins=60
real(8), dimension(nkbins) :: kbinc
real(8), dimension(nkbins+1) :: kbinb
real(8), parameter :: massUNIT= 1.0d10
character*3, parameter :: assign='CIC'
character*3,parameter :: ffttype='mes'
real(8) :: box
real :: rhobar
integer,parameter :: chunk=10
integer,parameter :: nthr=10
integer ,parameter :: doSN=0
integer ,parameter :: doRed=0
integer, parameter :: NAngBins=1
real(8), parameter :: pi=3.14159265d0
	
		
contains

!///////////////////////////////////////////////////////////////////////////////
!______________________________________GENBINK__________________________________
!///////////////////////////////////////////////////////////////////////////////
subroutine genbink(kmin,kmax,bintype)
implicit none
integer :: i
real(8) :: kmin,kmax,dklogbin,dklin
character*3 :: bintype

if (bintype=='log') then
dkLogBin=log10(kmax/kmin)/real(NkBins)
kBinB(1)=kmin
do i=2,NkBins+1
     kbinB(i)=log10(kmin)+(i-1)*dkLogBin
     kbinB(i)=10**kbinB(i)
     kbinC(i-1)=log10(kmin)+(real(i-1)-0.5)*dkLogBin
     kbinC(i-1)=10**kbinC(i-1)
enddo

else if (bintype=='lin') then

dklin=(kmax-kmin)/real(NkBins)
kBinB(1)=kmin
do i=2,NkBins+1
     kbinB(i)=kmin+(i-1)*dklin
     kbinC(i-1)=kmin+(real(i-1)-0.5)*dklin
enddo

else 
print*, 'Not a valid binning'
stop 
endif

end subroutine genbink

!///////////////////////////////////////////////////////////////////////////////
!______________________________________COMPROD__________________________________
!///////////////////////////////////////////////////////////////////////////////

subroutine comprod(in1,in2,out)
! computes the complex product of two grids
implicit none
real(8), dimension(:,:,:) :: in1,in2,out
integer :: i,j,k
real :: time1,time2
call cpu_time(time1)	
	out=0.
!$OMP PARALLEL DO SHARED(out,in1,in2) PRIVATE(k,j,i) NUM_THREADS(NTHR) SCHEDULE(DYNAMIC,CHUNK)
	do j=1,NCell	
  		do k=1,NCell
				do i=1,NCell+2,2
					out(i,j,k)=in1(i,j,k)*in2(i,j,k)+in1(i+1,j,k)*in2(i+1,j,k)
					out(i+1,j,k)=-in1(i,j,k)*in2(i+1,j,k)+in1(i+1,j,k)*in2(i,j,k)
				enddo
		enddo
 	enddo
!$OMP END PARALLEL DO

	call cpu_time(time2)
  time2=(time2-time1)
  write(*,'(a,f8.2)') ' Comprod done in', time2
  return	
end subroutine comprod

!///////////////////////////////////////////////////////////////////////////////
!______________________________________CICCORRECT______________________________
!///////////////////////////////////////////////////////////////////////////////

subroutine ciccorrect(GridArr,Ntot)
! Corrects for the CIC filter function
implicit none

real(8), dimension(:,:,:) :: GridArr
real(8) :: NTot,ncr,W,Wx,Wy,Wz,Wsx,Wsy,Wsz
integer :: i,j,k,nc,l,m,n
real(8) :: kx,ky,kz
real(8) :: factor,ws
real(8) :: expo
real(8), parameter :: pi=3.141592653589793d0

select case(assign)
case('NGP')
	expo=1.d0
case('CIC')
	expo=2.d0
case default
	expo=0.d0
	write( *, * ) 'Unspecified assignment Method'
end select


ncr=dble(NCell)
nc=NCell

factor=2.*pi/box

do i=1,nc+2,2
	if (i .gt. 1) then
		l=(i-1)/2
		Wx=(sin(pi*dble(l)/ncr)/(pi*dble(l)/ncr))
		Wsx=sin(pi*dble(l)/ncr)
	else
		l=1
		Wx=1.d0
		Wsx=1.0d0
	endif
	do j=1,nc	
		if (j .gt. nc/2+1) then
			m=j-nc-1
			Wy=(sin(pi*dble(m)/ncr)/(pi*dble(m)/ncr))
			Wsy=sin(pi*dble(m)/ncr)
		elseif (j .gt. 1) then
			m=j-1
			Wy=(sin(pi*dble(m)/ncr)/(pi*dble(m)/ncr))
			Wsy=sin(pi*dble(m)/ncr)
		else
			m=1
			Wy=1.d0
			Wsy=1.0d0
		endif
		do k=1,nc
			if (k .gt. nc/2+1) then
				n=k-nc-1
				Wz=(sin(pi*dble(n)/ncr)/(pi*dble(n)/ncr))
				Wsz=sin(pi*dble(n)/ncr)
			elseif (k .gt. 1) then
				n=k-1
				Wz=(sin(pi*dble(n)/ncr)/(pi*dble(n)/ncr))
				Wsz=sin(pi*dble(n)/ncr)
			else
				n=1
				Wz=1.d0
				Wsz=1.0d0
			endif
			W=Wx*Wy*Wz
			WS=(1.0d0-2.d0/3.d0*Wsx**2.0)*(1.0d0-2.d0/3.d0*Wsy**2.0)*(1.0d0-2.d0/3.d0*Wsz**2.0)
			if (W .ne. 0.) then
			if(doSN==1) then
				GridArr(i,j,k)=(GridArr(i,j,k)-1.0d0/Ntot*Ws)/W**expo
			else
				GridArr(i,j,k)=GridArr(i,j,k)/W**expo
			endif
				GridArr(i+1,j,k)=GridArr(i+1,j,k)/W**expo
			end if
			if (abs(W) < 0.0001) then
				print *, W,l,m,n
			endif
		enddo
	enddo
enddo
	write( *, * ) 'Finished assignment correction'
end subroutine ciccorrect



!///////////////////////////////////////////////////////////////////////////////
!______________________________________NORMALIZATION_____________________________
!///////////////////////////////////////////////////////////////////////////////
subroutine normalize(GridArr)
implicit none

	real(8), dimension(:,:,:) :: GridArr
	real(8) :: NTot
	real(8) :: ncr
	integer :: i,j,k,nc
	real(8) :: factor
	real :: time1,time2
	integer, parameter :: chunkloc=4
  call cpu_time(time1)				
	
	ncr=dble(NCell)
	nc=NCell
	NTot=0.d0


!$OMP PARALLEL DO SHARED(GridArr,nc) PRIVATE(k,j,i)  NUM_THREADS(chunk) SCHEDULE(DYNAMIC,CHUNKloc) &
!$OMP REDUCTION(+:NTot) 
	do j=1,nc	
		do k=1,nc
				do i=1,nc
				NTot=NTot+GridArr(i,j,k)
			enddo
		enddo
	enddo

!$OMP END PARALLEL DO



	write(*,'(a,F20.1)') ' Total Number assigned: ',Ntot
	
factor=1.d0/NTot*ncr*ncr*(ncr)
!$OMP PARALLEL DO SHARED(GridArr,NTot,factor,nc) PRIVATE(k,j,i)  NUM_THREADS(chunk) SCHEDULE(DYNAMIC,CHUNKloc)
	do j=1,nc	
		do k=1,nc
			do i=1,nc	
				GridArr(i,j,k)=GridArr(i,j,k)*factor-1.0d0
			enddo
		enddo
	enddo
!$OMP END PARALLEL DO
			

NTot=0.d0

!$OMP PARALLEL DO SHARED(GridArr,nc) PRIVATE(k,j,i)  NUM_THREADS(chunk) SCHEDULE(DYNAMIC,CHUNKloc) &
!$OMP REDUCTION(+:NTot) 		
	do j=1,nc	
		do k=1,nc
				do i=1,nc
				NTot=NTot+GridArr(i,j,k)
			enddo
		enddo
	enddo
!$OMP END PARALLEL DO
					


					
	write(*,'(a,ES8.1)') ' Check if sum=0: ',NTot
	call cpu_time(time2)
  time2=(time2-time1)
  write(*,'(a,f8.2)') ' Summation done in', time2
  return	
	
end subroutine normalize

!///////////////////////////////////////////////////////////////////////////////
!______________________________________CIC MASS__________________________________
!///////////////////////////////////////////////////////////////////////////////
subroutine cicmass(NPart,PosArr,GridArr)
	implicit none
	integer :: NC, NPart
	real(4), dimension(3,NPart) :: PosArr
	real(8), dimension(:,:,:) :: GridArr
  real(8) :: ncr
	real(8) :: factor
  integer :: i
  integer :: i1,i2,j1,j2,k1,k2
	real(8) :: x,y,z,dx1,dx2,dy1,dy2,dz1,dz2
	real :: time1,time2

  call cpu_time(time1)
	nc=ncell
	ncr=dble(ncell)
	factor=ncr/box
	
	select case(assign)
	
	case('NGP')
	
		print *,'NGP'
		do i=1,NPart
      		 
					
	x=mod(PosArr(1,i)*factor,ncr)
       y=mod(PosArr(2,i)*factor,ncr)
       z=mod(PosArr(3,i)*factor,ncr)
			 
			 i1=floor(x)+1
       j1=floor(y)+1
       k1=floor(z)+1
      
       GridArr(i1,j1,k1)=GridArr(i1,j1,k1)+1.0
  
	enddo
	
	case('CIC')
			print *,'CIC'
			
			do i=1,NPart
				x=mod(PosArr(1,i)*factor+ncr,ncr)
				y=mod(PosArr(2,i)*factor+ncr,ncr)
				z=mod(PosArr(3,i)*factor+ncr,ncr)
	 			 
	        i1=floor(x)+1
	        i2=mod(i1,nc)+1
	        dx1=real(i1)-x
	        dx2=1.0-dx1
	        
                j1=floor(y)+1
	        j2=mod(j1,nc)+1
		dy1=real(j1)-y
	        dy2=1.0-dy1
	        
                k1=floor(z)+1
	        k2=mod(k1,nc)+1
	        dz1=real(k1)-z
	        dz2=1.0-dz1
                
				GridArr(i1,j1,k1)=GridArr(i1,j1,k1)+dx1*dy1*dz1
				GridArr(i2,j1,k1)=GridArr(i2,j1,k1)+dx2*dy1*dz1
				GridArr(i1,j2,k1)=GridArr(i1,j2,k1)+dx1*dy2*dz1
				GridArr(i2,j2,k1)=GridArr(i2,j2,k1)+dx2*dy2*dz1
				GridArr(i1,j1,k2)=GridArr(i1,j1,k2)+dx1*dy1*dz2
				GridArr(i2,j1,k2)=GridArr(i2,j1,k2)+dx2*dy1*dz2
				GridArr(i1,j2,k2)=GridArr(i1,j2,k2)+dx1*dy2*dz2
				GridArr(i2,j2,k2)=GridArr(i2,j2,k2)+dx2*dy2*dz2
		enddo
                
	case default
		write( *, * ) 'Unspecified assignment Method'
	end select
  

	
	call cpu_time(time2)
  time2=(time2-time1)
  write(*,'(a,f8.2)') ' Assignment done in', time2
   return
	
end subroutine cicmass



!///////////////////////////////////////////////////////////////////////////////
!______________________________________Powerspectrum_____________________________
!///////////////////////////////////////////////////////////////////////////////

subroutine powerspectrum(delta,binCnt,binK,binP)
implicit none
real(8), dimension(:,:,:) :: delta
integer :: w
integer :: i,j,k
real(8) :: kr,kx,ky,kz,pow            
real :: time1,time2
real(8), dimension(nkbins) :: binP,binK
integer, dimension(nkbins) :: binCnt

call cpu_time(time1)

binCnt=0
binP=0.d0
binK=0.d0


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,i,kz,ky,kx,kr,pow,w) NUM_THREADS(chunk) &
!$OMP REDUCTION(+:binCnt,binP,binK)		
       do k=1,ncell
          if (k .lt. ncell/2+1) then
             kz=k-1
          else
             kz=k-1-ncell
          endif
          do j=1,ncell
             if (j .lt. ncell/2+1) then
                ky=j-1
             else
                ky=j-1-ncell
             endif
             do i=1,ncell+2,2
                kx=(i-1)/2
                kr=sqrt(kx**2+ky**2+kz**2)*2.0d0*pi/box
                if (kr .ne. 0) then

			pow=delta(i,j,k)
			do w=1,nkbins
				if ((kr>kbinb(w)) .and.(kr<=kbinb(w+1))) then
					binCnt(w)=binCnt(w)+1
					binP(w)=binP(w)+pow
					binK(w)=binK(w)+kr
				endif
			enddo
                endif
             enddo
          enddo
       enddo 
!$OMP END PARALLEL DO 


do w=1,nkbins
	if (binCnt(w) .gt. 0) then
		binP(w)=1.0d0/box**3.0d0*binP(w)/real(binCnt(w))
		binK(w)=binK(w)/real(binCnt(w))
	endif
enddo



call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Power spectrum calculated in', time2
return
end subroutine powerspectrum	

!///////////////////////////////////////////////////////////////////////////////
!_________________________________________FFT __________________________________
!///////////////////////////////////////////////////////////////////////////////
subroutine fft3d(a,b,c)
implicit none

	include 'fftw3.f'
	real(8), dimension(:,:,:):: a,b
	integer(8),save :: plan,iplan
  logical :: first_fft
	
 	character c
	real :: time1,time2

  call cpu_time(time1)				
					
	
					
	!create plan
  	data first_fft /.true./
	if (first_fft) then
  	first_fft=.false.
        print*,'creating plan'
		call dfftw_init_threads
		call dfftw_plan_with_nthreads(8)
		if(ffttype=='mes') then
			write(*,*) 'Measure'
			call dfftw_plan_dft_r2c_3d(plan, ncell, ncell, ncell, a, a,FFTW_MEASURE)
			call dfftw_plan_dft_c2r_3d(iplan, ncell, ncell, ncell, a, a,FFTW_MEASURE)
		else
			write(*,*) 'Estimate'
			call dfftw_plan_dft_r2c_3d(plan, ncell, ncell, ncell, a, a,FFTW_ESTIMATE)
			call dfftw_plan_dft_c2r_3d(iplan, ncell, ncell, ncell, a, a,FFTW_ESTIMATE)
		endif
	endif


	if (c .eq. 'f') then
		write(*,*) 'Forward FFT'
		call dfftw_execute_dft_r2c(plan,a,a)
		a=a/real(NCell)**(3.)*box**3.0d0
	else
   	write(*,*) 'Inverse FFT'
		call dfftw_execute_dft_c2r(iplan,a,a)
                a=a/box**3.0d0
	endif
    
		
	call cpu_time(time2)
  time2=(time2-time1)
  write(*,*) 'FFT Calculated in',time2
  return
end subroutine fft3d



end module eval