module eval
implicit none
	
integer :: NCell
integer, parameter :: NkBins=15
integer, parameter :: NmuBins=20
integer, parameter :: NrBins=60
real(8), dimension(nkbins) :: kbinc
real(8), dimension(nkbins+1) :: kbinb
real(8), dimension(nmubins) :: mubinc
real(8), dimension(nmubins+1) :: mubinb
real(8), parameter :: massUNIT= 1.0d10
character*3, parameter :: assign='CIC'
character*3,parameter :: ffttype='est'
real(8) :: box
real :: rhobar
integer,parameter :: chunk=10
integer,parameter :: nthr=16
integer,parameter :: nthrfftw=10
integer ,parameter :: doSN=0
integer ,parameter :: doRed=0
integer, parameter :: NAngBins=1
real(8), parameter :: pi=3.14159265d0
integer, parameter :: ithresh=512/2
		
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
!______________________________________GENBINK__________________________________
!///////////////////////////////////////////////////////////////////////////////
subroutine genbinmu()
implicit none
integer :: i
real(8) :: dmu
character*3 :: bintype

dmu=2.0d0/real(NmuBins)
muBinB(1)=-1.0d0
do i=2,NmuBins+1
     mubinB(i)=-1.0d0+(i-1)*dmu
     mubinC(i-1)=-1.0d0+(real(i-1)-0.5)*dmu
enddo


end subroutine genbinmu
!///////////////////////////////////////////////////////////////////////////////
!______________________________________GENBINK__________________________________
!///////////////////////////////////////////////////////////////////////////////
subroutine genbinr(kmin,kmax,bintype)
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

end subroutine genbinr

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


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,i,kz,ky,kx,kr,pow,w) NUM_THREADS(nthr) &
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
!				Bispectrum mmh kkk
!///////////////////////////////////////////////////////////////////////////////

subroutine bispectrum_threefield(deltaa,deltab,deltac,dosym,imax,binCnt,binB,binK)
implicit none 
integer :: imax,nc
real(8), dimension(:,:,:) :: deltaa,deltab,deltac
integer :: w1,w2,w3
integer :: i,i1,i2,i3,j1,j2,j3,k1,k2,k3
real(8) :: kr,kr1,kr2,kr3,kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,pow
real :: time1,time2
real(8), dimension(:,:,:) :: binB
real(8), dimension(:,:,:,:) :: binK
integer, dimension(:,:,:) :: binCnt
!integer :: OMP_GET_NUM_THREADS,nthr,TID,OMP_GET_THREAD_NUM
real(8) :: cyfac
complex(8) :: c1,c2,c3
logical :: dosym

print*,'Bispectrum Code Started'

call cpu_time(time1)
nc=ncell



binB=0.0d0
binCnt=0
binK=0.0d0
print*,imax
print*,'Starting Loops'
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k1,k2,k3,j1,j2,j3,i1,i2,i3,kz1,kz2,kz3,ky1,ky2,ky3,kx1,kx2,kx3,kr1,kr2,kr3,pow,w1,w2,w3,cyfac,c1,c2,c3) NUM_THREADS(nthr) &
!$OMP REDUCTION(+:binCnt,binB,binK)
do k1=1,nc
	!print*,'k1',k1
	!nthr=OMP_GET_NUM_THREADS()
	!TID = OMP_GET_THREAD_NUM()
	!print *,'Threads ',nthr,'ID ',tid
	
	if (k1 .lt. imax) then
		kz1=k1-1
	else if(k1 .gt. nc-imax) then
		kz1=k1-1-nc
	else
		cycle
	endif
	do j1=1,nc
		!print*,'j1',j1
		if (j1 .lt. imax) then
			ky1=j1-1
		else if (j1 .gt. nc-imax) then
		   ky1=j1-1-nc
		else
			cycle
		endif
		do i1=1,2*imax,2
			kx1=(i1-1)/2
			do k2=1,nc
				if (k2 .lt. imax) then
					kz2=k2-1
				else if(k2 .gt. nc-imax) then
					kz2=k2-1-nc
				else
					cycle
				endif
				do j2=1,nc
					!print*,'j2',j2
					if (j2 .lt. imax) then
					   ky2=j2-1
					else if (j2 .lt. nc-imax) then
					   ky2=j2-1-nc
					else
						cycle
					endif
					do i2=1,2*imax,2
					   kx2=(i2-1)/2
					
					if (sqrt((kx1+kx2)**2.0+(ky1+ky2)**2.0+(kz1+kz2)**2.0)<imax) then
					
						kx3=-kx1-kx2
						ky3=-ky1-ky2
						kz3=-kz1-kz2
						
						if (kx3<0) then
							i3=-2*kx3+1
							ky3=-ky3
							kz3=-kz3
							cyfac=-1.0d0
						else
							i3=2*kx3+1
							cyfac=1.0d0
						endif
						if (ky3<0) then	
							j3=ky3+nc+1
						else
							j3=ky3+1
						endif
						if (kz3<0) then
							k3=kz3+nc+1
						else
							k3=kz3+1
						end if
						
						
						
						c1=cmplx(deltaa(i1,j1,k1),deltaa(i1+1,j1,k1))
						c2=cmplx(deltab(i2,j2,k2),deltab(i2+1,j2,k2))
						c3=cmplx(deltac(i3,j3,k3),cyfac*deltac(i3+1,j3,k3))
						
						pow=real(c1*c2*c3)
						
						if (dosym .eqv. .true.) then
								c1=cmplx(deltab(i1,j1,k1),deltab(i1+1,j1,k1))
								c2=cmplx(deltaa(i2,j2,k2),deltaa(i2+1,j2,k2))
								c3=cmplx(deltac(i3,j3,k3),cyfac*deltac(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
								c1=cmplx(deltaa(i1,j1,k1),deltaa(i1+1,j1,k1))
								c2=cmplx(deltac(i2,j2,k2),deltac(i2+1,j2,k2))
								c3=cmplx(deltab(i3,j3,k3),cyfac*deltab(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
								c1=cmplx(deltac(i1,j1,k1),deltac(i1+1,j1,k1))
								c2=cmplx(deltab(i2,j2,k2),deltab(i2+1,j2,k2))
								c3=cmplx(deltaa(i3,j3,k3),cyfac*deltaa(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
	
								c1=cmplx(deltac(i1,j1,k1),deltac(i1+1,j1,k1))
								c2=cmplx(deltaa(i2,j2,k2),deltaa(i2+1,j2,k2))
								c3=cmplx(deltab(i3,j3,k3),cyfac*deltab(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
	
								c1=cmplx(deltab(i1,j1,k1),deltab(i1+1,j1,k1))
								c2=cmplx(deltac(i2,j2,k2),deltac(i2+1,j2,k2))
								c3=cmplx(deltaa(i3,j3,k3),cyfac*deltaa(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
			
								pow=pow/6.0d0
						endif
						
	
						
						
						
						kr1=sqrt(kx1**2.0+ky1**2.0+kz1**2.0)*2.0d0*pi/box
						kr2=sqrt(kx2**2.0+ky2**2.0+kz2**2.0)*2.0d0*pi/box
						kr3=sqrt(kx3**2.0+ky3**2.0+kz3**2.0)*2.0d0*pi/box
						
												
						if (kr1 > 0.0 .and. kr2 > 0.0 .and. kr3 > 0.0) then
							do w1=1,Nkbins
								do w2=1,Nkbins
									do w3=1,Nkbins
									   if ((kr1>=kbinb(w1)) .and.(kr1<kbinb(w1+1)) .and. (kr2>=kbinb(w2)) .and.(kr2<kbinb(w2+1)) .and. (kr3>=kbinb(w3)) .and.(kr3<kbinb(w3+1))) then
											!print*,kr1,kr2,kr3
											binCnt(w1,w2,w3)=binCnt(w1,w2,w3)+1
											binB(w1,w2,w3)=binB(w1,w2,w3)+pow
											binK(w1,w2,w3,1)=binK(w1,w2,w3,1)+kr1
											binK(w1,w2,w3,2)=binK(w1,w2,w3,2)+kr2
											binK(w1,w2,w3,3)=binK(w1,w2,w3,3)+kr3
									   endif
									enddo
								enddo
							enddo
						endif
						
					endif
					enddo

				enddo
			enddo 
		enddo
	enddo
enddo
!$OMP END PARALLEL DO 




do w1=1,Nkbins
	do w2=1,Nkbins
		do w3=1,Nkbins
			if (binCnt(w1,w2,w3) .gt. 0) then
				
				binB(w1,w2,w3)=1.0d0/box**3.0*binB(w1,w2,w3)/real(binCnt(w1,w2,w3))
				binK(w1,w2,w3,:)=binK(w1,w2,w3,:)/real(binCnt(w1,w2,w3))
			!print*,w1,w2,w3,binCnt(w1,w2,w3),binK(w1,w2,w3,:)
			endif
		enddo
	enddo
enddo
call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Bispectrum calculated in', time2
return
end subroutine bispectrum_threefield

subroutine bispectrum_mu_threefield(deltaa,deltab,deltac,dosym,imax,binCnt,binB,binK)
implicit none 
integer :: imax,nc
real(8), dimension(:,:,:) :: deltaa,deltab,deltac
integer :: w1,w2,w3
integer :: i,i1,i2,i3,j1,j2,j3,k1,k2,k3,x2,x3,mu
real(8) :: kr,kr1,kr2,kr3,kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,pow
real :: time1,time2
real(8), dimension(:,:,:) :: binB
real(8), dimension(:,:,:,:) :: binK
integer, dimension(:,:,:) :: binCnt
!integer :: OMP_GET_NUM_THREADS,nthr,TID,OMP_GET_THREAD_NUM
real(8) :: cyfac
complex(8) :: c1,c2,c3
logical :: dosym

print*,'Bispectrum Code Started'

call cpu_time(time1)
nc=ncell



binB=0.0d0
binCnt=0
binK=0.0d0
print*,imax
print*,'Starting Loops'
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(x2,x3,mu,k1,k2,k3,j1,j2,j3,i1,i2,i3,kz1,kz2,kz3,ky1,ky2,ky3,kx1,kx2,kx3,kr1,kr2,kr3,pow,w1,w2,w3,cyfac,c1,c2,c3) NUM_THREADS(nthr) &
!$OMP REDUCTION(+:binCnt,binB,binK)
do k1=1,nc
	!print*,'k1',k1
	!nthr=OMP_GET_NUM_THREADS()
	!TID = OMP_GET_THREAD_NUM()
	!print *,'Threads ',nthr,'ID ',tid
	
	if (k1 .lt. imax) then
		kz1=k1-1
	else if(k1 .gt. nc-imax) then
		kz1=k1-1-nc
	else
		cycle
	endif
	do j1=1,nc
		!print*,'j1',j1
		if (j1 .lt. imax) then
			ky1=j1-1
		else if (j1 .gt. nc-imax) then
		   ky1=j1-1-nc
		else
			cycle
		endif
		do i1=1,2*imax,2
			kx1=(i1-1)/2
			do k2=1,nc
				if (k2 .lt. imax) then
					kz2=k2-1
				else if(k2 .gt. nc-imax) then
					kz2=k2-1-nc
				else
					cycle
				endif
				do j2=1,nc
					!print*,'j2',j2
					if (j2 .lt. imax) then
					   ky2=j2-1
					else if (j2 .lt. nc-imax) then
					   ky2=j2-1-nc
					else
						cycle
					endif
					do i2=1,2*imax,2
					   kx2=(i2-1)/2
					
					if (sqrt((kx1+kx2)**2.0+(ky1+ky2)**2.0+(kz1+kz2)**2.0)<imax) then
					
						kx3=-kx1-kx2
						ky3=-ky1-ky2
						kz3=-kz1-kz2
						
						if (kx3<0) then
							i3=-2*kx3+1
							ky3=-ky3
							kz3=-kz3
							cyfac=-1.0d0
						else
							i3=2*kx3+1
							cyfac=1.0d0
						endif
						if (ky3<0) then	
							j3=ky3+nc+1
						else
							j3=ky3+1
						endif
						if (kz3<0) then
							k3=kz3+nc+1
						else
							k3=kz3+1
						end if
						
						
						
						c1=cmplx(deltaa(i1,j1,k1),deltaa(i1+1,j1,k1))
						c2=cmplx(deltab(i2,j2,k2),deltab(i2+1,j2,k2))
						c3=cmplx(deltac(i3,j3,k3),cyfac*deltac(i3+1,j3,k3))
						
						pow=real(c1*c2*c3)
						
						if (dosym .eqv. .true.) then
								c1=cmplx(deltab(i1,j1,k1),deltab(i1+1,j1,k1))
								c2=cmplx(deltaa(i2,j2,k2),deltaa(i2+1,j2,k2))
								c3=cmplx(deltac(i3,j3,k3),cyfac*deltac(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
								c1=cmplx(deltaa(i1,j1,k1),deltaa(i1+1,j1,k1))
								c2=cmplx(deltac(i2,j2,k2),deltac(i2+1,j2,k2))
								c3=cmplx(deltab(i3,j3,k3),cyfac*deltab(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
								c1=cmplx(deltac(i1,j1,k1),deltac(i1+1,j1,k1))
								c2=cmplx(deltab(i2,j2,k2),deltab(i2+1,j2,k2))
								c3=cmplx(deltaa(i3,j3,k3),cyfac*deltaa(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
	
								c1=cmplx(deltac(i1,j1,k1),deltac(i1+1,j1,k1))
								c2=cmplx(deltaa(i2,j2,k2),deltaa(i2+1,j2,k2))
								c3=cmplx(deltab(i3,j3,k3),cyfac*deltab(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
	
								c1=cmplx(deltab(i1,j1,k1),deltab(i1+1,j1,k1))
								c2=cmplx(deltac(i2,j2,k2),deltac(i2+1,j2,k2))
								c3=cmplx(deltaa(i3,j3,k3),cyfac*deltaa(i3+1,j3,k3))
								pow=pow+real(c1*c2*c3)
			
			
								pow=pow/6.0d0
						endif
						
	
						
						
						
						kr1=sqrt(kx1**2.0+ky1**2.0+kz1**2.0)*2.0d0*pi/box
						kr2=sqrt(kx2**2.0+ky2**2.0+kz2**2.0)*2.0d0*pi/box
						kr3=sqrt(kx3**2.0+ky3**2.0+kz3**2.0)*2.0d0*pi/box
						
												
						if (kr1 > 0.0 .and. kr2 > 0.0 .and. kr3 > 0.0) then
							x2=kr2/kr1;
							x3=kr3/kr1;
							mu=(x3**2.0-1.0-x2**2.0)/(2.0*x2)
							do w1=1,Nkbins
								do w2=1,Nkbins
									do w3=1,Nmubins
									   if ((kr1>=kbinb(w1)) .and.(kr1<kbinb(w1+1)) .and. (kr2>=kbinb(w2)) .and.(kr2<kbinb(w2+1)) .and. (mu>=MuBinB(w3)) .and.(mu<=MuBinB(w3+1))) then
											!print*,kr1,kr2,kr3
											binCnt(w1,w2,w3)=binCnt(w1,w2,w3)+1
											binB(w1,w2,w3)=binB(w1,w2,w3)+pow
											binK(w1,w2,w3,1)=binK(w1,w2,w3,1)+kr1
											binK(w1,w2,w3,2)=binK(w1,w2,w3,2)+kr2
											binK(w1,w2,w3,3)=binK(w1,w2,w3,3)+kr3
									   endif
									enddo
								enddo
							enddo
						endif
						
					endif
					enddo

				enddo
			enddo 
		enddo
	enddo
enddo
!$OMP END PARALLEL DO 




do w1=1,Nkbins
	do w2=1,Nkbins
		do w3=1,Nmubins
			if (binCnt(w1,w2,w3) .gt. 0) then
				
				binB(w1,w2,w3)=1.0d0/box**3.0*binB(w1,w2,w3)/real(binCnt(w1,w2,w3))
				binK(w1,w2,w3,:)=binK(w1,w2,w3,:)/real(binCnt(w1,w2,w3))
			!print*,w1,w2,w3,binCnt(w1,w2,w3),binK(w1,w2,w3,:)
			endif
		enddo
	enddo
enddo
call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Bispectrum calculated in', time2
return
end subroutine bispectrum_mu_threefield

!///////////////////////////////////////////////////////////////////////////////
!				Bispectrum mmh Mudependent
!///////////////////////////////////////////////////////////////////////////////
subroutine bispectrum_fft_equi(delta,radcnt,radcont)
implicit none
real(8) :: ncr  
integer :: nc
real(8), dimension(:,:,:) :: delta
integer :: w1,w2,w3
integer :: i,i1,i2,i3,j1,j2,j3,k1,k2,k3
real(8) :: kr,kr1,kr2,kr3,kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,pow,x2,x3,mu       
real :: time1,time2

real(8), dimension(NkBins) :: RadCont,RadContd,RadVar
real(8), dimension(:,:,:), allocatable :: deltaa,deltab,deltac,deltaad,deltabd,deltacd
integer, dimension(NkBins) :: RadCnt
integer :: OMP_GET_NUM_THREADS,nthr,TID,OMP_GET_THREAD_NUM
real(8) :: sumsqpower
real(8) :: cyfac
integer, parameter :: nmin=1



print*, 'start bispectrum'


call cpu_time(time1)
nc=ncell
ncr=ncell


RadCnt=0
RadCont=0.
RadContd=0.
RadVar=0.

sumsqpower=0.0d0





do w1=1,NkBins
	allocate(deltaa(Nc+2,Nc,Nc))
	deltaa=0.0d0
	allocate(deltaad(Nc+2,Nc,Nc))
	deltaad=0.0d0

	do k1=nmin,nc
	!print*,'k1',k1
	!nthr=OMP_GET_NUM_THREADS()
	!TID = OMP_GET_THREAD_NUM()
	!print *,'Threads ',nthr,'ID ',tid


	if (k1 .lt. ithresh) then
		    kz1=k1-1
		else if (k1 .gt. nc-ithresh) then
		    kz1=k1-1-nc
		else
		    cycle
		endif
		do j1=nmin,nc
		    !print*,'j1',j1
		    if (j1 .lt. ithresh) then
			ky1=j1-1
		    else if (j1 .gt. nc-ithresh) then
			ky1=j1-1-nc
		    else
			cycle
		    endif
		    do i1=1,ncell+2,2
		        kx1=(i1-1)/2
			kr1=sqrt(kx1**2.0+ky1**2.0+kz1**2.0)*2.0*pi/box
			if ((kr1>kBinB(w1)) .and.(kr1<=kBinB(w1+1))) then
			    deltaa(i1,j1,k1)=delta(i1,j1,k1)
			    deltaa(i1+1,j1,k1)=delta(i1+1,j1,k1)
	    		    deltaad(i1,j1,k1)=1.0d0
			    deltaad(i1+1,j1,k1)=0.0d0
			endif
			enddo
			enddo
			enddo
			call fft3d(deltaa,deltaa,'e')
			call fft3d(deltaad,deltaad,'e')
		
	 

	print*,'Done ffts'
	!$OMP PARALLEL DO SHARED(nc) PRIVATE(k1,j1,i1)  NUM_THREADS(nthr) SCHEDULE(DYNAMIC,chunk) &
	!$OMP REDUCTION(+:RadCont,RadContD) 
		do j1=1,nc	
		    do k1=1,nc
			do i1=1,nc
			    RadCont(w1)=RadCont(w1)+deltaa(i1,j1,k1)**3.0
			    RadContD(w1)=RadContD(w1)+deltaad(i1,j1,k1)**3.0
			enddo
		    enddo
		enddo
	!$OMP END PARALLEL DO

	if (RadContD(w1)>0.0d0) then
	    RadCont(w1)=RadCont(w1)/RadContD(w1)/box**3.0
	else
	    RadCont(w1)=0.0d0
	endif


	print*,'Done sum'

	!enddo
	deallocate(deltaa)
	deallocate(deltaad)
enddo


call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Bispectrum calculated in', time2
return
end subroutine bispectrum_fft_equi


!///////////////////////////////////////////////////////////////////////////////
!				Bispectrum mmh Mudependent
!///////////////////////////////////////////////////////////////////////////////
subroutine trispectrum_fft_equi(delta,radcnt,radcont)
implicit none
real(8) :: ncr  
integer :: nc
real(8), dimension(:,:,:) :: delta
integer :: w1,w2,w3
integer :: i,i1,i2,i3,j1,j2,j3,k1,k2,k3
real(8) :: kr,kr1,kr2,kr3,kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,pow,x2,x3,mu       
real :: time1,time2

real(8), dimension(NkBins) :: RadCont,RadContd,RadVar
real(8), dimension(:,:,:), allocatable :: deltaa,deltab,deltac,deltaad,deltabd,deltacd
integer, dimension(NkBins) :: RadCnt
integer :: OMP_GET_NUM_THREADS,nthr,TID,OMP_GET_THREAD_NUM
real(8) :: sumsqpower
real(8) :: cyfac
integer, parameter :: nmin=1



print*, 'start bispectrum'


call cpu_time(time1)
nc=ncell
ncr=ncell


RadCnt=0
RadCont=0.
RadContd=0.
RadVar=0.

sumsqpower=0.0d0





do w1=1,NkBins
	allocate(deltaa(Nc+2,Nc,Nc))
	deltaa=0.0d0
	allocate(deltaad(Nc+2,Nc,Nc))
	deltaad=0.0d0

	do k1=nmin,nc
	!print*,'k1',k1
	!nthr=OMP_GET_NUM_THREADS()
	!TID = OMP_GET_THREAD_NUM()
	!print *,'Threads ',nthr,'ID ',tid


	if (k1 .lt. ithresh) then
		    kz1=k1-1
		else if (k1 .gt. nc-ithresh) then
		    kz1=k1-1-nc
		else
		    cycle
		endif
		do j1=nmin,nc
		    !print*,'j1',j1
		    if (j1 .lt. ithresh) then
			ky1=j1-1
		    else if (j1 .gt. nc-ithresh) then
			ky1=j1-1-nc
		    else
			cycle
		    endif
		    do i1=1,ncell+2,2
		        kx1=(i1-1)/2
			kr1=sqrt(kx1**2.0+ky1**2.0+kz1**2.0)*2.0*pi/box
			if ((kr1>kBinB(w1)) .and.(kr1<=kBinB(w1+1))) then
			    deltaa(i1,j1,k1)=delta(i1,j1,k1)
			    deltaa(i1+1,j1,k1)=delta(i1+1,j1,k1)
	    		    deltaad(i1,j1,k1)=1.0d0
			    deltaad(i1+1,j1,k1)=0.0d0
			endif
			enddo
			enddo
			enddo
			call fft3d(deltaa,deltaa,'e')
			call fft3d(deltaad,deltaad,'e')
		
	 

	print*,'Done ffts'
	!$OMP PARALLEL DO SHARED(nc) PRIVATE(k1,j1,i1)  NUM_THREADS(nthr) SCHEDULE(DYNAMIC,chunk) &
	!$OMP REDUCTION(+:RadCont,RadContD) 
		do j1=1,nc	
		    do k1=1,nc
			do i1=1,nc
			    RadCont(w1)=RadCont(w1)+deltaa(i1,j1,k1)**4.0
			    RadContD(w1)=RadContD(w1)+deltaad(i1,j1,k1)**4.0
			enddo
		    enddo
		enddo
	!$OMP END PARALLEL DO

	if (RadContD(w1)>0.0d0) then
	    RadCont(w1)=RadCont(w1)/RadContD(w1)
	else
	    RadCont(w1)=0.0d0
	endif


	print*,'Done sum'

	!enddo
	deallocate(deltaa)
	deallocate(deltaad)
enddo


call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Bispectrum calculated in', time2
return
end subroutine trispectrum_fft_equi


!///////////////////////////////////////////////////////////////////////////////
!				Bispectrum mmh Mudependent
!///////////////////////////////////////////////////////////////////////////////
subroutine trispectrum_fft_equi_propa(delta,deltalin,radcnt,radcont)
implicit none
real(8) :: ncr  
integer :: nc
real(8), dimension(:,:,:) :: delta,deltalin
integer :: w1,w2,w3
integer :: i,i1,i2,i3,j1,j2,j3,k1,k2,k3
real(8) :: kr,kr1,kr2,kr3,kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,pow,x2,x3,mu       
real :: time1,time2

real(8), dimension(NkBins) :: RadCont,RadContd,RadVar
real(8), dimension(:,:,:), allocatable :: deltaa,deltab,deltac,deltaad,deltabd,deltacd
integer, dimension(NkBins) :: RadCnt
integer :: OMP_GET_NUM_THREADS,nthr,TID,OMP_GET_THREAD_NUM
real(8) :: sumsqpower
real(8) :: cyfac
integer, parameter :: nmin=1



print*, 'start bispectrum'


call cpu_time(time1)
nc=ncell
ncr=ncell


RadCnt=0
RadCont=0.
RadContd=0.
RadVar=0.

sumsqpower=0.0d0





do w1=1,NkBins
	allocate(deltaa(Nc+2,Nc,Nc))
	deltaa=0.0d0
	allocate(deltaa(Nc+2,Nc,Nc))
	deltab=0.0d0
	allocate(deltaad(Nc+2,Nc,Nc))
	deltaad=0.0d0

	do k1=nmin,nc
	!print*,'k1',k1
	!nthr=OMP_GET_NUM_THREADS()
	!TID = OMP_GET_THREAD_NUM()
	!print *,'Threads ',nthr,'ID ',tid


	if (k1 .lt. ithresh) then
		    kz1=k1-1
		else if (k1 .gt. nc-ithresh) then
		    kz1=k1-1-nc
		else
		    cycle
		endif
		do j1=nmin,nc
		    !print*,'j1',j1
		    if (j1 .lt. ithresh) then
			ky1=j1-1
		    else if (j1 .gt. nc-ithresh) then
			ky1=j1-1-nc
		    else
			cycle
		    endif
		    do i1=1,ncell+2,2
		        kx1=(i1-1)/2
			kr1=sqrt(kx1**2.0+ky1**2.0+kz1**2.0)*2.0*pi/box
			if ((kr1>kBinB(w1)) .and.(kr1<=kBinB(w1+1))) then
			    deltaa(i1,j1,k1)=delta(i1,j1,k1)
			    deltaa(i1+1,j1,k1)=delta(i1+1,j1,k1)
			    deltab(i1,j1,k1)=deltalin(i1,j1,k1)
			    deltab(i1+1,j1,k1)=deltalin(i1+1,j1,k1)
	    		    deltaad(i1,j1,k1)=1.0d0
			    deltaad(i1+1,j1,k1)=0.0d0
			endif
			enddo
			enddo
			enddo
			call fft3d(deltaa,deltaa,'e')
			call fft3d(deltaad,deltaad,'e')
		
	 

	print*,'Done ffts'
	!$OMP PARALLEL DO SHARED(nc) PRIVATE(k1,j1,i1)  NUM_THREADS(nthr) SCHEDULE(DYNAMIC,chunk) &
	!$OMP REDUCTION(+:RadCont,RadContD) 
		do j1=1,nc	
		    do k1=1,nc
			do i1=1,nc
			    RadCont(w1)=RadCont(w1)+deltaa(i1,j1,k1)**3.0*deltab(i1,j1,k1)
			    RadContD(w1)=RadContD(w1)+deltaad(i1,j1,k1)**4.0
			enddo
		    enddo
		enddo
	!$OMP END PARALLEL DO

	if (RadContD(w1)>0.0d0) then
	    RadCont(w1)=RadCont(w1)/RadContD(w1)
	else
	    RadCont(w1)=0.0d0
	endif


	print*,'Done sum'

	!enddo
	deallocate(deltaa)
	deallocate(deltaad)
enddo


call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Bispectrum calculated in', time2
return
end subroutine trispectrum_fft_equi_propa

!///////////////////////////////////////////////////////////////////////////////
!				Bispectrum mmh Mudependent
!///////////////////////////////////////////////////////////////////////////////
subroutine bispectrum_fft_kkk(delta,radcnt,radcont)
implicit none
real(8) :: ncr  
integer :: nc
real(8), dimension(:,:,:) :: delta
integer :: w1,w2,w3
integer :: i,i1,i2,i3,j1,j2,j3,k1,k2,k3
real(8) :: kr,kr1,kr2,kr3,kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,pow,x2,x3,mu       
real :: time1,time2

real(8), dimension(NkBins,NkBins,NkBins) :: RadCont,RadContd,RadVar
real(8), dimension(:,:,:), allocatable :: deltaa,deltab,deltac,deltaad,deltabd,deltacd
integer, dimension(NkBins,NkBins,NkBins) :: RadCnt
integer :: OMP_GET_NUM_THREADS,nthr,TID,OMP_GET_THREAD_NUM
real(8) :: sumsqpower
real(8) :: cyfac
integer, parameter :: nmin=1



print*, 'start bispectrum'


call cpu_time(time1)
nc=ncell
ncr=ncell




RadCnt=0
RadCont=0.
RadContd=0.
RadVar=0.

sumsqpower=0.0d0


do w1=1,NkBins
print*,'w1',w1
allocate(deltaa(Nc+2,Nc,Nc))
deltaa=0.0d0
allocate(deltaad(Nc+2,Nc,Nc))
deltaad=0.0d0

do k1=nmin,nc
!print*,'k1',k1
!nthr=OMP_GET_NUM_THREADS()
!TID = OMP_GET_THREAD_NUM()
!print *,'Threads ',nthr,'ID ',tid


if (k1 .lt. ithresh) then
            kz1=k1-1
        else if (k1 .gt. nc-ithresh) then
	    kz1=k1-1-nc
        else
	    cycle
        endif
        do j1=nmin,nc
	    !print*,'j1',j1
            if (j1 .lt. ithresh) then
		ky1=j1-1
            else if (j1 .gt. nc-ithresh) then!j2 until june 6th
		ky1=j1-1-nc
	    else
		cycle
            endif
	    do i1=1,ncell+2,2
                kx1=(i1-1)/2
		kr1=sqrt(kx1**2.0+ky1**2.0+kz1**2.0)*2.0*pi/box
		if ((kr1>=kBinB(w1)) .and.(kr1<=kBinB(w1+1))) then
		    deltaa(i1,j1,k1)=delta(i1,j1,k1)
		    deltaa(i1+1,j1,k1)=delta(i1+1,j1,k1)
    		    deltaad(i1,j1,k1)=1.0d0
		    deltaad(i1+1,j1,k1)=0.0d0
		endif
		enddo
		enddo
		enddo
		call fft3d(deltaa,deltaa,'e')
		call fft3d(deltaad,deltaad,'e')
		
 
 
do w2=1,w1!NRadBins
print*,'w2',w2
    allocate(deltab(Nc+2,Nc,Nc))
    deltab=0.0d0
    allocate(deltabd(Nc+2,Nc,Nc))
    deltabd=0.0d0
    do k2=nmin,nc
        if (k2 .lt. ithresh) then
	    kz2=k2-1
        else if (k2 .gt. nc-ithresh) then
            kz2=k2-1-nc
	else
	    cycle
        endif
        do j2=nmin,nc
            if (j2 .lt. ithresh) then
                ky2=j2-1
            else if (j2 .gt. nc-ithresh) then
                ky2=j2-1-nc
	    else
		cycle
            endif
            do i2=1,ncell+2,2
		kx2=(i2-1)/2
		
		
		kr2=sqrt(kx2**2.0+ky2**2.0+kz2**2.0)*2.0*pi/box
		if ((kr2>=kBinB(w2)) .and.(kr2<=kBinB(w2+1))) then
		    deltab(i2,j2,k2)=delta(i2,j2,k2)
		    deltab(i2+1,j2,k2)=delta(i2+1,j2,k2)
    		    deltabd(i2,j2,k2)=1.0d0
		    deltabd(i2+1,j2,k2)=0.0d0
		endif
		enddo
		enddo
		enddo
		call fft3d(deltab,deltab,'e')
		call fft3d(deltabd,deltabd,'e')


do w3=1,w2!NRadBins
if (kbinc(w3)<=(kbinc(w1)+kbinc(w2)) .and. kbinc(w3)>=abs(kbinc(w1)-kbinc(w2))) then
print*,'w1,w2,w3',w1,w2,w3
allocate(deltac(Nc+2,Nc,Nc))
deltac=0.0d0
allocate(deltacd(Nc+2,Nc,Nc))
deltacd=0.0d0

do k3=nmin,nc
        if (k3 .lt. ithresh) then
	    kz3=k3-1
        else if (k3 .gt. nc-ithresh) then
            kz3=k3-1-nc
	else
	    cycle
        endif
        do j3=nmin,nc
            if (j3 .lt. ithresh) then
                ky3=j3-1
            else if (j3 .gt. nc-ithresh) then
                ky3=j3-1-nc
	    else
		cycle
            endif
            do i3=1,ncell+2,2
		kx3=(i3-1)/2
		
		kr3=sqrt(kx3**2.0+ky3**2.0+kz3**2.0)*2.0*pi/box
		if ((kr3>=kBinB(w3)) .and.(kr3<=kBinB(w3+1))) then
		    deltac(i3,j3,k3)=delta(i3,j3,k3)
		    deltac(i3+1,j3,k3)=delta(i3+1,j3,k3)
		    deltacd(i3,j3,k3)=1.0d0
		    deltacd(i3+1,j3,k3)=0.0d0
		endif
		enddo
		enddo
		enddo
		call fft3d(deltac,deltac,'e')
		call fft3d(deltacd,deltacd,'e')

print*,'Done ffts'
!$OMP PARALLEL DO SHARED(nc) PRIVATE(k1,j1,i1)  NUM_THREADS(16) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:RadCont,RadContD) 
	do j1=1,nc	
	    do k1=1,nc
		do i1=1,nc
		    RadCont(w1,w2,w3)=RadCont(w1,w2,w3)+deltaa(i1,j1,k1)*deltab(i1,j1,k1)*deltac(i1,j1,k1)
		    RadContD(w1,w2,w3)=RadContD(w1,w2,w3)+deltaad(i1,j1,k1)*deltabd(i1,j1,k1)*deltacd(i1,j1,k1)
		enddo
	    enddo
	enddo
!$OMP END PARALLEL DO

if (RadContD(w1,w2,w3)>0.0d0) then
    RadCont(w1,w2,w3)=RadCont(w1,w2,w3)/RadContD(w1,w2,w3)
else
    RadCont(w1,w2,w3)=0.0d0
endif


print*,'Done sum'

deallocate(deltac)
deallocate(deltacd)


endif
enddo
deallocate(deltab)
deallocate(deltabd)
enddo
deallocate(deltaa)
deallocate(deltaad)
enddo


call cpu_time(time2)
time2=(time2-time1)
write(*,'(a,f8.2)') ' Bispectrum calculated in', time2
return
end subroutine bispectrum_fft_kkk

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
		call dfftw_plan_with_nthreads(nthrfftw)
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

!///////////////////////////////////////////////////////////////////////////////
!			Export HDF5 Bispectrum
!///////////////////////////////////////////////////////////////////////////////
subroutine bispecthdf5export(filename,dsetname,bcnt,b,kavg)
use hdf5
implicit none
integer :: error  
CHARACTER(*) :: filename ! File name
CHARACTER(*) :: dsetname    ! Dataset name

INTEGER(HID_T) :: file_id	! File identifier
INTEGER(HID_T) :: dset_id	! Dataset identifier
INTEGER(HID_T) :: dspace_id	! Dataspace identifier
REAL(8), DIMENSION(nkbins,nkbins,nkbins) :: b
REAL(8), DIMENSION(nkbins,nkbins,nkbins,3) :: kavg
INTEGER, DIMENSION(nkbins,nkbins,nkbins) :: bcnt
INTEGER(HSIZE_T), DIMENSION(3), parameter :: dimsb=(/nkbins,nkbins,nkbins/)! Dataset dimensions
INTEGER     ::   rankb = 3
INTEGER(HSIZE_T), DIMENSION(4), parameter :: dimsbk=(/nkbins,nkbins,nkbins,3/)! Dataset dimensions
INTEGER     ::   rankbk = 4  
INTEGER(HSIZE_T), DIMENSION(1), parameter :: dimsk=(/nkbins/)! Dataset dimensions
INTEGER     ::   rankk = 1  



CALL h5open_f(error)
	CALL h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

	CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, b, dimsb, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)

	CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname)//'_cnt', H5T_NATIVE_INTEGER, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, bcnt, dimsb, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)

	CALL h5screate_simple_f(rankbk, dimsbk, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname)//'_kavg', H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kavg, dimsbk, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)

	CALL h5screate_simple_f(rankk, dimsk, dspace_id, error)
		CALL h5dcreate_f(file_id,'kbinc', H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kbinc, dimsk, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)



	CALL h5fclose_f(file_id, error)
CALL h5close_f(error)
end subroutine bispecthdf5export

!///////////////////////////////////////////////////////////////////////////////
!			Export HDF5 Bispectrum Mu
!///////////////////////////////////////////////////////////////////////////////
subroutine bispectmuhdf5export(filename,dsetname,bcnt,b,kavg)
use hdf5
implicit none
integer :: error  
CHARACTER(*) :: filename ! File name
CHARACTER(*) :: dsetname    ! Dataset name

INTEGER(HID_T) :: file_id	! File identifier
INTEGER(HID_T) :: dset_id	! Dataset identifier
INTEGER(HID_T) :: dspace_id	! Dataspace identifier
REAL(8), DIMENSION(nkbins,nkbins,nmubins) :: b
REAL(8), DIMENSION(nkbins,nkbins,nmubins,3) :: kavg
INTEGER, DIMENSION(nkbins,nkbins,nmubins) :: bcnt
INTEGER(HSIZE_T), DIMENSION(3), parameter :: dimsb=(/nkbins,nkbins,nmubins/)! Dataset dimensions
INTEGER     ::   rankb = 3
INTEGER(HSIZE_T), DIMENSION(4), parameter :: dimsbk=(/nkbins,nkbins,nmubins,3/)! Dataset dimensions
INTEGER     ::   rankbk = 4  
INTEGER(HSIZE_T), DIMENSION(1), parameter :: dimsk=(/nkbins/)! Dataset dimensions
INTEGER     ::   rankk = 1  



CALL h5open_f(error)
	CALL h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

	CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, b, dimsb, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)

	CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname)//'_cnt', H5T_NATIVE_INTEGER, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, bcnt, dimsb, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)

	CALL h5screate_simple_f(rankbk, dimsbk, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname)//'_kavg', H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kavg, dimsbk, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)

	CALL h5screate_simple_f(rankk, dimsk, dspace_id, error)
		CALL h5dcreate_f(file_id,'kbinc', H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kbinc, dimsk, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)



	CALL h5fclose_f(file_id, error)
CALL h5close_f(error)
end subroutine bispectmuhdf5export


!///////////////////////////////////////////////////////////////////////////////
!			Export HDF5 Field
!///////////////////////////////////////////////////////////////////////////////
subroutine fieldhdf5export(filename,dsetname,b)
use hdf5
implicit none
integer :: error  
integer, parameter :: ncellh=128
CHARACTER(*) :: filename ! File name
CHARACTER(*) :: dsetname    ! Dataset name

INTEGER(HID_T) :: file_id	! File identifier
INTEGER(HID_T) :: dset_id	! Dataset identifier
INTEGER(HID_T) :: dspace_id	! Dataspace identifier
REAL(8), DIMENSION(ncellh,ncellh,ncellh) :: b
INTEGER(HSIZE_T), DIMENSION(3), parameter :: dimsb=(/ncellh,ncellh,ncellh/)! Dataset dimensions
INTEGER     ::   rankb = 3



CALL h5open_f(error)
	CALL h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

	CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
		CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
			CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, b, dimsb, error)
		CALL h5dclose_f(dset_id, error)
	CALL h5sclose_f(dspace_id, error)



	CALL h5fclose_f(file_id, error)
CALL h5close_f(error)
end subroutine fieldhdf5export

!///////////////////////////////////////////////////////////////////////////////
!				HDF5 Open
!///////////////////////////////////////////////////////////////////////////////
subroutine hdf5openfile(filename,file_id)
use hdf5
implicit none
integer :: error  
CHARACTER(*) :: filename 	!File name
INTEGER(HID_T) :: file_id 	!File identifier


CALL h5open_f(error)
CALL h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

end subroutine hdf5openfile

!///////////////////////////////////////////////////////////////////////////////
!				HDF5 Close
!///////////////////////////////////////////////////////////////////////////////
subroutine hdf5closefile(file_id)
use hdf5
implicit none
integer :: error  
INTEGER(HID_T) :: file_id 	!File identifier


CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

end subroutine hdf5closefile

!///////////////////////////////////////////////////////////////////////////////
!				HDF5 Write Vec
!///////////////////////////////////////////////////////////////////////////////
subroutine hdf5writevec(file_id,dsetname,vec)
use hdf5
implicit none
integer :: error  
CHARACTER(*) :: dsetname    ! Dataset name

INTEGER(HID_T) :: file_id	! File identifier
INTEGER(HID_T) :: dset_id	! Dataset identifier
INTEGER(HID_T) :: dspace_id	! Dataspace identifier
REAL(8), DIMENSION(nkbins) :: vec
INTEGER(HSIZE_T), DIMENSION(1), parameter :: dimsb=(/nkbins/)! Dataset dimensions
INTEGER     ::   rankb = 1

CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
	CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
		CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vec, dimsb, error)
	CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)


end subroutine hdf5writevec

!///////////////////////////////////////////////////////////////////////////////
!				HDF5 Write Vec Cnt
!///////////////////////////////////////////////////////////////////////////////
subroutine hdf5writeveccnt(file_id,dsetname,vec,veccnt)
use hdf5
implicit none
integer :: error  
CHARACTER(*) :: dsetname    ! Dataset name

INTEGER(HID_T) :: file_id	! File identifier
INTEGER(HID_T) :: dset_id	! Dataset identifier
INTEGER(HID_T) :: dspace_id	! Dataspace identifier
REAL(8), DIMENSION(nkbins) :: vec,veccnt
INTEGER(HSIZE_T), DIMENSION(1), parameter :: dimsb=(/nkbins/)! Dataset dimensions
INTEGER     ::   rankb = 1

CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
	CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
		CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vec, dimsb, error)
	CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)

CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
	CALL h5dcreate_f(file_id, trim(dsetname)//'_cnt', H5T_NATIVE_DOUBLE, dspace_id,dset_id, error)
		CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vec, dimsb, error)
	CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)

end subroutine hdf5writeveccnt
!///////////////////////////////////////////////////////////////////////////////
!				HDF5 Write Vec
!///////////////////////////////////////////////////////////////////////////////
subroutine hdf5writevec_int(file_id,dsetname,vec)
use hdf5
implicit none
integer :: error  
CHARACTER(*) :: dsetname    ! Dataset name

INTEGER(HID_T) :: file_id	! File identifier
INTEGER(HID_T) :: dset_id	! Dataset identifier
INTEGER(HID_T) :: dspace_id	! Dataspace identifier
INTEGER(4), DIMENSION(nkbins) :: vec
INTEGER(HSIZE_T), DIMENSION(1), parameter :: dimsb=(/nkbins/)! Dataset dimensions
INTEGER     ::   rankb = 1

CALL h5screate_simple_f(rankb, dimsb, dspace_id, error)
	CALL h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id,dset_id, error)
		CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec, dimsb, error)
	CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)


end subroutine hdf5writevec_int
end module eval
