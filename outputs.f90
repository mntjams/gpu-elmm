module OUTPUTS
 use PARAMETERS
 use BOUNDARIES
 !use BLASIUS
 use SCALARS
 use WALLMODELS
 use SMAGORINSKY
 use GEOMETRIC

 implicit none

 real(KND),dimension(:),allocatable:: profuavg,profuavg2,profvavg,profvavg2,profuuavg,profvvavg,profwwavg,&
                                      profU,profV,profuu,profvv,profww,proftauavg,proftau,proftausgs,proftausgsavg,&
                                      proftemp,proftempfl,proftempavg,proftempavg2,proftempflavg,&
                                      proftempflsgs,proftempflsgsavg,proftt,profttavg,&
                                      profuw,profuwavg,profuwsgs,profuwsgsavg,&
                                      profvw,profvwavg,profvwsgs,profvwsgsavg
 real(KND),allocatable:: Uavg(:,:,:),Vavg(:,:,:),Wavg(:,:,:),Pravg(:,:,:)
 real(TIM),allocatable,dimension(:):: times
 real(KND),allocatable,dimension(:):: Utime,Vtime,Wtime,Prtime,temptime,CDtime,CLtime,deltime,tke,dissip,dissip2
 real(KND),allocatable,dimension(:,:):: scaltime !which scalar, time
 real(KND),allocatable,dimension(:,:,:):: scalptime !which scalar, position, time
 logical:: outputframePr=.false.
 logical:: outputframeU=.true.
 logical:: outputframevort=.false.
 logical:: outputframelambda2=.false.
 logical:: outputframescalars=.false.
 logical:: outputframesumscalars=.true.
 logical:: outputframeT=.true.
contains

 subroutine OUTPUT(U,V,W,Pr)
 real(KND),dimension(-2:,-2:,-2:):: U,V,W
 real(KND),dimension(1:,1:,1:):: Pr
 real(KND),dimension(:,:,:),allocatable:: Upat,Vpat,Wpat,depos
 real(KND),dimension(1:Unx,1:Uny,1:Unz)::Uinterp,Uinterpdir
 real(KND),dimension(1:Vnx,1:Vny,1:Vnz)::Vinterp,Vinterpdir
 real(KND),dimension(1:Wnx,1:Wny,1:Wnz)::Winterp,Winterpdir
 character(70):: str
 character(8)::  scalname="scalar00"
 integer i,j,k,l,m
 real(KND) S,S2,Suu,Svv,Sww,Suv
 type(TIBPoint),pointer:: IBP
 type(WMPOINT),pointer:: WMP

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  open(11,file="Uaxis.txt")
  do j=0,Unx
   write (11,*) j,U(j,max(Uny/2,1),max(Unz/2,1))
  enddo
  close(11)
  open(11,file="Praxis.txt")
  do j=1,Prnx
   write (11,*) j,Pr(j,max(Uny/2,1),max(Unz/2,1))
  enddo
  close(11)
  open(11,file="Prtime.txt")
  do j=0,endstep
   write (11,*) times(j),Prtime(j)
  enddo
  close(11)  
  open(11,file="Utime.txt")
  do j=0,endstep
   write (11,*) times(j),Utime(j)
  enddo
  close(11)
  open(11,file="Vtime.txt")
  do j=0,endstep
   write (11,*) times(j),Vtime(j)
  enddo
  close(11)
  open(11,file="Wtime.txt")
  do j=0,endstep
   write (11,*) times(j),Wtime(j)
  enddo
  close(11)
  open(11,file="temptime.txt")
  do j=0,endstep
   write (11,*) times(j),temptime(j)
  enddo
  close(11)
  open(11,file="deltime.txt")
  do j=1,endstep
   write (11,*) times(j),deltime(j)
  enddo
  close(11)
  open(11,file="tke.txt")
  do j=0,endstep
   write (11,*) times(j),tke(j)
  enddo
  close(11)
  open(11,file="dissip.txt")
  do j=1,endstep
   write (11,*) times(j),dissip(j)
  enddo
  close(11)
  open(11,file="dissip2.txt")
  do j=1,endstep
   write (11,*) times(j),dissip2(j)
  enddo
  close(11)
  
  if (computescalars==1) then
  open(11,file="scaltottime.txt")
  do j=1,endstep
   write (11,*) times(j),scaltime(1,j)
  enddo
  close(11)
  endif
  if (computescalars>1) then
  open(11,file="scaltottime.txt")
  do j=1,endstep
   write (11,*) times(j),scaltime(:,j)
  enddo
  close(11)
  endif

  if (computescalars==1.and.tasktype==1) then
  open(11,file="scal1time.txt")
  do j=1,endstep
   write (11,*) times(j),scalptime(1,1,j),scalptime(1,2,j),scalptime(1,3,j)
  enddo
  close(11)
  endif

  if (computescalars==2.and.tasktype==1) then
  open(11,file="scal1time.txt")
  do j=1,endstep
   write (11,*) times(j),scalptime(1,1,j),scalptime(1,2,j),scalptime(1,3,j)
  enddo
  close(11)
  open(11,file="scal2time.txt")
  do j=1,endstep
   write (11,*) times(j),scalptime(2,1,j),scalptime(2,2,j),scalptime(2,3,j)
  enddo
  close(11)
  endif

  if (computescalars>0.and.tasktype==8) then
  open(11,file="scal1time.txt")
  do j=1,endstep
   write (11,*) times(j),sum(scalptime(:,1,j)),sum(scalptime(:,2,j)),sum(scalptime(:,3,j))
  enddo
  close(11)
  endif

  if (tasktype==4) then
  open(11,file="cavmidU.txt")
  write (11,*) 0.,0.
  do j=1,Uny
   write (11,*) yPr(j),U(Unx/2,j,Unz/2)
  enddo
  write (11,*) 1.,1.
  close(11)
  open(11,file="cavmidV.txt")
  write (11,*) 0.,0.
  do i=1,Vnx
   write (11,*) xPr(i),V(i,Vny/2,Vnz/2)
  enddo
  write (11,*) 1.,0.
  close(11)
  open(11,file="cavmidW.txt")
  !write (11,*) 0.,0.
  do i=0,Wnz-1
   write (11,*) zPr(i),W(Wnx/2,Wny/2,i)
  enddo
  !write (11,*) 1.,0.
  close(11)
  endif
  
  if ((tasktype==4.or.tasktype==7).and.averaging==1) then
  profU=0
   open(11,file="profU.txt")
   do j=1,Uny
    write (11,*) yPr(j),yPr(j)*meanustar*Re , profUavg(j),profUavg(j)/meanustar
   enddo
   close(11)
   open(11,file="proftau.txt")
   do j=1,Prny
    write (11,*) yPr(j),yPr(j)*meanustar*Re ,proftauavg(j)/(meanustar*meanustar)
   enddo
   close(11)
   
   open(11,file="profuu.txt")
   do j=1,Uny/2
    write (11,*) yPr(j) ,profuuavg(j)!/(1._KND*Unx*Unz*meanustar**2)
   enddo
   close(11)
   open(11,file="profvv.txt")
   do j=1,Vny/2
    write (11,*) yV(j) ,profvvavg(j)!/(1._KND*Vnx*Vnz*meanustar**2)
   enddo
   close(11)
   open(11,file="profww.txt")
   do j=1,Wny/2
    write (11,*) yPr(j) ,profwwavg(j)!/(1._KND*Wnx*Wnz*meanustar**2)
   enddo
   close(11)

   open(11,file="profunres.txt")
   do j=1,Prny
    Suu=0
    Svv=0
    Sww=0
    do k=1,Prnz
     do i=1,Prnx
      Suu=Suu+(Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/(meanustar*meanustar)
      Svv=Svv+(Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k))/dyPr(j))/(meanustar*meanustar)
      Sww=Sww+(Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/(meanustar*meanustar)
     enddo
    enddo
    Suu=Suu/(Prnx*Prnz)
    Svv=Svv/(Prnx*Prnz)
    Sww=Sww/(Prnx*Prnz)
    write(11,*) yPr(j), Suu,Svv,Sww
   enddo
   close(11)
   open(11,file="profuntau.txt")
   do j=1,Vny
    Suv=0
    do k=1,Prnz
     do i=1,Prnx
      Suv=Suv+(((Visc(i,j,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i+1,j+1,k))/4._KND)*&
          ((U(i,j+1,k)-U(i,j,k))/dyV(j) +(V(i+1,j,k)-V(i,j,k))/dxU(i)))
     enddo
    enddo
    Suv=abs(Suv/(Prnx*Prnz*meanustar*meanustar))
    write(11,*) yV(j), Suv, Suv+(proftauavg(j)/(meanustar*meanustar)+proftauavg(j+1)/(meanustar*meanustar))/2._KND
   enddo
   close(11)
   if (computescalars>0.and.tasktype==7) then
    open(11,file="profscal1.txt")
    S2=0
    S=0
    do k=1,Prnz
     do i=1,Prnx
      S2=S2+SCALAR(i,1,k,1)
     enddo
    enddo
    S2=S2/(Prnx*Prnz)
    do j=1,Prny
     S=0
     do k=1,Prnz
      do i=1,Prnx
       S=S+SCALAR(i,j,k,1)
      enddo
     enddo
     S=S/(Prnx*Prnz)
     write (11,*) yPr(j),yPr(j)*meanustar*Re , abs(S-S2)
    enddo 
    close(11)
    if (computescalars==2) then
      open(11,file="profscal2.txt")
      S2=0
      S=0
      do k=1,Prnz
       do i=1,Prnx
        S2=S2+SCALAR(i,1,k,2)
       enddo
      enddo
      S2=S2/(Prnx*Prnz)
      do j=1,Prny
       S=0
       do k=1,Prnz
        do i=1,Prnx
         S=S+SCALAR(i,j,k,2)
        enddo
       enddo
       S=S/(Prnx*Prnz)
       write (11,*) yPr(j),yPr(j)*meanustar*Re , abs(S-S2)
      enddo 
      close(11)
    endif
   endif 
  endif
  
  if (tasktype==8) then
!    if (buoyancy>0) then
!     call CBLPROFILES(U,V,W,profU,profV,proftau,temperature,proftemp,proftempfl)
!    else
!     call CBLPROFILES(U,V,W,profU,profV,proftau)
!    endif   
   open(11,file="profu.txt")
   do k=1,Unz
    write (11,*) zPr(k),profuavg(k)
   enddo
   close(11)
   open(11,file="profv.txt")
   do k=1,Vnz
    write (11,*) zPr(k),profvavg(k)
   enddo
   close(11)
   open(11,file="profuu.txt")
   do k=1,Unz
    write (11,*) zPr(k),profuuavg(k)
   enddo
   close(11)
   open(11,file="profvv.txt")
   do k=1,Vnz
    write (11,*) zPr(k),profvvavg(k)
   enddo
   close(11)
   open(11,file="profww.txt")
   do k=1,Wnz
    write (11,*) zW(k),profwwavg(k)
   enddo
   close(11)
   open(11,file="profuw.txt")
   do k=0,Prnz
    write (11,*) zW(k),profuwavg(k),profuwsgsavg(k)
   enddo
   close(11)
   open(11,file="profvw.txt")
   do k=0,Prnz
    write (11,*) zW(k),profvwavg(k),profvwsgsavg(k)
   enddo
   close(11)
   if (buoyancy>0) then
    open(11,file="proftemp.txt")
    do k=1,Prnz
     write (11,*) zPr(k),proftempavg(k)
    enddo
    close(11)
    open(11,file="proftempfl.txt")
    do k=0,Prnz
     write (11,*) zW(k),proftempflavg(k),proftempflsgsavg(k)
    enddo
    close(11)
    open(11,file="proftt.txt")
    do k=1,Prnz
     write (11,*) zPr(k),profttavg(k)
    enddo
    close(11)
    open(11,file="profRig.txt")
    do k=1,Prnz
     S=0
     do j=1,Prny
      do i=1,Prnx
       S=S+Rig(i,j,k,Uavg,Vavg,temperatureavg)
      enddo
     enddo
     S=S/(Prnx*Prny)
     write (11,*) zPr(k),S
    enddo
    close(11)
    open(11,file="profRf.txt")
    do k=1,Prnz
     S=0
     S2=0
     do j=1,Prny
      do i=1,Prnx
       S=S+(Uavg(i,j,k+1)+Uavg(i-1,j,k+1)-Uavg(i,j,k-1)-Uavg(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
       S2=S2+(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
      enddo
     enddo
     S=S/(Prnx*Prny)
     S2=S2/(Prnx*Prny)
     if (abs((profuwavg(k)+profuwsgsavg(k))*S+(profvwavg(k)+profvwsgsavg(k))*S2)>tiny(1._KND)) then
      S=(grav_acc/temperature_ref)*(proftempflavg(k)+proftempflsgsavg(k))/&
       ((profuwavg(k)+profuwsgsavg(k))*S+(profvwavg(k)+profvwsgsavg(k))*S2)
     else
         S=0
     endif
     write (11,*) zPr(k),S
    enddo
    close(11)
   endif
  endif
 
  open(11,file="out.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Prnx,Prny,Prnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Prnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Prny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Prnz)
  str="POINT_DATA"
  write (str(12:),*) Prnx*Prny*Prnz
  write (11,"(A)") str

  
  write (11,"(A)") "SCALARS p float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) Pr(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)

  if (buoyancy>0) then
    write (11,*) "SCALARS temperature float"
   write (11,*) "LOOKUP_TABLE default"
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       write (11,*) Temperature(i,j,k)
     enddo
    enddo
   enddo 
   write (11,*)
  endif
  
!  write (11,*) "SCALARS ptype float"
!  write (11,*) "LOOKUP_TABLE default"
!  do k=1,Prnz
!   do j=1,Prny
!    do i=1,Prnx
!      Write (11,*) Prtype(i,j,k)
!    enddo
!   enddo
!  enddo 
!  write (11,*)
   
!  write (11,*) "SCALARS div float"
!  write (11,*) "LOOKUP_TABLE default"
!  do k=1,Prnz
!   do j=1,Prny
!    do i=1,Prnx
!      Write (11,*) (U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
!    enddo
!   enddo
!  enddo 
!  write (11,*)

  write (11,*) "SCALARS lambda2 float"
  write (11,*) "LOOKUP_TABLE default"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) Lambda2(i,j,k,U,V,W)
    enddo
   enddo
  enddo 
  write (11,*)
   
!   write (11,*) "SCALARS visc float"
!   write (11,*) "LOOKUP_TABLE default"
!   do k=1,Prnz
!    do j=1,Prny
!     do i=1,Prnx
!       Write (11,*) Visc(i,j,k)
!     enddo
!    enddo
!   enddo 
!   write (11,*)

  write (11,"(A)") "VECTORS u float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
    enddo
   enddo
  enddo

  write (11,"(A)") "VECTORS vort float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                      -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                   (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                      -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                   (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                      -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin)
    enddo
   enddo
  enddo
  write (11,*)

  if (tasktype==6.and.averaging==1) then
  write (11,"(A)") "VECTORS udef float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND-(Uavg(i,j,k)+Uavg(i-1,j,k))/2._KND,&
      (V(i,j,k)+V(i,j-1,k))/2._KND-(Vavg(i,j,k)+Vavg(i,j-1,k))/2._KND,&
      (W(i,j,k)+W(i,j,k-1))/2._KND-(Wavg(i,j,k)+Wavg(i,j,k-1))/2._KND
    enddo
   enddo
  enddo
  write (11,*) 
  endif
  
  if (tasktype==3) then
  write (11,"(A)") "VECTORS du float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (U(i+1,j,k)-U(i-1,j,k))/(2*dxmin),(U(i,j+1,k)-U(i,j-1,k))/(2*dymin)&
       ,(U(i,j,k+1)-U(i,j,k-1))/(2*dzmin)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "VECTORS dv float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (V(i+1,j,k)-V(i-1,j,k))/(2*dxmin),(V(i,j+1,k)-V(i,j-1,k))/(2*dymin)&
       ,(V(i,j,k+1)-V(i,j,k-1))/(2*dzmin)
    enddo
   enddo
  enddo
  write (11,*)
  endif 
  close(11)

 if (computescalars>0) then
  open(11,file="scalars.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Prnx,Prny,Prnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Prnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Prny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Prnz)
  str="POINT_DATA"
  write (str(12:),*) Prnx*Prny*Prnz
  write (11,"(A)") str

  do l=1,computescalars
  write(scalname(7:8),"(I2.2)") l  
  write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) SCALAR(i,j,k,l)
    enddo
   enddo
  enddo
  write (11,*)
  enddo
  close(11)

  if (computedeposition>0) then
   allocate(depos(1:Prnx,1:Prny,computescalars))
   depos=0
   if (associated(FirstWMPoint)) then
    WMP => FirstWMPoint
    do
     if (allocated(WMP%depscalar)) then
      do i=1,computescalars
       depos(WMP%x,WMP%y,i)=depos(WMP%x,WMP%y,i)+WMP%depscalar(i)
      enddo
     endif
     if (associated(WMP%next)) then
      WMP=>WMP%next
     else
      exit
     endif
    enddo
   endif

  open(11,file="deposition.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Prnx,Prny,1
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Prnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Prny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') 1,"float"
  write (11,"(A)") str
  write (11,*) zW(0)
  str="POINT_DATA"
  write (str(12:),*) Prnx*Prny
  write (11,"(A)") str

  do l=1,computescalars
  write(scalname(7:8),"(I2.2)") l  
  write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
  write (11,"(A)") "LOOKUP_TABLE default"

   do j=1,Prny
    do i=1,Prnx
      write (11,*) depos(i,j,l)
    enddo
   enddo

  write (11,*)
  enddo
  close(11)
  deallocate(depos)
   
  endif
 endif  



  if (averaging==1) then
  open(11,file="avg.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Prnx,Prny,Prnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Prnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Prny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Prnz)
  str="POINT_DATA"
  write (str(12:),*) Prnx*Prny*Prnz
  write (11,"(A)") str

  
  write (11,"(A)") "SCALARS p float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) Pravg(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)

  write (11,*) "SCALARS ptype float"
  write (11,*) "LOOKUP_TABLE default"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) Prtype(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  
  if (buoyancy>0) then
    write (11,*) "SCALARS temperature float"
   write (11,*) "LOOKUP_TABLE default"
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       write (11,*) Temperatureavg(i,j,k)
     enddo
    enddo
   enddo 
   write (11,*)
  endif
   
  write (11,"(A)") "VECTORS vort float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (Wavg(i,j+1,k)-Wavg(i,j-1,k)+Wavg(i,j+1,k-1)-Wavg(i,j-1,k-1))/(4*dxmin)&
                      -(Vavg(i,j,k+1)-Vavg(i,j,k-1)+Vavg(i,j-1,k+1)-Vavg(i,j-1,k-1))/(4*dymin),&
                   (Uavg(i,j,k+1)-Uavg(i,j,k-1)+Uavg(i-1,j,k+1)-Uavg(i-1,j,k-1))/(4*dxmin)&
                      -(Wavg(i+1,j,k)-Wavg(i-1,j,k)+Wavg(i+1,j,k-1)-Wavg(i-1,j,k-1))/(4*dymin),&
                   (Vavg(i+1,j,k)-Vavg(i-1,j,k)+Vavg(i+1,j-1,k)-Vavg(i-1,j-1,k))/(4*dxmin)&
                      -(Uavg(i,j+1,k)-Uavg(i,j-1,k)+Uavg(i-1,j+1,k)-Uavg(i-1,j-1,k))/(4*dymin)
    enddo
   enddo
  enddo
  write (11,*)

  write (11,"(A)") "VECTORS u float"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) (Uavg(i,j,k)+Uavg(i-1,j,k))/2._KND,(Vavg(i,j,k)+Vavg(i,j-1,k))/2._KND,(Wavg(i,j,k)+Wavg(i,j,k-1))/2._KND
    enddo
   enddo
  enddo
  close(11)

  if (computescalars>0) then
  open(11,file="scalarsavg.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Prnx,Prny,Prnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Prnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Prny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Prnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Prnz)
  str="POINT_DATA"
  write (str(12:),*) Prnx*Prny*Prnz
  write (11,"(A)") str
  
  do l=1,computescalars
  write(scalname(7:8),"(I2.2)") l  
  write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
      write (11,*) SCALARavg(i,j,k,l)
    enddo
   enddo
  enddo
  write (11,*)
  enddo
  close(11)
 endif  

  if (tasktype==6) then
  allocate(Upat(1:tilenx,1:tileny,1:tilenz*3))
  allocate(Vpat(1:tilenx,1:tileny,1:tilenz*3))
  allocate(Wpat(1:tilenx,1:tileny,1:tilenz*3))
  Upat=0
  Vpat=0
  Wpat=0
   do m=1,patternny
    do l=1,patternnx
     do k=1,tilenz*3
      do j=1,tileny
       do i=1,tilenx
        Upat(i,j,k)=Upat(i,j,k)+(Uavg(i+(l-1)*tilenx,j+(m-1)*tileny,k)+&
                                  Uavg(i+(l-1)*tilenx-1,j+(m-1)*tileny,k))/(2*patternnx*patternny)
        Vpat(i,j,k)=Vpat(i,j,k)+(Vavg(i+(l-1)*tilenx,j+(m-1)*tileny,k)+&
                                 Vavg(i+(l-1)*tilenx,j+(m-1)*tileny-1,k))/(2*patternnx*patternny)
        Wpat(i,j,k)=Wpat(i,j,k)+(Wavg(i+(l-1)*tilenx,j+(m-1)*tileny,k)+&
                                 Wavg(i+(l-1)*tilenx,j+(m-1)*tileny,k-1))/(2*patternnx*patternny)
       enddo
      enddo
     enddo  
    enddo
   enddo
  open(11,file="pattern.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) tilenx*2,tileny*2,3*tilenz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') tilenx*2,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:tilenx*2)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') tileny*2,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:tileny*2)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') 3*tilenz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:3*tilenz)
  str="POINT_DATA"
  write (str(12:),*) 2*tilenx*2*tileny*3*tilenz
  write (11,"(A)") str
  write (11,"(A)") "VECTORS u float"
  do k=1,3*tilenz
   do m=1,tileny*2
    do l=1,tilenx*2
      i=1+MOD(l-1,tilenx)
      j=1+MOD(m-1,tileny)
      write (11,*) (Upat(i,j,k)+Upat(i-1,j,k))/2._KND,(Vpat(i,j,k)+Vpat(i,j-1,k))/2._KND,(Wpat(i,j,k)+Wpat(i,j,k-1))/2._KND
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "VECTORS vort float"
  do k=1,3*tilenz
   do m=1,tileny*2
    do l=1,tilenx*2
      i=1+MOD(l-1,tilenx)
      j=1+MOD(m-1,tileny)
      write (11,*) (Wpat(i,j+1,k)-Wpat(i,j-1,k)+Wpat(i,j+1,k-1)-Wpat(i,j-1,k-1))/(4*dxmin)&
                      -(Vpat(i,j,k+1)-Vpat(i,j,k-1)+Vpat(i,j-1,k+1)-Vpat(i,j-1,k-1))/(4*dymin),&
                   (Upat(i,j,k+1)-Upat(i,j,k-1)+Upat(i-1,j,k+1)-Upat(i-1,j,k-1))/(4*dxmin)&
                      -(Wpat(i+1,j,k)-Wpat(i-1,j,k)+Wpat(i+1,j,k-1)-Wpat(i-1,j,k-1))/(4*dymin),&
                   (Vpat(i+1,j,k)-Vpat(i-1,j,k)+Vpat(i+1,j-1,k)-Vpat(i-1,j-1,k))/(4*dxmin)&
                      -(Upat(i,j+1,k)-Upat(i,j-1,k)+Upat(i-1,j+1,k)-Upat(i-1,j-1,k))/(4*dymin)
    enddo
   enddo
  enddo
  close(11)
  !profile
  open(11,file="profU.txt")
  do k=1,Unz
   S=0
   m=0
   do j=1,Uny
    do i=1,Unx
     if (prtype(9,j,k)==0) then
       S=S+Uavg(i,j,k)
       m=m+1
     endif  
    enddo
   enddo
   S=S/m
   write (11,*) zPr(k),S
  enddo
  close(11)
  open(11,file="proftau.txt")
  do k=1,Prnz
   S=0
   m=0
    do j=1,Prny
     do i=1,Prnx
     if (prtype(i,j,k)==0) then
       S=S+(U(i,j,k)+U(i-1,j,k)/2._KND)*(W(i,j,k)+W(i,j,k-1))/2._KND
       m=m+1
     endif  
    enddo
   enddo
   S=S/m
   write (11,*) zPr(k),S
  enddo
  close(11)
  open(11,file="profuu.txt")
  do k=1,Unz
   S=0
   m=0
    do j=1,Uny
     do i=1,Unx
     if (prtype(i,j,k)==0) then
       S=S+(U(i,j,k)-Uavg(i,j,k))*(U(i,j,k)-Uavg(i,j,k))
       m=m+1
     endif  
    enddo
   enddo
   S=S/m
   S=sqrt(S)
   write (11,*) zPr(k),S
  enddo
  close(11)
  open(11,file="profww.txt")
  do k=1,Wnz
   S=0
   m=0
    do j=1,Wny
     do i=1,Wnx
     if (prtype(i,j,k)==0) then
       S=S+(W(i,j,k)-Wavg(i,j,k))*(W(i,j,k)-Wavg(i,j,k))
       m=m+1
     endif  
    enddo
   enddo
   S=S/m
   S=sqrt(S)
   write (11,*) zW(k),S
  enddo
  close(11)
  
  endif 

endif
  Uinterp=-1
  Uinterpdir=-1
  Vinterp=-1
  Vinterpdir=-1
  Winterp=-1
  Winterpdir=-1
  if (associated(FirstIBPoint)) then
   IBP => FirstIBPoint
   do
    i=IBP%x
    j=IBP%y
    k=IBP%z
    if (IBP%component==1) then
     Uinterp(i,j,k)=IBP%interp
    elseif (IBP%component==2) then
     Vinterp(i,j,k)=IBP%interp
    elseif (IBP%component==3) then
     Winterp(i,j,k)=IBP%interp
    endif
    if (IBP%component==1) then
     Uinterpdir(i,j,k)=IBP%interpdir
    elseif (IBP%component==2) then
     Vinterpdir(i,j,k)=IBP%interpdir
    elseif (IBP%component==3) then
     Winterpdir(i,j,k)=IBP%interpdir
    endif
    if (associated(IBP%next)) then
     IBP=>IBP%next
    else
     exit
    endif
   enddo
  endif

  open(11,file="U.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Unx,Uny,Unz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Unx,"float"
  write (11,"(A)") str
  write (11,*) xU(1:Unx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Uny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Uny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Unz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Unz)
  str="POINT_DATA"
  write (str(12:),*) Unx*Uny*Unz
  write (11,"(A)") str

  
  write (11,"(A)") "SCALARS U float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
      write (11,*) U(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "SCALARS Uinterp float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
      write (11,*) Uinterp(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "SCALARS Uinterpdir float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
      write (11,*) Uinterpdir(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)

  
  open(11,file="V.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Vnx,Vny,Vnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Vnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vny,"float"
  write (11,"(A)") str
  write (11,*) yV(1:Vny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Vnz)
  str="POINT_DATA"
  write (str(12:),*) Vnx*Vny*Vnz
  write (11,"(A)") str
  

  write (11,"(A)") "SCALARS V float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
      write (11,*) V(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "SCALARS Vinterp float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
      write (11,*) Vinterp(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "SCALARS Vinterpdir float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
      write (11,*) Vinterpdir(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)


  
  open(11,file="W.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Wnx,Wny,Wnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Wnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Wny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wnz,"float"
  write (11,"(A)") str
  write (11,*) zW(1:Wnz)
  str="POINT_DATA"
  write (str(12:),*) Wnx*Wny*Wnz
  write (11,"(A)") str
  

  write (11,"(A)") "SCALARS W float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
      write (11,*) W(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "SCALARS Winterp float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
      write (11,*) Winterp(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  write (11,"(A)") "SCALARS Winterpdir float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
      write (11,*) Winterpdir(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)

!  if (tasktype==3) then
!   OPEN(11,file="px.vtk")
!   write (11,"(A)") "# vtk DataFile Version 2.0"
!   write (11,"(A)") "CLMM output file"
!   write (11,"(A)") "ASCII"
!   write (11,"(A)") "DATASET RECTILINEAR_GRID"
!   str="DIMENSIONS"
!   write (str(12:),*) Prnx,Prny,Prnz
!   write (11,"(A)") str
!   str="X_COORDINATES"
!   write (str(15:),'(i5,2x,a)') Prnx,"float"
!   write (11,"(A)") str
!   write (11,*) xPr(1:Prnx)
!   str="Y_COORDINATES"
!   write (str(15:),'(i5,2x,a)') Prny,"float"
!   write (11,"(A)") str
!   write (11,*) yPr(1:Prny)
!   str="Z_COORDINATES"
!   write (str(15:),'(i5,2x,a)') Prnz,"float"
!   write (11,"(A)") str
!   write (11,*) zPr(1:Prnz)
!   str="POINT_DATA"
!   write (str(12:),*) Prnx*Prny*Prnz
!   write (11,"(A)") str
! 
!   
!   write (11,"(A)") "SCALARS p float"
!   write (11,"(A)") "LOOKUP_TABLE default"
!   do k=1,Prnz
!    do j=1,Prny
!     do i=1,Prnx
!       Write (11,*) (Pr(i+1,j,k)-Pr(i,j,k))/dxmin
!     enddo
!    enddo
!   enddo
!   write (11,*)
!   CLOSE(11)
!  endif 
  
  write (*,*) "saved"
  end subroutine OUTPUT
  


  subroutine FRAME(U,V,W,Pr,n)
  real(KND):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer n,i,j,k,l
  character(13):: fname
  character(70):: str
  character(8)::  scalname="scalar00"
  character,parameter:: lf=char(10)
  integer mini,maxi,minj,maxj,mink,maxk


  vtkformat=binaryvtk
  if (vtkformat==textvtk) then

   if (framedimension==3) then
   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time
   
   open(11,file=fname)
   write (11,"(A)") "# vtk DataFile Version 2.0"
   write (11,"(A)") "CLMM output file"
   write (11,"(A)") "ASCII"
   write (11,"(A)") "DATASET RECTILINEAR_GRID"
   str="DIMENSIONS"
   write (str(12:),*) Prnx,Prny,Prnz
   write (11,"(A)") str
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prnx,"float"
   write (11,"(A)") str
   write (11,*) xPr(1:Prnx)
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prny,"float"
   write (11,"(A)") str
   write (11,*) yPr(1:Prny)
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prnz,"float"
   write (11,"(A)") str
   write (11,*) zPr(1:Prnz)
   str="POINT_DATA"
   write (str(12:),*) Prnx*Prny*Prnz
   write (11,"(A)") str 

   if (outputframePr) then
    write (11,"(A)") "SCALARS p float"
    write (11,"(A)") "LOOKUP_TABLE default"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) Pr(i,j,k)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif
  
   if (outputframelambda2) then
    write (11,*) "SCALARS lambda2 float"
    write (11,*) "LOOKUP_TABLE default"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) Lambda2(i,j,k,U,V,W)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo 
    write (11,*)
   endif
  
   if (outputframescalars) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (11,*) SCALAR(i,j,k,l)
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
    enddo
   elseif (outputframesumscalars.and.computescalars>0) then
     write (11,"(A,1X,A,1X,A)") "SCALARS", "scalar" , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (11,*) SUM(SCALAR(i,j,k,:))
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
   endif

   if (outputframeT) then
    if (buoyancy>0) then
      write (11,*) "SCALARS temperature float"
      write (11,*) "LOOKUP_TABLE default"
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         if (Prtype(i,j,k)==0) then
          write (11,*) Temperature(i,j,k)
         else
          write (11,*) 0.
         endif
        enddo
       enddo
      enddo 
      write (11,*)
    endif
   endif
   
   if (outputframeU) then
    write (11,"(A)") "VECTORS u float"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (outputframevort) then
    write (11,"(A)") "VECTORS vort float"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                       -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                     (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                     (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                        -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin)
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
    close(11)
   endif

  
  else
   if (slicedir==1) then
     call Gridcoords(mini,minj,mink,slicex,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
     maxi=mini
     minj=1
     maxj=Prny
     mink=1
     maxk=Prnz
    elseif (slicedir==2) then
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,slicex,(zW(Prnz+1)+zW(0))/2._KND)
     maxj=minj
     mini=1
     maxi=Prnx
     mink=1
     maxk=Prnz
    else
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,(yV(Prny+1)+yV(0))/2._KND,slicex)
     maxk=mink
     mini=1
     maxi=Prnx
     minj=1
     maxj=Prny
    endif

   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time
   
   open(11,file=fname)
   write (11,"(A)") "# vtk DataFile Version 2.0"
   write (11,"(A)") "CLMM output file"
   write (11,"(A)") "ASCII"
   write (11,"(A)") "DATASET RECTILINEAR_GRID"
   str="DIMENSIONS"
   write (str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
   write (11,"(A)") str
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxi-mini+1,"float"
   write (11,"(A)") str
   write (11,*) xPr(mini:maxi)
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
   write (11,"(A)") str
   write (11,*) yPr(minj:maxj)
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
   write (11,"(A)") str
   write (11,*) zPr(mink:maxk)
   str="POINT_DATA"
   write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
   write (11,"(A)") str

   if (outputframePr) then
    write (11,"(A)") "SCALARS p float"
    write (11,"(A)") "LOOKUP_TABLE default"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) Pr(i,j,k)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (outputframelambda2) then
    write (11,*) "SCALARS lambda2 float"
    write (11,*) "LOOKUP_TABLE default"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) Lambda2(i,j,k,U,V,W)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo 
    write (11,*)
   endif
  
   if (outputframescalars) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (11,*) SCALAR(i,j,k,l)
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
    enddo
   elseif (outputframesumscalars.and.computescalars>0) then
     write (11,"(A,1X,A,1X,A)") "SCALARS", "scalar" , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (11,*) SUM(SCALAR(i,j,k,:))
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
   endif

   if (outputframeT) then
    if (buoyancy>0) then
     write (11,*) "SCALARS temperature float"
     write (11,*) "LOOKUP_TABLE default"
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (11,*) Temperature(i,j,k)
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo 
     write (11,*)
    endif
   endif

   if (outputframeU) then
    write (11,"(A)") "VECTORS u float"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (outputframevort) then
    write (11,"(A)") "VECTORS vort float"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                        -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                     (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                     (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                       -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin)
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
   endif  
   close(11)
   endif


  else  !Binary VTK format



   if (framedimension==3) then
   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time
   
   open(20,file=fname,access='stream',status='replace',form="unformatted",action="write")
   write (20) "# vtk DataFile Version 2.0",lf
   write (20) "CLMM output file",lf
   write (20) "BINARY",lf
   write (20) "DATASET RECTILINEAR_GRID",lf
   str="DIMENSIONS"
   write (str(12:),*) Prnx,Prny,Prnz
   write (20) str,lf
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prnx,"float"
   write (20) str,lf
   write (20) (real(xPr(i),SNG),i=1,Prnx),lf
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prny,"float"
   write (20) str,lf
   write (20) (real(yPr(j),SNG),j=1,Prny),lf
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prnz,"float"
   write (20) str,lf
   write (20) (real(zPr(k),SNG),k=1,Prnz),lf
   str="POINT_DATA"
   write (str(12:),*) Prnx*Prny*Prnz
   write (20) str,lf

   if (outputframePr) then
    write (20) "SCALARS p float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) real(Pr(1:Prnx,1:Prny,1:Prnz),SNG)
       else
        write (20) 0._SNG
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif
  
   if (outputframelambda2) then
    write (20) "SCALARS lambda2 float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) real(Lambda2(i,j,k,U,V,W),SNG)
       else
        write (20) 0._SNG
       endif
      enddo
     enddo
    enddo 
    write (20) lf
   endif
  
   if (outputframescalars) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (20) "SCALARS ", scalname , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (20) real(SCALAR(i,j,k,l),SNG)
        else
         write (20) 0._SNG
        endif
       enddo
      enddo
     enddo
     write (20) lf
    enddo
   elseif (outputframesumscalars.and.computescalars>0) then
     write (20) "SCALARS ", "scalar" , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (20) real(SUM(SCALAR(i,j,k,:)),SNG)
        else
         write (20) 0._SNG
        endif
       enddo
      enddo
     enddo
     write (20) lf
   endif

   if (outputframeT) then
    if (buoyancy>0) then
      write (20) "SCALARS temperature float",lf
      write (20) "LOOKUP_TABLE default",lf
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         if (Prtype(i,j,k)==0) then
          write (20) real(Temperature(i,j,k),SNG)
         else
          write (20) 0._SNG
         endif
        enddo
       enddo
      enddo 
      write (20) lf
    endif
   endif
   
   if (outputframeU) then
    write (20) "VECTORS u float",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) real((U(i,j,k)+U(i-1,j,k))/2._KND,SNG),real((V(i,j,k)+V(i,j-1,k))/2._KND,SNG)&
         ,real((W(i,j,k)+W(i,j,k-1))/2._KND,SNG)
       else
        write (20) 0._SNG,0._SNG,0._SNG
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (outputframevort) then
    write (20) "VECTORS vort float",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                       -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),SNG),&
                     real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),SNG),&
                     real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                        -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin),SNG)
       else
        write (20) 0._SNG,0._SNG,0._SNG
       endif
      enddo
     enddo
    enddo
    close(20)
   endif

  
  else
   if (slicedir==1) then
     call Gridcoords(mini,minj,mink,slicex,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
     maxi=mini
     minj=1
     maxj=Prny
     mink=1
     maxk=Prnz
    elseif (slicedir==2) then
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,slicex,(zW(Prnz+1)+zW(0))/2._KND)
     maxj=minj
     mini=1
     maxi=Prnx
     mink=1
     maxk=Prnz
    else
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,(yV(Prny+1)+yV(0))/2._KND,slicex)
     maxk=mink
     mini=1
     maxi=Prnx
     minj=1
     maxj=Prny
    endif

   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time
   
   open(20,file=fname,access='stream',status='replace',form="unformatted",action="write")
   write (20) "# vtk DataFile Version 2.0",lf
   write (20) "CLMM output file",lf
   write (20) "BINARY",lf
   write (20) "DATASET RECTILINEAR_GRID",lf
   str="DIMENSIONS"
   write (str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
   write (20) str,lf
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxi-mini+1,"float"
   write (20) str,lf
   write (20) (real(xPr(i),SNG),i=mini,maxi),lf
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
   write (20) str,lf
   write (20) (real(yPr(j),SNG),j=minj,maxj),lf
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
   write (20) str,lf
   write (20) (real(zPr(k),SNG),k=mink,maxk),lf
   str="POINT_DATA"
   write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
   write (20) str,lf

   if (outputframePr) then
    write (20) "SCALARS p float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) real(Pr(i,j,k),SNG)
       else
        write (20) 0._SNG
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (outputframelambda2) then
    write (20) "SCALARS lambda2 float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) real(Lambda2(i,j,k,U,V,W),SNG)
       else
        write (20) 0._SNG
       endif
      enddo
     enddo
    enddo 
    write (20) lf
   endif
  
   if (outputframescalars) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (20) "SCALARS ", scalname , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (20) real(SCALAR(i,j,k,l),SNG)
        else
         write (20) 0._SNG
        endif
       enddo
      enddo
     enddo
     write (20) lf
    enddo
   elseif (outputframesumscalars.and.computescalars>0) then
     write (20) "SCALARS ", "scalar" , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (20) real(SUM(SCALAR(i,j,k,:)),SNG)
        else
         write (20) 0._SNG
        endif
       enddo
      enddo
     enddo
     write (20) lf
   endif

   if (outputframeT) then
    if (buoyancy>0) then
     write (20) "SCALARS temperature float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (20) real(Temperature(i,j,k),SNG)
        else
         write (20) 0._SNG
        endif
       enddo
      enddo
     enddo 
     write (20) lf
    endif
   endif

   if (outputframeU) then
    write (20) "VECTORS u float",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) real((U(i,j,k)+U(i-1,j,k))/2._KND,SNG),real((V(i,j,k)+V(i,j-1,k))/2._KND,SNG)&
         ,real((W(i,j,k)+W(i,j,k-1))/2._KND,SNG)
       else
        write (20) 0._SNG,0._SNG,0._SNG
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (outputframevort) then
    write (20) "VECTORS vort float",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                        -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),SNG),&
                     real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),SNG),&
                     real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                       -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin),SNG)
       else
        write (20) 0._SNG,0._SNG,0._SNG
       endif
      enddo
     enddo
    enddo
   endif  
   close(20)
   endif

  endif
  endsubroutine FRAME

!  subroutine BLASPROF(U,V,W)
!  real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
!  real(KND) XI
!  integer j

!  write(*,*) "Saving U profile.."
!  write (*,*) "x=",xU(Unx-(Unx-iXup)/2)
  
!  OPEN(11,file="Uprof")
!  write(11,*) 0,0,0 
!  do j=0,Uny
!   XI=XIxy(xU(Unx-(Unx-iXup)/2),yPr(j))
!   write(11,*) XI,Ublas(xU(Unx-(Unx-iXup)/2),yPr(j)),U(Unx-(Unx-iXup)/2,j,Unz/2)
!  enddo
!  CLOSE(11)
!  write(*,*) "Saving V profile.."
!  write (*,*) "x=",xPr(Unx-(Unx-iXup)/2)

!  OPEN(11,file="Vprof")
!  write(11,*) 0,0,0 
!  do j=0,Vny
!   XI=XIxy(xPr(Unx-(Unx-iXup)/2),yV(j))
!   write(11,*) XI,Vblas(xPr(Unx-(Unx-iXup)/2),yV(j)),V(Unx-(Unx-iXup)/2,j,Vnz/2)
!  enddo
!  CLOSE(11)
!  do j=0,Vny
!  enddo

  
!  write (*,*) "Profiles saved"
!  endsubroutine



  subroutine PROFILES(U,V,W,profu,proftau,profuu, profvv, profww)
  real(KND),dimension(:):: profu,proftau,profuu,profvv,profww
  real(KND),dimension(-2:,-2:,-2:):: U,V,W
  real(KND):: S, S2
  integer i,j,k
  
   do j=1,Uny
    do k=1,Unz
     do i=1,Unx
      profU(j)=profU(j)+U(i,j,k)
     enddo
    enddo
    profU(j)=profU(j)/(Unx*Unz)
   enddo

   do j=1,Prny
    S=0
    S2=0
    do k=1,Prnz
     do i=1,Prnx
      S=S+((U(i,j,k)+U(i-1,j,k)/2._KND)-profUavg2(j))*(V(i,j,k)+V(i,j,k-1))/2._KND
      S2=S2+((W(i,j,k)+W(i-1,j,k)/2._KND)-profUavg2(j))*(V(i,j,k)+V(i,j-1,k))/2._KND
     enddo
    enddo
    proftau(j)=SQRT(S**2+s2**2)/(Prnz*Prnx)
   enddo

   do j=1,Uny
    S=0
    S2=0
    do k=1,Unz
     do i=1,Unx
      S=S+(U(i,j,k)-profUavg2(j))**2
     enddo
    enddo
    profuu(j)=S/(1._KND*Unx*Unz*meanustar**2)
   enddo

   do j=1,Vny
    S=0
    S2=0
    do k=1,Vnz
     do i=1,Vnx
      S=S+(V(i,j,k)*V(i,j,k))
     enddo
    enddo
    profvv(j)=S/(1._KND*Vnx*Vnz*meanustar**2)
   enddo

   do j=1,Wny
    S=0
    S2=0
    do k=1,Wnz
     do i=1,Wnx
      S=S+(W(i,j,k)*W(i,j,k))
     enddo
    enddo
    profww(j)=S/(1._KND*Wnx*Wnz*meanustar**2)
   enddo



  endsubroutine PROFILES


  subroutine CBLPROFILES(U,V,W,temperature)
  real(KND),dimension(-2:,-2:,-2:):: U,V,W
  real(KND),dimension(-1:,-1:,-1:),optional:: temperature
  real(KND):: S, S2
  real(KND) Str(1:3,1:3)
  real(KND),allocatable,save::fp(:),ht(:),gp(:)
  integer i,j,k,n
  integer,save:: called=0

   do k=0,Unz+1
    S=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)==0) then
       S=S+U(i,j,k)
       n=n+1
      endif
     enddo
    enddo
    profU(k)=S/n
   enddo
  
   do k=1,Vnz+1
    S=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)==0) then
       S=S+V(i,j,k)
       n=n+1
      endif
     enddo
    enddo
    profV(k)=S/n
   enddo
  
  if (present(temperature)) then
   do k=1,Prnz
    S=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k)==0) then
       S=S+temperature(i,j,k)
       n=n+1
      endif
     enddo
    enddo
    profTemp(k)=S/n
   enddo
  endif
   
   if (called==0) then
    allocate(fp(0:Prnx+1),ht(0:Prnz+1),gp(0:Prny+1))
    forall (i=0:Prnx+1)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
    forall (k=0:Prnz+1)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
    forall (j=0:Prny+1)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
   endif

   do k=0,Prnz
    S=0
    S2=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if ((Utype(i,j,k+1)==0.or.Utype(i,j,k)==0).and.(Wtype(i+1,j,k)==0.or.Wtype(i,j,k)==0)) then
       S=S+((ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))-((1-ht(k))*profU(k)+ht(k)*profU(k+1)))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))
       n=n+1
      endif
     enddo
    enddo
    profuw(k)=S/n
   enddo
   
   do k=0,Prnz
    S=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if ((Vtype(i,j,k+1)==0.or.Vtype(i,j,k)==0).and.(Wtype(i,j+1,k)==0.or.Wtype(i,j,k)==0)) then
       S=S+(ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k)-(ht(k)*profV(k+1)+(1-ht(k))*profV(k)))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))
       n=n+1
      endif
     enddo
    enddo
    profvw(k)=S/n
   enddo
   
   do k=0,Prnz
    S=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k+1)==0.or.Utype(i,j,k)==0) then
       S=S-0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)
       n=n+1
      endif
     enddo
    enddo
    profuwsgs(k)=S/n
   enddo
   
   do k=0,Prnz
    S=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k+1)==0.or.Vtype(i,j,k)==0) then
       S=S-0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)
       n=n+1
      endif
     enddo
    enddo
    profvwsgs(k)=S/n
   enddo
      
   do k=1,Unz
    S=0
    S2=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)==0) then
       S=S+(U(i,j,k)-profU(k))**2
       n=n+1
      endif
     enddo
    enddo
    profuu(k)=S/n
   enddo

   do k=1,Vnz
    S=0
    S2=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)==0) then
       S=S+(V(i,j,k)-profV(k))**2
       n=n+1
      endif
     enddo
    enddo
    profvv(k)=S/n
   enddo

   do k=0,Wnz
    S=0
    S2=0
    n=0
    do j=1,Wny
     do i=1,Wnx
      if (Wtype(i,j,k)==0) then
       S=S+(W(i,j,k))**2
       n=n+1
      endif
     enddo
    enddo
    profww(k)=S/n
   enddo

  if (present(temperature)) then
   do k=0,Prnz
    S=0
    S2=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
       S=S+0.5_KND*(temperature(i,j,k+1)+temperature(i,j,k))*(W(i,j,k))
       n=n+1
      endif
     enddo
    enddo
    proftempfl(k)=S/n
   enddo

   do k=1,Prnz
    S=0
    S2=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k)==0) then
       S=S+(temperature(i,j,k)-profTemp(k))**2
       n=n+1
      endif
     enddo
    enddo
    proftt(k)=S/n
   enddo

   do k=0,Prnz
    S=0
    S2=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
        S=S-(0.5_KND*(TDiff(i,j,k+1)+TDiff(i,j,k))*(temperature(i,j,k+1)-temperature(i,j,k)))/dzW(k)
        n=n+1
      endif
     enddo
    enddo
    proftempflsgs(k)=S/n
   enddo
  endif
  called=1
  endsubroutine CBLPROFILES

  
  real(KND) function TotKE(U,V,W)
  real(KND),dimension(-2:,-2:,-2:):: U,V,W
  real(KND) Um,Vm,Wm
  integer i,j,k
   TotKE=0
   Um=0
   Vm=0
   Wm=0
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      TotKE=TotKE+(((U(i-1,j,k)+U(i,j,k))/2._KND-Um)**2+&
                    ((V(i,j-1,k)+V(i,j,k))/2._KND-Vm)**2+&
                    ((W(i,j,k-1)+W(i,j,k))/2._KND-Wm)**2)
     enddo
    enddo
   enddo
   TotKE=TotKE*lx*lz*lz/2
   endfunction TotKE

   real(KND) function Vorticity(i,j,k,U,V,W)
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:):: U,V,W

    Vorticity=((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dymin)&
                      -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dzmin))**2+&
              ((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dzmin)&
                      -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dxmin))**2+&
              ((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                      -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin))**2
    Vorticity=Sqrt(Vorticity)
   endfunction Vorticity
   
 
   real(KND) function Lambda2(i,j,k,U,V,W)
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:):: U,V,W

   Lambda2=((U(i,j,k)-U(i-1,j,k))/dxmin)**2
   Lambda2=Lambda2+((V(i,j,k)-V(i,j-1,k))/dymin)**2
   Lambda2=Lambda2+((W(i,j,k)-W(i,j,k-1))/dzmin)**2
   Lambda2=Lambda2+((V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4*dxmin))**2
   Lambda2=Lambda2+((W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4*dxmin))**2
   Lambda2=Lambda2+((U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4*dymin))**2
   Lambda2=Lambda2+((W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4*dymin))**2
   Lambda2=Lambda2+((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4*dzmin))**2
   Lambda2=Lambda2+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4*dzmin))**2
   Lambda2=-Sqrt(Lambda2)
   Lambda2=Vorticity(i,j,k,U,V,W)+Lambda2
   

   endfunction Lambda2

   real(KND) function MOMTHICK(U)
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:):: U
   real(KND) p,S
    S=0
    do j=1,Uny
     p=0
     do k=1,Unz
      do i=1,Unx
       p=p+U(i,j,k)
      enddo
     enddo
     p=p/(Unx*Uny)
     S=S+(Uinlet-p)*(p+Uinlet)*dyPr(j)
    enddo
    MOMTHICK=S/(4._KND)
   endfunction MOMTHICK




  subroutine OUTPUTU2(U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: U,V,W
 integer i,j,k
 character(70):: str

  open(11,file="U2.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Unx,Uny,Unz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Unx,"float"
  write (11,"(A)") str
  write (11,*) xU(1:Unx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Uny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Uny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Unz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Unz)
  str="POINT_DATA"
  write (str(12:),*) Unx*Uny*Unz
  write (11,"(A)") str

  
  write (11,"(A)") "SCALARS U float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
      write (11,*) U(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)

  
  open(11,file="V2.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Vnx,Vny,Vnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Vnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vny,"float"
  write (11,"(A)") str
  write (11,*) yV(1:Vny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Vnz)
  str="POINT_DATA"
  write (str(12:),*) Vnx*Vny*Vnz
  write (11,"(A)") str
  

  write (11,"(A)") "SCALARS V float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
      write (11,*) V(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)


  
  open(11,file="W2.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Wnx,Wny,Wnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Wnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Wny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wnz,"float"
  write (11,"(A)") str
  write (11,*) zW(1:Wnz)
  str="POINT_DATA"
  write (str(12:),*) Wnx*Wny*Wnz
  write (11,"(A)") str
  

  write (11,"(A)") "SCALARS W float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
      write (11,*) W(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)
  endsubroutine OUTPUTU2


  subroutine OUTINLET(U,V,W)
  !for output of 2d data for use as an inilet condition later
  real(KND),dimension(-2:,-2:,-2:):: U,V,W
  integer i,j,k
  integer,save::fnum
  integer,save:: called=0

   if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
       .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
    if (called==0) then
     open(101,file="inletframeinfo.unf",form='unformatted',status='replace',action='write')
     write(101) Prny,Prnz  !for check of consistency of grids before use
     write(101) Vny
     write(101) Wnz
     write(101) dxPr(0)
     called=1
     fnum=0
    endif
    fnum=fnum+1
    write(101) time-timefram1
    call OUTINLETFRAME(U,V,W,fnum)
   elseif (time>timefram2+(timefram2-timefram1)/(frames-1).and.called==1) then
     close(101)
     called=2
   endif

  endsubroutine OUTINLET


  subroutine OUTINLETFRAME(U,V,W,n)
  real(KND):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer n,i,j,k,l
  character(12):: fname
  integer mini,maxi,minj,maxj,mink,maxk



  call Gridcoords(mini,minj,mink,slicex,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
  maxi=mini
  minj=1
  maxj=Prny
  mink=1
  maxk=Prnz


  fname(1:5)="frame"
  write(fname(6:8),"(I3.3)") n
  fname(9:12)=".unf"
  write(*,*) "Saving frame:",fname(1:6),"   time:",time
  
  open(11,file=fname,form='unformatted',access='sequential',status='replace',action='write')
  

  write(11) U(mini,1:Uny,1:Unz)
  write(11) V(mini,1:Vny,1:Vnz)
  write(11) W(mini,1:Wny,1:Wnz)
  if (buoyancy>0) then
       write(11) Temperature(mini,1:Prny,1:Prnz)
  endif
  close(11)

  endsubroutine OUTINLETFRAME



   pure real(KND) function TriLinInt(a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111)
   real(KND),intent(in):: a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111

    TriLinInt=   (1-a)*(1-b)*(1-c)*vel000+&
                 a*(1-b)*(1-c)*vel100+&
                 (1-a)*b*(1-c)*vel010+&
                 (1-a)*(1-b)*c*vel001+&
                 a*b*(1-c)*vel110+&
                 a*(1-b)*c*vel101+&
                 (1-a)*b*c*vel011+&
                 a*b*c*vel111

   endfunction TriLinInt


   subroutine GridCoordsU(xi,yj,zk,x,y,z)
   integer,intent(out):: xi,yj,zk
   real(KND),intent(in):: x,y,z
   integer i

   xi=Unx+1
   do i=1,Unx+1
    if (xPr(i+1)>=x) then
                  xi=i-1
                  exit
                 endif
   enddo

   yj=Vny
   do i=1,Vny
    if (yV(i)>=y) then
                  yj=i
                  exit
                 endif
   enddo
   zk=Wnz
   do i=1,Wnz
    if (zW(i)>=z) then
                  zk=i
                  exit
                 endif
   enddo
   endsubroutine GridCoordsU


   subroutine GridCoordsV(xi,yj,zk,x,y,z)
   integer,intent(out):: xi,yj,zk
   real(KND),intent(in):: x,y,z
   integer i

   xi=Unx
   do i=1,Unx
    if (xU(i)>=x) then
                  xi=i
                  exit
                 endif
   enddo

   yj=Vny+1
   do i=1,Vny+1
    if (yPr(i)>=y) then
                  yj=i-1
                  exit
                 endif
   enddo
   zk=Wnz
   do i=1,Wnz
    if (zW(i)>=z) then
                  zk=i
                  exit
                 endif
   enddo
   endsubroutine GridCoordsV


   subroutine GridCoordsW(xi,yj,zk,x,y,z)
   integer,intent(out):: xi,yj,zk
   real(KND),intent(in):: x,y,z
   integer i

   xi=Unx
   do i=1,Unx
    if (xU(i)>=x) then
                  xi=i
                  exit
                 endif
   enddo

   yj=Vny
   do i=1,Vny
    if (yV(i)>=y) then
                  yj=i
                  exit
                 endif
   enddo
   zk=Wnz+1
   do i=1,Wnz+1
    if (zPr(i)>=z) then
                  zk=i-1
                  exit
                 endif
   enddo
   endsubroutine GridCoordsW


   subroutine GridCoordsPr(xi,yj,zk,x,y,z)
   integer,intent(out):: xi,yj,zk
   real(KND),intent(in):: x,y,z
   integer i

   xi=Prnx
   do i=1,Prnx
    if (xU(i)>=x) then
                  xi=i
                  exit
                 endif
   enddo

   yj=Prny
   do i=1,Prny
    if (yV(i)>=y) then
                  yj=i
                  exit
                 endif
   enddo

   zk=Prnz
   do i=1,Prnz
    if (zW(i)>=z) then
                  zk=i
                  exit
                 endif
   enddo
   endsubroutine GridCoordsPr


end module OUTPUTS
