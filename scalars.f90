module SCALARS
use PARAMETERS
use WALLMODELS
use GEOMETRIC

implicit none
 
  private minmod,minmod2,minmod3,superbee,limiter,heav,GAMMAL,vanalbada,vanleer

 interface MINMOD
         module procedure MINMOD2, MINMOD3
  end interface

  real(KND),allocatable,dimension(:,:,:,:):: SCALAR,SCALARavg  !last index is number of scalar (because of paging)
  real(KND),dimension(:),allocatable:: partdiam,partrho,percdistrib !diameter of particles <=0 for gas
  

contains

  subroutine ADVSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  
   call PLMSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)

  endsubroutine ADVSCALAR

  pure subroutine CDSSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  integer nx,ny,nz,i,j,k
  real(KND) Ax,Ay,Az
        nx=Prnx
        ny=Prny
        nz=Prnz

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif

        Ax=0.5_KND*coef*dt/dxmin
        Ay=0.5_KND*coef*dt/dymin
        Az=0.5_KND*coef*dt/dzmin


        do k=1,nz
            do j=1,ny
                do i=1,nx
                    SCAL2(i,j,k)=SCAL2(i,j,k)- ((Az*(SCAL(i,j,k+1)+SCAL(i,j,k))*(W(i,j,k))&
                    -Az*(SCAL(i,j,k)+SCAL(i,j,k-1))*(W(i,j,k-1)))&
                    +(Ay*(SCAL(i,j+1,k)+SCAL(i,j,k))*(V(i,j,k))&
                    -Ay*(SCAL(i,j,k)+SCAL(i,j-1,k))*(V(i,j-1,k)))&
                    +(Ax*(SCAL(i+1,j,k)+SCAL(i,j,k))*(U(i,j,k))&
                    -Ax*(SCAL(i,j,k)+SCAL(i-1,j,k))*(U(i-1,j,k))))
                enddo
            enddo
        enddo
  end subroutine CDSSCALAR

  pure subroutine UDSSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  integer nx,ny,nz,i,j,k
  real(KND) Ax,Ay,Az
        nx=Prnx
        ny=Prny
        nz=Prnz

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif
        
        Ax=coef*dt/dxmin
        Ay=coef*dt/dymin
        Az=coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if (U(i,j,k)>0) then
                     SCAL2(i,j,k)=SCAL2(i,j,k)-Ax*SCAL(i,j,k)*(U(i,j,k))
                    else
                     SCAL2(i,j,k)=SCAL2(i,j,k)-Ax*SCAL(i+1,j,k)*(U(i,j,k))
                    endif
                    if (U(i-1,j,k)>0) then
                     SCAL2(i,j,k)=SCAL2(i,j,k)+Ax*SCAL(i-1,j,k)*(U(i-1,j,k))
                    else
                     SCAL2(i,j,k)=SCAL2(i,j,k)+Ax*SCAL(i,j,k)*(U(i-1,j,k))
                    endif

                     if (V(i,j,k)>0) then
                      SCAL2(i,j,k)=SCAL2(i,j,k)-Ay*SCAL(i,j,k)*(V(i,j,k))
                     else
                      SCAL2(i,j,k)=SCAL2(i,j,k)-Ay*SCAL(i,j+1,k)*(V(i,j,k))
                     endif
                     if (V(i,j-1,k)>0) then
                     SCAL2(i,j,k)=SCAL2(i,j,k)+Ay*SCAL(i,j-1,k)*(V(i,j-1,k))
                    else
                     SCAL2(i,j,k)=SCAL2(i,j,k)+Ay*SCAL(i,j,k)*(V(i,j-1,k))
                    endif


                    if (W(i,j,k)>0) then
                     SCAL2(i,j,k)=SCAL2(i,j,k)-Az*SCAL(i,j,k)*(W(i,j,k))
                    else
                     SCAL2(i,j,k)=SCAL2(i,j,k)-Az*SCAL(i,j,k+1)*(W(i,j,k))
                    endif
                    if (W(i,j,k-1)>0) then
                     SCAL2(i,j,k)=SCAL2(i,j,k)+Az*SCAL(i,j,k-1)*(W(i,j,k-1))
                    else
                     SCAL2(i,j,k)=SCAL2(i,j,k)+Az*SCAL(i,j,k)*(W(i,j,k-1))
                    endif

                enddo
            enddo
        enddo
  end subroutine UDSSCALAR



  pure subroutine QUICKSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  real(KND)::Scal3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer nx,ny,nz,i,j,k
  real(KND) Ax,Ay,Az,F

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif
        
   Ax=coef*dt/dxmin
   Ay=coef*dt/dymin
   Az=coef*dt/dzmin
   SCAL2=SCAL
   SCAL3=SCAL
   do k=1,Prnz
    do j=1,Prny
     do i=0,Prnx
       if (U(i,j,k)>0) then
        F=Ax*U(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k)+(3._KND/8._KND)*SCAL3(i-1,j,k)-(1._KND/8._KND)*SCAL3(i-1,j,k))
       else
        F=Ax*U(i,j,k)*((6._KND/8._KND)*SCAL3(i+1,j,k)+(3._KND/8._KND)*SCAL3(i,j,k)-(1._KND/8._KND)*SCAL3(i+2,j,k))
       endif
       SCAL2(i,j,k)=SCAL2(i,j,k)-F
       SCAL2(i+1,j,k)=SCAL2(i+1,j,k)+F
     enddo
    enddo
   enddo
       
   if (sctype==1) then
    call BOUND_Temp(SCAL2)
   else 
    call BOUND_PASSSCALAR(SCAL2)
   endif
   SCAL3=SCAL2
   do k=1,Prnz
    do j=0,Prny
     do i=1,Prnx
       if (V(i,j,k)>0) then
        F=Ay*V(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k)+(3._KND/8._KND)*SCAL3(i,j-1,k)-(1._KND/8._KND)*SCAL3(i,j-1,k))
       else
        F=Ay*V(i,j,k)*((6._KND/8._KND)*SCAL3(i,j+1,k)+(3._KND/8._KND)*SCAL3(i,j,k)-(1._KND/8._KND)*SCAL3(i,j+1,k))
       endif
       SCAL2(i,j,k)=SCAL2(i,j,k)-F
       SCAL2(i,j+1,k)=SCAL2(i,j+1,k)+F
     enddo
    enddo
   enddo

   if (sctype==1) then
    call BOUND_Temp(SCAL2)
   else 
    call BOUND_PASSSCALAR(SCAL2)
   endif
   SCAL3=SCAL2
   do k=0,Prnz
    do j=1,Prny
     do i=1,Prnx
       if (W(i,j,k)>0) then
        F=Az*W(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k)+(3._KND/8._KND)*SCAL3(i,j,k-1)-(1._KND/8._KND)*SCAL3(i,j,k-1))
       else
        F=Az*W(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k+1)+(3._KND/8._KND)*SCAL3(i,j,k)-(1._KND/8._KND)*SCAL3(i,j,k+1))
       endif
       SCAL2(i,j,k)=SCAL2(i,j,k)-F
       SCAL2(i,j,k+1)=SCAL2(i,j,k+1)+F
     enddo
    enddo
   enddo
   SCAL2=SCAL2-SCAL
  end subroutine QUICKSCALAR





  subroutine PLMSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  integer i,j,k
  real(KND) SL,SR,FLUX,A
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2):: SLOPE
  logical,save::direction=.true.

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif
  
  A=coef*dt
  SCAL2=SCAL2+SCAL

  if (direction) then
   SLOPE=0
   do i=0,Prnx+1
    do j=0,Prny+1
     do k=0,Prnz+1
      SR=(SCAL2(i+1,j,k)-SCAL2(i,j,k))/dxU(i)
      SL=(SCAL2(i,j,k)-SCAL2(i-1,j,k))/dxU(i-1)
      SLOPE(i,j,k)=LIMITER(SL,SR)*dxPr(i)
     enddo
    enddo
   enddo


   do i=0,Prnx
    do j=1,Prny
     do k=1,Prnz
      if (U(i,j,k)>0) then
       FLUX=U(i,j,k)*SCAL2(i,j,k)+U(i,j,k)*(1-U(i,j,k)*A/dxPr(i))*SLOPE(i,j,k)/2._KND
      else 
       FLUX=U(i,j,k)*SCAL2(i+1,j,k)-U(i,j,k)*(1+U(i,j,k)*A/dxPr(i+1))*SLOPE(i+1,j,k)/2._KND
      endif
      SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dxPr(i)
      SCAL2(i+1,j,k)=SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
     enddo
    enddo
   enddo    

   SLOPE=0
   do i=0,Prnx+1
    do j=0,Prny+1
     do k=0,Prnz+1
      SR=(SCAL2(i,j+1,k)-SCAL2(i,j,k))/dyV(j)
      SL=(SCAL2(i,j,k)-SCAL2(i,j-1,k))/dyV(j-1)
      SLOPE(i,j,k)=LIMITER(SL,SR)*dyPr(j)
     enddo
    enddo
   enddo


   do i=1,Prnx
    do j=0,Prny
     do k=1,Prnz
      if (V(i,j,k)>0) then
       FLUX=V(i,j,k)*SCAL2(i,j,k)+V(i,j,k)*(1-V(i,j,k)*A/dyPr(j))*SLOPE(i,j,k)/2._KND
      else 
       FLUX=V(i,j,k)*SCAL2(i,j+1,k)-V(i,j,k)*(1+V(i,j,k)*A/dyPr(j+1))*SLOPE(i,j+1,k)/2._KND
      endif
      SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dyPr(j)
      SCAL2(i,j+1,k)=SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
     enddo
    enddo
   enddo    

   SLOPE=0
   do i=0,Prnx+1
    do j=0,Prny+1
     do k=0,Prnz+1
      SR=(SCAL2(i,j,k+1)-SCAL2(i,j,k))/dzW(k)
      SL=(SCAL2(i,j,k)-SCAL2(i,j,k-1))/dzW(k-1)
      SLOPE(i,j,k)=LIMITER(SL,SR)*dzPr(k)
     enddo
    enddo
   enddo

   do i=1,Prnx
    do j=1,Prny
     do k=0,Prnz
      if (W(i,j,k)>0) then
       FLUX=W(i,j,k)*SCAL2(i,j,k)+W(i,j,k)*(1-W(i,j,k)*A/dzPr(k))*SLOPE(i,j,k)/2._KND
      else 
       FLUX=W(i,j,k)*SCAL2(i,j,k+1)-W(i,j,k)*(1+W(i,j,k)*A/dzPr(k+1))*SLOPE(i,j,k+1)/2._KND
      endif
      SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dzPr(k)
      SCAL2(i,j,k+1)=SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
     enddo
    enddo
   enddo
  else
   SLOPE=0
   do i=0,Prnx+1
    do j=0,Prny+1
     do k=0,Prnz+1
      SR=(SCAL2(i,j,k+1)-SCAL2(i,j,k))/dzW(k)
      SL=(SCAL2(i,j,k)-SCAL2(i,j,k-1))/dzW(k-1)
      SLOPE(i,j,k)=LIMITER(SL,SR)*dzPr(k)
     enddo
    enddo
   enddo

   do i=1,Prnx
    do j=1,Prny
     do k=0,Prnz
      if (W(i,j,k)>0) then
       FLUX=W(i,j,k)*SCAL2(i,j,k)+W(i,j,k)*(1-W(i,j,k)*A/dzPr(k))*SLOPE(i,j,k)/2._KND
      else 
       FLUX=W(i,j,k)*SCAL2(i,j,k+1)-W(i,j,k)*(1+W(i,j,k)*A/dzPr(k+1))*SLOPE(i,j,k+1)/2._KND
      endif
      SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dzPr(k)
      SCAL2(i,j,k+1)=SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
     enddo
    enddo
   enddo

   SLOPE=0
   do i=0,Prnx+1
    do j=0,Prny+1
     do k=0,Prnz+1
      SR=(SCAL2(i,j+1,k)-SCAL2(i,j,k))/dyV(j)
      SL=(SCAL2(i,j,k)-SCAL2(i,j-1,k))/dyV(j-1)
      SLOPE(i,j,k)=LIMITER(SL,SR)*dyPr(j)
     enddo
    enddo
   enddo


   do i=1,Prnx
    do j=0,Prny
     do k=1,Prnz
      if (V(i,j,k)>0) then
       FLUX=V(i,j,k)*SCAL2(i,j,k)+V(i,j,k)*(1-V(i,j,k)*A/dyPr(j))*SLOPE(i,j,k)/2._KND
      else 
       FLUX=V(i,j,k)*SCAL2(i,j+1,k)-V(i,j,k)*(1+V(i,j,k)*A/dyPr(j+1))*SLOPE(i,j+1,k)/2._KND
      endif
      SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dyPr(j)
      SCAL2(i,j+1,k)=SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
     enddo
    enddo
   enddo    

   SLOPE=0
   do i=0,Prnx+1
    do j=0,Prny+1
     do k=0,Prnz+1
      SR=(SCAL2(i+1,j,k)-SCAL2(i,j,k))/dxU(i)
      SL=(SCAL2(i,j,k)-SCAL2(i-1,j,k))/dxU(i-1)
      SLOPE(i,j,k)=LIMITER(SL,SR)*dxPr(i)
     enddo
    enddo
   enddo


   do i=0,Prnx
    do j=1,Prny
     do k=1,Prnz
      if (U(i,j,k)>0) then
       FLUX=U(i,j,k)*SCAL2(i,j,k)+U(i,j,k)*(1-U(i,j,k)*A/dxPr(i))*SLOPE(i,j,k)/2._KND
      else 
       FLUX=U(i,j,k)*SCAL2(i+1,j,k)-U(i,j,k)*(1+U(i,j,k)*A/dxPr(i+1))*SLOPE(i+1,j,k)/2._KND
      endif
      SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dxPr(i)
      SCAL2(i+1,j,k)=SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
     enddo
    enddo
   enddo    
  endif

  SCAL2=SCAL2-SCAL
  endsubroutine PLMSCALAR



  pure subroutine PLMSCALARNOSPLIT(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  integer i,j,k
  real(KND) SL,SR,FLUX,A
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2):: SLOPE

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif
  
  A=coef*dt

  SLOPE=0
  do i=0,Prnx+1
   do j=0,Prny+1
    do k=0,Prnz+1
     SR=(SCAL(i+1,j,k)-SCAL(i,j,k))/dxU(i)
     SL=(SCAL(i,j,k)-SCAL(i-1,j,k))/dxU(i-1)
     SLOPE(i,j,k)=LIMITER(SL,SR)*dxPr(i)
    enddo
   enddo
  enddo


  do i=0,Prnx
   do j=1,Prny
    do k=1,Prnz
     if (U(i,j,k)>0) then
      FLUX=U(i,j,k)*SCAL(i,j,k)+U(i,j,k)*(1-U(i,j,k)*A/dxPr(i))*SLOPE(i,j,k)/2._KND
     else 
      FLUX=U(i,j,k)*SCAL(i+1,j,k)-U(i,j,k)*(1+U(i,j,k)*A/dxPr(i+1))*SLOPE(i+1,j,k)/2._KND
     endif
     SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dxPr(i)
     SCAL2(i+1,j,k)=SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo    

  SLOPE=0
  do i=0,Prnx+1
   do j=0,Prny+1
    do k=0,Prnz+1
     SR=(SCAL(i,j+1,k)-SCAL(i,j,k))/dyV(j)
     SL=(SCAL(i,j,k)-SCAL(i,j-1,k))/dyV(j-1)
     SLOPE(i,j,k)=LIMITER(SL,SR)*dyPr(j)
    enddo
   enddo
  enddo


  do i=1,Prnx
   do j=0,Prny
    do k=1,Prnz
     if (V(i,j,k)>0) then
      FLUX=V(i,j,k)*SCAL(i,j,k)+V(i,j,k)*(1-V(i,j,k)*A/dyPr(j))*SLOPE(i,j,k)/2._KND
     else 
      FLUX=V(i,j,k)*SCAL(i,j+1,k)-V(i,j,k)*(1+V(i,j,k)*A/dyPr(j+1))*SLOPE(i,j+1,k)/2._KND
     endif
     SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dyPr(j)
     SCAL2(i,j+1,k)=SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo    

  SLOPE=0
  do i=0,Prnx+1
   do j=0,Prny+1
    do k=0,Prnz+1
     SR=(SCAL(i,j,k+1)-SCAL(i,j,k))/dzW(k)
     SL=(SCAL(i,j,k)-SCAL(i,j,k-1))/dzW(k-1)
     SLOPE(i,j,k)=LIMITER(SL,SR)*dzPr(k)
    enddo
   enddo
  enddo

  do i=1,Prnx
   do j=1,Prny
    do k=0,Prnz
     if (W(i,j,k)>0) then
      FLUX=W(i,j,k)*SCAL(i,j,k)+W(i,j,k)*(1-W(i,j,k)*A/dzPr(k))*SLOPE(i,j,k)/2._KND
     else 
      FLUX=W(i,j,k)*SCAL(i,j,k+1)-W(i,j,k)*(1+W(i,j,k)*A/dzPr(k+1))*SLOPE(i,j,k+1)/2._KND
     endif
     SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dzPr(k)
     SCAL2(i,j,k+1)=SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  endsubroutine PLMSCALARNOSPLIT




  subroutine KAPPASCALAR(SCAL2,SCAL,U,V,W,sctype,coef) !Kappa scheme with flux limiter
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in):: sctype
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2):: SLOPE
  real(KND),parameter::eps=1e-8

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif
  
  A=coef*dt
  SLOPE=0
  do i=0,Prnx
   do j=1,Prny
    do k=1,Prnz
     if (U(i,j,k)>0) then
      SR=(SCAL(i+1,j,k)-SCAL(i,j,k))!/dxU(i)
      SL=(SCAL(i,j,k)-SCAL(i-1,j,k))!/dxU(i-1)
     else
      SR=(SCAL(i,j,k)-SCAL(i+1,j,k))!/dxU(i)
      SL=(SCAL(i+1,j,k)-SCAL(i+2,j,k))!/dxU(i-1)
     endif
!      if (abs(SR+eps*sign(1._KND,SL))>1E5*abs(SL+eps*sign(1._KND,SL)))&
!        write(*,*) SR,SL,SR+eps*sign(1._KND,SL),SL+eps*sign(1._KND,SL)
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo


  do i=0,Prnx
   do j=1,Prny
    do k=1,Prnz
     if (U(i,j,k)>0) then
      FLUX=U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._KND)
     endif
     SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dxPr(i)
     SCAL2(i+1,j,k)=SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo    


  SLOPE=0
  do i=1,Prnx
   do j=0,Prny
    do k=1,Prnz
     if (V(i,j,k)>0) then
      SR=(SCAL(i,j+1,k)-SCAL(i,j,k))!/dxU(i)
      SL=(SCAL(i,j,k)-SCAL(i,j-1,k))!/dxU(i-1)
     else
      SR=(SCAL(i,j,k)-SCAL(i,j+1,k))!/dxU(i)
      SL=(SCAL(i,j+1,k)-SCAL(i,j+2,k))!/dxU(i-1)
     endif
!      if (abs(SR+eps*sign(1._KND,SL))>1E5*abs(SL+eps*sign(1._KND,SL)))&
!          write(*,*) SR,SL,SR+eps*sign(1._KND,SL),SL+eps*sign(1._KND,SL)
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo


  do i=1,Prnx
   do j=0,Prny
    do k=1,Prnz
     if (V(i,j,k)>0) then
      FLUX=V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._KND)
     endif

     SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dyPr(j)
     SCAL2(i,j+1,k)=SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo    


  SLOPE=0
  do i=1,Prnx
   do j=1,Prny
    do k=0,Prnz
    if (W(i,j,k)>0) then
      SR=(SCAL(i,j,k+1)-SCAL(i,j,k))!/dxU(i)
      SL=(SCAL(i,j,k)-SCAL(i,j,k-1))!/dxU(i-1)
     else
      SR=(SCAL(i,j,k)-SCAL(i,j,k+1))!/dxU(i)
      SL=(SCAL(i,j,k+1)-SCAL(i,j,k+2))!/dxU(i-1)
     endif
!      if (abs(SR+eps*sign(1._KND,SL))>1E5*abs(SL+eps*sign(1._KND,SL)))&
!        write(*,*) SR,SL,SR+eps*sign(1._KND,SL),SL+eps*sign(1._KND,SL)
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo


  do i=1,Prnx
   do j=1,Prny
    do k=0,Prnz
     if (W(i,j,k)>0) then
      FLUX=W(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j,k-1))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=W(i,j,k)*(SCAL(i,j,k+1)+(SCAL(i,j,k+1)-SCAL(i,j,k+2))*SLOPE(i,j,k)/2._KND)
     endif

     SCAL2(i,j,k)=SCAL2(i,j,k)-A*FLUX/dzPr(k)
     SCAL2(i,j,k+1)=SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  endsubroutine KAPPASCALAR




  subroutine DIFFSCALAR(SCAL2,SCAL,sctype,coef)
  real(KND),intent(inout)::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in):: coef
  integer,intent(in):: sctype
  integer nx,ny,nz,i,j,k,l
  real(KND) p,S,Ap(-1:Prnx+1,-1:Prny+1,-1:Prnz+1),Scal3(-1:Prnx+1,-1:Prny+1,-1:Prnz+1)
  type(TScalFlIBPoint),pointer::SFlIBP

   nx=Prnx
   ny=Prny
   nz=Prnz


  if (Re>0) then
   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif

   do k=1,Prnz  !initital value using forward Euler
    do j=1,Prny
     do i=1,Prnx
      SCAL3(i,j,k)=(((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL(i+1,j,k)-SCAL(i,j,k))/dxU(i)-&
        (TDiff(i,j,k)+TDiff(i-1,j,k))*(SCAL(i,j,k)-SCAL(i-1,j,k))/dxU(i-1))/(dxPr(i))+&
       ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL(i,j+1,k)-SCAL(i,j,k))/dyV(j)-&
        (TDiff(i,j,k)+TDiff(i,j-1,k))*(SCAL(i,j,k)-SCAL(i,j-1,k))/dyV(j-1))/(dyPr(j))+&
       ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL(i,j,k+1)-SCAL(i,j,k))/dzW(k)-&
        (TDiff(i,j,k)+TDiff(i,j,k-1))*(SCAL(i,j,k)-SCAL(i,j,k-1))/dzW(k-1))/(dzPr(k)))
     enddo
    enddo
   enddo

   Ap(1,1,1)=dt*coef
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       SCAL2(i,j,k)=SCAL(i,j,k)+Ap(1,1,1)*SCAL3(i,j,k)
     enddo
    enddo
   enddo

   call ScalFlSOURC(SCAL2,sctype)
   if (associated(FirstScalFlIBPoint)) then
    SFlIBP => FirstScalFlIBPoint
    do
     i=SFlIBP%x
     j=SFlIBP%y
     k=SFlIBP%z
     SCAL2(i,j,k)=SCAL2(i,j,k)+SFlIBP%ScalSourc*dt
     SCAL3(i,j,k)=SCAL3(i,j,k)+SFlIBP%ScalSourc*dt
     if (associated(SFlIBP%next)) then
      SFlIBP=>SFlIBP%next
     else
      exit
     endif
    enddo
   endif

    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       Ap(i,j,k)=1._KND/(1._KND/(coef*dt)+(((TDiff(i+1,j,k)+TDiff(i,j,k))/dxU(i)+&
                            (TDiff(i,j,k)+TDiff(i-1,j,k))/dxU(i-1))/(4._KND*dxPr(i))+&
                            ((TDiff(i,j+1,k)+TDiff(i,j,k))/dyV(j)+&
                            (TDiff(i,j,k)+TDiff(i,j-1,k))/dyV(j-1))/(4._KND*dyPr(j))+&
                            ((TDiff(i,j,k+1)+TDiff(i,j,k))/dzW(k)+&
                            (TDiff(i,j,k)+TDiff(i,j,k-1))/dzW(k-1))/(4._KND*dzPr(k))))
      enddo
     enddo
    enddo

   do l=1,maxCNiter
   S=0
   if (sctype==1) then
    call BOUND_Temp(SCAL2)
   else 
    call BOUND_PASSSCALAR(SCAL2)
   endif
    !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:S)
    !$OMP DO
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (REDBLACKPr(i,j,k)) then
        p=(SCAL(i,j,k)/(coef*dt))+(SCAL3(i,j,k)+&
         ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k))/dxU(i)-&
          (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k))/dxU(i-1))/(dxPr(i))+&
         ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k))/dyV(j)-&
          (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k))/dyV(j-1))/(dyPr(j))+&
         ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1))/dzW(k)-&
          (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1))/dzW(k-1))/(dzPr(k))&
         )/4._KND
         p=p*Ap(i,j,k)
         S=max(S,abs(p-SCAL2(i,j,k)))
         SCAL2(i,j,k)=p
       endif
      enddo
     enddo
    enddo
    !$OMP ENDDO
    !$OMP DO
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (.not.REDBLACKPr(i,j,k)) then
        p=(SCAL(i,j,k)/(coef*dt))+(SCAL3(i,j,k)+&
         ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k))/dxU(i)-&
          (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k))/dxU(i-1))/(dxPr(i))+&
         ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k))/dyV(j)-&
          (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k))/dyV(j-1))/(dyPr(j))+&
         ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1))/dzW(k)-&
          (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1))/dzW(k-1))/(dzPr(k))&
         )/4._KND
         p=p*Ap(i,j,k)
         S=max(S,abs(p-SCAL2(i,j,k)))
         SCAL2(i,j,k)=p
       endif
      enddo
     enddo
    enddo
    !$OMP ENDDO
   !$OMP ENDPARALLEL
      write (*,*) "CNscalar ",l,S
     if (S<=epsCN) exit
   enddo
  else
   SCAL2=SCAL
   call ScalFlSOURC(SCAL2,sctype)
   if (associated(FirstScalFlIBPoint)) then
    SFlIBP => FirstScalFlIBPoint
    do
     i=SFlIBP%x
     j=SFlIBP%y
     k=SFlIBP%z
     SCAL2(i,j,k)=SCAL2(i,j,k)+SFlIBP%ScalSourc*dt

     if (associated(SFlIBP%next)) then
      SFlIBP=>SFlIBP%next
     else
      exit
     endif
    enddo
   endif
  endif
  endsubroutine DIFFSCALAR


  pure subroutine BOUND_PASSSCALAR(SCAL)
  real(KND),intent(inout):: SCAL(-1:,-1:,-1:)
  integer i,j,k,nx,ny,nz
   nx=Prnx
   ny=Prny
   nz=Prnz


   if (tasktype==7) then
     do k=1,nz
      do j=1,ny                      !Periodic BC
       SCAL(0,j,k)=SCAL(nx,j,k)
      enddo
     enddo
   elseif (BtypeW==NOSLIP) then
    do k=1,nz
     do j=1,ny                      !No flux on noslip
      SCAL(0,j,k)=SCAL(1,j,k)
      SCAL(-1,j,k)=SCAL(1,j,k)
     enddo
    enddo
   else
    do k=1,nz
     do j=1,ny                      !Other BCs
      SCAL(0,j,k)=0
      SCAL(-1,j,k)=0
     enddo
    enddo
   endif

   if (tasktype==7) then
     do k=1,nz
      do j=1,ny                      !Periodic BC
       SCAL(nx+1,j,k)=SCAL(1,j,k)
      enddo
     enddo
   elseif (BtypeE==NOSLIP) then
    do k=1,nz
     do j=1,ny                      !No flux on noslip
      SCAL(nx+1,j,k)=SCAL(nx,j,k)
      SCAL(nx+2,j,k)=SCAL(nx,j,k)
     enddo
    enddo
   else
    do k=1,nz
     do j=1,ny                      !Other BCs
      SCAL(nx+1,j,k)=SCAL(nx,j,k)
      SCAL(nx+2,j,k)=SCAL(nx,j,k)
     enddo
    enddo
   endif

   if (BtypeS==PERIODIC) then
    do k=1,nz
     do i=-1,nx+2                      !Periodic BC
      SCAL(i,0,k)=SCAL(i,ny,k)
      SCAL(i,-1,k)=SCAL(i,ny-1,k)
     enddo
    enddo
   elseif (BtypeS==NOSLIP) then
    do k=1,nz
     do i=-1,nx+2                      !No flux on noslip
      SCAL(i,0,k)=SCAL(i,1,k)
      SCAL(i,-1,k)=SCAL(i,1,k)
     enddo
    enddo
   else
    do k=1,nz
     do i=-1,nx+2                      !Other BCs
      SCAL(i,0,k)=SCAL(i,1,k)
      SCAL(i,-1,k)=SCAL(i,1,k)
     enddo
    enddo
   endif 

   if (BtypeN==PERIODIC) then
    do k=1,nz
     do i=-1,nx+2                      !Periodic BC
      SCAL(i,ny+1,k)=SCAL(i,1,k)
      SCAL(i,ny+2,k)=SCAL(i,2,k)
     enddo
    enddo
   elseif (BtypeN==NOSLIP) then
    do k=1,nz
     do i=-1,nx+2                      !No flux on noslip
      SCAL(i,ny+1,k)=SCAL(i,ny,k)
      SCAL(i,ny+2,k)=SCAL(i,ny,k)
     enddo
    enddo
   else
    do k=1,nz
     do i=-1,nx+2                      !Other BCs
      SCAL(i,ny+1,k)=SCAL(i,ny,k)
      SCAL(i,ny+2,k)=SCAL(i,ny,k)
     enddo
    enddo
   endif 

   if (BtypeB==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2                      !Periodic BC
      SCAL(i,j,0)=SCAL(i,j,nz)
      SCAL(i,j,-1)=SCAL(i,j,nz-1)
     enddo
    enddo
   elseif (BtypeB==NOSLIP) then
    do j=-1,ny+2
     do i=-1,nx+2                      !No flux on noslip
      SCAL(i,j,0)=SCAL(i,j,1)
      SCAL(i,j,-1)=SCAL(i,j,1)
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2                      !Other BCs
      SCAL(i,j,0)=SCAL(i,j,1)
      SCAL(i,j,-1)=SCAL(i,j,1)
     enddo
    enddo
   endif 

   if (BtypeT==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2                      !Periodic BC
      SCAL(i,j,nz+1)=SCAL(i,j,1)
      SCAL(i,j,nz+2)=SCAL(i,j,2)
     enddo
    enddo
   elseif (BtypeT==NOSLIP) then
    do j=-1,ny+2
     do i=-1,nx+2                      !No flux on noslip
      SCAL(i,j,nz+1)=SCAL(i,j,nz)
      SCAL(i,j,nz+2)=SCAL(i,j,nz)
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2                      !Other BCs
      SCAL(i,j,nz+1)=SCAL(i,j,nz)
      SCAL(i,j,nz+2)=SCAL(i,j,nz)
     enddo
    enddo
   endif 
  end subroutine BOUND_PASSSCALAR


  



  pure subroutine BOUND_TEMP(Theta)
  real(KND),intent(inout):: theta(-1:,-1:,-1:)
  integer i,j,k,nx,ny,nz
   nx=Prnx
   ny=Prny
   nz=Prnz
!    if (TBtypeW==DIRICHLET) then
!      do k=1,nz
!       do j=1,ny                     
!        theta(0,j,k)=WsideTemp!-(theta(1,j,k)-WsideTemp)
!        theta(-1,j,k)=WsideTemp!-(theta(2,j,k)-WsideTemp)
!       enddo
!      enddo
   if (TBtypeW==DIRICHLET) then
     do k=1,nz
      do j=1,ny                     
       theta(0,j,k)=Tempin(j,k)!-(theta(1,j,k)-WsideTemp)
       theta(-1,j,k)=Tempin(j,k)!-(theta(2,j,k)-WsideTemp)
      enddo
     enddo
   else if (TBtypeW==PERIODIC) then
    do k=1,nz
     do j=1,ny
      theta(0,j,k)=theta(nx,j,k)
      theta(-1,j,k)=theta(nx-1,j,k)
     enddo
    enddo
   else if (TBtypeW==NEUMANN) then
    do k=1,nz
     do j=1,ny                     
      theta(0,j,k)=theta(1,j,k)-WsideTemp*dxU(0)
      theta(-1,j,k)=theta(1,j,k)-WsideTemp*(dxU(0)+dxU(-1))
     enddo
    enddo
   else if (TBtypeW==CONSTFLUX) then
    do k=1,nz
     do j=1,ny                      
!       theta(0,j,k)=(theta(1,j,k)*(U(0,j,k)/2._KND-(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+WsideTemp)/&
!                                   (U(0,j,k)/2._KND+(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
      theta(0,j,k)=(theta(1,j,k)*((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+WsideTemp)/&
                                  ((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
      theta(-1,j,k)=theta(0,j,k)-(theta(1,j,k)-theta(0,j,k))
     enddo
    enddo
   else    
    do k=1,nz
     do j=1,ny                      
      theta(0,j,k)=theta(1,j,k)-WsideTemp*dxU(0)
      theta(-1,j,k)=theta(1,j,k)-WsideTemp*(dxU(0)+dxU(-1))
     enddo
    enddo
   endif

   if (TBtypeE==DIRICHLET) then
    do k=1,nz
     do j=1,ny
      theta(nx+1,j,k)=EsideTemp!-(theta(nx,j,k)-EsideTemp)
      theta(nx+2,j,k)=EsideTemp!-(theta(nx-1,j,k)-EsideTemp)
     enddo
    enddo
   else if (TBtypeE==PERIODIC) then
    do k=1,nz
     do j=1,ny
      theta(nx+1,j,k)=theta(1,j,k)
      theta(nx+2,j,k)=theta(2,j,k)
     enddo
    enddo
   else if (TBtypeE==NEUMANN) then
    do k=1,nz
     do j=1,ny
      theta(nx+1,j,k)=theta(nx,j,k)+EsideTemp*dxU(nx+1)
      theta(nx+2,j,k)=theta(nx,j,k)+EsideTemp*(dxU(nx+1)+dxU(nx+2))
     enddo
    enddo
   else if (TBtypeE==CONSTFLUX) then
    do k=1,nz
     do j=1,ny
!       theta(nx+1,j,k)=(theta(nx,j,k)*(U(nx,j,k)/2._KND+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+EsideTemp)/&
!        (-U(nx,j,k)/2._KND+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
      theta(nx+1,j,k)=(theta(nx,j,k)*((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+EsideTemp)/&
       ((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
      theta(nx+2,j,k)=theta(nx+1,j,k)-(theta(nx,j,k)-theta(nx+1,j,k))
     enddo
    enddo   
   else
    do k=1,nz
     do j=1,ny
      theta(nx+1,j,k)=theta(nx,j,k)+EsideTemp*dxU(nx+1)
      theta(nx+2,j,k)=theta(nx,j,k)+EsideTemp*(dxU(nx+1)+dxU(nx+2))
     enddo
    enddo
   endif

  
   if (TBtypeS==DIRICHLET) then
    do k=1,nz
     do i=-1,nx+2
      theta(i,0,k)=SsideTemp!-(theta(i,1,k)-SsideTemp)
      theta(i,-1,k)=SsideTemp!-(theta(i,2,k)-SsideTemp)
     enddo
    enddo
   else if (TBtypeS==PERIODIC) then
    do k=1,nz
     do i=-1,nx+2
      theta(i,0,k)=theta(i,ny,k)
      theta(i,-1,k)=theta(i,ny-1,k)
     enddo
    enddo
   else if (TBtypeS==NEUMANN) then
    do k=1,nz
     do i=-1,nx+2
      theta(i,0,k)=theta(i,1,k)-SsideTemp*dyV(0)
      theta(i,-1,k)=theta(i,1,k)-SsideTemp*(dyV(0)+dyV(-1))
     enddo
    enddo
   else if (TBtypeS==CONSTFLUX) then
    do k=1,nz
     do i=-1,nx+2
!       theta(i,0,k)=(theta(i,1,k)*(V(i,0,k)/2._KND-(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+SsideTemp)/&
!                                   (V(i,0,k)/2._KND+(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
      theta(i,0,k)=(theta(i,1,k)*((TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+SsideTemp)/&
                                  ((TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
      theta(i,-1,k)=theta(i,0,k)-(theta(i,1,k)-theta(i,0,k))
     enddo
    enddo
   else    
    do k=1,nz
     do i=-1,nx+2
      theta(i,0,k)=theta(i,1,k)-SsideTemp*dyV(0)
      theta(i,-1,k)=theta(i,1,k)-SsideTemp*(dyV(0)+dyV(-1))
     enddo
    enddo
   endif

   if (TBtypeN==DIRICHLET) then
    do k=1,nz
     do i=-1,nx+2
       theta(i,ny+1,k)=NsideTemp!-(theta(i,ny,k)-NsideTemp)
       theta(i,ny+2,k)=NsideTemp!-(theta(i,ny-1,k)-NsideTemp)
      enddo
     enddo
   else if (TBtypeN==PERIODIC) then
    do k=1,nz
     do i=-1,nx+2
      theta(i,ny+1,k)=theta(i,1,k)
      theta(i,ny+2,k)=theta(i,2,k)
     enddo
    enddo
   else if (TBtypeN==NEUMANN) then
    do k=1,nz
     do i=-1,nx+2
      theta(i,ny+1,k)=theta(i,ny,k)+NsideTemp*dyV(ny+1)
      theta(i,ny+2,k)=theta(i,ny,k)+NsideTemp*(dyV(ny+1)+dyV(ny+2))
     enddo
    enddo
   else if (TBtypeN==CONSTFLUX) then
    do k=1,nz
     do i=-1,nx+2
!       theta(i,ny+1,k)=(theta(i,ny,k)*(V(i,ny,k)/2._KND+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+NsideTemp)/&
!        (-V(i,ny,k)/2._KND+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
      theta(i,ny+1,k)=(theta(i,ny,k)*((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+NsideTemp)/&
       ((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
      theta(i,ny+2,k)=theta(i,ny+1,k)-(theta(i,ny,k)-theta(i,ny+1,k))
     enddo
    enddo   
   else
    do k=1,nz
     do i=-1,ny+2
      theta(i,ny+1,k)=theta(i,ny,k)+NsideTemp*dyV(ny+1)
      theta(i,ny+2,k)=theta(i,ny,k)+NsideTemp*(dyV(ny+1)+dyV(ny+2))
     enddo
    enddo
   endif



!    if (TBtypeB==DIRICHLET)  then
!     do j=-1,ny+2
!      do i=-1,nx+2
!       theta(i,j,0)=BsideTemp!-(theta(i,j,1)-BsideTemp)
!       theta(i,j,-1)=BsideTemp!-(theta(i,j,2)-BsideTemp)
!      enddo
!     enddo
!    else 
   if (TBtypeB==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0)=theta(i,j,nz)
      theta(i,j,-1)=theta(i,j,nz-1)
     enddo
    enddo
   else if (TBtypeB==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0)=theta(i,j,1)-BsideTemp*dzW(0)
      theta(i,j,-1)=theta(i,j,1)-BsideTemp*(dzW(0)+dzW(-1))
     enddo
    enddo
   else if (TBtypeB==CONSTFLUX.or.TBtypeB==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
      if (abs(BsideTFLArr(i,j))<tiny(1._KND).and.TDiff(i,j,1)<1.1_KND/(Prandtl*Re).and.TBtypeB==DIRICHLET) then
       theta(i,j,0)=BsideTArr(i,j)
      else
       theta(i,j,0)=theta(i,j,1)+BsideTFLArr(i,j)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._KND))
      endif
      theta(i,j,-1)=theta(i,j,0)-(theta(i,j,1)-theta(i,j,0))
     enddo
    enddo
   else    
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0)=theta(i,j,1)-BsideTemp*dzW(0)
      theta(i,j,-1)=theta(i,j,1)-BsideTemp*(dzW(0)+dzW(-1))
     enddo
    enddo
   endif

   if (TBtypeT==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
       theta(i,j,nz+1)=TsideTemp!-(theta(i,j,nz)-TsideTemp)
       theta(i,j,nz+2)=TsideTemp!-(theta(i,j,nz-1)-TsideTemp)
      enddo
     enddo
   else if (TBtypeT==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1)=theta(i,j,1)
      theta(i,j,nz+2)=theta(i,j,2)
     enddo
    enddo
   else if (TBtypeT==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1)=theta(i,j,nz)+TsideTemp*dzW(nz+1)
      theta(i,j,nz+2)=theta(i,j,nz)+TsideTemp*(dzW(nz+1)+dzW(nz+2))
     enddo
    enddo
   else if (TBtypeT==CONSTFLUX) then
    do j=-1,ny+2
     do i=-1,nx+2
!       theta(i,j,nz+1)=(theta(i,j,nz)*(W(i,j,nz)/2._KND+(TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2*dzW(nz)))+TsideTemp)/&
!        (-W(i,j,nz)/2._KND+(TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2*dzW(nz)))
      theta(i,j,nz+1)=theta(i,j,nz)-TsideTemp*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._KND))
      theta(i,j,nz+2)=theta(i,j,nz+1)-(theta(i,j,nz)-theta(i,j,nz+1))
     enddo
    enddo   
   else
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1)=theta(i,j,nz)+TsideTemp*dzW(nz+1)
      theta(i,j,nz+2)=theta(i,j,nz)+TsideTemp*(dzW(nz+1)+dzW(nz+2))
     enddo
    enddo
   endif
  endsubroutine BOUND_TEMP





  pure subroutine BOUND_Visc(Nu)
  real(KND),intent(inout):: Nu(-1:,-1:,-1:)
  integer i,j,k,nx,ny,nz

  nx=Prnx
  ny=Prny
  nz=Prnz

  if (BtypeE==PERIODIC) then
   do k=1,nz
    do j=1,ny
      Nu(0,j,k)=Nu(nx,j,k)
      Nu(nx+1,j,k)=Nu(1,j,k)
    enddo
   enddo   
  else  
   do k=1,nz
    do j=1,ny
      Nu(0,j,k)=Nu(1,j,k)
      Nu(nx+1,j,k)=Nu(nx,j,k)
    enddo
   enddo
  endif

  if (BtypeN==PERIODIC) then
   do k=1,nz
    do i=1,nx
      Nu(i,0,k)=Nu(i,ny,k)
      Nu(i,ny+1,k)=Nu(i,1,k)
    enddo
   enddo   
  else  
   do k=1,nz
    do i=1,nx
      Nu(i,0,k)=Nu(i,1,k)
      Nu(i,ny+1,k)=Nu(i,ny,k)
    enddo
   enddo
  endif 

  if (BtypeT==PERIODIC) then
   do j=1,ny
    do i=1,nx
      Nu(i,j,0)=Nu(i,j,nz)
      Nu(i,j,nz+1)=Nu(i,j,1)
    enddo
   enddo
  else
   do j=1,ny
    do i=1,nx
      Nu(i,j,0)=Nu(i,j,1)
      Nu(i,j,nz+1)=Nu(i,j,nz)
    enddo
   enddo
  endif 
  endsubroutine BOUND_Visc
  
















  subroutine ScalFlSOURC(SCAL,sctype)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
   real(KND),intent(inout):: SCAL(-1:,-1:,-1:)
   integer,intent(in):: sctype
   real(KND) intscal,intTDiff
   type(TScalFlIBPoint),pointer:: SIBP
   integer xi,yj,zk,dirx,diry,dirz,n1,i

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   elseif (sctype==3) then
    call BOUND_Visc(SCAL)
   else 
    call BOUND_PASSSCALAR(SCAL)
   endif

  if (associated(FirstScalFlIBPoint)) then
   SIBP => FirstScalFlIBPoint
   do
    xi=SIBP%x
    yj=SIBP%y
    zk=SIBP%z
    intscal=0
    intTDiff=0
    n1=0
    do i=1,SIBP%interp
     intscal=intscal+SIBP%intcoef(i)*SCAL(SIBP%intpointi(i),SIBP%intpointj(i),SIBP%intpointk(i))
     if (sctype/=3.and.SIBP%intcoef(i)/=0._KND) then
                                 intTDiff=intTDiff+TDiff(SIBP%intpointi(i),SIBP%intpointj(i),SIBP%intpointk(i))
                                 n1=n1+1
     endif
    enddo

    if (sctype/=3) then
     intTDiff=intTDiff+TDiff(xi,yj,zk)
     intTDiff=intTDiff/(n1+1)

     if (intTDiff>0)  intscal=intscal+SIBP%Flux*SIBP%dist/intTDiff
    endif

    SIBP%ScalSourc=(intscal-SCAL(xi,yj,zk))/dt
    if (associated(SIBP%next)) then
     SIBP=>SIBP%next
    else
     exit
    endif
   enddo
  endif
  end subroutine ScalFlSOURC


    












  pure real(DBL) function AirDensity(press,temp)
  real(KND),intent(in):: press,temp
  AirDensity=press/(287.05_DBL*temp)
  endfunction AirDensity


  pure real(DBL) function AirDynVisc(temp)
  real(KND),intent(in):: temp
  AirDynVisc=1.85e-5_DBL
  endfunction AirDynVisc



  pure real(DBL) function CorrFactor(dp,press,temp)
  real(KND),intent(in):: dp,press,temp
  real(DBL) l
  l=MeanFreePath(press,temp)
  CorrFactor=1+(2*l/dp)*(1.257_KND+0.4_KND*exp(-0.55_KND*dp/l))
  endfunction CorrFactor


  pure real(DBL) function MeanFreePath(press,temp)
  real(KND),intent(in):: press,temp
  MeanFreePath=2.24e-5_DBL*temp/press
  endfunction MeanFreePath


  pure real(DBL) function BrownDiffusivity(dp,press,temp)
  real(KND),intent(in):: dp,press,temp
  real(DBL) C
  C=Corrfactor(dp,press,temp)
  BrownDiffusivity=BoltzC*temp*C/(3._KND*pi*AirDynVisc(temp)*dp)
  endfunction BrownDiffusivity


  pure real(KND) function BrownEff(dp,press,temp)
  real(KND),intent(in):: dp,press,temp
  real(KND) Sc
  Sc=(AirDynVisc(temp)/AirDensity(press,temp))/BrownDiffusivity(dp,press,temp)
  BrownEff=Sc**(-0.54_KND)
  endfunction BrownEff


  pure real(KND) function ImpactEff(dp,press,temp,ustar,visc)
  real(KND),intent(in):: dp,press,temp,ustar,visc
  real(KND) St
  St=SedimVelocity2(dp,press,temp)*ustar**2/visc
  ImpactEff=St**2/(400+St**2)
  endfunction ImpactEff


  pure real(KND) function SedimVelocity2(dp,press,temp)
  real(KND),intent(in):: dp,temp,press
  real(KND) C
  C=CorrFactor(dp,press,temp)
  SedimVelocity2=AirDensity(press,temp)*dp**2*9.81_KND*C/(18._KND*AirDynVisc(temp))
  endfunction SedimVelocity2

  pure real(DBL) function SedimVelocity(dp,rhop,press,temp)
  real(KND),intent(in):: dp,rhop,temp,press
  real(DBL) C,us,rho,mu
  rho=AirDensity(press,temp)
  mu=AirDynVisc(temp)
  C=CorrFactor(dp,press,temp)
  us=1+(0.42_DBL*C**2*rho*rhop/(108*mu**2))*dp**3*(1-rho/rhop)*9.81_DBL
  us=sqrt(us)
  us=(12._DBL*mu/(0.42_DBL*C*rho*dp))*(us-1._DBL)
  SedimVelocity=us
  endfunction SedimVelocity


  pure real(KND) function DyerH(zL)
  real(KND),intent(in):: zL
  if (zL>=0) then
   DyerH=1._KND+5._KND*zl
  else
   DyerH=1._KND/sqrt(1._KND-16._KND*zl)
  endif
  endfunction DyerH



  pure real(KND) function AerResist(z,z0,zL,ustar,visc)
  real(KND),intent(in):: z,z0,zL,ustar,visc
  real(KND),parameter:: yplcrit=11.225

  if (z>z0.and.z0>0) then
   AerResist=(log(z/z0)-DyerH(zl))/(0.4_KND*ustar)
  else
   if ((z*ustar/visc)<yplcrit) then
     AerResist=(z*ustar/visc)
   else
    AerResist=(log(abs(ustar*z/visc))/0.4_KND+5.2_KND)/ustar
   endif
  endif
  endfunction AerResist

  pure real(KND) function DepositionVelocity3(dp,rhop,press,temp,z,z0,zL,ustar) !EMRAS recommended values
  real(KND),intent(in):: dp,rhop,press,temp,z,z0,zL,ustar

   if (dp<1e-6) then
    DepositionVelocity3=0.5e-4
   elseif (dp<2e-6) then
    DepositionVelocity3=1.5e-4
   elseif (dp<10e-6) then
    DepositionVelocity3=10e-4
   else
    DepositionVelocity3=80e-4
   endif
  endfunction DepositionVelocity3



  pure real(KND) function DepositionVelocity2(dp,press,temp,z,z0,zL,ustar)
  real(KND),intent(in):: dp,press,temp,z,z0,zL,ustar
  real(KND) SurfResist,St,visc

   visc=AirDynVisc(temp)/AirDensity(press,temp)
   St=SedimVelocity2(dp,press,temp)*ustar**2/(visc)
   SurfResist=1._KND/(3._KND*ustar*(BrownEff(dp,press,temp)+ImpactEff(dp,press,temp,ustar,visc)))
   DepositionVelocity2=SedimVelocity2(dp,press,temp)+1._KND/(AerResist(z,z0,zL,ustar,visc)+SurfResist)
  endfunction DepositionVelocity2


  pure real(KND) function DepositionVelocity(dp,rhop,press,temp,z,z0,zL,ustar) !Kharchenko
  real(KND),intent(in):: dp,rhop,press,temp,z,z0,zL,ustar
  real(DBL) visc,us,Intz,Intexp,BD,tp
  real(DBL),parameter:: zexp=0.01

   us=SedimVelocity(dp,rhop,press,temp)
   visc=AirDynVisc(temp)/AirDensity(press,temp)
   tp=(us/9.81_KND)*ustar**2/visc
   BD=BrownDiffusivity(dp,press,temp)

   if (zl>=0) then
    Intz=-log(z/zexp+6*(zl-zexp*zl/z))/Karman
   else
    Intz=((sqrt(1-9*zl)-1)*(sqrt(1-9*zexp*zl/z)+1))
    Intz=Intz/((sqrt(1-9*zl)+1)*(sqrt(1-9*zexp*zl/z)-1))
    Intz=-Intz/Karman
   endif

   Intexp=-367.8_KND
   Intexp=Intexp+16.4*log(visc/BD)
   Intexp=Intexp-0.73*log(100*dp)*log(1e4*BD)-0.5*(log(100*dp))**2
   Intexp=Intexp+0.13*log(0.03/z0)
   Intexp=Intexp+0.25*log(0.2/ustar)*(1-0.2*log(0.03/z0))
   Intexp=Intexp-0.03*log(tp)*log(0.03/z0)
   Intexp=Intexp-32.7*log(100*dp)
   Intexp=-exp(Intexp)

   DepositionVelocity=us/(1._KND-exp((us/ustar)*(Intexp+Intz)))
  endfunction DepositionVelocity


  pure real(KND) function DepositionFlux(WMP,conc,partdiam,rhop)
  type(WMPoint),intent(in):: WMP
  real(KND),intent(in):: conc,partdiam,rhop
  real(KND):: press,temp,depvel

   press=101300
   temp=temperature_ref
   if ((.2_KND*WMP%distz)**2>(WMP%distx)**2+(WMP%disty)**2.and.WMP%distz>0) then
    depvel=DepositionVelocity3(partdiam,rhop,press,temp,WMP%distz,WMP%z0,0._KND,WMP%ustar)
   else
    depvel=0
   endif
   DepositionFlux=abs(depvel)*abs(conc)

  endfunction DepositionFlux

  subroutine Deposition(SCAL,coef)
  real(KND),dimension(-1:,-1:,-1:,1:),intent(inout):: SCAL
  real(KND),intent(in):: coef
  type(WMPoint),pointer:: WMP 
  integer i
  real(KND) deptmp
   if (associated(FirstWMPoint)) then
   WMP => FirstWMPoint
   do
    if (allocated(WMP%depscalar)) then    
    if (partdistrib>0) then
     do i=1,partdistrib
      deptmp=abs(DepositionFlux(WMP,SCAL(WMP%x,WMP%y,WMP%z,1)*percdistrib(i),partdiam(i),partrho(i)))&
          *coef*dt*dxPr(WMP%x)*dyPr(WMP%y)
      WMP%depscalar(1)=WMP%depscalar(1)+deptmp/(dxPr(WMP%x)*dyPr(WMP%y)*dzPr(WMP%z))
      SCAL(WMP%x,WMP%y,WMP%z,1)=SCAL(WMP%x,WMP%y,WMP%z,1)-deptmp/(dxPr(WMP%x)*dyPr(WMP%y)*dzPr(WMP%z))
     enddo
    else
     do i=1,computescalars
      deptmp=abs(DepositionFlux(WMP,SCAL(WMP%x,WMP%y,WMP%z,i),partdiam(i),partrho(i)))&
          *coef*dt*dxPr(WMP%x)*dyPr(WMP%y)
      WMP%depscalar(i)=WMP%depscalar(i)+deptmp/(dxPr(WMP%x)*dyPr(WMP%y)*dzPr(WMP%z))
      SCAL(WMP%x,WMP%y,WMP%z,i)=SCAL(WMP%x,WMP%y,WMP%z,i)-deptmp/(dxPr(WMP%x)*dyPr(WMP%y)*dzPr(WMP%z))
     enddo
    endif
    endif
    if (associated(WMP%next)) then
     WMP=>WMP%next
    else
     exit
    endif
   enddo
   endif
  endsubroutine Deposition


  subroutine Gravsettling(SCAL,coef)
  real(KND),dimension(-1:,-1:,-1:,1:):: SCAL
  integer i,j,k,l
  real(KND),dimension(Prnx,Prny,Prnz)::flux
  real(KND):: coef,press,temp,us

  press=101300
  temp=temperature_ref
  if (partdistrib==0) then
   do l=1,computescalars
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       us=SedimVelocity(partdiam(l),partrho(l),press,temp)
       flux(i,j,k)=us*SCAL(i,j,k+1,l)*coef*dt*dxPr(i)*dyPr(j)
      enddo
     enddo
    enddo
    do k=1,Prnz-1
     do j=1,Prny
      do i=1,Prnx
       SCAL(i,j,k+1,l)=SCAL(i,j,k+1,l)-flux(i,j,k)/(dxPr(i)*dyPr(j)*dzPr(k+1))
       SCAL(i,j,k,l)=SCAL(i,j,k,l)+flux(i,j,k)/(dxPr(i)*dyPr(j)*dzPr(k))
      enddo
     enddo
    enddo
   enddo
  endif
  endsubroutine Gravsettling  


  pure real(KND) function Prt(i,j,k,U,V,temperature)
  integer,intent(in):: i,j,k
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V
  real(KND),dimension(-1:,-1:,-1:),intent(in):: temperature
  
!   if (buoyancy>0) then
!    Prt=0.8_KND+min(max(3._KND*Rig(i,j,k,U,V,temperature),0._KND),sqrt(huge(1._KND)))
!   else
   Prt=debugparam!0.6_KND
!   endif
  endfunction Prt


  pure real(KND) function Rig(i,j,k,U,V,temperature)
  integer,intent(in):: i,j,k
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V
  real(KND),dimension(-1:,-1:,-1:),intent(in):: temperature
  real(KND) num,denom

  num=(grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
  denom=((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
  denom=denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
  if (abs(denom)>tiny(1._KND)*100) then
   Rig=num/denom
  else
   Rig=0
  endif
  endfunction Rig



  elemental real(KND) function FLUXLIMITER(r)
  real(KND),intent(in):: r
  real(KND) K

  K=(1+2._KND*r)/3._KND  !assuming kappa=1/3 scheme (Koren, 1993; Hundsdorfer, 1994)
  if (limparam>0) then
   FLUXLIMITER=max(0._KND,min(2._KND*r,min(limparam,K)))
  elseif (limparam<0) then
   FLUXLIMITER=K
  else
   FLUXLIMITER=0
  endif
  endfunction FLUXLIMITER
  

  elemental real(KND) function LIMITER(a,b)
  real(KND),intent(in):: a,b
  real(KND) R
  real(KND),parameter:: epsil=0.00001_KND
  
  if (limitertype==minmodlim) then
   LIMITER=MINMOD(a,b)
  elseif (limitertype==extminmodlim) then
   LIMITER=MINMOD(limparam*a,limparam*b,(a+b)/2._KND)
  elseif (limitertype==gammalim) then
   if (abs(b)>epsil.and.abs(a)>epsil) then
    R=a/b
    LIMITER=b*GAMMAL(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==vanalbadalim) then
   if (abs(b)>epsil) then
    R=a/b
    LIMITER=b*VANALBADA(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==vanleerlim) then
   if (abs(b)>epsil) then
    R=a/b
    LIMITER=b*VANLEER(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==superbeelim) then
   if (abs(b)>epsil) then
    R=a/b
    LIMITER=b*SUPERBEE(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==0) then
   LIMITER=0
  else
   LIMITER=(a+b)/2._KND
  endif
  endfunction LIMITER

  elemental real(KND) function HEAV(x)
  real(KND),intent(in):: x
  HEAV=(1._KND+SIGN(1._KND,x))/2._KND
  endfunction HEAV

  elemental real(KND) function GAMMAL(R)
  real(KND),intent(in):: R
  GAMMAL=((1-limparam)/limparam)*R*(HEAV(R) - HEAV(r-(limparam/(1-limparam)))) + HEAV(r-(limparam/(1-limparam)))
  endfunction GAMMAL

  elemental real(KND) function VANALBADA(R)
  real(KND),intent(in):: R
  VANALBADA=(R+R*R)/(1+R*R)
  endfunction VANALBADA

  elemental real(KND) function VANLEER(R)
  real(KND),intent(in):: R
  VANLEER=(R+abs(R))/(1+abs(R))
  endfunction VANLEER

  elemental real(KND) function SUPERBEE(R)
  real(KND),intent(in):: R
  SUPERBEE=MAX(0._KND,MAX(MIN(2*R,1._KND),MIN(R,2._KND)))
  endfunction SUPERBEE

  elemental real(KND) function MINMOD2(a,b)
  real(KND),intent(in):: a,b
      MINMOD2=(SIGN(1._KND,a)+SIGN(1._KND,b))*MIN(ABS(a),ABS(b))/2._KND
  endfunction MINMOD2

  elemental real(KND) function MINMOD3(a,b,c)
  real(KND),intent(in):: a,b,c
    if ((a>0).and.(b>0).and.(c>0)) then  
        MINMOD3=MIN(a,b,c)
    elseif ((a<0).and.(b<0).and.(c<0)) then
        MINMOD3=MAX(a,b,c)
    else
        MINMOD3=0
    endif
  endfunction MINMOD3

end module SCALARS
