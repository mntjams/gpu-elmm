
  !$hmpp <tsteps> PrePoisson codelet
  subroutine PrePoisson_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,correctcompatibility,&
                              Btype,sideU,&
                              dt2,dxmin,dymin,dzmin,&
                              Uin,Vin,Win,U,V,W,RHS,uncompatibility,divergence)

    implicit none

#include "hmpp-include.f90"


    integer,intent(in)      :: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,correctcompatibility
    integer,intent(in)      :: Btype(6)
    real(KND),intent(in)    :: sideU(3,6)
    real(KND),intent(in)    :: dxmin,dymin,dzmin
    real(KND),dimension(-2:Prny+3,-2:Prnz+3),intent(in)       :: Uin,Vin,Win
    real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
    real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
    real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    real(KND),intent(out)   :: RHS(0:Prnx+1,0:Prny+1,0:Prnz+1)
    real(KND),intent(out)   :: uncompatibility,divergence
    real(KND),intent(in)    :: dt2
    real(KND) :: S,S2
    integer   :: i,j,k


    call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz,&
                         Btype,sideU,&
                         Uin,U,0)
    call BoundU_GPU(2,Vnx,Vny,Vnz,Prny,Prnz,&
                         Btype,sideU,&
                         Vin,V,0)
    call BoundU_GPU(3,Wnx,Wny,Wnz,Prny,Prnz,&
                         Btype,sideU,&
                         Win,W,0)


    if (correctcompatibility==1) then
      S=0
      !$hmppcg permute (k,i,j)
      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify(k,i), reduce(+:S)
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         S = S + (-((U(i,j,k)-U(i-1,j,k))/(dxmin)+(V(i,j,k)-V(i,j-1,k))/(dymin)+(W(i,j,k)-W(i,j,k-1))&
                      /(dzmin)))
        end do
       end do
      end do

      S=S*dxmin/(Prny*Prnz)

      uncompatibility = S

      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify(k,j)
      do k=-2,Unz+3
        do j=-2,Uny+3
          U(Unx+1,j,k)=U(Unx+1,j,k)+S
        end do
      end do
    end if

    S=0
    S2=0

    !$hmppcg permute (k,i,j)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i), reduce(+:S2)
    do k=1,Prnz            !divergence of U -> RHS
     do j=1,Prny
      do i=1,Prnx
           RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxmin)+(V(i,j,k)-V(i,j-1,k))/(dymin)+(W(i,j,k)-W(i,j,k-1))/(dzmin)
           S2=S2+abs(RHS(i,j,k))
           RHS(i,j,k)=RHS(i,j,k)/(dt2)
      end do
     end do
    end do

    divergence = S2/(Prnx*Prny*Prnz)

   end subroutine PrePoisson_GPU




  !$hmpp <tsteps> PostPoisson codelet
  subroutine PostPoisson_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                              Btype,sideU,&
                              dt2,dt3,dxmin,dymin,dzmin,&
                              Uin,Vin,Win,U,V,W,Pr,Phi)

    implicit none

#include "hmpp-include.f90"


    integer,intent(in)      :: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz
    integer,intent(in)      :: Btype(6)
    real(KND),intent(in)    :: sideU(3,6)
    real(KND),intent(in)    :: dxmin,dymin,dzmin
    real(KND),dimension(-2:Prny+3,-2:Prnz+3),intent(in) :: Uin,Vin,Win
    real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
    real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
    real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    real(KND),intent(inout) :: Pr(1:Unx+1,1:Vny+1,1:Wnz+1)
    real(KND),intent(inout) :: Phi(0:Prnx+1,0:Prny+1,0:Prnz+1)
    real(KND),intent(in)    :: dt2,dt3
    real(KND) :: Phiref,Au,Av,Aw,dxmin2,dymin2,dzmin2
    integer   :: i,j,k



    Phiref=0!Phi(Prnx/2,Prny/2,Prnz/2)

    !$hmppcg permute (k,i,j)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i)
    do k=0,Prnz+1
     do j=0,Prny+1
      do i=0,Prnx+1
       Phi(i,j,k) = Phi(i,j,k) - Phiref
      end do
     end do
    end do

    Au = dt2/dxmin
    Av = dt2/dymin
    Aw = dt2/dzmin

    dxmin2 = dxmin**2
    dymin2 = dymin**2
    dzmin2 = dzmin**2


    !$hmppcg permute (k,i,j)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
        U(i,j,k)=U(i,j,k)-Au*(Phi(i+1,j,k)-Phi(i,j,k))
      end do
     end do
    end do
    !$hmppcg permute (k,i,j)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i)
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
        V(i,j,k)=V(i,j,k)-Av*(Phi(i,j+1,k)-Phi(i,j,k))
      end do
     end do
    end do
    !$hmppcg permute (k,i,j)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i)
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
        W(i,j,k)=W(i,j,k)-Aw*(Phi(i,j,k+1)-Phi(i,j,k))
      end do
     end do
    end do

    !$hmppcg permute (k,i,j)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i)
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
        Pr(i,j,k)=Pr(i,j,k)+Phi(i,j,k)-dt3*(((Phi(i+1,j,k)-Phi(i,j,k))-(Phi(i,j,k)-Phi(i-1,j,k)))/dxmin2+&
                                            ((Phi(i,j+1,k)-Phi(i,j,k))-(Phi(i,j,k)-Phi(i,j-1,k)))/dymin2+&
                                            ((Phi(i,j,k+1)-Phi(i,j,k))-(Phi(i,j,k)-Phi(i,j,k-1)))/dzmin2)
      end do
     end do
    end do


    call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz,&
                         Btype,sideU,&
                         Uin,U,0)
    call BoundU_GPU(2,Vnx,Vny,Vnz,Prny,Prnz,&
                         Btype,sideU,&
                         Vin,V,0)
    call BoundU_GPU(3,Wnx,Wny,Wnz,Prny,Prnz,&
                         Btype,sideU,&
                         Win,W,0)

    call BoundPr_GPU(Unx,Vny,Wnz,Prnx,Prny,Prnz,Btype,Pr)


   end subroutine PostPoisson_GPU

