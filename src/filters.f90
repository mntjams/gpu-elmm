module Filters

  use Parameters

  implicit none

  private

  public filtertype, filter_ratios, Filter

  integer filtertype

  real(knd), parameter :: filter_ratios(0:1) = [ 1.0_knd, 2._knd]

 
  contains

  
    subroutine Filter3ord(U,Utype,dir)
      !3D filter with 3 vanishing moments width 2 delta x
      !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
      real(knd),dimension(-2:,-2:,-2:),intent(inout)  :: U
      integer,dimension(-2:,-2:,-2:),intent(inout)  :: Utype
      integer,intent(in) :: dir
      real(knd),parameter :: w(-3:3) = (/ -1._knd/32, 0._knd, 9._knd/32, 1._knd/2, 9._knd/32, 0._knd, -1._knd/32 /)
      real(knd) S,q
      integer i,j,k,ii,jj,kk
      integer mini,minj,mink
      integer maxi,maxj,maxk

      real(knd) :: tmp(-2:ubound(U,dir))
      
      mini = lbound(U,1)+3
      maxi = ubound(U,1)-3
      minj = lbound(U,2)+3
      maxj = ubound(U,2)-3
      mink = lbound(U,3)+3
      maxk = ubound(U,3)-3

      !The condition is intentionally ==0( not <=0)
      !  to avoid points nearest to the boundary.
      
      if (dir==1) then
        !$omp parallel do private(i,j,k,ii,S,q,tmp) shared(U,Utype,mini,maxi,minj,maxj,mink,maxk) schedule(runtime)
        do k = mink,maxk
         do j = minj,maxj
           tmp = U(:,j,k)
           do i = 1,maxi
             if (Utype(i,j,k) == 0) then
               S = 0
               q = 0
               do ii = i-3,i+3
                 if (Utype(ii,j,k)<=0) then
                   S = S + w(ii-i) * tmp(ii)
                   q = q + w(ii-i)
                 end if
               end do
               if (abs(q)>0.74) U(i,j,k) = S / q
             end if
           end do
         end do
        end do
        !$omp end parallel do
      else if (dir==2) then
        !$omp parallel do private(i,j,k,jj,S,q,tmp) shared(U,Utype,mini,maxi,minj,maxj,mink,maxk) schedule(runtime)
        do k = mink,maxk
         do i = mini,maxi
           tmp = U(i,:,k)
           do j = 1,maxj
             if (Utype(i,j,k) == 0) then
               S = 0
               q = 0
               do jj = j-3,j+3
                 if (Utype(i,jj,k)<=0) then
                   S = S + w(jj-j) * tmp(jj)
                   q = q + w(jj-j)
                 end if
               end do
               if (abs(q)>.74) U(i,j,k) = S / q
             end if
           end do
         end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i,j,k,kk,S,q,tmp) shared(U,Utype,mini,maxi,minj,maxj,mink,maxk) schedule(runtime)
        do j = minj,maxj
         do i = mini,maxi
           tmp = U(i,j,:)
           do k = 1,maxk
             if (Utype(i,j,k) == 0) then               
               S = 0
               q = 0
               do kk = k-3,k+3
                 if (Utype(i,j,kk)<=0) then
                   S = S + w(kk-k) * tmp(kk)
                   q = q + w(kk-k)
                 end if
               end do
               if (abs(q)>.74) U(i,j,k) = S / q
             end if
           end do
         end do
        end do
        !$omp end parallel do
      end if
    end subroutine Filter3ord



    subroutine Filter(U,Utype)    !Calls a selected filter  U2 <- filt(U1)
      real(knd),dimension(:,:,:),intent(inout) :: U
      integer,dimension(-2:,-2:,-2:),intent(inout)  :: Utype

      if (filtertype==1) then
        if (Prnx==1) then
          call Filter3ord(U,Utype,2)
          call Filter3ord(U,Utype,3)
        else if (Prny==1) then
          call Filter3ord(U,Utype,1)
          call Filter3ord(U,Utype,3)
        else if (Prnz==1) then
          call Filter3ord(U,Utype,3)
          call Filter3ord(U,Utype,3)
        else
          call Filter3ord(U,Utype,1)
          call Filter3ord(U,Utype,2)
          call Filter3ord(U,Utype,3)
        end if
      end if

    end subroutine Filter

end module Filters
