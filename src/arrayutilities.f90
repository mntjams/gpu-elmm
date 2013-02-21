! module MyBLAS
!   !Warning, all supplied arrays must be contiguous! This attribute not used due to compatibility with Oracle Solrais Studio 12.3
! 
!   use iso_c_binding
! 
!   interface
!     subroutine saxpy_(n,a,px,incx,py,incy) bind(C,name="saxpy_")
!       import
!       integer,intent(in)      :: n
!       real(c_float),intent(in) :: a
!       integer(c_intptr_t),value       :: px, py
!       integer,intent(in)      :: incx, incy
!     end subroutine
! 
!     subroutine daxpy_(n,a,px,incx,py,incy) bind(C,name="daxpy_")
!       import
!       integer,intent(in)        :: n
!       real(c_double),intent(in) :: a
!       integer(c_intptr_t),value         :: px, py
!       integer,intent(in)        :: incx, incy
!     end subroutine
!   end interface
! 
!   interface axpy
!     module procedure saxpy1D
!     module procedure daxpy1D
!     module procedure saxpy2D
!     module procedure daxpy2D
!     module procedure saxpy3D
!     module procedure daxpy3D
!   end interface
! 
!   contains
! 
!     subroutine saxpy1D(a, x, y)
!       real(c_float),target :: a, x(:), y(:)
!       integer(c_intptr_t) :: px, py
!       !y = y + a * x
!       px = loc(x(1))
!       py = loc(y(1))
!       call saxpy_(size(x),a,px,1,py,1)
!     end subroutine
! 
!     subroutine daxpy1D(a, x, y)
!       real(c_double),target :: a, x(:), y(:)
!       integer(c_intptr_t) :: px, py
!       !y = y + a * x
!       px = loc(x(1))
!       py = loc(y(1))
!       call daxpy_(size(x),a,px,1,py,1)
!     end subroutine
! 
!     subroutine saxpy2D(a, x, y)
!       real(c_float),target :: a, x(:,:), y(:,:)
!        integer(c_intptr_t) :: px, py
!      !y = y + a * x
!       px = loc(x(1,1))
!       py = loc(y(1,1))
!       call saxpy_(size(x),a,px,1,py,1)
!     end subroutine
! 
!     subroutine daxpy2D(a, x, y)
!       real(c_double),target :: a, x(:,:), y(:,:)
!       integer(c_intptr_t) :: px, py
!       !y = y + a * x
!       px = loc(x(1,1))
!       py = loc(y(1,1))
!       call daxpy_(size(x),a,px,1,py,1)
!     end subroutine
! 
!     subroutine saxpy3D(a, x, y)
!       real(c_float),target :: a, x(:,:,:), y(:,:,:)
!       integer(c_intptr_t) :: px, py
!       !y = y + a * x
!       px = loc(x(1,1,1))
!       py = loc(y(1,1,1))
!       call saxpy_(size(x),a,px,1,py,1)
!     end subroutine
! 
!     subroutine daxpy3D(a, x, y)
!       real(c_double),target :: a, x(:,:,:), y(:,:,:)
!       integer(c_intptr_t) :: px, py
!       !y = y + a * x
!       px = loc(x(1,1,1))
!       py = loc(y(1,1,1))
!       call daxpy_(size(x),a,px,1,py,1)
!     end subroutine
! 
! end module MyBLAS



module ArrayUtilities
  use Kinds
!   use MyBLAS

  implicit none

  private
  public exchange_alloc, assign, add, &
         set, multiply, add_multiplied, &
         multiply_and_add_scalar, reciprocal

  interface exchange_alloc
    module procedure exchange_alloc_1D
    module procedure exchange_alloc_2D
    module procedure exchange_alloc_3D
  end interface

  interface assign
    module procedure assign_1D
    module procedure assign_2D
    module procedure assign_3D
  end interface

  interface add
    module procedure add_1D
    module procedure add_2D
    module procedure add_3D
  end interface

  interface set
    module procedure set_1D
    module procedure set_2D
    module procedure set_3D
  end interface

  interface multiply
    module procedure multiply_1D
    module procedure multiply_2D
    module procedure multiply_3D
  end interface

  interface add_multiplied
    module procedure add_multiplied_1D
    module procedure add_multiplied_2D
    module procedure add_multiplied_3D
  end interface

  interface multiply_and_add_scalar
    module procedure multiply_and_add_scalar_1D
    module procedure multiply_and_add_scalar_2D
    module procedure multiply_and_add_scalar_3D
  end interface

  interface reciprocal
    module procedure reciprocal_1D
    module procedure reciprocal_2D
    module procedure reciprocal_3D
  end interface

  contains

    subroutine exchange_alloc_1D(A,B)
      real(KND),allocatable,intent(inout) :: A(:),B(:)
      real(KND),allocatable :: tmp(:)

      call move_alloc(A,tmp)
      call move_alloc(B,A)
      call move_alloc(tmp,B)
    end subroutine
    
    subroutine exchange_alloc_2D(A,B)
      real(KND),allocatable,intent(inout) :: A(:,:),B(:,:)
      real(KND),allocatable :: tmp(:,:)

      call move_alloc(A,tmp)
      call move_alloc(B,A)
      call move_alloc(tmp,B)
    end subroutine
    
    subroutine exchange_alloc_3D(A,B)
      real(KND),allocatable,intent(inout) :: A(:,:,:),B(:,:,:)
      real(KND),allocatable :: tmp(:,:,:)

      call move_alloc(A,tmp)
      call move_alloc(B,A)
      call move_alloc(tmp,B)
    end subroutine


    subroutine assign_1D(to,from)
      real(KND),intent(out) :: to(:)
      real(KND),intent(in)  :: from(:)
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = from(i)
      end do
      !$omp end parallel do
    end subroutine

    subroutine assign_2D(to,from)
      real(KND),intent(out) :: to(:,:)
      real(KND),intent(in)  :: from(:,:)
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = from(:,j)
      end do
      !$omp end parallel do
    end subroutine

    subroutine assign_3D(to,from)
      real(KND),intent(out) :: to(:,:,:)
      real(KND),intent(in)  :: from(:,:,:)
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = from(:,:,k)
      end do
      !$omp end parallel do
    end subroutine


    subroutine add_1D(to,from)
      real(KND),intent(inout) :: to(:)
      real(KND),intent(in)  :: from(:)
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = to(i) + from(i)
      end do
      !$omp end parallel do
    end subroutine

    subroutine add_2D(to,from)
      real(KND),intent(inout) :: to(:,:)
      real(KND),intent(in)  :: from(:,:)
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = to(:,j) + from(:,j)
      end do
      !$omp end parallel do
    end subroutine

    subroutine add_3D(to,from)
      real(KND),intent(inout) :: to(:,:,:)
      real(KND),intent(in)  :: from(:,:,:)
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = to(:,:,k) + from(:,:,k)
      end do
      !$omp end parallel do
!       call axpy(1._KND,from,to)
    end subroutine


    subroutine set_1D(to,from)
      real(KND),intent(out) :: to(:)
      real(KND),intent(in)  :: from
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = from
      end do
      !$omp end parallel do
    end subroutine

    subroutine set_2D(to,from)
      real(KND),intent(out) :: to(:,:)
      real(KND),intent(in)  :: from
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = from
      end do
      !$omp end parallel do
    end subroutine

    subroutine set_3D(to,from)
      real(KND),intent(out) :: to(:,:,:)
      real(KND),intent(in)  :: from
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = from
      end do
      !$omp end parallel do
    end subroutine
    

    subroutine multiply_1D(to,a)
      real(KND),intent(inout) :: to(:)
      real(KND),intent(in)  :: a
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = to(i) * a
      end do
      !$omp end parallel do
    end subroutine

    subroutine multiply_2D(to,a)
      real(KND),intent(inout) :: to(:,:)
      real(KND),intent(in)  :: a
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = to(:,j) * a
      end do
      !$omp end parallel do
    end subroutine

    subroutine multiply_3D(to,a)
      real(KND),intent(inout) :: to(:,:,:)
      real(KND),intent(in)  :: a
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = to(:,:,k) * a
!         call sscal(size(to,2)*size(to,1),a,to,1)
      end do
      !$omp end parallel do
    end subroutine


    subroutine add_multiplied_1D(to,from,a)
      real(KND),intent(inout) :: to(:)
      real(KND),intent(in)  :: from(:)
      real(KND),intent(in)  :: a
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = to(i) + from(i) * a
      end do
      !$omp end parallel do
    end subroutine

    subroutine add_multiplied_2D(to,from,a)
      real(KND),intent(inout) :: to(:,:)
      real(KND),intent(in)  :: from(:,:)
      real(KND),intent(in)  :: a
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = to(:,j) + from(:,j) * a
      end do
      !$omp end parallel do
    end subroutine

    subroutine add_multiplied_3D(to,from,a)
      real(KND),intent(inout) :: to(:,:,:)
      real(KND),intent(in)  :: from(:,:,:)
      real(KND),intent(in)  :: a
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = to(:,:,k) +from(:,:,k) * a
      end do
      !$omp end parallel do
!        call axpy(a,from,to)
    end subroutine


    subroutine multiply_and_add_scalar_1D(to,b,a)
      ! to = b * to + a
      real(KND),intent(inout) :: to(:)
      real(KND),intent(in)  :: b
      real(KND),intent(in)  :: a
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = b * to(i) + a
      end do
      !$omp end parallel do
    end subroutine

    subroutine multiply_and_add_scalar_2D(to,b,a)
      ! to = b * to + a
      real(KND),intent(inout) :: to(:,:)
      real(KND),intent(in)  :: b
      real(KND),intent(in)  :: a
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = b * to(:,j) + a
      end do
      !$omp end parallel do
    end subroutine

    subroutine multiply_and_add_scalar_3D(to,b,a)
      ! to = b * to + a
      real(KND),intent(inout) :: to(:,:,:)
      real(KND),intent(in)  :: b
      real(KND),intent(in)  :: a
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = b * to(:,:,k) + a
      end do
      !$omp end parallel do
    end subroutine


    subroutine reciprocal_1D(to,a)
      ! to = a / to
      real(KND),intent(inout) :: to(:)
      real(KND),intent(in)  :: a
      integer i

      !$omp parallel do
      do i=1,size(to)
        to(i) = a / to(i)
      end do
      !$omp end parallel do
    end subroutine

    subroutine reciprocal_2D(to,a)
      ! to = a / to
      real(KND),intent(inout) :: to(:,:)
      real(KND),intent(in)  :: a
      integer j

      !$omp parallel do
      do j=1,size(to,2)
        to(:,j) = a / to(:,j)
      end do
      !$omp end parallel do
    end subroutine

    subroutine reciprocal_3D(to,a)
      ! to = a / to
      real(KND),intent(inout) :: to(:,:,:)
      real(KND),intent(in)  :: a
      integer k

      !$omp parallel do
      do k=1,size(to,3)
        to(:,:,k) = a / to(:,:,k)
      end do
      !$omp end parallel do
    end subroutine


end module ArrayUtilities
