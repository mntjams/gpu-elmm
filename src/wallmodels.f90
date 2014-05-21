module Sort
  use iso_c_binding

  implicit none

  interface
    subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
      import
      type(c_ptr),value       :: array
      integer(c_size_t),value :: elem_count
      integer(c_size_t),value :: elem_size
      type(c_funptr),value    :: compare !int(*compare)(const void *, const void *)
    end subroutine qsort !standard C library qsort
  end interface
end module Sort

module WMPoint_types
  use Kinds

  type WMpoint   !points in which we apply wall model

    integer   :: xi
    integer   :: yj
    integer   :: zk

    real(knd) :: distx
    real(knd) :: disty
    real(knd) :: distz

    real(knd) :: area_factor ![m^-1] area of the solid wall divided by the volume of the cell(xi,yj,zk)

    real(knd) :: z0 = 0
    real(knd) :: ustar = 1

    real(knd) :: temperature = 0
    real(knd) :: temperature_flux = 0

    real(knd) :: moisture_flux = 0

    real(knd) :: div = 0 !divergence when zeroing out trans-boundary velocities

    real(knd) :: wallu = 0
    real(knd) :: wallv = 0
    real(knd) :: wallw = 0

    real(knd) :: albedo = 0.1 !for shortwave radiation - asphalt
    real(knd) :: emissivity = 0.9 !for longwave radiation - asphalt
    
    real(knd) :: evaporative_fraction = 0

    real(knd),allocatable:: depscalar(:)

  end type WMpoint


  type WMpointUVW   !points in which we apply wall model
  
    real(knd) :: flux

    integer   :: xi
    integer   :: yj
    integer   :: zk

    real(knd) :: distx
    real(knd) :: disty
    real(knd) :: distz

    real(knd) :: z0 = 0
    real(knd) :: ustar = 1

    real(knd) :: wallu = 0
    real(knd) :: wallv = 0
    real(knd) :: wallw = 0

    real(knd) :: temperature = 0
    real(knd) :: temperature_flux = 0

  end type WMpointUVW

end module

module WMPointLists
  use WMPoint_types
#define TYPEPARAM type(WMPoint)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module

module WMPointUVWLists
  use WMPoint_types
#define TYPEPARAM type(WMPointUVW)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module


module Wallmodels
  use Parameters
  use WMPoint_types
  use WMPointLists, only: WMPointList => list
  use WMPointUVWLists, only: WMPointUVWList => list

  implicit none

  private
  public WMPoint, WMpointUVW, &
         AddWMPoint, AddWMPointUVW, &
         MoveWMPointsToArray, GetOutsideBoundariesWM, InitWMMasks, &
         ComputeViscsWM, ComputeUVWFluxesWM, DivergenceWM, &
         InitTempFl, GroundDeposition, GroundUstar, wallmodeltype
#ifdef __HMPP
  public hmppWMpoint,WMPtoHMPP
#endif

  integer, parameter, public :: MINUSX = 1, PLUSX = 2, MINUSY = 3, PLUSY = 4, MINUSZ = 5, PLUSZ = 6

  type(WMPointUVWList) :: UxmWMPsL, UxpWMPsL, UymWMPsL, UypWMPsL, UzmWMPsL, UzpWMPsL, &
                          VxmWMPsL, VxpWMPsL, VymWMPsL, VypWMPsL, VzmWMPsL, VzpWMPsL, &
                          WxmWMPsL, WxpWMPsL, WymWMPsL, WypWMPsL, WzmWMPsL, WzpWMPsL

  type(WMPointList) :: WMPointsList

  type(WMPointUVW), dimension(:), allocatable, public :: &
      UxmWMpoints, UxpWMpoints, UymWMpoints, UypWMpoints, UzmWMpoints, UzpWMPoints, &
      VxmWMpoints, VxpWMpoints, VymWMpoints, VypWMpoints, VzmWMpoints, VzpWMpoints, &
      WxmWMpoints, WxpWMpoints, WymWMpoints, WypWMpoints, WzmWMpoints, WzpWMpoints

  type(WMPoint), dimension(:), allocatable, public :: WMPoints

  logical(slog),dimension(:,:,:),allocatable,public :: Scflx_mask, Scfly_mask, Scflz_mask

  logical(slog),dimension(:,:,:),allocatable,public :: Uflx_mask, Ufly_mask, Uflz_mask
  logical(slog),dimension(:,:,:),allocatable,public :: Vflx_mask, Vfly_mask, Vflz_mask
  logical(slog),dimension(:,:,:),allocatable,public :: Wflx_mask, Wfly_mask, Wflz_mask

  integer :: wallmodeltype

contains

#ifdef __HMPP
  elemental subroutine WMPtoHMPP(ToWMP,FromWMP)
    type(hmppWMPoint),intent(out) :: ToWMP
    type(WMPoint),intent(in)  :: FromWMP

    ToWMP%xi = FromWMP%xi
    ToWMP%yj = FromWMP%yj
    ToWMP%zk = FromWMP%zk

    ToWMP%distx = FromWMP%distx
    ToWMP%disty = FromWMP%disty
    ToWMP%distz = FromWMP%distz

    ToWMP%z0     = FromWMP%z0
    ToWMP%ustar  = FromWMP%ustar
    ToWMP%temperature   = FromWMP%temperature
    ToWMP%temperature_flux = FromWMP%temperature_flux

  end subroutine WMPtoHMPP
#endif


  elemental subroutine ComputeAreaFactor(p)
    type(WMpoint),intent(inout):: p
    real(knd) :: out_norm(3)
    !FIXME: Not accurate!!!
    if (all([-p%distx, -p%disty, -p%distz]==0)) return !HACK!!!
    out_norm = [-p%distx, -p%disty, -p%distz] / &
                     hypot(hypot(p%distx,p%disty),p%distz)
    p%area_factor = 1._knd/dot_product([dxmin,dymin,dzmin],abs(out_norm))
  end subroutine  



  subroutine AddWMPoint(WMP)
    type(WMPoint), intent(in) :: WMP

    call WMPointsList%add(WMP)

  end subroutine


  subroutine AddWMPointUVW(WMP, component, direction)
    type(WMPointUVW), intent(in) :: WMP
    integer, intent(in) :: component, direction

    select case (component)
      case (1)
        call add(WMP, direction, UxmWMPsL, UxpWMPsL, UymWMPsL, UypWMPsL, UzmWMPsL, UzpWMPsL)
      case (2)
        call add(WMP, direction, VxmWMPsL, VxpWMPsL, VymWMPsL, VypWMPsL, VzmWMPsL, VzpWMPsL)
      case (3)
        call add(WMP, direction, WxmWMPsL, WxpWMPsL, WymWMPsL, WypWMPsL, WzmWMPsL, WzpWMPsL)
    end select

  contains
    subroutine add(p, dir, xm, xp, ym, yp, zm, zp)
      type(WMPointUVW), intent(in) :: p
      integer, intent(in) :: dir
      type(WMPointUVWList), intent(inout) :: xm, xp, ym, yp, zm, zp
      select case (dir)
        case (MINUSX)
          call xm%add(p)
        case (PLUSX)
          call xp%add(p)
        case (MINUSY)
          call ym%add(p)
        case (PLUSY)
          call yp%add(p)
        case (MINUSZ)
          call zm%add(p)
        case (PLUSZ)
          call zp%add(p)
      end select
    end subroutine
  end subroutine




  subroutine MoveWMPointsToArray

    call MoveWMPointsToArrayPr
    call MoveWMPointsToArrayUVW
  end subroutine




  subroutine MoveWMPointsToArrayPr
    type(WMPoint),pointer :: p
    integer :: i

    allocate(WMPoints(WMPointsList%len()))

    if (size(WMPoints)>0) then

      call WMPointsList%iter_restart

      do i = 1, size(WMPoints)
        p => WMPointsList%iter_next()
        if (associated(p)) then
          WMPoints(i) = p
        else
          call error_stop("Assert error, pointer not associated. File "//__FILE__//" line ",__LINE__)
        end if
      end do

      call RemoveDuplicateWMPoints(WMPoints)

      call ComputeAreaFactor(WMPoints)

    end if

    call WMPointsList%finalize

  end subroutine MoveWMPointsToArrayPr



  subroutine MoveWMPointsToArrayUVW
    type(WMPointUVW),pointer :: p
    integer :: i

    call helper(UxmWMPsL, UxmWMpoints)
    call helper(UxpWMPsL, UxpWMpoints)
    call helper(UymWMPsL, UymWMpoints)
    call helper(UypWMPsL, UypWMpoints)
    call helper(UzmWMPsL, UzmWMpoints)
    call helper(UzpWMPsL, UzpWMpoints)

    call helper(VxmWMPsL, VxmWMpoints)
    call helper(VxpWMPsL, VxpWMpoints)
    call helper(VymWMPsL, VymWMpoints)
    call helper(VypWMPsL, VypWMpoints)
    call helper(VzmWMPsL, VzmWMpoints)
    call helper(VzpWMPsL, VzpWMpoints)

    call helper(WxmWMPsL, WxmWMpoints)
    call helper(WxpWMPsL, WxpWMpoints)
    call helper(WymWMPsL, WymWMpoints)
    call helper(WypWMPsL, WypWMpoints)
    call helper(WzmWMPsL, WzmWMpoints)
    call helper(WzpWMPsL, WzpWMpoints)

  contains
    subroutine helper(l, arr)
      type(WMPointUVWList), intent(inout) :: l
      type(WMPointUVW), allocatable, intent(out) :: arr(:)

      allocate(arr(l%len()))

      if (size(arr)>0) then

        call l%iter_restart

        do i = 1, size(arr)
          p => l%iter_next()
          if (associated(p)) then
            arr(i) = p
          else
            write(*,*) "Assert error, pointer not associated. File ",__FILE__," line ",__LINE__
            call error_stop
          end if
        end do

      end if

      call l%finalize
    end subroutine
  end subroutine MoveWMPointsToArrayUVW



  subroutine RemoveDuplicateWMPoints(WMPoints)
    use iso_c_binding
    use Sort
    type(WMPoint),allocatable,dimension(:),target,intent(inout)  :: WMPoints
    type(WMPoint),allocatable,dimension(:) :: TMP
    integer i,n

    !Choose the one closer to a wall. If of the same distance, choose the later one.
    !For wider compatibility we do not use MOLD= or SOURCE= in allocate.

    if (size(WMPoints)>1) then

      allocate(TMP(size(WMPoints)))

      call qsort(c_loc(WMPoints(1)), &
                 size(WMPoints,kind = c_size_t), &
                 storage_size(WMPoints,c_size_t)/storage_size(c_char_'a',c_size_t), &
  !                c_sizeof(WMPoints), &  F08+TS29113 seem to require C interoperable variable as argument.
                 c_funloc(CompareWMPoints))

      TMP(1) = WMPoints(1)
      n = 1
      do i = 2,size(WMPoints)
          if (WMPoints(i-1)%xi/=WMPoints(i)%xi .or. &
              WMPoints(i-1)%yj/=WMPoints(i)%yj .or. &
              WMPoints(i-1)%zk/=WMPoints(i)%zk) then !if the point is not duplicate of previous one

                n = n+1
                TMP(n) = WMPoints(i)

          end if
      end do

      WMPoints = TMP(1:n)

   end if

  end subroutine RemoveDuplicateWMPoints






  function CompareWMPoints(Aptr,Bptr) bind(C,name="CompareWMPoints") result(res)
    use iso_c_binding
    integer(c_int)         :: res
    type(c_ptr),value :: Aptr,Bptr
    type(WMPoint),pointer  :: A,B

    call c_f_pointer(Aptr,A)
    call c_f_pointer(Bptr,B)

    if ((A%xi+(A%yj-1)*Prnx+(A%zk-1)*Prnx*Prny) < (B%xi+(B%yj-1)*Prnx+(B%zk-1)*Prnx*Prny)) then
      res = -1_c_int
    else if ((A%xi+(A%yj-1)*Prnx+(A%zk-1)*Prnx*Prny) > (B%xi+(B%yj-1)*Prnx+(B%zk-1)*Prnx*Prny)) then
      res =  1_c_int
    else if (A%distx**2+A%disty**2+A%distz**2 < B%distx**2+B%disty**2+B%distz**2) then
      res = -1_c_int
    else if (A%distx**2+A%disty**2+A%distz**2 > B%distx**2+B%disty**2+B%distz**2) then
      res =  1_c_int
    else
      res =  0_c_int
    end if

  end function CompareWMPoints









 
  subroutine GetOutsideBoundariesWM(nscalars)
    integer, intent(in) :: nscalars

    call GetOutsideBoundariesWM_Pr(nscalars)

    call GetOutsideBoundariesWM_UVW
  end subroutine






  subroutine GetOutsideBoundariesWM_Pr(nscalars)
    integer, intent(in) :: nscalars
    integer       :: i,j,k
    type(WMPoint) :: WMP

    allocate(WMP%depscalar(nscalars))
    WMP%depscalar = 0
    WMP%evaporative_fraction = 0.1

    !type==0 below, therefore we need 
    ! type of free bounderies not to be -1
    if (Btype(We)==NOSLIP) then
      do k = 1,Prnz
       do j = 1,Prny
         if (Prtype(1,j,k)==0) then
           WMP%xi = 1
           WMP%yj = j
           WMP%zk = k
           WMP%distx = xU(0) - xPr(1)
           WMP%disty = 0
           WMP%distz = 0
           WMP%ustar = 1

           WMP%z0 = z0W
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(Ea)==NOSLIP) then
      do k = 1,Prnz
       do j = 1,Prny
         if (Prtype(Prnx,j,k)==0) then
           WMP%xi = Prnx
           WMP%yj = j
           WMP%zk = k
           WMP%distx = xU(Unx+1) - xPr(Prnx)
           WMP%disty = 0
           WMP%distz = 0
           WMP%ustar = 1

           WMP%z0 = z0E
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(So)==NOSLIP.or.(Btype(So)==DIRICHLET.and.sideU(2,So)==0)) then
      do k = 1,Prnz
       do i = 1,Prnx
         if (Prtype(i,1,k)==0) then
           WMP%xi = i
           WMP%yj = 1
           WMP%zk = k
           WMP%distx = 0
           WMP%disty = yV(0) - yPr(1)
           WMP%distz = 0
           WMP%ustar = 1

           if (Btype(So)==DIRICHLET) then
             WMP%wallu = sideU(1,So)
             WMP%wallv = 0
             WMP%wallw = sideU(3,So)
           end if

           WMP%z0 = z0S
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(No)==NOSLIP.or.(Btype(No)==DIRICHLET.and.sideU(2,No)==0)) then
      do k = 1,Prnz
       do i = 1,Prnx
         if (Prtype(i,Prny,k)==0) then
           WMP%xi = i
           WMP%yj = Prny
           WMP%zk = k
           WMP%distx = 0
           WMP%disty = yV(Vny+1) - yPr(Prny)
           WMP%distz = 0
           WMP%ustar = 1

           if (Btype(No)==DIRICHLET) then
             WMP%wallu = sideU(1,No)
             WMP%wallv = 0
             WMP%wallw = sideU(3,No)
           end if

           WMP%z0 = z0N
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

    if (Btype(Bo)==NOSLIP.or.(Btype(Bo)==DIRICHLET.and.sideU(3,Bo)==0)) then
      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,1)==0) then

           WMP%xi = i
           WMP%yj = j
           WMP%zk = 1
           WMP%distx = 0
           WMP%disty = 0
           WMP%distz = zW(0) - zPr(1)
           WMP%ustar = 1

           if (Btype(Bo)==DIRICHLET) then
             WMP%wallu = sideU(1,Bo)
             WMP%wallv = sideU(2,Bo)
             WMP%wallw = 0
           end if

           WMP%z0 = z0B

           if (TempBtype(Bo)==CONSTFLUX.or.TempBtype(Bo)==WALL_FLUX) then
             WMP%temperature_flux = sideTemp(Bo)
           end if

           if (TempBtype(Bo)==DIRICHLET.or.TempBtype(Bo)==WALL_DIRICHLET) then
             WMP%temperature = sideTemp(Bo)
           end if

           call AddWMPoint(WMP)

         end if
       end do
      end do
    end if

    if (Btype(To)==NOSLIP.or.(Btype(To)==DIRICHLET.and.sideU(3,To)==0)) then

      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,Prnz)==0) then
           WMP%xi = i
           WMP%yj = j
           WMP%zk = Prnz
           WMP%distx = 0
           WMP%disty = 0
           WMP%distz = zW(Wnz+1) - zPr(Prnz)
           WMP%ustar = 1

           if (Btype(To)==DIRICHLET) then
             WMP%wallu = sideU(1,To)
             WMP%wallv = sideU(2,To)
             WMP%wallw = 0
           end if

           WMP%z0 = z0T
           call AddWMPoint(WMP)
         end if
       end do
      end do
    end if

  end subroutine GetOutsideBoundariesWM_Pr


  subroutine GetOutsideBoundariesWM_UVW

    call helper(1, Unx, Uny, Unz, xU(-2:), yPr, zPr, Utype)

    call helper(2, Vnx, Vny, Vnz, xPr, yV(-2:), zPr, Vtype)

    call helper(3, Wnx, Wny, Wnz, xPr, yPr, zW(-2:), Wtype)

  contains

    subroutine helper(component, nx, ny, nz, x, y, z, Xtype)
      integer, intent(in) :: component, nx, ny, nz
      real(knd), intent(in) :: x(-2:), y(-2:), z(-2:)
      integer, intent(in) :: Xtype(-2:,-2:,-2:)
      integer       :: i,j,k
      type(WMPointUVW) :: p

      !type==0 below, therefore we need 
      ! type of free bounderies not to be -1
      if (Btype(We)==NOSLIP) then
        do k = 1,nz
         do j = 1,ny
           if (Xtype(1,j,k)==0.and.Xtype(0,j,k)<=0) then
             p%xi = 1
             p%yj = j
             p%zk = k
             p%distx = xU(0) - x(1)
             p%disty = 0
             p%distz = 0
             p%ustar = 1

             p%z0 = z0W
             call AddWMPointUVW(p, component, We)
           end if
         end do
        end do
      end if

      if (Btype(Ea)==NOSLIP) then
        do k = 1,nz
         do j = 1,ny
           if (Xtype(nx,j,k)==0.and.Xtype(nx+1,j,k)<=0) then
             p%xi = nx
             p%yj = j
             p%zk = k
             p%distx = xU(Unx+1) - x(nx)
             p%disty = 0
             p%distz = 0
             p%ustar = 1

             p%z0 = z0E
             call AddWMPointUVW(p, component, Ea)
           end if
         end do
        end do
      end if

      if (Btype(So)==NOSLIP.or.(Btype(So)==DIRICHLET.and.sideU(2,So)==0)) then
        do k = 1,nz
         do i = 1,nx
           if (Xtype(i,1,k)==0.and.Xtype(i,0,k)<=0) then
             p%xi = i
             p%yj = 1
             p%zk = k
             p%distx = 0
             p%disty = yV(0) - y(1)
             p%distz = 0
             p%ustar = 1

             if (Btype(So)==DIRICHLET) then
               p%wallu = sideU(1,So)
               p%wallv = 0
               p%wallw = sideU(3,So)
             end if

             p%z0 = z0S
             call AddWMPointUVW(p, component, So)
           end if
         end do
        end do
      end if

      if (Btype(No)==NOSLIP.or.(Btype(No)==DIRICHLET.and.sideU(2,No)==0)) then
        do k = 1,nz
         do i = 1,nx
           if (Xtype(i,ny,k)==0.and.Xtype(i,ny+1,k)<=0) then
             p%xi = i
             p%yj = ny
             p%zk = k
             p%distx = 0
             p%disty = yV(Vny+1) - y(ny)
             p%distz = 0
             p%ustar = 1

             if (Btype(No)==DIRICHLET) then
               p%wallu = sideU(1,No)
               p%wallv = 0
               p%wallw = sideU(3,No)
             end if

             p%z0 = z0N
             call AddWMPointUVW(p, component, No)
           end if
         end do
        end do
      end if

      if (Btype(Bo)==NOSLIP.or.(Btype(Bo)==DIRICHLET.and.sideU(3,Bo)==0)) then
        do j = 1,ny
         do i = 1,nx
           if (Xtype(i,j,1)==0.and.Xtype(i,j,0)<=0) then

             p%xi = i
             p%yj = j
             p%zk = 1
             p%distx = 0
             p%disty = 0
             p%distz = zW(0) - z(1)
             p%ustar = 1

             if (Btype(Bo)==DIRICHLET) then
               p%wallu = sideU(1,Bo)
               p%wallv = sideU(2,Bo)
               p%wallw = 0
             end if

             p%z0 = z0B

             call AddWMPointUVW(p, component, Bo)

           end if
         end do
        end do
      end if

      if (Btype(To)==NOSLIP.or.(Btype(To)==DIRICHLET.and.sideU(3,To)==0)) then

        do j = 1,ny
         do i = 1,nx
           if (Xtype(i,j,nz)==0.and.Xtype(i,j,nz+1)<=0) then
             p%xi = i
             p%yj = j
             p%zk = nz
             p%distx = 0
             p%disty = 0
             p%distz = zW(Wnz+1) - z(nz)
             p%ustar = 1

             if (Btype(To)==DIRICHLET) then
               p%wallu = sideU(1,To)
               p%wallv = sideU(2,To)
               p%wallw = 0
             end if

             p%z0 = z0T
             call AddWMPointUVW(p, component, To)
           end if
         end do
        end do
      end if

    end subroutine helper

  end subroutine GetOutsideBoundariesWM_UVW




  
  
  subroutine InitWMMasks
    integer :: i, j, k

    allocate(Scflx_mask(0:Prnx,Prny,Prnz))
    allocate(Scfly_mask(Prnx,0:Prny,Prnz))
    allocate(Scflz_mask(Prnx,Prny,0:Prnz))

    Scflx_mask = .true.
    Scfly_mask = .true.
    Scflz_mask = .true.

    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          if (Prtype(i,j,k)>0.or.Prtype(i+1,j,k)>0) Scflx_mask(i,j,k) = .false.
          if (Prtype(i,j,k)>0.or.Prtype(i,j+1,k)>0) Scfly_mask(i,j,k) = .false.
          if (Prtype(i,j,k)>0.or.Prtype(i,j,k+1)>0) Scflz_mask(i,j,k) = .false.
        end do
      end do
    end do

    call set_masks(Uflx_mask, Ufly_mask, Uflz_mask, &
                   UxmWMpoints, UxpWMpoints, UymWMpoints, UypWMpoints, UzmWMpoints, UzpWMpoints, Unx, Uny, Unz, Utype)

    call set_masks(Vflx_mask, Vfly_mask, Vflz_mask, &
                   VxmWMpoints, VxpWMpoints, VymWMpoints, VypWMpoints, VzmWMpoints, VzpWMpoints, Vnx, Vny, Vnz, Vtype)

    call set_masks(Wflx_mask, Wfly_mask, Wflz_mask, &
                   WxmWMpoints, WxpWMpoints, WymWMpoints, WypWMpoints, WzmWMpoints, WzpWMpoints, Wnx, Wny, Wnz, Wtype)

  contains
  
    subroutine set_masks(flx, fly, flz, xm, xp, ym, yp, zm, zp, nx, ny, nz, Xtype)
      logical(slog), dimension(:,:,:), allocatable, intent(out) :: flx, fly, flz
      type(WMpointUVW),  dimension(:), intent(in) :: xm, xp, ym, yp, zm, zp
      integer, intent(in) :: nx, ny, nz
      integer, intent(in) :: Xtype(-2:,-2:,-2:)
      integer :: i, j, k

      allocate(flx(nx+1,ny,nz))
      allocate(fly(nx,ny+1,nz))
      allocate(flz(nx,ny,nz+1))
  
      flx = .true.
      fly = .true.
      flz = .true.


      do i = 1, size(xm)
        associate(p => xm(i))
            flx(p%xi,p%yj,p%zk) = .false.
        end associate
      end do
      do i = 1, size(xp)
        associate(p => xp(i))
            flx(p%xi+1,p%yj,p%zk) = .false.
        end associate
      end do
      do i = 1, size(ym)
        associate(p => ym(i))
            fly(p%xi,p%yj,p%zk) = .false.
        end associate
      end do
      do i = 1, size(yp)
        associate(p => yp(i))
            fly(p%xi,p%yj+1,p%zk) = .false.
        end associate
      end do
      do i = 1, size(zm)
        associate(p => zm(i))
            flz(p%xi,p%yj,p%zk) = .false.
        end associate
      end do
      do i = 1, size(zp)
        associate(p => zp(i))
            flz(p%xi,p%yj,p%zk+1) = .false.
        end associate
      end do

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if (Xtype(i,j,k)>0.and.Xtype(i-1,j,k)>0) flx(i,j,k) = .false.
            if (Xtype(i,j,k)>0.and.Xtype(i,j-1,k)>0) fly(i,j,k) = .false.
            if (Xtype(i,j,k)>0.and.Xtype(i,j,k-1)>0) flz(i,j,k) = .false.
          end do
        end do
      end do
    end subroutine

  end subroutine InitWMMasks



  subroutine DivergenceWM(U, V, W)
    real(knd), contiguous, intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer :: i, xi, yj, zk
    real(knd) :: div

    !$omp parallel do private(i, xi, yj, zk, div)
    do i=1,size(WMPoints)
      xi = WMPoints(i)%xi
      yj = WMPoints(i)%yj
      zk = WMPoints(i)%zk

      div = 0
      if (Prtype(xi+1,yj,zk)<=0) div = div + U(xi,yj,zk) / dxmin
      if (Prtype(xi-1,yj,zk)<=0) div = div - U(xi-1,yj,zk) / dxmin
      if (Prtype(xi,yj+1,zk)<=0) div = div + V(xi,yj,zk) / dymin
      if (Prtype(xi,yj-1,zk)<=0) div = div - V(xi,yj-1,zk) / dymin
      if (Prtype(xi,yj,zk+1)<=0) div = div + W(xi,yj,zk) / dzmin
      if (Prtype(xi,yj,zk-1)<=0) div = div - W(xi,yj,zk-1) / dzmin

      WMPoints(i)%div = div
    end do
    !$omp end parallel do

  end subroutine





















  pure subroutine WMFlatUstar(ustar,vel,dist)
    real(knd),intent(inout) :: ustar
    real(knd),intent(in) :: vel,dist
    real(knd),parameter :: eps = 1e-4_knd
    real(knd),parameter :: yplcrit = 11.225_knd
    real(knd) :: ustar2,ustar_lam
    integer i

    ustar_lam = sqrt(vel/(dist*Re))

    if ((dist*ustar_lam*Re)<yplcrit) then

      ustar = ustar_lam

    else   !turbulent region
      i = 0

      do
        i = i+1
        ustar2 = ustar

        if ((dist*ustar2*Re)<yplcrit) then
          ustar = sqrt(vel/(dist*Re))
        else
          ustar = vel/(log(abs(ustar2*dist*Re))/0.41_knd+5.2_knd)
        end if

        if  (abs(ustar-ustar2)/abs(ustar)<eps) exit

        if (i>=50) then
                    ustar = 0
                    exit
        end if

      end do

    end if

  end subroutine WMFlatUstar


  subroutine WMFlatStress(ustar,distvect,uvect,walluvect,tan_vel,distance,tan_vect)
    real(knd), intent(inout) :: ustar
    real(knd), intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd), optional, intent(out)   :: tan_vel, distance, tan_vect(3)
    real(knd) :: vel, dist, tan_vect_loc(3)

    call vel_and_dist(tan_vect_loc, dist, uvect, walluvect, distvect)

    vel = norm2(tan_vect_loc)

    if (vel/=0) then

      call WMFlatUstar(ustar,vel,dist)

    end if

    if (ustar<0) ustar = 0

    if (present(tan_vel)) tan_vel = vel
    if (present(distance)) distance = dist
    if (present(tan_vect)) tan_vect = tan_vect_loc

  end subroutine WMFlatStress


  subroutine WMFlatVisc(visc,ustar,distvect,uvect,walluvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd) :: vel, dist

    call WMFlatStress(ustar,distvect,uvect,walluvect,vel,dist)

    if (vel>0) then

     if (dist*ustar*Re>1) then
       visc = ustar**2 * dist/vel
     else if (Re>0) then
       visc = 1._knd/Re
     else
       visc = 0
     end if

    else if (Re>0) then

      visc = 1._knd/Re

    else

      visc = 0

    end if

  end subroutine WMFlatVisc




  pure subroutine WMRoughUstar(ustar,vel,dist,z0)
    real(knd),intent(inout) :: ustar
    real(knd),intent(in) :: vel,dist,z0
    real(knd),parameter  :: eps = 1e-4_knd
    real(knd),parameter  :: yplcrit = 11.225_knd
    real(knd) :: dist_plus

    dist_plus = sqrt(dist * Re * vel)

    !under z0 the whole concept of rougness parameter breaks down
    !if in the laminar region for flat boundary layer also treat as flat and possibly laminar
    if (dist<=z0 .or. dist_plus<yplcrit) then
      if (Re>0) then
        call WMFlatUstar(ustar,vel,dist)
      else
        ustar = 0
      end if
    else 
      ustar = vel * 0.41_knd / log(dist/z0)
    end if

  end subroutine WMRoughUstar


  pure subroutine WMRoughStress(ustar,z0,distvect,uvect,walluvect,tan_vel,distance,tan_vect)
    real(knd), intent(inout) :: ustar
    real(knd), intent(in)    :: z0
    real(knd), intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd), optional, intent(out)   :: tan_vel, distance, tan_vect(3)
    real(knd) :: vel, dist, tan_vect_loc(3)

    call vel_and_dist(tan_vect_loc, dist, uvect, walluvect, distvect)

    vel = norm2(tan_vect_loc)

    if (vel/=0) then

      call WMRoughUstar(ustar,vel,dist,z0)

    end if

    if (ustar<0) ustar = 0

    if (present(tan_vel)) tan_vel = vel
    if (present(distance)) distance = dist
    if (present(tan_vect)) tan_vect = tan_vect_loc

  end subroutine WMRoughStress


  pure subroutine WMRoughVisc(visc,ustar,z0,distvect,uvect,walluvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0
    real(knd),intent(in)    :: distvect(3), uvect(3), walluvect(3)
    real(knd) vect(3), vel, dist

    call WMRoughStress(ustar,z0,distvect,uvect,walluvect,vel,dist)

    if (vel>0 .and. ustar**2 * dist/vel>1._knd/Re) then
      visc = ustar**2 * dist/vel
    else if (Re>0) then
      visc = 1._knd/Re
    else
      visc = 0
    end if

  end subroutine WMRoughVisc


  pure subroutine vel_and_dist(tan_vect, dist, uvect, walluvect, distvect)
    real(knd), intent(out) :: tan_vect(3),  dist
    real(knd), intent(in)  :: uvect(3), walluvect(3), distvect(3)

    dist = norm2(distvect)

    tan_vect = uvect - walluvect

    tan_vect = tan_vect - dot_product(tan_vect,distvect) * distvect / dist**2  !tangential part
  end subroutine





  pure subroutine WMFlatPrGradUstar(ustar,vel,prgrad,dist)
    real(knd),intent(out) :: ustar
    real(knd),intent(in) :: vel,prgrad,dist
    real(knd),parameter :: eps = 1e-3_knd
    real(knd),parameter :: yplcrit = 11.225_knd
    real(knd),parameter :: ustar_div = 1000
    real(knd) :: ustar1,ustar2,ustar_lam
    integer i
    integer,parameter :: maxiter = 30

    !u/u_* = z * u_* / nu  +  dp/dx * z**2 / (2 * u_*)
    ustar_lam = sqrt(abs((1._knd/(dist*Re)) * (vel - dist**2 * prgrad/2) ))

    if ((dist*ustar_lam*Re)<yplcrit) then

      ustar = ustar_lam

    else   !turbulent region


      i = 0

      ustar1 = ustar_lam

      do
        i = i+1

        ustar2 = newguess(ustar1)
        if (ustar2 < 0) then
          if (ustar1>max(100*vel,10*Uinlet).or.i>=30) then
            ustar = ustar_lam
            exit
          end if
          ustar1 = ustar1 * 10
        else if (ustar2<tiny(1._knd)) then
          ustar = ustar_lam
          exit
        else if (abs(ustar1-ustar2)/abs(ustar1)<eps) then
          ustar = ustar2
          exit
        else if (i<20) then
          ustar1 = ustar2
        else
          if (abs(ustar1-ustar2)/abs(ustar1)<0.1) then
            ustar = ustar2
          else
            ustar = ustar_lam
          end if
          exit
        end if
      end do

    end if !laminar/turbulent

    contains

      pure function newguess(ustar)
         !linearize the function by letting the ustar in log constant
         !  and solve the quadratic equation for the larger root
         real(knd) newguess
         real(knd),intent(in) :: ustar
         real(knd) a,b,c,D
         !u/u_* = dp/dx * z / (k (u_*)**2)  + (1/k) * ln(z * u_* / nu) + B

         a = log(ustar*dist*Re)/0.41_knd + 5.2_knd
         b = - vel
         c = prgrad * dist / 0.41
         !function to find root of is f(ustar) = a*ustar**2 + b*ustar + c
         D = b**2 - 4*a*c

         if (D<0) then !solution does not exist
           newguess = -1
         else
           newguess = (-b+sqrt(D))/(2*a)
         end if
      end function

  end subroutine WMFlatPrGradUstar




  pure subroutine WMFlatPrGradVisc(visc,ustar,distvect,uvect,walluvect,prgradvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: distvect(3),uvect(3),walluvect(3),prgradvect(3)
    real(knd) vect(3),vel,dist,prgrad

    dist = sqrt(sum(distvect**2))

    vect = uvect - walluvect

    vect = vect - dot_product(vect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    if (vel>=tiny(1._knd)) then
      !in the same direction as tangential velocity vector
      prgrad = dot_product( prgradvect , vect ) / vel

    else
      !tangential to the wall
      prgrad = sqrt(sum((prgradvect - dot_product(prgradvect,distvect) * distvect / dist**2)**2))

    end if

    call WMFlatPrGradUstar(ustar,vel,prgrad,dist)

    if (ustar<0) ustar = 0

    if (vel>0) then

     if (dist*ustar*Re>1) then
       visc = ustar**2 * dist/vel
     else if (Re>0) then
       visc = 1._knd/Re
     else
       visc = 0
     end if

    else if (Re>0) then

      visc = 1._knd/Re

    else

      visc = 0

    end if

  end subroutine WMFlatPrGradVisc



















  pure real(knd) function PsiM_MO(zeta)
    real(knd),intent(in):: zeta
    real(knd) x

    if (zeta<0) then
     x = (1-15._knd*zeta)**(1/4._knd)
     PsiM_MO = log(((1+x**2)/2._knd)*((1+x)/2._knd)**2)-2._knd*atan(x)+pi/2
    else
     PsiM_MO = -4.8_knd*zeta !GABLS recommend ation
    end if
  end function PsiM_MO


  pure real(knd) function PsiH_MO(zeta)
    real(knd),intent(in):: zeta
    real(knd) x

    if (zeta<0) then
     x = (1-15._knd*zeta)**(1/4._knd)
     PsiH_MO = 2._knd*log((1+x**2)/2._knd)
    else
     PsiH_MO = -7.8_knd*zeta !GABLS recommend ation
    end if
  end function PsiH_MO


  pure real(knd) function Obukhov_zL(ustar,temperature_flux,tempref,g,z)
    real(knd),intent(in):: ustar,temperature_flux,tempref,g,z

    Obukhov_zL = z*(0.4_knd*(g/tempref)*temperature_flux)/(-ustar**3)
  end function Obukhov_zL



  pure subroutine WM_MO_FLUX_ustar(vel,dist,ustar,z0,temperature_flux,Re,temperature_ref,grav_acc)
    implicit none

    real(knd),intent(inout) :: ustar
    real(knd),parameter  :: eps = 1e-3
    real(knd),parameter  :: yplcrit = 11.225_knd
    real(knd),intent(in) :: vel,dist,z0,temperature_flux
    real(knd),intent(in) :: Re,temperature_ref,grav_acc
    real(knd) ustar2,zL,zL2,Psi
    integer i

    if (dist<=z0) then

     if (Re>0) then

      if ((dist*ustar*Re)<yplcrit) then
        ustar = sqrt(vel/(dist*Re))
      else
        ustar = vel/(log(abs(ustar*dist*Re))/0.4_knd+5.2_knd)
      end if
     else
       ustar = 0
     end if

    else

     i = 0
     zL = 0
     Psi = 0

     do
      i = i+1
      ustar2 = ustar

      ustar = ustar+(max(vel*0.4_knd/(log(max((dist/z0)-Psi,1E-5))),0._knd)-ustar)/2

      if (ustar<1E-4) then
       zL = -10000
      else
       zL2 = Obukhov_zL(ustar,temperature_flux,temperature_ref,grav_acc,dist)
       zL = zL+(zL2-zL)/2
      end if

      Psi = PsiM_MO(zL)

      if  (abs(ustar-ustar2)/max(abs(ustar),1.e-3_knd)<eps) exit

      if (i>=50) then
                  ustar = 0
                  exit
      end if

     end do

    end if

  end subroutine WM_MO_FLUX_ustar



  pure subroutine WM_MO_DIRICHLET_ustar_tfl(ustar,temperature_flux,vel,dist,z0,tempdif)
    real(knd),intent(inout) :: ustar,temperature_flux
    real(knd),intent(in) :: vel,dist,z0,tempdif
    real(knd),parameter :: eps = 1e-3
    real(knd),parameter :: yplcrit = 11.225
    real(knd) :: zL,zL0,Rib
    integer i

    call WMRoughUstar(ustar,vel,dist,z0)

    if (dist<=z0) then

      if (Re>0) then
       if ((dist*ustar*Re)<yplcrit) then
         ustar = sqrt(vel/(dist*Re))
       else
         ustar = vel/(log(abs(ustar*dist*Re))/0.4_knd+5.2_knd)
       end if
      else
       ustar = 0
      end if

    else

      Rib = -grav_acc*dist*tempdif/(temperature_ref*max(vel**2,1E-6_knd))
      zL = 0
      i = 0
      if (Rib>0.34_knd) then
                          ustar = 0
                          temperature_flux = 0
      else
        do
          i = i+1
          zL0 = zL
          zL = Rib*(log(dist/z0)-PsiM_MO(zl))**2/(log(dist/z0)-PsiH_MO(zl))
          if  (abs(zL-zL0)/max(abs(zL),1.e-3_knd)<eps) exit
          if (i>=50.or.zL>100) exit
        end do

        if (i>=50.or.zL>100) then
          ustar = 0
          temperature_flux = 0
        else
          ustar = vel*0.4_knd/(log(dist/z0)-PsiM_MO(zL))
          temperature_flux = 0.4_knd*ustar*tempdif/(log(dist/z0)-PsiH_MO(zL))
        end if
      end if
    end if
  end subroutine WM_MO_DIRICHLET_ustar_tfl







  pure subroutine WM_MO_FLUX(visc,ustar,temperature_flux,z0,distvect,uvect)
    real(knd),intent(out) :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0,temperature_flux
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - dot_product(uvect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    if (vel/=0) then
      call WM_MO_FLUX_ustar(vel,dist,ustar,z0,temperature_flux,Re,temperature_ref,grav_acc)
      if (ustar<0) ustar = 0
    end if

    if (vel>0.and.ustar*ustar*dist/vel>1._knd/Re) then
      visc = ustar*ustar*dist/vel
    else if (Re>0) then
      visc = 1._knd/Re
    else
      visc = 0
    end if

  end subroutine WM_MO_FLUX



  pure subroutine WM_MO_DIRICHLET(visc,ustar,temperature_flux,z0,tempdif,distvect,uvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar,temperature_flux
    real(knd),intent(in)    :: z0
    real(knd),intent(in)    :: tempdif ! temperature difference surface - nearest point
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - dot_product(uvect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    call WM_MO_DIRICHLET_ustar_tfl(ustar,temperature_flux,vel,dist,z0,tempdif)

    if (ustar<0) ustar = 0

    if ((vel>0.or.temperature_flux/=0)) then
      visc = max(ustar*ustar*dist/vel,1._knd/Re)
    else if (Re>0) then
      visc = 1._knd/Re
    else
      visc = 0
    end if

  end subroutine WM_MO_DIRICHLET



  pure subroutine BOUND_temperature_flux(Nu)
    real(knd),intent(inout):: Nu(-1:,-1:)
    integer i,j,nx,ny

    nx = Prnx
    ny = Prny

    if (Btype(Ea)==PERIODIC) then
      do j = 1,ny
        Nu(0,j) = Nu(nx,j)
        Nu(nx+1,j) = Nu(1,j)
      end do
    else
      do j = 1,ny
        Nu(0,j) = Nu(1,j)
        Nu(nx+1,j) = Nu(nx,j)
      end do
    end if

    if (Btype(No)==PERIODIC) then
      do i = 1,nx
        Nu(i,0) = Nu(i,ny)
        Nu(i,ny+1) = Nu(i,1)
      end do
    else
      do i = 1,nx
        Nu(i,0) = Nu(i,1)
        Nu(i,ny+1) = Nu(i,ny)
      end do
    end if
  end subroutine BOUND_temperature_flux





  subroutine InitTempFL(Temperature)
    real(knd),intent(in) :: Temperature(-1:,-1:,-1:)
    integer i

    if (enable_buoyancy.and.TempBtype(Bo)==DIRICHLET) then
      do i = 1,size(WMPoints)
        if (WMPoints(i)%zk==1) WMPoints(i)%temperature_flux = -TDiff(WMPoints(i)%xi,WMPoints(i)%yj,1)*&
                       (temperature(WMPoints(i)%xi,WMPoints(i)%yj,1) - temperature(WMPoints(i)%xi,WMPoints(i)%yj,0))
      end do
    end if

  end subroutine InitTempFL


  pure recursive subroutine WallPrGradient(prgrad,i,j,k,Pr,Prtype)
    real(knd),intent(out) :: prgrad(3)
    integer,intent(in)    :: i,j,k
    real(knd),intent(in)  :: Pr(1:,1:,1:)
    integer,intent(in)    :: Prtype(0:,0:,0:)
    integer n

    prgrad = 0

    n=0
    if (Prtype(i+1,j,k)>0 .and. i<Prnx) then
      prgrad(1) = prgrad(1) + (Pr(i+1,j,k) - Pr(i,j,k))/(dxU(i))
      n = n + 1
    end if
    if (Prtype(i-1,j,k)>0 .and. i>1) then
      prgrad(1) = prgrad(1) + (Pr(i,j,k) - Pr(i-1,j,k))/(dxU(i-1))
      n = n + 1
    end if
    if (n>0) prgrad(1) = prgrad(1)/n

    n=0
    if (Prtype(i,j+1,k)>0 .and. j<Prny) then
      prgrad(2) = prgrad(2) + (Pr(i,j+1,k) - Pr(i,j,k))/(dyV(j))
      n = n + 1
    end if
    if (Prtype(i,j-1,k)>0 .and. j>1) then
      prgrad(2) = prgrad(2) + (Pr(i,j,k) - Pr(i,j-1,k))/(dyV(j-1))
      n = n + 1
    end if
    if (n>0) prgrad(2) = prgrad(2)/n

    n=0
    if (Prtype(i,j,k+1)>0 .and. k<Prnz) then
      prgrad(3) = prgrad(3) + (Pr(i,j,k+1) - Pr(i,j,k))/(dzW(k))
      n = n + 1
    end if
    if (Prtype(i,j,k-1)>0 .and. k>1) then
      prgrad(3) = prgrad(3) + (Pr(i,j,k) - Pr(i,j,k-1))/(dzW(k-1))
      n = n + 1
    end if
    if (n>0) prgrad(3) = prgrad(3)/n

    prgrad = prgrad + [prgradientx,prgradienty,0._knd]

  end subroutine WallPrGradient



  subroutine ComputeViscsWM(U,V,W,Pr,Temperature)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(knd),dimension(1:,1:,1:),   intent(in) :: Pr
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Temperature
    integer i,j,xi,yj,zk
    real(knd) tdif
    real(knd) dist(3), vel(3), wallvel(3), prgrad(3)
    real(knd) visc


    !$omp parallel do private(i,xi,yj,zk,tdif, vel, wallvel, dist, prgrad, visc)
    do i = 1,size(WMPoints)

      xi = WMPoints(i)%xi
      yj = WMPoints(i)%yj
      zk = WMPoints(i)%zk

      dist(:) = [WMPoints(i)%distx, WMPoints(i)%disty, WMPoints(i)%distz]

      if (all(dist==0))then
        write(*,*) "ijk",xi,yj,zk
        write(*,*) "dist",dist
        call error_stop("Error, WM point can not be exactly on the wall!")
      end if

      vel(1) = (U(xi,yj,zk)+U(xi-1,yj,zk))/2._knd
      vel(2) = (V(xi,yj,zk)+V(xi,yj-1,zk))/2._knd
      vel(3) = (W(xi,yj,zk)+W(xi,yj,zk-1))/2._knd

      wallvel(:) = [WMPoints(i)%wallu, WMPoints(i)%wallv, WMPoints(i)%wallw]


      if (WMPoints(i)%z0>0) then

         if (enable_buoyancy .and. WMPoints(i)%temperature>0) then

           tdif = WMPoints(i)%temperature - Temperature(xi,yj,zk)

           call WM_MO_DIRICHLET(visc, WMPoints(i)%ustar, &
                                WMPoints(i)%temperature_flux, &
                                WMPoints(i)%z0, tdif, dist, vel)

         else if (enable_buoyancy .and. WMPoints(i)%temperature_flux>0) then

           call WM_MO_FLUX(visc, WMPoints(i)%ustar, WMPoints(i)%temperature_flux, &
                           WMPoints(i)%z0, dist, vel)

         else

           call WMRoughVisc(visc, &
                            WMPoints(i)%ustar, WMPoints(i)%z0, &
                            dist, vel, wallvel)
         end if

       else

         if (Re<=0) then
           call error_stop("The wall model requires positive viscosity or roughness length.")
         end if

         if (wallmodeltype == 2) then

           call WallPrGradient(prgrad,xi,yj,zk,Pr,Prtype)

           call WMFlatPrGradVisc(visc, &
                           WMPoints(i)%ustar, &
                           dist, vel, wallvel, prgrad)

         else

           call WMFlatVisc(visc, &
                           WMPoints(i)%ustar, &
                           dist, vel, wallvel)

         end if
       end if

       Viscosity(xi,yj,zk) = visc

    end do
    !$omp end parallel do

    if (enable_buoyancy.and. TempBtype(Bo)==DIRICHLET) call Bound_temperature_flux(BsideTFLArr)


  end subroutine ComputeViscsWM






  subroutine ComputeUVWFluxesWM(U,V,W,Pr,Temperature)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(knd),dimension(1:,1:,1:),   intent(in) :: Pr
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Temperature




    if (enable_buoyancy.and. TempBtype(Bo)==DIRICHLET) call Bound_temperature_flux(BsideTFLArr)

    call fluxes(UxmWMpoints, 1, MINUSX)
    call fluxes(UxpWMpoints, 1, PLUSX)
    call fluxes(UymWMpoints, 1, MINUSY)
    call fluxes(UypWMpoints, 1, PLUSY)
    call fluxes(UzmWMpoints, 1, MINUSZ)
    call fluxes(UzpWMpoints, 1, PLUSZ)

    call fluxes(VxmWMpoints, 2, MINUSX)
    call fluxes(VxpWMpoints, 2, PLUSX)
    call fluxes(VymWMpoints, 2, MINUSY)
    call fluxes(VypWMpoints, 2, PLUSY)
    call fluxes(VzmWMpoints, 2, MINUSZ)
    call fluxes(VzpWMpoints, 2, PLUSZ)

    call fluxes(WxmWMpoints, 3, MINUSX)
    call fluxes(WxpWMpoints, 3, PLUSX)
    call fluxes(WymWMpoints, 3, MINUSY)
    call fluxes(WypWMpoints, 3, PLUSY)
    call fluxes(WzmWMpoints, 3, MINUSZ)
    call fluxes(WzpWMpoints, 3, PLUSZ)


  contains

    subroutine fluxes(points, component, direction)
      type(WMPointUVW), intent(inout), target :: points(:)
      integer, intent(in) :: component, direction
      integer i,j,xi,yj,zk
      real(knd) dist(3), vel(3), wallvel(3), tan_vect(3), mag
      type(WMPointUVW), pointer :: p
      real(knd), parameter :: eps = 0.0001_knd

      !$omp parallel do private(i,xi,yj,zk,dist,vel,wallvel,tan_vect,p,mag)
      do i = 1,size(points)

!         associate (p => points(i))  NOTE: ASSOCIATE supported only in OpenMP 4
          p => points(i)
          xi = p%xi
          yj = p%yj
          zk = p%zk

          dist = [p%distx, p%disty, p%distz]

          if (all(dist==0))then
            write(*,*) "ijk",xi,yj,zk
            write(*,*) "dist",dist
            call error_stop("Error, WM point can not be exactly on the wall!")
          end if

          vel = local_velocity(U,V,W,component,xi,yj,zk)

          wallvel = [p%wallu, p%wallv, p%wallw]


          if (p%z0>0) then
!TODO: The temperature flux is always 0 now in momentum points, so no stability effect here!
             if (enable_buoyancy .and. p%temperature>0) then

               call error_stop("Not implemented!")

             else if (enable_buoyancy .and. p%temperature_flux>0) then

               call error_stop("Not implemented!")

             else

               call WMRoughStress(p%ustar, p%z0, &
                                  dist, vel, wallvel, tan_vect = tan_vect)
             end if

          else

             if (Re<=0) then
               call error_stop("The wall model requires positive viscosity or roughness length.")
             end if


             call WMFlatStress(p%ustar, &
                               dist, vel, wallvel, tan_vect = tan_vect)
           end if

           mag = norm2(tan_vect)

           if (mag>eps) then

             p%flux = tan_vect(component) / mag * p%ustar**2

             if (mod(direction,2)==1) p%flux = - p%flux

           else

             p%flux = 0

           end if

!         end associate

      end do
      !$omp end parallel do

    end subroutine

  end subroutine ComputeUVWFluxesWM


  pure function local_velocity(U,V,W,component,xi,yj,zk) result(vel)
    real(knd) :: vel(3)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    integer, intent(in) :: component, xi, yj, zk

    select case (component)
      case (1)
        vel(1) =   U(xi, yj, zk)
        vel(2) = ( V(xi-1,yj,  zk) + &
                   V(xi,  yj,  zk) + &
                   V(xi-1,yj+1,zk) + &
                   V(xi,  yj+1,zk) &
                 ) / 4._knd
        vel(3) = ( W(xi-1,yj, zk  ) + &
                   W(xi,  yj, zk  ) + &
                   W(xi-1,yj, zk+1) + &
                   W(xi,  yj, zk+1) &
                 ) / 4._knd
      case (2)
        vel(1) = ( U(xi  ,yj-1,zk) + &
                   U(xi,  yj,  zk) + &
                   U(xi+1,yj-1,zk) + &
                   U(xi+1,yj,  zk) &
                 ) / 4._knd
        vel(2) =   V(xi, yj, zk)
        vel(3) = ( W(xi, yj-1,zk  ) + &
                   W(xi, yj,  zk  ) + &
                   W(xi, yj-1,zk+1) + &
                   W(xi, yj,  zk+1) &
                 ) / 4._knd
      case (3)
        vel(1) = ( U(xi  ,yj, zk-1) + &
                   U(xi,  yj, zk  ) + &
                   U(xi+1,yj, zk-1) + &
                   U(xi+1,yj, zk  ) &
                 ) / 4._knd
        vel(2) = ( V(xi, yj,  zk-1) + &
                   V(xi, yj,  zk  ) + &
                   V(xi, yj+1,zk-1) + &
                   V(xi, yj+1,zk  ) &
                 ) / 4._knd
        vel(3) =   W(xi, yj, zk)
    end select
  end function


  pure real(knd) function GroundUstar()
    if (any(WMPoints%zk == 1)) then
      GroundUstar = sum(WMPoints%ustar, mask = (WMPoints%zk == 1)) / count(WMPoints%zk == 1)
    else
      GroundUstar = 0
    end if
  end function GroundUstar




  pure real(knd) function TotalUstar()
    if (size(WMPoints) > 0) then
      TotalUstar = sum(WMPoints%ustar) / size(WMPoints%zk)
    else
      TotalUstar = 0
    end if
  end function TotalUstar





  pure function GroundDeposition() result(depos)
    real(knd), dimension(:) :: depos(1:Prnx,1:Prny,num_of_scalars)

    integer :: i, j

    depos = 0

    do j = 1, size(WMPoints)

      if (allocated(WMPoints(j)%depscalar)) then

        do i = 1, num_of_scalars
          depos(WMPoints(j)%xi,WMPoints(j)%yj,i) = depos(WMPoints(j)%xi,WMPoints(j)%yj,i) + &
                                                   WMPoints(j)%depscalar(i)
        end do

      end if

    end do
  end function GroundDeposition



  pure real(knd) function SurfTemperature(x,y,t)
   real(knd),intent(in):: x,y
   real(TIM),intent(in):: t

   SurfTemperature = sideTemp(Bo)  ! Needs bet to allow time evolution somehow
  end function


 end module Wallmodels