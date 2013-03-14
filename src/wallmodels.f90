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



module Wallmodels
use PARAMETERS


implicit none

 private
 public WMPoint, AddWMPoint, FirstWMPoint, ComputeViscsWM, MoveWMPointsToArray, GetOutsideBoundariesWM,&
        InitTempFl, GroundDeposition, GroundUstar, WMPoints, ListLength, wallmodeltype
#ifdef __HMPP
 public hmppWMpoint,WMPtoHMPP
#endif

 type WMpoint   !points in which we apply wall model

   integer   :: xi
   integer   :: yj
   integer   :: zk

   real(knd) :: distx
   real(knd) :: disty
   real(knd) :: distz

   real(knd) :: z0 = 0
   real(knd) :: ustar = 1
   real(knd) :: temp = 0
   real(knd) :: tempfl = 0
   real(knd) :: wallu = 0
   real(knd) :: wallv = 0
   real(knd) :: wallw = 0

   real(knd),allocatable:: depscalar(:)

   type(WMpoint),pointer:: next=>null()

 end type WMpoint

 type(WMpoint),pointer::FirstWMPoint=>null(), LastWMPoint=>null()

 type(WMPoint),dimension(:),allocatable :: WMPoints

 interface ListLength
   module procedure WMPoint_ListLength
 end interface

 interface assignment(=)
   module procedure WMPtoWMP
 end interface

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
    ToWMP%temp   = FromWMP%temp
    ToWMP%tempfl = FromWMP%tempfl

  end subroutine WMPtoHMPP
#endif

  subroutine WMPtoWMP(ToWMP,FromWMP)
    type(WMPoint),intent(out) :: ToWMP
    type(WMPoint),intent(in)  :: FromWMP

    ToWMP%xi = FromWMP%xi
    ToWMP%yj = FromWMP%yj
    ToWMP%zk = FromWMP%zk

    ToWMP%distx = FromWMP%distx
    ToWMP%disty = FromWMP%disty
    ToWMP%distz = FromWMP%distz

    ToWMP%z0     = FromWMP%z0
    ToWMP%ustar  = FromWMP%ustar
    ToWMP%temp   = FromWMP%temp
    ToWMP%tempfl = FromWMP%tempfl
    ToWMP%wallu  = FromWMP%wallu
    ToWMP%wallv  = FromWMP%wallv
    ToWMP%wallw  = FromWMP%wallw

    allocate(ToWMP%depscalar(size(FromWMP%depscalar)))

    ToWMP%depscalar = FromWMP%depscalar

  end subroutine WMPtoWMP


  subroutine AddWMpoint(WMP)
    type(WMpoint),intent(in):: WMP

    if (.not.associated(LastWMPoint)) then

     allocate(FirstWMPoint)

     if (computedeposition>0) then
       allocate( FirstWMPoint%depscalar(size(WMP%depscalar)) )
       FirstWMPoint%depscalar = 0
     end if

     FirstWMPoint = WMP
     LastWMPoint=>FirstWMPoint

    else

     allocate(LastWMPoint%next)

     if (computedeposition>0) then
       allocate( LastWMPoint%next%depscalar(size(WMP%depscalar)) )
       LastWMPoint%next%depscalar = 0
     end if

     LastWMPoint%next = WMP
     LastWMPoint=>LastWMPoint%next

    end if
  end subroutine AddWMPoint


  subroutine WMPoint_DeallocateList(WMP)
    type(WMPoint),pointer,intent(inout) :: WMP
    type(WMPoint),pointer :: Aux,Aux2

    if (.not.associated(WMP)) return

    Aux => WMP

    do
      if (allocated(Aux%depscalar)) deallocate(Aux%depscalar)
      if (associated(Aux%next)) then
        Aux => Aux%next
      else
        exit
      end if
    end do


    Aux =>WMP

    do

      if (associated(Aux%next)) then
        Aux2 => Aux%next
      else
        Aux2 => null()
      end if

      deallocate(Aux)

      if (associated(Aux2)) then
        Aux => Aux2
      else
        exit
      end if

    end do

    WMP => null()

  end subroutine WMPoint_DeallocateList



  subroutine MoveWMPointsToArray
    type(WMPoint),pointer :: CurrentWMPoint
    integer :: i

    allocate(WMPoints(WMPoint_ListLength(FirstWMPoint)))

    if (size(WMPoints)>0) then

      CurrentWMPoint => FirstWMPoint
      i = 0
      do
       if (associated(CurrentWMPoint)) then
         i = i + 1
         WMPoints(i) = CurrentWMPoint
       else
         exit
       end if
       CurrentWMPoint => CurrentWMPoint%next
      end do

      call RemoveDuplicateWMPoints(WMPoints)

    end if

    call WMPoint_DeallocateList(FirstWMPoint)

  end subroutine MoveWMPointsToArray



  function WMPoint_ListLength(WMP) result(nWMP)
    integer :: nWMP
    type(WMPoint),pointer :: WMP
    type(WMPoint),pointer :: CurrentWMPoint

    nWMP = 0

    if (associated(WMP)) then

      CurrentWMPoint => WMP

      do
        nWMP = nWMP + 1

        if (associated(CurrentWMPoint%next)) then
          CurrentWMPoint => CurrentWMPoint%next
        else
          exit
        end if
      end do

    end if

  end function WMPoint_ListLength


  subroutine RemoveDuplicateWMPoints(WMPoints)
    use iso_c_binding
    use Sort
    type(WMPoint),allocatable,dimension(:),target,intent(inout)  :: WMPoints !Requires f95TS
    type(WMPoint),allocatable,dimension(:) :: TMP
    integer i,n

    !Choose the one closer to a wall. If of the same distance, choose the later one.
    !For wider compatibility we do not use MOLD= or SOURCE= in allocate.

    allocate(TMP(size(WMPoints)))

    call qsort(c_loc(WMPoints(1)),&
               size(WMPoints,kind = c_size_t),&
               int(storage_size(WMPoints)/storage_size('a'),c_size_t),& !c_sizeof not supported by Solaris Studio 12.3
               c_funloc(CompareWMPoints))

    TMP(1) = WMPoints(1)
    n = 1
    do i = 2,size(WMPoints)
        if (WMPoints(i-1)%xi/=WMPoints(i)%xi .or.&
            WMPoints(i-1)%yj/=WMPoints(i)%yj .or.&
            WMPoints(i-1)%zk/=WMPoints(i)%zk) then !if the point is not duplicate of previous one

              n = n+1
              TMP(n) = WMPoints(i)

        end if
    end do

    deallocate(WMPoints) !would not be needed in Fortran 2003
    allocate(WMPoints(n))!

    WMPoints = TMP(1:n)

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
      res =  0_C_int
    end if

  end function CompareWMPoints

















  subroutine GetOutsideBoundariesWM(nscalars)
    integer, intent(in) :: nscalars
    integer       :: i,j,k
    type(WMPoint) :: WMP

    allocate(WMP%depscalar(nscalars))
    WMP%depscalar = 0

    if (Btype(We)==NOSLIP) then
      do k = 1,Prnz
       do j = 1,Prny
         if (Prtype(1,j,k)<=0) then
           WMP%xi = 1
           WMP%yj = j
           WMP%zk = k
           WMP%distx = (xPr(1)-xU(0))
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
         if (Prtype(Prnx,j,k)<=0) then
           WMP%xi = Prnx
           WMP%yj = j
           WMP%zk = k
           WMP%distx = (xPr(Prnx)-xU(Unx+1))
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
         if (Prtype(i,1,k)<=0) then
           WMP%xi = i
           WMP%yj = 1
           WMP%zk = k
           WMP%distx = 0
           WMP%disty = (yPr(1)-yV(0))
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
         if (Prtype(i,Prny,k)<=0) then
           WMP%xi = i
           WMP%yj = Prny
           WMP%zk = k
           WMP%distx = 0
           WMP%disty = (yPr(Prny)-yV(Vny+1))
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
         if (Prtype(i,j,1)<=0) then

           WMP%xi = i
           WMP%yj = j
           WMP%zk = 1
           WMP%distx = 0
           WMP%disty = 0
           WMP%distz = (zPr(1)-zW(0))
           WMP%ustar = 1

           if (Btype(Bo)==DIRICHLET) then
             WMP%wallu = sideU(1,Bo)
             WMP%wallv = sideU(2,Bo)
             WMP%wallw = 0
           end if

           WMP%z0 = z0B

           if (TempBtype(Bo)==CONSTFLUX) then
             WMP%tempfl = sideTemp(Bo)
           else
             WMP%temp = 0
           end if

           if (TempBtype(Bo)==DIRICHLET) then
             WMP%temp = sideTemp(Bo)
           end if

           call AddWMPoint(WMP)

         end if
       end do
      end do
    end if

    if (Btype(To)==NOSLIP.or.(Btype(To)==DIRICHLET.and.sideU(3,To)==0)) then

      do j = 1,Prny
       do i = 1,Prnx
         if (Prtype(i,j,Prnz)<=0) then
           WMP%xi = i
           WMP%yj = j
           WMP%zk = Prnz
           WMP%distx = 0
           WMP%disty = 0
           WMP%distz = (zPr(Prnz)-zW(Wnz+1))
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

  end subroutine GetOutsideBoundariesWM































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


  pure subroutine WMFlatVisc(visc,ustar,distvect,uvect,walluvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: distvect(3),uvect(3),walluvect(3)
    real(knd) vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - walluvect

    vect = vect - dot_product(vect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    if (vel/=0) then

      call WMFlatUstar(ustar,vel,dist)

    end if

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

  end subroutine WMFlatVisc




  pure subroutine WMRoughUstar(ustar,vel,dist,z0)
    real(knd),intent(inout) :: ustar
    real(knd),intent(in) :: vel,dist,z0
    real(knd),parameter  :: eps = 1e-4_knd
    real(knd),parameter  :: yplcrit = 11.225_knd

    if (dist<=z0) then
     if (Re>0) then
      call WMFlatUstar(ustar,vel,dist)
     else
      ustar = 0
     end if
    else
      ustar = vel * 0.41_knd / log(dist/z0)
    end if

  end subroutine WMRoughUstar


  pure subroutine WMRoughVisc(visc,ustar,z0,distvect,uvect,walluvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0
    real(knd),intent(in)    :: distvect(3),uvect(3),walluvect(3)
    real(knd) vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - walluvect

    vect = vect - dot_product(vect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    if (vel/=0) then

      call WMRoughUstar(ustar,vel,dist,z0)

    end if

    if (ustar<0) ustar = 0

    if (vel>0 .and. ustar**2 * dist/vel>1._knd/Re) then
      visc = ustar**2 * dist/vel
    else if (Re>0) then
      visc = 1._knd/Re
    else
      visc = 0
    end if

  end subroutine WMRoughVisc





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


  pure real(knd) function Obukhov_zL(ustar,tempfl,tempref,g,z)
    real(knd),intent(in):: ustar,tempfl,tempref,g,z

    Obukhov_zL = z*(0.4_knd*(g/tempref)*tempfl)/(-ustar**3)
  end function Obukhov_zL



  pure subroutine WM_MO_FLUX_ustar(vel,dist,ustar,z0,tempflux,Re,temperature_ref,grav_acc)
    implicit none

    real(knd),intent(inout) :: ustar
    real(knd),parameter  :: eps = 1e-3
    real(knd),parameter  :: yplcrit = 11.225_knd
    real(knd),intent(in) :: vel,dist,z0,tempflux
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
       zL2 = Obukhov_zL(ustar,tempflux,temperature_ref,grav_acc,dist)
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



  pure subroutine WM_MO_DIRICHLET_ustar_tfl(ustar,tempflux,vel,dist,z0,tempdif)
    real(knd),intent(inout) :: ustar,tempflux
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
                          tempflux = 0
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
          tempflux = 0
        else
          ustar = vel*0.4_knd/(log(dist/z0)-PsiM_MO(zL))
          tempflux = 0.4_knd*ustar*tempdif/(log(dist/z0)-PsiH_MO(zL))
        end if
      end if
    end if
  end subroutine WM_MO_DIRICHLET_ustar_tfl







  pure subroutine WM_MO_FLUX(visc,ustar,tempfl,z0,distvect,uvect)
    real(knd),intent(out) :: visc
    real(knd),intent(inout) :: ustar
    real(knd),intent(in)    :: z0,tempfl
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - dot_product(uvect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    if (vel/=0) then
      call WM_MO_FLUX_ustar(vel,dist,ustar,z0,tempfl,Re,temperature_ref,grav_acc)
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



  pure subroutine WM_MO_DIRICHLET(visc,ustar,tempfl,z0,tempdif,distvect,uvect)
    real(knd),intent(out)   :: visc
    real(knd),intent(inout) :: ustar,tempfl
    real(knd),intent(in)    :: z0
    real(knd),intent(in)    :: tempdif ! temperature difference surface - nearest point
    real(knd),intent(in)    :: distvect(3),uvect(3)
    real(knd) vect(3),vel,dist

    dist = sqrt(sum(distvect**2))

    vect = uvect - dot_product(uvect,distvect) * distvect / dist**2  !tangential part

    vel = sqrt(sum(vect**2))

    call WM_MO_DIRICHLET_ustar_tfl(ustar,tempfl,vel,dist,z0,tempdif)

    if (ustar<0) ustar = 0

    if ((vel>0.or.tempfl/=0)) then
      visc = max(ustar*ustar*dist/vel,1._knd/Re)
    else if (Re>0) then
      visc = 1._knd/Re
    else
      visc = 0
    end if

  end subroutine WM_MO_DIRICHLET



  pure subroutine BOUND_tempfl(Nu)
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
  end subroutine BOUND_tempfl





  subroutine InitTempFL(Temperature)
    real(knd),intent(in) :: Temperature(-1:,-1:,-1:)
    integer i

    if (enable_buoyancy==1.and.TempBtype(Bo)==DIRICHLET) then
      do i = 1,size(WMPoints)
        if (WMPoints(i)%zk==1) WMPoints(i)%tempfl = -TDiff(WMPoints(i)%xi,WMPoints(i)%yj,1)*&
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


    if (enable_buoyancy==1.and.TempBtype(Bo)==DIRICHLET) then
      !$omp parallel private(i,j)
      !$omp do
      do j = 1,Prny
        do i = 1,Prnx
          BsideTArr(i,j) = SurfTemperature(xPr(i),yPr(j), time+dt/2._TIM)
        end do
       end do
      !$omp end do
      !$omp end parallel
    end if



    !$omp parallel do private(i,xi,yj,zk,tdif, vel, wallvel, dist, prgrad)
    do i = 1,size(WMPoints)

      xi = WMPoints(i)%xi
      yj = WMPoints(i)%yj
      zk = WMPoints(i)%zk

      dist(:) = [WMPoints(i)%distx, WMPoints(i)%disty, WMPoints(i)%distz]

      vel(1) = (U(xi,yj,zk)+U(xi-1,yj,zk))/2._knd
      vel(2) = (V(xi,yj,zk)+V(xi,yj-1,zk))/2._knd
      vel(3) = (W(xi,yj,zk)+W(xi,yj,zk-1))/2._knd

      wallvel(:) = [WMPoints(i)%wallu, WMPoints(i)%wallv, WMPoints(i)%wallw]


      if (WMPoints(i)%z0>0) then

         if (enable_buoyancy==1 .and. TempBtype(Bo)==CONSTFLUX) then

           call WM_MO_FLUX(Visc(xi, yj, zk), WMPoints(i)%ustar, WMPoints(i)%tempfl,&
                           WMPoints(i)%z0, dist, vel)

         else if (enable_buoyancy==1 .and. TempBtype(Bo)==DIRICHLET) then

           if (WMPoints(i)%zk==1) WMPoints(i)%temp = BsideTArr(WMPoints(i)%xi,WMPoints(i)%yj)

           tdif = WMPoints(i)%temp - Temperature(xi,yj,zk)

           call WM_MO_DIRICHLET(Visc(xi, yj, zk), WMPoints(i)%ustar, WMPoints(i)%tempfl,&
                                WMPoints(i)%z0, tdif, dist, vel)

           if (zk==1) BsideTFLArr(xi,yj) = WMPoints(i)%tempfl

         else

           call WMRoughVisc(Visc(xi, yj, zk),&
                            WMPoints(i)%ustar, WMPoints(i)%z0,&
                            dist, vel, wallvel)
         end if

       else

         if (Re<=0) then
           stop "The wall model requires positive viscosity or roughness length."
         end if

         if (wallmodeltype == 2) then

           call WallPrGradient(prgrad,xi,yj,zk,Pr,Prtype)

           call WMFlatPrGradVisc(Visc(xi, yj, zk),&
                           WMPoints(i)%ustar,&
                           dist, vel, wallvel, prgrad)

         else


           call WMFlatVisc(Visc(xi, yj, zk),&
                           WMPoints(i)%ustar,&
                           dist, vel, wallvel)

         end if
       end if

    end do
    !$omp end parallel do

    if (enable_buoyancy==1.and. TempBtype(Bo)==DIRICHLET) call Bound_tempfl(BsideTFLArr)


  end subroutine ComputeViscsWM






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
          depos(WMPoints(j)%xi,WMPoints(j)%yj,i) = depos(WMPoints(j)%xi,WMPoints(j)%yj,i) + WMPoints(j)%depscalar(i)
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