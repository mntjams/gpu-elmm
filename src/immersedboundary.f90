module ImmersedBoundaryWM
  use Parameters
  use SolidBodies
  use GeometricShapes, only: Terrain, Ray
  use WallModels
  
  implicit none

  private
  public GetSolidBodiesWM, GetSolidBodiesWM_UVW, GetWMFluxes
 
contains
  
  subroutine GetSolidBodiesWM
    type(WMPoint)            :: WMP
    type(SolidBody),pointer :: CurrentSB => null()
    integer                  :: neighbours(3,6)
    real(knd)     :: dist,nearx,neary,nearz
    integer       :: i,j,k,m,n,o,p
    integer       :: nb

    allocate(WMP%depscalar(num_of_scalars))

    !six triplets [1,0,0], [-1,0,0], [0,1,0],...
    neighbours = 0
    neighbours(1,1) =  1
    neighbours(1,2) = -1
    neighbours(2,3) =  1
    neighbours(2,4) = -1
    neighbours(3,5) =  1
    neighbours(3,6) = -1

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)<=0) then
        dist = huge(dist)
        nb = 0
         
        do p=1,6
           m=neighbours(1,p)
           n=neighbours(2,p)
           o=neighbours(3,p)

           if ((Prtype(i+m,j+n,k+o)>0).and.Prtype(i+m,j+n,k+o)/=nb.and.(sum(abs([m,n,o]))==1)) then
             call SetCurrentSB(CurrentSB,Prtype(i+m,j+n,k+o))
             call CurrentSB%Closest(nearx,neary,nearz,xPr(i),yPr(j),zPr(k))
             if (sqrt((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)<dist) then
              dist = sqrt((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)
              nb = Prtype(i+m,j+n,k+o)
             end if
           end if
        end do

        if (nb>0) then
                    
          call SetCurrentSB(CurrentSB,nb)
          
          WMP%xi = i
          WMP%yj = j
          WMP%zk = k
          WMP%distx = nearx-xPr(i)
          WMP%disty = neary-yPr(j)
          WMP%distz = nearz-zPr(k)
          WMP%ustar = 1

          select type (geomshape => CurrentSB%GeometricShape)
            type is (Terrain)
              if (geomshape%PrPoints(i,j)%rough) then
                WMP%z0 = geomshape%PrPoints(i,j)%z0
              else
                WMP%z0 = 0
              end if
            class default
              if (CurrentSB%rough) then
                WMP%z0 = CurrentSB%z0
              else
                WMP%z0 = 0
              end if
          end select

          call AddWMPoint(WMP)
          
        end if
       end if
      end do
     end do
    end do

    
  end subroutine GetSolidBodiesWM
  
  
  
  
  subroutine GetSolidBodiesWM_UVW
    integer                  :: neighbours(3,MINUSX:PLUSZ)
    real(knd), target        :: r_neighbours(3,MINUSX:PLUSZ)

    !six triplets [1,0,0], [-1,0,0], [0,1,0],...
    neighbours = 0
    neighbours(1, MINUSX) = -1
    neighbours(1, PLUSX)  =  1
    neighbours(2, MINUSY) = -1
    neighbours(2, PLUSY)  =  1
    neighbours(3, MINUSZ) = -1
    neighbours(3, PLUSZ)  =  1
    
    r_neighbours = real(neighbours, knd)
    
    call helper(1, Unx, Uny, Unz, xU(-2:), yPr, zPr, Utype)

    call helper(2, Vnx, Vny, Vnz, xPr, yV(-2:), zPr, Vtype)

    call helper(3, Wnx, Wny, Wnz, xPr, yPr, zW(-2:), Wtype)

  contains
  
    subroutine helper(component, nx, ny, nz, x, y, z, Xtype)
      integer, intent(in) :: component, nx, ny, nz
      real(knd), intent(in) :: x(-2:), y(-2:), z(-2:)
      integer, intent(in) :: Xtype(-2:,-2:,-2:)
      type(WMPointUVW) :: p
      type(SolidBody), pointer :: SB
      real(knd), pointer, contiguous :: dirvec(:)
      real(knd)     :: nearx, neary, nearz, t, distvec(3)
      integer       :: i, j, k, m, n, o, dir

      do k = 1, nz
       do j = 1, ny
        do i = 1, nx
        
          if (Xtype(i,j,k) < 0) then

            do dir = MINUSX, PLUSZ

              dirvec => r_neighbours(:,dir)

              m = neighbours(1,dir)
              n = neighbours(2,dir)
              o = neighbours(3,dir)
              
              if ((Xtype(i+m,j+n,k+o)>0)) then
              
                call SetCurrentSB(SB, Xtype(i+m,j+n,k+o))
                
                call SB%Closest(nearx ,neary, nearz, x(i), y(j), z(k))
                
                distvec = [nearx - x(i), neary - y(j), nearz - z(k)]
                
                if (.not. right_direction(distvec, dirvec)) then
                  t = SB%ClosestOnLineOut( x(i+m), y(j+n), z(k+o), &
                                                   x(i),   y(j),   z(k) )
                  nearx = xU(i+m) + t * ( x(i) - x(i+m) )
                  neary = yPr(j+n) + t * ( y(j) - y(j+n) )
                  nearz = zPr(k+o) + t * ( z(k) - z(k+o) )
                  distvec = [nearx - x(i), neary - y(j), nearz - z(k)]
                end if
                
                p%xi = i
                p%yj = j
                p%zk = k
                
                p%distx = distvec(1)
                p%disty = distvec(2)
                p%distz = distvec(3)
                
                p%ustar = 1
                
                select type (geomshape => SB%GeometricShape)
                  type is (Terrain)
                  
                    select case (component)
                      case (1)
                        if (geomshape%UPoints(i,j)%rough) then
                          p%z0 = geomshape%UPoints(i,j)%z0
                        else
                          p%z0 = 0
                        end if
                      case (2)
                        if (geomshape%VPoints(i,j)%rough) then
                          p%z0 = geomshape%VPoints(i,j)%z0
                        else
                          p%z0 = 0
                        end if
                      case (3)
                        if (geomshape%PrPoints(i,j)%rough) then
                          p%z0 = geomshape%PrPoints(i,j)%z0
                        else
                          p%z0 = 0
                        end if
                    end select
                    
                  class default
                  
                    if (SB%rough) then
                      p%z0 = SB%z0
                    else
                      p%z0 = 0
                    end if
                    
                end select

                call AddWMPointUVW(p, component, dir)
              
              end if

            end do
             
          end if

        end do
       end do
      end do
    end subroutine
  
    logical function right_direction(a, b)
      !checks if the angle between two vectors is small enough
      real(knd), intent(in) :: a(3), b(3)
      right_direction = ( dot_product(a, b) / (norm2(a) * norm2(b)) ) > 0.4
    end function

  end subroutine GetSolidBodiesWM_UVW
  
  
  subroutine GetWMFluxes
    use SolarRadiation
    use PhysicalProperties
    use VolumeSources
    
    real(knd) :: inc_radiation_flux, radiation_balance, angle_to_sun, &
                 total_heat_flux, sensible_heat_flux, latent_heat_flux
    real(knd) :: out_norm(3), xr, yr, zr, svf
    
    type(Ray) :: r
    
    integer :: i
    
real(knd),allocatable :: svf_array(:,:,:)   
allocate(svf_array(1:Prnx,1:Prny,1));svf_array = 0
    
    enable_radiation = .true.
    
    if (enable_buoyancy .and. enable_radiation) then
      do i=1,size(WMPoints)
      
        associate(p => WMPoints(i))
          out_norm = [-p%distx, -p%disty, -p%distz] / &
                     hypot(hypot(p%distx,p%disty),p%distz)
      
          angle_to_sun = max(0._knd, dot_product(out_norm, vector_to_sun))


          xr = xPr(p%xi)+p%distx*0.9_knd
          yr = yPr(p%yj)+p%disty*0.9_knd
          zr = zPr(p%zk)+p%distz*0.9_knd
          if (angle_to_sun>0._knd) then
            r = Ray([xr, yr, zr], &
                    vector_to_sun)
                     
            if ((SolidBodiesList%any(DoIntersect)) .or. &
                (SourceBodiesList%any(DoIntersectPlants))) then
              angle_to_sun = 0
            end if
          end if

          !temporary
          svf = sky_view_factor(xr,yr,zr)
      
          inc_radiation_flux = angle_to_sun * solar_direct_flux() + &
                               solar_diffuse_flux()*svf + &
                               in_lw_radiation()

          radiation_balance = inc_radiation_flux * (1-p%albedo) - &
                              out_lw_radiation(p%emmissivity, temperature_ref) * svf

          total_heat_flux = radiation_balance - &
                            (0.1 * radiation_balance) !crude guess of the storage flux

          if (enable_moisture) then
            latent_heat_flux = total_heat_flux * p%evaporative_fraction
            sensible_heat_flux = total_heat_flux - latent_heat_flux
            p%moisture_flux = latent_heat_flux / (rho_air_ref * Lv_water_ref)
          else
            sensible_heat_flux = total_heat_flux
          end if
          
          p%temperature_flux = sensible_heat_flux / (rho_air_ref * Cp_air_ref)

        end associate  
        
      end do
      
      call SaveFluxes
    end if
    
    
    contains
    
      real(knd) function sky_view_factor(x,y,z)
        real(knd),intent(in) :: x,y,z
        integer :: i,nfree

        
        nfree = 0
        
        do i=1,svf_nrays
          r = Ray([x,y,z],svf_vecs(:,i))
        
          if (.not.((SolidBodiesList%any(DoIntersect)) .or. &
                    (SourceBodiesList%any(DoIntersectPlants)))) then
            nfree = nfree + 1
          end if
        end do
        
        sky_view_factor = real(nfree,knd) / svf_nrays
        
      end function
  
      logical function DoIntersect(item) result(res)
        class(*) :: item
        
        select type (item)
          class is (SolidBody)
            res = item%IntersectsRay(r)
          class default
            stop "Non-SolidBody i the list."
        end select
      end function
      
      logical function DoIntersectPlants(item) result(res)
        class(*) :: item
        
        select type (item)
          class is (PlantBody)
            res = item%IntersectsRay(r)
          class is (VolumeSourceBody)
            res = .false.
          class default
            stop "Non-VolumeSourceBody i the list."
        end select
      end function
      
      subroutine SaveFluxes
        use VTKArray
        real(knd),allocatable :: temperature_flux(:,:,:)
        allocate(temperature_flux(1:Prnx,1:Prny,1:Prnz))
        
        temperature_flux = 0
        
        do i=1,size(WMPoints)
          associate(p => WMPoints(i))
            temperature_flux(p%xi,p%yj,p%zk) = p%temperature_flux
          end associate
        end do
        
        call VtkArraySimple("output/tempfl.vtk",temperature_flux)
        
        call VtkArraySimple("output/svf.vtk",svf_array)
        
      end subroutine
      
  end subroutine GetWMFluxes

end module ImmersedBoundaryWM


















module IBPoint_types
  use Kinds
  
  implicit none

  type TInterpolationPoint
    integer   :: xi                !coordinates of the interpolation points
    integer   :: yj
    integer   :: zk
    real(knd) :: coef              !interpolation coefficients for the interpolation points
  endtype TInterpolationPoint

  type :: TVelIBPoint
    integer   :: component      !1..U, 2..V, 3..W
    integer   :: xi             !coordinates of the grid point
    integer   :: yj
    integer   :: zk
    real(knd) :: distx          !vector to the nearest boundary point
    real(knd) :: disty
    real(knd) :: distz
    integer   :: dirx           !integer form of the above vector (~sign(distx))
    integer   :: diry
    integer   :: dirz
    integer   :: interp         !kind of interpolation 0.. none (boundarypoint), 1..linear, 2..bilinear, 3..trilinear
    integer   :: interpdir      !direction of interpolation in the linear case, for bilinear it is normal direction to the interpolation plane
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
 !  contains
 !    procedure Create         => TVelIBPoint_Create
  end type TVelIBPoint



  type :: TScalFlIBPoint
    integer                :: xi                       !coordinates of the grid point
    integer                :: yj
    integer                :: zk
    real(knd)              :: dist                     !distance to the boundary
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    integer                :: interp                   !kind of interpolation 1.. none (1 point outside), 2..linear, 4..bilinear  other values not allowed
    real(knd)              :: temperature_flux = 0      !desired temperature flux
    real(knd)              :: moisture_flux = 0      !desired temperature flux
    integer :: n_WMPs = 0 !number of associated wall model points feeding the  !  contains
 !    procedure Create         => TScalFlIBPoint_Create
  end type TScalFlIBPoint
end module IBPoint_types


module VelIBPoint_list
  use IBPoint_types
#define TYPEPARAM type(TVelIBPoint)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module VelIBPoint_list


module ScalFlIBPoint_list
  use IBPoint_types
#define TYPEPARAM type(TScalFlIBPoint)
#include "list-inc-def.f90"
contains
#include "list-inc-proc.f90"
#undef TYPEPARAM
end module ScalFlIBPoint_list






module ImmersedBoundary

  use Parameters
  use IBPoint_types
  use VelIBPoint_list, only: VelIBPointList => list
  use ScalFlIBPoint_list, only: ScalFlIBPointList => list
  use SolidBodies
  use ImmersedBoundaryWM

  implicit none

  private

  public TIBPoint, TIBPoint_Interpolate, TIBPoint_MomentumSource, TIBPoint_ScalFlSource, TIBPoint_Viscosity, &
         UIBPoints, VIBPoints, WIBPoints, ScalFlIBPoints, &
         GetSolidBodiesBC, InitIBPFluxes!, SetIBPFluxes
         !InitSolidBodies imported from SolidBodies




  interface Create
    module procedure TVelIBPoint_Create
    module procedure TScalFlIBPoint_Create
  end interface Create



  type TIBPoint
    integer   :: xi
    integer   :: yj
    integer   :: zk
    real(knd) :: dist
    real(knd) :: temperature_flux = 0
    real(knd) :: moisture_flux = 0
    integer   :: interp
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    integer :: n_WMPs = 0 !number of associated wall model points feeding the fluxes
 !   contains
 !     procedure Interpolate      => TIBPoint_Interpolate
 !     procedure InterpolateTDiff => TIBPoint_Interpolate_TDiff
 !     procedure ScalFlSource     => TIBPoint_ScalFlSource
 !     procedure MomentumSource   => TIBPoint_MomentumSource
 !     procedure Viscosity        => TIBPoint_Viscosity
  end type TIBPoint

  type(VelIBPointList)  :: UIBPointsList, VIBPointsList, WIBPointsList
  type(ScalFlIBPointList)  :: ScalFlIBPointsList

  type(TIBPoint),dimension(:),allocatable,save :: UIBPoints, VIBPoints, WIBPoints
  type(TIBPoint),dimension(:),allocatable,save :: ScalFlIBPoints

  interface assignment (=)
    module procedure VelIBPtoIBP
  end interface

  interface assignment (=)
    module procedure ScalFlIBPtoIBP
  end interface
  
  interface IB_interpolation_coefs
    module procedure IB_interpolation_coefs_1st_order
  end interface

  !number of least square interpolation points
  !must be higher than the interpolation polynomial order
  integer, parameter :: n_ls = 2
   
contains


  subroutine VelIBPtoIBP(IBP,VelIBP)
    type(TIBPoint),intent(out)    :: IBP
    type(TVelIBPoint),intent(in)  :: VelIBP

    IBP%xi = VelIBP%xi
    IBP%yj = VelIBP%yj
    IBP%zk = VelIBP%zk
    IBP%interp = VelIBP%interp

    allocate(IBP%IntPoints(size(VelIBP%IntPoints)))

    IBP%IntPoints = VelIBP%IntPoints
  end subroutine VelIBPtoIBP

  subroutine ScalFlIBPtoIBP(IBP,ScalFlIBP)
    type(TIBPoint),intent(out)    :: IBP
    type(TScalFlIBPoint),intent(in)  :: ScalFlIBP

    IBP%xi = ScalFlIBP%xi
    IBP%yj = ScalFlIBP%yj
    IBP%zk = ScalFlIBP%zk
    IBP%dist = ScalFlIBP%dist
    IBP%temperature_flux = ScalFlIBP%temperature_flux
    IBP%moisture_flux = ScalFlIBP%moisture_flux
    IBP%interp = ScalFlIBP%interp

    allocate(IBP%IntPoints(size(ScalFlIBP%IntPoints)))

    IBP%IntPoints = ScalFlIBP%IntPoints
    
    IBP%n_WMPs = ScalFlIBP%n_WMPs
  end subroutine ScalFlIBPtoIBP


  
  pure recursive function TIBPoint_Interpolate(IBP,U,lb) result(Uint)
    real(knd) :: Uint
    type(TIBPoint),intent(in) :: IBP
    integer,intent(in) :: lb
    real(knd),dimension(lb:,lb:,lb:),intent(in) :: U
    integer i

    Uint = 0

    do i=1,IBP%interp
      Uint = Uint + IBP%IntPoints(i)%coef * U(IBP%IntPoints(i)%xi,&
                                              IBP%IntPoints(i)%yj,&
                                              IBP%IntPoints(i)%zk)
    end do
  end function TIBPoint_Interpolate


  pure recursive function TIBPoint_InterpolateTDiff(IBP,U) result(Uint)
    real(knd) :: Uint
    type(TIBPoint),intent(in) :: IBP
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: U
    integer i,n

    n = 0
    Uint = 0

    do i=1,IBP%interp
      if (abs(IBP%IntPoints(i)%coef-0._knd)>epsilon(1._knd)) then
        Uint = Uint + U(IBP%IntPoints(i)%xi,&
                        IBP%IntPoints(i)%yj,&
                        IBP%IntPoints(i)%zk)
        n = n + 1
      end if
    end do

    Uint = Uint + U(IBP%xi,IBP%yj,IBP%zk)
    Uint = Uint / (n+1)
  end function TIBPoint_InterpolateTDiff

  
  pure recursive function TIBPoint_ScalFlSource(IBP,Scalar,sctype)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
    real(knd) :: src

    type(TIBPoint),intent(in) :: IBP
    real(knd),intent(in)      :: Scalar(-1:,-1:,-1:)
    integer,intent(in)        :: sctype

    real(knd) intscal,intTDiff

    intscal = TIBPoint_Interpolate(IBP,Scalar,-1)

    if (sctype==ScalarTypeTemperature) then
      intTDiff = TIBPoint_InterpolateTDiff(IBP,TDiff)

      if (intTDiff>0)  intscal = intscal + IBP%temperature_flux * IBP%dist / intTDiff
    else if (sctype==ScalarTypeMoisture) then
      intTDiff = TIBPoint_InterpolateTDiff(IBP,TDiff)

      if (intTDiff>0)  intscal = intscal + IBP%moisture_flux * IBP%dist / intTDiff
    end if

    src = (intscal - Scalar(IBP%xi,IBP%yj,IBP%zk)) / dt

  end function TIBPoint_ScalFlSource


  
  pure recursive function TIBPoint_Viscosity(IBP,Viscosity)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
    real(knd) :: src

    type(TIBPoint),intent(in) :: IBP
    real(knd),intent(in)      :: Viscosity(-1:,-1:,-1:)

    src = TIBPoint_Interpolate(IBP,Viscosity,-1)

  end function TIBPoint_Viscosity


  
  pure recursive function TIBPoint_MomentumSource(IBP,U) result(src)
    real(knd) :: src

    type(TIBpoint),intent(in) :: IBP
    real(knd),intent(in)                      :: U(-2:,-2:,-2:)

    src = (TIBPoint_Interpolate(IBP,U,-2) - U(IBP%xi,IBP%yj,IBP%zk)) / dt

  end function TIBPoint_MomentumSource


  recursive subroutine TVelIBPoint_Create(IBP,xi,yj,zk,xU,yU,zU,Utype,component)
    type(TVelIBPoint),intent(out)               :: IBP
    integer,intent(in)                          :: xi,yj,zk
    real(knd),dimension(-2:),intent(in)         :: xU,yU,zU
    integer,dimension(-2:,-2:,-2:),intent(in)   :: Utype
    integer,intent(in)                          :: component

    type(SolidBody),pointer :: SB
    integer :: dirx,diry,dirz,n1,n2,nx,ny,nz
    real(knd) :: x,y,z
    real(knd) :: x2,y2,z2
    real(knd) :: xnear,ynear,znear
    real(knd) :: tx,ty,tz
    logical :: freexm,freeym,freezm
    logical :: freexp,freeyp,freezp
    logical :: free1xm,free1ym,free1zm
    logical :: free1xp,free1yp,free1zp
    real(knd) :: distxm, distym, distzm
    real(knd) :: distxp, distyp, distzp
    integer :: i

    !real coordinates of the IB forcing point
    x = xU(xi)
    y = yU(yj)
    z = zU(zk)
    call SetCurrentSB(SB,Utype(xi,yj,zk))
    
    IBP%component = component
    
    !integer grid coordinates
    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk
! call null_point; return
    if (.not. SB%Inside(x,y,z,0._knd)) then
      call null_point
      return
    end if

    call SB%ClosestOut(xnear,ynear,znear,x,y,z)

    if (hypot(xnear,hypot(ynear,znear))<=0.1_knd*(dxmin*dymin*dzmin)**(1._knd/3)) then
      call null_point
      return
    end if

    freexm = all([ ( Utype(xi-i,yj  ,zk  )<=0, i = 1, n_ls ) ])
    freeym = all([ ( Utype(xi  ,yj-i,zk  )<=0, i = 1, n_ls ) ])
    freezm = all([ ( Utype(xi  ,yj  ,zk-i)<=0, i = 1, n_ls ) ])
    freexp = all([ ( Utype(xi+i,yj  ,zk  )<=0, i = 1, n_ls ) ])
    freeyp = all([ ( Utype(xi  ,yj+i,zk  )<=0, i = 1, n_ls ) ])
    freezp = all([ ( Utype(xi  ,yj  ,zk+i)<=0, i = 1, n_ls ) ])
    
    free1xm = Utype(xi-1,yj  ,zk  )<=0
    free1ym = Utype(xi  ,yj-1,zk  )<=0
    free1zm = Utype(xi  ,yj  ,zk-1)<=0
    free1xp = Utype(xi+1,yj  ,zk  )<=0
    free1yp = Utype(xi  ,yj+1,zk  )<=0
    free1zp = Utype(xi  ,yj  ,zk+1)<=0
    
    if ((free1xm.and.free1xp) .or. (free1ym.and.free1yp) .or. (free1zm.and.free1zp)) then
      call null_point
      return
    end if  

    if ( (.not.freexm) .and. (.not.freexp) ) then
      dirx = 0
    else if (freexm) then
      dirx = -1
    else
      dirx = 1
    end if
      
    if ( (.not.freeym) .and. (.not.freeyp) ) then
      diry = 0
    else if (freeym) then
      diry = -1
    else
      diry = 1
    end if
      
    if ( (.not.freezm) .and. (.not.freezp) ) then
      dirz = 0
    else if (freezm) then
      dirz = -1
    else
      dirz = 1
    end if
      
    tx = 0
    ty = 0
    tz = 0
      
    if (dirx/=0) then
      x2 = xU(xi+dirx)
      tx = SB%ClosestOnLineOut(x,y,z,x2,y,z)
      
      if (tx>1 .or. tx<0.05) then
        dirx = 0
        tx = 0
      end if
      
      IBP%distx = tx * (xU(xi+dirx) - xU(xi))
    end if
      
    if (diry/=0) then
      y2 = yU(yj+diry)
      ty = SB%ClosestOnLineOut(x,y,z,x,y2,z)
      
      if (ty>1 .or. ty<0.05) then
        diry = 0
        ty = 0
      end if
      
      IBP%disty = ty * (yU(yj+diry) - yU(yj))
    end if
      
    if (dirz/=0) then
      z2 = zU(zk+dirz)
      tz = SB%ClosestOnLineOut(x,y,z,x,y,z2)
      
      if (tz>1 .or. tz<0.05) then
        dirz = 0
        tz = 0
      end if
      
      IBP%distz = tz * (zU(zk+dirz) - zU(zk))
    end if
    
    IBP%dirx = dirx
    IBP%diry = diry
    IBP%dirz = dirz

    IBP%IntPoints = InterpolationPoints(IBP,xU,yU,zU)
     
    IBP%interp = size(IBP%IntPoints)
    
  contains
  
    subroutine null_point
      IBP%interp = 0
      allocate(IBP%IntPoints(0))
    end subroutine
    
  end subroutine TVelIBPoint_Create


  recursive subroutine TVelIBPoint_Create_old(IBP,xi,yj,zk,xU,yU,zU,Utype,component)
    type(TVelIBPoint),intent(out)               :: IBP
    integer,intent(in)                          :: xi,yj,zk
    real(knd),dimension(-2:),intent(in)         :: xU,yU,zU
    integer,dimension(-2:,-2:,-2:),intent(in)   :: Utype
    integer,intent(in)                          :: component

    type(SolidBody),pointer :: SB
    integer dirx,diry,dirz,n1,n2,nx,ny,nz
    real(knd) x,y,z,xnear,ynear,znear,t
    logical free100,free010,free001
    real(knd) x2,y2,z2
    integer :: i

    x = xU(xi)                                !real coordinates of the IB forcing point
    y = yU(yj)
    z = zU(zk)
    call SetCurrentSB(SB,Utype(xi,yj,zk))
    call SB%ClosestOut(xnear,ynear,znear,x,y,z)

    IBP%component = component
    IBP%xi = xi                                !integer grid coordinates
    IBP%yj = yj
    IBP%zk = zk
    IBP%distx = xnear-x                       !real distance to the boundary in the x,y,z direction
    IBP%disty = ynear-y
    IBP%distz = znear-z
    IBP%dirx = nint(sign(1.0_knd,IBP%distx))  !integer value denoting direction to the boundary
    IBP%diry = nint(sign(1.0_knd,IBP%disty))
    IBP%dirz = nint(sign(1.0_knd,IBP%distz))

    dirx = abs(IBP%dirx)                      !local temporary variable with abs(dir)
    diry = abs(IBP%diry)
    dirz = abs(IBP%dirz)

    if (.not.SB%Inside(xU(xi),yPr(yj),zPr(zk)))  then  !For now, if actually outside the body, set an artificial boundary
      IBP%interp = 0                                   !point here. In future we can use another interpolation.
      IBP%interpdir = 0

      allocate(IBP%IntPoints(0))

      return

    end if


    if (abs(IBP%distx)<(xU(xi+1)-xU(xi-1))/1000._knd) then      !if too close to the boundary, set the distance to 0
      IBP%distx = 0
      dirx = 0
      IBP%dirx = 0
    end if

    if (abs(IBP%disty)<(yU(yj+1)-yU(yj-1))/1000._knd) then
      IBP%disty = 0
      diry = 0
      IBP%diry = 0
    end if

    if (abs(IBP%distz)<(zU(zk+1)-zU(zk-1))/1000._knd) then
      IBP%distz = 0
      dirz = 0
      IBP%dirz = 0
    end if


    if (dirx==0) then                                               !If there is a cell free of solid bodies in the direction
      free100=.false.
    else
      free100 = all([ ( Utype(xi+IBP%dirx*i,yj,zk)<=0, i=1,3 ) ])
    end if

    if (diry==0) then
      free010=.false.
    else
      free010 = all([ ( Utype(xi,yj+IBP%diry*i,zk)<=0, i=1,3 ) ])
    end if

    if (dirz==0) then
      free001=.false.
    else
      free001 = all([ ( Utype(xi,yj,zk+IBP%dirz*i)<=0, i=1,3 ) ])
    end if

    nx = 0
    ny = 0
    nz = 0
    n1 = 0                                                          ! n1 number of neighbouring cells free of solid bodies
    n2 = 0                                                          ! n2 number of nonzero coordinate directions to boundary

    if (Utype(xi+1,yj,zk)<=0) then
      n1 = n1+1
      nx = nx+1
    end if

    if (Utype(xi-1,yj,zk)<=0) then
      n1 = n1+1
      nx = nx+1
    end if

    if (Utype(xi,yj+1,zk)<=0) then
      n1 = n1+1
      ny = ny+1
    end if

    if (Utype(xi,yj-1,zk)<=0) then
      n1 = n1+1
      ny = ny+1
    end if

    if (Utype(xi,yj,zk+1)<=0) then
      n1 = n1+1
      nz = nz+1
    end if

    if (Utype(xi,yj,zk-1)<=0) then
      n1 = n1+1
      nz = nz+1
    end if

    if (dirx/=0) n2 = n2+1
    if (diry/=0) n2 = n2+1
    if (dirz/=0) n2 = n2+1

    if (n1>n2) then   ! If too many free directions, treat as directly on the boundary,
      IBP%interp = 0  ! because we are probably at some edge.
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%dirz = 0
    end if

      !At least in one direction free  cells on both sides.
      !Unresolvably small object or a thin wall -> treat as directly on the boundary.
    if (nx>1.or.ny>1.or.nz>1) then
      IBP%interp = 0
      IBP%interpdir = 0
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%dirz = 0
    end if

    
    if ((dirx==0.and.diry==0.and.dirz==0).or.((.not.free001).and.(.not.free010).and.(.not.free100))) then
                              !If no free direction, the boundary is here.
      IBP%interp = 0
      IBP%interpdir = 0
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%dirz = 0

    elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then !Three free directions,

      IBP%interp = 9                                                                           !use trilinear interpolation.
      IBP%interpdir = 0

    elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
      !Two free directions y and z use bilinear interpolation normal to x.
      IBP%interp = 6
      IBP%interpdir = 1

      if (dirx==1) then
        !If dirx /= 0 then forget the x component of dist vector        
        y2 = yU(yj+IBP%diry)
        z2 = zU(zk+IBP%dirz)
        ! and find an intersection of the new vector with the boundary
        t = SB%ClosestOnLineOut(x,y,z,x,y2,z2)
        
        if (t>=2) then
          IBP%interp = 0
          IBP%interpdir = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = 0
          dirx = 0
          diry = 0
          dirz = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%dirz = 0
        else
          dirx = 0
          IBP%dirx = 0
          IBP%distx = 0
          IBP%disty = (y2-y)*t
          IBP%distz = (z2-z)*t
        end if
      end if

    elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then  !the same normal to y

      IBP%interp = 6
      IBP%interpdir = 2

      if (diry==1) then
        x2 = xU(xi+IBP%dirx)
        z2 = zU(zk+IBP%dirz)
        
        t = SB%ClosestOnLineOut(x,y,z,x2,y,z2)
        
        if (t>=2) then
          IBP%interp = 0
          IBP%interpdir = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = 0
          dirx = 0
          diry = 0
          dirz = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%dirz = 0
        else
          diry = 0
          IBP%diry = 0
          IBP%distx = (x2-x)*t
          IBP%disty = 0
          IBP%distz = (z2-z)*t
        end if
      end if

    elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then  !the same normal to z

      IBP%interp = 6
      IBP%interpdir = 3

      if (dirz==1) then
        x2 = xU(xi+IBP%dirx)
        y2 = yU(yj+IBP%diry)
        
        t = SB%ClosestOnLineOut(x,y,z,x2,y2,z)
        
        if (t>=2) then
          IBP%interp = 0
          IBP%interpdir = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = 0
          dirx = 0
          diry = 0
          dirz = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%dirz = 0
        else
          dirz = 0
          IBP%dirz = 0
          IBP%distx = (x2-x)*t
          IBP%disty = (y2-y)*t
          IBP%distz = 0
        end if
      end if

    elseif (free100.and.dirx==1) then  !Only one free direction, use linear interpolation in direction x.

      IBP%interp = 3
      IBP%interpdir = 1

      if (diry==1.or.dirz==1) then                   !If other dir components nonzero, delete them and
        x2 = xU(xi+IBP%dirx)
        
        t = SB%ClosestOnLineOut(x,y,z,x2,y,z)  !find an intersection of the new vector with the boundary

        if (t>=2) then
          IBP%interp = 0
          IBP%interpdir = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = 0
          dirx = 0
          diry = 0
          dirz = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%dirz = 0
        else
          diry = 0
          dirz = 0
          IBP%diry = 0
          IBP%dirz = 0
          IBP%distx = t*(x2-x)
          IBP%disty = 0
          IBP%distz = 0
        end if
      end if

    elseif (free010.and.diry==1) then               !the same in y

      IBP%interp = 3
      IBP%interpdir = 2

      if (dirx==1.or.dirz==1) then
        y2 = yU(yj+IBP%diry)
        
        t = SB%ClosestOnLineOut(x,y,z,x,y2,z)
        
        if (t>=2) then
          IBP%interp = 0
          IBP%interpdir = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = 0
          dirx = 0
          diry = 0
          dirz = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%dirz = 0
        else
          dirx = 0
          dirz = 0
          IBP%dirx = 0
          IBP%dirz = 0
          IBP%distx = 0
          IBP%disty = t*(y2-y)
          IBP%distz = 0
        end if
      end if

    elseif (free001.and.dirz==1) then               !the same in z

      IBP%interp = 3
      IBP%interpdir = 3

      if (dirx==1.or.diry==1) then
        z2 = zU(zk+IBP%dirz)
        
        t = SB%ClosestOnLineOut(x,y,z,x,y,z2)
        
        if (t>=2) then
          IBP%interp = 0
          IBP%interpdir = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = 0
          dirx = 0
          diry = 0
          dirz = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%dirz = 0
        else
          dirx = 0
          diry = 0
          IBP%dirx = 0
          IBP%diry = 0
          IBP%distx = 0
          IBP%disty = 0
          IBP%distz = (z2-z)
        end if
      end if

    else
      ! We should have find some interpolation and not come here!
      write(*,*) "Assert error line",__LINE__,__FILE__
      write(*,*) "component",component
      write(*,*) "xi,yj,zk",xi,yj,zk
      write(*,*) "xyz",x,y,z
      write(*,*) "xyznear",xnear,ynear,znear
      write(*,*) "interp",IBP%interp
      write(*,*) "interpdir",IBP%interpdir
      write(*,*) "distxyz",IBP%distx,IBP%disty,IBP%distz
      write(*,*) "Utype",Utype(xi,yj,zk)
      write(*,*) "dirxyz",IBP%dirx,IBP%diry,IBP%dirz
      write(*,*) "free",free100,free010,free001
      stop
    end if
    
    if (abs(IBP%distx)>dxmin*2.or. abs(IBP%disty)>dymin*2.or. abs(IBP%distz)>dzmin*2 ) then
      write(*,*) "Assert error line",__LINE__,__FILE__
      write(*,*) "component",component
      write(*,*) "xi,yj,zk",xi,yj,zk
      write(*,*) "xyz",x,y,z
      write(*,*) "xyznear",xnear,ynear,znear
      write(*,*) "interp",IBP%interp
      write(*,*) "interpdir",IBP%interpdir
      write(*,*) "distxyz",IBP%distx,IBP%disty,IBP%distz
      write(*,*) "Utype",Utype(xi,yj,zk)
      write(*,*) "dirxyz",IBP%dirx,IBP%diry,IBP%dirz
      write(*,*) "free",free100,free010,free001
      stop
    end if

    allocate(IBP%IntPoints(IBP%interp))

    call TVelIBPoint_InterpolationCoefs(IBP,xU,yU,zU)

    if (any(IBP%IntPoints%coef>1E-5)) then
      write(*,*) "Assert error line",__LINE__,__FILE__
      write(*,*) "point coefs",IBP%IntPoints%coef
      write(*,*) "component",component
      write(*,*) "xi,yj,zk",xi,yj,zk
      write(*,*) "xyz",x,y,z
      write(*,*) "xyznear",xnear,ynear,znear
      write(*,*) "interp",IBP%interp
      write(*,*) "interpdir",IBP%interpdir
      write(*,*) "distxyz",IBP%distx,IBP%disty,IBP%distz
      write(*,*) "Utype",Utype(xi,yj,zk)
      write(*,*) "dirxyz",IBP%dirx,IBP%diry,IBP%dirz
      write(*,*) "free",free100,free010,free001
      stop
    end if

  end subroutine TVelIBPoint_Create_old

  recursive function InterpolationPoints(IBP,xU,yU,zU) result(res)
    type(TInterpolationPoint),allocatable :: res(:)
    type(TVelIBpoint), intent(in)         :: IBP
    real(knd),dimension(-2:), intent(in)  :: xU,yU,zU
    
    type(TInterpolationPoint),target      :: tmp(3*n_ls)
    type(TInterpolationPoint), pointer    :: t(:)
    real(knd) :: xr, yr, zr, x(0:n_ls), y(0:n_ls), z(0:n_ls)
    real(knd) :: b(3), d(3), c
    integer :: i,xi,yj,zk,dirx,diry,dirz
    
!     interface IB_interpolation_coefs
!       module procedure IB_interpolation_coefs_1st_order
!     end interface
    

    xi = IBP%xi
    yj = IBP%yj
    zk = IBP%zk
    dirx = IBP%dirx
    diry = IBP%diry
    dirz = IBP%dirz
    
    !distance from the wall
    d = 1
    
    if (dirx/=0) then
     
      do i = 1, n_ls
        tmp(i)%xi = xi + i*dirx
        tmp(i)%yj = yj
        tmp(i)%zk = zk
      end do
      
      xr = xU(xi) + IBP%distx
      x  = [ ( xU(xi+i*dirx), i = 0, n_ls ) ]
      tmp(1:n_ls)%coef = IB_interpolation_coefs(xr,x)
      
      d(1)  = abs(x(0)-xr)
        
    end if

    if (diry/=0) then

      do i = 1, n_ls
        tmp(n_ls + i)%xi = xi
        tmp(n_ls + i)%yj = yj + i*diry
        tmp(n_ls + i)%zk = zk
      end do
      
      yr = yU(yj) + IBP%disty
      y  = [ ( yU(yj+i*diry) , i = 0, n_ls ) ]
      tmp(n_ls+1:2*n_ls)%coef = IB_interpolation_coefs(yr,y)
      
      d(2) = abs(y(0)-yr)
      
    end if
    
    if (dirz/=0) then
     
      do i = 1, n_ls
        tmp(2*n_ls + i)%xi = xi
        tmp(2*n_ls + i)%yj = yj
        tmp(2*n_ls + i)%zk = zk + i*dirz
      end do

      zr = zU(zk) + IBP%distz
      z  = [ ( zU(zk+i*dirz), i = 0, n_ls )  ]
      tmp(2*n_ls+1:3*n_ls)%coef = IB_interpolation_coefs(zr,z)
      
      d(3) = abs(z(0)-zr)
      
    end if

    !beta 1,2,3 in eq. 17-19 in Peller et al., doi:1.1002/fld.1227
    b = 0
    
    if (dirx/=0) then
      b(1) = product(d)/(d(1))**2
    end if
    
    if (diry/=0) then
      b(2) = product(d)/(d(2))**2
    end if
    
    if (dirz/=0) then
      b(3) = product(d)/(d(3))**2
    end if
    
    !sum of betas
    c = sum(b)

    allocate(res(0))
    
    !NOTE: when well supported by compilers use associate (cf. http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56386)
    if (dirx/=0) then
      t=>tmp(1:n_ls)
      t%coef = t%coef * b(1)/c
      res = [res, t]
    end if

    if (diry/=0) then
      t=>tmp(n_ls+1:2*n_ls)
      t%coef = t%coef * b(2)/c
      res = [res, t]
    end if

    if (dirz/=0) then
      t=>tmp(2*n_ls+1:3*n_ls)
      t%coef = t%coef * b(3)/c
      res = [res, t]
    end if

  end function

  recursive subroutine TVelIBPoint_InterpolationCoefs(IBP,xU,yU,zU)
    type(TVelIBpoint),intent(inout)     :: IBP
    real(knd),dimension(-2:),intent(in) :: xU,yU,zU
    real(knd) xr, yr, zr, x(0:3), y(0:3), z(0:3)
    real(knd) b1, b2, b3, c
    integer xi,yj,zk,dirx,diry,dirz
    
!     interface IB_interpolation_coefs
!       module procedure IB_interpolation_coefs_1st_order
!     end interface
    

    xi = IBP%xi
    yj = IBP%yj
    zk = IBP%zk

    dirx = IBP%dirx
    diry = IBP%diry
    dirz = IBP%dirz

    if (IBP%interp==3) then


        IBP%IntPoints(1)%xi = xi+dirx
        IBP%IntPoints(1)%yj = yj+diry
        IBP%IntPoints(1)%zk = zk+dirz

        IBP%IntPoints(2)%xi = xi+2*dirx
        IBP%IntPoints(2)%yj = yj+2*diry
        IBP%IntPoints(2)%zk = zk+2*dirz

        IBP%IntPoints(3)%xi = xi+3*dirx
        IBP%IntPoints(3)%yj = yj+3*diry
        IBP%IntPoints(3)%zk = zk+3*dirz

        if (IBP%interpdir==1) then
          xr = xU(xi) + IBP%distx
          x  = [ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  ]
        elseif (IBP%interpdir==2) then
          xr = yU(yj) + IBP%disty
          x  = [ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  ]
        else
          xr = zU(zk) + IBP%distz
          x  = [ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  ]
        end if

        IBP%IntPoints%coef = IB_interpolation_coefs(xr,x)

    elseif (IBP%interp==6) then


        if (IBP%interpdir==1) then

          IBP%IntPoints(1)%xi = xi
          IBP%IntPoints(1)%yj = yj + diry
          IBP%IntPoints(1)%zk = zk

          IBP%IntPoints(2)%xi = xi
          IBP%IntPoints(2)%yj = yj + 2*diry
          IBP%IntPoints(2)%zk = zk

          IBP%IntPoints(3)%xi = xi
          IBP%IntPoints(3)%yj = yj + 3*diry
          IBP%IntPoints(3)%zk = zk

          IBP%IntPoints(4)%xi = xi
          IBP%IntPoints(4)%yj = yj
          IBP%IntPoints(4)%zk = zk + dirz

          IBP%IntPoints(5)%xi = xi
          IBP%IntPoints(5)%yj = yj
          IBP%IntPoints(5)%zk = zk + 2*dirz

          IBP%IntPoints(6)%xi = xi
          IBP%IntPoints(6)%yj = yj
          IBP%IntPoints(6)%zk = zk + 3*dirz

          xr = yU(yj) + IBP%disty
          x  = [ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  ]

          yr = zU(zk) + IBP%distz
          y  = [ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  ]


        elseif (IBP%interpdir==2) then

          IBP%IntPoints(1)%xi = xi
          IBP%IntPoints(1)%yj = yj
          IBP%IntPoints(1)%zk = zk + dirz

          IBP%IntPoints(2)%xi = xi
          IBP%IntPoints(2)%yj = yj
          IBP%IntPoints(2)%zk = zk + 2*dirz

          IBP%IntPoints(3)%xi = xi
          IBP%IntPoints(3)%yj = yj
          IBP%IntPoints(3)%zk = zk + 3*dirz

          IBP%IntPoints(4)%xi = xi + dirx
          IBP%IntPoints(4)%yj = yj
          IBP%IntPoints(4)%zk = zk

          IBP%IntPoints(5)%xi = xi + 2*dirx
          IBP%IntPoints(5)%yj = yj
          IBP%IntPoints(5)%zk = zk

          IBP%IntPoints(6)%xi = xi + 3*dirx
          IBP%IntPoints(6)%yj = yj
          IBP%IntPoints(6)%zk = zk

          xr = zU(zk) + IBP%distz
          x  = [ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  ]

          yr = xU(xi) + IBP%distx
          y  = [ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  ]

        else

          IBP%IntPoints(1)%xi = xi + dirx
          IBP%IntPoints(1)%yj = yj
          IBP%IntPoints(1)%zk = zk

          IBP%IntPoints(2)%xi = xi + 2*dirx
          IBP%IntPoints(2)%yj = yj
          IBP%IntPoints(2)%zk = zk

          IBP%IntPoints(3)%xi = xi + 3*dirx
          IBP%IntPoints(3)%yj = yj
          IBP%IntPoints(3)%zk = zk

          IBP%IntPoints(4)%xi = xi
          IBP%IntPoints(4)%yj = yj + diry
          IBP%IntPoints(4)%zk = zk

          IBP%IntPoints(5)%xi = xi
          IBP%IntPoints(5)%yj = yj + 2*diry
          IBP%IntPoints(5)%zk = zk

          IBP%IntPoints(6)%xi = xi
          IBP%IntPoints(6)%yj = yj + 3*diry
          IBP%IntPoints(6)%zk = zk

          xr = xU(xi) + IBP%distx
          x  = [ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  ]

          yr = yU(yj) + IBP%disty
          y  = [ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  ]

        end if

        IBP%IntPoints(1:3)%coef = IB_interpolation_coefs(xr,x)

        IBP%IntPoints(4:6)%coef = IB_interpolation_coefs(yr,y)

        !beta 1 and  2 in eq. 17-19 in Peller et al., doi:1.1002/fld.1227
        b1 = abs(y(0)-yr)/abs(x(0)-xr)
        b2 = 1/b1
        !sum of betas
        c = b1 + b2

        IBP%IntPoints(1:3)%coef = IBP%IntPoints(1:3)%coef * b1/c

        IBP%IntPoints(4:6)%coef = IBP%IntPoints(4:6)%coef * b2/c

    elseif (IBP%interp==9) then

        IBP%IntPoints(1)%xi = xi + dirx
        IBP%IntPoints(1)%yj = yj
        IBP%IntPoints(1)%zk = zk

        IBP%IntPoints(2)%xi = xi + 2*dirx
        IBP%IntPoints(2)%yj = yj
        IBP%IntPoints(2)%zk = zk

        IBP%IntPoints(3)%xi = xi + 3*dirx
        IBP%IntPoints(3)%yj = yj
        IBP%IntPoints(3)%zk = zk

        IBP%IntPoints(4)%xi = xi
        IBP%IntPoints(4)%yj = yj + diry
        IBP%IntPoints(4)%zk = zk

        IBP%IntPoints(5)%xi = xi
        IBP%IntPoints(5)%yj = yj + 2*diry
        IBP%IntPoints(5)%zk = zk

        IBP%IntPoints(6)%xi = xi
        IBP%IntPoints(6)%yj = yj + 3*diry
        IBP%IntPoints(6)%zk = zk

        IBP%IntPoints(7)%xi = xi
        IBP%IntPoints(7)%yj = yj
        IBP%IntPoints(7)%zk = zk + dirz

        IBP%IntPoints(8)%xi = xi
        IBP%IntPoints(8)%yj = yj
        IBP%IntPoints(8)%zk = zk + 2*dirz

        IBP%IntPoints(9)%xi = xi
        IBP%IntPoints(9)%yj = yj
        IBP%IntPoints(9)%zk = zk + 3*dirz

        xr = xU(xi) + IBP%distx
        x  = [ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  ]

        yr = yU(yj) + IBP%disty
        y  = [ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  ]

        zr = zU(zk) + IBP%distz
        z  = [ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  ]

        IBP%IntPoints(1:3)%coef = IB_interpolation_coefs(xr,x)

        IBP%IntPoints(4:6)%coef = IB_interpolation_coefs(yr,y)

        IBP%IntPoints(7:9)%coef = IB_interpolation_coefs(zr,z)

        !beta 1 and  2 in eq. 17-19 in Peller et al., doi:1.1002/fld.1227
        b1 = abs(y(0)-yr)*abs(z(0)-zr)/abs(x(0)-xr)
        b2 = abs(z(0)-zr)*abs(x(0)-xr)/abs(y(0)-yr)
        b3 = abs(x(0)-xr)*abs(y(0)-yr)/abs(z(0)-zr)
        !sum of betas
        c = b1 + b2 + b3

        IBP%IntPoints(1:3)%coef = IBP%IntPoints(1:3)%coef * b1/c

        IBP%IntPoints(4:6)%coef = IBP%IntPoints(4:6)%coef * b2/c

        IBP%IntPoints(7:9)%coef = IBP%IntPoints(7:9)%coef * b3/c

    else if (IBP%interp/=0) then
        write(*,*) "Unknown interpolation",__FILE__,__LINE__
        stop
    end if

  end subroutine TVelIBPoint_InterpolationCoefs



  pure function IB_interpolation_coefs_2nd_order(xr,x) result(Coefs)
    real(knd),intent(in)  :: xr,x(0:)
    real(knd) :: Coefs(size(x)-1)
    real(knd) :: A1, A2, A4
    integer   :: n

    n = size(x)-1
    if ( n<3 .or. size(x)<=n ) then
      Coefs = 0*x(1:)
    else
      A1 = sum( x(1:) - xr )**2
      A2 = sum( (x(1:)**2 - xr**2) * (x(1:) - xr) )
      A4 = sum( x(1:)**2 - xr**2 )**2
      
      Coefs = ( ( A2 * (x(1:)**2 - xr**2) - A4 * (x(1:) - xr) ) * (x(0)    - xr)&
            +   ( A2 * (x(1:) - xr) - A1 * (x(1:)**2 - xr**2) ) * (x(0)**2 - xr**2) )&
            / (A2**2 - A1*A4)
    end if

  end function


   function IB_interpolation_coefs_1st_order(xr,x) result(Coefs)
    real(knd),intent(in)  :: xr,x(0:)
    real(knd) :: Coefs(size(x)-1)
    real(knd) :: A
    integer   :: n

    n = size(x)-1
    if ( n<2 .or. size(x)<=n ) then
      Coefs = 0*x(1:)
    else
      A = sum( (x(1:) - xr)**2 )

      Coefs = (x(1:) - xr) * (x(0) - xr) / A
    end if

  end function


  subroutine TScalFlIBPoint_Create(IBP,xi,yj,zk)
    type(TScalFlIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk               !grid coordinates of the forcing point
    type(SolidBody),pointer :: SB
    integer dirx,diry,dirz,dirx2,diry2,dirz2,nfreedirs,ndirs,i
    real(knd) x,y,z,xnear,ynear,znear,distx,disty,distz,t,tx,ty,tz
    logical freep00,free0p0,free00p,freem00,free0m0,free00m

    x = xPr(xi)                                   !physical coordinates of the forcing point
    y = yPr(yj)
    z = zPr(zk)
    call SetCurrentSB(SB,Prtype(xi,yj,zk))
    call SB%ClosestOut(xnear,ynear,znear,x,y,z)

    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk

    distx = xnear-x                               !distance to the nearest point on the boundary
    disty = ynear-y
    distz = znear-z
    dirx = nint(sign(1.0_knd,distx))              !direction to the boundary point
    diry = nint(sign(1.0_knd,disty))
    dirz = nint(sign(1.0_knd,distz))


    if (abs(distx)<(xPr(xi+1)-xPr(xi-1))/100._knd) then  !If very close to the boundary, set the boundary here
      distx = 0
      dirx = 0
    end if

    if (abs(disty)<(yPr(yj+1)-yPr(yj-1))/100._knd) then
      disty = 0
      diry = 0
    end if

    if (abs(distz)<(zW(zk+1)-zW(zk-1))/100._knd) then
      distz = 0
      dirz = 0
    end if


    freep00 = (Prtype(xi+1,yj,zk)<=0)   !logicals denoting if the cell in plus x direction is free of SB
    freem00 = (Prtype(xi-1,yj,zk)<=0)
    free0p0 = (Prtype(xi,yj+1,zk)<=0)
    free0m0 = (Prtype(xi,yj-1,zk)<=0)
    free00p = (Prtype(xi,yj,zk+1)<=0)
    free00m = (Prtype(xi,yj,zk-1)<=0)



    if (.not.(freep00.or.freem00)) then
      dirx = 0
      distx = 0
    end if
    if (.not.(free0p0.or.free0m0)) then
      diry = 0
      disty = 0
    end if
    if (.not.(free00p.or.free00m)) then
      dirz = 0
      distz = 0
    end if

    nfreedirs = 0
    ndirs = 0

    if (freep00) nfreedirs = nfreedirs+1
    if (free0p0) nfreedirs = nfreedirs+1
    if (free00p) nfreedirs = nfreedirs+1
    if (freem00) nfreedirs = nfreedirs+1
    if (free0m0) nfreedirs = nfreedirs+1
    if (free00m) nfreedirs = nfreedirs+1

    if (dirx/=0) ndirs = ndirs+1
    if (diry/=0) ndirs = ndirs+1
    if (dirz/=0) ndirs = ndirs+1

    if (nfreedirs>ndirs.or.(freep00.and.freem00).or.(free0p0.and.free0m0).or.(free00p.and.free00m)) then
      IBP%interp = 0    !If more free spaces than directions, or free space in both oposite directions,
      distx = 0         !Assume is an unresolvably small feature and treat as on the boundary.
      disty = 0
      distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      dirx = 0
      diry = 0
      dirz = 0
    end if


    if (dirx==0.and.diry==0.and.dirz==0) then   !If dir>0 treat as on the boundary.

      dirx2 = 0
      diry2 = 0
      dirz2 = 0
      if (freep00) dirx2 = dirx2+1     !Find a free neighbouring cell and interpolate from there.
      if (freem00) dirx2 = dirx2-1
      if (free0p0) diry2 = diry2+1
      if (free0m0) diry2 = diry2-1
      if (free00p) dirz2 = dirz2+1
      if (free00m) dirz2 = dirz2-1

      if (dirx2==0.and.diry2==0.and.dirz2==0) then  !probably a very thin body, just average free points around

        IBP%interp = nfreedirs

        allocate(IBP%IntPoints(nfreedirs))

        i = 1
        IBP%dist = 0
        if (freep00) then
          IBP%IntPoints(i) = TInterpolationPoint(xi + 1, yj, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + xPr(xi+1)-xPr(xi)
        end if
        if (freem00) then
          IBP%IntPoints(i) = TInterpolationPoint(xi - 1, yj, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + xPr(xi)-xPr(xi-1)
        end if
        if (free0p0) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj + 1, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + yPr(yj+1)-yPr(yj)
        end if
        if (free0m0) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj - 1, zk, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + yPr(yj)-yPr(yj-1)
        end if
        if (free00p) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj, zk + 1, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + zPr(zk+1)-zPr(zk)
        end if
        if (free00m) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj, zk - 1, 1._knd / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + zPr(zk)-zPr(zk-1)
        end if
        IBP%dist = IBP%dist / nfreedirs

      else  !some edge

        allocate(IBP%IntPoints(1))

        IBP%IntPoints(1)%xi = IBP%xi+dirx2
        IBP%IntPoints(1)%yj = IBP%yj+diry2
        IBP%IntPoints(1)%zk = IBP%zk+dirz2
        IBP%IntPoints%coef = 1._knd
        IBP%interp = 1
        IBP%dist = sqrt((x-xPr(IBP%xi+dirx2))**2+(y-yPr(IBP%yj+diry2))**2+(z-zPr(IBP%zk+dirz2))**2)

      end if

    elseif (nfreedirs==1.or.ndirs==1) then     !Only one free point outside, interpolate from there.

      allocate(IBP%IntPoints(1))

      IBP%IntPoints(1)%xi = IBP%xi+dirx
      IBP%IntPoints(1)%yj = IBP%yj+diry
      IBP%IntPoints(1)%zk = IBP%zk+dirz
      IBP%IntPoints%coef = 1._knd
      IBP%interp = 1
      IBP%dist = sqrt((x-xPr(IBP%xi+dirx))**2+(y-yPr(IBP%yj+diry))**2+(z-zPr(IBP%zk+dirz))**2)

    elseif (nfreedirs==2.or.ndirs==2) then    !Two free directions, use bilinear interpolation in the plane contaning these two neigbours.

      allocate(IBP%IntPoints(2))

      if (dirx==0) then     !plane yz

        if (abs(disty/distz)<abs(yPr(yj+diry)-y)/abs(zPr(zk+dirz)-z)) then  !Which gridline does the line from this point
          t = (zPr(zk+dirz)-z)/distz                                           !to the boundary intersect? Then use the points
          IBP%IntPoints(1)%xi = IBP%xi
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk+dirz
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
        else
          t = (yPr(yj+diry)-y)/disty
          IBP%IntPoints(1)%xi = IBP%xi
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj+diry
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
        end if

      elseif (diry==0) then     !plane xz

        if (abs(distx/distz)<abs(xPr(xi+dirx)-x)/abs(zPr(zk+dirz)-z)) then
          t = (zPr(zk+dirz)-z)/distz
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk+dirz
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
        else
          t = (xPr(xi+dirx)-x)/distx
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
          IBP%IntPoints(2)%xi = IBP%xi+dirx
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
        end if

      else                 !plane xy

        if (abs(distx/disty)<abs(xPr(xi+dirx)-x)/abs(yPr(yj+diry)-y)) then
          t = (yPr(yj+diry)-y)/disty
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk
          IBP%IntPoints(1)%coef = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj+diry
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
        else
          t = (xPr(xi+dirx)-x)/distx
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk
          IBP%IntPoints(1)%coef = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
          IBP%IntPoints(2)%xi = IBP%xi+dirx
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
        end if
      end if

    else                         !more than two free directions

      tx = (xPr(xi+dirx)-x)/distx   !a vector pointing from a free point to this forcing point.
      ty = (yPr(yj+diry)-y)/disty
      tz = (zPr(zk+dirz)-z)/distz

      IBP%interp = 4

      allocate(IBP%IntPoints(4))

      if (tx<=ty.and.tx<=tz) then
        !coordinates of interpolation point are therefore
        !xPr(xi+dirx) = x+tx*distx
        !y+tx*disty
        !z+tx*distz
        IBP%dist = sqrt((distx*tx)**2+(disty*tx)**2+(distz*tx)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                  abs(distz*tx)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(2)%xi = IBP%xi+dirx
        IBP%IntPoints(2)%yj = IBP%yj
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                  abs(distz*tx)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj+diry
        IBP%IntPoints(3)%zk = IBP%zk
        IBP%IntPoints(3)%coef = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                  (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
        IBP%IntPoints(4)%xi = IBP%xi+dirx
        IBP%IntPoints(4)%yj = IBP%yj
        IBP%IntPoints(4)%zk = IBP%zk
        IBP%IntPoints(4)%coef = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                  (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
      elseif (ty<=tx.and.ty<=tz) then
        !the same with ty
        IBP%dist = sqrt((distx*ty)**2+(disty*ty)**2+(distz*ty)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                  abs(distz*ty)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(2)%xi = IBP%xi
        IBP%IntPoints(2)%yj = IBP%yj+diry
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                  abs(distz*ty)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj+diry
        IBP%IntPoints(3)%zk = IBP%zk
        IBP%IntPoints(3)%coef = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                  (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
        IBP%IntPoints(4)%xi = IBP%xi
        IBP%IntPoints(4)%yj = IBP%yj+diry
        IBP%IntPoints(4)%zk = IBP%zk
        IBP%IntPoints(4)%coef = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                  (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
      else
        !the same with tz
        IBP%dist = sqrt((distx*tz)**2+(disty*tz)**2+(distz*tz)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                  abs(disty*tz)/abs(yPr(yj+diry)-y)
        IBP%IntPoints(2)%xi = IBP%xi
        IBP%IntPoints(2)%yj = IBP%yj+diry
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                  abs(disty*tz)/abs(yPr(yj+diry)-y)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj
        IBP%IntPoints(3)%zk = IBP%zk+dirz
        IBP%IntPoints(3)%coef = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                  (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
        IBP%IntPoints(4)%xi = IBP%xi
        IBP%IntPoints(4)%yj = IBP%yj
        IBP%IntPoints(4)%zk = IBP%zk+dirz
        IBP%IntPoints(4)%coef = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                  (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
      end if

    end if

    IBP%temperature_flux = SB%temperature_flux

  end subroutine TScalFlIBPoint_Create


  subroutine MoveIBPointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
    integer :: i

    allocate(UIBPoints(UIBPointsList%Len()))
    allocate(VIBPoints(VIBPointsList%Len()))
    allocate(WIBPoints(WIBPointsList%Len()))
    allocate(ScalFlIBPoints(ScalFlIBPointsList%Len()))

    i = 0
    call UIBPointsList%for_each(CopyVelIBPoint)
    i = 0
    call VIBPointsList%for_each(CopyVelIBPoint)
    i = 0
    call WIBPointsList%for_each(CopyVelIBPoint)
    i = 0
    call ScalFlIBPointsList%for_each(CopyScalFlIBPoint)

    contains

      subroutine CopyVelIBPoint(CurrentIBPoint)
        type(TVelIBPoint) :: CurrentIBPoint

        i = i + 1

        if (CurrentIBPoint%component==1) then
          UIBPoints(i) = CurrentIBPoint
        elseif (CurrentIBPoint%component==2) then
          VIBPoints(i) = CurrentIBPoint
        else
          WIBPoints(i) = CurrentIBPoint
        end if

      end subroutine

      subroutine CopyScalFlIBPoint(CurrentIBPoint)
        type(TScalFlIBPoint) :: CurrentIBPoint

        i = i + 1

        ScalFlIBPoints(i) = CurrentIBPoint

      end subroutine

  end subroutine MoveIBPointsToArray


  subroutine AuxNeighbours(Xtype,nx,ny,nz,s)
    !helper procedure to FindNeighbouringCells
    integer,intent(in)    :: nx,ny,nz,s
    integer,intent(inout) :: Xtype(s:,s:,s:)
    integer i,j,k
!     if (Btype(We)==DIRICHLET.or.Btype(We)==NOSLIP) Xtype(:0,:,:) = -2
!     if (Btype(Ea)==DIRICHLET.or.Btype(Ea)==NOSLIP) Xtype(nx+1:,:,:) = -2
!     if (Btype(So)==DIRICHLET.or.Btype(So)==NOSLIP) Xtype(:,:0,:) = -2
!     if (Btype(No)==DIRICHLET.or.Btype(No)==NOSLIP) Xtype(:,ny+1:,:) = -2
!     if (Btype(Bo)==DIRICHLET.or.Btype(Bo)==NOSLIP) Xtype(:,:,:0) = -2
!     if (Btype(To)==DIRICHLET.or.Btype(To)==NOSLIP) Xtype(:,:,nz+1:) = -2

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (Xtype(i,j,k)==0.and. &
              (any(Xtype(i-1:i+1,j-1:j+1,k-1:k+1)>0) )) &!.or. &
!                any(Xtype(i-1:i+1,j-1:j+1,k-1:k+1)<-1) ) ) &
                                                        Xtype(i,j,k) = -1
        end do
      end do
    end do

  end subroutine

  subroutine FindNeighbouringCells
    !sets type of the cells closest to the solid body to -1
    call AuxNeighbours(Prtype,Prnx,Prny,Prnz,0)
    call AuxNeighbours(Utype,Unx,Uny,Unz,-2)
    call AuxNeighbours(Vtype,Vnx,Vny,Vnz,-2)
    call AuxNeighbours(Wtype,Wnx,Wny,Wnz,-2)
  end subroutine FindNeighbouringCells

  
  integer function FindScalarIBCellIndex(xi,yj,zk) result(res)
    !binary search of the immersed boundary cell index with the prescribed coordinates
    integer,intent(in) :: xi,yj,zk
    integer :: idx,d,u,i,idx_d,idx_u,idx_i
    
    idx = f_ijk(xi,yj,zk)
    
    d = 1
    u = size(ScalFlIBPoints)
    
    idx_d = f(d)
    if (idx_d==idx) then
      res = d
      return
    end if
    
    idx_u = f(u)
    if (idx_u==idx) then
      res = u
      return
    end if
    
    res = -1
    
    do while (u-d>1)
      i = (d+u)/2
      idx_i = f(i)
      if (idx_i == idx) then
        res = i
        return
      else if (idx_i<idx) then
        d = i
      else
        u = i
      end if
    end do
    
    contains
      
      pure integer function f_ijk(i,j,k)
        integer,intent(in) :: i,j,k
        f_ijk = (i-1) + Prnx * (j-1) + Prnx * Prny * (k-1)
      end function
      pure integer function f(n)
        integer,intent(in) :: n
        f = (ScalFlIBPoints(n)%xi-1) + &
            Prnx * (ScalFlIBPoints(n)%yj-1) + &
            Prnx * Prny * (ScalFlIBPoints(n)%zk-1)
      end function
  end function FindScalarIBCellIndex
  
!   
!   subroutine BindWMstoIBPs
!     use ArrayUtilities, only: add_element
!     use WallModels, only: WMPoints
!     integer :: neighbours(3,6)
!     integer :: i, j
!     integer :: xn, yn, zn, idx
!     
!     neighbours = 0
!     neighbours(1,1) =  1
!     neighbours(1,2) = -1
!     neighbours(2,3) =  1
!     neighbours(2,4) = -1
!     neighbours(3,5) =  1
!     neighbours(3,6) = -1
!     
!     do i=1,size(WMPoints)
!       associate (p => WMPoints(i)) !gfortran 4.8 bug prevents more items here
!         allocate(p%bound_IBPs(0))
! 
!         do j=1,6
!           xn = p%xi+neighbours(1,j)
!           yn = p%yj+neighbours(2,j)
!           zn = p%zk+neighbours(3,j)
!           if (Prtype(xn,yn,zn)>0.and. &
!               xn>0 .and. xn<=Prnx .and. &
!               yn>0 .and. yn<=Prny .and. &
!               zn>0 .and. zn<=Prnz) then
! 
!             idx = FindScalarIBCellIndex(xn, &
!                                         yn, &
!                                         zn)
!             if (idx>0) then
!               ScalFlIBPoints(idx)%n_WMPs = ScalFlIBPoints(idx)%n_WMPs + 1
!               call add_element(p%bound_IBPs, idx)
!             else
!               stop "This point does not have an IBP!"
!             end if
!           end if
!         end do
!         
!       end associate
!     end do
!   end subroutine BindWMstoIBPs
  
  
  subroutine InitIBPFluxes
  
    !find the corresponding immersed boundary points to wall model points
!     call BindWMstoIBPs
    !compute initial temperature and moisture fluxes
    call GetWMFluxes
  end subroutine
  
!   subroutine SetIBPFluxes
!     use WallModels, only: WMPoints
!     integer i,j
!          
!     ScalFlIBPoints%temperature_flux = 0
! 
!     do i=1,size(WMPoints)
!       associate (p => WMPoints(i))
!         if (enable_buoyancy) then
!           do j=1,size(p%bound_IBPs)
!             associate (IBP => ScalFlIBPoints(p%bound_IBPs(j)))
!               IBP%temperature_flux = IBP%temperature_flux + p%temperature_flux/IBP%n_WMPs
!             end associate
!           end do
!         end if
!         
!         if (enable_moisture) then
!           do j=1,size(p%bound_IBPs)
!             associate (IBP => ScalFlIBPoints(p%bound_IBPs(j)))
!               IBP%moisture_flux = IBP%moisture_flux + p%moisture_flux/IBP%n_WMPs
!             end associate
!           end do
!         end if
!         
!         if (p%zk==1) then
!           if (enable_buoyancy) BsideTFlArr(p%xi,p%yj) = p%temperature_flux
!           if (enable_moisture) BsideMFlArr(p%xi,p%yj) = p%moisture_flux
!         end if
!       end associate
!     end do
!   
!   end subroutine SetIBPFluxes
!   
  


  subroutine InitImBoundaries
    type(TVelIBPoint) IBP
    type(TScalFlIBPoint) SIBP
    integer i,j,k


    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)>0) then
        if (Utype(i+1,j,k)<=0.or.Utype(i-1,j,k)<=0.or.Utype(i,j+1,k)<=0&
          .or.Utype(i,j-1,k)<=0.or.Utype(i,j,k+1)<=0.or.Utype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xU(-2:),yPr(-2:),zPr(-2:),Utype,1)
            call  UIBPointsList%add(IBP)
        end if
       end if
      end do
     end do
    end do

    do k = 1,Vnz
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)>0) then
        if (Vtype(i+1,j,k)<=0.or.Vtype(i-1,j,k)<=0.or.Vtype(i,j+1,k)<=0&
          .or.Vtype(i,j-1,k)<=0.or.Vtype(i,j,k+1)<=0.or.Vtype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xPr(-2:),yV(-2:),zPr(-2:),Vtype,2)
            call  VIBPointsList%add(IBP)
        end if
       end if
      end do
     end do
    end do

    do k = 1,Wnz
     do j = 1,Wny
      do i = 1,Wnx
       if (Wtype(i,j,k)>0) then
        if (Wtype(i+1,j,k)<=0.or.Wtype(i-1,j,k)<=0.or.Wtype(i,j+1,k)<=0&
          .or.Wtype(i,j-1,k)<=0.or.Wtype(i,j,k+1)<=0.or.Wtype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xPr(-2:),yPr(-2:),zW(-2:),Wtype,3)
            call  WIBPointsList%add(IBP)
        end if
       end if
      end do
     end do
    end do

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)>0) then
        if (Prtype(i+1,j,k)<=0.or.Prtype(i-1,j,k)<=0.or.Prtype(i,j+1,k)<=0&
          .or.Prtype(i,j-1,k)<=0.or.Prtype(i,j,k+1)<=0.or.Prtype(i,j,k-1)<=0)  then
            call  Create(SIBP,i,j,k)
            call  ScalFlIBPointsList%add(SIBP)
        end if
       end if
      end do
     end do
    end do

    call MoveIBPointsToArray

    call UIBPointsList%finalize
    call VIBPointsList%finalize
    call WIBPointsList%finalize
    call ScalFlIBPointsList%finalize

  end subroutine InitImBoundaries


  subroutine GetSolidBodiesBC

    call FindInsideCells

    call FindNeighbouringCells

    call InitImBoundaries

    call GetSolidBodiesWM

    call GetSolidBodiesWM_UVW

  end subroutine GetSolidBodiesBC


end module ImmersedBoundary
