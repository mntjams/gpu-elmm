    integer :: xi, dii, dij, dik
    integer :: prev

    if (enable_multiple_domains) then

      if (domain_index==1) then

        !send buffers should be counted from 1, there may be more of them from several nested domains
        allocate(domain_bc_send_buffers(2*domain_nyims(2)*domain_nzims(2)))! + &
!                                         2*domain_nxims(2)*domain_nzims(2) + &
!                                         domain_nxims(2)*domain_nyims(2)))

        prev = 0

        do dik = 1, domain_nzims(2)
          do dij = 1, domain_nyims(2)
            call parent_buffer(num=dij+(dik-1)*domain_nyims(2), &
                               dir=We, domain=2, im=[1,dij,dik])
         end do
        end do

        prev = prev + domain_nyims(2)*domain_nzims(2)

        do dik = 1, domain_nzims(2)
          do dij = 1, domain_nyims(2)
            call parent_buffer(num=dij+(dik-1)*domain_nyims(2) + prev, &
                               dir=Ea, domain=2, im=[1,dij,dik])
         end do
        end do

!         prev = prev + domain_nyims(2)*domain_nzims(2)
! 
!         do dik = 1, domain_nzims(2)
!           do dii = 1, domain_nxims(2)
!             call parent_buffer(num=dii+(dik-1)*domain_nxims(2) + prev, &
!                                dir=So, domain=2, im=[dii,1,dik])
!          end do
!         end do
! 
!         prev = prev + domain_nxims(2)*domain_nzims(2)
! 
!         do dik = 1, domain_nzims(2)
!           do dii = 1, domain_nxims(2)
!             call parent_buffer(num=dii+(dik-1)*domain_nxims(2) + prev, &
!                                dir=No, domain=2, im=[dii,domain_nyims(2),dik])
!          end do
!         end do
! 
!         prev = prev + domain_nxims(2)*domain_nzims(2)
! 
!         do dij = 1, domain_nyims(2)
!           do dii = 1, domain_nxims(2)
!             call parent_buffer(num=dii+(dij-1)*domain_nxims(2) + prev, &
!                                dir=To, domain=2, im=[dii,dij,domain_nzims(2)])
! 
!          end do
!         end do

      else if (domain_index==2) then

        allocate(domain_bc_recv_buffers(We:To))

        if (iim==1) then
          Btype(We) = BC_DOMAIN_COPY
          TempBtype(We) = BC_DOMAIN_COPY
          MoistBtype(We) = BC_DOMAIN_COPY
          ScalBtype(We) = BC_DOMAIN_COPY

          call child_buffer(dir=We, domain=1)
        endif

        if (iim==nxims) then
          Btype(Ea) = BC_DOMAIN_COPY
          TempBtype(Ea) = BC_DOMAIN_COPY
          MoistBtype(Ea) = BC_DOMAIN_COPY
          ScalBtype(Ea) = BC_DOMAIN_COPY

          call child_buffer(dir=Ea, domain=1)
        end if
! 
!         if (jim==1) then
!           Btype(So) = BC_DOMAIN_COPY
!           TempBtype(So) = BC_DOMAIN_COPY
!           MoistBtype(So) = BC_DOMAIN_COPY
!           ScalBtype(So) = BC_DOMAIN_COPY
! 
!           call child_buffer(dir=So, domain=1)
!         endif
! 
!         if (jim==nyims) then
!           Btype(No) = BC_DOMAIN_COPY
!           TempBtype(No) = BC_DOMAIN_COPY
!           MoistBtype(No) = BC_DOMAIN_COPY
!           ScalBtype(No) = BC_DOMAIN_COPY
! 
!           call child_buffer(dir=No, domain=1)
!         end if
! 
!         if (kim==nzims) then
!           Btype(To) = BC_DOMAIN_COPY
!           TempBtype(To) = BC_DOMAIN_COPY
!           MoistBtype(To) = BC_DOMAIN_COPY
!           ScalBtype(To) = BC_DOMAIN_COPY
! 
!           call child_buffer(dir=To, domain=1)
!         end if

      end if


    end if !enable_multiple_domains


contains

    function x_coord(x) result(res)
      integer :: res
      real(knd), intent(in) :: x
      res = min( max(nint( (x - im_xmin)/dxmin ),0) , Unx+1)
    end function

    function y_coord(y) result(res)
      integer :: res
      real(knd), intent(in) :: y
      res = min( max(nint( (y - im_ymin)/dymin ),0) , Vny+1)
    end function

    function z_coord(z) result(res)
      integer :: res
      real(knd), intent(in) :: z
      res = min( max(nint( (z - im_zmin)/dzmin ),0) , Wnz+1)
    end function

    subroutine parent_buffer(num, dir, domain, im)
      integer, intent(in) :: num, dir, domain, im(3)

      real(knd) :: cxmax, cxmin, cymax, cymin, czmax, czmin
      integer :: cxi1, cxi2, cyj1, cyj2, czk1, czk2
      integer :: pos

      associate (b=>domain_bc_send_buffers(num))

        b%comm = world_comm
        b%remote_rank = domain_ranks_grid(domain)%arr(im(1),im(2),im(3))
        b%remote_domain = domain
        b%direction = dir
        b%enabled = .true.

        !get the child image extent
        cxmin = domain_grids(domain)%xmins(im(1))
        cxmax = domain_grids(domain)%xmaxs(im(1))
        cymin = domain_grids(domain)%ymins(im(2))
        cymax = domain_grids(domain)%ymaxs(im(2))
        czmin = domain_grids(domain)%zmins(im(3))
        czmax = domain_grids(domain)%zmaxs(im(3))

        !find the child grid indexes
        cxi1 = x_coord(cxmin)
        cxi2 = x_coord(cxmax)
        cyj1 = y_coord(cymin)
        cyj2 = y_coord(cymax)
        czk1 = z_coord(czmin)
        czk2 = z_coord(czmax)

        select case (dir)
          case (We)
            pos = cxi1
          case (Ea)
            pos = cxi2
          case (So)
            pos = cyj1
          case (No)
            pos = cyj2
          case (Bo)
            pos = czk1
          case (To)
            pos = czk2
        end select

        if (dir==We .or. dir==Ea) then

          if (dir==We) then
            Ui1 = pos-2
            Ui2 = pos+2

            Vi1 = pos-2
            Vi2 = pos+3

            Wi1 = pos-2
            Wi2 = pos+3

            Pri1 = pos-2
            Pri2 = pos+3
          else if (dir==Ea) then
            Ui1 = pos-2
            Ui2 = pos+2

            Vi1 = pos-2
            Vi2 = pos+3

            Wi1 = pos-2
            Wi2 = pos+3

            Pri1 = pos-2
            Pri2 = pos+2
          end if

          Uj1 = cyj1-2
          Uj2 = cyj2+3
          Uk1 = czk1-2
          Uk2 = czk2+3

          Vj1 = cyj1-2
          if (domain_grids(domain)%internal_or_periodic_No(im(2))) then
            Vj2 = cyj2+3
          else
            Vj2 = cyj2+2
          end if
          Vk1 = czk1-2
          Vk2 = czk2+3

          Wj1 = cyj1-2
          Wj2 = cyj2+3
          Wk1 = czk1-2
          if (domain_grids(domain)%internal_or_periodic_To(im(3))) then
            Wk2 = czk2+3
          else
            Wk2 = czk2+2
          end if

          Prj1 = cyj1-1
          Prj2 = cyj2+2
          Prk1 = czk1-1
          Prk2 = czk2+2

        else if (dir==So .or. dir==No) then

          if (dir==So) then
            Uj1 = pos-2
            Uj2 = pos+2

            Vj1 = pos-2
            Vj2 = pos+2

            Wj1 = pos-2
            Wj2 = pos+2

            Prj1 = pos-1
            Prj2 = pos+2
          else if (dir==No) then
            Uj1 = pos-1
            Uj2 = pos+3

            Vj1 = pos-2
            Vj2 = pos+2

            Wj1 = pos-1
            Wj2 = pos+3

            Prj1 = pos-1
            Prj2 = pos+2
          end if

          Ui1 = cxi1-2
          if (domain_grids(domain)%internal_or_periodic_Ea(im(1))) then
            Ui2 = cxi2+3
          else
            Ui2 = cxi2+2
          end if
          Uk1 = czk1-2
          Uk2 = czk2+3

          Vi1 = cxi1-2
          Vi2 = cxi2+3
          Vk1 = czk1-2
          Vk2 = czk2+3

          Wi1 = cxi1-2
          Wi2 = cxi2+3
          Wk1 = czk1-2
          if (domain_grids(domain)%internal_or_periodic_To(im(3))) then
            Wk2 = czk2+3
          else
            Wk2 = czk2+2
          end if

          Pri1 = cxi1-1
          Pri2 = cxi2+2
          Prk1 = czk1-1
          Prk2 = czk2+2

        else

          if (dir==Bo) then
            Uk1 = pos-2
            Uk2 = pos+2

            Vk1 = pos-2
            Vk2 = pos+2

            Wk1 = pos-2
            Wk2 = pos+2

            Prk1 = pos-1
            Prk2 = pos+2
          else if (dir==To) then
            Uk1 = pos-1
            Uk2 = pos+3

            Vk1 = pos-1
            Vk2 = pos+3

            Wk1 = pos-2
            Wk2 = pos+2

            Prk1 = pos-1
            Prk2 = pos+2
          end if

          Ui1 = cxi1-2
          if (domain_grids(domain)%internal_or_periodic_Ea(im(1))) then
            Ui2 = cxi2+3
          else
            Ui2 = cxi2+2
          end if
          Uj1 = cyj1-2
          Uj2 = cyj2+3

          Vi1 = cxi1-2
          Vi2 = cxi2+3
          Vj1 = cyj1-2
          if (domain_grids(domain)%internal_or_periodic_No(im(2))) then
            Vj2 = cyj2+3
          else
            Vj2 = cyj2+2
          end if

          Wi1 = cxi1-2
          Wi2 = cxi2+3
          Wj1 = cyj1-2
          Wj2 = cyj2+3

          Pri1 = cxi1-1
          Pri2 = cxi2+2
          Prj1 = cyj1-1
          Prj2 = cyj2+2

        end if

        b%position = pos

        b%Ui1 = Ui1 
        b%Ui2 = Ui2 
        b%Uj1 = Uj1 
        b%Uj2 = Uj2 
        b%Uk1 = Uk1 
        b%Uk2 = Uk2 

        b%Vi1 = Vi1 
        b%Vi2 = Vi2 
        b%Vj1 = Vj1 
        b%Vj2 = Vj2 
        b%Vk1 = Vk1 
        b%Vk2 = Vk2 

        b%Wi1 = Wi1 
        b%Wi2 = Wi2 
        b%Wj1 = Wj1 
        b%Wj2 = Wj2 
        b%Wk1 = Wk1 
        b%Wk2 = Wk2 

        b%Pri1 = Pri1
        b%Pri2 = Pri2
        b%Prj1 = Prj1
        b%Prj2 = Prj2
        b%Prk1 = Prk1
        b%Prk2 = Prk2


        allocate(b%U(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%V(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%W(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        allocate(b%dU_dt(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%dV_dt(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%dW_dt(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        if (enable_buoyancy) then
          allocate(b%Temperature(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dTemperature_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if
        if (enable_moisture) then
          allocate(b%Moisture(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dMoisture_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if

      end associate

    end subroutine parent_buffer


    subroutine child_buffer(dir, domain)
      integer, intent(in) :: dir, domain
      integer :: i, width

      associate(b => domain_bc_recv_buffers(dir))

        b%comm = world_comm
        b%remote_rank = domain_ranks_grid(1)%arr(1,1,1) !arr(domain_nxims(1),jim,kim)        
        b%remote_domain = domain
        b%direction = dir
        b%enabled = .true.

        b%spatial_ratio = domain_spatial_ratio
        b%time_step_ratio = domain_time_step_ratio

        width = b%spatial_ratio * 2

        select case  (dir)
          case (We)
            b%position = 0

            Ui1 = -2
            Ui2 = width
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = width
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = width
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = width
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = +2
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = +3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = +3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = -2
            b%r_Pri2 = +3
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (Ea)
            b%position = Prnx

            Ui1 = Unx-width+1
            Ui2 = Prnx+2
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = Vnx-width+1
            Vi2 = Prnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = Wnx-width+1
            Wi2 = Prnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = Prnx-width+1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = Prnx/b%spatial_ratio-2
            b%r_Ui2 = Prnx/b%spatial_ratio+2
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = Prnx/b%spatial_ratio-2
            b%r_Vi2 = Prnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = Prnx/b%spatial_ratio-2
            b%r_Wi2 = Prnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = Prnx/b%spatial_ratio-2
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (So)
            b%position = 0

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = +2
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = +2
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = +2
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = +2
            Prk1 = -1
            Prk2 = Prnz+2

          case (No)
            b%position = Prny

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = Prny-1
            Uj2 = Prny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = Prny-2
            Vj2 = Prny+2
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = Prny-1
            Wj2 = Prny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = Prny-1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

          case (Bo)
            b%position = 0
          case (To)
            b%position = Prnz

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = Prnz-1
            Uk2 = Prnz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = Prnz-1
            Vk2 = Prnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = Prnz-2
            Wk2 = Prnz+2

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = Prnz-1
            Prk2 = Prnz+2

        end select    


        b%Ui1 = Ui1 
        b%Ui2 = Ui2 
        b%Uj1 = Uj1 
        b%Uj2 = Uj2 
        b%Uk1 = Uk1 
        b%Uk2 = Uk2 

        b%Vi1 = Vi1 
        b%Vi2 = Vi2 
        b%Vj1 = Vj1 
        b%Vj2 = Vj2 
        b%Vk1 = Vk1 
        b%Vk2 = Vk2 

        b%Wi1 = Wi1 
        b%Wi2 = Wi2 
        b%Wj1 = Wj1 
        b%Wj2 = Wj2 
        b%Wk1 = Wk1 
        b%Wk2 = Wk2 

        b%Pri1 = Pri1
        b%Pri2 = Pri2
        b%Prj1 = Prj1
        b%Prj2 = Prj2
        b%Prk1 = Prk1
        b%Prk2 = Prk2


        b%bUi1 = Ui1 
        b%bUi2 = Ui2 
        b%bUj1 = Uj1 
        b%bUj2 = Uj2 
        b%bUk1 = Uk1 
        b%bUk2 = Uk2 

        b%bVi1 = Vi1 
        b%bVi2 = Vi2 
        b%bVj1 = Vj1 
        b%bVj2 = Vj2 
        b%bVk1 = Vk1 
        b%bVk2 = Vk2 

        b%bWi1 = Wi1 
        b%bWi2 = Wi2 
        b%bWj1 = Wj1 
        b%bWj2 = Wj2 
        b%bWk1 = Wk1 
        b%bWk2 = Wk2 

        b%bPri1 = Pri1
        b%bPri2 = Pri2
        b%bPrj1 = Prj1
        b%bPrj2 = Prj2
        b%bPrk1 = Prk1
        b%bPrk2 = Prk2

        select case  (dir)
          case (We)
            b%bUi2 = Ui2 - width 
            b%bVi2 = Vi2 - width 
            b%bWi2 = Wi2 - width 
            b%bPri2 = Pri2 - width 
          case (Ea)
            b%bUi1 = Ui1 + width 
            b%bVi1 = Vi1 + width 
            b%bWi1 = Wi1 + width 
            b%bPri1 = Pri1 + width 
          case (So)
            b%bUj2 = Uj2 - width 
            b%bVj2 = Vj2 - width 
            b%bWj2 = Wj2 - width 
            b%bPrj2 = Prj2 - width 
          case (No)
            b%bUj1 = Uj1 + width 
            b%bVj1 = Vj1 + width 
            b%bWj1 = Wj1 + width 
            b%bPrj1 = Prj1 + width 
          case (Bo)
          case (To)
            b%bUk1 = Uk1 + width 
            b%bVk1 = Vk1 + width 
            b%bWk1 = Wk1 + width 
            b%bPrk1 = Prk1 + width 
        end select

         
        b%r_x0 = im_xmin
        b%r_y0 = im_ymin
        b%r_z0 = im_zmin
        b%r_dx = domain_grids(domain)%dx
        b%r_dy = domain_grids(domain)%dy
        b%r_dz = domain_grids(domain)%dz
        
        allocate(b%r_xU(b%r_Ui1:b%r_Ui2))
        allocate(b%r_yV(b%r_Vj1:b%r_Vj2))
        allocate(b%r_zW(b%r_Wk1:b%r_Wk2))
        allocate(b%r_x(b%r_Vi1:b%r_Vi2))
        allocate(b%r_y(b%r_Uj1:b%r_Uj2))
        allocate(b%r_z(b%r_Uk1:b%r_Uk2))
        b%r_xU = [(b%r_x0 + i*b%r_dx , i = b%r_Ui1,b%r_Ui2)]
        b%r_yV = [(b%r_y0 + i*b%r_dy , i = b%r_Vj1,b%r_Vj2)]
        b%r_zW = [(b%r_z0 + i*b%r_dz , i = b%r_Wk1,b%r_Wk2)]
        b%r_x = [(b%r_x0 + (i-0.5_knd)*b%r_dx , i = b%r_Vi1,b%r_Vi2)]
        b%r_y = [(b%r_y0 + (i-0.5_knd)*b%r_dy , i = b%r_Uj1,b%r_Uj2)]
        b%r_z = [(b%r_z0 + (i-0.5_knd)*b%r_dz , i = b%r_Vk1,b%r_Vk2)]


        allocate(b%U(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%V(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%W(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        allocate(b%dU_dt(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
        allocate(b%dV_dt(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
        allocate(b%dW_dt(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

        if (enable_buoyancy) then
          allocate(b%Temperature(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dTemperature_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if
        if (enable_moisture) then
          allocate(b%Moisture(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          allocate(b%dMoisture_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
        end if

        allocate(b%r_U(b%r_Ui1:b%r_Ui2,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        allocate(b%r_V(b%r_Vi1:b%r_Vi2,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        allocate(b%r_W(b%r_Wi1:b%r_Wi2,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        allocate(b%r_dU_dt(b%r_Ui1:b%r_Ui2,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        allocate(b%r_dV_dt(b%r_Vi1:b%r_Vi2,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        allocate(b%r_dW_dt(b%r_Wi1:b%r_Wi2,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        if (enable_buoyancy) then
          allocate(b%r_Temperature(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
          allocate(b%r_dTemperature_dt(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
        end if
        if (enable_moisture) then
          allocate(b%r_Moisture(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
          allocate(b%r_dMoisture_dt(b%r_Pri1:b%r_Pri2,b%r_Prj1:b%r_Prj2,b%r_Prk1:b%r_Prk2))
        end if

      end associate

    end subroutine child_buffer
