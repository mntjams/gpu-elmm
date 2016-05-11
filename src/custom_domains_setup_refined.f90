    integer :: xi, dii, djj, dkk
    integer :: prev
    integer :: child_domain, i_child_domain, side
    integer :: n_send_buffers
    integer :: err
    integer, allocatable :: requests(:)
    integer :: request

    if (enable_multiple_domains) then

      allocate(requests(0))

      allocate(domain_parent_buffers(size(child_domains)))

      !send buffers should be counted from 1, there may be more of them from several nested domains
      n_send_buffers = 0
      do i_child_domain = 1, size(child_domains)
        child_domain = child_domains(i_child_domain)

        domain_parent_buffers(i_child_domain)%remote_domain = child_domain

        allocate(domain_parent_buffers(i_child_domain)%bs(domain_nxims(child_domain), &
                                                          domain_nyims(child_domain), &
                                                          domain_nzims(child_domain)))

        do dkk = 1, domain_nzims(child_domain)
          do djj = 1, domain_nyims(child_domain)
            do dii = 1, domain_nxims(child_domain)
              call create_domain_parent_buffer(domain_parent_buffers(i_child_domain)%bs(dii,djj,dkk), &
                                               child_domain, [dii, djj, dkk])
            end do
          end do
        end do

        do side = We, Ea
          if (domain_is_boundary_nested(side,child_domain)) &
            n_send_buffers = n_send_buffers + domain_nyims(child_domain)*domain_nzims(child_domain)
        end do 
        do side = So, No
          if (domain_is_boundary_nested(side,child_domain)) &
            n_send_buffers = n_send_buffers + domain_nxims(child_domain)*domain_nzims(child_domain)
        end do 
        do side = Bo, To
          if (domain_is_boundary_nested(side,child_domain)) &
            n_send_buffers = n_send_buffers + domain_nxims(child_domain)*domain_nyims(child_domain)
        end do 
      end do

      allocate(domain_bc_send_buffers(n_send_buffers))

      prev = 0
      do i_child_domain = 1, size(child_domains)

        child_domain = child_domains(i_child_domain)

        if (domain_is_boundary_nested(We,child_domain)) then   
          do dkk = 1, domain_nzims(child_domain)
            do djj = 1, domain_nyims(child_domain)
              call create_boundary_parent_buffer(num=djj+(dkk-1)*domain_nyims(child_domain) + prev, &
                                                 dir=We, domain=child_domain, im=[1,djj,dkk])
           end do
          end do

          prev = prev + domain_nyims(child_domain)*domain_nzims(child_domain)
        end if

        if (domain_is_boundary_nested(Ea,child_domain)) then   
          do dkk = 1, domain_nzims(child_domain)
            do djj = 1, domain_nyims(child_domain)
              call create_boundary_parent_buffer(num=djj+(dkk-1)*domain_nyims(child_domain) + prev, &
                                                 dir=Ea, domain=child_domain, im=[1,djj,dkk])
           end do
          end do

          prev = prev + domain_nyims(child_domain)*domain_nzims(child_domain)
        end if

        if (domain_is_boundary_nested(So,child_domain)) then   
          do dkk = 1, domain_nzims(child_domain)
            do dii = 1, domain_nxims(child_domain)
              call create_boundary_parent_buffer(num=dii+(dkk-1)*domain_nxims(child_domain) + prev, &
                                                 dir=So, domain=child_domain, im=[dii,1,dkk])
           end do
          end do
        end if

        if (domain_is_boundary_nested(No,child_domain)) then   
          prev = prev + domain_nxims(child_domain)*domain_nzims(child_domain)

          do dkk = 1, domain_nzims(child_domain)
            do dii = 1, domain_nxims(child_domain)
              call create_boundary_parent_buffer(num=dii+(dkk-1)*domain_nxims(child_domain) + prev, &
                                                 dir=No, domain=child_domain, im=[dii,domain_nyims(child_domain),dkk])
           end do
          end do

          prev = prev + domain_nxims(child_domain)*domain_nzims(child_domain)
        end if

       if (domain_is_boundary_nested(Bo,child_domain)) then   
           do djj = 1, domain_nyims(child_domain)
            do dii = 1, domain_nxims(child_domain)
              call create_boundary_parent_buffer(num=dii+(djj-1)*domain_nxims(child_domain) + prev, &
                                                 dir=Bo, domain=child_domain, im=[dii,djj,1])

           end do
          end do

          prev = prev + domain_nxims(child_domain)*domain_nyims(child_domain)
        end if

        if (domain_is_boundary_nested(To,child_domain)) then   
          do djj = 1, domain_nyims(child_domain)
            do dii = 1, domain_nxims(child_domain)
              call create_boundary_parent_buffer(num=dii+(djj-1)*domain_nxims(child_domain) + prev, &
                                                 dir=To, domain=child_domain, im=[dii,djj,domain_nzims(child_domain)])

           end do
          end do

          prev = prev + domain_nxims(child_domain)*domain_nyims(child_domain)
        end if

      end do


      if (parent_domain>0) then
        allocate(domain_bc_recv_buffers(We:To))

        allocate(domain_child_buffer)

        call create_domain_child_buffer(domain_child_buffer)

        do side = We, To
          if (is_domain_boundary_nested(side).and.is_boundary_domain_boundary(side)) then
            Btype(side) = BC_DOMAIN_COPY
            TempBtype(side) = BC_DOMAIN_COPY
            MoistBtype(side) = BC_DOMAIN_COPY
            ScalBtype(side) = BC_DOMAIN_COPY

            call create_boundary_child_buffer(dir=side, domain=parent_domain)
          endif
       end do

      end if

      call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, err)

      if (parent_domain>0) then
        if (domain_child_buffer%exchange_pr_gradient_x) &
          enable_pr_gradient_x_uniform = .true.

        if (domain_child_buffer%exchange_pr_gradient_y) &
          enable_pr_gradient_y_uniform = .true.
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

    subroutine create_boundary_parent_buffer(num, dir, domain, im)
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
            Uj2 = pos+3

            Vj1 = pos-2
            Vj2 = pos+2

            Wj1 = pos-2
            Wj2 = pos+3

            Prj1 = pos-2
            Prj2 = pos+3
          else if (dir==No) then
            Uj1 = pos-2
            Uj2 = pos+3

            Vj1 = pos-2
            Vj2 = pos+2

            Wj1 = pos-2
            Wj2 = pos+3

            Prj1 = pos-2
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
            Uk2 = pos+3

            Vk1 = pos-2
            Vk2 = pos+3

            Wk1 = pos-2
            Wk2 = pos+2

            Prk1 = pos-2
            Prk2 = pos+3
          else if (dir==To) then
            Uk1 = pos-2
            Uk2 = pos+3

            Vk1 = pos-2
            Vk2 = pos+3

            Wk1 = pos-2
            Wk2 = pos+2

            Prk1 = pos-2
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

    end subroutine create_boundary_parent_buffer


    subroutine create_boundary_child_buffer(dir, domain)
      integer, intent(in) :: dir, domain
      integer :: i, width

      associate(b => domain_bc_recv_buffers(dir))

        b%comm = world_comm
        b%remote_rank = domain_ranks_grid(domain)%arr(parent_image(1), &
                                                      parent_image(2), &
                                                      parent_image(3))       
        b%remote_domain = domain
        b%direction = dir
        b%enabled = .true.

        b%relaxation = .true.

        b%interp_order = 2

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
            Uj2 = width
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = width
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = width
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = width
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = -2
            b%r_Uj2 = +3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = +2
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = +3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = -2
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = +3
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (No)
            b%position = Prny

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = Uny-width+1
            Uj2 = Prny+3
            Uk1 = -2
            Uk2 = Unz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = Vny-width+1
            Vj2 = Prny+2
            Vk1 = -2
            Vk2 = Vnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = Wny-width+1
            Wj2 = Prny+3
            Wk1 = -2
            Wk2 = Wnz+3

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = Prny-width+1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = Prny/b%spatial_ratio-2
            b%r_Uj2 = Prny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = Unz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = Prny/b%spatial_ratio-2
            b%r_Vj2 = Prny/b%spatial_ratio+2
            b%r_Vk1 = -2
            b%r_Vk2 = Vnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = Prny/b%spatial_ratio-2
            b%r_Wj2 = Prny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = Wnz/b%spatial_ratio+3
            
            b%r_Pri1 = -1
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = Prny/b%spatial_ratio-2
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = Prnz/b%spatial_ratio+2

          case (Bo)
            b%position = 0

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = -2
            Uk2 = width

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = -2
            Vk2 = width

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = -2
            Wk2 = width

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = -1
            Prk2 = width

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = -2
            b%r_Uk2 = +3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = -2
            b%r_Vk2 = +3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = -2
            b%r_Wk2 = +2
            
            b%r_Pri1 = -2
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = -1
            b%r_Prk2 = +3

          case (To)
            b%position = Prnz

            Ui1 = -2
            Ui2 = Unx+3
            Uj1 = -2
            Uj2 = Uny+3
            Uk1 = Unz-width+1
            Uk2 = Prnz+3

            Vi1 = -2
            Vi2 = Vnx+3
            Vj1 = -2
            Vj2 = Vny+3
            Vk1 = Vnz-width+1
            Vk2 = Prnz+3

            Wi1 = -2
            Wi2 = Wnx+3
            Wj1 = -2
            Wj2 = Wny+3
            Wk1 = Wnz-width+1
            Wk2 = Prnz+2

            Pri1 = -1
            Pri2 = Prnx+2
            Prj1 = -1
            Prj2 = Prny+2
            Prk1 = Prnz-width+1
            Prk2 = Prnz+2

            b%r_Ui1 = -2
            b%r_Ui2 = Unx/b%spatial_ratio+3
            b%r_Uj1 = -2
            b%r_Uj2 = Uny/b%spatial_ratio+3
            b%r_Uk1 = Prnz/b%spatial_ratio-2
            b%r_Uk2 = Prnz/b%spatial_ratio+3
            
            b%r_Vi1 = -2
            b%r_Vi2 = Vnx/b%spatial_ratio+3
            b%r_Vj1 = -2
            b%r_Vj2 = Vny/b%spatial_ratio+3
            b%r_Vk1 = Prnz/b%spatial_ratio-2
            b%r_Vk2 = Prnz/b%spatial_ratio+3
            
            b%r_Wi1 = -2
            b%r_Wi2 = Wnx/b%spatial_ratio+3
            b%r_Wj1 = -2
            b%r_Wj2 = Wny/b%spatial_ratio+3
            b%r_Wk1 = Prnz/b%spatial_ratio-2
            b%r_Wk2 = Prnz/b%spatial_ratio+2
            
            b%r_Pri1 = -1
            b%r_Pri2 = Prnx/b%spatial_ratio+2
            b%r_Prj1 = -1
            b%r_Prj2 = Prny/b%spatial_ratio+2
            b%r_Prk1 = Prnz/b%spatial_ratio-2
            b%r_Prk2 = Prnz/b%spatial_ratio+2

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
            b%bUk2 = Uk2 - width 
            b%bVk2 = Vk2 - width 
            b%bWk2 = Wk2 - width 
            b%bPrk2 = Prk2 - width 
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

        if (has_domain_boundary_turbulence_generator(dir)) then
          b%turb_generator_enabled = .true.
          b%relaxation = .false.

          allocate(b%turb_generator)

          select case (dir)
            case (We:Ea)
              allocate(b%U_turb(Uj1:Uj2,Uk1:Uk2))
              allocate(b%V_turb(Vj1:Vj2,Vk1:Vk2))
              allocate(b%W_turb(Wj1:Wj2,Wk1:Wk2))
              allocate(b%turb_generator%sgs_tke(1:Prny, 1:Prnz))
            case (So:No)
              allocate(b%U_turb(Ui1:Ui2,Uk1:Uk2))
              allocate(b%V_turb(Vi1:Vi2,Vk1:Vk2))
              allocate(b%W_turb(Wi1:Wi2,Wk1:Wk2))
              allocate(b%turb_generator%sgs_tke(1:Prnx, 1:Prnz))
            case (Bo:To)
              allocate(b%U_turb(Ui1:Ui2,Uj1:Uj2))
              allocate(b%V_turb(Vi1:Vi2,Vj1:Vj2))
              allocate(b%W_turb(Wi1:Wi2,Wj1:Wj2))
              allocate(b%turb_generator%sgs_tke(1:Prnx, 1:Prny))
          end select

          b%turb_generator%sgs_tke = 0.01
          b%turb_generator%L_y = dymin * b%spatial_ratio / 2
          b%turb_generator%L_z = dzmin * b%spatial_ratio / 2
          b%turb_generator%T_lag = time_stepping%dt_constant * b%time_step_ratio * 2
          call b%turb_generator%init()
        end if

        if (dir==Ea) then
          b%relax_factor = 10
        end if
      end associate

    end subroutine create_boundary_child_buffer



    subroutine create_domain_parent_buffer(b, child_domain, im)
      use Strings, only: itoa
      type(dom_parent_buffer), intent(out) :: b
      integer, intent(in) :: child_domain
      integer, intent(in) :: im(3)
      real(knd) :: cxmax, cxmin, cymax, cymin, czmax, czmin
      integer :: cxi1, cxi2, cyj1, cyj2, czk1, czk2

      b%comm = world_comm
      b%remote_image = im
      b%remote_rank = domain_ranks_grid(child_domain)%arr(im(1),im(2),im(3))

      !get the child image extent
      cxmin = domain_grids(child_domain)%xmins(im(1))
      cxmax = domain_grids(child_domain)%xmaxs(im(1))
      cymin = domain_grids(child_domain)%ymins(im(2))
      cymax = domain_grids(child_domain)%ymaxs(im(2))
      czmin = domain_grids(child_domain)%zmins(im(3))
      czmax = domain_grids(child_domain)%zmaxs(im(3))

      !find the child grid indexes
      cxi1 = x_coord(cxmin)
      cxi2 = x_coord(cxmax)
      cyj1 = y_coord(cymin)
      cyj2 = y_coord(cymax)
      czk1 = z_coord(czmin)
      czk2 = z_coord(czmax)

      !check that the grids are aligned
      if (.not.(abs(cxi1*dxmin + im_xmin - cxmin) < dxmin/100)) &
        call error_stop("The west boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned.")
      if (.not.(abs(cxi2*dxmin + im_xmin - cxmax) < dxmin/100)) &
        call error_stop("The east boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned.")
      if (.not.(abs(cyj1*dymin + im_ymin - cymin) < dymin/100)) &
        call error_stop("The south boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned.")
      if (.not.(abs(cyj2*dymin + im_ymin - cymax) < dymin/100)) &
        call error_stop("The north boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned.")
      if (.not.(abs(czk1*dzmin + im_zmin - czmin) < dzmin/100)) &
        call error_stop("The bottom boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned.")
      if (.not.(abs(czk2*dzmin + im_zmin - czmax) < dzmin/100)) &
        call error_stop("The top boundary of domains "//itoa(domain_index)// &
                        " and "//itoa(child_domain)//" is not aligned.")

      !corresponds to index 0 on the child grid
      b%i1 = cxi1
      !we need index Prnx+1 on the child grid
      b%i2 = cxi2 + 1

      b%j1 = cyj1
      b%j2 = cyj2 + 1

      b%k1 = czk1
      b%k2 = czk2 + 1




      b%exchange_pr_gradient_x = enable_fixed_flow_rate .and. flow_rate_x_fixed
      b%exchange_pr_gradient_y = enable_fixed_flow_rate .and. flow_rate_y_fixed

      call MPI_ISend(b%exchange_pr_gradient_x, 1, MPI_LOGICAL, &
                     b%remote_rank, 1, b%comm, &
                     request, err)
      requests = [requests, request]

      call MPI_ISend(b%exchange_pr_gradient_y, 1, MPI_LOGICAL, &
                     b%remote_rank, 2, b%comm, &
                     request, err)
      requests = [requests, request]
    end subroutine



    subroutine create_domain_child_buffer(b)
      type(dom_child_buffer), intent(out) :: b
      integer :: i, j, k

      b%comm = world_comm
      b%remote_image = parent_image
      b%remote_rank = domain_ranks_grid(parent_domain)%arr(parent_image(1), &
                                                           parent_image(2), &
                                                           parent_image(3))
      b%spatial_ratio = domain_spatial_ratio
      b%time_step_ratio = domain_time_step_ratio

      b%r_i1 = 0
      b%r_i2 = Prnx/b%spatial_ratio + 1
      b%r_j1 = 0
      b%r_j2 = Prny/b%spatial_ratio + 1
      b%r_k1 = 0
      b%r_k2 = Prnz/b%spatial_ratio + 1

      b%r_x0 = im_xmin
      b%r_y0 = im_ymin
      b%r_z0 = im_zmin
      b%r_dx = domain_grids(parent_domain)%dx
      b%r_dy = domain_grids(parent_domain)%dy
      b%r_dz = domain_grids(parent_domain)%dz
      
      allocate(b%r_xU(b%r_i1:b%r_i2))
      allocate(b%r_yV(b%r_j1:b%r_j2))
      allocate(b%r_zW(b%r_k1:b%r_k2))
      allocate(b%r_x(b%r_i1:b%r_i2))
      allocate(b%r_y(b%r_j1:b%r_j2))
      allocate(b%r_z(b%r_k1:b%r_k2))
      b%r_xU(:) = [(b%r_x0 + i*b%r_dx , i = b%r_i1,b%r_i2)]
      b%r_yV(:) = [(b%r_y0 + i*b%r_dy , i = b%r_j1,b%r_j2)]
      b%r_zW(:) = [(b%r_z0 + i*b%r_dz , i = b%r_k1,b%r_k2)]
      b%r_x(:) = [(b%r_x0 + (i-0.5_knd)*b%r_dx , i = b%r_i1,b%r_i2)]
      b%r_y(:) = [(b%r_y0 + (i-0.5_knd)*b%r_dy , i = b%r_j1,b%r_j2)]
      b%r_z(:) = [(b%r_z0 + (i-0.5_knd)*b%r_dz , i = b%r_k1,b%r_k2)]

      call MPI_IRecv(b%exchange_pr_gradient_x, 1, MPI_LOGICAL, &
                     b%remote_rank, 1, b%comm, &
                     request, err)
      requests = [requests, request]

      call MPI_IRecv(b%exchange_pr_gradient_y, 1, MPI_LOGICAL, &
                     b%remote_rank, 2, b%comm, &
                     request, err)
      requests = [requests, request]
    end subroutine
