!1D, 2D, 3D domain decomposition via recursive coordinate bisection,
!maximazing the area of the sub-domains
!Besides re-ordered array of indexes (dd_a), the script returns
!a 2D array of (rank, [ghost cell, rank it belongs to])
!
!Platon Karpov
module rcbmod
  implicit none
  type slimer
     integer :: rank ! impi
     integer, allocatable :: ighost(:)
  end type slimer

contains

  recursive function bisect(a_arg,nmpi,counter_arg,x,y,z,index_s,index_e, counts)
    implicit none
    integer :: nmpi, counter, x, y, z, index_s, index_e, i, new_size
    integer, intent(in) :: a_arg(:)
    integer :: counter_arg
    integer, dimension(size(a_arg)) :: a
    integer, dimension(z*y*x) :: dd_a
    integer, dimension(x,y,z) :: a_dum
    integer, dimension(:), allocatable :: bisect
    integer, dimension(:,:,:), allocatable :: d1, d2
    integer, dimension(:), allocatable :: counts
    integer, dimension(nmpi) :: counts_dum
    counter = counter_arg
    a(:) = a_arg(:)

    !check if ran out of nmpi
    if (counter < nmpi) then
       !print*, a
       !reshape 1D input array to given domain dimensions
       a_dum = reshape(a(index_s:index_e),(/x,y,z/))

       !print*, a_dum(:,:,1)

       !multiply counter by 2, since it is unique to each recursive pass
       counter = counter * 2

       !if x >= y and x >= z, then cut through x
       if (size(a_dum(:,1,1)) >= size(a_dum(1,:,1)) .and. size(a_dum(:,1,1)) >= size(a_dum(1,1,:)))  then
          allocate (d1(x/2, y, z))
          allocate (d2(x - x/2, y, z))
          d1 = a_dum(1 : size(a_dum(:,1,1))/2,:,:)
          d2 = a_dum(size(a_dum(:,1,1))/2+1 : size(a_dum(:,1,1)),:,:)
          dd_a = restruct(d1,d2,counter,nmpi, counts)

          !if y >= x and y >= z, then cut through y
       elseif (size(a_dum(1,:,1)) >= size(a_dum(:,1,1)) .and. size(a_dum(1,:,1)) >= size(a_dum(1,1,:))) then
          allocate (d1(x, y/2, z))
          allocate (d2(x, y - y/2, z))
          d1 = a_dum(:,1 : size(a_dum(1,:,1))/2,:)
          d2 = a_dum(:,size(a_dum(1,:,1))/2+1 : size(a_dum(1,:,1)),:)
          dd_a = restruct(d1,d2,counter,nmpi,counts)
          !if x < y < z, then cut through z
       else
          allocate (d1(x, y, z/2))
          allocate (d2(x, y, z-z/2))
          d1 = a_dum(:,:,1 : size(a_dum(1,1,:))/2)
          d2 = a_dum(:,:,size(a_dum(1,1,:))/2+1 : size(a_dum(1,1,:)))
          dd_a = restruct(d1,d2,counter,nmpi,counts)

       endif

       !check if we have the smallest division,
       !add it's length to displ to later cut dd_a
       if (counter > nmpi/2) then
          counts_dum(:) = (/((1), i=1,nmpi)/)
          counts_dum(1:size(counts)) = counts(:)
          !print*, counts
          new_size = size(counts)+2
          deallocate (counts)
          allocate (counts(new_size))
          !print*, 'new-size -->>>>>', new_size

          counts(:) = counts_dum(1:size(counts))

          counts(new_size-1) = size(d1)
          counts(new_size) = size(d2)

       endif

       allocate (bisect(size(dd_a)))
       bisect(:) = dd_a(:)
       return

    else
       !if cannot cut anymore, return the input 'a'
       allocate (bisect(size(a)))
       bisect(:) = a(:)
       return
    endif

  end function bisect

  recursive function restruct(d1,d2,counter_arg,nmpi, counts)
    implicit none
    integer :: nmpi, counter, x1, x2, y1, y2, z1, z2
    integer, intent(in) :: d1(:,:,:), d2(:,:,:)
    integer :: counter_arg
    integer, dimension(:), allocatable :: a, d1_dum, d2_dum, restruct, dd1, dd2
    integer, dimension(:), allocatable :: counts

    counter = counter_arg

    !get dimensions of each of the 2 cuts
    x1 = size(d1(:,1,1))
    y1 = size(d1(1,:,1))
    z1 = size(d1(1,1,:))
    x2 = size(d2(:,1,1))
    y2 = size(d2(1,:,1))
    z2 = size(d2(1,1,:))

    allocate (d1_dum(x1*y1*z1))
    allocate (d2_dum(x2*y2*z2))
    allocate (a(x1*y1*z1+x2*y2*z2))
    !print*, d1
    !reshape back to 1D form
    d1_dum = reshape(d1,(/x1*y1*z1/))
    d2_dum = reshape(d2,(/x2*y2*z2/))
    !print*, x1,y1,z1
    !join them together; now 'a' is re-ordered
    a(1 : x1*y1*z1) = d1_dum(:)
    a(x1*y1*z1+1 : x1*y1*z1+x2*y2*z2) = d2_dum(:)

    !input 'a' to cut recursively
    allocate (dd1(x1*y1*z1))
    allocate (dd2(x2*y2*z2))
    dd1 = bisect(a,nmpi,counter,x1,y1,z1,1,x1*y1*z1, counts)
    dd2 = bisect(a,nmpi,counter,x2,y2,z2,x1*y1*z1+1,x1*y1*z1+x2*y2*z2, counts)

    !if cannot cut anymore, then dd1=dd2, hence we only need to return dd1
    if (counter > nmpi/2) then
       allocate (restruct(size(dd1)))
       restruct(1:size(dd1)) = dd1(1:size(dd1))
       return

       !otherwise will combine and return dd1 + dd2
    else
       allocate (restruct(size(dd1)+size(dd2)))
       restruct(1:size(dd1)) = dd1(:)
       restruct(size(dd1)+1:size(dd1)+size(dd2)) = dd2(:)
       return
    endif

  end function restruct

  subroutine check_neighbor(ighost,num_ghosts, dd_a, a_dum, counter, x_d, y_d, z_d, &
       & side, outer, inner, nmpi, counts, displs, box, x,y,z)

    integer :: escaped, look_i, look_n, x_d, y_d, z_d, index, side
    integer :: nmpi, counter, num_ghosts, i, j, outer, inner, rank, x, y, z
    integer :: one, two
    integer, dimension(nmpi) :: counts, displs
    integer :: box(x,y,z)
    integer, dimension(x_d*y_d*z_d) :: dd_a
    integer :: a_dum(x_d,y_d,z_d)
    integer :: ighost(num_ghosts)

    do j=1, outer
       do i=1, inner

          escaped = -1

          !check if the cells are on the domain boundary (side)
          !which means the ghost cells are out of the domain (-1)
          if(side==1) then
             look_i = box(i,j,1)
             look_n = look_i-x_d*y_d
             yolo1: do one = 1, y_d
                do two = 1, x_d
                   if(look_i==a_dum(two,one,1)) then
                      escaped = 1
                      exit yolo1
                   endif
                enddo
             enddo yolo1
          elseif(side==2) then
             look_i = box(i,j,z)
             look_n = look_i+x_d*y_d
             yolo2: do one = 1, y_d
                do two = 1, x_d
                   if(look_i==a_dum(two,one,z_d)) then
                      escaped = 1
                      exit yolo2
                   endif
                enddo
             enddo yolo2
          elseif(side==3) then
             look_i = box(i,1,j)
             look_n = look_i-y_d
             yolo3: do one = 1, z_d
                do two = 1, x_d
                   if(look_i==a_dum(two,1,one)) then
                      escaped = 1
                      exit yolo3
                   endif
                enddo
             enddo yolo3
          elseif(side==4) then
             look_i = box(i,y,j)
             look_n = look_i+y_d
             yolo4: do one = 1, z_d
                do two = 1, x_d
                   if(look_i==a_dum(two,y_d,one)) then
                      escaped = 1
                      exit yolo4
                   endif
                enddo
             enddo yolo4
          elseif(side==5) then
             look_i = box(1,i,j)
             look_n = look_i-1
             yolo5: do one = 1, z_d
                do two = 1, y_d
                   if(look_i==a_dum(1,two,one)) then
                      escaped = 1
                      exit yolo5
                   endif
                enddo
             enddo yolo5
          elseif(side==6) then
             look_i = box(x,i,j)
             look_n = look_i+1
             yolo6: do one = 1, z_d
                do two = 1, y_d
                   if(look_i==a_dum(x_d,two,one)) then
                      escaped = 1
                      exit yolo6
                   endif
                enddo
             enddo yolo6
          endif
          !put ghost cell index into dictionary
          ighost(counter-1) = look_n

          !if particles escaped, assign rank=-1, hence out of boundary
          if(escaped == 1) then
             ighost(counter-1) = -1
             rank = -1
             !else, find the rank, which has the ghost cell
          else
             do index = 1, size(dd_a)
                if(look_n==dd_a(index)) then
                   do rank = 1, nmpi
                      if(index>=displs(rank) .and. index<=displs(rank)+counts(rank)-1) exit
                   enddo
                   exit
                endif
             enddo
          endif

          !put rank #, which contains the ghost cell index into dictionary
          ighost(counter) = rank

          counter = counter + 2
       enddo
    enddo
    return

  end subroutine check_neighbor

  subroutine ghost_busters(dd_a, x_d, y_d, z_d, counts, displs, nmpi, sizes, ghosts)
    implicit none
    integer :: x_d, y_d, z_d, x,y,z,i, n, num_ghosts, counter
    integer :: nmpi
    integer, dimension(x_d*y_d*z_d) :: dd_a, a
    integer :: a_dum(x_d,y_d,z_d)
    integer, dimension(nmpi) :: counts, displs
    integer, allocatable :: box(:,:,:)
    integer :: sizes(nmpi*3)
    type(slimer), intent(out),allocatable :: ghosts(:)

    !create a 3D version of the domain
    a(:) = (/((i), i=1,size(a))/)
    a_dum = reshape(a(:),(/x_d,y_d,z_d/))

    allocate(ghosts(nmpi))

    !iterate through all of the ranks
    do n = 1, nmpi
       ghosts(n)%rank = n

       !get the size of the specific rank (box)
       x = sizes(n*3-2)
       y = sizes(n*3-1)
       z = sizes(n*3)

       allocate(box(x,y,z))
       box = reshape(dd_a(displs(n):displs(n)+counts(n)-1),(/x,y,z/))

       num_ghosts = 2*(2*x*y + 2*x*z + 2*y*z)

       allocate(ghosts(n)%ighost(num_ghosts))

       counter = 2

       !get ghosts of all 6 sides
       call check_neighbor(ghosts(n)%ighost,num_ghosts, dd_a, a_dum, counter, &
            & x_d, y_d, z_d, 1, y, x, nmpi, counts, displs, box, x,y,z)

       call check_neighbor(ghosts(n)%ighost,num_ghosts, dd_a, a_dum, counter, &
            & x_d, y_d, z_d, 2, y, x, nmpi, counts, displs, box, x,y,z)

       call check_neighbor(ghosts(n)%ighost,num_ghosts, dd_a, a_dum, counter, &
            & x_d, y_d, z_d, 3, z, x, nmpi, counts, displs, box, x,y,z)

       call check_neighbor(ghosts(n)%ighost,num_ghosts, dd_a, a_dum, counter, &
            & x_d, y_d, z_d, 4, z, x, nmpi, counts, displs, box, x,y,z)

       call check_neighbor(ghosts(n)%ighost,num_ghosts, dd_a, a_dum, counter, &
            & x_d, y_d, z_d, 5, z, y, nmpi, counts, displs, box, x,y,z)

       call check_neighbor(ghosts(n)%ighost,num_ghosts, dd_a, a_dum, counter, &
            & x_d, y_d, z_d, 6, z, y, nmpi, counts, displs, box, x,y,z)

       deallocate(box)
    enddo

    return

  end subroutine ghost_busters


  subroutine dimensions(dd_a, x, y, z, counts, displs, nmpi, sizes)
    implicit none
    integer :: x,y,z,i, j, k, dim_x, dim_y, dim_z
    integer :: nmpi
    integer, dimension(z*y*x) :: dd_a
    integer :: sizes(nmpi*3)
    integer, dimension(nmpi) :: counts, displs

    do i=1, nmpi
       dim_x = 0
       dim_y = 0

       !calculate dimension 'x' of the rank
       do j=displs(i), displs(i)+counts(i)-1
          if(dd_a(j)==x*y*z) then
             dim_x = dim_x + 1
             exit
          endif

          if (dd_a(j)<dd_a(displs(i))+x) then
             dim_x = dim_x + 1
          else
             exit
          endif
          if (dd_a(j+1)-dd_a(j)/=1) exit
       enddo

       !get absolute 'z' boundary of the rank to calculate 'y'
       do k=1, z
          if(x*y*(k-1)<dd_a(displs(i)) .and. dd_a(displs(i))<=x*y*k) exit
       enddo

       !calculate dimension 'y' of the rank
       do j=displs(i), displs(i)+counts(i)-1, dim_x
          if(dd_a(j)>k*x*y) exit
          dim_y = dim_y + 1
       enddo

       !calculate dimension 'z' of the rank
       dim_z = counts(i)/(dim_x*dim_y)

       sizes(i*3-2) = dim_x
       sizes(i*3-1) = dim_y
       sizes(i*3) = dim_z
    enddo

  end subroutine dimensions
end module rcbmod
