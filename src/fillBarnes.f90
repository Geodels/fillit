!*****************************************************************************
! Copyright 2018 Tristan Salles
!
! fillIT is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or any later version.
!
! fillIT is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with fillIT.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************

subroutine fillInitialise(Z, type, step, m, n)
!*****************************************************************************
! This function initialises structured mesh parameters.

  use queues
  implicit none

  integer :: m, n
  integer,intent(in) :: type
  integer,intent(in) :: step
  real(kind=8),intent(in) :: Z(m,n)

  cartOn = .True.
  call cartesian_mesh_parameters(Z, m, n, type, step)

  return

end subroutine fillInitialise

subroutine fillinitialise_unst_init(Z, boundary, cells_nodes, cells_edges, edges_nodes, &
                                              type, m, n, o, p)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m, n, p, o
  integer,intent(in) :: type
  real(kind=8),intent(in) :: Z(m)
  integer,intent(in) :: cells_nodes(n,3)
  integer,intent(in) :: cells_edges(n,3)
  integer,intent(in) :: edges_nodes(o,2)
  integer,intent(in) :: boundary(p)

  cartOn = .False.
  call unstructured_mesh_parameters_init(Z, m, n, o, p, type, boundary, &
                                                                      cells_nodes, cells_edges, edges_nodes)

  return

end subroutine fillinitialise_unst_init

subroutine fillinitialise_unst_fast(Z, boundary, ngbIDs, ngbNb, type, m, p)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m, p
  integer,intent(in) :: type
  real(kind=8),intent(in) :: Z(m)
  integer,intent(in) :: ngbIDs(m,12)
  integer,intent(in) :: ngbNb(m)
  integer,intent(in) :: boundary(p)
  integer :: mID(m)

  cartOn = .False.
  mID = 1
  call unstructured_mesh_parameters_initfast(Z, m, p, type, boundary, &
                                                                              ngbIDs, ngbNb, mID)

  return

end subroutine fillinitialise_unst_fast

subroutine fillinitialise_unst(Z, type, m)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m
  integer,intent(in) :: type
  real(kind=8),intent(in) :: Z(m)

  cartOn = .False.
  call unstructured_mesh_parameters(Z, type, m)

  return

end subroutine fillinitialise_unst

subroutine cellconnect(nbcell, meshcells)
!*****************************************************************************
! This function defines cell connectivity for cartesian regular mesh

  use queues
  implicit none

  integer,intent(in) :: nbcell
  integer,intent(out) :: meshcells(nbcell,4)

  integer :: k, i, j, p

  k = 0
  p = 1
  do j = 1,ncols-1
    do i = 1,nrows-1
      meshcells(p,1) = k+i-1
      meshcells(p,2) = k+i
      meshcells(p,3) = k+i+nrows
      meshcells(p,4) = k+i+nrows-1
      p = p + 1
    enddo
    k = k + nrows
  enddo

  return

end subroutine cellconnect

subroutine fillpit_struct(m, n, Filled)
!*****************************************************************************
! This function implements priority-flood depression filling algorithm from Barnes 2014
! Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and
! Watershed-Labeling Algorithm for Digital Elevation Models". Computers & Geosciences.
! Vol 62, Jan 2014, pp 117–127.

  use queues
  implicit none

  integer,intent(in) :: m, n
  real(kind=8),intent(out) :: Filled(m,n)

  integer :: c, k, nc

  type (node)  :: ptID

  ! Perform pit filling using priority flood algorithm variant from Barnes 2014
  do while(priorityqueue%n >0 .or. plainqueue%n > 0)
    if( plainqueue%n > 0 )then
      ptID = plainqueue%pop()
    else
      ptID = priorityqueue%PQpop()
    endif
    c = ptID%id
    do k = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,k)
      if(.not.Flag(nc))then
        Flag(nc) = .True.
        if(Fill(nc) <= Fill(c))then
          Fill(nc) = Fill(c)
          call plainqueue%push(Fill(nc), nc)
        else
          call priorityqueue%PQpush(Fill(nc), nc)
        endif
      endif
    enddo
  enddo

  if(cartOn)then
    do k = 1, ntot
      Filled(meshIDs(k,1),meshIDs(k,2)) = Fill(k)
    enddo
  endif

  return

end subroutine fillpit_struct

subroutine fillpit_eps_struct(m, n, eps, Filled)
!*****************************************************************************
! This function implements priority-flood depression filling with epsilon algorithm from
! Barnes 2014. Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and
! Watershed-Labeling Algorithm for Digital Elevation Models". Computers & Geosciences.
! Vol 62, Jan 2014, pp 117–127.

  use queues
  implicit none

  integer,intent(in) :: m, n
  real(kind=8),intent(in) :: eps
  real(kind=8),intent(out) :: Filled(m,n)

  integer :: k, c, nc

  type (pnodeT)  :: ptID

  ! Perform pit filling using priority flood algorithm variant from Barnes 2014
  ! Here we use only one priority total queue as the plain queue doesn't ensure
  ! a consistent pit+epsilon filling in our case... not sure why?
  do while(totalqueue%n >0)
    ptID = totalqueue%PQpopT()
    c = ptID%id
    do k = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,k)
      if(.not.Flag(nc))then
        Flag(nc) = .True.
        Fill(nc) = max(Fill(nc),Fill(c)+eps)
        call totalqueue%PQpushT(Fill(nc), cellnb, nc)
        cellnb = cellnb + 1
      endif
    enddo
  enddo

  if(cartOn)then
    do k = 1, ntot
      Filled(meshIDs(k,1),meshIDs(k,2)) = Fill(k)
    enddo
  endif

  return

end subroutine fillpit_eps_struct

subroutine fillpit_unstruct(m, Filled)
!*****************************************************************************
! This function implements priority-flood depression filling algorithm from Barnes 2014
! Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and
! Watershed-Labeling Algorithm for Digital Elevation Models". Computers & Geosciences.
! Vol 62, Jan 2014, pp 117–127.

  use queues
  implicit none

  integer,intent(in) :: m
  real(kind=8),intent(out) :: Filled(m)

  integer :: c, k, nc

  type (node)  :: ptID

  ! Perform pit filling using priority flood algorithm variant from Barnes 2014
  do while(priorityqueue%n >0 .or. plainqueue%n > 0)
    if( plainqueue%n > 0 )then
      ptID = plainqueue%pop()
    else
      ptID = priorityqueue%PQpop()
    endif
    c = ptID%id
    do k = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,k)
      if(.not.Flag(nc))then
        Flag(nc) = .True.
        if(Fill(nc) <= Fill(c))then
          Fill(nc) = Fill(c)
          call plainqueue%push(Fill(nc), nc)
        else
          call priorityqueue%PQpush(Fill(nc), nc)
        endif
      endif
    enddo
  enddo

  Filled = Fill

  return

end subroutine fillpit_unstruct

subroutine fillpit_eps_unstruct(m, eps, Filled)
!*****************************************************************************
! This function implements priority-flood depression filling with epsilon algorithm from
! Barnes 2014. Barnes, Lehman, Mulla. "Priority-Flood: An Optimal Depression-Filling and
! Watershed-Labeling Algorithm for Digital Elevation Models". Computers & Geosciences.
! Vol 62, Jan 2014, pp 117–127.

  use queues
  implicit none

  integer,intent(in) :: m
  real(kind=8),intent(in) :: eps
  real(kind=8),intent(out) :: Filled(m)

  integer :: k, c, nc

  type (pnodeT)  :: ptID

  ! Perform pit filling using priority flood algorithm variant from Barnes 2014
  ! Here we use only one priority total queue as the plain queue doesn't ensure
  ! a consistent pit+epsilon filling in our case... not sure why?
  do while(totalqueue%n >0)
    ptID = totalqueue%PQpopT()
    c = ptID%id
    do k = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,k)
      if(.not.Flag(nc))then
        Flag(nc) = .True.
        Fill(nc) = max(Fill(nc),Fill(c)+eps)
        call totalqueue%PQpushT(Fill(nc), cellnb, nc)
        cellnb = cellnb + 1
      endif
    enddo
  enddo

  Filled = Fill

  return

end subroutine fillpit_eps_unstruct
