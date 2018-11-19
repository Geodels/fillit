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

module queues
! This module is inspired by the fortran implementation of the priority queue
! algorithm defined in:
!                      https://rosettacode.org/wiki/Priority_queue

  implicit none

  !*****************************************************************************
  !                             QUEUES PARAMETERS DECLARATION
  !*****************************************************************************
  ! Queue node definition: index and elevation
  type node
    integer :: id
    real(kind=8) :: Z
  end type

  ! Watershed node definition: index of connected watershed and their respective elevation
  type wnode
    integer :: id
    integer :: w1
    integer :: w2
    real(kind=8) :: Z
  end type

  ! Total priority queue node definition: index, nb and elevation
  type pnodeT
    integer :: id
    integer :: nb
    real(kind=8) :: Z
  end type

  ! Definition of plain queue (no priority)
  type queue
    type(node), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: pop
    procedure :: top
    procedure :: push
  end type

  ! Definition of priority queue (1 priority: elevation)
  type pqueue
    type(node), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: PQpop
    procedure :: PQtop
    procedure :: PQpush
    procedure :: siftdown
  end type

  ! Definition of watersheds graph
  type wgraph
    type(wnode), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: wpop
    procedure :: wpush
  end type

  ! Definition of watershed priority queue  (2 priorities: elevation - first and nb - second)
  type pqueueT
    type(pnodeT), allocatable :: buf(:)
    integer :: n = 0
  contains
    procedure :: PQpopT
    procedure :: PQtopT
    procedure :: PQpushT
    procedure :: siftdownT
  end type

  !*****************************************************************************
  !                       PIT FILLING PARAMETERS DECLARATION
  !*****************************************************************************

  !  Fortran allocated variables used to defined series of mesh variables
  ! for which dimensions are unknown initially...
  integer,dimension(:,:),allocatable :: depressionPts     ! used in Zhou's algorithm to define
                                                                                          ! spill nodes and depression pit position

  ! Cartesian regular grid parameters
  integer :: ncols          ! number of cols
  integer :: nrows        ! number of rows
  integer :: cellnb        ! visited cell number
  integer :: ntot           ! total number of nodes in the DEM
  integer :: nborder    ! total number of nodes on the border of the DEM
  integer :: wlabel       ! index of watershed
  integer :: nids
  logical :: cartOn
  integer :: gridextent(4)
  integer,dimension(:),allocatable :: uextent

  logical,dimension(:),allocatable :: Flag                              ! visited node ID flag
  integer,dimension(:,:),allocatable :: labels                          ! label to specify pit and slope nodes
  integer,dimension(:,:),allocatable :: meshIDs                   ! given node ID its i,j position
  integer,dimension(:),allocatable :: meshuIDs                  ! Unstructured grid inside mesh IDs
  integer,dimension(:),allocatable :: meshBorder               ! node ID of DEM border
  integer,dimension(:,:),allocatable :: meshNgbhs              ! neighbors IDs for a given node
  integer,dimension(:),allocatable :: meshNgbhsNb          ! number of neighbors for a given node

  integer,dimension(:),allocatable :: inIDs                           ! Unstructured grid communication node IDs inside mesh
  integer,dimension(:),allocatable :: outIDs                        ! Unstructured grid communication node IDs outside mesh

  real(kind=8),dimension(:),allocatable :: Fill                       ! FIlled elevation
  real(kind=8),dimension(:,:),allocatable :: XYcoords          ! XY coordinates

  type (wgraph) :: graphW
  type(wgraph) :: graphTile

  type (queue) :: plainqueue
  type (queue) :: traceQueue
  type (queue) :: potentialQueue
  type (queue) :: depressionQueue

  type (pqueueT) :: totalqueue
  type (pqueue) :: priorityqueue

contains

  !*****************************************************************************
  !                                     PLAIN QUEUE FUNCTIONS
  !*****************************************************************************
  function top(this) result (res)
  !*****************************************************************************
  ! This function returns the top value in a plain queue

    class(queue) :: this
    type(node)   :: res

    res = this%buf(1)

  end function top

  function pop(this) result (res)
  !*****************************************************************************
  ! This function pops first values in a plain queue

    class(queue) :: this
    type(node)   :: res

    res = this%buf(1)
    this%buf(1) = this%buf(this%n)
    this%n = this%n - 1

  end function pop

  subroutine push(this, Z, id)
  !*****************************************************************************
  ! This function pushes new values in a plain queue

    class(queue), intent(inout) :: this
    real(kind=8) :: Z
    integer  :: id
    type(node)  :: x
    type(node), allocatable  :: tmp(:)

    x%Z = Z
    x%id = id
    this%n = this%n +1

    if (.not.allocated(this%buf)) allocate(this%buf(1))
    if (size(this%buf)<this%n) then
      allocate(tmp(2*size(this%buf)))
      tmp(1:this%n-1) = this%buf
      call move_alloc(tmp, this%buf)
    end if

    this%buf(this%n) = x

  end subroutine push

  !*****************************************************************************
  !                              PRIORITY QUEUE FUNCTIONS
  !*****************************************************************************
  subroutine siftdown(this, a)
  !*****************************************************************************
  ! This function sort the queue based on increasing elevations priority

    class (pqueue)  :: this
    integer :: a, parent, child

    associate (x => this%buf)
    parent = a

    do while(parent*2 <= this%n)
      child = parent*2
      if (child + 1 <= this%n) then
        if (x(child+1)%Z < x(child)%Z ) then
          child = child +1
        end if
      end if

      if (x(parent)%Z > x(child)%Z) then
        x([child, parent]) = x([parent, child])
        parent = child
      else
        exit
      end if
    end do
    end associate

  end subroutine siftdown

  function PQtop(this) result (res)
  !*****************************************************************************
  ! This function returns the top value in a priority queue

    class(pqueue) :: this
    type(node)   :: res

    res = this%buf(1)

  end function PQtop

  function PQpop(this) result (res)
  !*****************************************************************************
  ! This function pops first values in a priority queue

    class(pqueue) :: this
    type(node)   :: res

    res = this%buf(1)
    this%buf(1) = this%buf(this%n)
    this%n = this%n - 1

    call this%siftdown(1)

  end function PQpop

  subroutine PQpush(this, Z, id)
  !*****************************************************************************
  ! This function pushes new values in a priority queue

    class(pqueue), intent(inout) :: this
    real(kind=8) :: Z
    integer  :: id

    type(node)  :: x
    type(node), allocatable  :: tmp(:)
    integer :: ii

    x%Z = Z
    x%id = id
    this%n = this%n +1

    if (.not.allocated(this%buf)) allocate(this%buf(1))
    if (size(this%buf)<this%n) then
      allocate(tmp(2*size(this%buf)))
      tmp(1:this%n-1) = this%buf
      call move_alloc(tmp, this%buf)
    end if

    this%buf(this%n) = x
    ii = this%n
    do
      ii = ii / 2
      if (ii==0) exit
      call this%siftdown(ii)
    end do

  end subroutine PQpush

  !*****************************************************************************
  !                            TOTAL PRIORITY QUEUE FUNCTIONS
  !*****************************************************************************
  subroutine siftdownT(this, a)
  !*****************************************************************************
  ! This function sort the queue based on increasing elevations first and nb in second

    class (pqueueT)  :: this
    integer :: a, parent, child

    associate (x => this%buf)
    parent = a

    do while(parent*2 <= this%n)
      child = parent*2
      if (child + 1 <= this%n) then
        if (x(child+1)%Z < x(child)%Z ) then
          child = child +1
        elseif (x(child+1)%Z == x(child)%Z ) then
          if (x(child+1)%nb < x(child)%nb ) then
            child = child +1
          end if
        end if
      end if

      if (x(parent)%Z > x(child)%Z) then
        x([child, parent]) = x([parent, child])
        parent = child
      elseif (x(parent)%Z == x(child)%Z) then
        if (x(parent)%nb > x(child)%nb) then
          x([child, parent]) = x([parent, child])
          parent = child
        else
          exit
        endif
      else
        exit
      end if
    end do
    end associate

  end subroutine siftdownT

  function PQtopT(this) result (res)
  !*****************************************************************************
  ! This function returns the top value in a total priority queue

    class(pqueueT) :: this
    type(pnodeT)   :: res

    res = this%buf(1)

  end function PQtopT

  function PQpopT(this) result (res)
  !*****************************************************************************
  ! This function pops first values in a total priority queue

    class(pqueueT) :: this
    type(pnodeT)   :: res

    res = this%buf(1)
    this%buf(1) = this%buf(this%n)
    this%n = this%n - 1
    call this%siftdownT(1)

  end function PQpopT

  subroutine PQpushT(this, Z, nb, id)
  !*****************************************************************************
  ! This function pushes new values in a total priority queue

    class(pqueueT), intent(inout) :: this
    real(kind=8) :: Z
    integer  :: id
    integer  :: nb
    type(pnodeT)  :: x
    type(pnodeT), allocatable  :: tmp(:)
    integer :: ii

    x%Z = Z
    x%nb = nb
    x%id = id
    this%n = this%n +1

    if (.not.allocated(this%buf)) allocate(this%buf(1))
    if (size(this%buf)<this%n) then
      allocate(tmp(2*size(this%buf)))
      tmp(1:this%n-1) = this%buf
      call move_alloc(tmp, this%buf)
    end if

    this%buf(this%n) = x
    ii = this%n
    do
      ii = ii / 2
      if (ii==0) exit
      call this%siftdownT(ii)
    end do

  end subroutine PQpushT

  !*****************************************************************************
  !                                    GRAPH FUNCTIONS
  !*****************************************************************************
  function wpop(this) result (res)
  !*****************************************************************************
  ! This function pops first values in the watershed graph

    class(wgraph) :: this
    type(wnode)   :: res

    res = this%buf(1)
    this%buf(1) = this%buf(this%n)
    this%n = this%n - 1

  end function wpop

  subroutine wpush(this, w1, w2, Z, id)
  !*****************************************************************************
  ! This function pushes new values in the  watershed graph

    class(wgraph), intent(inout) :: this
    real(kind=8) :: Z
    integer  :: w1, w2, id, k
    type(wnode)  :: x
    type(wnode), allocatable  :: tmp(:)
    logical :: add

    if (.not.allocated(this%buf)) allocate(this%buf(1))

    add = .True.
    lp: do k = 1, this%n
      x = this%buf(k)
      if(w1 == x%w1 .and. w2 == x%w2)then
        if(Z < x%Z)then
          x%Z = Z
          x%id = id
          this%buf(k) = x
        endif
        add = .False.
        exit lp
      endif
    enddo lp

    if(add .or. this%n == 0)then
      x%Z = Z
      x%id = id
      x%w1 = w1
      x%w2 = w2
      this%n = this%n+1
      if (size(this%buf)<this%n) then
        allocate(tmp(2*size(this%buf)))
        tmp(1:this%n-1) = this%buf
        call move_alloc(tmp, this%buf)
      end if
      this%buf(this%n) = x
    endif

  end subroutine wpush

  !*****************************************************************************
  !                                     MESHING FUNCTIONS
  !*****************************************************************************
  subroutine cartesian_mesh_parameters(Z, m, n, type, step)
  !*****************************************************************************
  ! This function defines nodes and neighbors connectivity for cartesian regular mesh

    implicit none

    integer :: m,n
    integer,intent(in) :: type
    integer,intent(in) :: step
    real(kind=8),intent(in) :: Z(m,n)

    integer :: i, j, p, ip, jp, nbij, nbord, nbngbhs, id

    integer, dimension(8) :: ik=(/1, 1, 0, -1, -1,  -1,  0, 1/)
    integer, dimension(8) :: jk=(/0, 1,  1, 1,  0,  -1, -1, -1/)

    if(step == 1)then
      if(allocated(Fill)) deallocate(Fill)
      if(allocated(Flag)) deallocate(Flag)
      if(allocated(meshIDs)) deallocate(meshIDs)
      if(allocated(meshNgbhs)) deallocate(meshNgbhs)
      if(allocated(meshBorder)) deallocate(meshBorder)
      if(allocated(meshNgbhsNb)) deallocate(meshNgbhsNb)
      ncols = m
      nrows = n
      ntot = m*n
      nborder = (m-1+n-1)*2
      allocate(Fill(ntot))
      allocate(Flag(ntot))
      allocate(meshIDs(ntot,2))
      allocate(meshNgbhs(ntot,8))
      allocate(meshNgbhsNb(ntot))
      allocate(meshBorder(nborder))
    endif

    if(type==2)then
      ! Zhou labelling IDs
      if(step == 1)then
        if(allocated(labels)) deallocate(labels)
        if(allocated(meshuIDs)) deallocate(meshuIDs)
        allocate(labels(ntot,2))
        allocate(meshuIDs(ntot))
      endif
      meshuIDs = 1
      labels(:,1) = -1000
      labels(:,2) = 0
    endif

    if(step==1)then
      nbij = 1
      nbord = 1
      meshNgbhs = 0
      do j = 1, nrows
        do i = 1, ncols
          meshIDs(nbij,1) = i
          meshIDs(nbij,2) = j
          Fill(nbij) = Z(i,j)
          nbngbhs = 0
          if(i == 1 .or. i == ncols .or. j == 1 .or. j == nrows)then
            meshBorder(nbord) = nbij
            nbord = nbord + 1
          endif
          do p = 1, 8
            ip = i + ik(p)
            jp = j + jk(p)
            if(ip>=1 .and. ip<=ncols .and. jp>=1 .and. jp<=nrows)then
              nbngbhs = nbngbhs + 1
              id = (jp-1)*ncols + ip
              meshNgbhs(nbij,nbngbhs) = id
            endif
          enddo
          meshNgbhsNb(nbij) = nbngbhs
          nbij = nbij+1
        enddo
      enddo
    else
      do p = 1, ntot
        Fill(p) = Z(meshIDs(p,1),meshIDs(p,2))
      enddo
    endif

    ! Push edges to priority queue
    cellnb = 0
    Flag = .False.
    do p = 1, nborder
      id = meshBorder(p)
      if(type == 1)then
        call totalqueue%PQpushT(Fill(id), cellnb, id)
      else
        call priorityqueue%PQpush(Fill(id), id)
        if(type == 2) labels(id,1) = 0
      endif
      Flag(id) = .True.
      cellnb = cellnb + 1
    enddo

    return

  end subroutine cartesian_mesh_parameters

  subroutine tesselate(cells_nodes, cells_edges, edges_nodes, n, m, o)
  !*****************************************************************************
  ! Compute for a specific triangulation the characteristics of each node

    implicit none

    integer :: m, n, o
    integer, intent(in) :: cells_nodes(n, 3)
    integer, intent(in) :: cells_edges(n,3)
    integer, intent(in) :: edges_nodes(o, 2)

    integer :: i, n1, n2, k, l, p, eid
    integer :: nid(2), nc(3), edge(ntot, 12)
    integer :: edgeNb(3), edges(3,2), cell_ids(ntot, 12)

    logical :: inside

    cell_ids = -1
    edge = -1
    meshNgbhs = -1
    meshNgbhsNb = 0

    ! Find all cells surrounding a given vertice
    do i = 1, n
      nc = cells_nodes(i,1:3)+1
      do p = 1, 3
        inside = .False.
        lp: do k = 1, 12
          if( cell_ids(nc(p),k) == i-1 )then
            exit lp
          elseif( cell_ids(nc(p),k) == -1 )then
            inside = .True.
            exit lp
          endif
        enddo lp
        if( inside )then
          cell_ids(nc(p),k)  = i-1
        endif
      enddo
    enddo

    ! Find all edges connected to a given vertice
    do i = 1, o
      n1 = edges_nodes(i,1)+1
      n2 = edges_nodes(i,2)+1
      inside = .False.
      lp0: do k = 1, 12
        if(edge(n1,k) == i-1)then
          exit lp0
        elseif(edge(n1,k) == -1)then
          inside = .True.
          exit lp0
        endif
      enddo lp0
      if( inside )then
        edge(n1,k)  = i-1
        meshNgbhsNb(n1) = meshNgbhsNb(n1) + 1
      endif
      inside = .False.
      lp1: do k = 1, 12
        if(edge(n2,k) == i-1)then
          exit lp1
        elseif(edge(n2,k) == -1)then
          inside = .True.
          exit lp1
        endif
      enddo lp1
      if( inside )then
        edge(n2,k)  = i-1
        meshNgbhsNb(n2) = meshNgbhsNb(n2) + 1
      endif
    enddo
    do k = 1, ntot
      ! Get triangulation node IDs
      l = 0
      do eid = 1, meshNgbhsNb(k)
        nid = edges_nodes(edge(k,eid)+1,1:2)
        if( nid(1) == k-1)then
          l = l + 1
          meshNgbhs(k,l) = nid(2)+1
        else
          l = l + 1
          meshNgbhs(k,l) = nid(1)+1
        endif
      enddo
    enddo

  end subroutine tesselate

  subroutine unstructured_mesh_parameters_init(Z, m, n, o, p, type, boundary, cells_nodes, &
                                                                            cells_edges, edges_nodes)
  !*****************************************************************************
  ! This function defines nodes and neighbors connectivity for unstructured mesh

    implicit none

    integer :: m,n,p,o
    integer,intent(in) :: type
    real(kind=8),intent(in) :: Z(m)
    integer,intent(in) :: cells_nodes(n,3)
    integer,intent(in) :: cells_edges(n,3)
    integer,intent(in) :: edges_nodes(o,2)
    integer,intent(in) :: boundary(p)

    integer :: i, id

    if(allocated(Fill)) deallocate(Fill)
    if(allocated(Flag)) deallocate(Flag)
    if(allocated(meshuIDs)) deallocate(meshuIDs)
    if(allocated(meshNgbhs)) deallocate(meshNgbhs)
    if(allocated(meshBorder)) deallocate(meshBorder)
    if(allocated(meshNgbhsNb)) deallocate(meshNgbhsNb)

    ntot = m
    nborder = p

    allocate(Fill(ntot))
    allocate(Flag(ntot))
    allocate(meshuIDs(ntot))
    allocate(meshNgbhs(ntot,12))
    allocate(meshNgbhsNb(ntot))
    allocate(meshBorder(nborder))

    meshuIDs = 1
    meshBorder = boundary
    Fill = Z
    call tesselate(cells_nodes, cells_edges, edges_nodes, n, m, o)

    if(type==2)then
      ! Zhou labelling IDs
      if(allocated(labels)) deallocate(labels)
      allocate(labels(ntot,2))
      labels(:,1) = -1000
      labels(:,2) = 0
    endif

    ! Push edges to priority queue
    cellnb = 0
    Flag = .False.
    do i = 1, nborder
      id = meshBorder(i) + 1
      if(type == 1)then
        call totalqueue%PQpushT(Fill(id), cellnb, id)
      else
        call priorityqueue%PQpush(Fill(id), id)
        if(type == 2) labels(id,1) = 0
      endif
      Flag(id) = .True.
      cellnb = cellnb + 1
    enddo

    return

  end subroutine unstructured_mesh_parameters_init

  subroutine unstructured_mesh_parameters_initfast(Z, m, p, type, boundary, ngbIDs, ngbNb, IDs)
  !*****************************************************************************
  ! This function defines nodes and neighbors connectivity for unstructured mesh

    implicit none

    integer :: m,p
    integer,intent(in) :: type
    real(kind=8),intent(in) :: Z(m)
    integer,intent(in) :: ngbIDs(m,12)
    integer,intent(in) :: ngbNb(m)
    integer,intent(in) :: boundary(p)
    integer,intent(in) :: IDs(m)

    integer :: i, id

    if(allocated(Fill)) deallocate(Fill)
    if(allocated(Flag)) deallocate(Flag)
    if(allocated(meshuIDs)) deallocate(meshuIDs)
    if(allocated(meshNgbhs)) deallocate(meshNgbhs)
    if(allocated(meshBorder)) deallocate(meshBorder)
    if(allocated(meshNgbhsNb)) deallocate(meshNgbhsNb)

    ntot = m
    nborder = p

    allocate(Fill(ntot))
    allocate(Flag(ntot))
    allocate(meshuIDs(ntot))
    allocate(meshNgbhs(ntot,12))
    allocate(meshNgbhsNb(ntot))
    allocate(meshBorder(nborder))

    Fill = Z
    meshuIDs = IDs
    meshNgbhs = ngbIDs+1
    meshNgbhsNb = ngbNb
    meshBorder = boundary

    if(type==2)then
      ! Zhou labelling IDs
      if(allocated(labels)) deallocate(labels)
      allocate(labels(ntot,2))
      labels(:,1) = -1000
      labels(:,2) = 0
    endif

    ! Push edges to priority queue
    cellnb = 0
    Flag = .False.
    do i = 1, nborder
      id = meshBorder(i) + 1
      if(type == 1)then
        call totalqueue%PQpushT(Fill(id), cellnb, id)
      else
        call priorityqueue%PQpush(Fill(id), id)
        if(type == 2) labels(id,1) = 0
      endif
      Flag(id) = .True.
      cellnb = cellnb + 1
    enddo

    return

  end subroutine unstructured_mesh_parameters_initfast

  subroutine unstructured_mesh_parameters(Z, type, m)
  !*****************************************************************************
  ! This function defines nodes and neighbors connectivity for unstructured mesh

    implicit none

    integer :: m
    integer,intent(in) :: type
    real(kind=8),intent(in) :: Z(m)

    integer :: i, id

    Fill = Z

    if(type==2)then
      ! Zhou labelling IDs
      labels(:,1) = -1000
      labels(:,2) = 0
    endif

    ! Push edges to priority queue
    cellnb = 0
    Flag = .False.
    do i = 1, nborder
      id = meshBorder(i) + 1
      if(type == 1)then
        call totalqueue%PQpushT(Fill(id), cellnb, id)
      else
        call priorityqueue%PQpush(Fill(id), id)
        if(type == 2) labels(id,1) = 0
      endif
      Flag(id) = .True.
      cellnb = cellnb + 1
    enddo

    return

  end subroutine unstructured_mesh_parameters

end module queues
