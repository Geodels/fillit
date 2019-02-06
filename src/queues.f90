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

  !*****************************************************************************
  !                       PIT FILLING PARAMETERS DECLARATION
  !*****************************************************************************

  !  Fortran allocated variables used to defined series of mesh variables
  ! for which dimensions are unknown initially...
  integer,dimension(:,:),allocatable :: depressionPts     ! used in Zhou's algorithm to define
                                                                                          ! spill nodes and depression pit position

  ! Cartesian regular grid parameters
  integer :: ntot           ! total number of nodes in the DEM
  integer :: gtot           ! total number of nodes in the DEM
  integer :: nborder    ! total number of nodes on the border of the DEM
  integer :: gborder    ! total number of nodes on the border of the DEM
  integer :: wlabel       ! index of watershed
  integer :: nids
  integer,dimension(:),allocatable :: uextent

  logical,dimension(:),allocatable :: Flag                              ! visited node ID flag
  integer,dimension(:,:),allocatable :: labels                          ! label to specify pit and slope nodes
  integer,dimension(:,:),allocatable :: meshIDs                   ! given node ID its i,j position
  integer,dimension(:),allocatable :: meshuIDs                  ! Unstructured grid inside mesh IDs
  integer,dimension(:),allocatable :: meshBorder               ! node ID of DEM border
  integer,dimension(:),allocatable :: mBorders               ! node ID of DEM border
  integer,dimension(:,:),allocatable :: meshNgbhs              ! neighbors IDs for a given node
  integer,dimension(:),allocatable :: meshNgbhsNb          ! number of neighbors for a given node
  integer,dimension(:),allocatable :: oBorders               ! node ID of DEM border
  integer,dimension(:),allocatable :: eBorders               ! node ID of DEM border
  integer,dimension(:),allocatable :: idpit

  integer,dimension(:,:),allocatable :: GmeshNgbhs              ! neighbors IDs for a given node
  integer,dimension(:),allocatable :: GmeshNgbhsNb          ! number of neighbors for a given node
  integer,dimension(:),allocatable :: GmeshBounds          ! number of neighbors for a given node
  real(kind=8),dimension(:),allocatable :: GmeshArea
  real(kind=8),dimension(:),allocatable :: gpitVol
  real(kind=8),dimension(:),allocatable :: zmin          

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
  subroutine mesh_parameters(Z, seaIDs, m, p, n, boundary, ngbIDs, ngbNb, IDs)
  !*****************************************************************************
  ! This function defines nodes and neighbors connectivity for unstructured mesh

    implicit none

    integer :: m, p, n
    real(kind=8),intent(in) :: Z(m)
    integer,intent(in) :: ngbIDs(m,12)
    integer,intent(in) :: ngbNb(m)
    integer,intent(in) :: boundary(p)
    integer,intent(in) :: IDs(m)
    integer,intent(in) :: seaIDs(n)

    integer :: i, id

    if(allocated(Fill)) deallocate(Fill)
    if(allocated(Flag)) deallocate(Flag)
    if(allocated(meshuIDs)) deallocate(meshuIDs)
    if(allocated(meshNgbhs)) deallocate(meshNgbhs)
    if(allocated(meshBorder)) deallocate(meshBorder)
    if(allocated(mBorders)) deallocate(mBorders)
    if(allocated(meshNgbhsNb)) deallocate(meshNgbhsNb)

    ntot = m
    nborder = p

    allocate(Fill(ntot))
    allocate(Flag(ntot))
    allocate(meshuIDs(ntot))
    allocate(meshNgbhs(ntot,12))
    allocate(meshNgbhsNb(ntot))
    allocate(meshBorder(nborder))
    allocate(mBorders(ntot))

    Fill = Z
    meshuIDs = IDs
    meshNgbhs = ngbIDs+1
    meshNgbhsNb = ngbNb
    meshBorder = boundary
    mBorders = 0

    ! Zhou labelling IDs
    if(allocated(labels)) deallocate(labels)
    allocate(labels(ntot,2))
    labels(:,1) = -1000
    labels(:,2) = 0

    ! Push edges to priority queue
    Flag = .False.
    do i = 1, nborder
      id = meshBorder(i) + 1
      mBorders(id) = 1
      call priorityqueue%PQpush(Fill(id), id)
      labels(id,1) = 0
      Flag(id) = .True.
    enddo

    ! Push deep marine vertices to priority queue
    do i = 1, n
      id = seaIDs(i) + 1
      if(id>0)then
        call priorityqueue%PQpush(Fill(id), id)
        labels(id,1) = 0
        Flag(id) = .True.
      endif
    enddo

    return

  end subroutine mesh_parameters

  subroutine mesh_parameters_fast(Z, seaIDs, m, n)
  !*****************************************************************************
  ! This function defines nodes and neighbors connectivity for unstructured mesh

    implicit none

    integer :: m, n
    real(kind=8),intent(in) :: Z(m)
    integer,intent(in) :: seaIDs(n)

    type (node)  :: ptID


    ! type (pqueue) :: ptyqueue
    ! type (node)  :: ptID
    integer :: i, id, k
    real(kind=8) :: h
    logical :: done(ntot)

    Fill = Z
    ! Zhou labelling IDs
    labels(:,1) = -1000
    labels(:,2) = 0

    ! Push edges to priority queue
    Flag = .False.
    do i = 1, nborder
      id = meshBorder(i) + 1
      call priorityqueue%PQpush(Fill(id), id)
      labels(id,1) = 0
      Flag(id) = .True.
    enddo

    ! Push deep marine vertices to priority queue
    do i = 1, n
      id = seaIDs(i) + 1
      if(id>0)then
        call priorityqueue%PQpush(Fill(id), id)
        labels(id,1) = 0
        Flag(id) = .True.
      endif
    enddo

    return

  end subroutine mesh_parameters_fast

end module queues
