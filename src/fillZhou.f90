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

subroutine fillInitialise(Z, extent, step, m, n)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m, n
  integer,intent(in) :: extent(4)
  integer,intent(in) :: step
  real(kind=8),intent(in) :: Z(m,n)

  cartOn = .True.
  wlabel = 0
  gridextent = extent
  call cartesian_mesh_parameters(Z, m, n, 2, step)

  return

end subroutine fillInitialise

subroutine fillinitialise_unst_init(coords, boundary, cells_nodes, cells_edges, edges_nodes, &
                                              extent, m, n, o, p)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m, n, p, o
  integer,intent(in) :: extent(m)
  real(kind=8),intent(in) :: coords(m,3)
  integer,intent(in) :: cells_nodes(n,3)
  integer,intent(in) :: cells_edges(n,3)
  integer,intent(in) :: edges_nodes(o,2)
  integer,intent(in) :: boundary(p)

  cartOn = .False.
  wlabel = 0

  if(allocated(XYcoords)) deallocate(XYcoords)
  allocate(XYcoords(m,2))
  if(allocated(uextent)) deallocate(uextent)
  allocate(uextent(m))
  uextent = extent
  XYcoords = coords(:,1:2)
  call unstructured_mesh_parameters_init(coords(:,3), m, n, o, p, 2, boundary, &
                                            cells_nodes, cells_edges, edges_nodes)

  return

end subroutine fillinitialise_unst_init

subroutine fillinitialise_unst_fast(coords, boundary, ngbIDs, ngbNb, &
                                                        mIDs, extent, m, p)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m, p
  integer,intent(in) :: extent(m)
  real(kind=8),intent(in) :: coords(m,3)
  integer,intent(in) :: ngbIDs(m,12)
  integer,intent(in) :: ngbNb(m)
  integer,intent(in) :: boundary(p)
  integer,intent(in) :: mIDs(m)

  cartOn = .False.
  wlabel = 0

  if(allocated(XYcoords)) deallocate(XYcoords)
  allocate(XYcoords(m,2))
  if(allocated(uextent)) deallocate(uextent)
  allocate(uextent(m))
  uextent = extent
  XYcoords = coords(:,1:2)
  call unstructured_mesh_parameters_initfast(coords(:,3), m, p, 2, boundary, &
                                                            ngbIDs, ngbNb, mIDs)

  return

end subroutine fillinitialise_unst_fast

subroutine fillinitialise_unst(Z, tp, m)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m
  integer,intent(in) :: tp
  real(kind=8),intent(in) :: Z(m)

  cartOn = .False.
  if(tp == 2)then
    wlabel = 0
  elseif(tp == 3)then
    if(.not. allocated(flowdir)) allocate(flowdir(m))
  endif
  call unstructured_mesh_parameters(Z, tp, m)

  return

end subroutine fillinitialise_unst

subroutine combineregtiles( elev, watershed, ext, newgraph, graphnb, m, n)
!*****************************************************************************
! Combine tiles along each edges based on watershed numbers and elevations

  use queues
  implicit none

  integer :: m, n
  real( kind=8 ), intent(in) :: elev(m,n)
  integer, intent(in) :: watershed(m,n)
  integer, intent(in) :: ext(4)

  integer,intent(out) :: graphnb
  real(kind=8),intent(out) :: newgraph((m+n)*2,3)

  integer :: i, ii, p, lb1, lb2

  type(wnode)  :: wID
  real(kind=8) :: eo

  integer, dimension(3) :: ik=(/-1, 0, 1/)

  ! Local edges
  do i = ext(3), ext(4)
    do p = 1, 3
      ii = i + ik(p)
      if( ii > 0 .and. ii <= m)then
        ! East
        if(watershed(i,n-1) .ne. watershed(ii,n))then
          eo = max(elev(i,n-1),elev(ii,n))
          if(watershed(i,n-1)<watershed(ii,n))then
            lb1 = watershed(i,n-1)
            lb2 = watershed(ii,n)
          else
            lb2 = watershed(i,n-1)
            lb1 = watershed(ii,n)
          endif
          call graphTile%wpush(lb1, lb2, eo, 0)
        endif
        ! West
        if(watershed(i,2) .ne.  watershed(ii,1))then
          eo = max(elev(i,2),elev(ii,1))
          if(watershed(i,2)<watershed(ii,1))then
            lb1 = watershed(i,2)
            lb2 = watershed(ii,1)
          else
            lb2 = watershed(i,2)
            lb1 = watershed(ii,1)
          endif
          call graphTile%wpush(lb1, lb2, eo, 0)
        endif
      endif
    enddo
  enddo

  do i = ext(1), ext(2)
    do p = 1, 3
      ii = i + ik(p)
      if( ii > 0 .and. ii <= n)then
        ! North
        if(watershed(m-1,i) .ne.  watershed(m,ii))then
          eo = max(elev(m-1,i),elev(m,ii))
          if(watershed(m-1,i)<watershed(m,ii))then
            lb1 = watershed(m-1,i)
            lb2 = watershed(m,ii)
          else
            lb2 = watershed(m-1,i)
            lb1 = watershed(m,ii)
          endif
          call graphTile%wpush(lb1, lb2, eo, 0)
        endif
        ! South
        if(watershed(2,i) .ne.  watershed(1,ii))then
          eo = max(elev(2,i),elev(1,ii))
          if(watershed(2,i)<watershed(1,ii))then
            lb1 = watershed(2,i)
            lb2 = watershed(1,ii)
          else
            lb2 = watershed(2,i)
            lb1 = watershed(1,ii)
          endif
          call graphTile%wpush(lb1, lb2, eo, 0)
        endif
      endif
    enddo
  enddo

  i = 1
  graphnb = graphTile%n
  do while(graphTile%n >0)
    wID = graphTile%wpop()
    newgraph(i,1) = wID%w1
    newgraph(i,2) = wID%w2
    newgraph(i,3) = wID%Z
    i = i+1
  enddo

  return

end subroutine combineregtiles

subroutine combine_unstgrids( elev, watershed, ins, outs, newgraph, graphnb, m, n)
!*****************************************************************************
! Combine unstructured grids along each edges based on watershed numbers and elevations

  use queues
  implicit none

  integer :: m, n
  real( kind=8 ), intent(in) :: elev(m)
  integer, intent(in) :: watershed(m)
  integer, intent(in) :: ins(n)
  integer, intent(in) :: outs(m)

  integer,intent(out) :: graphnb
  real(kind=8),intent(out) :: newgraph((m+n)*2,4)

  integer :: i, c, p, nc, lb1, lb2

  type(wnode)  :: wID
  real(kind=8) :: eo

  nids = n
  if(allocated(inIDs)) deallocate(inIDs)
  if(allocated(outIDs)) deallocate(outIDs)
  allocate(inIDs(nids))
  allocate(outIDs(ntot))
  inIDs = ins
  outIDs = outs

  ! Local edges
  do i = 1, nids
    c = inIDs(i)+1
    do p = 1, meshNgbhsNb(c)
        nc = meshNgbhs(c,p)
        if(outIDs(nc) > 0)then
          if(watershed(c) .ne. watershed(nc))then
            eo = max(elev(c),elev(nc))
            if(watershed(c)<watershed(nc))then
              lb1 = watershed(c)
              lb2 = watershed(nc)
            else
              lb2 = watershed(nc)
              lb1 = watershed(c)
            endif
            if(elev(c)>elev(nc))then
              call graphTile%wpush(lb1, lb2, eo, c)
            else
              call graphTile%wpush(lb1, lb2, eo, nc)
            endif
          endif
        endif
    enddo
  enddo

  i = 1
  graphnb = graphTile%n
  do while(graphTile%n >0)
    wID = graphTile%wpop()
    newgraph(i,1) = wID%w1
    newgraph(i,2) = wID%w2
    newgraph(i,3) = wID%Z
    newgraph(i,4) = wID%id
    i = i+1
  enddo

  return

end subroutine combine_unstgrids

subroutine combine_unstgrids_fast( n, elev, watershed, newgraph, graphnb, m)
!*****************************************************************************
! Combine unstructured grids along each edges based on watershed numbers and elevations

  use queues
  implicit none

  integer :: m
  integer, intent(in)  :: n
  real( kind=8 ), intent(in) :: elev(m)
  integer, intent(in) :: watershed(m)

  integer,intent(out) :: graphnb
  real(kind=8),intent(out) :: newgraph((m+n)*2,4)

  integer :: i, c, p, nc, lb1, lb2

  type(wnode)  :: wID
  real(kind=8) :: eo

  ! Local edges
  do i = 1, nids
    c = inIDs(i)+1
    do p = 1, meshNgbhsNb(c)
        nc = meshNgbhs(c,p)
        if( outIDs(nc) > 0)then
          if(watershed(c) .ne. watershed(nc))then
            eo = max(elev(c),elev(nc))
            if(watershed(c)<watershed(nc))then
              lb1 = watershed(c)
              lb2 = watershed(nc)
            else
              lb2 = watershed(nc)
              lb1 = watershed(c)
            endif
            if(elev(c)>elev(nc))then
              call graphTile%wpush(lb1, lb2, eo, c)
            else
              call graphTile%wpush(lb1, lb2, eo, nc)
            endif
          endif
        endif
    enddo
  enddo

  i = 1
  graphnb = graphTile%n
  do while(graphTile%n >0)
    wID = graphTile%wpop()
    newgraph(i,1) = wID%w1
    newgraph(i,2) = wID%w2
    newgraph(i,3) = wID%Z
    newgraph(i,4) = wID%id
    i = i+1
  enddo

  return

end subroutine combine_unstgrids_fast

subroutine graphfill(nb, cgraph, maxnghbs, nelev, m)
!*****************************************************************************
! This function returns filled graph based on priority flood algorithm.

  use queues
  implicit none

  integer :: m
  integer, intent(in) :: nb, maxnghbs
  real(kind=8),intent(in) :: cgraph(m,3)

  real(kind=8),intent(out) :: nelev(nb)

  integer :: k, c, nc, n1, n2
  logical :: inFlag(nb)
  integer :: ngbNb(nb)
  type (node)  :: ptID

  integer :: ngbhArr(nb,maxnghbs)
  real(kind=8) :: spill(nb,maxnghbs)

  ! Initialise graph as a mesh
  inFlag = .False.
  ngbNb = 0
  ngbhArr = -1
  nelev = 1.e8
  nelev(1) = -1.e8
  do k = 1, m
    n1 = int(cgraph(k,1))+1
    n2 = int(cgraph(k,2))+1
    ngbNb(n1) = ngbNb(n1)+1
    ngbNb(n2) = ngbNb(n2)+1
    ngbhArr(n1,ngbNb(n1)) = n2
    ngbhArr(n2,ngbNb(n2)) = n1
    spill(n1,ngbNb(n1)) = cgraph(k,3)
    spill(n2,ngbNb(n2)) = cgraph(k,3)
  enddo

  ! Perform pit filling using priority flood algorithm from Barnes 2014
  inFlag = .False.
  call priorityqueue%PQpush(nelev(1), 1)

  do while(priorityqueue%n > 0)
    ptID = priorityqueue%PQpop()
    c = ptID%id
    if(.not.inFlag(c))then
      nelev(c) = ptID%Z
      inFlag(c) = .True.
      do k = 1, ngbNb(c)
        nc = ngbhArr(c,k)
        if(nc>0)then
          if(.not.inFlag(nc))then
            call priorityqueue%PQpush(max(spill(c,k),nelev(c)), nc)
          endif
        endif
      enddo
    endif
  enddo

  return

end subroutine graphfill

subroutine watershedsmeet(c, nc)
!*****************************************************************************
! This function finds the connectivities between different watersheds

  use queues
  implicit none

  integer :: c, nc
  integer :: lb1, lb2
  real(kind=8) :: oelev

  if(nc == -1)then
    call graphW%wpush(labels(c,2), 0, Fill(c), c)
    return
  endif

  if(labels(nc,2) == 0) return
  if(labels(nc,2) == labels(c,2)) return

  oelev = max(Fill(c),Fill(nc))

  if(labels(nc,2)<labels(c,2))then
    lb1 = labels(nc,2)
    lb2 = labels(c,2)
  else
    lb1 = labels(c,2)
    lb2 = labels(nc,2)
  endif

  if( Fill(c)>Fill(nc) )then
    call graphW%wpush(lb1, lb2, oelev, c)
  else
    call graphW%wpush(lb1, lb2, oelev, nc)
  endif


  return

end subroutine watershedsmeet

subroutine processPitOnePass(nlabel)
!*****************************************************************************
! This function implements processPit from Zhou 2016 - one pass variant

  use queues
  implicit none

  integer :: nlabel
  integer :: c, p, nc

  real(kind=8) :: iSpill,spill

  type (node)  :: ptID
  ! type (pnodeT)  :: ptL

  do while(depressionQueue%n >0)
    ptID = depressionQueue%pop()
    c = ptID%id
    spill = ptID%Z
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(meshuIDs(nc)>0)then
        call watershedsmeet(c,nc)
        if(.not.Flag(nc))then
          iSpill = FIll(nc)
          labels(nc,2) = labels(c,2)
          if( iSpill >= spill )then
            ! Slope cell
            Flag(nc) = .True.
            labels(nc,1) = -1
            call traceQueue%push(iSpill, nc)
          else
            ! Depression cell
            labels(nc,1) = nlabel
            Fill(nc) = spill
            Flag(nc) = .True.
            call depressionQueue%push(Fill(nc), nc)
          endif
        endif
      endif
    enddo
  enddo

  return

end subroutine processPitOnePass

subroutine processTraceQueueOnePass(nlabel)
!*****************************************************************************
! This function implements processTraceQueue from Zhou 2016 - one pass variant

  use queues
  implicit none

  integer :: nlabel
  logical :: bInPQ, isBoundary
  type (node)  :: ptID
  integer :: c, p, k, nc, kc

  real(kind=8) :: iSpill, spill,kSpill

  do while(traceQueue%n >0)
    ptID = traceQueue%pop()
    c = ptID%id
    spill = ptID%Z
    bInPQ = .False.
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(meshuIDs(nc)>0)then
        call watershedsmeet(c,nc)
        if(.not.Flag(nc))then
          iSpill = Fill(nc)
          if(iSpill > spill)then
            Flag(nc) = .True.
            labels(nc,1) = -1
            call traceQueue%push(iSpill, nc)
            labels(nc,2) = labels(c,2)
          else
            if(.not.bInPQ)then
              isBoundary = .True.
              lp: do k = 1, meshNgbhsNb(nc)
                kc = meshNgbhs(nc,k)
                if(meshuIDs(kc)>0)then
                  kSpill = Fill(kc)
                  if(Flag(kc) .and. kSpill < iSpill)then
                    isBoundary = .False.
                    exit lp
                  endif
                endif
              enddo lp
              if(isBoundary)then
                labels(c,1) = -nlabel
                call priorityQueue%PQpush(spill, c)
                bInPQ = .True.
              endif
            endif
          endif
        endif
      endif
    enddo
  enddo

  return

end subroutine processTraceQueueOnePass

subroutine fillpit_struct(m, n, Filled, pitLabel, watershedLabel, graphN)
!*****************************************************************************
! This function implements priority-flood depression filling algorithm from Zhou 2016
! Zhou, Sun, Fu. "An efficient variant of the Priority-Flood algorithm for filling
! depressions in raster digital elevation models". Computers & Geosciences.
! Vol 90, Feb 2016, pp 87–96.

  use queues
  implicit none

  integer,intent(in) :: m, n
  integer,intent(out) :: pitLabel(m,n)
  integer,intent(out) :: watershedLabel(m,n)
  real(kind=8),intent(out) :: Filled(m,n)
  integer,intent(out) :: graphN

  logical :: labelexist
  type (node)  :: ptID
  integer :: nlabel, c, p, nc, i, j
  real(kind=8) :: iSpill, spill

  ! Perform pit filling using priority flood algorithm one-pass variant from Zhou 2016
  nlabel = 1
  do while(priorityQueue%n >0)
    ptID = priorityQueue%PQpop()
    c = ptID%id
    spill = ptID%Z

    if(labels(c,2) == 0)then
      labelexist = .False.
      lp: do p = 1, meshNgbhsNb(c)
        nc = meshNgbhs(c,p)
        if(labels(nc,2) > 0 .and. Fill(nc) <= Fill(c))then
          labels(c,2) = labels(nc,2)
          labelexist = .True.
          exit lp
        endif
      enddo lp
      if(.not.labelexist)then
        wlabel = wlabel + 1
        labels(c,2) = wlabel
      endif
    endif

    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      call watershedsmeet(c,nc)
      if(.not.Flag(nc))then
        iSpill = Fill(nc)
        labels(nc,2) = labels(c,2)
        if(iSpill < spill)then
          ! Depression cell
          labels(nc,1) = nlabel
          Fill(nc) = spill
          Flag(nc) = .True.
          call depressionQueue%push(Fill(nc), nc)
          call processPitOnePass(nlabel)
          nlabel = nlabel + 1
        else
          ! Slope cell
          labels(nc,1) = -1
          Flag(nc) = .True.
          call traceQueue%push(iSpill, nc)
        endif
       call processTraceQueueOnePass(nlabel)
      endif
    enddo
  enddo

  if(cartOn)then
    ! West edge
    if(gridextent(3) == 1)then
      j = 1
      do i = 1, ncols
        p = (j-1)*ncols+i
        call watershedsmeet(p,-1)
      enddo
    endif
    ! East edge
    if(gridextent(4) == 1)then
      j = nrows
      do i = 1, ncols
        p = (j-1)*ncols+i
        call watershedsmeet(p,-1)
      enddo
    endif
    ! North edge
    if(gridextent(1) == 1)then
      i = 1
      do j = 1, nrows
        p = (j-1)*ncols+i
        call watershedsmeet(p,-1)
      enddo
    endif
    ! South edge
    if(gridextent(2) == 1)then
      i = ncols
      do j = 1, nrows
        p = (j-1)*ncols+i
        call watershedsmeet(p,-1)
      enddo
    endif
  endif

  do p = 1, ntot
    Filled(meshIDs(p,1),meshIDs(p,2)) = Fill(p)
    pitLabel(meshIDs(p,1),meshIDs(p,2)) = labels(p,1)
    watershedLabel(meshIDs(p,1),meshIDs(p,2)) = labels(p,2)
  enddo

  graphN = graphW%n

  return

end subroutine fillpit_struct

subroutine fillpit_unstruct(m, Filled, pitLabel, watershedLabel, graphN)
!*****************************************************************************
! This function implements priority-flood depression filling algorithm from Zhou 2016
! Zhou, Sun, Fu. "An efficient variant of the Priority-Flood algorithm for filling
! depressions in raster digital elevation models". Computers & Geosciences.
! Vol 90, Feb 2016, pp 87–96.

  use queues
  implicit none

  integer,intent(in) :: m
  integer,intent(out) :: pitLabel(m)
  integer,intent(out) :: watershedLabel(m)
  real(kind=8),intent(out) :: Filled(m)
  integer,intent(out) :: graphN

  logical :: labelexist
  type (node)  :: ptID
  integer :: nlabel, c, p, nc
  real(kind=8) :: iSpill, spill

  ! Perform pit filling using priority flood algorithm one-pass variant from Zhou 2016
  nlabel = 1
  do while(priorityQueue%n >0)
    ptID = priorityQueue%PQpop()
    c = ptID%id
    spill = ptID%Z

    if(labels(c,2) == 0)then
      labelexist = .False.
      lp: do p = 1, meshNgbhsNb(c)
        nc = meshNgbhs(c,p)
        if(meshuIDs(nc)>0 .and. labels(nc,2) > 0 .and. Fill(nc) <= Fill(c))then
          labels(c,2) = labels(nc,2)
          labelexist = .True.
          exit lp
        endif
      enddo lp
      if(.not.labelexist)then
        wlabel = wlabel + 1
        labels(c,2) = wlabel
      endif
    endif

    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(meshuIDs(nc)>0)then
        call watershedsmeet(c,nc)
        if(.not.Flag(nc))then
          iSpill = Fill(nc)
          if(iSpill < spill)then
            wlabel = wlabel + 1
            ! Depression cell
            labels(nc,1) = nlabel
            labels(nc,2) = wlabel
            Fill(nc) = spill
            Flag(nc) = .True.
            call depressionQueue%push(Fill(nc), nc)
            call processPitOnePass(nlabel)
            nlabel = nlabel + 1
          else
            ! Slope cell
            labels(nc,2) = labels(c,2)
            labels(nc,1) = -1
            Flag(nc) = .True.
            call traceQueue%push(iSpill, nc)
          endif
         call processTraceQueueOnePass(nlabel)
        endif
      endif
    enddo
  enddo

  do p = 1, ntot
    if(meshuIDs(p)>0 .and. uextent(p)==1)  call watershedsmeet(p,-1)
  enddo

  Filled = Fill
  pitLabel = labels(:,1)
  watershedLabel = labels(:,2)

  graphN = graphW%n

  return

end subroutine fillpit_unstruct

subroutine global_flowdir(gflow, nmax, sfd, nb)
!*****************************************************************************
! This function implements flow direction and path on the edges of a given grid.

  use queues
  implicit none

  integer :: nb
  integer :: nmax
  real(kind=8),intent(in) :: gflow(nb,4)
  integer,intent(out) :: sfd(nmax,2)

  type (pqueue) :: pQ
  type (node)  :: ptID

  integer :: id, c, nc, p
  logical :: Done(nmax)
  integer :: d8(nmax)
  integer :: ngbh(nb,nb),nbNgbh(nb),idg(nmax)
  real(kind=8) :: z(nmax),border(nmax)
  ! integer :: uid(nmax)

  ! Create neighborhood
  nbNgbh = 0
  idg = -1
  do id = 1, nb
    c = int(gflow(id,1))+1
    if(idg(c) < 0) idg(c) = id-1
    nc = int(gflow(id,2))+1
    z(c) = gflow(id,3)
    border(c) = gflow(id,4)
    if(nc .ne. c)then
      nbNgbh(c) = nbNgbh(c)+1
      nbNgbh(nc) = nbNgbh(nc)+1
      ngbh(c,nbNgbh(c)) = nc
      ngbh(nc,nbNgbh(nc)) = c
    endif
  enddo

  do id = 1, nmax
    write(*,*)'id',id-1,'connect',ngbh(id,1:nbNgbh(id))-1,'border',int(border(id))
  enddo


  ! Define borders
  Done = .False.
  d8 = -1
  sfd = -1
  do id = 1, nmax
    if(border(id)>0.)then
      call pQ%PQpush(z(id), id)
      d8(id) = id
      Done(id) = .True.
      ! write(*,*)'id',id-1
    endif
  enddo

  ! Perform flow direction as a depression-carving operation from Barnes 2014
  do while(pQ%n >0)
    ptID = pQ%PQpop()
    c = ptID%id
    ! write(*,*)'entry',c-1
    do p = 1, nbNgbh(c)
      nc = ngbh(c,p)
      if(.not.Done(nc))then
        ! write(*,*)'  subject',nc-1,'point to ',c-1
        Done(nc) = .True.
        call pQ%PQpush(z(nc), nc)
        d8(nc) = c
      endif
    enddo
  enddo
  d8 = d8-1

  sfd(:,1) = d8
  sfd(:,2) = idg

  return

end subroutine global_flowdir

subroutine flowdir_unstruct(m, sfd, path, opath, bpath)
!*****************************************************************************
! This function implements flow direction and path on the edges of a given grid.

  use queues
  implicit none

  integer,intent(in) :: m
  integer,intent(out) :: sfd(m)
  integer,intent(out) :: path(m)
  integer,intent(out) :: opath(m)
  integer,intent(out) :: bpath(m)

  ! type (queue) :: Q
  type (node)  :: ptID

  integer :: c, p, nc, c0, cn, i
  real(kind=8) :: zmin
  logical :: keepLooping

  sfd = flowdir
  ! Perform flow direction as a depression-carving operation from Barnes 2014
  do while(priorityqueue%n >0)
    ptID = priorityqueue%PQpop()
    c = ptID%id
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(.not.Flag(nc) .and. outIDs(nc)==0 .and. oBorders(nc)==0 )then
        Flag(nc) = .True.
        call priorityqueue%PQpush(Fill(nc), nc)
        sfd(nc) = c
      endif
    enddo
  enddo

  ! Compute flow accumulation
  ! D = 0
  ! FA = 0
  ! do c = 1, ntot
  !   c0 = sfd(c)
  !   if(c0>0)then
  !     if(sfd(c0)>0) D(c0) = D(c0)+1
  !   endif
  ! enddo
  ! do c = 1, ntot
  !   if(D(c) == 0 .and. sfd(c)>0) call Q%push(Fill(c), c)
  ! enddo
  ! do while(Q%n > 0)
  !   ptID = Q%pop()
  !   c = ptID%id
  !   FA(c) = FA(c)+1
  !   c0 = sfd(c)
  !   if(c0>0 .and. sfd(c0)>0)then
  !     FA(c0) = FA(c0) + FA(c)
  !     D(c0) = D(c0) - 1
  !     if( D(c0) == 0)then
  !       call Q%push(Fill(c0), c0)
  !     endif
  !   endif
  ! enddo
  ! write(*,*)'ddf'
  ! Determine flow path on the perimeter based on local edges
  path = -2
  do i = 1, nborder
    c = meshBorder(i)+1
    keepLooping = .True.
    c0 = c
    if(uextent(c)==0 .or. eBorders(c)==1 )then
      ! if(nborder==16)write(*,*)'init',c-1,sfd(c)-1,nborder
      if(sfd(c) .ne. c)then
        if(uextent(sfd(c))==0 .and. mBorders(sfd(c)) == 1)then
          keepLooping = .False.
          path(c) = sfd(c)
          ! if(nborder==16)write(*,*)'   change',c-1,sfd(c)-1
        endif
      endif
    endif
    do while(keepLooping)
      ! Terminate flow as we hit the global mesh boundary  (gbounds)
      if(sfd(c) == c)then !if(uextent(c) == 1 .and. sfd(c) == c)then
         path(c0) = c    ! FlowTerminate
         keepLooping = .False.
      endif
      if(keepLooping)then
        cn = sfd(c)
        if(cn <= 0)then ! We go out the partition grid points
          if(c == c0)then
            path(c0) = -1 ! FlowExternal # Should not pass there I think...
          else
            path(c0) = c
          endif
          keepLooping = .False.
        endif
        if(outIDs(cn)>0)then
          path(c0) = c
          keepLooping = .False.
        endif
        c = cn
      endif
    enddo
  enddo
  ! write(*,*)'ddfd2'

  opath = -2
  bpath = -2
  do i = 1, nborder
    c = meshBorder(i)+1
    zmin = Fill(sfd(c))
    if(uextent(c)==0 .and. path(c)==c)then
      do p = 1, meshNgbhsNb(c)
        nc = meshNgbhs(c,p)
        if(outIDs(nc) > 0 .and. Fill(nc)<=zmin)then
          opath(c) = nc
        endif
      enddo
      ! lp: do p = 1, meshNgbhsNb(c)
      !   nc = meshNgbhs(c,p)
      !   if(uextent(nc) == 1)then
      !     bpath(c) = nc
      !     exit lp
      !   endif
      ! enddo lp
    endif
    cn = sfd(c)
  enddo

  ! do i = 1, nborder
  !   c = meshBorder(i)+1
  !   cn = sfd(c)
  !   if(outIDs(cn)>0)then
  !     opath(c) = cn
  !   endif
  !   if(uextent(c) == 0)then
  !     lp: do p = 1, meshNgbhsNb(c)
  !       nc = meshNgbhs(c,p)
  !       if(uextent(nc) == 1)then
  !         bpath(c) = nc
  !         exit lp
  !       endif
  !     enddo lp
  !   endif
  ! enddo

  sfd = sfd-1


  ! do i = 1, nids
  !   c0 = inIDs(i)+1
  !   keepLooping = .True.
  !   c = c0
  !   do while(keepLooping)
  !     ! Terminate flow as we hit the global mesh boundary  (gbounds)
  !     if(uextent(c) > 0)then
  !        path(c0) = -1
  !        keepLooping = .False.
  !     endif
  !     cn = sfd(c)+1
  !     if(uextent(cn) > 0)then
  !       path(c0) = -1
  !       keepLooping = .False.
  !     ! Local points that will be updated by the neighboring partition (idComm)
  !     elseif(outIDs(nc) > 0 )then
  !       path(c0) = c-1
  !       keepLooping = .False.
  !     elseif(cn-1 == sfd(cn)) then
  !       path(c0) = c-1
  !       keepLooping = .False.
  !     endif
  !     c = cn
  !   enddo
  ! enddo

  return

end subroutine flowdir_unstruct

subroutine flowdir_combined(sfd, nsfd, FA, m)
!*****************************************************************************
! This function implements flow direction and path on the edges of a given grid.

  use queues
  implicit none

  integer :: m
  integer,intent(in) :: sfd(m)
  integer,intent(out) :: nsfd(m)
  integer,intent(out) :: FA(m)

  type (queue) :: Q
  type (node)  :: ptID

  integer :: c, p, nc, c0, id
  integer :: D(m)


  Flag = .False.
  nsfd = -2
  do id = 1, ntot
    if(sfd(id)>-1)then
      call priorityqueue%PQpush(Fill(id), id)
      nsfd(id) = sfd(id)+1
      Flag(id) = .True.
    elseif(uextent(id)>0)then
      call priorityqueue%PQpush(Fill(id), id)
      nsfd(id) = id
      Flag(id) = .True.
    elseif(outIDs(id)>0)then
      call priorityqueue%PQpush(Fill(id), id)
      nsfd(id) = -1
      Flag(id) = .True.
    elseif(oBorders(id)>0)then
      nsfd(id) = -3
      Flag(id) = .True.
    endif
  enddo

  ! done = .False.
  ! ptID = priorityqueue%PQpop()
  ! id = ptID%id
  ! if(flowdir(id) < 0)flowdir(id) = id
  ! done(id) = .True.
  !
  ! do while(priorityqueue%n>0)
  !   ptID = priorityqueue%PQpop()
  !   id = ptID%id
  !   if(outIDs(id)==0 .and. oBorders(id)==0)then
  !     do k = 1, meshNgbhsNb(id)
  !       i = meshNgbhs(id,k)
  !       if( done(i) .and. outIDs(i)==0)then
  !         flowdir(id) = i
  !       endif
  !     enddo
  !     done(id) = .True.
  !   endif
  ! enddo
  !
  ! do id = 1, ntot
  !   if(Flag(id))then
  !     call priorityqueue%PQpush(Fill(id), id)
  !   endif
  ! enddo

  ! Perform flow direction as a depression-carving operation from Barnes 2014
  do while(priorityqueue%n >0)
    ptID = priorityqueue%PQpop()
    c = ptID%id
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(.not.Flag(nc))then
        Flag(nc) = .True.
        call priorityqueue%PQpush(Fill(nc), nc)
        nsfd(nc) = c
      endif
    enddo
  enddo

  ! Compute flow accumulation
  D = 0
  FA = 0
  do c = 1, ntot
    c0 = nsfd(c)
    if(c0>0)then
      if(nsfd(c0)>0) D(c0) = D(c0)+1
    endif
  enddo
  do c = 1, ntot
    if(D(c) == 0 .and. nsfd(c)>0) call Q%push(Fill(c), c)
  enddo
  do while(Q%n > 0)
    ptID = Q%pop()
    c = ptID%id
    FA(c) = FA(c)+1
    c0 = nsfd(c)
    if(c0>0 .and. nsfd(c0)>0)then
      FA(c0) = FA(c0) + FA(c)
      D(c0) = D(c0) - 1
      if( D(c0) == 0)then
        call Q%push(Fill(c0), c0)
      endif
    endif
  enddo
  nsfd = nsfd-1

  return

end subroutine flowdir_combined

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

subroutine spillPts(graphnb, newwgraph)
!*****************************************************************************
! This function returns pit/depression graph.

  use queues
  implicit none

  integer,intent(in) :: graphnb
  real(kind=8),intent(out) :: newwgraph(graphnb,4)
  type(wnode)  :: wID
  integer :: p

  p = 1
  do while(graphW%n >0)
    wID = graphW%wpop()
    newwgraph(p,1) = wID%w1
    newwgraph(p,2) = wID%w2
    newwgraph(p,3) = wID%Z
    newwgraph(p,4) = wID%id
    p = p+1
  enddo

  return

end subroutine  spillPts

subroutine params_unstPit(zi, zf, area, depIDs, totpit, pitNb, pitVol, spillPts, m)
!*****************************************************************************
! Extract for unstructured grid pit information: unique label and volume

use queues
implicit none

integer :: m
integer, intent(in) :: totpit
real(kind=8), intent(in) :: area(m)
real(kind=8), intent(in) :: zf(m)
real(kind=8), intent(in) :: zi(m)
integer, intent(in) :: depIDs(m)

integer, intent(out) :: pitNb(m)
real(kind=8), intent(out) :: pitVol(totpit)
real(kind=8), intent(out) :: spillPts(totpit,2)

real(kind=8) :: z0
integer :: depID(m)
integer :: pitLabel(totpit)
integer :: p, k, k1, k2, kk, k0

pitLabel = -1
depID = depIDs
spillPts = -1

do k = 1,m
  if(meshuIDs(k)>0)then
    do kk = 1, meshNgbhsNb(k)
      p = meshNgbhs(k,kk)

      if(depID(k)>0 .and. depID(p)>0)then
        if(depID(k)>depID(p))then
          if(pitLabel(depID(k))>0)then
            pitLabel(depID(k)) = min(pitLabel(depID(k)),depID(p))
          else
            pitLabel(depID(k)) = depID(p)
          endif
        elseif(depID(k)<depID(p))then
          if(pitLabel(depID(p))>0)then
            pitLabel(depID(p)) = min(pitLabel(depID(p)),depID(k))
          else
            pitLabel(depID(p)) = depID(k)
          endif
        endif
      endif

      if(depID(k)>0 .and. depID(p) .ne. depID(k))then
        if(zf(p) <= zf(k))then
          k0 = int(spillPts(depID(k),1))
          z0 = spillPts(depID(k),2)
          if(k0>0)then
            if(z0>zf(p))then
              spillPts(depID(k),1) = p
              spillPts(depID(k),2) = zf(p)
            endif
          else
            spillPts(depID(k),1) = p
            spillPts(depID(k),2) = zf(p)
          endif
        endif
      endif

    enddo
  endif
enddo

do k = 1, m
  if(depID(k)>0 .and. meshuIDs(k)>0)then
    kk = pitLabel(depID(k))
    if(kk>0) depID(k) = kk
  endif
enddo

pitLabel = -1

! Local edges
do k0 = 1, nids
  k = inIDs(k0)+1
  do p = 1, meshNgbhsNb(k)
      kk = meshNgbhs(k,p)
      if(outIDs(kk) > 0)then
        k1 = depID(k)
        k2 = depID(kk)
        if( k2>0 .and. k1>0)then
          if(k2<k1) pitLabel(k1) = k2
          if(k2>k1) pitLabel(k2) = k1
        endif
      endif
  enddo
enddo

pitNb =  -1
pitVol = 0.

do k = 1, m
  pitNb(k) = depID(k)
  if( depID(k)>0)then
    if(pitLabel(depID(k))>0)then
      pitNb(k) = pitLabel(depID(k))
    else
      pitNb(k) = depID(k)
    endif
    if(meshuIDs(k)>0) pitVol(pitNb(k)) = pitVol(pitNb(k))+(zf(k)-zi(k))*area(k)
  endif
enddo

return

end subroutine params_unstPit

subroutine params_regPit(zi, zf, area, depIDs, totpit, pitNb, pitVol, spillPts, m, n)
!*****************************************************************************
! Extract for regular grid pit information: unique label and volume

  use queues
  implicit none

  integer :: m, n
  integer, intent(in) :: totpit
  real(kind=8), intent(in) :: area
  real(kind=8), intent(in) :: zf(m,n)
  real(kind=8), intent(in) :: zi(m,n)
  integer, intent(in) :: depIDs(m,n)

  integer, intent(out) :: pitNb(m,n)
  real(kind=8), intent(out) :: pitVol(totpit)
  real(kind=8), intent(out) :: spillPts(totpit,3)

  real(kind=8) :: z0
  integer :: depID(m,n)
  integer :: pitLabel(totpit)
  integer :: p, k, p1, p2, k1, k2, r, kk, pp, p0, k0

  integer, dimension(3) :: ik=(/-1, 0, 1/)

  pitLabel = -1
  depID = depIDs
  spillPts = -1
  do k = 1,ntot
    k1 = meshIDs(k,1)
    k2 = meshIDs(k,2)
    do p = 1, meshNgbhsNb(k)
      pp = meshNgbhs(k,p)
      p1 = meshIDs(pp,1)
      p2 = meshIDs(pp,2)
      if(depID(k1,k2)>-1 .and. depID(p1,p2)>-1)then
        if(depID(k1,k2)>depID(p1,p2))then
          if(pitLabel(depID(k1,k2))>-1)then
            pitLabel(depID(k1,k2)) = min(pitLabel(depID(k1,k2)),depID(p1,p2))
          else
            pitLabel(depID(k1,k2)) = depID(p1,p2)
          endif
        else
          if(pitLabel(depID(p1,p2))>-1)then
            pitLabel(depID(p1,p2)) = min(pitLabel(depID(p1,p2)),depID(k1,k2))
          else
            pitLabel(depID(p1,p2)) = depID(k1,k2)
          endif
        endif
      endif
      if(depID(k1,k2)>-1 .and. depID(p1,p2) .ne. depID(k1,k2))then
        if(zf(p1,p2) <= zf(k1,k2))then
          k0 = int(spillPts(depID(k1,k2),1))
          p0 = int(spillPts(depID(k1,k2),2))
          z0 = spillPts(depID(k1,k2),3)
          if(k0>0 .and. p0>0)then
            if(z0>zf(p1,p2))then
              spillPts(depID(k1,k2),1) = p1
              spillPts(depID(k1,k2),2) = p2
              spillPts(depID(k1,k2),3) = zf(p1,p2)
            endif
          else
            spillPts(depID(k1,k2),1) = p1
            spillPts(depID(k1,k2),2) = p2
            spillPts(depID(k1,k2),3) = zf(p1,p2)
          endif
        endif
      endif
    enddo
  enddo

  do k = 1, m
    do p = 1, n
      if(depID(k,p)>0)then
          kk = pitLabel(depID(k,p))
          if(kk>0) depID(k,p) = kk
      endif
    enddo
  enddo

  pitLabel = -1

  do p = 2, n-1
    do r = 1, 3
      pp = p + ik(r)
      if( pp > 0 .and. pp <= n)then
        k1 = depID(2,p)
        k2 = depID(1,pp)
        if( k2>0 .and. k1>0)then
          if(k2<k1) pitLabel(k1) = k2
          if(k2>k1) pitLabel(k2) = k1
        endif
        k1 = depID(m-1,p)
        k2 = depID(m,pp)
        if( k2>0 .and. k1>0)then
          if(k2<k1) pitLabel(k1) = k2
          if(k2>k1) pitLabel(k2) = k1
        endif
      endif
    enddo
  enddo

  do k = 2, m-1
    do r = 1, 3
      kk = k + ik(r)
      if( kk > 0 .and. kk <= m)then
        p1 = depID(k,2)
        p2 = depID(kk,1)
        if( p1>0  .and. p2>0)then
          if(p1<p2) pitLabel(p2) = p1
          if(p1>p2) pitLabel(p1) = p2
        endif
        p1 = depID(k,n-1)
        p2 = depID(kk,n)
        if( p1>0  .and. p2>0)then
          if(p1<p2) pitLabel(p2) = p1
          if(p1>p2) pitLabel(p1) = p2
        endif
      endif
    enddo
  enddo

  pitNb =  -1
  pitVol = 0.

  do k = 1, m
    do p = 1, n
      pitNb(k,p) = depID(k,p)
      if( depID(k,p)>-1)then
        if(pitLabel(depID(k,p))>-1)then
          pitNb(k,p) = pitLabel(depID(k,p))
        else
          pitNb(k,p) = depID(k,p)
        endif
        if(k>1 .and. k<m .and. p>1 .and. p<n)then
          pitVol(pitNb(k,p)) = pitVol(pitNb(k,p))+(zf(k,p)-zi(k,p))*area
        endif
      endif
    enddo
  enddo

  return

end subroutine params_regPit
