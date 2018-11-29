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


!*****************************************************************************
!*****************************************************************************
!         INITIALISATION OF UNSTRUCTURED MESH PARAMETERS
!*****************************************************************************
!*****************************************************************************


subroutine build_mesh(coords, boundary, cells_nodes, cells_edges, &
                               edges_nodes, extent, m, n, o, p)
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

end subroutine build_mesh

subroutine escape_grid(coords, boundary, ngbIDs, ngbNb, &
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

end subroutine escape_grid

subroutine escape_grid_fast(Z, tp, m)
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

end subroutine escape_grid_fast

!*****************************************************************************
!*****************************************************************************
!                   PERFORM EDGES COMBINATIONS IN PARALLEL MODE
!*****************************************************************************
!*****************************************************************************

subroutine combine_edges( elev, watershed, ins, outs, newgraph, graphnb, m, n)
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
    newgraph(i,4) = wID%id-1
    i = i+1
  enddo

  return

end subroutine combine_edges

subroutine combine_edges_fast(n, elev, watershed, newgraph, graphnb, m)
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
    newgraph(i,4) = wID%id-1
    i = i+1
  enddo

  return

end subroutine combine_edges_fast

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

subroutine processPitOnePass()
!*****************************************************************************
! This function implements processPit from Zhou 2016 - one pass variant

  use queues
  implicit none

  integer :: c, p, nc

  real(kind=8) :: iSpill,spill

  type (node)  :: ptID

  do while(depressionQueue%n >0)
    ptID = depressionQueue%pop()
    c = ptID%id
    spill = ptID%Z
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(meshuIDs(nc)>0)then
        call watershedsmeet(c,nc)
        if(.not.Flag(nc))then
          iSpill = Fill(nc)
          labels(nc,2) = labels(c,2)
          if( iSpill > spill )then
            ! Slope cell
            Flag(nc) = .True.
            call traceQueue%push(iSpill, nc)
          else
            ! Depression cell
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

subroutine processTraceQueueOnePass()
!*****************************************************************************
! This function implements processTraceQueue from Zhou 2016 - one pass variant

  use queues
  implicit none

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

subroutine fillpit(m, Filled, watershedLabel, graphN)
!*****************************************************************************
! This function implements priority-flood depression filling algorithm from Zhou 2016
! Zhou, Sun, Fu. "An efficient variant of the Priority-Flood algorithm for filling
! depressions in raster digital elevation models". Computers & Geosciences.
! Vol 90, Feb 2016, pp 87–96.

  use queues
  implicit none

  integer,intent(in) :: m
  integer,intent(out) :: watershedLabel(m)
  real(kind=8),intent(out) :: Filled(m)
  integer,intent(out) :: graphN

  logical :: labelexist
  type (node)  :: ptID
  integer :: c, p, nc
  real(kind=8) :: iSpill, spill

  ! Perform pit filling using priority flood algorithm one-pass variant from Zhou 2016
  do while(priorityQueue%n >0)
    ptID = priorityqueue%PQpop()
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
            labels(nc,2) = wlabel
            Fill(nc) = spill
            Flag(nc) = .True.
            call depressionQueue%push(Fill(nc), nc)
            call processPitOnePass()
          else
            ! Slope cell
            labels(nc,2) = labels(c,2)
            Flag(nc) = .True.
            call traceQueue%push(iSpill, nc)
          endif
         call processTraceQueueOnePass() !wlabel)
        endif
      endif
    enddo
  enddo

  do p = 1, ntot
    if(meshuIDs(p)>0 .and. uextent(p)==1)  call watershedsmeet(p,-1)
  enddo

  Filled = Fill
  watershedLabel = labels(:,2)

  graphN = graphW%n

  return

end subroutine fillpit

subroutine get_spillover_nodes(graphnb, newwgraph)
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
    newwgraph(p,4) = wID%id-1
    p = p+1
  enddo

  return

end subroutine  get_spillover_nodes

!*****************************************************************************
!*****************************************************************************
!                   PERFORM DEPRESSION FILLING ON GLOBAL GRAPH
!*****************************************************************************
!*****************************************************************************

subroutine global_graph_fill(nb, cgraph, maxnghbs, nelev, spillrank, spillnodes, spillid, order, m)
!*****************************************************************************
! This function returns filled graph based on priority flood algorithm.

  use queues
  implicit none

  integer :: m
  integer, intent(in) :: nb, maxnghbs
  real(kind=8),intent(in) :: cgraph(m,5)

  integer,intent(out) :: spillrank(nb)
  integer,intent(out) :: spillnodes(nb)
  integer,intent(out) :: spillid(nb)
  integer,intent(out) :: order(nb)
  real(kind=8),intent(out) :: nelev(nb)

  integer :: k, c, nc, n1, n2, p
  logical :: inFlag(nb)
  integer :: ngbNb(nb)
  type (node)  :: ptID

  integer :: rank(nb,maxnghbs)
  integer :: ranknode(nb)
  integer :: tmp(nb)
  integer :: ngbhArr(nb,maxnghbs)
  integer :: spillnode(nb,maxnghbs)
  real(kind=8) :: spill(nb,maxnghbs)
  real(kind=8) :: spillz(nb)

  ! Initialise graph as a mesh
  inFlag = .False.
  ngbNb = 0
  ngbhArr = -1
  spillnodes = -1
  spillrank = -1
  nelev = 1.e8
  nelev(1) = -1.e8
  spillnode = -1
  spillz = 1.e8
  ranknode = -1
  spillid = -1
  rank = -1
  tmp = -1
  do k = 1, m
    n1 = int(cgraph(k,1))+1
    n2 = int(cgraph(k,2))+1
    ranknode(n1) = int(cgraph(k,5))
    ngbNb(n1) = ngbNb(n1)+1
    ngbNb(n2) = ngbNb(n2)+1
    ngbhArr(n1,ngbNb(n1)) = n2
    ngbhArr(n2,ngbNb(n2)) = n1
    spill(n1,ngbNb(n1)) = cgraph(k,3)
    spill(n2,ngbNb(n2)) = cgraph(k,3)
    spillnode(n1,ngbNb(n1)) = int(cgraph(k,4))
    spillnode(n2,ngbNb(n2)) = int(cgraph(k,4))
    rank(n1,ngbNb(n1)) = int(cgraph(k,5))
    rank(n2,ngbNb(n2)) = int(cgraph(k,5))
  enddo

  ! Perform pit filling using priority flood algorithm from Barnes 2014
  inFlag = .False.
  call priorityqueue%PQpush(nelev(1), 1)
  p = 0
  do while(priorityqueue%n > 0)
    ptID = priorityqueue%PQpop()
    c = ptID%id
    if(.not.inFlag(c))then
      p = p+1
      tmp(p) = c
      nelev(c) = ptID%Z
      inFlag(c) = .True.
      do k = 1, ngbNb(c)
        nc = ngbhArr(c,k)
        if(nc>0)then
          if(.not.inFlag(nc))then
            call priorityqueue%PQpush(max(spill(c,k),nelev(c)), nc)
            if(spillid(nc)>=0 .and. max(spill(c,k),nelev(c))<spillz(nc))then
              if(ranknode(c)==rank(c,k))then
                spillnodes(nc) = spillnode(c,k)
                spillid(nc) = c-1
                spillrank(nc) = rank(c,k)
                spillz(nc) = max(spill(c,k),nelev(c))
              endif
            elseif(spillid(nc)<0)then
              if(c==1)then
                spillnodes(nc) = spillnode(c,k)
                spillid(nc) = c-1
                spillrank(nc) = rank(c,k)
                spillz(nc) = max(spill(c,k),nelev(c))
              elseif(ranknode(c)==rank(c,k))then
                spillnodes(nc) = spillnode(c,k)
                spillid(nc) = c-1
                spillrank(nc) = rank(c,k)
                spillz(nc) = max(spill(c,k),nelev(c))
              endif
            endif
          endif
        endif
      enddo
    endif
  enddo

  order = -1
  do k =1, nb
    if(tmp(k)>-1)then
      order(tmp(k)) = k-2
    endif
  enddo

  return

end subroutine global_graph_fill

subroutine depression_info(zi, zf, area, depIDs, totpit, pitNb, pitVol, spillPts, m)
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

end subroutine depression_info
