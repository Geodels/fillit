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

subroutine fillinitialise_unst(Z, m)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m
  real(kind=8),intent(in) :: Z(m)

  cartOn = .False.
  wlabel = 0
  call unstructured_mesh_parameters(Z, 2, m)

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
          call graphTile%wpush(lb1, lb2, eo)
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
          call graphTile%wpush(lb1, lb2, eo)
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
          call graphTile%wpush(lb1, lb2, eo)
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
          call graphTile%wpush(lb1, lb2, eo)
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
  real(kind=8),intent(out) :: newgraph((m+n)*2,3)

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
            call graphTile%wpush(lb1, lb2, eo)
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
  real(kind=8),intent(out) :: newgraph((m+n)*2,3)

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
            call graphTile%wpush(lb1, lb2, eo)
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
        if(.not.inFlag(nc))then
          call priorityqueue%PQpush(max(spill(c,k),nelev(c)), nc)
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
    call graphW%wpush(labels(c,2), 0, Fill(c))
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

  call graphW%wpush(lb1, lb2, oelev)

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

subroutine spillPts(graphnb, wgraph)
!*****************************************************************************
! This function returns pit/depression graph.

  use queues
  implicit none

  integer,intent(in) :: graphnb
  real(kind=8),intent(out) :: wgraph(graphnb,3)
  type(wnode)  :: wID
  integer :: p

  p = 1
  do while(graphW%n >0)
    wID = graphW%wpop()
    wgraph(p,1) = wID%w1
    wgraph(p,2) = wID%w2
    wgraph(p,3) = wID%Z
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
