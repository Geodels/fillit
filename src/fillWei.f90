
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

subroutine fillInitialise(Z, m, n)
!*****************************************************************************
! This function initialises mesh parameters.

  use queues
  implicit none

  integer :: m, n
  real(kind=8),intent(in) :: Z(m,n)

  cartOn = .True.
  call cartesian_mesh_parameters(Z, m, n, 0, 1)

  return

end subroutine fillInitialise

subroutine processPit(eps)
!*****************************************************************************
! This function implements processPit from Wei 2018

  use queues
  implicit none

  integer :: c, p, nc
  real(kind=8) :: eps, iSpill, spill

  type (node)  :: pID

  do while(depressionQueue%n >0)
    pID = depressionQueue%pop()
    c = pID%id
    spill = pID%Z
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(.not.Flag(nc))then
        iSpill = FIll(nc)
        if( iSpill > spill)then
          ! Slope cell
          Flag(nc) = .True.
          call traceQueue%push(iSpill, nc)
        else
          ! Depression cell
          Fill(nc) = spill+eps
          Flag(nc) = .True.
          call depressionQueue%push(spill+eps, nc)
        endif
      endif
    enddo
  enddo

  return

end subroutine processPit

subroutine processTraceQueue()
!*****************************************************************************
! This function implements processTraceQueue from Wei 2018

  use queues
  implicit none

  type (node)  :: ptID
  integer :: c, p, k, nc, kc, i, j, ki, kj, pi, pj
  integer :: indexThreshold
  real(kind=8) :: iSpill, spill, kSpill
  logical :: HaveSpillPathOrLowerSpillOutlet

  logical :: Mask(5,5)

  ! index threshold, default to 2
  indexThreshold = 2

  do while(traceQueue%n >0)
    ptID= traceQueue%pop()
    c = ptID%id
    spill = ptID%Z
    i = meshIDs(c,1)
    j = meshIDs(c,2)
    Mask = .False.
    lp1: do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(.not.Flag(nc))then
        pi = meshIDs(nc,1)
        pj = meshIDs(nc,2)
        iSpill = Fill(nc)
        if(iSpill > spill)then
          Flag(nc) = .True.
          call traceQueue%push(iSpill, nc)
        else
          HaveSpillPathOrLowerSpillOutlet = .False.
          lp2 : do k = 1, meshNgbhsNb(nc)
            kc = meshNgbhs(nc,k)
            kSpill = Fill(kc)
            ki = meshIDs(kc,1)
            kj = meshIDs(kc,2)
            if( Mask(ki-i+2,kj-j+2) .or. (Flag(kc) .and. kSpill < spill) )then
              Mask(pi-i+2,pj-j+2) = .True.
              HaveSpillPathOrLowerSpillOutlet = .True.
              exit lp2
            endif
          enddo lp2
          if(.not. HaveSpillPathOrLowerSpillOutlet)then
            if(p<indexThreshold)then
              call potentialQueue%push(spill, c)
            else
              call priorityQueue%PQpush(spill, c)
            endif
            exit lp1
          endif
        endif
      endif
    enddo lp1
  enddo

  do while(potentialQueue%n >0)
    ptID = potentialQueue%pop()
    c = ptID%id
    lp3: do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(.not.Flag(nc))then
        call priorityQueue%PQpush(ptID%Z, c)
        exit lp3
      endif
    enddo lp3
  enddo

  return

end subroutine processTraceQueue

subroutine fillpit(m, n, eps, Filled)
!*****************************************************************************
! This function implements priority-flood depression filling algorithm from Wei 2018
! Wei, Zhou, Fu. "Efficient Priority-Flood depression filling in raster digital
! elevation models", International Journal of Digital Earth. Jan 2018,
! DOI: 10.1080/17538947.2018.1429503

  use queues
  implicit none

  integer,intent(in) :: m, n
  real(kind=8),intent(in) :: eps
  real(kind=8),intent(out) :: Filled(m,n)

  type (node)  :: ptID
  integer :: c, p, nc
  real(kind=8) :: iSpill, spill

  ! Perform pit filling using priority flood algorithm variant from Wei 2018
  do while(priorityQueue%n >0)
    ptID = priorityQueue%PQpop()
    c = ptID%id
    spill = ptID%Z
    do p = 1, meshNgbhsNb(c)
      nc = meshNgbhs(c,p)
      if(.not.Flag(nc))then
        iSpill = Fill(nc)
        if(iSpill <= spill)then
          ! Depression cell
          Fill(nc) = spill+eps
          Flag(nc) = .True.
          call depressionQueue%push(Fill(nc), nc)
          call processPit(eps)
        else
          ! Slope cell
          Flag(nc) = .True.
          call traceQueue%push(iSpill, nc)
        endif
        call processTraceQueue()
      endif
    enddo
  enddo

  do p = 1, ntot
    Filled(meshIDs(p,1),meshIDs(p,2)) = Fill(p)
  enddo

  return

end subroutine fillpit
