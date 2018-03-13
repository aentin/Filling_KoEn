! This module implements priority-queue for Fortran
! Example taken from: https://rosettacode.org/wiki/Priority_queue#Fortran
    
module priority_queue_mod
implicit none

type node
  real      :: priority ! Elevation of current node
  integer   :: c ! column index of current node
  integer   :: r ! row index of current node
end type
 
type queue
  type(node), allocatable :: buf(:)
  integer                 :: n = 0
contains
  procedure :: top
  procedure :: enqueue
  procedure :: siftdown
end type
 
contains
 
subroutine siftdown(this, a)
  class (queue)           :: this
  integer                 :: a, parent, child
  associate (x => this%buf)
  parent = a
  do while(parent*2 <= this%n)
    child = parent*2
    if (child + 1 <= this%n) then 
      !if (x(child+1)%priority > x(child)%priority ) then
      if (x(child+1)%priority < x(child)%priority) then
        child = child +1 
      end if
    end if
    !if (x(parent)%priority < x(child)%priority) then
    if (x(parent)%priority >= x(child)%priority) then
      x([child, parent]) = x([parent, child])
      parent = child
    else
      exit
    end if  
  end do      
  end associate
end subroutine
 
function top(this) result (res)
  class(queue) :: this
  type(node)   :: res
  res = this%buf(1)
  this%buf(1) = this%buf(this%n)
  this%n = this%n - 1
  call this%siftdown(1)
end function
 
subroutine enqueue(this, priority, c, r)
  class(queue), intent(inout) :: this
  real                        :: priority
  integer                     :: c, r
  type(node)                  :: x
  type(node), allocatable     :: tmp(:)
  integer                     :: i
  x%priority = priority
  x%c = c
  x%r = r
  this%n = this%n +1  
  if (.not.allocated(this%buf)) allocate(this%buf(1))
  if (size(this%buf)<this%n) then
    allocate(tmp(2*size(this%buf)))
    tmp(1:this%n-1) = this%buf
    call move_alloc(tmp, this%buf)
  end if
  this%buf(this%n) = x
  i = this%n
  do 
    i = i / 2
    if (i==0) exit
    call this%siftdown(i)
  end do
end subroutine
end module 
