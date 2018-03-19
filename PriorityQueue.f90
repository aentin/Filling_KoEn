! This module implements priority-queue for Fortran
! Example taken from: https://rosettacode.org/wiki/Priority_queue#Fortran
    
module priority_queue_mod
implicit none

type node
  real      :: priority ! Elevation of current node
  integer   :: priority2 ! Global addition order
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
  type(node)              :: temp
  integer                 :: a, parent, child, j
  associate (x => this%buf)
  parent = a
  do while(parent*2 <= this%n)
    child = parent*2
    if (child + 1 <= this%n) then 
      !if (x(child+1)%priority > x(child)%priority ) then
      if (x(child+1)%priority < x(child)%priority) then
          child = child +1
      else if (x(child+1)%priority == x(child)%priority) then
          if (x(child+1)%priority2 < x(child)%priority2) then
              child = child +1
          endif
      end if
    end if
    !if (x(parent)%priority < x(child)%priority) then
    if (x(parent)%priority > x(child)%priority) then
      x([child, parent]) = x([parent, child])
      parent = child
    else if (x(parent)%priority == x(child)%priority) then
        if (x(parent)%priority2 > x(child)%priority2) then
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
end subroutine
 
function top(this) result (res)
  class(queue) :: this
  type(node)   :: res
  res = this%buf(1)
  this%buf(1) = this%buf(this%n)
  this%n = this%n - 1
  call this%siftdown(1)
end function
 
subroutine enqueue(this, priority, priority2, c, r)
  class(queue), intent(inout) :: this
  real                        :: priority
  integer                     :: priority2, c, r
  type(node)                  :: x
  type(node), allocatable     :: tmp(:)
  integer                     :: i
  x%priority = priority
  x%priority2 = priority2
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
!
!subroutine heapsort(a)
! 
!   real, intent(in out) :: a(0:)
!   integer :: start, n, bottom
!   real :: temp
! 
!   n = size(a)
!   do start = (n - 2) / 2, 0, -1
!     call siftdown(a, start, n);
!   end do
! 
!   do bottom = n - 1, 1, -1
!     temp = a(0)
!     a(0) = a(bottom)
!     a(bottom) = temp;
!     call siftdown(a, 0, bottom)
!   end do
!end subroutine heapsort

