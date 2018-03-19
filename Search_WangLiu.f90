subroutine Wang_Liu2(Z,Z_flat,depr,out_list,q_out,Zmax,Nx,Ny,NODATA)
    
    use priority_queue_mod
    implicit none
    ! ����������, ������������ � ������������
    real :: Z(Nx,Ny) ! �������� ���
    real :: Z_flat(Nx,Ny) ! ��������� ���������� (��� ���������)
    integer :: depr(Nx,Ny) ! ������������� ���������� �� ������. 1 - ���������, 0 - �� ���������.
    integer :: out_list(Nx*Ny) ! ������ ����� ������ (����������, ���������� �� ����������� ������)
    integer :: q_out ! �������� ��� ������ ����� ������
    real :: Zmax ! ������������ ������ � �������� ������
    integer :: Nx,Ny ! ����������� �����
    real :: NODATA ! �������� ���� �������
    ! ���������� ���������� ������������
    type (queue) :: q ! ������� � �����������
    type (node)  :: x ! ������������ ���� ������� � �����������
    integer, allocatable :: mask(:,:) ! ������, � ������� ����� ���������� ���������������/����������������� ������. 
    ! 0 - �� �����������, 1 - � �������, 2 - �����������
    integer :: q2, q_max ! ����� ���������� ����� ��� ���������. 
    ! � ������ ������ ��������� ����� ���������� ����� � �������, ����� ��������� (�� ���� �������� "��� ������")
    integer :: percent_complete, percent_0
    integer :: c1,c2,r1,r2
    real :: Z1,Z2
    
    !! ���������� ���� 3�3 (������� �� ����������� "������" ������ ������� �������)
    !integer :: k
    !integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    !integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    
    ! ���������� ���� 3�3 (������� �������� ������ �������, ����� ������� �� ���������)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
! ������ ������
    ! ���������� ����������� ���������� � ��������.
    q_max = Nx*Ny ! ������� ������������� ���������� ����� ��� ���������
    Z_flat(1:Nx,1:Ny) = Z(1:Nx,1:Ny) ! ������ ��� ���������� "��� ���������"
    depr(1:Nx,1:Ny) = 0 ! ������� ����� ��������� (������, ����� ������)
    allocate(mask(Nx,Ny)) ! ����� �����������
    mask(1:Nx,1:Ny) = 0 
    
    ! ����������� ��������������� ������ �����
    q2 = 0
    call initialize_set()
    
    ! ����� ��������� ���������
    q_out = 0
    percent_0 = 0
    do while (q%n >0)
        x = q%top() ! ���������� ������� �������� �� �������
        Z1 = x%priority; c1 = x%c; r1 = x%r
        do k = 1,8 ! �������� ������� ������������ ��������
            c2 = c1+kx(k); r2 = r1 + ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
            if (Z(c2,r2) == nodata) cycle ! 
            
            ! �������� �������
            if (Z_flat(c2,r2) >= Z_flat(c1,r1) .and. depr(c1,r1) == 1 .and. depr(c2,r2) == 0) then 
                depr(c2,r2) = -1
                continue
            endif
            
            ! �������� �������������� ����� ������
            if (mask(c2,r2) == 1 .and. depr(c1,r1) == 1 .and. depr(c2,r2) == -1 .and. Z_flat(c2,r2) == Z_flat(c1,r1)) then
                depr(c2,r2) = -2
                q_out = q_out + 1
                out_list(q_out) = (c1 - 1) * Ny + r1
            endif
            
            if (mask(c2,r2) /= 0) cycle ! ������������ ������, ������� ��� ������������� � �������
            
            ! ���������� � �������
            q2 = q2 + 1
            call q%enqueue(Z(c2,r2), q2, c2, r2)          
            mask(c2,r2) = 1 ! ������������� �����
            
            ! �������� ���������
            if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
                Z_flat(c2,r2) = Z_flat(c1,r1)
                depr(c2,r2) = 1
                ! ���� �����, �� ������� �� ������, �� �������� ����������, ��������� � ��� ����� ������
                if( depr(c1,r1) /= 1) then 
                    depr(c1,r1) = -2
                    q_out = q_out + 1
                    out_list(q_out) = (c1 - 1) * Ny + r1
                endif
            endif
        enddo
            
        mask(c1,r1) = 2
            
        ! ����������� �������� ����������
        percent_complete = int(real(q2)/real(q_max) * 100)
        if (percent_complete > percent_0) then
            print *, "Searching for pits:", percent_complete, "% DEM scanned"
            percent_0 = percent_complete
        endif 
        
    enddo
    if(allocated(mask)) deallocate(mask)

contains

subroutine initialize_set()
! ���� �����:
! ��������������� ��� ������ ���. �� ��������� ���������, ��� ������ �� ������� ��������� � �������������� ������ (decision = .false.)
! ������, ���� ������ ��������� �� �������, ��� � �� ����� "��� ������" � ������� �������� �� ������������� (decision = .true.)
! ������ �� ���������� "��� ������" � ������ �� ���������, � ����� ��������� �� ������������
! � ����� �������� �����, ���� ������� �������������, ������ ����������� � �������������� ������

logical :: decision

! ������������� �� ������� ��� ������ �������
do c1 = 1,Nx
    do r1 = 1,Ny
        
        decision = .false.
        
        ! ���� � ������ ��� ������, ��������� � �� ������������ ��������
        if (Z(c1,r1) == nodata) then
            mask(c1,r1) = 2
            q_max = q_max - 1
            cycle
        endif
        
        ! ���� � ������ ���� ������, ��������������� � ��� ��������������
        do k = 1,8
            c2 = c1+kx(k); r2 = r1+ky(k) ! ���� � ������
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) then ! ���� ����� �� ���������� �����, ��������� ������� ������� ������ � ������
                decision = .true.
                exit
            endif
            if (Z(c2,r2) == nodata) then ! ���� ����� ����������, �� � �� "��� ������", ��������� ������� ������� ������ � ������
                decision = .true.
                exit
            endif
        enddo
        
        ! ���������� ������ � ������
        if (decision == .true.) then
            mask(c1,r1) = 1
            q2 = q2 + 1
            continue
            call q%enqueue(Z(c1,r1), q2, c1, r1)
        endif  
        
    enddo
enddo
end subroutine initialize_set

!subroutine mark_outlet(col,row)
!    
!    integer :: col, row
!    integer, allocatable :: out_list_temp
!
!    depr(col,row) = -2
!    q_out = q_out + 1
!    if (size(out_list)< q_out) then
!        allocate(out_list_temp(2*size(out_list)))
!        out_list_temp(1:q_out-1) = out_list
!        call move_alloc(out_list_temp,out_list)
!    endif
!    out_list(q_out) = (col - 1) * Ny + r1
!
!end subroutine mark_outlet

end subroutine Wang_Liu2