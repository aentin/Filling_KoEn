! ������������ FILLING: ���������� ��������� ��������� � �������������� ��������� A* � �������� ������������ ����� ����� �������
!   ������������ �������� ��� ��������������� � ���������������� ��������, ����������� ���������� ��� ���������� 
!   ������� ���������:
!   * Z - �������� ������� �����
!   * Nx, Ny - ���������� �����/�������� � �������
!   * Zmax - ������������ �������� ������ �� �������� ������� Z
!   * StepX, StepY - ��� ����� �� �����������
!   * NODATA - �������� ���� ������� �� �������� � ��������������� �������
!   �������� ���������:
!   * Zout - �������������� (�����������) ������� �����. ����������� ��������� ���������� ����������� ������� Z
!
! ������� ������������:
!   Wang-Liu � ������������ ��� ���������� ��������� ��������� ���� ���� �� ��������� �����-�� (��������� ���������� ������������ ��� ����������� ��������� � ����� ������)
!   dist_to_borders_Z_matrix � ������������ ��� ���������� ���������� �� ����� ��������� �� �����-������
!   matrix_Dijkstra � ��������������� ������������ ���������� ��������� � �������������� ��������� ��������. �������� �� �������� ����� ��������� � ���� �����
! ��� ������ ����������� ������������ ��������� � ����� CONTAINS
!
! ������� ������:
! 1) ���������� ��������� �� ��������� �����-�� (� ���������� �������������� ������������) � ������������� ����������� �� ���������������;
! 2) ���������� ������������� ������� ��� �������� �����
! 3) � ����� ��� ������ ���� ���������� � ����� ������:
!       �) ���������� ����� ����������� ��� ���� ���������� � ����� ������;
!       �) ���������� ����� ����������� � ������������� ��������;
! 4) ����� � ����� �������� �����

    subroutine FILLING(Z,Zout,Nx,Ny,Zmin,Zmax,StepX,StepY,NODATA)
    implicit none
    
    ! ����������, ������������ � ���������
    real :: Z(Nx,Ny) ! �������� ������
    real :: Zout(Nx,Ny) ! ��������� ����������
    real(4) :: Zmin, Zmax ! ����������� � ������������ �������� ������ �� �������� ������
    integer :: Nx,Ny ! ����������� �������
    real(8) :: StepX,StepY ! ��� �����
    real :: NODATA ! �������� ���� �������
    
    ! ���������� ���������� ������������
    real, allocatable :: Zout_temp(:,:) ! ������������� ������� �����
    real, allocatable :: Z_flat(:,:) ! ��� ������ ��������� ���������� ����� ������, ����������� ���� ����������
    integer, allocatable :: depr(:,:) ! �������, ������ ������� �������� ���������� � ����������, �������� � ������ ������
    integer, allocatable :: outlets_list(:) ! ������ ����� ������
    integer :: number_of_outlets
    ! ��������� ��� ������
    integer :: c1,c2,r1,r2,i,j,q

    ! ��� ������ � ��������
    character(25) :: DATE, ZEIT, ZONE
    integer :: time_values(8)
    
    integer :: percent_complete, percent_0
    
! ������ ������
    
    ! 1. ���������� ��������� ���� ���������� (� ����� ����������� ��������������� ���������)
    allocate(Z_flat(Nx,Ny))
    allocate(depr(Nx,Ny))
    allocate(outlets_list(Nx*Ny))
    outlets_list(1:Nx*Ny) = 0
    
    CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
    write(*,'(A,i2,A,i2,a,i2)') "Wang-Liu started  at:    ", time_values(5), ":", time_values(6), ":", time_values(7)
    
    call Wang_Liu2(Z,Z_flat,depr,outlets_list,number_of_outlets,Zmax,Nx,Ny,NODATA)
    
    CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
    write(*,'(A,i2,A,i2,a,i2)') "Wang-Liu finished at:    ", time_values(5), ":", time_values(6), ":", time_values(7)
    
    print *, "Found", number_of_outlets, "outlet points"
    
    ! �� ������������ ��������: ��������� ������, ������� ��������� (�������������), ������ �����-�������
    
    ! 2. ����� ������
    !call find_borders()
    ! � �������������� �������: -1 - �������, 0 - �� �������
    
    ! 3. ���������� �������������� ������� (������ ����� ��������� ���������� �� �������� ���� �������
    do i = 1,Nx
        do j = 1,Ny
            if (depr(i,j) == 0) then
                Zout(i,j) = Z(i,j)
            else
                Zout(i,j) = NODATA
            endif
        enddo
    enddo
    
    ! 4. ����������, ���������� ���������. ���������������, �� ����� ���� ���������� � ����� ������
    CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
    write(*,'(A,i2,A,i2,a,i2)') "Filling started at:    ", time_values(5), ":", time_values(6), ":", time_values(7)
    
    !!!! ���������� �� "�������" ���������������� ������ ��������� ���������
    !Zout(1:Nx,1:Ny) = depr(1:Nx,1:Ny)
    !return
    
    percent_0 = 0
    allocate(Zout_temp(Nx,Ny)) ! ������������� ������� ��� ������ ���������� ����������
    do q = 1, number_of_outlets
        if (outlets_list(q) == 0) exit ! �������, ����� ��������� �� ������ ����� ������
        call one_to_two(c1,r1, outlets_list(q)) ! �������� ������� ��������� ����� ������
        if (depr(c1,r1) /= -2) cycle ! ������������ ����� ������, ���� ��� ��� �� ����� ������
        call filling_ASTAR(Z,depr,c1,r1,Z(c1,r1),Zout_temp,Nx,Ny,Zmax,StepX,StepY,NODATA)
        do i = 1,Nx
            do j = 1,Ny
                ! ���� ������ � ��������������� ������ ������, ��� �� ����������, ��������� �������� �� ��������� ������� �� ����������
                ! �� ������� ���������� �������, ��� ������ ��������� ��������� ��������� (NODATA)
                ! �������, ���� � ������� ��� �������� ������, ��� ������ ����������� �� ������������
                Select case (depr(i,j))
                case (0) ! �����, �� ����������� � ���������, ������� ��� ������
                    cycle
                case (-1) ! �������
                    if (Zout(i,j) < Zout_temp(i,j) .or. Zout(i,j) == NODATA) Zout(i,j) = Zout_temp(i,j)
                case (-2) ! �����
                    if (Zout(i,j) > Zout_temp(i,j) .or. Zout(i,j) == NODATA) Zout(i,j) = Zout_temp(i,j)
                case default ! ���������
                    if (Zout(i,j) > Zout_temp(i,j) .or. Zout(i,j) == NODATA) Zout(i,j) = Zout_temp(i,j)
                end select
            enddo
        enddo
        
        ! ����������� �������� ����������
        percent_complete = int(real(q)/real(number_of_outlets) * 100)
        if (percent_complete > percent_0) then
            print *, "Filling:", percent_complete, "% completed"
            percent_0 = percent_complete
        endif 
        
    enddo
    
    CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
    write(*,'(A,i2,A,i2,a,i2)') "Filling finished at:    ", time_values(5), ":", time_values(6), ":", time_values(7)
    
    ! �� ������ ����� ������� Zout, ���������� ����������� ���. ��� ������������� ������� �������
    if(allocated(Z_flat)) deallocate(Z_flat)
    if(allocated(depr)) deallocate(depr)
    if(allocated(Zout_temp)) deallocate(Zout_temp)
    if(allocated(outlets_list)) deallocate(outlets_list)
    
    contains

    ! ������������ ������ �����-������
    subroutine find_borders()
    
    ! ���������� ���� 3�3 (������� �� ����������� "������" ������ ������� �������)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    

    ! ����� �� ������. ��������������� ��� ������ �������. ���� ���� �� ���� �� ������� ��������� � ���������, �� ������ ��������� ��������
    do c1 = 1, Nx
        do r1 = 1 ,Ny
            if (depr(c1,r1) /= 0) cycle
            do k = 1,8
                c2 = c1+kx(k); r2 = r1+ky(k)
                if (c2<1 .or. c2> Nx .or. r2 <1 .or. r2> Ny) cycle
                if (depr(c2,r2) > 0) depr(c1,r1) = -1
            enddo
        enddo
    enddo

    end subroutine find_borders
    
integer function two_to_one(column,row)
integer :: row,column
two_to_one = (column-1)*Ny+row
end function two_to_one

subroutine one_to_two(column,row,index)
integer :: column,row,index
row = mod(index,Ny)
if (row == 0) then ! ���� �� ������ ����� ��������, ������ ����� ������������ �����������
    row = Ny
    index = index - Ny
endif
column = index/Ny + 1    
end subroutine one_to_two

    end subroutine filling