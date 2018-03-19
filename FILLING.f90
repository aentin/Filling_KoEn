! Подпрограмма FILLING: заполнение локальных понижений с использованием алгоритма A* и линейной интерполяции высот вдоль ломаных
!   Подпрограмма включает все предварительные и подготовительные операции, необходимые собственно для заполнения 
!   Входные параметры:
!   * Z - исходная матрица высот
!   * Nx, Ny - количество строк/столбцов в матрице
!   * Zmax - максимальное значение высоты на исходной матрице Z
!   * StepX, StepY - шаг сетки по координатам
!   * NODATA - значение «нет данных» на исходной и результирующией матрице
!   Выходные параметры:
!   * Zout - результирующая (заполненная) матрица высот. Размерность полностью аналогична размерности матрицы Z
!
! Внешние подпрограммы:
!   Wang-Liu — подпрограмма для заполнения локальных понижений «под бокс» по алгоритму Ванга-Лю (результат заполнения используется для определения понижений и точек выхода)
!   dist_to_borders_Z_matrix — подпрограмма для вычисления расстояний от ячеек понижения до ячеек-границ
!   matrix_Dijkstra — непосредственно подпрограмма заполнения понижений с использованием алгоритма Дейкстры. Работает по принципу «одно понижение — один выход»
! Все прочие необходимые подпрограммы размещены в блоке CONTAINS
!
! Порядок работы:
! 1) Заполнение понижений по алгоритму Ванга-Лю (с получением горизонтальных поверхностей) и одновременное определение их местонахождений;
! 2) Подготовка промежуточной матрицы для переноса высот
! 3) В цикле для каждой пары «понижение — точка выхода»:
!       а) построение новой поверхности для пары «понижение — точка выхода»;
!       б) совмещение новой поверхности с промежуточной матрицей;
! 4) Выход с новой матрицей высот

    subroutine FILLING(Z,Zout,Nx,Ny,Zmin,Zmax,StepX,StepY,NODATA)
    implicit none
    
    ! Переменные, передаваемые в программу
    real :: Z(Nx,Ny) ! Исходная модель
    real :: Zout(Nx,Ny) ! Результат заполнения
    real(4) :: Zmin, Zmax ! Минимальное и максимальное значение высоты на исходной модели
    integer :: Nx,Ny ! Размерность матрицы
    real(8) :: StepX,StepY ! Шаг сетки
    real :: NODATA ! Значение «нет данных»
    
    ! Внутренние переменные подпрограммы
    real, allocatable :: Zout_temp(:,:) ! Промежуточная матрица высот
    real, allocatable :: Z_flat(:,:) ! Для работы программы заполнения нужна модель, заполненная «под плоскость»
    integer, allocatable :: depr(:,:) ! Матрица, ячейки которой содержат информацию о понижениях, границах и точках выхода
    integer, allocatable :: outlets_list(:) ! Список точек выхода
    integer :: number_of_outlets
    ! Итераторы для циклов
    integer :: c1,c2,r1,r2,i,j,q

    ! Для работы с временем
    character(25) :: DATE, ZEIT, ZONE
    integer :: time_values(8)
    
    integer :: percent_complete, percent_0
    
! НАЧАЛО РАБОТЫ
    
    ! 1. Заполнение понижений «под плоскость» (и сразу определение местонахождения понижений)
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
    
    ! Из подпрограммы переданы: «плоская» модель, матрица понижений (целочисленная), список ячеек-выходов
    
    ! 2. Поиск границ
    !call find_borders()
    ! В результирующей матрице: -1 - граница, 0 - не граница
    
    ! 3. Подготовка результирующей матрицы (высоты ячеек понижений заменяются на значение «нет данных»
    do i = 1,Nx
        do j = 1,Ny
            if (depr(i,j) == 0) then
                Zout(i,j) = Z(i,j)
            else
                Zout(i,j) = NODATA
            endif
        enddo
    enddo
    
    ! 4. Собственно, заполнение понижений. Последовательно, по одной паре «понижение — точка выхода»
    CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
    write(*,'(A,i2,A,i2,a,i2)') "Filling started at:    ", time_values(5), ":", time_values(6), ":", time_values(7)
    
    !!!! Инструкция по "выбросу" индексированного растра локальных понижений
    !Zout(1:Nx,1:Ny) = depr(1:Nx,1:Ny)
    !return
    
    percent_0 = 0
    allocate(Zout_temp(Nx,Ny)) ! Промежуточная матрица для одного результата заполнения
    do q = 1, number_of_outlets
        if (outlets_list(q) == 0) exit ! Выходим, когда добрались до пустой части списка
        call one_to_two(c1,r1, outlets_list(q)) ! Получаем индексы очередной точки выхода
        if (depr(c1,r1) /= -2) cycle ! Прокручиваем точку выхода, если она уже не точка выхода
        call filling_ASTAR(Z,depr,c1,r1,Z(c1,r1),Zout_temp,Nx,Ny,Zmax,StepX,StepY,NODATA)
        do i = 1,Nx
            do j = 1,Ny
                ! Если высота в «результирующей» модели больше, чем во «временной», переносим значение из временной матрицы на постоянную
                ! Мы заранее заготовили матрицу, где высоты локальных понижений устранены (NODATA)
                ! Поэтому, если в матрице нет значения высоты, оно просто переносится из подпрограммы
                Select case (depr(i,j))
                case (0) ! Точка, не относящаяся к понижению, границе или выходу
                    cycle
                case (-1) ! Граница
                    if (Zout(i,j) < Zout_temp(i,j) .or. Zout(i,j) == NODATA) Zout(i,j) = Zout_temp(i,j)
                case (-2) ! Выход
                    if (Zout(i,j) > Zout_temp(i,j) .or. Zout(i,j) == NODATA) Zout(i,j) = Zout_temp(i,j)
                case default ! Понижение
                    if (Zout(i,j) > Zout_temp(i,j) .or. Zout(i,j) == NODATA) Zout(i,j) = Zout_temp(i,j)
                end select
            enddo
        enddo
        
        ! Отображение процента выполнения
        percent_complete = int(real(q)/real(number_of_outlets) * 100)
        if (percent_complete > percent_0) then
            print *, "Filling:", percent_complete, "% completed"
            percent_0 = percent_complete
        endif 
        
    enddo
    
    CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
    write(*,'(A,i2,A,i2,a,i2)') "Filling finished at:    ", time_values(5), ":", time_values(6), ":", time_values(7)
    
    ! На выходе имеем матрицу Zout, содержащую заполненную ЦМР. Все промежуточные матрицы удаляем
    if(allocated(Z_flat)) deallocate(Z_flat)
    if(allocated(depr)) deallocate(depr)
    if(allocated(Zout_temp)) deallocate(Zout_temp)
    if(allocated(outlets_list)) deallocate(outlets_list)
    
    contains

    ! Подпрограмма поиска ячеек-границ
    subroutine find_borders()
    
    ! Скользящее окно 3х3 (поворот от направления "вправо" против часовой стрелки)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    

    ! Здесь всё просто. Просматриваются все ячейки матрицы. Если хотя бы один из соседей относится к понижению, то ячейка считается границей
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
if (row == 0) then ! Если не ввести такую поправку, индекс будет возвращаться неправильно
    row = Ny
    index = index - Ny
endif
column = index/Ny + 1    
end subroutine one_to_two

    end subroutine filling