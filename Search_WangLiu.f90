subroutine Wang_Liu2(Z,Z_flat,depr,out_list,q_out,Zmax,Nx,Ny,NODATA)
    
    use priority_queue_mod
    implicit none
    ! Переменные, передаваемые в подпрограмму
    real :: Z(Nx,Ny) ! Исходная ЦМР
    real :: Z_flat(Nx,Ny) ! Результат заполнения (под плоскость)
    integer :: depr(Nx,Ny) ! Идентификатор плоскостей на модели. 1 - плоскость, 0 - не плоскость.
    integer :: out_list(Nx*Ny) ! Список точек выхода (одномерный, сортировка по возрастанию высоты)
    integer :: q_out ! Итератор для списка точек выхода
    real :: Zmax ! Максимальная высота в пределах модели
    integer :: Nx,Ny ! Размерность сетки
    real :: NODATA ! Значение «нет данных»
    ! Внутренние переменные подпрограммы
    type (queue) :: q ! Очередь с приоритетом
    type (node)  :: x ! Элементарный узел очереди с приоритетом
    integer, allocatable :: mask(:,:) ! Массив, в котором будет отмечаться просмотренность/непросмотренность ячейки. 
    ! 0 - не просмотрена, 1 - в очереди, 2 - просмотрена
    integer :: q2, q_max ! Общее количество точек для обработки. 
    ! В начале работы программы равно количеству ячеек в матрице, затем снижается (за счёт значений "нет данных")
    integer :: percent_complete, percent_0
    integer :: c1,c2,r1,r2
    real :: Z1,Z2
    
    !! Скользящее окно 3х3 (поворот от направления "вправо" против часовой стрелки)
    !integer :: k
    !integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    !integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    
    ! Скользящее окно 3х3 (сначала просмотр прямых соседей, потом соседей по диагонали)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
! Начало работы
    ! Подготовка необходимых переменных и массивов.
    q_max = Nx*Ny ! Счётчик максимального количества точек для обработки
    Z_flat(1:Nx,1:Ny) = Z(1:Nx,1:Ny) ! Модель для заполнения "под плоскость"
    depr(1:Nx,1:Ny) = 0 ! Матрица меток понижений (границ, точек выхода)
    allocate(mask(Nx,Ny)) ! Метки процессинга
    mask(1:Nx,1:Ny) = 0 
    
    ! Определение первоначального списка точек
    q2 = 0
    call initialize_set()
    
    ! Поиск локальных понижений
    q_out = 0
    percent_0 = 0
    do while (q%n >0)
        x = q%top() ! Извлечение первого элемента из очереди
        Z1 = x%priority; c1 = x%c; r1 = x%r
        do k = 1,8 ! Просмотр соседей извлечённого элемента
            c2 = c1+kx(k); r2 = r1 + ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
            if (Z(c2,r2) == nodata) cycle ! 
            
            ! Проверка границы
            if (Z_flat(c2,r2) >= Z_flat(c1,r1) .and. depr(c1,r1) == 1 .and. depr(c2,r2) == 0) then 
                depr(c2,r2) = -1
                continue
            endif
            
            ! Проверка дополнительных точек выхода
            if (mask(c2,r2) == 1 .and. depr(c1,r1) == 1 .and. depr(c2,r2) == -1 .and. Z_flat(c2,r2) == Z_flat(c1,r1)) then
                depr(c2,r2) = -2
                q_out = q_out + 1
                out_list(q_out) = (c1 - 1) * Ny + r1
            endif
            
            if (mask(c2,r2) /= 0) cycle ! Пропускаются соседи, которые уже «засветились» в очереди
            
            ! Добавление в очередь
            q2 = q2 + 1
            call q%enqueue(Z(c2,r2), q2, c2, r2)          
            mask(c2,r2) = 1 ! Устанавливаем маску
            
            ! Проверка понижения
            if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
                Z_flat(c2,r2) = Z_flat(c1,r1)
                depr(c2,r2) = 1
                ! Если точка, из которой мы пришли, не является понижением, маркируем её как точку выхода
                if( depr(c1,r1) /= 1) then 
                    depr(c1,r1) = -2
                    q_out = q_out + 1
                    out_list(q_out) = (c1 - 1) * Ny + r1
                endif
            endif
        enddo
            
        mask(c1,r1) = 2
            
        ! Отображение процента выполнения
        percent_complete = int(real(q2)/real(q_max) * 100)
        if (percent_complete > percent_0) then
            print *, "Searching for pits:", percent_complete, "% DEM scanned"
            percent_0 = percent_complete
        endif 
        
    enddo
    if(allocated(mask)) deallocate(mask)

contains

subroutine initialize_set()
! Идея такая:
! Просматриваются все ячейки ЦМР. По умолчанию считается, что ячейку не следует добавлять в первоначальный список (decision = .false.)
! Однако, если ячейка находится на границе, или у неё сосед "нет данных" — решение меняется на положительное (decision = .true.)
! Ячейки со значениями "нет данных" в список не заносятся, а сразу убираются из рассмотрения
! В конце главного цикла, если решение положительное, ячейка добавляется в первоначальный список

logical :: decision

! Просматриваем по очереди все ячейки матрицы
do c1 = 1,Nx
    do r1 = 1,Ny
        
        decision = .false.
        
        ! Если в ячейке нет данных, исключаем её из рассмотрения насовсем
        if (Z(c1,r1) == nodata) then
            mask(c1,r1) = 2
            q_max = q_max - 1
            cycle
        endif
        
        ! Если в ячейке есть данные, присматриваемся к ней повнимательнее
        do k = 1,8
            c2 = c1+kx(k); r2 = r1+ky(k) ! Идем к соседу
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) then ! Если сосед не существует вовсе, принимаем решение занести ячейку в список
                decision = .true.
                exit
            endif
            if (Z(c2,r2) == nodata) then ! Если сосед существует, но в нём "нет данных", принимаем решение занести ячейку в список
                decision = .true.
                exit
            endif
        enddo
        
        ! Добавление ячейки в список
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