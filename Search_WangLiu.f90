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

    !! Внутренние переменные подпрограммы
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
    
    type (queue) :: q
    type (node)  :: x
    
! Начало работы
    ! Подготовка необходимых переменных и массивов.
    q_max = Nx*Ny
    Z_flat(1:Nx,1:Ny) = Z(1:Nx,1:Ny)
    
    depr(1:Nx,1:Ny) = 0
    allocate(mask(Nx,Ny))
    mask(1:Nx,1:Ny) = 0
    
    percent_0 = 0
    
    ! Определение первоначального списка точек
    q2 = 0
    call initialize_set()
    
    ! Обработка ЦМР
    q_out = 0
    do while (q%n >0)
        ! Извлекаем первый элемент из списка
        x = q%top()
        Z1 = x%priority; c1 = x%c; r1 = x%r
        do k = 1,8
            c2 = c1+kx(k); r2 = r1 + ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle ! Если вышли за границы модели, прокручиваем
            if (mask(c2,r2) /= 0) cycle ! Пропускаются соседи, у которых значение mask не равно нулю (просмотренные ранее, Nodata, etc.)
            call q%enqueue(Z(c2,r2), c2, r2)! Добавляем в очередь
            mask(c2,r2) = 1 ! Устанавливаем маску
            
            ! Проверяем, есть ли в соседе понижение
            if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
                Z_flat(c2,r2) = Z_flat(c1,r1)
                depr(c2,r2) = 1
                ! Если точка, из которой мы пришли, не является понижением, маркируем её как точку выхода
                if( depr(c1,r1) == 0) then 
                    depr(c1,r1) = -2
                    q_out = q_out + 1
                    out_list(q_out) = (c1 - 1) * Ny + r1
                endif
            endif
        enddo
            
        mask(c1,r1) = 2
            
        ! Отображение процента выполнения
        q2 = q2 + 1
        percent_complete = int(real(q2)/real(q_max) * 100)
        if (percent_complete > percent_0) then
            print *, "Searching for pits:", percent_complete, "% completed"
            percent_0 = percent_complete
        endif 
            
        
            
    !    if (q2 > q_max) exit
    !    call one_to_two(c1,r1,list_cells(q2)) ! Выбираем ячейку для просмотра на текущем шаге
    !
    !    do k = 1,8 ! Идём в гости к соседям
    !        c2 = c1+kx(k); r2 = r1+ky(k)
    !        if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
    !        if (mask(c2,r2) /= 0) cycle ! Если у этого соседа мы уже были, переходим на следующую итерацию
    !        
    !        ! Добавляем новонайденного соседа в список
    !        call add_cell(c2,r2,q2)
    !                   
    !        ! Проверяем, есть ли в соседе понижение
    !        if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
    !        Z_flat(c2,r2) = Z_flat(c1,r1)
    !        depr(c2,r2) = 1
    !        
    !        if( depr(c1,r1) == 0) then ! Если точка, из которой мы пришли, не является понижением, маркируем её как точку выхода
    !            depr(c1,r1) = -2
    !            q_out = q_out + 1
    !            out_list(q_out) = two_to_one(c1,r1)
    !        endif
    !        
    !        endif
    !    enddo
    !    mask(c1,r1) = 2
    !    q2 = q2+1
    !    ! Отображение процента выполнения
    !    percent_complete = int(real(q2)/real(q_max) * 100)
    !    if (percent_complete > percent_0) then
    !        print *, "Searching for pits:", percent_complete, "% completed"
    !        percent_0 = percent_complete
    !    endif 
    enddo
    if(allocated(mask)) deallocate(mask)

contains

subroutine initialize_set()
! Идея такая:
! Просматриваются все ячейки ЦМР. По умолчанию считается, что ячейку не следует добавлять в первоначальный список (decision = .false.)
! Однако, если ячейка находится на границе, или у неё сосед "нет данных" — решение меняется на положительное (decision = .true.)
! Ячейки со значениями "нет данных" в список не заносятся, а сразу убираются из рассмотрения
! В конце подпрограммы, если решение положительное, ячейка добавляется в первоначальный список

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
            call q%enqueue(Z(c1,r1), c1, r1)
        endif     
    enddo
enddo
end subroutine initialize_set

end subroutine Wang_Liu2