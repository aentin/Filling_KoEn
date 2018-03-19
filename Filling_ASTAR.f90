! Программа предназначена для заполнения локальных понижений через своеобразную реализацию алгоритма А*
    
subroutine filling_ASTAR(input_model,depr,cout,rout,Zo,output_model,Nx,Ny,Zmax,StepX,StepY,nodata)
    
    use priority_queue_mod

    implicit none
    
    ! Переменные, передаваемые в подпрограму
    real :: input_model(Nx,Ny) ! Исходная матрица высот
    real :: output_model(Nx,Ny) ! Матрица для записи результата заполнения
    integer :: depr(Nx,Ny) ! Матрица понижений, границ и точек выхода
    integer :: Nx,Ny ! Размерность матриц
    real(8) :: StepX,StepY ! Шаг сетки (в предположении, что шаги по X и Y могут различаться)
    real :: Zmax ! Максимальная высота в пределах исходной модели. Должна быть прописана в свойствах файла
    integer :: cout, rout  ! Колонка и столбец ячейки выхода
    real :: Zo ! Высота ячейки-выхода
    real :: nodata
    logical :: connect_depressions
    
    ! Переменные, работающие в программе
    real, allocatable :: waylength(:,:) ! Массив, в который записывается длина пути. Это не то же самое, что стоимость пути
    integer, allocatable :: way_prev(:,:) ! В этот массив мы пишем, из какой точки мы возвращаемся
    integer, allocatable  :: flag_status(:,:) ! Массив, отражающий «состояние» ячейки: 
    ! 0 — нет сведений, 
    ! 1 — относится к понижению, обработана, 
    ! 2 — относится к границе, поставлена в очередь на «обратный ход»,
    ! 3 — относится к понижению, поставлен в очередь на обработку
    ! -1 — выход из понижения
    integer, allocatable :: list_cells(:) ! Общий список ячеек
    real, allocatable :: list_Z(:) ! Список высот ячеек. ! Параллелен списку list_cells. Необходим для удобного сравнения высоты (не надо каждый раз обращаться к исходной матрице)
    integer, allocatable :: list_borders(:) ! Общий список граничных ячеек. Нужен, чтобы не итерировать каждый раз всю матрицу при построении обратных путей
    
    integer :: c1,c2,c3,r1,r2,r3 ! Переменные для индексов в матрице
    real :: length
    real :: dZ ! Разность высот (используется при расчёте новой высоты точки
    
    ! Первая переменная соответствует индексу, вторая маркирует "просмотренность" (1 - в очереди, 0 - просмотрено)
    integer :: q1, q2, q3 ! Счётчик общего количества точек и индекс текущей точки в списке
    integer :: s,s2 ! Переменные для конвертации строки-столбца в одно целое число (и последующей записи этого числа в список)
    integer :: i,j,k ! Итераторы для циклов
    integer :: i1, j1
    
    ! Переменные для работы с алгоритмом A*
    type (queue) :: queue_points, queue_borders
    type (node)  :: x
    
    
    real(8) :: Z0,Z1,Z2,Z3,Z_temp ! Высота точки истечения, временная новая высота заполняемой точки
    
    ! Переменные для присвоения высоты тем ячейкам, через которые не проходит кратчайший путь
    real :: dist1, dist2
    
    ! Скользящее окно 3х3 (сначала просмотр прямых соседей, потом соседей по диагонали)
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
! НАЧАЛО РАБОТЫ
    
    allocate(waylength(Nx,Ny))
    allocate(way_prev(Nx,Ny))
    allocate(flag_status(Nx,Ny))
    allocate(list_cells(Nx*Ny))
    allocate(list_Z(Nx*Ny))
    allocate(list_borders(Nx*Ny))
    
    connect_depressions = .TRUE.
    
    ! Заполняем выходной массив
    do i = 1,Nx
        do j = 1,Ny
            if (depr(i,j) > 0) then
                output_model(i,j) = nodata
            else
                output_model(i,j) = input_model(i,j)
            endif
        enddo
    enddo
    
    ! Записываем индексы ячейки-выхода в переменные, которые будут использованы в цикле
    c1 = cout; r1 = rout; 
    
    ! Подготовительные операции
    list_cells(1:Nx*Ny) = 0 ! Общий список ячеек
    list_Z(1:Nx*Ny) = nodata ! Общий список высот ячеек
    list_borders(1:Nx*Ny) = 0 ! Общий список точек выхода
    waylength(1:Nx,1:Ny) = nodata

    ! Ставим первую точку в очередь
    flag_status(1:Nx,1:Ny) = 0
    q1 = 1; q2 = 1; q3 = 0
    list_cells(1) = two_to_one(c1,r1); list_Z(1) = input_model(c1,r1); flag_status(c1,r1) = -1
    ! Переменная Z1 хранит вторую высоту, используемую в интерполяции
    Z1 = Zmax ! Конкретное значение высоты будет найдено в процессе работы программы
    waylength(c1,r1) = 0.0
    ! Подготовительный этап на этом закончен
    
    ! Процедура поиска кратчайших путей по алгоритму AT
    do
        if (q2 > 1 .and. list_cells(q2) == 0) exit ! Завершаем поиск кратчайших путей, когда в очереди не остаётся ячеек
        ! То есть на очередном месте в списке-очереди — пустая запись
                
        call one_to_two(c1,r1,list_cells(q2)) ! Выбираем ячейку для просмотра на текущем шаге
        !!!x = queue_points%top()
        
        if (depr(c1,r1) == 0 .or. depr(c1,r1) == -1) cycle ! Пропускаем ячейку из очереди, если она не относится к понижению
        
        do k = 1,8 ! Идём в гости к соседям
            c2 = c1+kx(k); r2 = r1+ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
            if (depr(c2,r2) == 0) cycle ! Если сосед не относится к понижению, прокручиваем его
            if (flag_status(c2,r2) == 2 .or. flag_status(c2,r2) == 3) cycle ! Если сосед уже был обработан ранее, тоже прокручиваем его
            
            ! Вычисляем расстояние для соседа
            waylength(c2,r2) = waylength(c1,r1) + sqrt( ( (c2-c1)*StepX )**2 + ( (r2-r1)*StepY )**2 )
            ! Отмечаем путь возврата для соседа
            way_prev(c2,r2) = two_to_one(c1,r1)
            ! Отмечаем соседа как требующего обработки
            flag_status(c2,r2) = 3
            ! Добавляем новонайденного соседа в список
            Select case (depr(c2,r2)) ! Обработка соседей в зависимости от типа
            case(-2) ! Сосед — точка выхода
                ! Обработка встреченных точек выхода будет зависеть от параметра connect_depressions
                ! Если TRUE, система понижений обрабатывается как единое целое
                ! Если FALSE, каждое понижение будет обрабатываться отдельно
                Select case (connect_depressions)
                case (.true.)
                    if (input_model(c2,r2) > Zo) then
                        output_model(c2,r2) = nodata
                        depr(c2,r2) = 1
                        flag_status(c2,r2) = 3
                        call add_cell_1(c2,r2,q1,q2)
                    endif
                    !!!call q%enqueue(Z(c2,r2), c2, r2)! Добавляем в очередь
                    
                case (.false.)
                    flag_status(c2,r2) = 2
                    if (Z1 > input_model(c2,r2) .and. input_model(c2,r2) > input_model(cout,rout)) then
                    Z1 = input_model(c2,r2)
                    c3 = c2; r3 = r2
                    endif
                end select
            case (-1) ! Сосед — граница
                ! Вставить инструкцию записи в список
                flag_status(c2,r2) = 2
                q3 = q3 + 1
                list_borders(q3) = two_to_one(c2,r2)
                ! Инструкция поиска минимальной высоты среди граничных ячеек. По моему скромному мнению, она не нужна
                if (Z1 > input_model(c2,r2) .and. input_model(c2,r2) > input_model(cout,rout)) then
                    Z1 = input_model(c2,r2)
                    c3 = c2; r3 = r2
                endif
            
            case (0) ! Сосед — просто сосед
                continue
            case default ! Сосед — просто сосед
                flag_status(c2,r2) = 3
                call add_cell_1(c2,r2,q1,q2)  
            end select
        
        enddo
        q2 = q2+1
    enddo
       
    ! Ищем высоты для интерполяции (Z0 и Z1) 
    !! Пытаемся уменьшить высоту точки выхода
    !!call decreasing_outlet()
    Z0 = output_model(cout,rout)
    ! Пытаемся «приподнять» граничные ячейки.
    !!call increasing_border(Z1)
    
    ! К настоящему моменту имеем: предварительно «приподнятую» временную модель; грид расстояний, грид "соединений" (каждая ячейка с предыдущей)
    ! Имееем также «верхнюю» (Z1) и «нижнюю» (Z0) высоты интерполяции
    
    ! Теперь пойдём в обратном направлении (от границ к выходу). Просматривать будем все ячейки границ (flag_status == 2) подряд
    c1 = 0;r1 = 0; c2 = 0; r2 = 0
    q2 = 0
            
    do q2 = 1, q3
        call one_to_two(c1,r1,list_borders(q2))
            
        if (flag_status(c1,r1) /= 2) cycle ! Если точка не относится к границе, прокручиваем её
        ! Пришли в точку, из которой надо восстанавливать линию.            
        dZ = input_model(c1,r1) - Z0
        !dZ = Z1 - Z0
        length = waylength(c1,r1)
        do
            s = way_prev(c1,r1)
            call one_to_two(c2,r2,s)
            if (c2 == cout .and. r2 == rout) exit ! Выходим, если последующая точка является точкой выхода
            Z_temp = Z0 + (dZ/length)*waylength(c2,r2)
            if (output_model(c2,r2) == nodata .or. output_model(c2,r2) > Z_temp) output_model(c2,r2) = Z_temp
            flag_status(c2,r2) = 1
            c1 = c2; r1 = r2
        enddo      
    enddo
    
    flag_status(cout,rout) = 1
    
    ! Заполнение точек, через которые не проходит кратчайший путь
    ! Заполнение производится путём экстраполяции. Используются: высота точки выхода (Z0) и высота ближайшей_уже_обработанной точки (Z2)
    do i = 1, q1
        call one_to_two(c3,r3,list_cells(i))
        if (flag_status(c3,r3) /= 3) cycle
        dist1 = waylength(c3,r3) ! Полный путь от обрабатываемой ячейки до выхода
        c1 = c3; r1 = r3
        ! Теперь ищем ближайшего уже обработанного соседа вниз по течению   
        do
            s = way_prev(c1,r1)
            call one_to_two(c2,r2,s)
            if (flag_status(c2,r2) == -1) then ! Если, двигаясь обратно, мы дошли до ячейки-выхода, то...
                Z2 = output_model(c2,r2) + 0.1
                dist2 = dist1
            else if (flag_status(c2,r2) == 1) then ! Если, двигаясь обратно, мы дошли до точки с уже известной высотой, то...
                Z2 = output_model(c2,r2) ! Устанавливаем новое значение «базовой» высоты (равное высоте заполненной ранее точки)
                dist2 = waylength(c2,r2) ! И запоминаем расстояние от незаполненной точки до выхода
                exit ! И выходим из цикла
            else
                c1 = c2; r1 = r2 ! Если нет, то смотрим дальше
            endif
        enddo
        
        Z3 = Z0 + (Z2-Z0)/dist2*dist1       
        output_model(c3,r3) = Z3
        
        continue
        
    enddo

    ! Собственно, всё. Матрица output_model передаётся обратно в главную программу

    if(allocated(waylength)) deallocate(waylength)
    if(allocated(flag_status)) deallocate(flag_status)
    if(allocated(way_prev)) deallocate(way_prev)
    if(allocated(list_cells)) deallocate(list_cells)
    if(allocated(list_Z)) deallocate(list_Z)
    if(allocated(list_borders)) deallocate(list_borders)

contains

! Подпрограмма добавления ячейки

subroutine add_cell_1(c,r,q1,q2)

integer :: c,r,q1,q2 ! Столбец и строка ячейки, добавляемой в список, а также номер текущего элемента в списке
real(4) :: Z_index
integer :: index
integer :: q
integer :: new_element ! Индекс, под которым будет записан новый элемент. По умолчанию — в конец списка
integer :: first, last
integer :: ncase
    
    index = two_to_one(c,r)

    q1 = q1+1
    new_element = q1 ! Собственно, вот оно, "по умолчанию"
    
    if (q1 > 1) then ! Последующие операции проводятся только тогда, когда в списке больше одного элемента
        if (input_model(c,r) < list_Z(q1-1)) then 
            ! Если новый элемент меньше последнего, запускаем процедуру сдвига.
            ! (Если новый элемент больше последнего, сразу ставим его в конец)
            
            first = q2+1; last = q1-1
            ! Ищем место элемента в списке. Сортировка кучей (точнее, её подобие)
            if (last /= first) then
                do while (last /= first)
                if (input_model(c,r) < list_Z( (last-first)/2 + first) ) then
                    last = (last-first)/2 + first
                else
                    first = (last-first)/2 + first + 1
                endif
                enddo
            endif
            new_element = first
        call move_cells_in_queue_1(new_element,q1)
        endif   
    endif        
    
    ! Записываем новый элемент в оба списка в ту позицию, которая была определена (q, она же new_element)
    list_cells(new_element) = index
    list_Z(new_element) = input_model(c,r)
end subroutine add_cell_1

subroutine move_cells_in_queue_1(first_cell,last_cell)
integer :: number
integer :: first_cell,last_cell

do number = last_cell, first_cell+1, -1
    list_cells(number) = list_cells(number-1)
    list_Z(number) = list_Z(number-1)
enddo

end subroutine move_cells_in_queue_1

! Подпрограмма снижения высоты точки выхода
subroutine decreasing_outlet()
    real(4) :: Z0_temp
    real(4),parameter :: threshold = 0.0001
    
    Z0_temp = output_model(cout,rout)
    k = 0
    do k = 1,8
        c2 = cout+kx(k); r2 = rout+ky(k)
        
        if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
        if (depr(c2,r2) > 0) cycle ! Не рассматриваем ячейки из понижений
        if (output_model(c2,r2) + threshold < Z0_temp) then
            Z0_temp = output_model(c2,r2) + threshold
        endif
    enddo
    if (output_model(cout,rout) > Z0_temp) then
        output_model(cout,rout) = Z0_temp 
        print *, "Outlet point height decreased from ", input_model(cout,rout), " to ", Z0_temp
    else
        continue
        print *, "Outlet point height did not decreased"
    endif

end subroutine decreasing_outlet

! Подпрограмма увеличения граничной точки
subroutine increasing_border(Z1)
    real, allocatable :: matrix_Z1_temp(:,:)
    real(4),parameter :: threshold = 1.0
    real(4) :: Z1
    real(4) :: Z1_temp
    integer :: c,r
    
    ! Подпрограмма работает в три этапа:
    ! 1) Для каждой точки границы, от которой будет идти интерполяция (flag_status = 2), ищется теоретически максимальное значение высоты
    ! Если точка является точкой выхода (не для этого понижения!), то её значение не увеличивается (максимальное значение равно собственному значению)
    ! 2) Среди всех максимальных высот (уже увеличенных!) ищется минимальная;
    ! 3) Минимальная высота принимается за высоту интерполяции, все высоты, которые ниже её (но выше исходной высоты точки выхода!), приподнимаются до её значения    
    
    allocate(matrix_Z1_temp(Nx,Ny)) ! Матрица максимальных высот
    matrix_Z1_temp(1:Nx,1:Ny) = nodata
    Z1_temp = Zmax ! Предполагаемое новое, повышенное значение высоты
    c = 0; r = 0
    ! 1. Определяем «максимальные» высоты точек (для каждой ячейки)
    do c1 = 1,Nx
        do r1 = 1,Ny
            if (flag_status(c1,r1) /= 2) cycle
            matrix_Z1_temp(c1,r1) = input_model(c1,r1)
            if (depr(c1,r1) /= -1) cycle ! Для всех точек, которые «не просто границы», сохраняем существующую высоту
            ! Все точки, которые «просто границы», пытаемся приподнять. Просматриваем для этой цели всех соседей
            do k = 1,8
                c2 = c1+kx(k); r2 = r1+ky(k)
                if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
                if (depr(c2,r2) > 0) cycle ! Соседей из понижений пропускаем
                if (input_model(c2,r2) < input_model(c1,r1) .or. input_model(c2,r2) == nodata) cycle ! И вообще более низких соседей пропускаем
                if (matrix_Z1_temp(c1,r1) < (input_model(c2,r2) - threshold) ) then
                    matrix_Z1_temp(c1,r1) = input_model(c2,r2) - threshold
                endif
            enddo
            ! 2. Ищем «минимум среди максимумов»
            if (Z1_temp > matrix_Z1_temp(c1,r1)) then
                Z1_temp = matrix_Z1_temp(c1,r1); c = c1; r = r1
            endif
        enddo
    enddo
    
    ! Проверка, увеличены ли высоты
    if (Z1_temp == input_model(c,r)) then
        c = 0; r = 0
    endif
    
    Select case (c+r)
    case(0) ! Высоты не увеличены
        print *, "Border point height did not increased"
    case default ! Высоты увеличены
        Z1 = Z1_temp
        do c1 = 1,Nx
            do r1 = 1,Ny
                if (flag_status(c1,r1) /= 2) cycle
                if (depr(c1,r1) /= -1) cycle ! Всё, что не есть граница (не точка выхода!), приподнято не будет
                if (output_model(c1,r1) > Z1) cycle ! Если ячейка уже выше Z1, нечего её увеличивать
                output_model(c1,r1) = Z1 ! Увеличиваем все более низкие высоты до максимума.
            enddo
        enddo
        print *, "Border point height increased from ", input_model(c,r), " to ", Z1
    end select
    
    deallocate(matrix_Z1_temp)
    
end subroutine increasing_border

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

end subroutine filling_ASTAR
