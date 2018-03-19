program FLOW

implicit none

! Переменные, описывающие ЦМР
real, allocatable :: Z(:,:),ZX(:),ZY(:),Fx(:,:),Fy(:,:),Slope(:,:)
real(8)           :: Xmin,Ymin,Xmax,Ymax,StepX,StepY,SizeX,SizeY
integer           :: Nx,Ny
real(8)           :: floweps
real(4)           :: Zmin, Zmax


type tGRDHeader
    sequence
    character  Signature*4 ! DSBB
    integer(2) Nx;
    integer(2) Ny;
    real(8)    Xmin;
    real(8)    Xmax;
    real(8)    Ymin;
    real(8)    Ymax;
    real(8)    Zmin;
    real(8)    Zmax;
end type

integer(4), parameter  :: GRD_NODATA_HEX  = Z'7EFFFFEE' ! "Нет данных" на сёрферном языке
character*4, parameter :: GRD_SIGNATURE = 'DSBB'
integer, parameter     :: GRD_FMT = 0
real, parameter        :: GRD_NODATA = transfer(GRD_NODATA_HEX,0.0)

! объявление переменной производного типа (GRD_Header - переменная, tGRDHeader - тип, объявленный выше)
type(tGRDHeader)       :: GRD_Header 

real, allocatable   :: Zout(:,:) ! массивы для обработки и вывода данных

real                    :: nodata

integer                 :: i,j,ii,jj,ic,k,kk,k1,k2,row,column     ! переменные для разного рода циклов
integer                 :: ier           ! код ошибки
character(128)          :: directory, inputfile, outputfile  ! Директория, имя входного и выходного файла. Все имеют по 128 символов

! Для работы с временем
real     :: current_time, previous_time


! НАЧАЛО РАБОТЫ
write(*,'(A)') 'Program FILLING searches closed local depressions and fills them with linear interpolation' 
write(*,'(A)') 'To fill a model, write its full filename (like D:\Projects\models\etc\dem.grd) in the field below and press Enter'
write(*,'(A)') "Be sure if your file is in Surfer grid (*.grd) format with Surfer's NODATA value, 'cause the program does not currently support other formats"
write(*,'(A)') 'The input is not case-sensitive, e.g. DEM.GRD = dem.grd'
write(*,'(A)') 'And, unfortunately, the program cannot read Cyrillic pathways'

!directory = "C:\DATA\TEST\"              ! Директория, в которой будем работать ! ИСПРАВИТЬ!!!
write(*,'(A\)') 'Enter input grid file name (with extension): '
!read(*,'(A)') inputfile
!inputfile = "D:\DATA\FILLING\ex.grd"
inputfile = "D:\DATA\FILLING\khoper.grd"

! Чтение файла

open(file = inputfile,unit=20,status='old',form='binary')
    read(20) GRD_Header
    if(GRD_header.signature == GRD_SIGNATURE) then  ! Проверка соответствия формата
        Nx = GRD_Header%Nx 
        Ny = GRD_Header%Ny
        allocate(Z(Nx,Ny),stat=ier); 
        if(ier/=0) goto 9003 ! Прекращение работы программы в случае отсутствия нужного файла
        do j=1,Ny
            read(20) (Z(i,j),i=1,Nx)
        enddo
    endif
close(20)
! Cообщаем об окончании чтения файла
print *, ' Reading completed '

! Рассчитываем недостающие переменные
Xmin = GRD_Header%Xmin
Xmax = GRD_Header%Xmax
Ymin = GRD_Header%Ymin
Ymax = GRD_Header%Ymax
Zmax = GRD_Header%Zmax
Zmin = GRD_Header%Zmin
StepX = (Xmax-Xmin) / (Nx - 1)
StepY = (Ymax-Ymin) / (Ny - 1)
nodata = GRD_NODATA

! Рабочий блок
    
allocate (Zout(Nx,Ny))
Zout(1:Nx,1:Ny) = nodata
call filling(Z,Zout,Nx,Ny,Zmin,Zmax,StepX,StepY,GRD_NODATA)
print *, "Ready to drop filled DEM"
call save_grd(Zout,Nx,Ny,Xmin,Ymin,Xmax,Ymax)

! Переход к концу программы
goto 9000

! Возможные ошибки
    9001 print *,' *** ERROR *** : memory allocation'; goto 9000
    9002 print *,' *** ERROR *** : illegal input file format'; goto 9000
    9003 print *,' *** ERROR *** : input grid file does not exist'; goto 9000
        
! Продолжаем после получения сообщения об ошибке
9000 Continue

call sleep (10) ! УДАЛИТЬ!

! Освободить память
if(allocated(Z)) deallocate(Z)
if(allocated(Zx)) deallocate(Zx)
if(allocated(Zy)) deallocate(Zy)


contains ! Всопомогательные подпрограммы для выполнения отдельных операций

subroutine save_grd(Zout,Nx,Ny,Xmin,Ymin,Xmax,Ymax) ! Сохранение grd-файла

type tGRDHeader
    sequence
    character  Signature*4 ! DSBB
    integer(2) Nx;
    integer(2) Ny;
    real(8)    Xmin;
    real(8)    Xmax;
    real(8)    Ymin;
    real(8)    Ymax;
    real(8)    Zmin;
    real(8)    Zmax;
end type

real :: Zout(Nx,Ny)
real(8) :: Xmin,Ymin,Xmax,Ymax
integer :: Nx, Ny
real(4) :: Zmin, Zmax
type(tGRDHeader) :: GRD_Header 
character(128)    :: outputfile, directory

GRD_Header%Signature = 'DSBB'
GRD_Header%Nx = Nx
GRD_Header%Ny = Ny
GRD_Header%Xmin = Xmin
GRD_Header%Xmax = Xmax
GRD_Header%Ymin = Ymin
GRD_Header%Ymax = Ymax
GRD_Header%Zmin = minval(Zout(1:Nx,1:Ny))
GRD_Header%Zmax = maxval(Zout(1:Nx,1:Ny))

!directory = "C:\DATA\TEST\"              ! Директория, в которой будем работать ! ИСПРАВИТЬ!!!
write(*,'(A\)') 'Enter output grid filename (with extension): '
!read(*,'(A)') outputfile
outputfile = "D:\DATA\FILLING\out.grd"
!outputfile = trim(directory) // trim(outputfile) !Исправить!

    open(file = outputfile, unit=20,status='unknown',form='binary')
    write(20) GRD_Header
    do j=1,Ny
        write(20) Zout(1:Nx,j)
    enddo
    close (20)
    print *, ' Writing completed '  
end subroutine save_grd

! Работа с временем

real function timestamp()

character(25) :: DATE, ZEIT, ZONE
integer :: time_values(8)
real :: timesec1, timesec2

CALL DATE_AND_TIME(DATE,ZEIT,ZONE,time_values)
timestamp = 3600*time_values(5)+ 60*time_values(6) + time_values(7) + 0.001 * time_values(8)

end function timestamp

end program FLOW
