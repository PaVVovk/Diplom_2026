module bolsig_data
    implicit none

    ! Тип данных для хранения результатов
    type :: bolsig_data_t
        integer :: n_points                         ! количество точек данных
        real, allocatable :: E_over_N(:)            ! E/N (Td)
        real, allocatable :: mean_energy(:)         ! A1 - средняя энергия (eV)
        real, allocatable :: mobility_N(:)          ! A2 - подвижность * N (1/m/V/s)
        real, allocatable :: diffusion_N(:)         ! A6 - диффузия * N (1/m/s)
        real, allocatable :: energy_mobility_N(:)   ! A11 - энерг. подвижность * N
        real, allocatable :: energy_diffusion_N(:)  ! A12 - энерг. диффузия * N
        real, allocatable :: ionization_freq_N(:)   ! A16 - частота ионизации * N (m3/s)
        character(len=256) :: source_file           ! исходный файл данных
    end type bolsig_data_t

contains

    !----------------------------------------------------------------------
    ! Основная функция чтения файла
    !----------------------------------------------------------------------
    function read_bolsig_file(filename) result(b_data)
        implicit none
        character(len=*), intent(in) :: filename
        type(bolsig_data_t) :: b_data
        integer :: file_unit, io_status, i, n_rows
        character(len=2000) :: line
        logical :: found_table
        real :: dummy

        ! Инициализация
        b_data%n_points = 0
        b_data%source_file = filename
        found_table = .false.

        ! Открываем файл
        open(newunit=file_unit, file=filename, status='old', action='read', iostat=io_status)
        if (io_status /= 0) then
            print *, 'ERROR: Cannot open file ', trim(filename)
            return
        end if

        ! ================================================================
        ! ПЕРВЫЙ ПРОХОД: Чтение шапки и подсчёт количества строк в таблице
        ! ================================================================
        n_rows = 0

        do
            read(file_unit, '(a)', iostat=io_status) line
            if (io_status /= 0) exit

            ! Ищем начало таблицы
            if (index(line, 'R#') > 0 .and. index(line, 'E/N') > 0) then
                found_table = .true.
                exit
            end if
        end do

        if (.not. found_table) then
            print *, 'ERROR: Table header not found in ', trim(filename)
            close(file_unit)
            return
        end if

        ! Подсчитываем количество строк данных
        do
            read(file_unit, *, iostat=io_status) dummy
            if (io_status /= 0) exit
            n_rows = n_rows + 1
        end do

        if (n_rows == 0) then
            print *, 'WARNING: No data rows in ', trim(filename)
            close(file_unit)
            return
        end if

        ! ================================================================
        ! ВТОРОЙ ПРОХОД: Выделение памяти и чтение данных
        ! ================================================================

        ! Выделяем память под массивы
        b_data%n_points = n_rows
        allocate(b_data%E_over_N(n_rows))
        allocate(b_data%mean_energy(n_rows))
        allocate(b_data%mobility_N(n_rows))
        allocate(b_data%diffusion_N(n_rows))
        allocate(b_data%energy_mobility_N(n_rows))
        allocate(b_data%energy_diffusion_N(n_rows))
        allocate(b_data%ionization_freq_N(n_rows))

        ! Возвращаемся в начало файла
        rewind(file_unit)

        ! Пропускаем строки до заголовка таблицы
        do
            read(file_unit, '(a)', iostat=io_status) line
            if (io_status /= 0) exit
            if (index(line, 'R#') > 0 .and. index(line, 'E/N') > 0) exit
        end do

        ! Читаем данные
        do i = 1, n_rows
            read(file_unit, *, iostat=io_status) dummy, &
                b_data%E_over_N(i),              &  ! колонка 2
                b_data%mean_energy(i),           &  ! колонка 3
                b_data%mobility_N(i),            &  ! колонка 4
                b_data%diffusion_N(i),           &  ! колонка 5
                b_data%energy_mobility_N(i),     &  ! колонка 6
                b_data%energy_diffusion_N(i),    &  ! колонка 7
                dummy, dummy,                  &  ! колонки 8, 9 (A13, A14)
                b_data%ionization_freq_N(i)         ! колонка 10 (A16)

            if (io_status /= 0) then
                print *, 'ERROR reading row', i, 'in', trim(filename)
                exit
            end if
        end do

        close(file_unit)
    end function read_bolsig_file

    !----------------------------------------------------------------------
    ! Подпрограмма для освобождения памяти
    !----------------------------------------------------------------------
    subroutine free_bolsig_data(b_data)
        implicit none
        type(bolsig_data_t), intent(inout) :: b_data

        if (allocated(b_data%E_over_N)) deallocate(b_data%E_over_N)
        if (allocated(b_data%mean_energy)) deallocate(b_data%mean_energy)
        if (allocated(b_data%mobility_N)) deallocate(b_data%mobility_N)
        if (allocated(b_data%diffusion_N)) deallocate(b_data%diffusion_N)
        if (allocated(b_data%energy_mobility_N)) deallocate(b_data%energy_mobility_N)
        if (allocated(b_data%energy_diffusion_N)) deallocate(b_data%energy_diffusion_N)
        if (allocated(b_data%ionization_freq_N)) deallocate(b_data%ionization_freq_N)

        b_data%n_points = 0
    end subroutine free_bolsig_data

    !----------------------------------------------------------------------
    ! Функция для интерполяции значения при заданном E/N
    !----------------------------------------------------------------------
!    function interpolate_value(b_data, target_EN, value_type) result(res)
!        implicit none
!        type(bolsig_data_t), intent(in) :: b_data
!        real, intent(in) :: target_EN
!        character(len=*), intent(in) :: value_type
!        real :: res
!        integer :: i
!        real, pointer :: values(:)
!
!        ! Выбираем нужный массив значений
!        select case(trim(value_type))
!            case('mean_energy', 'energy')
!                values => b_data%mean_energy
!            case('mobility')
!                values => b_data%mobility_N
!            case('diffusion')
!                values => b_data%diffusion_N
!            case('energy_mobility')
!                values => b_data%energy_mobility_N
!            case('energy_diffusion')
!                values => b_data%energy_diffusion_N
!            case('ionization')
!                values => b_data%ionization_freq_N
!            case default
!                print *, 'ERROR: Unknown value_type: ', trim(value_type)
!                result = 0.0
!                return
!        end select
!    end function interpolate_value


    !----------------------------------------------------------------------
    ! Подпрограмма для вывода информации о данных
    !----------------------------------------------------------------------
    subroutine print_bolsig_info(b_data)
        implicit none
        type(bolsig_data_t), intent(in) :: b_data
        integer :: i

        print *, ''
        print *, '========================================'
        print *, 'BOLSIG+ Data Information'
        print *, '========================================'
        print *, 'Source file:    ', trim(b_data%source_file)
        print *, 'Number of points:', b_data%n_points
        print *, 'Data range:'
        print '(a, f10.2, a, f10.2, a)', '  E/N:        ', &
            b_data%E_over_N(1), ' -', b_data%E_over_N(b_data%n_points), ' Td'
        print '(a, f10.3, a, f10.3, a)', '  Mean energy:', &
            minval(b_data%mean_energy), ' -', maxval(b_data%mean_energy), ' eV'
        print '(a, es10.3, a, es10.3, a)', '  Mobility:   ', &
            minval(b_data%mobility_N), ' -', maxval(b_data%mobility_N), ' 1/m/V/s'
        print *, '========================================'
        print *, ''

    end subroutine print_bolsig_info

end module bolsig_data
