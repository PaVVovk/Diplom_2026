program concentration
implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer :: M, k, step_count, t_counter
    integer, parameter :: T_count = 50
    real(dp) :: r0, delta, h_first, sigma, tau, E_divide_N
    real(dp), allocatable :: r_k(:), r_half_k(:), h_k(:), h_half_k(:)
    real(dp) :: l_e, l_i, D_e, D_i, k_e, k_i, nu_ion
    real(dp) :: p, N, t_max, t, t0

    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: k_b = 1.380649E-23_dp !Дж/К
    real(dp), parameter :: e = 1.602176634E-19_dp !Кл
    real(dp), parameter :: gamma_e = 0.7104_dp
    real(dp), parameter :: gamma_i = 0.7104_dp
    real(dp), parameter :: T_gas = 300_dp  !К
    real(dp), parameter :: beta_ei = 1E-18 ! м^3/c

    real(dp), allocatable :: n_e_i(:), n_i_i(:), E_r_i(:)
    real(dp), allocatable :: n_e_final(:), n_i_final(:), E_r_final(:)
    real(dp) :: t_out(0:T_count)
    logical :: repeat_flag

    character(len=*), parameter :: bolsig_name = "Data"
    character(len=50) :: folder_name

    folder_name = "Data_Fortran\" // data_string()
    call create_folder(folder_name)
    call create_folder(subf("\n_e"))
    call create_folder(subf("\n_i"))
    call create_folder(subf("\j_e"))
    call create_folder(subf("\j_i"))
    call create_folder(subf("\potential"))

    ! Ввод параметров
    print *, 'Enter the number of intervals M:'
    read *, M
    print *, 'Enter the radius of the tube r0:'
    read *, r0
    print *, 'Enter the DELTA thickening parameter:'
    read *, delta
    !print *, 'Enter a numeric parameter sigma'
    !read *, sigma
    sigma = 0.5_dp
    print *, 'Enter the pressure in pascals'
    read *, p


    ! Выделение памяти
    allocate(r_k(0:M))
    allocate(r_half_k(0:M-1))
    allocate(h_k(0:M-1))
    allocate(h_half_k(1:M-1))
    allocate(n_e_i(0:M))
    allocate(n_i_i(0:M))
    allocate(E_r_i(0:M))
    allocate(n_e_final(0:M))
    allocate(n_i_final(0:M))
    allocate(E_r_final(0:M))

    !Значения коэффициентов
    N = p/(k_b*T_gas) ! 1/м^3
    print *, 'N = ', N
    nu_ion = 0.8093E-16_dp*N !1/c
    E_divide_N = 100   !Тд
    k_e = 0.7890E+24_dp/N  ! м^2/(В*с)
    k_i = ((0.286_dp + 0.669_dp*EXP(-E_divide_N/179.5_dp) + 0.679_dp*EXP(-E_divide_N*1305_dp)))*1e-4 ! м^2/(В*с)
    D_e = 0.6932E+25_dp/N  ! м^2/c
    D_i = k_i*T/e  !м^2/c
    l_e = 1/(N*1E-20_dp)
    l_i = 1/(N*1E-20_dp)
    n_i_i = nu_ion/beta_ei
    n_e_i = n_i_i
    E_r_i = 0.0_dp
    t = 0.0_dp
    t_max = 1e-4_dp
    step_count = 0


    ! Генерация основной сетки
    r_k(0) = 0.0_dp
    if(delta == 0) then
        h_first = r0/M
    else
        h_first = (delta * r0) / ((1.0_dp + delta)**M - 1.0_dp)
    end if
    do k = 0, M-1
        h_k(k) = h_first * (1.0_dp + delta)**(k)
        r_k(k+1) = r_k(k) + h_k(k)
        r_half_k(k) = (r_k(k) + r_k(k+1)) / 2.0_dp
    end do

    do k = 1, M-1
        h_half_k(k) = (h_k(k-1) + h_k(k)) / 2.0_dp
    end do
      ! Вывод результатов
    print *, '========================================'
    print *, 'Grid parameters:'
    print *, 'M =', M
    print *, 'r0 =', r0
    print *, 'delta =', delta
    print *, 'Initial step h_first =', h_first
    print *, '========================================'
    print *, 'Main mesh (r_k):'
    do k = 0, M
        print '(I4, F12.6)', k, r_k(k)
    end do

    print *, '========================================'
    print *, 'Main grid steps h_k:'
    do k = 0, M-1
        print '(I4, F12.6)', k, h_k(k)
    end do

    print *, '========================================'
    print *, 'Auxiliary grid (r_{k+1/2}):'
    do k = 0, M-1
        print '(I4, F12.6)', k, r_half_k(k)
    end do

    print *, '========================================'
    print *, 'Auxiliary grid steps h_half_k:'
    do k = 1, M-1
        print '(I4, F12.6)', k, h_half_k(k)
    end do

    ! Проверка последнего узла
    print *, '========================================'
    print *, 'Examination: r(M) =', r_k(M), 'must be =', r0

    print *, 'Min step h_k:', minval(h_k)
    if (minval(h_k) < 1e-30) then
        print *, 'ATTENTION: very small step!'
    end if
    print *, 'n_e_i = ', n_e_i
    print *, 'n_i_i = ', n_i_i
    !Заполняем массив t_out

    t0 = t_max / 10 ** (T_count / 10)
    do k = 0, T_count
        t_out(k) = t0 * 10 ** (k / 10)
        print *, t_out(k)
    end do

    tau = minval(h_k**2/(4*D_e))

    t_counter = 0
    do while(t < t_max)
        repeat_flag = .false.
        step_count = step_count + 1

        if (mod(step_count, 10) == 0) then
            print *, 'Step:', step_count, 't =', t, 'tau =', tau
        end if
        call solve_iterations(n_e_i, n_i_i, E_r_i, &
                           n_e_final, n_i_final, E_r_final, &
                           repeat_flag, tau)
        do while(repeat_flag .eqv. .true.)
            call solve_iterations(n_e_i, n_i_i, E_r_i, &
                           n_e_final, n_i_final, E_r_final, &
                           repeat_flag, tau)
        end do

        if (t > t_out(t_counter)) then
            call write_parameters(n_e_final, n_i_final, n_e_i, n_i_i, E_r_i, t, t_counter)
            t_counter = t_counter + 1
        end if
        !print *, 'n_e_final = ', n_e_final(M/2)
        !print *, 'n_i_final = ', n_i_final(M/2)
        !print *, 'E_r_final = ', E_r_final(M/2)
        n_e_i = n_e_final
        n_i_i = n_i_final
        E_r_i = E_r_final
        t = t + tau
    end do



    print *, '========================================'
    do k = 0, M
       print *, 'Electron concentration n_e =', n_e_i(k)
    end do
    do k = 0, M
       print *, 'Ion concentration n_i =', n_i_i(k)
    end do
    do k = 0, M
       print *, 'Field E_r =', E_r_i(k)
    end do

    deallocate(r_k, r_half_k, h_k, h_half_k, n_e_i, n_i_i, &
               E_r_i, n_e_final, n_i_final, E_r_final)
contains

    !Прогонка для электронов

    subroutine solve_electron_continuity(n_e_i, n_i_im, E_r_im, n_e_m1)
    implicit none
    real(dp), intent(in) :: n_e_i(0:M)
    real(dp), intent(in) :: n_i_im(0:M), E_r_im(0:M)
    real(dp), intent(out) :: n_e_m1(0:M)

    integer :: k
    real(dp) :: A(0:M), B(0:M), C(0:M), F(0:M)
    real(dp) :: alpha(0:M+1), beta(0:M+1)
    real(dp) :: y(0:M)

    real(dp) :: a_k, b_k, u_k, w_k


    B(0) = 4.0_dp*D_e/(h_k(0)**2) + 2.0_dp*k_e*E_r_im(1)/h_k(0)
    C(0) = 1.0_dp/(sigma*tau) - nu_ion + beta_ei * n_i_im(0) + 4.0_dp * D_e / (h_k(0)**2)
    F(0) = n_e_i(0)/tau + (1.0_dp - sigma) * (B(0)*n_e_i(1) - (C(0) - 1.0_dp/(sigma*tau))*n_e_i(0))


    do k = 1, M-1
        a_k = (r_half_k(k-1)/r_k(k)) * (1.0_dp/h_half_k(k))
        b_k = (r_half_k(k+1)/r_k(k)) * (1.0_dp/h_half_k(k))
        u_k = (b_k - a_k) / 2.0_dp
        w_k = b_k / h_k(k) + a_k / h_k(k-1)
        A(k) = a_k * (D_e/h_k(k-1) - k_e*E_r_im(k-1)/2.0_dp)
        B(k) = b_k * (D_e/h_k(k) + k_e*E_r_im(k+1)/2.0_dp)
        C(k) = 1.0_dp/(sigma*tau) - nu_ion + beta_ei * n_i_im(k) + w_k * D_e - u_k * k_e * E_r_im(k)
        F(k) = n_e_i(k)/tau + (1.0_dp - sigma) * (A(k)*n_e_i(k-1) - (C(k) - 1.0_dp/(sigma*tau))*n_e_i(k) + B(k)*n_e_i(k+1))
    end do

    !A(M) = gamma_e*l_e/h_k(M-1)
    !C(M) = 1.0_dp + gamma_e * l_e / h_k(M-1)
    A(M) = 0.0_dp
    C(M) = 1.0_dp
    B(M) = 0.0_dp
    F(M) = 0.0_dp

    ! Метод прогонки

    ! Прямой ход
    alpha(1) = B(0) / C(0)
    beta(1) = F(0) / C(0)
    do k = 1, M-1
    !Сохранить знаменатель в отдельную перем
        alpha(k+1) = B(k) / (C(k) - alpha(k) * A(k))
    end do

    do k = 1, M
        beta(k+1) = (beta(k) * A(k) + F(k)) / (C(k) - alpha(k) * A(k))
    end do
    ! Обратный ход
    y(M) = beta(M+1)
    do k = M-1, 0, -1
        y(k) = alpha(k+1) * y(k+1) + beta(k+1)
    end do

    do k = 0, M
        n_e_m1(k) = y(k) / sigma
    end do

    end subroutine solve_electron_continuity

    !Прогонка для ионов

    subroutine solve_ion_continuity(n_i_i, n_e_im, E_r_im, n_i_m1)
    real(dp), intent(in) :: n_i_i(0:M)
    real(dp), intent(in) :: n_e_im(0:M), E_r_im(0:M)
    real(dp), intent(out) :: n_i_m1(0:M)

    integer :: k
    real(dp) :: A(0:M), B(0:M), C(0:M), F(0:M)
    real(dp) :: alpha(0:M+1), beta(0:M+1)
    real(dp) :: y(0:M)

    real(dp) :: a_k, b_k, u_k, w_k


    B(0) = 4.0_dp*D_i/(h_k(0)**2) - 2.0_dp*k_i*E_r_im(1)/h_k(0)
    C(0) = 1.0_dp/(sigma*tau) + beta_ei*n_e_im(0) + 4.0_dp*D_i/(h_k(0)**2)
    F(0) = n_i_i(0)/tau + nu_ion*n_e_im(0) + &
           (1-sigma)*(B(0)*n_i_i(1) - (C(0) - 1/(sigma*tau))*n_i_i(0))
    do k = 1, M-1
        a_k = (r_half_k(k-1)/r_k(k)) * (1.0_dp/h_half_k(k))
        b_k = (r_half_k(k+1)/r_k(k)) * (1.0_dp/h_half_k(k))
        u_k = (b_k - a_k) / 2.0_dp
        w_k = b_k / h_k(k) + a_k / h_k(k-1)
        A(k) = a_k*(D_i/h_k(k-1) + k_i*E_r_im(k-1)/2)
        B(k) = b_k*(D_i/h_k(k) - k_i*E_r_im(k+1)/2)
        C(k) = 1/(sigma*tau) + beta_ei*n_e_im(k) + w_k*D_i + u_k*k_i*E_r_im(k)
        F(k) = n_i_i(k)/tau + nu_ion*n_e_im(k) + &
               (1 - sigma)*(A(k)*n_i_i(k-1) - &
               (C(k) - 1/(sigma* tau))*n_i_i(k) + B(k)*n_i_i(k+1))
    end do
    !A(M) = gamma_i * l_i / h_k(M-1)
    !C(M) = 1.0_dp + gamma_i * l_i / h_k(M-1)
    A(M) = 0.0_dp
    C(M) = 1.0_dp
    B(M) = 0.0_dp
    F(M) = 0.0_dp

    ! Прямой ход
    alpha(1) = B(0) / C(0)
    beta(1) = F(0) / C(0)
    do k = 1, M-1
        alpha(k+1) = B(k) / (C(k) - alpha(k) * A(k))
    end do

    do k = 1, M
        beta(k+1) = (beta(k) * A(k) + F(k)) / (C(k) - alpha(k) * A(k))
    end do
    ! Обратный ход
    y(M) = beta(M+1)
    do k = M-1, 0, -1
        y(k) = alpha(k+1) * y(k+1) + beta(k+1)
    end do

    do k = 0, M
        n_i_m1(k) = y(k) / sigma
    end do

    end subroutine solve_ion_continuity

    !Решение уравнения Пуассона
    subroutine solve_poisson_trapezoidal(n_e_m1, n_i_m1, E_r_m1)
    implicit none
    real(dp), intent(in) :: n_e_m1(0:M), n_i_m1(0:M)
    real(dp), intent(out) :: E_r_m1(0:M)

    integer :: k
    real(dp) :: F(1:M), const
    const = 2.0_dp * pi * e
    E_r_m1(0) = 0
    E_r_m1(1) = const*h_k(0)*((n_i_m1(1)-n_e_m1(1)) + &
                 (n_i_m1(0) - n_e_m1(0)))
    F(1) = r_k(1)*E_r_m1(1)
    do k = 2, M
        F(k) = F(k-1) + const*h_k(k-1)*(r_k(k-1)*(n_i_m1(k-1)-n_e_m1(k-1))+ &
               r_k(k)*(n_i_m1(k) - n_e_m1(k)))
        E_r_m1(k) = F(k)/r_k(k)
    end do
    end subroutine solve_poisson_trapezoidal

    !Итерации

    subroutine solve_iterations(n_e_i, n_i_i, E_r_i, &
                                n_e_final, n_i_final, E_r_final, &
                                repeat_flag, tau)

    ! Входные параметры

    real(dp), intent(in) :: n_e_i(0:M), n_i_i(0:M), E_r_i(0:M)
    logical, intent(inout) :: repeat_flag
    real(dp), intent(inout) :: tau

    ! Выходные параметры

    real(dp), intent(out) :: n_e_final(0:M), n_i_final(0:M), E_r_final(0:M)

    !Локальные переменные

    integer :: iter, k, j
    real(dp) :: error_ne, error_ni, error_E
    real(dp) :: n_e_m(0:M), n_i_m(0:M), E_r_m(0:M)
    real(dp) :: n_e_im(0:M), n_i_im(0:M), E_r_im(0:M)
    real(dp) :: n_e_m1(0:M), n_i_m1(0:M), E_r_m1(0:M)
    real(dp) :: max_error, tol

    n_e_m = n_e_i
    n_i_m = n_i_i
    E_r_m = E_r_i

    tol = 10E-8
    iter = 10


    do k = 1, iter
        max_error = 0.0_dp
        do j = 0, M
            n_e_im(j) = sigma * n_e_m(j) + (1.0_dp - sigma) * n_e_i(j)
            n_i_im(j) = sigma * n_i_m(j) + (1.0_dp - sigma) * n_i_i(j)
            E_r_im(j) = sigma * E_r_m(j) + (1.0_dp - sigma) * E_r_i(j)
        end do
        call solve_electron_continuity(n_e_i, n_i_im, E_r_im, n_e_m1)
        call solve_ion_continuity(n_i_i, n_e_im, E_r_im, n_i_m1)
        call solve_poisson_trapezoidal(n_e_m1, n_i_m1, E_r_m1)

         do j = 0, M
            error_ne = abs(n_e_m1(j) - n_e_m(j)) / abs(n_e_m1(j))
            error_ni = abs(n_i_m1(j) - n_i_m(j)) / abs(n_i_m1(j))
            error_E =  abs(E_r_m1(j) - E_r_m(j)) / abs(E_r_m1(j))
            max_error = max(error_ne, error_ni, error_E, max_error)
        end do
        !do j = 0, M
            !print *, 'n_e_m1= ', n_e_m1(j)
            !print *, 'n_i_m1= ', n_i_m1(j)
            !print *, 'E_r_m1= ', E_r_m1(j)
            !print *, 'error_ne= ', error_ne
            !print *, 'error_ni= ', error_ni
            !print *, 'error_E= ',  error_E
            !print *, 'max_error= ', max_error
        !end do
        do j = 0, M
            n_e_m(j) = n_e_m1(j)
            n_i_m(j) = n_i_m1(j)
            E_r_m(j) = E_r_m1(j)
        end do
        if (max_error < tol) then
            n_e_final = n_e_m1
            n_i_final = n_i_m1
            E_r_final = E_r_m1
            exit
        end if
        if (k == iter) then
            tau = tau/2
            print *, 'tau = ', tau
            repeat_flag = .true.
        else
            repeat_flag = .false.
        end if
    end do
    end subroutine solve_iterations

    real(dp) function df_dr(f1, f2, delta)
        real(dp) :: f1, f2, delta
        df_dr = (f2 - f1)/delta
        !df_dr = (delta2*(f2 - f1) / delta1  + delta1*(f3 - f2) / delta2) / ( delta1 + delta2)
    end function df_dr

    subroutine create_folder(folder_name)
        character(len=*), intent(in) :: folder_name
        character(len=10), parameter :: MKDIR = 'mkdir '
        integer :: folder_status

        folder_status = system(MKDIR // folder_name)
        if (folder_status /= 0) then
            print *, 'mkdir: failed to create folder!'
        end if
    end subroutine

    subroutine write_parameters(n_e_final, n_i_final, n_e_i, n_i_i, E_r_i, time, ind)

        integer, intent(in) :: ind
        real(dp), intent(in) :: time
        real(dp), intent(in) :: n_e_final(0:M), n_i_final(0:M)
        real(dp), intent(in) :: n_e_i(0:M), n_i_i(0:M), E_r_i(0:M)
        real(dp) :: potential(0:M), j_e(0:M), j_i(0:M), difference_e(1:M-1), difference_i(1:M-1)
        character(len=150) :: filepath
        character(len=40) :: filename
        integer :: file_code, i

        character(len=20) :: parameters(1:5)
        data parameters /"potential","j_e","j_i","n_e","n_i"/


        do i = 1, 5
            filename = ""
            file_code = 10 + i
            write(filename, "(A,A,I0)") trim(parameters(i)), "_", ind
            print *, filename

            filepath = trim(subf("\" // parameters(i))) // "\" // &
                       trim(filename) // ".txt"
            print *, filepath
            open(file_code, file = filepath, status = "new")
            write(file_code, "(A,E16.8,A)") "Time = ", time, " seconds"
        end do

        potential(0) = 0

        do k = 1, M
            potential(k) = potential(k-1) + (E_r_i(k) + E_r_i(k-1)) * h_k(k-1) / 2.0_dp
        end do

        j_e(0) = 0
        j_i(0) = 0
        do k = 1, M-1
            difference_e(k) = n_e_final(k) - n_e_i(k)
            if (difference_e(k) > 0) then
                j_e(k) = - D_e * df_dr(n_e_i(k), n_e_i(k+1), h_k(k)) &
                - k_e * n_e_i(k) * E_r_i(k)
            else
                j_e(k) = - D_e * df_dr(n_e_i(k-1), n_e_i(k), h_k(k-1)) &
                - k_e * n_e_i(k) * E_r_i(k)
            end if
        end do
        do k = 1, M - 1
            difference_i(k) = n_i_final(k) - n_i_i(k)
            if(difference_i(k) > 0) then
                j_i(k) = - D_i * df_dr(n_i_i(k), n_i_i(k+1), h_k(k)) &
                + k_i * n_i_i(k) * E_r_i(k)
            else
                 j_i(k) = - D_i * df_dr(n_i_i(k-1), n_i_i(k), h_k(k-1)) &
                 + k_i * n_i_i(k) * E_r_i(k)
            end if
        end do

        j_e(M) = -D_e * (n_e_i(M) - n_e_i(M-1)) / h_k(M-1) - k_e * n_e_i(M) * E_r_i(M)
        j_i(M) = -D_i * (n_i_i(M) - n_i_i(M-1)) / h_k(M-1) + k_i * n_i_i(M) * E_r_i(M)

        do k = 0, M
            write (11, '(F20.6,A,E20.6)') r_k(k)," ", potential(k)
            write (12, '(F20.6,A,E20.6)') r_k(k)," ", j_e(k)
            write (13, '(F20.6,A,E20.6)') r_k(k)," ", j_i(k)
            write (14, '(F20.6,A,E20.6)') r_k(k)," ", n_e_i(k)
            write (15, '(F20.6,A,E20.6)') r_k(k)," ", n_i_i(k)
        end do

        close(11)
        close(12)
        close(13)
        close(14)
        close(15)

    end subroutine

    character(len=2) function str(num)
        integer, intent(in) :: num
        if (num < 10) then
            write(str, "(A,I0)") "0", num
        else
            write(str, "(I0)") num
        end if
    end function str

    character(len=30) function data_string()
        integer :: dt(8)
        call date_and_time(values=dt)
        !character(len=2) :: day, month, minute, second
        write (data_string, "(A,A,A,A,A,I0,A,A,A,A,A,A)") "Data_", &
                            str(dt(3)), "-", str(dt(2)), "-", dt(1), "_", &
                            str(dt(5)), "-", str(dt(6)), "-", str(dt(7))
    end function data_string

    character(len=100) function subf(name)
        character(len=*), intent(in) :: name
        subf = trim(folder_name) // "\" // adjustl(name)
    end function subf

    !subroutine read_bolsig()

    !end subroutine

end program concentration
