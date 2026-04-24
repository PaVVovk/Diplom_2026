    !Итерации
    subroutine solve_iterations(tau,repeat_flag)
    use variables
    implicit none
    ! Входные параметры
    logical, intent(inout) :: repeat_flag
    real(dp), intent(inout) :: tau
    ! Выходные параметры - передаются через module
    !Локальные переменные
    integer :: k, j
    real(dp) :: error_ne, error_ni!, error_E
    real(dp) :: n_e_m(0:M), n_i_m(0:M), E_r_m(0:M)
    real(dp) :: n_e_im(0:M), n_i_im(0:M), E_r_im(0:M)
    real(dp) :: n_e_m1(0:M), n_i_m1(0:M), E_r_m1(0:M)
    real(dp) :: max_error, tol, ne_max,ni_max,E_max

    n_e_m = n_e_i
    n_i_m = n_i_i
    E_r_m = E_r_i

    tol = 1.0E-4

    do k = 1, k_maxit
        max_error = 0.0_dp
        do j = 0, M
            n_e_im(j) = sigma * n_e_m(j) + sigm1 * n_e_i(j)
            n_i_im(j) = sigma * n_i_m(j) + sigm1 * n_i_i(j)
            E_r_im(j) = sigma * E_r_m(j) + sigm1 * E_r_i(j)
        end do
        call solve_electron_continuity(tau,n_e_m1,n_i_im,E_r_im)
!n_e_im(j) = sigma * n_e_m(j) + sigm1 * n_e_i(j)
        call solve_ion_continuity(tau,n_i_m1, n_e_im, E_r_im)
        call solve_poisson_trapezoidal(n_e_m1, n_i_m1, E_r_m1)
        ne_max = maxval(abs(n_e_m1)); ni_max = maxval(abs(n_i_m1)); E_max = maxval(abs(E_r_m1))
!write (*,*)
         do j = 0, M
            error_ne = abs(n_e_m1(j) - n_e_m(j)) / ne_max
            error_ni = abs(n_i_m1(j) - n_i_m(j)) / ni_max
!           error_E =  abs(E_r_m1(j) - E_r_m(j)) / E_max
            max_error = max(error_ne, error_ni, max_error)!, error_E)!)
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

        if (max_error < tol) then
            n_e_final = n_e_m1
            n_i_final = n_i_m1
            E_r_final = E_r_m1
            if (k < k_maxit / 2) tau = tau*2.0_dp
            repeat_flag = .false.
            exit
        end if

        do j = 0, M
            n_e_m(j) = n_e_m1(j)
            n_i_m(j) = n_i_m1(j)
            E_r_m(j) = E_r_m1(j)
        end do
    end do
        if (k == k_maxit) then
            tau = tau/2.0_dp
            repeat_flag = .true.
        end if
    end subroutine solve_iterations
