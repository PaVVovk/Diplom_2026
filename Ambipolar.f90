subroutine ambipolar(D_amb, nu_amb)
    use variables
    implicit none
    real(dp), intent(in) :: D_amb,nu_amb
    real(dp) :: bessj0,x,xxx,ttt
    integer k
!-----------------------------------------------------
    n_e_i(0) = ne0; n_i_i(0) = ne0
    xxx = sqrt(nu_amb/D_amb)
    do k = 1,M
        x = r_k(k)*xxx
        ttt = ne0*bessj0(x)
        n_e_i(k) = ttt; n_i_i(k) = ttt
        E_r_i(k-1) = 0.0_dp
    end do
return
end
