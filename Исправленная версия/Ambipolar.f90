subroutine ambipolar(D_amb,R,ne0)
    use variables
    implicit none
    real(dp), intent(in) :: D_amb,R,ne0
    real (8) bessj0,x,lam_01,xxx,ttt
    integer k
!-----------------------------------------------------
    n_e_i(0) = ne0; n_i_i(0) = ne0; E_r_i(0) = 0.0_dp
    lam_01 = 2.4048d0
    xxx = sqrt(nu_ion/D_amb)
    write (33,'(1p,E14.6,3(a,E14.6))') ne0,ch,xxx,ch,nu_ion,ch,D_amb
    do k = 1,M
        x = r_k(k)*xxx
        ttt = ne0*bessj0(x)
        n_e_i(k) = ttt; n_i_i(k) = ttt
        E_r_i(k) = 0.0_dp
        write (33,'(1p,E14.6,4(a,E14.6))') r_k(k),ch,n_e_i(k),ch,n_i_i(k),ch,x,ch,bessj0(x)
    end do
    close (33)
return
end
