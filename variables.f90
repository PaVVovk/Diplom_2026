module variables
	implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer :: M,k_maxit,iter
    integer :: mode, bc_mode
    real(dp) :: sigma, sigm1, const,tol,T_e, ne0
    real(dp), allocatable :: r_k(:), r_half_k(:), h_k(:), h_half_k(:)
    real(dp) :: EN,D_e,D_i,k_e,k_i,nu_ion,p,N,beta_ei,N_L,E_mean
    real(dp) :: I, U_D, R, L
    real(dp), parameter :: k_b = 1.380649E-23_dp, pi = 3.14159265_dp
    real(dp), parameter :: e = 1.602176634E-19_dp, epsilon0 = 8.85E-12_dp, m_e = 9.1e-31_dp, mu = 4.0e-2_dp
    real(dp), parameter :: gamma_e = 0.7104_dp 
    real(dp), parameter :: gamma_i = 0.7104_dp
    real(dp), parameter :: T_gas = 300.0_dp
    real(dp), allocatable :: n_e_i(:), n_i_i(:), E_r_i(:)
    real(dp), allocatable :: n_e_final(:), n_i_final(:), E_r_final(:)
    real(dp), allocatable :: a_km(:),b_km(:),w_km(:)
    character ch    !fav
    character*6 ch_nt0
end module variables
