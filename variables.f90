module variables
    public
    integer, parameter :: dp = kind(1.0d0)
    integer :: M, k_maxit
    real(dp), parameter :: sigma = 0.5_dp
    real(dp), parameter :: sigm1 = 1 - sigma
    real(dp), allocatable :: r_k(:), r_half_k(:), h_k(:), h_half_k(:)
    real(dp) :: l_e, l_i, D_e, D_i, k_e, k_i, nu_ion, p, N, beta_ei, N_L, ne0
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: k_b = 1.380649E-23_dp
    real(dp), parameter :: e = 1.602176634E-19_dp
    real(dp), parameter :: const = 2 * pi * e
    real(dp), parameter :: gamma_e = 0.7104_dp
    real(dp), parameter :: gamma_i = 0.7104_dp
    real(dp), parameter :: T_gas = 300.0_dp
    real(dp), allocatable :: n_e_i(:), n_i_i(:), E_r_i(:)
    real(dp), allocatable :: n_e_final(:), n_i_final(:), E_r_final(:)
    real(dp), allocatable :: a_km(:),b_km(:),u_km(:),w_km(:)
    character ch    !fav
    character*6 ch_nt0
    character(len=50) :: folder_name

end module variables
