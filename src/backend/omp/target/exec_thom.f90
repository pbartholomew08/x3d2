module m_omptgt_exec_thom

  use m_common, only: dp
  use m_tdsops, only: tdsops_t

  use m_omptgt_kernels_thom, only: der_univ_thom, der_univ_thom_per

  implicit none

  private
  public :: exec_thom_tds_compact

contains

  subroutine exec_thom_tds_compact(du, u, tdsops, n_groups)

    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    type(tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: n_groups

    integer :: k

    if (tdsops%periodic) then
      call exec_thom_tds_per_(du, u, &
                              tdsops%n_tds, &
                              tdsops%coeffs, tdsops%alpha, &
                              tdsops%thom_f, tdsops%thom_s, tdsops%thom_w, tdsops%thom_p, &
                              tdsops%stretch, &
                              n_groups &
                              )
    else
!!!       call exec_thom_tds_compact_(du, u, &
!!!                                   tdsops%n_tds, tdsops%n_rhs, &
!!!                                   tdsops%coeffs_s, tdsops%coeffs_e, tdsops%coeffs, &
!!!                                   tdsops%thom_f, tdsops%thom_s, tdsops%thom_w, &
!!!                                   tdsops%stretch, &
!!!                                   n_groups &
!!!                                   )
       call der_univ_thom(du, u, &
                          tdsops%n_tds, tdsops%n_rhs, &
                          tdsops%coeffs_s, tdsops%coeffs_e, tdsops%coeffs, &
                          tdsops%thom_f, tdsops%thom_s, tdsops%thom_w, &
                          tdsops%stretch, &
                          n_groups &
                          )
                          
    end if

  end subroutine exec_thom_tds_compact

  subroutine exec_thom_tds_compact_(du, u, &
                                    n_tds, n_rhs, &
                                    coeffs_s, coeffs_e, coeffs, &
                                    thom_f, thom_s, thom_w, &
                                    stretch, &
                                    n_groups &
                                    )
    
    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    integer, intent(in) :: n_tds, n_rhs
    real(dp), intent(in), dimension(:, :) :: coeffs_s, coeffs_e ! start/end
    real(dp), intent(in), dimension(:) :: coeffs
    real(dp), intent(in), dimension(:) :: thom_f, thom_s, thom_w
    real(dp), intent(in), dimension(:) :: stretch
    integer, intent(in) :: n_groups

!!!    !$omp target data map(to:u, coeffs_s, coeffs_e, coeffs, thom_f, thom_s, thom_w, stretch) map(from:du)
!!!    !$omp target teams distribute
!!!    do k = 1, n_groups
!!!      call der_univ_thom( &
!!!        du(:, :, k), u(:, :, k), n_tds, n_rhs, &
!!!        coeffs_s, coeffs_e, coeffs, &
!!!        thom_f, thom_s, thom_w, &
!!!        stretch &
!!!        )
!!!    end do
!!!    !$omp end target teams distribute
!!!    !$omp end target data

    du = 42.0_dp
    print *, "+++:", minval(du), maxval(du)
    !$omp target data map(to:u, coeffs_s, coeffs_e, coeffs, thom_f, thom_s, thom_w, stretch) map(from:du)
    call der_univ_thom( &
      du, u, n_tds, n_rhs, &
      coeffs_s, coeffs_e, coeffs, &
      thom_f, thom_s, thom_w, &
      stretch, n_groups &
      )
    !$omp end target data
    print *, "+++:", minval(du), maxval(du)

  end subroutine

  subroutine exec_thom_tds_per_(du, u, &
                                n_tds, &
                                coeffs, alpha, &
                                thom_f, thom_s, thom_w, thom_p, &
                                stretch, &
                                n_groups &
                                )

    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u
    integer, intent(in) :: n_tds
    real(dp), intent(in), dimension(:) :: coeffs
    real(dp), intent(in) :: alpha
    real(dp), intent(in), dimension(:) :: thom_f, thom_s, thom_w, thom_p
    real(dp), intent(in), dimension(:) :: stretch
    integer, intent(in) :: n_groups

    integer :: k

    !$omp target data map(to:u, coeffs, thom_f, thom_s, thom_w, thom_p, stretch) map(from:du)
    !$omp target teams distribute
    do k = 1, n_groups
      call der_univ_thom_per( &
        du(:, :, k), u(:, :, k), n_tds, coeffs, alpha, &
        thom_f, thom_s, thom_w, thom_p, &
        stretch &
        )
    end do
    !$omp end target teams distribute
    !$omp end target data

  end subroutine

end module m_omptgt_exec_thom
