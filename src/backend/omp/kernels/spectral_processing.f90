module m_omp_spectral
  use m_common, only: dp
  implicit none

contains

  subroutine process_spectral_000( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Post-process div U* in spectral space for all periodic BCs.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Grid size in spectral space
    integer, intent(in) :: nx_spec, ny_spec, nz_spec
    !> Offsets in the permuted pencils in spectral space
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, ix, iy, iz
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    !$omp parallel do private(div_r, div_c, ix, iy, iz, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          ! normalisation
          div_r = real(div_u(i, j, k), kind=dp)/nx/ny/nz
          div_c = aimag(div_u(i, j, k))/nx/ny/nz

          ix = i + x_sp_st
          iy = j + y_sp_st
          iz = k + z_sp_st

          ! post-process forward
          ! post-process in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) + tmp_c*az(iz)
          div_c = tmp_c*bz(iz) - tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! post-process in y
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*by(iy) + tmp_c*ay(iy)
          div_c = tmp_c*by(iy) - tmp_r*ay(iy)
          if (iy > ny/2 + 1) div_r = -div_r
          if (iy > ny/2 + 1) div_c = -div_c

          ! post-process in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
          div_c = tmp_c*bx(ix) - tmp_r*ax(ix)

          ! Solve Poisson
          tmp_r = real(waves(i, j, k), kind=dp)
          tmp_c = aimag(waves(i, j, k))
          if ((tmp_r < 1.e-16_dp) .or. (tmp_c < 1.e-16_dp)) then
            div_r = 0._dp; div_c = 0._dp
          else
            div_r = -div_r/tmp_r
            div_c = -div_c/tmp_c
          end if

          ! post-process backward
          ! post-process in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) - tmp_c*az(iz)
          div_c = -tmp_c*bz(iz) - tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! post-process in y
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*by(iy) + tmp_c*ay(iy)
          div_c = tmp_c*by(iy) - tmp_r*ay(iy)
          if (iy > ny/2 + 1) div_r = -div_r
          if (iy > ny/2 + 1) div_c = -div_c

          ! post-process in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
          div_c = -tmp_c*bx(ix) + tmp_r*ax(ix)

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_000

  subroutine process_spectral_010( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Post-process div U* in spectral space, for non-periodic BC in y-dir.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Grid size in spectral space
    integer, intent(in) :: nx_spec, ny_spec, nz_spec
    !> Offsets in the permuted pencils in spectral space
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, ix, iy, iz, iy_r
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c

    !$omp parallel do private(div_r, div_c, ix, iz, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          ix = i + x_sp_st
          iz = k + z_sp_st

          ! normalisation
          div_r = real(div_u(i, j, k), kind=dp)/nx/ny/nz
          div_c = aimag(div_u(i, j, k))/nx/ny/nz

          ! postprocess in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) + tmp_c*az(iz)
          div_c = tmp_c*bz(iz) - tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! postprocess in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
          div_c = tmp_c*bx(ix) - tmp_r*ax(ix)
          if (ix > nx/2 + 1) div_r = -div_r
          if (ix > nx/2 + 1) div_c = -div_c

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(div_r, div_c, iy, iy_r, l_r, l_c, r_r, r_c) collapse(3)
    do k = 1, nz_spec
      do j = 2, ny_spec/2 + 1
        do i = 1, nx_spec
          iy = j + y_sp_st
          iy_r = ny_spec - j + 2 + y_sp_st

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, ny_spec - j + 2, k), kind=dp)
          r_c = aimag(div_u(i, ny_spec - j + 2, k))

          ! update the entry
          div_u(i, j, k) = 0.5_dp*cmplx( & !&
            l_r*by(iy) + l_c*ay(iy) + r_r*by(iy) - r_c*ay(iy), &
            -l_r*ay(iy) + l_c*by(iy) + r_r*ay(iy) + r_c*by(iy), kind=dp &
            )
          div_u(i, ny_spec - j + 2, k) = 0.5_dp*cmplx( & !&
            r_r*by(iy_r) + r_c*ay(iy_r) + l_r*by(iy_r) - l_c*ay(iy_r), &
            -r_r*ay(iy_r) + r_c*by(iy_r) + l_r*ay(iy_r) + l_c*by(iy_r), &
           kind=dp &
           )
        end do
      end do
    end do
    !$omp end parallel do

    ! Solve Poisson
    !$omp parallel do private(div_r, div_c, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))

          tmp_r = real(waves(i, j, k), kind=dp)
          tmp_c = aimag(waves(i, j, k))
          if (abs(tmp_r) < 1.e-16_dp) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < 1.e-16_dp) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
          if (i == nx/2 + 1 .and. k == nz/2 + 1) div_u(i, j, k) = 0._dp
        end do
      end do
    end do
    !$omp end parallel do

    ! post-process backward
    !$omp parallel do private(div_r, div_c, iy, iy_r, l_r, l_c, r_r, r_c) collapse(3)
    do k = 1, nz_spec
      do j = 2, ny_spec/2 + 1
        do i = 1, nx_spec
          iy = j + y_sp_st
          iy_r = ny_spec - j + 2 + y_sp_st

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, ny_spec - j + 2, k), kind=dp)
          r_c = aimag(div_u(i, ny_spec - j + 2, k))

          ! update the entry
          div_u(i, j, k) = cmplx( & !&
            l_r*by(iy) - l_c*ay(iy) + r_r*ay(iy) + r_c*by(iy), &
            l_r*ay(iy) + l_c*by(iy) - r_r*by(iy) + r_c*ay(iy), kind=dp &
            )
          div_u(i, ny_spec - j + 2, k) = cmplx( & !&
            r_r*by(iy_r) - r_c*ay(iy_r) + l_r*ay(iy_r) + l_c*by(iy_r), &
            r_r*ay(iy_r) + r_c*by(iy_r) - l_r*by(iy_r) + l_c*ay(iy_r), &
            kind=dp &
            )
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(div_r, div_c, ix, iz, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          ix = i + x_sp_st
          iz = k + z_sp_st

          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))

          ! post-process in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) - tmp_c*az(iz)
          div_c = tmp_c*bz(iz) + tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! post-process in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) - tmp_c*ax(ix)
          div_c = tmp_c*bx(ix) + tmp_r*ax(ix)
          if (ix > nx/2 + 1) div_r = -div_r
          if (ix > nx/2 + 1) div_c = -div_c

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_010

end module m_omp_spectral
