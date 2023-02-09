module m_stencil_definitions
  use m_stencil, only: stencil
  implicit none

  enum, bind(c)
     enumerator :: ecompact6
     enumerator :: edirichlet
  end enum

  type(stencil), parameter :: compact6_order_1 = &
       & stencil( &
       & nodes = [-2, -1, 1, 2], &
       & coeffs = [- 1. / 36., -7. / 9., + 7. / 9., + 1. / 36.], &
       & upper = 1. / 3., &
       & lower = 1. / 3. &
       &)

  type(stencil), parameter :: compact6_order_2 = &
       & stencil( &
       & nodes = [-2, -1, 0, 1, 2], &
       & coeffs = [3. / 44., 12. / 11., &
         & - 2. * (12. / 11. + 3. / 44.), &
         & 12. / 11., 3. / 44.], &
       & upper = 2. / 11., &
       & lower = 2. / 11. &
       &)

  type(stencil) :: dirichlet_order_1(2)
  type(stencil) :: dirichlet_order_2(2)
  type(stencil) :: compact6(2)

contains

  dirichlet(1, 1) = stencil_type( &
       & order = 1, &
       & nodes = [0, 1, 2, 3], &
       & coeffs = [-5. / 2., 2., 0.5, 0.], &
       & lower = 0., upper = 2. &
       & )
  dirichlet(2, 1) = stencil_type( &
       & order = 1, &
       & nodes = [-1, 0, 1, 2], &
       & coeffs = [-3. / 4., 0., 3. / 4., 0.], &
       & lower = 1. / 4., upper = 1. / 4. &
       & )

  dirichlet(1, 2) = stencil_type( &
       & order = 2, &
       & nodes = [0, 1, 2, 3], &
       & coeffs = [13., -27., 15., -1.], &
       & lower = 0., upper = 11. &
       & )
  dirichlet(2, 2) = stencil_type( &
       & order = 2, &
       & nodes = [-1, 0, 1, 2], &
       & coeffs = [6. / 5., -12. / 5., 6. / 5., 0.], &
       & lower = 1. / 10., upper = 1. / 10. &
       & )

  pure function get_stencil(key, order, side) result(s)
    integer, intent(in) :: key, order
    character(*), optional :: side
    type(stencil) :: s
    if (.not. any([1, 2] == order)) then
       error stop "order must be 1 or 2"
    end if
    if present(side) then
       if (.not. any(["left", "right"] == side)) then
          error stop "'side' argument must be 'left' or 'right'"
       end if
    end if

    if (key == ecompact6) then
       s = compact6(order)
    else
       error stop "Unknown key"
    end if
  end function get_stencil

  pure function get_boundary_stencils(key, order) result(s)
    integer, intent(in) :: key, order
    type(stencil) :: s(2)
    if (.not. any([1, 2] == order)) then
       error stop "order must be 1 or 2"
    end if
    if (key == edirichlet) then
       s = dirichlet(:, order)
    else
       error stop "Unknown key"
    end if
    if (side == "right") then
       s = s%flip()
    end if
  end function get_stencil
end module m_stencil_definitions