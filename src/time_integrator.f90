module m_time_integrator
  use m_allocator, only: allocator_t, field_t, flist_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X

  implicit none

  private adams_bashforth

  type :: time_intg_t
    integer :: istep, order, nvars, nolds
    real(dp) :: coeffs(4, 4)
    type(flist_t), allocatable :: olds(:, :)
    type(flist_t), allocatable :: curr(:)
    type(flist_t), allocatable :: deriv(:)
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
  contains
    procedure :: step
    procedure :: adams_bashforth
  end type time_intg_t

  interface time_intg_t
    module procedure init
  end interface time_intg_t

contains

  function init(backend, allocator, order, nvars)
    implicit none

    type(time_intg_t) :: init
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    integer, intent(in), optional :: order
    integer, intent(in), optional :: nvars

    integer :: i, j

    ! initialize Adams-Bashforth coefficients
    init%coeffs(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    init%coeffs(:, 2) = [1.5_dp, -0.5_dp, 0.0_dp, 0.0_dp]
    init%coeffs(:, 3) = &
      [23._dp/12._dp, -4._dp/3._dp, 5._dp/12._dp, 0.0_dp]
    init%coeffs(:, 4) = &
      [55._dp/24._dp, -59._dp/24._dp, 37._dp/24._dp, -3._dp/8._dp]

    ! set variables
    init%backend => backend
    init%allocator => allocator

    if (present(order)) then
      init%order = order
    else
      init%order = 1
    end if

    if (present(nvars)) then
      init%nvars = nvars
    else
      init%nvars = 3
    end if

    init%istep = 1
    init%nolds = init%order - 1

    ! allocate memory
    allocate (init%olds(init%nvars, init%nolds))
    allocate (init%curr(init%nvars))
    allocate (init%deriv(init%nvars))

    ! Request all the storage for old timesteps
    do i = 1, init%nvars
      do j = 1, init%nolds
        init%olds(i, j)%ptr => allocator%get_block(DIR_X)
      end do
    end do

  end function init

  subroutine step(self, u, v, w, du, dv, dw, dt)
    implicit none

    class(time_intg_t), intent(inout) :: self
    class(field_t), target, intent(inout) :: u, v, w
    class(field_t), target, intent(in) :: du, dv, dw

    real(dp), intent(in) :: dt

    ! assign pointer to variables
    self%curr(1)%ptr => u
    self%curr(2)%ptr => v
    self%curr(3)%ptr => w

    ! assign pointer to variables
    self%deriv(1)%ptr => du
    self%deriv(2)%ptr => dv
    self%deriv(3)%ptr => dw

    call self%adams_bashforth(dt)

    ! increment step counter
    self%istep = self%istep + 1
  end subroutine step

  subroutine adams_bashforth(self, dt)
    class(time_intg_t), intent(inout) :: self
    real(dp), intent(in) :: dt

    integer :: i, j
    integer :: order

    order = min(self%istep, self%order)
    do i = 1, self%nvars
      ! update solution
      call self%backend%vecadd(self%coeffs(1, order)*dt, &
                               self%deriv(i)%ptr, &
                               1._dp, self%curr(i)%ptr)
      do j = 2, order
        call self%backend%vecadd(self%coeffs(j, order)*dt, &
                                 self%olds(i, j - 1)%ptr, &
                                 1._dp, self%curr(i)%ptr)
      end do

      ! rotate pointers
      if (order < self%order) then
        ! for startup
        if (self%istep > 1) then
          call rotate(self%olds(i, :), order)
        end if
      else
        ! after startup
        if (self%order > 2) then
          call rotate(self%olds(i, :), order - 1)
        end if
      end if

      ! update olds(1) with new derivative
      if (self%order > 1) then
        call self%backend%vecadd(1.0_dp, self%deriv(i)%ptr, 0._dp, &
                                 self%olds(i, 1)%ptr)
      end if
    end do

  end subroutine adams_bashforth

  subroutine rotate(sol, n)
    type(flist_t), intent(inout) :: sol(:)
    integer, intent(in) :: n

    integer :: i
    class(field_t), pointer :: ptr

    ! rotate pointer
    ptr => sol(n)%ptr
    do i = n, 2, -1
      sol(i)%ptr => sol(i - 1)%ptr
    end do
    sol(1)%ptr => ptr

  end subroutine rotate
end module m_time_integrator
