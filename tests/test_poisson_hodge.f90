!!! tests/test_poisson_hodge.f90
!!
!! Tests the Poisson solver by enforcement of the Hodge decomposition.

program test_hodge

  use mpi
  
  use m_common, only: dp, pi

  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t

  implicit none

  class(mesh_t), allocatable :: mesh
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend

  integer, dimension(3) :: dims_global
  real(dp), dimension(3) :: L_global
  integer, dimension(3) :: nproc_dir
  character(len=20), dimension(2) :: BC_x, BC_y, BC_z

  logical :: use_2decomp

  integer :: nproc, nrank
  integer :: ierr

  ! Initialise variables and arrays
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
  use_2decomp = .false.
#else
  use_2decomp = .true.
#endif

  ! Global number of cells in each direction
  dims_global = [64, 32, 128]

  ! Global domain dimensions
  L_global = [2*pi, 2*pi, 2*pi]

  ! Domain decomposition in each direction
  nproc_dir = [1, 1, nproc]

  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z, use_2decomp=use_2decomp)

  allocator => initialise_allocator(mesh)
  backend => initialise_backend(allocator)

contains

  !! Initialises the generic allocator pointer
  function initialise_allocator(mesh) result(allocator)
#ifdef CUDA
    use m_cuda_allocator, only: cuda_allocator_t
    use m_cuda_common, only: SZ
#else
    use m_omp_common, only: SZ
#endif
    class(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), pointer :: allocator
#ifdef CUDA
    type(cuda_allocator_t), target :: cuda_allocator
    
    cuda_allocator = cuda_allocator_t(mesh, SZ)
    allocator => cuda_allocator
    print *, 'CUDA allocator instantiated'
#else
    type(allocator_t), target :: omp_allocator
    
    omp_allocator = allocator_t(mesh, SZ)
    allocator => omp_allocator
#endif
  end function initialise_allocator

  !! Initialises the generic backend pointer
  function initialise_backend(allocator) result(backend)
#ifdef CUDA
    use m_cuda_backend, only: cuda_backend_t
#else
    use m_omp_backend, only: omp_backend_t
#endif
    class(allocator_t), pointer, intent(in) :: allocator
    class(base_backend_t), pointer :: backend
#ifdef CUDA
    type(cuda_backend_t), target :: cuda_backend

    cuda_backend = cuda_backend_t(mesh, allocator)
    backend => cuda_backend
    print *, 'CUDA backend instantiated'
#else
    type(omp_backend_t), target :: omp_backend

    omp_backend = omp_backend_t(mesh, allocator)
    backend => omp_backend
#endif
  end function initialise_backend
  
end program test_hodge
