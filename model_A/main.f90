!-------------------------------------------------------------------------------
!   MPI Finite Difference Phase Field Code of Allen-Cahn Equation.
!
!   Ported from the serial reference:
!   https://github.com/Shahid718/Programming-Phase-field-in-Fortran/blob/main/model_A/allen_cahn/fd_ac.f90
!   Domain decomposition is along y. NOTE: the free-energy derivative (dfdphi)
!   uses the correct double-well form 2*A*phi*(1-phi)*(1-2*phi), i.e. the
!   derivative of f(phi) = A*phi^2*(1-phi)^2 -- matching the Cahn-Hilliard
!   code's potential. The original repo's dfdphi has an extra (1-phi) factor,
!   which does not match that potential's derivative, so it was not carried over.
!
!   To compile and run:
!             mpif90 -std=gnu main.f90 -o main
!             mpirun -np 8 ./main
!-------------------------------------------------------------------------------

program allen_cahn_mpi
    use mpi
    use iso_fortran_env, only: real64, dp=>real64      ! fortran intrinsic module 
    implicit none

    !=============================================================================
    !                    PARAMETERS  (match the serial reference)
    !=============================================================================
    ! simulation cell parameters
    integer, parameter :: Nx = 128
    integer, parameter :: Ny = 128
    real(dp), parameter :: dx = 2.0d0
    real(dp), parameter :: dy = 2.0d0

    ! time integration parameters
    integer, parameter :: nsteps = 2000       
    integer, parameter :: nprint = 1000                ! changed to 1000
    real(dp), parameter :: dt = 0.01d0

    ! material specific parameters
    real(dp), parameter :: phi_0 = 0.5d0
    real(dp), parameter :: mobility = 1.0d0
    real(dp), parameter :: grad_coef = 1.0d0

    ! microstructure parameters
    real(dp), parameter :: noise = 0.02d0
    real(dp), parameter :: A = 1.0d0

    ! MPI variables
    integer :: rank, nprocs, ierr
    integer :: status(MPI_STATUS_SIZE)

    ! MPI domain decomposition (along y)
    integer :: local_Ny, start_y, end_y
    integer :: down, up

    ! Arrays (allocatable for local domain with halo rows 0 and local_Ny+1)
    real(dp), allocatable, dimension(:,:) :: phi, phi_new
    real(dp), allocatable, dimension(:,:) :: dfdphi, lap_phi, dummy_phi

    ! Gather variables
    real(dp), allocatable :: full_phi(:, :)
    real(dp), allocatable :: recv_buf(:, :)
    real(dp), allocatable :: send_buf(:, :)
    integer :: p_local_Ny, p_start, p

    ! Local variables
    integer :: i, j, tsteps, im, ip, jp, jm
    real(dp) :: r
    real(dp) :: mpi_start_time, mpi_end_time

    !=============================================================================
    ! Initialize MPI
    !=============================================================================
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    !=============================================================================
    ! Domain decomposition along y-direction
    !=============================================================================
    local_Ny = Ny / nprocs
    start_y  = rank * local_Ny + 1
    end_y    = merge(Ny, start_y + local_Ny - 1, rank == nprocs - 1)
    local_Ny = end_y - start_y + 1

    !=============================================================================
    ! Allocate local arrays with halo rows 0 and local_Ny+1
    !=============================================================================
    allocate(phi(1:Nx,0:local_Ny+1))
    allocate(phi_new(1:Nx,0:local_Ny+1))
    allocate(dfdphi(1:Nx,0:local_Ny+1))
    allocate(lap_phi(1:Nx,0:local_Ny+1))
    allocate(dummy_phi(1:Nx,0:local_Ny+1))

    phi = 0.0d0
    phi_new = 0.0d0
    dfdphi = 0.0d0
    lap_phi = 0.0d0
    dummy_phi = 0.0d0

    !=============================================================================
    ! Initial microstructure (only local domain; matches serial reference)
    !=============================================================================
    do j = 1, local_Ny
        do i = 1, Nx
            call random_number(r)
            phi(i, j) = phi_0 + noise * (0.5d0 - r)
        end do
    end do
 
    !=============================================================================
    ! Save Data from Each Rank
    !============================================================================= 
    ! if (rank == 0) then
        ! open(unit=10, file='ac_rank0.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: ac_rank0.dat'
    ! end if

    ! if (rank == 1) then
        ! open(unit=10, file='ac_rank1.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: ac_rank1.dat'
    ! end if

    ! if (rank == 2) then
        ! open(unit=10, file='ac_rank2.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: ac_rank2.dat'
    ! end if
    
    ! if (rank == 3) then
        ! open(unit=10, file='ac_rank3.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: ac_rank3.dat'
    ! end if

    !=============================================================================
    ! Start timing
    !=============================================================================
    mpi_start_time = MPI_Wtime()

    !=============================================================================
    ! Main evolution loop
    !=============================================================================
    TIME_LOOP: do tsteps = 1, nsteps

        down = rank - 1
        up   = rank + 1
        if (rank == 0)          down = nprocs - 1
        if (rank == nprocs - 1) up   = 0

        ! Exchange halo rows (periodic in y across ranks)
        call MPI_Sendrecv( &
            phi(1,1), Nx, MPI_DOUBLE_PRECISION, down, 0, &
            phi(1,local_Ny+1), Nx, MPI_DOUBLE_PRECISION, up, 0, &
            MPI_COMM_WORLD, status, ierr)

        call MPI_Sendrecv( &
            phi(1,local_Ny), Nx, MPI_DOUBLE_PRECISION, up, 1, &
            phi(1,0), Nx, MPI_DOUBLE_PRECISION, down, 1, &
            MPI_COMM_WORLD, status, ierr)

        !=========================================================================
        ! Free energy derivative, Laplacian, and time integration
        ! (formula matches the serial reference exactly)
        !=========================================================================
        do j = 1, local_Ny
            do i = 1, Nx

                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1

                ! Periodic boundary in x
                if (ip > Nx) ip = 1
                if (im < 1)  im = Nx

                dfdphi(i,j) = A * ( 2.0d0 * phi(i,j) * (1.0d0 - phi(i,j)) * &
                                    (1.0d0 - phi(i,j)) - 2.0d0 * phi(i,j) * &
                                    phi(i,j) * (1.0d0 - phi(i,j)) )

                lap_phi(i,j) = ( phi(ip,j) + phi(im,j) + phi(i,jm) + phi(i,jp) - &
                                 4.0d0 * phi(i,j) ) / (dx * dy)

                dummy_phi(i,j) = dfdphi(i,j) - grad_coef * lap_phi(i,j)

                phi_new(i,j) = phi(i,j) - dt * mobility * dummy_phi(i,j)

                if (phi_new(i,j) >= 0.99999d0) phi_new(i,j) = 0.99999d0
                if (phi_new(i,j) <  0.00001d0) phi_new(i,j) = 0.00001d0

            end do
        end do

        phi(1:Nx,1:local_Ny) = phi_new(1:Nx,1:local_Ny)

        if (mod(tsteps, nprint) == 0 .and. rank == 0) &
            write(*, '(A, I0)') 'Done steps = ', tsteps

    end do TIME_LOOP

    !=============================================================================
    ! End timing
    !=============================================================================
    mpi_end_time = MPI_Wtime()

    if (rank == 0) then
        write(*, '(A)') '---------------------------------'
        write(*, '(A, F10.3, A)') '  MPI Time    = ', mpi_end_time - mpi_start_time, ' seconds.'
    end if

    !=============================================================================
    ! Save Data from Each Rank
    !=============================================================================
    
    ! if (rank == 0) then
        ! open(unit=10, file='rank0.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: rank0.dat'
    ! end if

    ! if (rank == 1) then
        ! open(unit=10, file='rank1.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: rank1.dat'
    ! end if

    ! if (rank == 2) then
        ! open(unit=10, file='rank2.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: rank2.dat'
    ! end if
    
    ! if (rank == 3) then
        ! open(unit=10, file='rank3.dat', status='replace', action='write')
        ! do j = 1, local_Ny
            ! write(10, *) (phi(i, j), i = 1, Nx)
        ! end do
        ! close(10)
        ! write(*, '(A)') 'Results written to: rank3.dat'
    ! end if

    !=============================================================================
    ! Gather results to main and write ac.dat (full grid)
    !=============================================================================
    if (rank == 0) then
        allocate(full_phi(1:Nx, 1:Ny))

        do j = 1, local_Ny
            do i = 1, Nx
                full_phi(i, j) = phi(i, j)
            end do
        end do

        do p = 1, nprocs - 1
            p_local_Ny = Ny / nprocs
            p_start = p * p_local_Ny + 1
            if (p == nprocs - 1) p_local_Ny = Ny - p_start + 1

            allocate(recv_buf(1:Nx, 1:p_local_Ny))
            call MPI_Recv(recv_buf, p_local_Ny * Nx, MPI_DOUBLE_PRECISION, p, 0, &
                          MPI_COMM_WORLD, status, ierr)

            do j = 1, p_local_Ny
                do i = 1, Nx
                    full_phi(i, p_start + j - 1) = recv_buf(i, j)
                end do
            end do

            deallocate(recv_buf)
        end do

        open(unit=10, file='ac.dat', status='replace', action='write')
        do i = 1, Nx
            write(10, *) (full_phi(i, j), j = 1, Ny)
        end do
        close(10)
        write(*, '(A)') 'Results written to: ac.dat'

        deallocate(full_phi)
    else
        allocate(send_buf(1:Nx, 1:local_Ny))

        do j = 1, local_Ny
            do i = 1, Nx
                send_buf(i, j) = phi(i, j)
            end do
        end do

        call MPI_Send(send_buf, local_Ny * Nx, MPI_DOUBLE_PRECISION, 0, 0, &
                      MPI_COMM_WORLD, ierr)

        deallocate(send_buf)
    end if

    deallocate(phi, phi_new, dfdphi, lap_phi, dummy_phi)

    call MPI_Finalize(ierr)

end program allen_cahn_mpi