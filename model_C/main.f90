!-------------------------------------------------------------------------------
!   MPI Finite Difference Phase Field Code for Solving Model C
!   (Cahn-Hilliard + Allen-Cahn coupled equations)
!
!   To compile and run:
!             mpif90 -std=gnu main.f90 -o main
!             mpirun -np 8 ./main
!-------------------------------------------------------------------------------

program fd_ch_ac_mpi
    use mpi
    use iso_fortran_env, only: real64, dp=>real64
    implicit none

    !=============================================================================
    !                    PARAMETERS
    !=============================================================================
    ! simulation cell parameters
    integer, parameter :: Nx = 128
    integer, parameter :: Ny = 128
    real(dp), parameter :: dx = 1.0
    real(dp), parameter :: dy = 1.0

    ! time integration parameters
    integer, parameter :: nsteps = 10000
    integer, parameter :: nprint = 2000
    real(dp), parameter :: dt = 0.03d0

    ! material specific parameters
    real(dp), parameter :: A = 1.0d0
    real(dp), parameter :: B = 1.0d0
    real(dp), parameter :: D = 1.0d0
    real(dp), parameter :: mobility_con = 0.5d0
    real(dp), parameter :: mobility_phi = 0.5d0
    real(dp), parameter :: grad_coef_con = 1.5d0
    real(dp), parameter :: grad_coef_phi = 1.5d0
    real(dp), parameter :: radius = 10.0d0

    ! MPI variables
    integer :: rank, size, ierr
    integer :: status(MPI_STATUS_SIZE)

    ! MPI Domain decomposition
    integer :: local_ny, start_y, end_y
    integer :: down, up

    ! Arrays (allocatable for local domain with halo)
    real(dp), allocatable, dimension(:,:) :: con, phi
    real(dp), allocatable, dimension(:,:) :: dfdcon, dfdphi
    real(dp), allocatable, dimension(:,:) :: lap_con, lap_phi
    real(dp), allocatable, dimension(:,:) :: dummy_con, lap_dummy, phi_dummy

    ! MPI Variables for gathering
    real(dp), allocatable :: full_con(:, :), full_phi(:, :)
    real(dp), allocatable :: recv_buf(:, :)
    real(dp), allocatable :: send_buf(:, :)
    integer :: p_local_Ny, p_start, p

    ! Local variables
    integer :: i, j, tsteps, im, ip, jp, jm
    real(dp) :: mpi_start_time, mpi_end_time

    !=============================================================================
    ! Initialize MPI
    !=============================================================================
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    !=============================================================================
    ! Domain decomposition along y-direction
    !=============================================================================  
    local_ny = Ny / size
    start_y  = rank * local_ny + 1
    end_y    = merge(Ny, start_y + local_ny - 1, rank == size - 1)
    local_ny = end_y - start_y + 1

    !=============================================================================      
    ! Allocate arrays with 1-based indexing for both dimensions
    !=============================================================================  
    allocate(con(1:Nx,0:local_ny+1))
    allocate(phi(1:Nx,0:local_ny+1))
    allocate(dfdcon(1:Nx,0:local_ny+1))
    allocate(dfdphi(1:Nx,0:local_ny+1))
    allocate(lap_con(1:Nx,0:local_ny+1))
    allocate(lap_phi(1:Nx,0:local_ny+1))
    allocate(dummy_con(1:Nx,0:local_ny+1))
    allocate(lap_dummy(1:Nx,0:local_ny+1))
    allocate(phi_dummy(1:Nx,0:local_ny+1))

    !=============================================================================  
    ! Initialize arrays to zero
    !=============================================================================
    con = 0.0d0
    phi = 0.0d0
    dfdcon = 0.0d0
    dfdphi = 0.0d0
    lap_con = 0.0d0
    lap_phi = 0.0d0
    dummy_con = 0.0d0
    lap_dummy = 0.0d0
    phi_dummy = 0.0d0

    !=============================================================================    
    ! Initialize microstructure (only local domain)
    !=============================================================================
    do j = 1, local_ny
        do i = 1, Nx
            if ((i - Nx/2)*(i - Nx/2) + (j + start_y - 1 - Ny/2)* &
                (j + start_y - 1 - Ny/2) < radius**2) then
                con(i, j) = 1.0d0
                phi(i, j) = 1.0d0
            else
                con(i, j) = 0.02d0
                phi(i, j) = 0.0d0
            endif
        end do
    end do

    !=============================================================================    
    ! Initialize halos (will be updated via MPI)
    !=============================================================================  
    do i = 1, Nx
        con(i,0) = 0.0d0
        con(i,local_ny+1) = 0.0d0
        phi(i,0) = 0.0d0
        phi(i,local_ny+1) = 0.0d0
    end do

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

        if (rank == 0)      down = size - 1
        if (rank == size-1) up   = 0

        ! Exchange halos for con
        call MPI_Sendrecv( &
            con(1,1), Nx, MPI_DOUBLE_PRECISION, down, 0, &
            con(1,local_ny+1), Nx, MPI_DOUBLE_PRECISION, up, 0, &
            MPI_COMM_WORLD, status, ierr)

        call MPI_Sendrecv( &
            con(1,local_ny), Nx, MPI_DOUBLE_PRECISION, up, 1, &
            con(1,0), Nx, MPI_DOUBLE_PRECISION, down, 1, &
            MPI_COMM_WORLD, status, ierr)

        ! Exchange halos for phi
        call MPI_Sendrecv( &
            phi(1,1), Nx, MPI_DOUBLE_PRECISION, down, 0, &
            phi(1,local_ny+1), Nx, MPI_DOUBLE_PRECISION, up, 0, &
            MPI_COMM_WORLD, status, ierr)

        call MPI_Sendrecv( &
            phi(1,local_ny), Nx, MPI_DOUBLE_PRECISION, up, 1, &
            phi(1,0), Nx, MPI_DOUBLE_PRECISION, down, 1, &
            MPI_COMM_WORLD, status, ierr)

        !=========================================================================
        ! Derivative wrt concentration and phi, and first laplacian for con
        !=========================================================================
        do j = 1, local_ny
            do i = 1, Nx

                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1

                ! Periodic boundary in x
                if (ip > Nx) ip = 1
                if (im < 1)  im = Nx

                ! derivative wrt concentration and phi
                dfdcon(i,j) = 2*A*con(i,j)*(1-( phi(i,j)**3*(10 - 15*phi(i,j) + &
                              6*phi(i,j)**2 ))) - 2*B*(1 - con(i,j))* &
                              ( phi(i,j)**3*(10 - 15*phi(i,j) + 6*phi(i,j)**2 ) )

                dfdphi(i,j) = -A*con(i,j)*con(i,j)*(3*phi(i,j)**2*(10 - &
                              15*phi(i,j) + 6*phi(i,j)**2 ) + phi(i,j)**3* &
                              (12*phi(i,j) - 15)) + 2*B*(1 - con(i,j))* &
                              (1 - con(i,j))*(3*phi(i,j)**2*(10 - 15*phi(i,j) + &
                              6*phi(i,j)**2 ) + phi(i,j)**3*(12*phi(i,j) - 15)) + &
                              2*D*phi(i,j)*(1 - phi(i,j))*(1 - 2*phi(i,j))

                ! concentration laplacian
                lap_con(i,j) = ( con(ip,j) + con(im,j) + con(i,jp) + con(i,jm) - &
                                 4.0d0*con(i,j) ) / (dx*dy)
                dummy_con(i,j) = dfdcon(i,j) - grad_coef_con*lap_con(i,j)

            end do
        end do

        ! Exchange halos for dummy_con
        call MPI_Sendrecv( &
            dummy_con(1,1), Nx, MPI_DOUBLE_PRECISION, down, 0, &
            dummy_con(1,local_ny+1), Nx, MPI_DOUBLE_PRECISION, up, 0, &
            MPI_COMM_WORLD, status, ierr)

        call MPI_Sendrecv( &
            dummy_con(1,local_ny), Nx, MPI_DOUBLE_PRECISION, up, 1, &
            dummy_con(1,0), Nx, MPI_DOUBLE_PRECISION, down, 1, &
            MPI_COMM_WORLD, status, ierr)

        !========================================================================= 
        ! Second laplacian for con, laplacian for phi, and time integration
        !=========================================================================
        do j = 1, local_ny
            do i = 1, Nx

                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1

                ! Periodic boundary in x
                if (ip > Nx) ip = 1
                if (im < 1)  im = Nx

                ! Laplacian of dummy
                lap_dummy(i,j) = ( dummy_con(ip,j) + dummy_con(im,j) + &
                                   dummy_con(i,jp) + dummy_con(i,jm) - &
                                   4.0d0*dummy_con(i,j) ) / (dx*dy)

                ! phi laplacian
                lap_phi(i,j) = ( phi(ip,j) + phi(im,j) + phi(i,jp) + phi(i,jm) - &
                                 4.0d0*phi(i,j) ) / (dx*dy)
                phi_dummy(i,j) = dfdphi(i,j) - grad_coef_phi*lap_phi(i,j)

                ! time integration
                con(i,j) = con(i,j) + dt*mobility_con*lap_dummy(i,j)
                phi(i,j) = phi(i,j) - dt*mobility_phi*phi_dummy(i,j)

            end do
        end do

        !=========================================================================
        ! For small deviations
        !=========================================================================
        do j = 1, local_ny
            do i = 1, Nx
                if (phi(i,j) >= 0.99999d0) phi(i,j) = 0.99999d0
                if (phi(i,j) < 0.00001d0)  phi(i,j) = 0.00001d0
            end do
        end do

        ! Print progress (only main)
        if (mod(tsteps, nprint) == 0 .and. rank == 0) write(*, '(A, I0)') 'Done steps = ', tsteps

    end do TIME_LOOP

    !========================================================================= 
    ! output utilities
    !=========================================================================
    ! End timing
    mpi_end_time = MPI_Wtime()

    ! Print timing (only main)
    if (rank == 0) then
        write(*, '(A)') '---------------------------------'
        write(*, '(A, F10.3, A)') '  MPI Time    = ', mpi_end_time - mpi_start_time, ' seconds.'
    end if

    ! Gather results to main
    if (rank == 0) then
        ! main: allocate full arrays
        allocate(full_con(Nx, 1:Ny))
        allocate(full_phi(Nx, 1:Ny))

        ! Copy local data
        do j = 1, local_ny
            do i = 1, Nx
                full_con(i, j) = con(i, j)
                full_phi(i, j) = phi(i, j)
            end do
        end do

        ! Receive from other processes
        do p = 1, size - 1
            p_local_Ny = Ny / size
            p_start = p * p_local_Ny + 1
            if (p == size - 1) p_local_Ny = Ny - p_start + 1

            allocate(recv_buf(1:Nx,1:p_local_Ny))
            call MPI_Recv(recv_buf, p_local_Ny * Nx, MPI_DOUBLE_PRECISION, p, 0, &
                          MPI_COMM_WORLD, status, ierr)

            do j = 1, p_local_Ny
                do i = 1, Nx
                    full_con(i, p_start + j - 1) = recv_buf(i, j)
                end do
            end do

            deallocate(recv_buf)

            allocate(recv_buf(1:Nx,1:p_local_Ny))
            call MPI_Recv(recv_buf, p_local_Ny * Nx, MPI_DOUBLE_PRECISION, p, 1, &
                          MPI_COMM_WORLD, status, ierr)

            do j = 1, p_local_Ny
                do i = 1, Nx
                    full_phi(i, p_start + j - 1) = recv_buf(i, j)
                end do
            end do

            deallocate(recv_buf)
        end do

        ! Write results to file
        open(unit=12, file='ch_ac.dat', status='replace', action='write')
        do i = 1, Nx
            write(12, *) (full_con(i, j), j = 1, Ny)
        end do
        close(12)
        write(*, '(A)') 'Results written to: ch_ac.dat'

        deallocate(full_con, full_phi)

    else
        ! Non-main processes: send local data
        allocate(send_buf(1:Nx, 1:local_ny))

        do j = 1, local_ny
            do i = 1, Nx
                send_buf(i, j) = con(i, j)
            end do
        end do

        call MPI_Send(send_buf, local_ny * Nx, MPI_DOUBLE_PRECISION, 0, 0, &
                      MPI_COMM_WORLD, ierr)

        do j = 1, local_ny
            do i = 1, Nx
                send_buf(i, j) = phi(i, j)
            end do
        end do

        call MPI_Send(send_buf, local_ny * Nx, MPI_DOUBLE_PRECISION, 0, 1, &
                      MPI_COMM_WORLD, ierr)

        deallocate(send_buf)
    end if

    ! Clean up
    deallocate(con, phi, dfdcon, dfdphi, lap_con, lap_phi)
    deallocate(dummy_con, lap_dummy, phi_dummy)

    ! Finalize MPI
    call MPI_Finalize(ierr)

end program fd_ch_ac_mpi