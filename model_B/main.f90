!-------------------------------------------------------------------------------
!   MPI Finite Difference Phase Field Code of Cahn-Hilliard Equation.
!
!   To compile and run:
!             mpif90 -std=gnu main.f90 -fcheck=all -fbacktrace -Wall -o main
!             mpirun -np 8 ./main
!-------------------------------------------------------------------------------

program cahn_hilliard_mpi
    use mpi
    use iso_fortran_env, only: real64, dp=>real64
    implicit none
    
    !=============================================================================
    !                    PARAMETERS
    !=============================================================================
    ! simulation cell parameters
    integer, parameter :: Nx = 1024
    integer, parameter :: Ny = 1024
    real(dp), parameter :: dx = 1.0
    real(dp), parameter :: dy = 1.0
    
    ! time integration parameters
    integer, parameter :: nsteps = 2000
    integer, parameter :: nprint = 1000
    real(dp), parameter :: dt = 0.01d0
    
    ! material specific parameters
    real(dp), parameter :: con_0 = 0.4d0
    real(dp), parameter :: mobility = 1.0d0
    real(dp), parameter :: grad_coef = 0.5d0
    
    ! microstructure parameters
    real(dp), parameter :: noise = 0.02d0
    real(dp), parameter :: A = 1.0d0
    
    ! MPI variables
    integer :: rank, size, ierr
    integer :: status(MPI_STATUS_SIZE)
    
    ! MPI Domain decomposition
    integer :: local_ny, start_y, end_y
    integer :: down, up
    
    ! Arrays (allocatable for local domain with halo)
    real(dp), allocatable, dimension(:,:) :: con, con_new
    real(dp), allocatable, dimension(:,:) :: dfdcon, lap_con
    real(dp), allocatable, dimension(:,:) :: dummy_con, lap_dummy
    
    ! MPI Variables for gathering
    real(dp), allocatable :: full_con(:, :)
    real(dp), allocatable :: recv_buf(:, :)
    real(dp), allocatable :: send_buf(:, :)
    integer :: p_local_Ny, p_start, p
    
    ! Local variables
    integer :: i, j, tsteps, im, ip, jp, jm
    real(dp) :: c, new_val, r
    real(dp) :: mpi_start_time, mpi_end_time !MPI TIME VARIABLES

    !=============================================================================
    ! Initialize MPI
    !=============================================================================
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    !=============================================================================
    ! Domain decomposition along y-direction
    !=============================================================================  
    local_Ny = Ny / size
    start_y  = rank * local_Ny + 1
    end_y    = merge(Ny, start_y + local_Ny - 1, rank == size - 1)
    local_Ny = end_y - start_y + 1
    !=============================================================================      
    ! Allocate arrays with 1-based indexing for both dimensions
    !=============================================================================  
    allocate(con(1:Nx,0:local_Ny+1))
    allocate(con_new(1:Nx,0:local_Ny+1))
    allocate(dfdcon(1:Nx,0:local_Ny+1))
    allocate(lap_con(1:Nx,0:local_Ny+1))
    allocate(dummy_con(1:Nx,0:local_Ny+1))
    allocate(lap_dummy(1:Nx,0:local_Ny+1))
    !=============================================================================  
    ! Initialize arrays to zero
    !=============================================================================
    con = 0.0d0
    con_new = 0.0d0
    dfdcon = 0.0d0
    lap_con = 0.0d0
    dummy_con = 0.0d0
    lap_dummy = 0.0d0
    !=============================================================================    
    ! Initialize with noise (only local domain)
    !=============================================================================
    do j = 1, local_Ny
        do i = 1, Nx
            call random_number(r)
            con(i, j) = con_0 + noise * (0.5d0 - r)
        end do
    end do
    !=============================================================================    
    ! Initialize halos (will be updated via MPI)
    !=============================================================================  
    do i = 1, Nx
        con(i,0) = 0.0d0
        con(i,local_Ny+1) = 0.0d0
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

        if(rank==0)      down = size-1
        if(rank==size-1) up   = 0
               
        call MPI_Sendrecv( &
            con(1,1), Nx, MPI_DOUBLE_PRECISION, down, 0, &
            con(1,local_Ny+1), Nx, MPI_DOUBLE_PRECISION, up, 0, &
            MPI_COMM_WORLD,status,ierr)

        call MPI_Sendrecv( &
            con(1,local_Ny), Nx, MPI_DOUBLE_PRECISION, up, 1, &
            con(1,0), Nx, MPI_DOUBLE_PRECISION, down, 1, &
            MPI_COMM_WORLD,status,ierr)
        !=========================================================================
        ! first laplacian 
        !=========================================================================
        do j = 1, local_Ny
            do i = 1, Nx

                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1
                
                ! Periodic boundary in y
                if (ip > Nx) ip = 1
                if (im < 1)  im = Nx
                
                ! Free energy derivative
                c = con(i, j)
                dfdcon(i, j) = A * (2.0d0 * c * (1.0d0 - c) * (1.0d0 - c) - &
                                    2.0d0 * c * c * (1.0d0 - c))
                
                lap_con(i,j) = ( con(ip,j) + con(im,j) + con(i,jp) + con(i,jm) - &
                                4.0d0*c ) / (dx*dy)
                
                dummy_con(i, j) = dfdcon(i, j) - grad_coef * lap_con(i, j)
            end do
        end do
        
        ! Exchange halos for dummy_con
        call MPI_Sendrecv( &
                dummy_con(1,1), Nx, MPI_DOUBLE_PRECISION, down, 0, &
                dummy_con(1,local_Ny+1), Nx, MPI_DOUBLE_PRECISION, up, 0, &
                MPI_COMM_WORLD,status,ierr)

        call MPI_Sendrecv( &
                dummy_con(1,local_Ny), Nx, MPI_DOUBLE_PRECISION, up, 1, &
                dummy_con(1,0), Nx, MPI_DOUBLE_PRECISION, down, 1, &
                MPI_COMM_WORLD,status,ierr)
 
        !========================================================================= 
        ! second laplacian & time integration
        !=========================================================================
        do j = 1, local_Ny
            do i = 1, Nx
            
                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1
                
                ! Periodic boundary in y 
                if (ip > Nx) ip = 1
                if (im < 1)  im = Nx
                
                ! Laplacian of dummy
                lap_dummy(i, j) = (dummy_con(ip, j) + dummy_con(im, j) + &
                                   dummy_con(i, jp) + dummy_con(i, jm) - &
                                   4.0d0 * dummy_con(i, j)) / (dx * dy)
                
                ! Time integration
                new_val = con(i, j) + dt * mobility * lap_dummy(i, j)
                con_new(i, j) = min(0.99999d0, max(0.00001d0, new_val))
            end do
        end do
        
        ! Swap arrays
        con = con_new
        
        !Print progress (only main)
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
        ! main: allocate full array
        allocate(full_con(Nx, 1:Ny))
        
        ! Copy local data
        do j = 1, local_Ny
            do i = 1, Nx
                full_con(i, j) = con(i, j)
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
        end do
        
        ! Write results to file
        open(unit=10, file='ch.dat', status='replace', action='write')
        do i = 1, min(64, Nx)
            do j = 1, min(64, Ny)
                write(10, '(F12.6)', advance='no') full_con(i, j)
                if (j < min(100, Ny)) write(10, '(A)', advance='no') ' '
            end do
            write(10, *)
        end do
        close(10)
        write(*, '(A)') 'Results written to: ch.dat'
        
        deallocate(full_con)
        
    else
        ! Non-main processes: send local data
        allocate(send_buf(1:Nx, 1:local_Ny))
        
        do j = 1, local_Ny
            do i = 1, Nx
                send_buf(i, j) = con(i, j)
            end do
        end do
        
        call MPI_Send(send_buf, local_Ny * Nx, MPI_DOUBLE_PRECISION, 0, 0, &
                      MPI_COMM_WORLD, ierr)
        
        deallocate(send_buf)
    end if
    
    ! Clean up
    deallocate(con, con_new, dfdcon, lap_con, dummy_con, lap_dummy)
    
    ! Finalize MPI
    call MPI_Finalize(ierr)
    
end program cahn_hilliard_mpi