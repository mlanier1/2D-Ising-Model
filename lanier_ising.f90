program lanier_ising
    implicit none
    double precision, allocatable, dimension (:,:) :: State, spin
    double precision :: T, deltaE, deltaM, chance, Rx, Ry, temp
    double precision :: M, E, TotalE, TotalM, TotalEsq, TotalMsq
    double precision :: mu, k, j, mB, B, curie
    integer :: numtemps, bigNum, bigNum2, it, i, l, N, N2, ix, iy
    integer :: size, seed
    double precision :: heatCapacity, magSus, stdDev
    double precision :: sumE, sumE2, sumM, sumM2, avgE, avgE2, avgM, avgM2
    integer :: x, y, ri, rj
    double precision :: DE
    character(len=50) :: energyTempFile, magTempFile, energyBFieldFile, magBFieldFile, heatTempFile

    ! Parameters for the Ising Model
    N = 10 ! Lattice dimension (N x N2)
    N2 = 10 ! Lattice dimension (N x N2)
    mu = 4.89 * 9.274E-24 ! Magnetic moment - 4.89 * Bohr magneton [J/T]
    numtemps = 100 ! Number of temperature points to simulate
    deltaE = 16 ! Temperature step size
    ! bigNum = 100000000 NOTE: Larger values improve accuracy but take longer.
    bigNum = 1000 ! Number of Monte Carlo steps
    k = 1.380649E-23 ! Boltzmann Constant [J/K]
    B = 1 ! External magnetic field
    j = (3.0 / 16.0) * k * 1043.0 ! Exchange interaction constant (strength of spin coupling)
    ! NOTE: j=0 for paramagnetic.

    ! Parameters for the Extended Ising Model
    size = 25 ! Lattice dimensionion (size x size)
    mB = 2.74009994d-24 ! Magnetic moment for extended Ising model (different from mu)
    curie = 1043.2 ! Curie temperature [K]
    ! bigNum2 = 500000 NOTE: Larger values improve accuracy but take longer.
    bigNum2 = 500 ! Number of Monte Carlo steps

    ! Initialize random seed
    call random_seed()

    ! Allocate and initialize Ising Model state
    allocate(State(0:N + 1, 0:N2 + 1)) ! Dynamic allocation of a 2D array. Use extra rows and columns for boundary conditions
    State = 0 ! Initialized to known state

    ! Allocate and initialize Extended Ising Model state
    allocate(spin(size, size))

    ! Output file paths
    energyTempFile = 'IdealEnergyVsTemperature.txt'
    magTempFile = 'IdealMagnetizationVsTemperature.txt'
    energyBFieldFile = 'IdealEnergyVsMagneticField.txt'
    magBFieldFile = 'IdealMagnetizationVsMagneticField.txt'
    heatTempFile = 'IdealSpecificHeatVsTemperature.txt'

    ! Open output files
    open(10, file=energyTempFile, status='replace')
    open(20, file=magTempFile, status='replace')
    open(30, file=energyBFieldFile, status='replace')
    open(40, file=magBFieldFile, status='replace')
    open(50, file=heatTempFile, status='replace')

    ! Initial conditions for single temperature step
    do it = 1, numtemps
        T = it * deltaE ! Current temperature T incremented by deltaE for each iteration

        ! Initialize Ising Model state
        do i = 1, N
            do l = 1, N2
                State(i, l) = 1 ! Set all spins to +1 (up-spin)
            end do
        end do

        ! Calculate energy contribution from interaction between spins and external magnetic field B
        E = -N * N2 * B * mu ! mu = (look up value)
        E = E - (4 * j * 2) - (2 * (N - 2) * j * 3) - (2 * (N2 - 2) * j * 3) - ((N - 2) * (N2 - 2) * j * 4)
        ! Terms:
        ! 4 * j * 2 = edges with 2 neighbors
        ! 2 * (N - 2) * j * 3 = interior edges with 3 neighbors
        ! 2 * (N2 - 2) * j * 3 = similar contribution for the vertical direction
        ! ((N - 2) * (N2 - 2) * j * 4) = interior spins with 4 neighbors

        ! Initialize total magnetization M of the lattice
        M = N * N2 ! NOTE: Since all spins are up-spin initially, the magnetization is simply the number of spins.
        ! Initialize accumlators for energy and magnetization sums
        TotalE = E 
        TotalM = M
        TotalEsq = E**2
        TotalMsq = M**2

        ! Perform Metropolis Hastings algorithm for Ising Model
        do i = 1, bigNum
            call random_number(Rx) ! Generates random number between 0 and 1
            ix = INT(Rx * N) + 1 ! Selects a random x-coordinate on the lattice
            call random_number(Ry) ! Generates another random number
            iy = INT(Ry * N2) + 1 ! Selects a random y-coordinate on the lattice
            State(ix, iy) = -1 * State(ix, iy) ! Flip the spin at (ix, iy) NOTE: Flip is temporary until the energy change deltaE is evaluated.

            ! Calculate energy change caused by flipping the spin at (ix, iy)
            deltaE = -2 * B * mu * State(ix, iy) - 4 * j * (State(ix, iy) * State(ix - 1, iy) + State(ix, iy) * State(ix + 1, iy) &
                + State(ix, iy) * State(ix, iy - 1) + State(ix, iy) * State(ix, iy + 1))
            ! Updates the magneization as the flipped spin changes from +1 to -1 or vice versa.
                deltaM = 2 * State(ix, iy)

            call random_number(chance) ! Generate a random number for acceptance criteria
            if (chance < exp(-deltaE / (k * T))) then
                ! NOTE: If the energy decreases, the flip is always accepted.
                E = E + deltaE ! Accept the flip: update the total energy
                M = M + deltaM ! Accept the flip: update the total magnetization
            else
                State(ix, iy) = -1 * State(ix, iy) ! Reject the flip: revert the spin
            end if

            ! Updates statistical accumulators for energy and magnetization
            sumE = TotalE + E
            sumM = TotalM + M
            TotalEsq = TotalEsq + E**2
            TotalMsq = TotalMsq + M**2
        end do

        ! Outputs average energy per spin for the current temperature
        write(10, *) T, TotalE / (bigNum + 1)

        ! Perform Metropolis Hastings algorithm for Extended Ising Model
        E = 0.0 ! Reset energy 
        M = 0.0 ! Reset magnetization

        ! Randomly initialize the lattice spins for the extended model
        do x = 1, size
            do y = 1, size
                call random_number(temp) ! Generate random value
                if (temp .le. 0.5) then
                    spin(x, y) = -1 ! Set spin to -1
                else
                    spin(x, y) = 1 ! Set spin to +1
                end if
            end do
        end do

        ! Calculate initial energy and magnetization of extended model
        do x = 1, size
            do y = 1, size
                E = E - j / 2.0 * (spin(mod(x, size) + 1, y) + spin(x, mod(y, size) + 1) + spin(mod(x + size - 2, size) + 1, y) &
                    + spin(x, mod(y + size - 2, size) + 1)) * spin(x, y) - mB * B * spin(x, y)
                M = M + spin(x, y)
            end do
        end do

        ! Perform Metropolis Hastings algorithm for the extended model
        ! NOTE: First loop - setting initial averages
        ! Allows system to reach thermal equilibrium at the current temperature T
        ! Stabilizes the system before starting the measurement phase
        do x = 1, bigNum2
            call random_number(temp)
            ri = floor(temp * size + 1) ! Select random spin (row)
            call random_number(temp)
            rj = floor(temp * size + 1) ! Select random spin (column)
            spin(ri, rj) = -spin(ri, rj) ! Temporarily flip the spin

            DE = -2 * j * (spin(mod(ri, size) + 1, rj) + spin(ri, mod(rj, size) + 1) + spin(mod(ri + size - 2, size) + 1, rj) &
                + spin(ri, mod(rj + size - 2, size) + 1)) * spin(ri, rj) - mB * B * spin(ri, rj)

            if (DE < 0) then ! Always accept if energy decreases
                E = E + DE
                M = M + 2 * spin(ri, rj)
            else
                call random_number(temp)
                if (temp < exp(-DE / (k * T))) then
                    E = E + DE ! Accept with probability exp(-DE / (k * T))
                    M = M + 2 * spin(ri, rj)
                else
                    spin(ri, rj) = -spin(ri, rj) ! Reject: revert spin
                end if
            end if
        end do

        sumE = 0.0
        sumE2 = 0.0
        sumM = 0.0
        sumM2 = 0.0

        ! NOTE: Second loop - collecting data
        ! Collects statistical data for the energy and magnetization after the system reached equilibrium
        do x = 1, bigNum2
            call random_number(temp)
            ri = floor(temp * size + 1)
            call random_number(temp)
            rj = floor(temp * size + 1)
            spin(ri, rj) = -spin(ri, rj)

            DE = -2 * j * (spin(mod(ri, size) + 1, rj) + spin(ri, mod(rj, size) + 1) + spin(mod(ri + size - 2, size) + 1, rj) &
                + spin(ri, mod(rj + size - 2, size) + 1)) * spin(ri, rj) - mB * B * spin(ri, rj)

            if (DE < 0) then
                E = E + DE
                M = M + 2 * spin(ri, rj)
            else
                call random_number(temp)
                if (temp < exp(-DE / (k * T))) then
                    E = E + DE
                    M = M + 2 * spin(ri, rj)
                else
                    spin(ri, rj) = -spin(ri, rj)
                    sumE = sumE + E
                    sumE2 = sumE2 + E**2
                    sumM = sumM + M
                    sumM2 = sumM2 + M**2
                end if
            end if
        end do

        ! Calculates averages and thermodynamic quantities
        avgE = 0.0
        avgE2 = 0.0
        avgM = 0.0
        avgM2 = 0.0
        stdDev = 0
        avgE = sumE / bigNum2
        avgE2 = sumE2 / bigNum2
        avgM = sumM / bigNum2
        avgM2 = sumM2 / bigNum2

        heatCapacity = (1.0 / (k * T**2)) * (avgE2 - avgE**2) ! Specific heat [J/K]
        magSus = (avgM2 - avgM**2) / (k * T) ! Magnetic susceptibility [1/J]

        ! Write results for Extended Ising Model to files
        write(20, *) T, avgE
        write(30, *) B, avgE
        write(40, *) B, abs(avgM)
        write(50, *) T, heatCapacity

    end do

    ! Close output files
    close(10)
    close(20)
    close(30)
    close(40)
    close(50)

end program lanier_ising