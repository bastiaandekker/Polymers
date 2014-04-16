! **********************************************************************************
! B. Dekker - ICCP - March 2014
! Main program of PERM algorithm. Generates and anlyzes an ensemble of polymers.
! **********************************************************************************

program Rosenbluth
    use helpers ! For initializing and writing the results
    use plot ! For plotting polymers
    implicit none
    
    real(8), dimension(:,:), allocatable :: polymerPos, r2End2End, weights
    integer, dimension(:), allocatable :: polymerCount, cloneCount, pruneCount
    real(8) :: temperature, sigma, sigma6, sigma12, eps, pi, polymerStartWeight
    integer :: numBeads, numPolymerStarts, numTheta
    integer :: iPolymer, maxPolymerCount
    logical :: usePerm, plotPolymer


    ! Get user defined parameters. See rosenbluth.f90 for further info.
    call getParameters(numBeads, numPolymerStarts, numTheta, eps, sigma, &
                                temperature, usePerm, plotPolymer)
    ! Initialize all arrays and other variables
    call Initialize()


    ! Polymers are grown by calling the recursive routine AddBead
    ! The beads positions are stored in 'polymerPos'
    polymerPos(1,:) = (/0.d0, 0.d0/)
    polymerPos(2,:) = (/1.d0, 0.d0/)
    DO iPolymer=1, numPolymerStarts
        polymerStartWeight = 1.d0
        call AddBead(polymerPos, 3, polymerStartWeight)
    end do
    

    ! Results are saved to an .txt file
    call WriteResults(usePerm, polymerCount, weights, r2End2End, numBeads, &
                        numPolymerStarts)
    
    contains

! **********************************************************************************
    ! Main recursive routine that grows the polymers, and implements PERM
    recursive subroutine AddBead(polymerPos, iNewBead, polymerWeight)
        real(8), intent(inout) :: polymerPos(numBeads,2), polymerWeight
        integer, intent(in) :: iNewBead
        real(8) :: thetaWeightsSum, rand, newPolymerWeight, &
                    permUplim, permLowLim, averageWeight, weight3Beads
        
        ! Get the position of the new bead
        polymerPos(iNewBead, :) = GetNewPosition(polymerPos, iNewBead, &
                                    thetaWeightsSum)
        
        ! Store the end-to-end distance squared of the newly created polymer
        polymerCount(iNewBead) = polymerCount(iNewBead)+1 
        r2End2End(iNewBead, polymerCount(iNewBead)) = &
            sum(PolymerPos(iNewBead,:)*PolymerPos(iNewBead,:))
        
        ! Calculate and store the weight of the polymer
        polymerWeight = polymerWeight*thetaWeightsSum
        if (usePerm) polymerWeight = polymerWeight/(0.75*numTheta) ! Scale for PERM
        weights(iNewBead, polymerCount(iNewBead)) = polymerWeight
        
        ! Start plotting finished polymer if required
        if (plotPolymer .and. iNewBead .eq. numBeads &
              .and. polymerCount(numBeads) > int(numPolymerStarts/3)) then
            call drawNicePolymer(polymerPos)
        end if
        
        ! Add next bead according to selected alogrithm
        if (usePerm .and. iNewBead < numBeads) then ! PERM
            ! Calculate limits for pruning and cloning
            averageWeight = sum(weights(iNewBead,1:polymerCount(iNewBead)))
            weight3Beads = sum(weights(3,1:polymerCount(3)))
            permUpLim = 2.0*averageWeight/weight3Beads
            permLowLim = 1.2*averageWeight/weight3Beads
            if (polymerWeight > permUpLim) then ! Clone
                if (polymerCount(iNewBead) >= maxPolymerCount) print *, &
                     'Increase maxPolymerCount'
                cloneCount(iNewBead) = cloneCount(iNewBead)+1
                newPolymerWeight = 0.5d0*polymerWeight
                call AddBead(polymerPos, iNewBead+1, newPolymerWeight)
                newPolymerWeight = 0.5d0*polymerWeight
                call AddBead(polymerPos, iNewBead+1, newPolymerWeight)
            else if (polymerWeight < permLowLim) then ! Prune
                call random_number(rand)
                if (rand < 0.5d0) then
                    pruneCount(iNewBead) = pruneCount(iNewBead)+1
                    newPolymerWeight = 2.d0*polymerWeight
                    call AddBead(polymerPos, iNewBead+1, newPolymerWeight)
                end if
            else
                call AddBead(polymerPos, iNewBead+1, polymerWeight)
            end if
        else if ((.not. usePerm) .and. iNewBead < numBeads) then ! No PERM
            call AddBead(polymerPos, iNewBead+1, polymerWeight)
        end if
    end subroutine AddBead

! **********************************************************************************
    ! Subroutine of addbead(). Returns the position of the new bead
    function GetNewPosition(polymerPos, iNewBead, thetaWeightsSum)
        real(8), intent(in) :: polymerPos(:,:)
        integer, intent(in) :: iNewBead
        real(8), intent(out) :: thetaWeightsSum
        real(8), dimension(2) :: GetNewposition
        real(8) thetas(numTheta), thetaPos(numTheta, 2), thetaWeights(numTheta), &
                thetaPrbs(numTheta), thetaStart, rand, probSum
        integer :: iTheta, iNewTheta = 1       
        
        ! Generate random angles for new positions and calculate probabilities
        call random_number(rand)       
        thetaStart = rand*2*pi
        do iTheta = 1, numTheta
            thetas(iTheta) = thetaStart + 2*pi*(dble(iTheta-1))/dble(numTheta)
            thetaPos(iTheta, :) = polymerPos(iNewBead-1,:) + &
                                    (/cos(thetas(iTheta)), sin(thetas(iTheta))/)
            thetaWeights(iTheta) = exp(-PotentialNewBead(thetaPos(iTheta, :), &
                                    iNewBead, polymerPos)/temperature)
        end do
        thetaWeightsSum = sum(thetaWeights)
        thetaPrbs = thetaWeights/thetaWeightsSum
        
        ! Pick an angle according to weights
        call random_number(rand)
        probSum = 0
        do iNewTheta = 1, numTheta
            probSum = probSum + thetaPrbs(iNewTheta)
            if (rand < probSum) exit
        end do
        if (iNewTheta > numTheta) print *, 'Error, no angles available!'
        
        ! Pass the new position
        GetNewposition = thetaPos(iNewTheta,:)
    end function GetNewPosition

! **********************************************************************************
    ! Subroutine of GetNewPosition(). Returns the potential energy associated with
    ! a potential location for a bead.
    real(8) function PotentialNewBead(newPos, iNewBead, polymerPos)
        real(8) :: polymerPos(numBeads, 2), NewPos(2), r(2), r2Inv
        integer :: iBead, iNewBead
        PotentialNewBead = 0.d0
        do iBead=1, iNewBead-2
          r = PolymerPos(iBead,:)-newPos(:)
          r2Inv = 1.d0/sum(r*r)          
          PotentialNewBead = PotentialNewBead + &
            4.d0*eps*(sigma12*r2Inv**6 - sigma6*r2Inv**3)
        end do
        if (PotentialNewBead>500._8) then
            PotentialNewBead=500._8
        end if
    end function PotentialNewBead

! **********************************************************************************
    ! Main initilizing routine
    subroutine Initialize
        ! Calculate constants
        pi = 4.d0*DATAN(1.d0)
        sigma6 = Sigma**6
        sigma12 = Sigma**12
        maxPolymerCount = 10*numPolymerStarts ! Important for max arraysize required
        
        ! Allocate arrays
        allocate(polymerPos(numBeads,2)) 
        allocate(r2End2End(numBeads, maxPolymerCount), &
                weights(numBeads, maxPolymerCount))
        allocate(polymerCount(numBeads), pruneCount(numBeads), cloneCount(numBeads)) 
        
        ! Initilize random generator and integers
        call random_seed
        polymerCount = 0
        pruneCount = 0
        cloneCount = 0
    end subroutine Initialize

end program Rosenbluth
! **********************************************************************************
! **********************************************************************************
