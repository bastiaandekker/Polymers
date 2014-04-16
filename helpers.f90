! **********************************************************************************
! Module containing helper functions for rosenbluth.f90
! **********************************************************************************
module helpers
  implicit none
  private

  public WriteResults, GetParameters

contains
! **********************************************************************************   
    ! Reads 'rosenbluth.params'
    subroutine GetParameters(numBeads,numPolymerStarts,numTheta,eps,sigma,&
                                temperature,usePerm,plotPolymer)
        real(8), intent(out) :: eps,sigma,temperature
        integer, intent(out) :: numBeads,numPolymerStarts,numTheta
        logical, intent(out) :: usePerm,plotPolymer

        open(12, file="rosenbluth.params")
        read(12,*) numBeads         ! Max bead length
        read(12,*) numPolymerStarts ! Number times the recursive algoritm is called
        read(12,*) numTheta         ! Number of angles in the algorthm
        read(12,*) eps              ! Epsilon    
        read(12,*) sigma            ! Sigma
        read(12,*) temperature      ! Temperature
        read(12,*) usePerm          ! Should PERM be used?
        read(12,*) plotPolymer      ! Should polymers be plotted
        close(12)
    end subroutine GetParameters
    
! **********************************************************************************
    ! Calculates the final results and writes them to a .txt file
    subroutine WriteResults(usePerm, polymerCount, weights, r2End2End, numBeads, &
                            numPolymerStarts)
        logical, intent(in) :: usePerm
        integer, intent(in) :: polymerCount(:), numBeads, &
                                numPolymerStarts
        real(8), intent(in) :: weights(:, :), r2End2End(:,:)
        
        integer iBead
        character(40) :: strFileName
        character(25) :: strAlgroritm
        
        if (usePerm) strAlgroritm = 'PERM'
        if (.not. usePerm) strAlgroritm = 'Rosenbluth'
        write(strFileName, '(a, a,I0, a, I0, a)') &
            trim(strAlgroritm), '_Beads', numBeads, '_Pol', numPolymerStarts, '.txt'
        write(strFileName, '(a,a)') trim(strAlgroritm), 'Calc.txt'
       
        open(15, file=strFileName)
        do iBead = 3, numBeads! todo: smaller file?-> datablocking
            write (15, *) iBead, polymerCount(iBead), &
                                  weights(iBead,1:maxval(polymerCount)) 
            write (15, *) iBead, polymerCount(iBead), &
                                r2End2End(iBead,1:maxval(polymerCount))
        end do
    
        close(15)
        write( *, '(a, a)') 'Results saved in ', strFileName
        write(*,*)
    end subroutine WriteResults

end module helpers
! **********************************************************************************
! **********************************************************************************
