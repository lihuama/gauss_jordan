SUBMODULE (simul_mod) simul_exec
    CONTAINS
    MODULE PROCEDURE simul
    
    REAL, PARAMETER :: EPSILON = 1.0E-6  ! A "small" number for comparison
                                     ! when determining singular eqns 

! Data dictionary: declare local variable types & definitions
    INTEGER :: n 
    REAL :: factor                       ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
    INTEGER :: irow                      ! Number of the equation currently
                                     ! being processed
    INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
    INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
    INTEGER :: kcol                      ! Index over all columns of eqn
    REAL :: temp                         ! Scratch value

    REAL, ALLOCATABLE, DIMENSION(:,:) :: a0      
    REAL, ALLOCATABLE, DIMENSION(:) :: b0        

    n = SIZE(b)
    ALLOCATE(c(n))
    a0 = a
    b0 = b
    mainloop: DO irow = 1, n
! Find peak pivot for column irow in rows irow to n
        ipeak = irow
        max_pivot: DO jrow = irow+1, n
            IF (ABS(a0(jrow,irow)) > ABS(a0(ipeak,irow))) THEN
                ipeak = jrow
            END IF
        END DO max_pivot
! Check for singular equations.  
        singular: IF ( ABS(a0(ipeak,irow)) < EPSILON ) THEN
            error = 1
            RETURN
        END IF singular
! Otherwise, if ipeak /= irow, swap equations irow & ipeak
        swap_eqn: IF ( ipeak /= irow ) THEN
            DO kcol = 1, n
                temp          = a0(ipeak,kcol)
                a0(ipeak,kcol) = a0(irow,kcol)
                a0(irow,kcol)  = temp 
            END DO
            temp     = b0(ipeak)
            b0(ipeak) = b0(irow)
            b0(irow)  = temp
        END IF swap_eqn
! Multiply equation irow by -a(jrow,irow)/a(irow,irow),  
! and add it to Eqn jrow (for all eqns except irow itself).
        eliminate: DO jrow = 1, n
            IF ( jrow /= irow ) THEN
                factor = -a0(jrow,irow)/a0(irow,irow)
                DO kcol = 1, n
                    a0(jrow,kcol) = a0(irow,kcol)*factor + a0(jrow,kcol)
                END DO
                b0(jrow) = b0(irow)*factor + b0(jrow)
            END IF
        END DO eliminate
    END DO mainloop
! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
    divide: DO irow = 1, n
        c(irow)      = b0(irow) / a0(irow,irow)
    END DO divide
! Set error flag to 0 and return.
    error = 0
    END PROCEDURE simul
    
END SUBMODULE simul_exec