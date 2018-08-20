MODULE simul_mod
    INTERFACE
    
    MODULE SUBROUTINE simul(a,b,c,error)
    IMPLICIT NONE
    REAL, ALLOCATABLE, DIMENSION(:,:) :: a      ! Array of coefficients (n x n).
    REAL, ALLOCATABLE, DIMENSION(:) :: b        ! Input: Right-hand side of eqns.
    REAL, ALLOCATABLE, DIMENSION(:) :: c        ! Output: Solution vector.
    INTEGER :: error        ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations
    INTENT(IN) :: a,b
    INTENT(OUT) :: c,error
    END SUBROUTINE
    
    END INTERFACE
END MODULE
