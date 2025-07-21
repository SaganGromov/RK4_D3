! filepath: /home/sagan/RK4_D3/ode_solver_rk4.f90
! =============================================================================
! Programa para resolver um sistema de EDOs usando Runge-Kutta de 4ª ordem,
! agora utilizando precisão dupla (real(8)).
!
! Mudança Final: O arquivo 'order_study.dat' agora contém apenas
!              o passo (h) e o valor final de y(t).
! =============================================================================
PROGRAM ode_solver_rk4
    IMPLICIT NONE

    ! --- Declaração de Variáveis ---
    REAL(8), PARAMETER :: y0 = 0.0d0, z0 = -acos(-1.0d0), t0 = 0.0d0
    REAL(8), PARAMETER :: t_final = 1.0d0
    INTEGER, PARAMETER :: num_refinements = 8
    
    REAL(8) :: h, t, y, z
    REAL(8) :: h_values(num_refinements)
    REAL(8) :: y_final(num_refinements)
    INTEGER :: i, j, n_steps

    ! ==========================================================================
    ! Parte 1: Resolver o sistema com um passo fixo
    ! ==========================================================================
    h = 0.1d0
    t = t0
    y = y0
    z = z0
    n_steps = NINT(t_final / h)

    OPEN(UNIT=10, FILE='solution.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(10, '(A8, A18, A18)') 't', 'y(t)', 'z(t)'
    WRITE(10, '(F8.4, 2(F18.10))') t, y, z

    DO i = 1, n_steps
        CALL rk4_step(t, y, z, h)
        t = t0 + i * h
        WRITE(10, '(F8.4, 2(F18.10))') t, y, z
    END DO
    CLOSE(10)

    PRINT *, '== Trajetória da Solução =='
    PRINT *, 'Trajetória da solução salva em "solution.dat".'
    PRINT *, ''

    ! ==========================================================================
    ! Parte 2: Gerar dados de h vs. y(t_final)
    ! ==========================================================================
    
    ! --- Calcular y(t_final) para uma faixa de passos decrescentes ---
    h = 0.2d0 ! Passo inicial para o estudo
    DO i = 1, num_refinements
        h_values(i) = h
        t = t0
        y = y0
        z = z0
        n_steps = NINT(t_final / h)

        DO j = 1, n_steps
            CALL rk4_step(t, y, z, h)
        END DO

        ! Armazenar o valor final de y
        y_final(i) = dabs(y - sin(acos(-1.0d0))) ! Ajuste para o valor esperado

        h = h / 2.0d0 ! Reduzir o passo pela metade
    END DO

    ! --- Gravar os dados de h vs. y(t_final) em um arquivo ---
    OPEN(UNIT=11, FILE='order_study.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(11, '(A15, A25)') '# h', 'y(t=1.0)'
    DO i = 1, num_refinements
        WRITE(11, '(E15.6, F25.15)') h_values(i), y_final(i)
    END DO
    CLOSE(11)

    PRINT *, '== Estudo de Convergência de y(t) =='
    PRINT *, 'Dados de h vs. y(t) salvos em "order_study.dat".'

CONTAINS

    ! --- Função para y' = f(t,y,z) ---
    REAL(8) FUNCTION f_func(t, y, z)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: t, y, z
        f_func = z
    END FUNCTION f_func

    ! --- Função para z' = g(t,y,z) ---
    REAL(8) FUNCTION g_func(t, y, z)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: t, y, z
        g_func = -acos(-1.0d0)**2 * y
    END FUNCTION g_func

    ! --- Sub-rotina para um passo do método Runge-Kutta de 4ª ordem ---
    SUBROUTINE rk4_step(t_in, y_inout, z_inout, h)
        IMPLICIT NONE
        REAL(8), INTENT(IN)    :: t_in, h
        REAL(8), INTENT(INOUT) :: y_inout, z_inout
        REAL(8) :: k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z

        k1y = f_func(t_in, y_inout, z_inout)
        k1z = g_func(t_in, y_inout, z_inout)
        
        k2y = f_func(t_in + h/2.0d0, y_inout + h/2.0d0 * k1y, z_inout + h/2.0d0 * k1z)
        k2z = g_func(t_in + h/2.0d0, y_inout + h/2.0d0 * k1y, z_inout + h/2.0d0 * k1z)

        k3y = f_func(t_in + h/2.0d0, y_inout + h/2.0d0 * k2y, z_inout + h/2.0d0 * k2z)
        k3z = g_func(t_in + h/2.0d0, y_inout + h/2.0d0 * k2y, z_inout + h/2.0d0 * k2z)

        k4y = f_func(t_in + h, y_inout + h * k3y, z_inout + h * k3z)
        k4z = g_func(t_in + h, y_inout + h * k3y, z_inout + h * k3z)

        y_inout = y_inout + (h/6.0d0) * (k1y + 2.0d0*k2y + 2.0d0*k3y + k4y)
        z_inout = z_inout + (h/6.0d0) * (k1z + 2.0d0*k2z + 2.0d0*k3z + k4z)
    END SUBROUTINE rk4_step

END PROGRAM ode_solver_rk4