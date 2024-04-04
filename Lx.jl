## Eigenfunctions and spectra of Lₓ 
## Gianni V. Vinci
module Lx

function exp_mat2_full_scalar(X, t)
    s = 0.5 * (X[1,1] + X[2,2])
    det_X_minus_sI = (X[1,1]-s) * (X[2,2]-s)  -  X[1,2] * X[2,1]
    q = sqrt(-det_X_minus_sI)
    if(abs(q) < 1e-15)
        sinh_qt_over_q = t
    else
        sinh_qt_over_q = sinh(q*t) / q
    end
    
    cosh_qt = cosh(q*t)
    cosh_q_t_minus_s_sinh_qt_over_qt = cosh_qt - s*sinh_qt_over_q
    exp_st = exp(s*t)
    exp_Xt= complex(zeros(2, 2))
    
    # abbreviations for the case of exp_Xt referencing the same array as X
    E_00 = exp_st * (cosh_q_t_minus_s_sinh_qt_over_qt   +   sinh_qt_over_q * X[1,1])
    E_01 = exp_st * (sinh_qt_over_q * X[1,2])
    E_10 = exp_st * (sinh_qt_over_q * X[2,1])
    E_11 = exp_st * (cosh_q_t_minus_s_sinh_qt_over_qt   +   sinh_qt_over_q * X[2,2])
    
    exp_Xt[1,1] = E_00
    exp_Xt[1,2] = E_01
    exp_Xt[2,1] = E_10
    exp_Xt[2,2] = E_11
    return exp_Xt
end

function defineSpectralProblem(F,σ,xMax,xMin,Nx)
    a(x)=-F(x)
    b(x)=(σ^2)/2
    δb(x)=0.0

    Xgrid=range(xMin,stop=xMax,length=Nx) 
    Nx = size(Xgrid)[1]

    function ϕ_and_q(init_q_r, λ)       
   
        ϕ=complex(zeros(size(Xgrid)))
        q=complex(zeros(size(Xgrid)))

        #Boundary condition:
        ϕ[Nx] = init_q_r + 0.0im 
        q[Nx] = 0.0 + 0.0im 



        #Intialize A matrix
        A = complex(zeros(2, 2))
        A[1, 2] = λ
        # grid spacing, allowing non-uniform grids
        dX = Xgrid[3]-Xgrid[2]

        for k=Nx-1:-1:1 # k= Nx-1 , ... , 1

            
            xVal=(Xgrid[k+1]+Xgrid[k])/2.0

            A[2, 1]  =  1/b(xVal)
            A[2, 2] = (a(xVal)+δb(xVal))/b(xVal)
           
            exp_A_dX=exp_mat2_full_scalar(A, dX) 
            
            q[k] = exp_A_dX[1,1]*q[k+1] + exp_A_dX[1,2]*ϕ[k+1]
            ϕ[k] = exp_A_dX[2,1]*q[k+1] + exp_A_dX[2,2]*ϕ[k+1]

            
        end   
        return ϕ,q,dX,Xgrid

    end

    function ψ_and_δψ(init_ψ_r, λ)    
   

        δψ=complex(zeros(size(Xgrid)))
        ψ=complex(zeros(size(Xgrid)))
        #Boundary condition:
        δψ[1] = 0.0 + 0.0im # left boundary condition of adjoint operator
        ψ[1] = init_ψ_r + 0.0im # linear eq = arbitrary complex scaling
    
    
    
        #Intialize A matrix
        A = complex(zeros(2, 2))
        A[1, 2] = 1.0 
    
    
        for k=2:Nx
            # grid spacing, allowing non-uniform grids
            dX = Xgrid[k]-Xgrid[k-1]
            xVal=(Xgrid[k]+Xgrid[k-1])/2.0
    
            A[2, 1]  =  λ/b(xVal)
            A[2, 2] = a(xVal)/b(xVal)
    
     
            exp_A_dX=exp_mat2_full_scalar(A, dX) 
            ψ[k] = exp_A_dX[1,1]*ψ[k-1] + exp_A_dX[1,2]*δψ[k-1]
            δψ[k] = exp_A_dX[2,1]*ψ[k-1] + exp_A_dX[2,2]*δψ[k-1]
    
        end   
        return ψ,δψ
    end

    return  ϕ_and_q,ψ_and_δψ
end






end