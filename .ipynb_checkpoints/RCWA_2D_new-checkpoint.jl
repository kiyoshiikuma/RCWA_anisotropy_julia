using LinearAlgebra

include("K_calc.jl")
include("homogeneous_eigen.jl")
include("S_calc.jl")
include("P_Q_calc.jl")
include("layer_eigen.jl")
include("Redheffer_star.jl")
include("Pol_Cin_calc.jl")



function RCWA2D_new(λ, θ, ϕ, Λx, Λy, d, E, E_1, N, M, Φ)

    #degree to rad
    θ = θ * π / 180
    ϕ = ϕ * π / 180
    Φ = Φ * π / 180

    #polalization amprittude

    # Φ = 0 p波
    # Φ = 90 s波
    
    pₑ = sin(Φ)
    pₘ = cos(Φ)

    # unit vector
    n̂ = [0, 0, -1]  
    n̂ₘ = [0, 1, 0]    

    NM = (2 * N + 1) * (2 * M + 1)

    # vacuum wave number
    k₀ = 2π / λ

    # storage
    S_matrix = [] 
    kz_matrix = []

    mᵣ = 1
    eᵣ = 1
    #n_i = √(eᵣ * mᵣ)
    n_i = 1

    # incident wave number
    kx_in = n_i * sin(θ) * cos(ϕ)
    ky_in = n_i * sin(θ) * sin(ϕ)  
    #kz_in = √(eᵣ * 1 - kx_in^2 - ky_in^2)
    kz_in = n_i * cos(θ)

    # Wave number Matrix Calculation (K_calc.jl)
    Kx, Ky = K_matrix(kx_in, ky_in, k₀, Λx, Λy, N, M);

    # W,V,Kz of 0th layer Calculation (homogeneous_eigen.jl)
    W₀, V₀, Kz₀ = homogeneous_layer(Kx, Ky, 1)

    #--------------------------------------------------------------#

    # W,V,Kz of reflection layer Calculation (homogeneous_eigen.jl)
    Wᵣ, Vᵣ, kzᵣ = homogeneous_layer(Kx, Ky, eᵣ);

    #kz_matrix.append(kzᵣ)
    append!(kz_matrix, kzᵣ) 

    #ここ一応書き換え
    #Ar, Br = A_B_matrices(W₀, Wᵣ, V₀, Vᵣ);
    Ar, Br = A_B_matrices_isotropy(W₀, Wᵣ, V₀, Vᵣ);

    # Calculation of the scattering matrix on the reflective side
    S_ref, Sr_dict = S_R(Ar, Br)  
    push!(S_matrix, S_ref)
    Sg = Sr_dict

    #--------------------------------------------------------------#

    # The scattering matrix is obtained for each layer,
    # and the final scattering matrix is obtained by performing the Redheffer star product.

    
    @views @inbounds for i in 1:length(E)
        
        # ith layer material parameters
        E_matrix = E[i]
        E_matrix_1 = E_1[i]
        

        #U = UR[i]
        mu_conv = Matrix(I, NM, NM)

        # longitudinal k_vector
        #P_Q_kz_new(Kx, Ky, E_matrix_1, E_matrix, mu_conv)
        #E_matrix_1 は Q の方に代入される
        
        P, Q, kzl = P_Q_kz_new(Kx, Ky, E_matrix_1, E_matrix, mu_conv)
        push!(kz_matrix, kzl)
        Gamma_squared = P * Q
        
        # E-field modes that can propagate in the medium, these are well-conditioned
        W_i, lambda_matrix = eigen_W(Gamma_squared)
        V_i = eigen_V(Q, W_i, lambda_matrix)

        # now define A and B, slightly worse conditioned than W and V
        A, B = A_B_matrices(W_i, W₀, V_i, V₀)  # ORDER HERE MATTERS A LOT because W_i is not diagonal

        # calculate scattering matrix
        Li = d[i]
        
        S_layer, Sl_dict = S_layers(A, B, Li, k₀, lambda_matrix)
        push!(S_matrix, S_layer)
        
        # update global scattering matrix using redheffer star
        Sg_matrix, Sg = RedhefferStar(Sg, Sl_dict)
        
    end

    
    #--------------------------------------------------------------#
    m_t = 1
    e_t = 1

    # W,V,Kz of reflection layer Calculation (homogeneous_eigen.jl)
    Wt, Vt, kz_trans = homogeneous_layer(Kx, Ky, e_t)

    #ここも書き換えた
    #At, Bt = A_B_matrices(W₀, Wt, V₀, Vt)
    At, Bt = A_B_matrices_isotropy(W₀, Wt, V₀, Vt)

    # Calculation of the scattering matrix on the reflective side
    ST, ST_dict = S_T(At, Bt, N)
    push!(S_matrix, ST)

    Sg_matrix, Sg = RedhefferStar(Sg, ST_dict)

    #--------------------------------------------------------------#


    K_in_vector = n_i * [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]

    E_in, cin, Polarization = initial_conditions(K_in_vector, θ, n̂, pₑ, pₘ, N, M)

    cin = Wᵣ \ cin
    
    reflected = Wᵣ * Sg["S11"] * cin
    transmitted = Wt * Sg["S21"] * cin

    rx = reflected[1:NM, :]  
    ry = reflected[NM+1:end, :]
    tx = transmitted[1:NM, :]
    ty = transmitted[NM+1:end, :]

    rz = -kzᵣ\ (Kx * rx + Ky * ry)
    tz = -kz_trans \ (Kx * tx + Ky * ty)
    
    # 
    r_sq = abs2.(rx) + abs2.(ry) + abs2.(rz)
    t_sq = abs2.(tx) + abs2.(ty) + abs2.(tz)
    t_sq_p = abs2.(tx) + abs2.(tz)
    
    R = real(kzᵣ) * r_sq * pinv(real(kz_in))
    T = real(kz_trans) * t_sq / real(kz_in)
    T_p = real(kz_trans) * t_sq_p / real(kz_in)

    return sum(R), sum(T) ,sum(T_p), T_p
    
end