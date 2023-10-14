using LinearAlgebra

function A(W_layer, Wg, V_layer, Vg)

    A = W_layer \ Wg + V_layer \ Vg
    
    return A
end

function B(W_layer, Wg, V_layer, Vg)

    B = W_layer \ Wg - V_layer \ Vg
    
    return B
end

function A_B_matrices(W_layer, Wg, V_layer, Vg)
    
    a = A(W_layer, Wg, V_layer, Vg)
    b = B(W_layer, Wg, V_layer, Vg)
    
    return a, b
end

function A_B_matrices_isotropy(W_layer, Wg, V_layer, Vg)
    
    N = size(W_layer,1)
    a = 2 .* Matrix(I,N,N);
    b = 0 .* Matrix(I,N,N);
    
    return a, b
end

function S_layers(A, B, Li, k0, modes)
    
    X_i = Diagonal(exp.(-diag(modes) * Li * k0))
    term1 = A - X_i * B * (A \ X_i) * B
    S11 = term1 \ (X_i * B * (A \ X_i) * A - B)
    S12 = term1 \ ((X_i) * (A - B * (A \ B)))
    S22 = S11
    S21 = S12
    S_dict = Dict("S11" => S11, "S22" => S22, "S12" => S12, "S21" => S21)
    S = [S11 S12; S21 S22]
    
    return S, S_dict
end

function S_R(Ar, Br)

    S11 = -Ar \ Br
    S12 = 2 * pinv(Ar)
    S21 = 0.5 * (Ar - Br / Ar * Br)
    S22 = Br * pinv(Ar)

    S_dict = Dict("S11" => S11, "S22" => S22, "S12" => S12, "S21" => S21)
    S = [S11 S12; S21 S22]
    
    return S, S_dict
end

function S_T(At, Bt, N)

    In = Matrix(I,N,N)

    S11 = Bt / At
    S21 = 2 * inv(At)
    S12 = 0.5 * (At - Bt * (At \ Bt))
    S22 = -At \ Bt
    S_dict = Dict("S11" => S11, "S22" => S22, "S12" => S12, "S21" => S21)
    S = [S11 S12; S21 S22]
    
    return S, S_dict
end