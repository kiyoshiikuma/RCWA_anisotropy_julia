using LinearAlgebra

# solve the eigen eq to get the W which is eigen matrix for E
function eigen_W(Gamma_squared)
    
    Lambda, W = eigen(Gamma_squared) 
    lambda_squared_matrix = Diagonal(Lambda)
    lambda_matrix = .√(complex.(lambda_squared_matrix))
    return W, lambda_matrix
    
end

# calc V which is eigen matrix for H
function eigen_V(Q, W, lambda_matrix)

    V = Q * W / lambda_matrix 

    return V
end