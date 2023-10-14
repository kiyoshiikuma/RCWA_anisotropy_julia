using LinearAlgebra

function eigen_W(Gamma_squared)
    
    Lambda, W = eigen(Gamma_squared) 
    lambda_squared_matrix = Diagonal(Lambda)
    lambda_matrix = .√(complex.(lambda_squared_matrix))
    return W, lambda_matrix

end

function eigen_V(Q, W, lambda_matrix)
    
    V = Q * W / lambda_matrix 

    return V
end