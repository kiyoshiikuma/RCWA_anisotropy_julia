using LinearAlgebra

function K_matrix(kx_in, ky_in, k0, Λx, Λy, N_x, N_y)
    kx = kx_in .- 2π*collect(-N_x:N_x)/(k0*Λx)
    ky = ky_in .- 2π*collect(-N_y:N_y)/(k0*Λy)
    
    K̄x = repeat(kx, inner = (1, size(kx,1)))
    K̄y = transpose(repeat(ky, inner = (1,size(ky,1))))

    Kx = Diagonal(K̄x[:])
    Ky = Diagonal(K̄y[:])

    return Kx, Ky
    
end