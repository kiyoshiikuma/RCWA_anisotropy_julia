using LinearAlgebra
using BlockDiagonals

function homogeneous_layer(Kx, Ky, e_r)

    m_r = 1

    N = size(Kx,1)
    Ī = Matrix{ComplexF64}(I, N, N)
    P = (e_r^-1) * [Kx * Ky e_r * m_r * Ī - Kx^2; Ky^2 - m_r * e_r * Ī -Ky * Kx]
    Q = (e_r/m_r) * P
    W = Matrix{ComplexF64}(I, 2*N, 2*N)


    arg = m_r * e_r * Ī - Kx.^2 - Ky.^2
    arg = complex.(arg)
    Kz = conj.(.√(arg))
    
    #==
    
    Kz = sqrt.(Complex.(Ī - Kx.^2 - Ky.^2))
	     Kz[imag.(Kz).<0].*=-1
    ==#
    
    eigenvalues = BlockDiagonal([im*Kz, im*Kz])
    
    V = Q / eigenvalues
    #V = Q * pinv(eigenvalues)
    
    return W, V, Kz
    
end