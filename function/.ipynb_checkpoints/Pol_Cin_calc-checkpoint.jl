using LinearAlgebra

function delta_vector(P, Q)

    fourier_grid = zeros(P, Q)
    fourier_grid[div(P, 2) + 1, div(Q, 2) + 1] = 1
    vector = fourier_grid[:]
    return reshape(vector, (1, length(vector)))
    
end

function initial_conditions(K_inc_vector, theta, normal_vector, pte, ptm, P, Q)

    if theta != 0
        ate_vector = cross(normal_vector, K_inc_vector)
        ate_vector /= norm(ate_vector)
    else
        ate_vector = [0, 1, 0]
    end

    atm_vector = cross(K_inc_vector, ate_vector)
    atm_vector /= norm(atm_vector)

    Polarization = pte * ate_vector + ptm * atm_vector
    E_inc = Polarization
    Polarization = vec(Polarization)
    delta = delta_vector(2P+1, 2Q+1)

    esrc = [Polarization[1]*delta Polarization[2]*delta]
    esrc = reshape(esrc, (length(esrc), 1))

    return E_inc, esrc, Polarization
end