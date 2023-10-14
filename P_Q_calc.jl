using LinearAlgebra

function Q_matrix(Kx, Ky, e_conv, mu_conv)
    
    # μ = 1
    a_1 = Kx * Ky
    a_2 = e_conv - Kx * Kx
    a_3 = Ky * Ky - e_conv
    a_4 = -Ky * Kx
    
    return [a_1 a_2; a_3 a_4]
    
end

function P_matrix(Kx, Ky, e_conv, mu_conv)

    b_1 = Kx / e_conv * Ky
    b_2 = mu_conv - Kx / e_conv * Kx
    b_3 = Ky / e_conv * Ky -mu_conv
    b_4 = -Ky / e_conv * Kx
    
    P = [b_1  b_2; b_3  b_4]
    
    return P

end

function P_Q_kz(Kx, Ky, e_conv, mu_conv)


    argument = e_conv - Kx.^2 - Ky.^2
    Kz = sqrt.(Complex.(conj.(argument)))

    
    #==
    Kz = sqrt.(Complex.(e_conv - Kx.^2 - Ky.^2))
	     Kz[imag.(Kz).<0].*=-1
    ==#
    
    q = Q_matrix(Kx, Ky, e_conv, mu_conv)
    p = P_matrix(Kx, Ky, e_conv, mu_conv)

    return p, q, Kz
end

function P_Q_kz_new(Kx, Ky, e_conv_1, e_conv, mu_conv)


    argument = e_conv - Kx.^2 - Ky.^2
    Kz = sqrt.(Complex.(conj.(argument)))
    
    
    #==
    Kz = sqrt.(Complex.(e_conv - Kx.^2 - Ky.^2))
	     Kz[imag.(Kz).<0].*=-1
    ==#
    
    q = Q_matrix(Kx, Ky, e_conv_1, mu_conv)
    p = P_matrix(Kx, Ky, e_conv, mu_conv)

    return p, q, Kz
end
