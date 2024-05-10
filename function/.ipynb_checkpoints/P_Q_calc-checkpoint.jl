using LinearAlgebra

function Q_matrix(Kx, Ky, e_conv_x_to_y, e_conv_y_to_x, mu_conv)
    
    # μ = 1
    a_1 = Kx * Ky
    a_2 = e_conv_y_to_x - Kx * Kx
    a_3 = Ky * Ky - e_conv_x_to_y
    a_4 = -Ky * Kx
    
    Q = [a_1 a_2; a_3 a_4]
    
    return Q
    
end

function P_matrix(Kx, Ky, e_conv, mu_conv)

    b_1 = Kx / e_conv * Ky
    b_2 = mu_conv - Kx / e_conv * Kx
    b_3 = Ky / e_conv * Ky -mu_conv
    b_4 = -Ky / e_conv * Kx
    
    P = [b_1  b_2; b_3  b_4]
    
    return P

end

function P_Q_kz_new(Kx, Ky, e_conv, e_conv_x_to_y, e_conv_y_to_x, mu_conv)

    argument = e_conv - Kx.^2 - Ky.^2
    Kz = sqrt.(Complex.(conj.(argument)))
    
    #==
    Kz = sqrt.(Complex.(e_conv - Kx.^2 - Ky.^2))
	     Kz[imag.(Kz).<0].*=-1
    ==#
    
    q = Q_matrix(Kx, Ky, e_conv_x_to_y, e_conv_y_to_x, mu_conv)
    p = P_matrix(Kx, Ky, e_conv, mu_conv)

    return p, q, Kz
end

#============================================================================================================#

# 複屈折物質　の　P, Q マトリックス (ノートの式と符号逆転してる)

function Q_matrix_bire(Kx, Ky, exx_conv_x_to_y, exy_conv_x_to_y, eyy_conv_y_to_x, eyx_conv_y_to_x, mu_conv)
    
    # μ = 1
    a_1 = eyx_conv_y_to_x + Kx * Ky
    a_2 = eyy_conv_y_to_x - Kx * Kx
    a_3 = Ky * Ky - exx_conv_x_to_y
    a_4 = -Ky * Kx - exy_conv_x_to_y
    
    Q = [a_1 a_2; a_3 a_4]
    
    return Q
    
end

function P_matrix_bire(Kx, Ky, e_conv, mu_conv)

    b_1 = Kx / e_conv * Ky
    b_2 = mu_conv - Kx / e_conv * Kx
    b_3 = Ky / e_conv * Ky -mu_conv
    b_4 = -Ky / e_conv * Kx
    
    P = [b_1  b_2; b_3  b_4]
    
    return P

end

function P_Q_kz_bire(Kx, Ky, ezz_conv, exx_conv_x_to_y, exy_conv_x_to_y, eyy_conv_y_to_x, eyx_conv_y_to_x, mu_conv)

    argument = ezz_conv - Kx.^2 - Ky.^2
    Kz = sqrt.(Complex.(conj.(argument)))
    
    #==
    Kz = sqrt.(Complex.(e_conv - Kx.^2 - Ky.^2))
	     Kz[imag.(Kz).<0].*=-1
    ==#
    
    q = Q_matrix_bire(Kx, Ky, exx_conv_x_to_y, exy_conv_x_to_y, eyy_conv_y_to_x, eyx_conv_y_to_x, mu_conv)
    p = P_matrix_bire(Kx, Ky, ezz_conv, mu_conv)

    return p, q, Kz
end