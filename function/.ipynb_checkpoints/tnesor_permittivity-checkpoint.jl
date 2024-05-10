include("/Users/ikumakiyoshi/Library/Mobile Documents/com~apple~CloudDocs/study_optics/Jupyter/2D_RCWA_julia_debug_5/function/fft_permittivity.jl")

function rot_perm(E, alpha)

    alpha = alpha * (pi/180)

    rot_mat_left = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1]

    rot_mat_right = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1]

    E_rot = rot_mat_left * E.^2 * rot_mat_right
    
    return E_rot
end

function tensor_perm_distribution(E, alpha, Nx, Ny, x_start, x_end, y_start, y_end, angle)

    E_rot = rot_perm(E, alpha)

    E_xx = Fourier_permittivity(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[1, 1], angle)

    E_xy = Fourier_permittivity(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[1, 2], angle)

    E_yy = Fourier_permittivity(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[2, 2], angle)

    E_yx = Fourier_permittivity(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[2, 1], angle)


    return E_xx, E_xy, E_yy, E_yx

end

function Toeplotz_tesor_perm_homogenous(E, alpha, P, Q, angle)

    E_rot = rot_perm(E, alpha)

    Toe_Ezz = E_rot[3, 3] * Matrix(I, (2*P + 1)*(2*Q + 1), (2*P + 1)*(2*Q + 1))
    Toe_Exx = E_rot[1, 1] * Matrix(I, (2*P + 1)*(2*Q + 1), (2*P + 1)*(2*Q + 1))
    Toe_Eyy = E_rot[2, 2] * Matrix(I, (2*P + 1)*(2*Q + 1), (2*P + 1)*(2*Q + 1))
    Toe_Exy = E_rot[1, 2] * Matrix(I, (2*P + 1)*(2*Q + 1), (2*P + 1)*(2*Q + 1))
    Toe_Eyx = E_rot[2, 1] * Matrix(I, (2*P + 1)*(2*Q + 1), (2*P + 1)*(2*Q + 1))

    return Toe_Exx, Toe_Exy, Toe_Eyy, Toe_Eyx, Toe_Ezz

end

function Toeplotz_tesor_perm(E, alpha, Nx, Ny, x_start, x_end, y_start, y_end, P, Q, angle)

    E_rot = rot_perm(E, alpha)

    E_xx = Toeplitz_matrix_x_to_y(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[1, 1], angle, P)

    E_xy = Toeplitz_matrix_x_to_y(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[1, 2], angle, P)

    E_yy = Toeplitz_matrix_y_to_x(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[2, 2], angle, P)

    E_yx = Toeplitz_matrix_y_to_x(Nx, Ny, x_start, x_end, y_start, y_end, E_rot[2, 1], angle, P)


    return E_xx, E_xy, E_yy, E_yx

end