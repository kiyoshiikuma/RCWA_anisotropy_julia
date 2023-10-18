using PyCall
using PyPlot
using ToeplitzMatrices, LinearAlgebra
using FFTW

scipy_ndimage = pyimport("scipy.ndimage")

# 1 / εpq making inverse of permittivity distribution 

function Fourier_permittivity_new_y(Nx, Ny, x_start, x_end, y_start, y_end, n, angle)
    
    A = ones(Float64, Nx, Ny) 

    A[y_start + 1:y_end, x_start + 1:x_end] .= n
    #A[x_start:x_end, y_start:y_end] .= n

    # Calculate the center of the rectangle
    center_x = div(x_start + x_end, 2)
    center_y = div(y_start + y_end, 2)

    # Rotate rectangle around its center using scipy's rotate function
    A_rotated = scipy_ndimage.rotate(A, angle, reshape=false, order=1, cval=1)

    # permittivity figure
    PyPlot.imshow(A_rotated, cmap="viridis", origin="lower")
    PyPlot.colorbar(label="permittivity")
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    PyPlot.grid(true, linestyle="dotted", color="white")
    PyPlot.title("Permittivity distribution")
    PyPlot.show()

    # making 1 / εpq

    A_rotated_inverse = 1 ./ A_rotated;
    
    return A_rotated_inverse
end


# Make the Toeplitz matrix
# Return 
function Make_Toeplitz_matrix_y(Nx, Ny, x_start, x_end, y_start, y_end, n, P, angle)

    F_r = Fourier_permittivity_new_y(Nx, Ny, x_start, x_end, y_start, y_end, n, angle);

    nn = Int(P * 2 + 1)

    nnn = nn * nn

    Af_a = zeros(Complex{Float64}, nnn, Ny)

    for i in 1 : Ny
        ## do fft yに沿っって
        Af = (1 / Ny) * fftshift(fft(F_r[:, i]'))

        Af_r = reverse(Af)

        #print(size(Af))
    
        o_0 = Int((length(Af) + 2) / 2)

        # o_0 is 0 order amplitude
        # N is 2 times of max order for Fourier 
        # P is max order for Fourier 

        N = Int(2 * P)

        # make the Toeplitz matrix center part is o_0
    
        Af_matrix = Toeplitz(Af_r[o_0 - 1 : o_0 - 1 + N], Af[o_0 : o_0 + N])

        #print(size(Af_matrix))

        II = Matrix{Complex{Float64}}(I, nn, nn)

        # calculate inverse of Toeplitz matrix
    
        Af_matrix = II / Af_matrix

        # 逆行列を 1行　 に並べる  nn * nn 要素数

        Af_a[:, i] = vec(Af_matrix)
    
    end

    size_A = size(Af_a)
    
    A_aa = zeros(Complex{Float64}, nnn, Nx)
    
    for i in 1 : nnn
        
        A_aa[i, :] = (1 / Nx) * fftshift(fft(Af_a[i, :]))
    
    end
        
    o_0 = Int((2 + size(A_aa)[2]) / 2)

    N = Int(2 * P)

    A_aa = A_aa[: , o_0　 - N  : o_0　 + N];

    nn = Int(P * 2 + 1)

    nnn = nn * nn

    # Toeplitz　行列作成 gamma_mnを作っていく
    nn_1 = Int(4*P + 1)

    #center
    c_0 = Int((nn_1 + 1) / 2)

    gamma_values = zeros(Complex{Float64}, nnn, nn, nn)
    
    for i in 1:nnn

        # i 行目の要素等を取得
        gamma = A_aa[i, c_0:nn_1]
    
        gamma_re = reverse(A_aa[i, 1:c_0])
    
        gamma_i = Toeplitz(gamma_re, gamma)
    
        eval(Meta.parse("gamma_$i = $gamma_i"))

        gamma_values[i, :, :] = gamma_i 
        
    end

    gamma_matrix = [gamma_values[i, :, :] for i in 1:nnn];
    
    return gamma_matrix
end


# make Toeplitz block matrix

function gamma_block_matrix(N::Int, matrices::Vararg{AbstractMatrix{T}}) where T
    
    if length(matrices) % N != 0
        throw(ArgumentError("The number of matrices must be a multiple of N"))
    end

    num_blocks = length(matrices) ÷ N  # ブロックの数を計算
    block_matrix = Matrix{T}(undef, size(matrices[1], 1) * num_blocks, size(matrices[1], 2) * N)

    for i in 1:num_blocks
        block_start = (i - 1) * N + 1
        block_end = i * N
        block = hcat(matrices[block_start:block_end]...)
        block_matrix[(i - 1) * size(block, 1) + 1:i * size(block, 1), :] = block
    end

    return block_matrix
end

create_block_matrix(matrices::Vector{Matrix{ComplexF64}}) = create_block_matrix(matrices...)

function Toeplitz_matrix_y_to_x(Nx, Ny, x_start, x_end, y_start, y_end, n, P, angle)

    # Make_Toeplitz_matrix function
    
    nn = Int(P * 2 + 1)
    nn_1 = Int(4 * P + 1)

    gamma_matrix_1 = Make_Toeplitz_matrix_y(Nx, Ny, x_start, x_end, y_start, y_end, n, P, angle)

    matrices = [gamma_matrix_1[i] for i in 1:length(gamma_matrix_1)]

    gamma_block = gamma_block_matrix(nn, matrices...)
    
    return gamma_block
    
end



        