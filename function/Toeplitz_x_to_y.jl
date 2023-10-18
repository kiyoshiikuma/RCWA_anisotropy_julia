using PyCall
using PyPlot
using ToeplitzMatrices, LinearAlgebra
using FFTW

scipy_ndimage = pyimport("scipy.ndimage")
np = pyimport("numpy")

# 1 / εpq making inverse of permittivity distribution 

function Fourier_permittivity_new(Nx, Ny, x_start, x_end, y_start, y_end, n, angle)
    
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
function Make_Toeplitz_matrix(Nx, Ny, x_start, x_end, y_start, y_end, n, angle, P)

    F_r = Fourier_permittivity_new(Nx, Ny, x_start, x_end, y_start, y_end, n, angle);

    nn = Int(P * 2 + 1)

    nnn = nn * nn

    Af_a = zeros(Complex{Float64}, Nx, nnn)

    for i in 1 : Nx
        
        ## do fft along x
        Af = (1 / Nx) * fftshift(fft(F_r[i, :]))

        Af_r = reverse(Af)
    
        o_0 = Int((length(Af) + 2) / 2)

        # o_0 is 0 order amplitude
        # N is 2 times of max order for Fourier 
        # P is max order for Fourier 

        N = Int(2 * P)
    
        Af_matrix = Toeplitz(Af_r[o_0 - 1 : o_0 - 1 + N], Af[o_0 : o_0 + N])

        II = Matrix{Complex{Float64}}(I, nn, nn)
    
        Af_matrix = II / Af_matrix

        # Inverse matrix in a row nn * nn Number of elements

        Af_a[i, :] = vec(Af_matrix')'

    end

    size_A = size(Af_a)

    #A_aa = fftshift(fft(Af_a, 1))/Ny

    A_aa = zeros(Complex{Float64}, Ny, nnn)
    
    for j in 1 : nnn
        
        A_aa[:, j] = (1 / Ny) * fftshift(fft(Af_a[:, j]))
    
    end

    o_0 = Int((2 + size(A_aa)[1]) / 2)

    N = Int(2 * P)

    # -2N ~ 2N Fourier along y

    A_aa = A_aa[o_0　 - N  : o_0　 + N , :];

    return A_aa
    
end

# function beta is 

function beta(array, cols)

    output_matrix = reshape(array, cols, :)

    output_matrix = output_matrix'

    return output_matrix
end


# make Toeplitz block matrix

function create_block_matrix(matrices::Vararg{AbstractMatrix{T}}) where T
    # Get the dimensions of the first input matrix
    if isempty(matrices)
        throw(ArgumentError("At least one input matrix is required."))
    end
    matrix_size = size(matrices[1])

    num_matrices = length(matrices)

    # Check that all input matrices have the same dimensions
    for matrix in matrices
        if size(matrix) != matrix_size
            throw(ArgumentError("Input matrices must have the same shape."))
        end
    end

    # Calculate the size of the block matrix
    block_matrix = zeros(T, matrix_size[1] * num_matrices, matrix_size[2] * num_matrices)

    #== miss
    # Construct the block matrix
    for i in 1:num_matrices
        for j in 1:num_matrices
            block_matrix[(i-1)*matrix_size[1]+1:i*matrix_size[1], (j-1)*matrix_size[2]+1:j*matrix_size[2]] = matrices[(i+j-2) % num_matrices + 1]
        end
    end
    ==#

    # Construct the Toeplitz for block matrix　　　　　　
    for i in 1:num_matrices
        for j in 1:num_matrices
            if i == 1
                # For the first row of blocks, use (i+j-2) % num_matrices to access matrices
                block_matrix[(i-1)*matrix_size[1]+1:i*matrix_size[1], (j-1)*matrix_size[2]+1:j*matrix_size[2]] = matrices[(i+j-2) % num_matrices + 1]
                
            else
                # For the subsequent rows, the access index should be updated
                block_matrix[(1 + num_matrices - i)*matrix_size[1]+1:(num_matrices - (i - 2))*matrix_size[1], (j-1)*matrix_size[2]+1:j*matrix_size[2]] = matrices[(i+j-2) % num_matrices + 1]
            end
        end
    end

    return block_matrix
end

# Return Toeplitz beta

create_block_matrix(matrices::Vector{Matrix{ComplexF64}}) = create_block_matrix(matrices...)

function Toeplitz_matrix_x_to_y(Nx, Ny, x_start, x_end, y_start, y_end, n, angle, P)

    # Make_Toeplitz_matrix function
    A_aa = Make_Toeplitz_matrix(Nx, Ny, x_start, x_end, y_start, y_end, n, angle, P)
    
    nn = Int(P * 2 + 1)
    nn_1 = Int(4 * P + 1)
    #nn_2 = Int((P * 2 + 1))^2
    
    beta_values = zeros(Complex{Float64}, nn_1, nn, nn)
    
    for i in 1:nn_1
        # beta function
        beta_i = beta(A_aa[i, :], nn)
        eval(Meta.parse("beta_$i = $beta_i"))
        beta_values[i, :, :] = beta_i 
    end


    
    beta_matrix = Vector{Matrix{ComplexF64}}()

    
    for i in (nn_1 + 1) ÷ 2:-1:1
        push!(beta_matrix, beta_values[i, :, :])
    end

    #===================
        
    for i in (nn_1 + 1) ÷ 2 + 1:nn_1
        push!(beta_matrix, beta_values[i, :, :])
    end
    ====================#

    for i in 1:(nn_1 + 1) ÷ 2 - 1
        
        push!(beta_matrix, beta_values[nn_1 - i + 1, :, :])
    
    end
    

    #beta_matrix = [beta_values[i, :, :] for i in 1:nn_1]

    #print(beta_matrix)

    #beta_block = bt.create_beta_n_matrix(beta_values, nn_1)
    
    beta_block = create_block_matrix(beta_matrix)

    N_0 = size(beta_1)[1]

    N_1 = N_0^2

    beta_block = beta_block[1:N_1, 1:N_1]
    
    return beta_block
    
end



        