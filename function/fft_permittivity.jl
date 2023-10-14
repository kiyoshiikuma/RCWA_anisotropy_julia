using PyCall
using PyPlot

scipy_ndimage = pyimport("scipy.ndimage")

function Fourier_permittivity(Nx, Ny, x_start, x_end, y_start, y_end, n, angle)
    
    A = ones(Complex{Float64}, Nx, Ny) 

    A[y_start + 1:y_end, x_start + 1:x_end] .= n
    #A[x_start:x_end, y_start:y_end] .= n

    # Calculate the center of the rectangle
    center_x = div(x_start + x_end, 2)
    center_y = div(y_start + y_end, 2)

    # Rotate rectangle around its center using scipy's rotate function
    A_rotated = scipy_ndimage.rotate(A, angle, reshape=false, order=1, cval=1)

    # permittivity figure
    PyPlot.imshow(real.(A_rotated), cmap="viridis", origin="lower")
    PyPlot.colorbar(label="permittivity")
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    PyPlot.grid(true, linestyle="dotted", color="white")
    PyPlot.title("Permittivity distribution")
    PyPlot.show()

    return A_rotated
end