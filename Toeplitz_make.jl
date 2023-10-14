using FFTW

function toeplitz(A, P, Q)
    N = size(A)

    NH = (2*P+1) * (2*Q+1)
    p = -P:P
    q = -Q:Q

    ## do fft
    Af = (1 / (N[1] * N[2])) * fftshift(fft(A))
    
    # central indices marking the (0,0) order
    p0 = div(N[2], 2)
    q0 = div(N[1], 2)

    C = zeros(ComplexF64, NH, NH)
    
    for qrow in 1:2*Q+1
        for prow in 1:2*P+1
            row = (qrow) * (2*P+1) + prow-2*P-1
            for qcol in 1:2*Q+1
                for pcol in 1:2*P+1
                    col = (qcol) * (2*P+1) + pcol-2*P-1
                    pfft = p[prow] - p[pcol]
                    qfft = q[qrow] - q[qcol]

                    C[row, col] = Af[q0 + pfft+1, p0 + qfft+1]
                
                end
            end
        end
    end

    return C
end
