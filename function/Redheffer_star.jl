using LinearAlgebra

function dict_to_matrix(S_dict)

    return [S_dict["S11"] S_dict["S12"]; S_dict["S21"] S_dict["S22"]]
end

function RedhefferStar(SA, SB)

    SA_11 = SA["S11"]
    SA_12 = SA["S12"]
    SA_21 = SA["S21"]
    SA_22 = SA["S22"]
    SB_11 = SB["S11"]
    SB_12 = SB["S12"]
    SB_21 = SB["S21"]
    SB_22 = SB["S22"]
    
    N = size(SA_11, 1) 
    In = Matrix{Float64}(I, N, N)

    D = In - SB_11 * SA_22
    F = In - SA_22 * SB_11

    SAB_11 = SA_11 + SA_12 / D * SB_11 * SA_21
    SAB_12 = SA_12 / D * SB_12
    SAB_21 = SB_21 / F * SA_21
    SAB_22 = SB_22 + SB_21 / F * SA_22 * SB_12

    SAB = [SAB_11 SAB_12; SAB_21 SAB_22]
    SAB_dict = Dict("S11" => SAB_11, "S22" => SAB_22, "S12" => SAB_12, "S21" => SAB_21)

    return SAB, SAB_dict
end