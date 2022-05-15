


function hamming_syndrome(bits)
    # https://www.youtube.com/watch?v=b3NxrZOu_CE&ab_channel=3Blue1Brown

    # Reduce by xor all indices of active bit
    return reduce(xor, [i for (i, bit) in enumerate(bits) if bit==1])

end

function hamming_decode(bits)

    decoded = []
    # @show bits

    for i in 3:length(bits) # skips power of 2 indexes, we know 1=2^0 and 2 = 2^1

        if !(isinteger(log2(i)))
            # append!(decoded, bits[i])
            decoded = [decoded; bits[i]]
        end

    end
    # @show decoded
    return Vector{Int}(decoded)

end


function majority_logic_decode(bits)

    decoded = zeros(length(bits[:,1])); i = 1;

    for row in eachrow(bits)
        decoded[i] = sum(row)/length(row); i += 1
    end
    
    return round.(decoded)

end

function hamming_decode(bits::String)

    decoded = "" #[]

    for i in 3:length(bits) # skips power of 2 indexes, we know 1=2^0 and 2 = 2^1

        if !(isinteger(log2(i)))
            # append!(decoded, bits[i])
            decoded = decoded * bits[i]
        end

    end

    return decoded

end


    