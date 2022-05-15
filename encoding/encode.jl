
include("decode.jl")

function hamming_encode(bits::String, n::Int, k::Int)

    # using 7,4 hamming code

    encoded = Int.(zeros(n))
    parity_bits = [2^(i-1) for i in 1:n-k]
    
    j = 1; # scope issue
    
    for i in 1:length(encoded)
        if !(i in parity_bits) # fill bits in non parity positions
            encoded[i] = parse(Int, bits[j])
            j += 1
        end

    end

    for bit in parity_bits
        
        # check parity positions and xor bit positions using hamming encode algorithm

        parity = [] 
        for i in 1:length(encoded)

            if i & bit != 0 && i != bit
                append!(parity, encoded[i])
            end
        end
        encoded[bit] = reduce(xor, parity)
    end

    return reduce(*, string.(encoded))
end

function hamming_encode(bits, n::Int, k::Int)


    encoded = Int.(zeros(n))
    parity_bits = [2^(i-1) for i in 1:n-k]
    
    j = 1; # scope issue
    
    for i in 1:length(encoded)
        if !(i in parity_bits) # fill bits in non parity positions
            encoded[i] = bits[j]
            j += 1
        end

    end

    for bit in parity_bits
        
        # check parity positions and xor bit positions using hamming encode algorithm

        parity = [] 
        for i in 1:length(encoded)

            if i & bit != 0 && i != bit
                append!(parity, encoded[i])
            end
        end
        encoded[bit] = reduce(xor, parity)
    end

    return encoded
end

function majority_logic_encode(bits, n::Int)

    encoded = Int.(zeros(length(bits), n)); i = 1;

    for b in bits # fill each row with one bit
        encoded[i,:] .= b; i += 1
    end

    return encoded

end
    