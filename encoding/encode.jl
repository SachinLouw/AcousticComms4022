

function hamming_encode(bits, n, k)

    # using 7,4 hamming code

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

bits = [rand(0:1) for i in 1:4]

println(bits)
println(hamming_encode(bits, 7, 4))
println(reduce(xor, [bit for bit in hamming_encode(bits, 7, 4) if bit==1]))
