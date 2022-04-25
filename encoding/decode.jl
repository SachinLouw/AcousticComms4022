


function hamming_syndrome(bits)
    # https://www.youtube.com/watch?v=b3NxrZOu_CE&ab_channel=3Blue1Brown

    # Reduce by xor all indices of active bit
    return reduce(xor, [i for (i, bit) in enumerate(bits) if bit==1])

end
