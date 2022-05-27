module QAM

    using Trapz

    function modulate(bit_stream, carriers, bits_to_send::Int, pulses_per_second::Int, fs::Int = 48000)
          
        modulated_signal = Array{Float64}(undef,0)
        start(i) = 1+Int(floor(fs/(bits_to_send*pulses_per_second)*(i-1)))
        stop(i) = Int(floor(fs/(bits_to_send*pulses_per_second)*(i)))
        
        for i in 1:bits_to_send:length(bit_stream)
            
            for j in 1:Int(2^bits_to_send/2)
                seq = parse(Int, join(bit_stream[i:i+bits_to_send-1]); base=2) # row_bits bit seq casted as base 2 int, join method converts array chunk to string
                if seq == j - 1
                    append!(modulated_signal, carriers[j][start(i):stop(i)] .+ carriers[j*2][start(i):stop(i)]);
                    break
                end
            end
            
        end     
        
        return modulated_signal
    
    end     


    function demodulate(signal, bits_to_send::Int, pulses_per_second::Int, f::Int, phases, col_bits::Int, row_bits::Int, fs::Int = 48000)

        start = 1
        slice = Int(floor(fs/(bits_to_send * pulses_per_second)))
        stop = Int(floor(start + slice))

        Δt = 1/fs
        t = 0:Δt:(length(signal) - 1)* Δt

        msg = zeros(row_bits, col_bits); msg_index = 1;
        i_ensembles = [];   q_ensembles = []

        t_array, signal_array = [], []
        while msg_index <= col_bits*row_bits
            t_slice = t[start:stop]
            y1_slice = signal[start:stop]

            for a in amplitudes
                i = y1_slice .* a*cos.(2*π*f.*t_slice)
                q = y1_slice .* a*sin.(2*π*f.*t_slice)

                append!(i_ensembles, trapz(t_slice, i))
                append!(q_ensembles, trapz(t_slice, q))
            end
            # println(ensembles)
            
            index = findmax(i_ensembles)[2]

            bit = reverse(digits(index - 1, base=2, pad = bits_to_send))
            i_ensembles = []
            msg[msg_index:msg_index + bits_to_send - 1] = bit[end-(bits_to_send - 1):end]
            
            index = 2*findmax(q_ensembles)[2]

            bit = reverse(digits(index - 1, base=2, pad = bits_to_send))
            q_ensembles = []
            msg[msg_index:msg_index + bits_to_send - 1] = bit[end-(bits_to_send - 1):end]

            msg_index = msg_index + bits_to_send
            
            start += slice
            stop += slice
                    
        end

        @show transpose(msg);
        return msg, t_array, signal_array

    end

end
