module PSK

    using Trapz

    function modulate(bit_stream, carriers, bits_to_send::Int, pulses_per_second::Int, fs::Int = 48000)
          
        modulated_signal = Array{Float64}(undef,0)
        start(i) = 1+Int(floor(fs/(bits_to_send*pulses_per_second)*(i-1)))
        stop(i) = Int(floor(fs/(bits_to_send*pulses_per_second)*(i)))
        
        for i in 1:bits_to_send:length(bit_stream)
            
            for j in 1:2^bits_to_send
                seq = parse(Int, join(bit_stream[i:i+bits_to_send-1]); base=2) # row_bits bit seq casted as base 2 int, join method converts array chunk to string
                if seq == j - 1
                    append!(modulated_signal, carriers[j][start(i):stop(i)]);
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
        ensembles = []

        t_array, signal_array = [], []
        while msg_index <= col_bits*row_bits
            t_slice = t[start:stop]
            y1_slice = signal[start:stop]

            for ϕ in phases
                y = y1_slice .* cos.(2*π*f.*t_slice .- ϕ)
                append!(ensembles, trapz(t_slice, y))
            end
            # println(ensembles)
            
            index = findmax(ensembles)[2]

            bit = reverse(digits(index - 1, base=2, pad = bits_to_send))
            ensembles = []
            msg[msg_index:msg_index + bits_to_send - 1] = bit[end-(bits_to_send - 1):end]
            msg_index = msg_index + bits_to_send
            
            start += slice
            stop += slice
                    
        end

        @show transpose(msg);
        return msg, t_array, signal_array

    end

end





















#== QPSK 
        b_odd = bit_stream[1:2:end]; b_even = bit_stream[2:2:end]
        xo = []; xe = []


        for j in 1:Int(length(bit_stream)/2)
            
            if b_odd[j] == 1
                xo = append!(xo, ones(1,Int(fs/pulses_per_second))) 
            else
                xo = append!(xo, -1*ones(1,Int(fs/pulses_per_second)))
            end
            if b_even[j] == 1
                xe = append!(xe, ones(1,Int(fs/pulses_per_second)))
            else
                xe = append!(xe, -1*ones(1,Int(fs/pulses_per_second)))
            end

        end

        f=500;# f signal frequency

        x1=cos.(2*pi*f.*t); y1=xo.*x1;
        x2=sin.(2*pi*f.*t); y2=xe.*x2;
        y = y1 .+ y2

        # TODO: iterate over array, split into subarrays of length/bits_to_send
        #       modulated_signal = sum of subarrays
        
    ==#
