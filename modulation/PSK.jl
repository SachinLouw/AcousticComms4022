module PSK

    using Trapz

    function modulate(bit_stream::String, bits_to_send::Int = 2, pulses_per_second::Int = 1, fs::Int = 48000)

        #==
        Δt = 1/fs; t = 0:Δt:(length(bit_stream)/bits_to_send-Δt)/pulses_per_second;
        
        # MPSK
      
        n_carriers = 2^bits_to_send; 
        carriers = []; phase = Array{Float64}(undef,0);
        
        f1=2000                # initial carrier value


        a = 1
        for n in 1:n_carriers

            ϕ = 2*π*(n-1)/n_carriers
            append!(carriers, [cos.(2*π*f1).*t - ϕ])
            append!(phase, f1+(n-1)*spacing)
        end
        
        modulated_signal = Array{Float64}(undef,0)
        start(i) = 1+Int(floor(fs/(bits_to_send*pulses_per_second)*(i-1)))
        stop(i) = Int(floor(fs/(bits_to_send*pulses_per_second)*(i)))
        
        for i in 1:bits_to_send:length(bit_stream)
            
            for j in 1:n_carriers
                seq = parse(Int, bit_stream[i:i+bits_to_send-1]; base=2) # n bit seq casted as base 2 int
                if seq == j - 1
                    append!(modulated_signal, carriers[j][start(i):stop(i)]);
                    break
                end
            end
            
        end        
        
        # TODO: iterate over array, split into subarrays of length/bits_to_send
        #       modulated_signal = sum of subarrays
        
        ==#
        
        # QPSK
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

        carrier=500;# carrier signal frequency

        x1=cos.(2*pi*carrier.*t); y1=xo.*x1;
        x2=sin.(2*pi*carrier.*t); y2=xe.*x2;
        y = y1 .+ y2

        return y

    end

    function demodulate(signal, bits_to_send::Int, pulses_per_second::Int, carrier::Int = 2000, fs::Int = 48000)

        start = 1
        slice = Int(floor(fs/pulses_per_second))
        end_ = Int(floor(start + slice))

        Δt = 1/fs
        t = 0:Δt:(length(signal) - 1)* Δt

        msg = []
        x = []; y = []

        while end_ <= length(signal)
            t_slice = t[start:end_]
            y1_slice = signal[start:end_]
        #     plot_td(t_slice, y1_slice, "Modulated output at receiver","Time","Amplitude");
            yx = y1_slice .* cos.(2*pi*carrier.*t_slice)
        #     plot_td(t_slice, yx, "Output * carrier","Time","Amplitude");
            avg1 = trapz(t_slice, yx)
            
            if avg1 > 0
                append!(msg,1)
            else
                append!(msg,0)
            end 
            
            yx = y1_slice .* sin.(2*pi*carrier.*t_slice)
        #     plot_td(t_slice, yx, "Output * quad carrier","Time","Amplitude");
            avg = trapz(t_slice, yx)
            append!(x, avg1); append!(y, avg)
        #     @show atan(avg, avg1)
            if avg > 0
                append!(msg,1)
            else
                append!(msg,0)
            end 
            start += slice
            end_ += slice
                    
        end

        @show msg;
    end

end