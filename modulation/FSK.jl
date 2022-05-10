module FSK

    export modulate, demodulate

    using Plots, FFTW

    function modulate(bit_stream, bits_to_send::Int, pulses_per_second::Int, fs::Int = 48000)

        Δt = 1/fs; t = 0:Δt:(length(bit_stream)/bits_to_send-Δt)/pulses_per_second;

        #TODO: move frequencies outside so demodulate can also use 
        n_carriers = 2^bits_to_send; 
        carriers = []; frequencies = Array{Int}(undef,0);
        f1=2000                # initial carrier value
        spacing = 2000         # spacing between frequencies

        for n in 1:n_carriers
            append!(carriers, [cos.(2*π*(f1+(n-1)*spacing).*t)])
            append!(frequencies, f1+(n-1)*spacing)
        end
        
        modulated_signal = Array{Float64}(undef,0)
        start(i) = 1+Int(floor(fs/(bits_to_send*pulses_per_second)*(i-1)))
        stop(i) = Int(floor(fs/(bits_to_send*pulses_per_second)*(i)))
        
        for i in 1:bits_to_send:length(bit_stream)
            
            for j in 1:n_carriers
                seq = parse(Int, join(bit_stream[i:i+bits_to_send-1]); base=2) # n bit seq casted as base 2 int, join method converts array chunk to string
                if seq == j - 1
                    append!(modulated_signal, carriers[j][start(i):stop(i)]);
                    break
                end
            end
            
        end

        return modulated_signal

    end

    function demodulate(signal, bits_to_send::Int, p::Int, spacing::Int, frequencies, ratio::Int, n::Int, fs::Int = 48000)

        start = 1
        slice = Int(floor(fs/(bits_to_send * p)))
        end_ = Int(floor(start + slice))
        
        Δt = 1/fs
        t = 0:Δt:(length(signal) - 1)* Δt
        
        msg = zeros(n, Int(ratio)); msg_index = 1;

        # t_array, signal_array = [], []
        
        while msg_index <= ratio*n # using < inst of <= to remove extra bit
             
            # global start; global end_; local stop; global msg;
            t_slice = t[start:end_]
            y1_slice = signal[start:end_]
        #     plot_td(t_slice, y1_slice, "Modulated output at receiver","Time","Amplitude");
        #     plot_fd(t_slice, y1_slice, "Modulated output at receiver","Time","Amplitude");
            # append!(t_array); append!(signal_array);

            Y1 = abs.(fft(y1_slice))
            
            N = length(y1_slice);
            Δf = 1/(N*Δt)  # spacing in frequency domain
        
            stop = Int(floor(length(Y1)/2))
        
            #create array of freq values stored in f_axis. First element maps to 0Hz
            f = (0:N-1)*Δf;    
            ind = Int(floor(findmax(Y1[1:stop])[2]))
            
            for j in frequencies
                if abs(f[ind] - j) <= spacing/2
                    # @show j
                    # @show indexin(j, frequencies)[1]
                    bit = reverse(digits(indexin(j, frequencies)[1] - 1, base=2, pad = bits_to_send))

                    msg[msg_index:msg_index + bits_to_send - 1] = bit[end-(bits_to_send - 1):end]
                    # @show msg[msg_index:msg_index + bits_to_send - 1]
                    # @show msg_index
                    msg_index = msg_index + bits_to_send
                    println(msg_index)
                    # break

                # elseif indexin(j, frequencies)[1] >= length(frequencies)
                #     msg = msg * "*" # could not determine bit or end of bits
                end
                
            end
        #     @show (f[ind] - f1)
            
            start += slice
            end_ += slice
            
        end
        
        @show transpose(msg);
        return msg, t_array, signal_array

    end


    # process input as strings

    function modulate(bit_stream::String, bits_to_send::Int, pulses_per_second::Int, fs::Int = 48000)

        Δt = 1/fs; t = 0:Δt:(length(bit_stream)/bits_to_send-Δt)/pulses_per_second;

        #TODO: move frequencies outside so demodulate can also use 
        n_carriers = 2^bits_to_send; 
        carriers = []; frequencies = Array{Int}(undef,0);
        f1=2000                # initial carrier value
        spacing = 2000         # spacing between frequencies

        for n in 1:n_carriers
            append!(carriers, [cos.(2*π*(f1+(n-1)*spacing).*t)])
            append!(frequencies, f1+(n-1)*spacing)
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

        return modulated_signal

    end


    # function demodulate(signal, bits_to_send::Int, p::Int, spacing::Int, frequencies, fs::Int = 48000)

    #     start = 1
    #     slice = Int(floor(fs/(bits_to_send * p)))
    #     end_ = Int(floor(start + slice))
        
    #     Δt = 1/fs
    #     t = 0:Δt:(length(signal) - 1)* Δt
        
    #     msg = "";

    #     t_array, signal_array = [], []
        
    #     while end_ < length(signal) # using < inst of <= to remove extra bit
             
    #         # global start; global end_; local stop; global msg;
    #         t_slice = t[start:end_]
    #         y1_slice = signal[start:end_]
    #     #     plot_td(t_slice, y1_slice, "Modulated output at receiver","Time","Amplitude");
    #     #     plot_fd(t_slice, y1_slice, "Modulated output at receiver","Time","Amplitude");
    #         append!(t_array); append!(signal_array);

    #         Y1 = abs.(fft(y1_slice))
            
    #         N = length(y1_slice);
    #         Δf = 1/(N*Δt)  # spacing in frequency domain
        
    #         stop = Int(floor(length(Y1)/2))
        
    #         #create array of freq values stored in f_axis. First element maps to 0Hz
    #         f = (0:N-1)*Δf;    
    #         ind = Int(floor(findmax(Y1[1:stop])[2]))
            
    #         for j in frequencies
    #             if abs(f[ind] - j) <= spacing/2
    #     #             @show j
    #     #             @show indexin(j, frequencies)[1]
    #                 bit = bitstring(indexin(j, frequencies)[1] - 1)
    #                 msg = msg * bit[end-(bits_to_send - 1):end]; break

    #             elseif indexin(j, frequencies)[1] >= length(frequencies)
    #                 msg = msg * "*" # could not determine bit or end of bits
    #             end
                
    #         end
    #     #     @show (f[ind] - f1)
            
    #         start += slice
    #         end_ += slice
            
    #     end
        
    #     @show msg;
    #     return msg, t_array, signal_array

    # end

end # module