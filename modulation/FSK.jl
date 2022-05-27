module FSK

    export modulate, demodulate

    include("../plotting/plotting.jl")
    using FFTW

    function modulate(bit_stream, carriers, bits_to_send::Int, pulses_per_second::Int, fs::Int = 48000)
        
        modulated_signal = Array{Float64}(undef,0)
        start(i) = 1+Int(floor(fs/(bits_to_send*pulses_per_second)*(i-1)))
        stop(i) = Int(floor(fs/(bits_to_send*pulses_per_second)*(i)))
        
        for i in 1:bits_to_send:length(bit_stream)
            
            for j in 1:2^bits_to_send
                seq = parse(Int, join(bit_stream[i:i+bits_to_send-1]); base=2) # n bit seq casted as base 2 int, join method converts array chunk to string
                if seq == j - 1
                    append!(modulated_signal, carriers[j][start(i):stop(i)]);
                    break
                end
            end
            
        end

        return modulated_signal

    end

    function demodulate(signal, bits_to_send::Int, pulses_per_second::Int, spacing::Int, frequencies, ratio::Int, n::Int, fs::Int = 48000)

        
        start = 1
        slice = Int(floor(fs/(bits_to_send * pulses_per_second)))
        end_ = Int(floor(start + slice))
        
        Δt = 1/fs
        t = 0:Δt:(length(signal) - 1)* Δt
        
        msg = zeros(n, ratio); msg_index = 1;

        t_array, signal_array = [], []
        # @show length(signal)
        while msg_index <= ratio*n

            # global start; global end_; local stop; global msg;
            t_slice = t[start:end_]
            y1_slice = signal[start:end_]
            # display(plot_td(t_slice, y1_slice, "Modulated output at receiver","Time","Amplitude"))
            # display(plot_fd(t_slice, y1_slice, "Modulated output Spectrum $(t[start])s","Time","Amplitude"))            
            # append!(t_array); append!(signal_array);

            Y1 = abs.(fft(y1_slice))
            stop = Int(floor(length(Y1)/2))
            N = length(y1_slice);
            Δf = 1/(N*Δt)  # spacing in frequency domain
            #create array of freq values stored in f_axis. First element maps to 0Hz
            f = (0:N-1)*Δf;    

            B = 500 # filter bandwidth in Hz
            rect(t) = (abs.(t).<=0.5)*1.0;
            H = 1 .- (rect(f/(2π*B)) + rect( (f .- 1/Δt)/(2π*B) ));

            # HPF to reject DC comp during recording
            Y1 = conj(Y1).*H

            ind = Int(floor(findmax(Y1[1:stop])[2]))
            # @show f[ind]
            for j in frequencies
                if abs(f[ind] - j) <= spacing/2
                    bit = reverse(digits(indexin(j, frequencies)[1] - 1, base=2, pad = bits_to_send))

                    msg[msg_index:msg_index + bits_to_send - 1] = bit[end-(bits_to_send - 1):end]
                    msg_index = msg_index + bits_to_send
                    break

                elseif j == frequencies[end] 
                    bit = reverse(digits(0, base=2, pad = bits_to_send))

                    msg[msg_index:msg_index + bits_to_send - 1] = bit[end-(bits_to_send - 1):end]
                    msg_index = msg_index + bits_to_send
                end
                
            end
        #     @show (f[ind] - f1)
            
            start += slice
            end_ += slice
            
        end
        
        # @show (msg);
        # @show transpose(msg);
        return msg, t_array, signal_array

    end

end