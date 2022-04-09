module PSK


    function modulate(bit_stream::String, bits_to_send::Int, pulses_per_second::Int, fs::Int = 48000)

        Δt = 1/fs; t = 0:Δt:(length(bit_stream)/bits_to_send-Δt)/pulses_per_second;
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
        a=1;# Amplitude scale for carrier signal

        f1=1000;# carrier signal frequency

        x1=a.*cos.(2*pi*f1.*t); y1=xo.*x1;
        x2=a.*sin.(2*pi*f1.*t); y2=xe.*x2;
        y = y1 .+ y2

        return y

    end

    function demodulate(signal, bits_to_send::Int, pulses_per_second::Int, spacing::Int, frequencies, fs::Int = 48000)

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
            yx = y1_slice .* x1[1:length(y1_slice)]
        #     plot_td(t_slice, yx, "Output * carrier","Time","Amplitude");
            avg1 = trapz(t_slice, yx)
            
            if avg1 > 0
                append!(msg,1)
            else
                append!(msg,0)
            end 
            
            yx = y1_slice .* x2[1:length(y1_slice)]
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

        @show bit_stream;
        @show msg;
    end

end