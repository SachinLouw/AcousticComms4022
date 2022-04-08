# acoustic communication script
# EEE4022F
# S Louw, 2022

include("modulation/FSK.jl")
using PortAudio


# inline functions
rect(t) = (abs.(t).<=0.5)*1.0;
chirp(t, T, f₀, K) = rect.(t/T).*cos.(2*pi*(f₀*t .+ 0.5*K*t.^2))

function add_chirp(signal, fs, p, len_bs)
    
    #define centre f, BW and T
    f₀ = 10e3; B = 10e3; T = 0.05;

    K = B/T #chirp rate
    Rmax = 10 #maximum range 
    c = 343
    t1 = 0:1/fs:(Rmax/c + T);

    #transmitted signal
    v_tx = cat(chirp(t1.-T/2, T, f₀, K), signal, dims=(1,1))
    # plot_td(t, v_tx, "Modulated output with chirp","Time","Amplitude");
    # plot_fd(t, v_tx, "Modulated output with chirp","Frequency","Amplitude");

    cap_offset = 3;
    v_tx = cat(zeros(Int(cap_offset*fs)), v_tx, dims = (1,1))

    return v_tx
end

function match_filter(recorded_signal)

    f₀ = 10e3; B = 10e3; T = 0.05;

    K = B/T #chirp rate
    Rmax = 10 #maximum range 
    c = 343
    t1 = 0:1/fs:(Rmax/c + T);

    chirp_length = length(chirp(t1.-T/2, T, f₀, K))
    padded_chirp = cat(chirp(t1.-T/2, T, f₀, K), zeros(length(recorded_signal) - chirp_length),dims=(1,1))

    H = conj(fft(padded_chirp))
    V = fft(recorded_signal)
    y4 = real(ifft(H .* V))

    match_offset = Int(floor(findmax(y4)[2] + chirp_length))
    
    return recorded_signal[match_offset:end]

end

function transmit(bit_stream::String, mod_technique::String, fs::Float64 = 48000.0, bits_to_send::Int = 1, pulses_per_second::Int = 1)

   if mod_technique == "FSK"
        modulated_signal = FSK.modulate(bit_stream, bits_to_send, pulses_per_second, fs)
        
    elseif mod_technique == "PSK"
        # PSK(bit_stream, bits_to_send, pulses_per_second)
    elseif mod_technique == "QAM"
        # QAM(bit_stream, bits_to_send, pulses_per_second) 
    end

    # modulated_signal_w_chirp = add_chirp(modulated_signal, fs, pulses_per_second, length(bit_stream))

    # PortAudioStream(0, 2) do stream
    #     write(stream, modulated_signal_w_chirp);
    # end
    demodulated_sig = FSK.demodulate(modulated_signal[1], bits_to_send, pulses_per_second, modulated_signal[2], modulated_signal[3], fs)
end

function receive(synched_signal, bits_to_send, pulses_per_second, fs)

    if mod_technique == "FSK"
        demodulated_sig = FSK.demodulate(synched_signal, bits_to_send, pulses_per_second, fs)
        
    elseif mod_technique == "PSK"
        # PSK(bit_stream, bits_to_send, pulses_per_second)
    elseif mod_technique == "QAM"
        # QAM(bit_stream, bits_to_send, pulses_per_second) 
    end
end

# start of program

println("----- Acoustic Communication thingy -------")

bits_to_send = 1
bit_stream = "";
for i in 1:10 #==Generate random bit stream==#
    global bit_stream
    bit_stream = bit_stream * string(rand(0:1))
end

@show bit_stream
transmit(bit_stream, "FSK")


println("----- End -------")