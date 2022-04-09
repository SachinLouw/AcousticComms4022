# acoustic communication script
# EEE4022F
# S Louw, 2022

include("modulation/FSK.jl")

using Plots, FFTW, SampledSignals, LibSndFile, PortAudio, .Threads
using FileIO: load, save, loadstreaming, savestreaming

# inline functions
rect(t) = (abs.(t).<=0.5)*1.0;
chirp(t, T, f₀, K) = rect.(t/T).*cos.(2*pi*(f₀*t .+ 0.5*K*t.^2))
CAP_OFFSET = 3;

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

    v_tx = cat(zeros(Int(CAP_OFFSET*fs)), v_tx, dims = (1,1))

    return v_tx
end

function match_filter(recorded_signal, fs)

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

function transmit(bit_stream::String, mod_technique::String, fs::Int = 48000, bits_to_send::Int = 1, pulses_per_second::Int = 1)

    if mod_technique == "FSK"
        modulated_signal = FSK.modulate(bit_stream, bits_to_send, pulses_per_second, fs)
        
    elseif mod_technique == "PSK"
        # PSK(bit_stream, bits_to_send, pulses_per_second)
    elseif mod_technique == "QAM"
        # QAM(bit_stream, bits_to_send, pulses_per_second) 
    end

    modulated_signal_w_chirp = add_chirp(modulated_signal[1], fs, pulses_per_second, length(bit_stream))

    PortAudioStream(0, 2) do stream
        write(stream, modulated_signal_w_chirp);
    end
    # demodulated_sig = FSK.demodulate(modulated_signal[1], bits_to_send, pulses_per_second, modulated_signal[2], modulated_signal[3], fs)
end

function receive(mod_technique::String, frequencies::Array{Int} = [], spacing::Int = 2000, fs::Int = 48000, bits_to_send::Int = 1, pulses_per_second::Int = 1)

    filename = "test/audio.wav"; time = (CAP_OFFSET + 40/pulses_per_second)
    recording = @async begin run(`ffmpeg -y -f alsa -i hw:CARD=PCH,DEV=0 -sample_rate $fs -t $time $filename`); end;

    wait(recording)
    recorded_signal = load(filename)[1][:,1] # recorded audio may be stereo (2-channels), only using one 
    
    synched_signal = match_filter(recorded_signal, fs)

    if mod_technique == "FSK"
        demodulated_sig = FSK.demodulate(synched_signal, bits_to_send, pulses_per_second, spacing, frequencies, fs)
        
    elseif mod_technique == "PSK"
        # PSK(bit_stream, bits_to_send, pulses_per_second)
    elseif mod_technique == "QAM"
        # QAM(bit_stream, bits_to_send, pulses_per_second) 
    end
end

# start of program

println("----- Acoustic Communication thingy -------")

# println("Choose modulation scheme"); 

# choices = "1) FSK \n2) PSK"

# println(choices); 

# modulation = readline()

println("Transmitting or Receiving?")

choices = "1) Transmitting \n2) Receiving"

println(choices); t_r = readline()

# println("Total number of bits?"); total_bits = readline()

# println("Number of bits to send within a pulse?"); total_bits = readline()
bits_to_send = 1
bit_stream = "";
for i in 1:20 #==Generate random bit stream==#
    global bit_stream
    bit_stream = bit_stream * string(rand(0:1))
end

n_carriers = 2^bits_to_send; 
frequencies = Array{Int}(undef,0);
f1=2000                # initial carrier value
spacing = 2000         # spacing between frequencies

for n in 1:n_carriers
    append!(frequencies, f1+(n-1)*spacing)
end

if t_r == "1"
    transmit(bit_stream, "FSK")
    @show bit_stream

elseif t_r == "2"
    receive("FSK", frequencies)
end


#== TODO: fix order of operation
1   transmit DONE
2   add_chirp DONE
3   record DONE
4   sync (apply match_filter and feed into demodulator) DONE
5   demodulate DONE
==#
println("----- End -------")