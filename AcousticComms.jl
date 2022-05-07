# acoustic communication script
# EEE4022F
# S Louw, 2022

include("modulation/FSK.jl")
include("modulation/PSK.jl")
include("encoding/encode.jl")
include("encoding/decode.jl")
include("plotting/plotting.jl")

using Plots, FFTW, SampledSignals, LibSndFile, PortAudio
# using Dash, DashCoreComponents, DashHtmlComponents
using FileIO: load, save, loadstreaming, savestreaming
plotly()

# constants
CAP_OFFSET = 3;
fs = 48000

# inline functions
rect(t) = (abs.(t).<=0.5)*1.0;
chirp(t, T, f₀, K) = rect.(t/T).*cos.(2*pi*(f₀*t .+ 0.5*K*t.^2))

function add_chirp(signal, fs, T = 0.05, f₀ = 10e3, B = 20e3)
    
    K = B/T #chirp rate
    Rmax = 10 #maximum range 
    c = 343
    t1 = 0:1/fs:(Rmax/c + T);

    #transmitted signal, adds zeros to cater for capacitive effect when recording
    v_tx = cat(zeros(Int(CAP_OFFSET*fs)), chirp(t1.-T/2, T, f₀, K), signal, dims = (1,1))

    return v_tx

end

function check_nulls()

    # create chirp 

    achirp = add_chirp(zeros(10), fs, 0.5)

    # start recording
    filename = "test/audio.wav"; time = (CAP_OFFSET + 5)
    recording = @async begin run(`ffmpeg -y -f alsa -i hw:CARD=PCH,DEV=0 -sample_rate $fs -t $time $filename`); end;
    # run(`ffmpeg -y -f alsa -i hw:CARD=PCH,DEV=0 -sample_rate $fs -t $time $filename`);

    # transmit
    PortAudioStream(0, 2) do stream
        write(stream, achirp);
    end

    wait(recording)

    # process recorded_signal
    recorded_signal = load(filename)[1][:,1] # recorded audio may be stereo (2-channels), only using one 
    
    synched_signal = match_filter(recorded_signal, fs)[1:length(achirp)]
        
    # get transfer function
    transfer_function = fft(synched_signal)./fft(achirp)

    # generate plots
    achirp_fig = plot_fd([], achirp, "Transmitted Chirp Spectrum")
    synched_signal_fig = plot_fd([], synched_signal, "Received Chirp Spectrum")
    transfer_function_fig  = plot_fd([], transfer_function, "Transfer function of Chirp Spectrum")

    return achirp_fig, synched_signal_fig, transfer_function_fig

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

function transmit(bit_stream::String, mod_technique::String, bits_to_send::Int = 1, pulses_per_second::Int = 1, fs::Int = 48000)

    if mod_technique == "FSK"
        modulated_signal = FSK.modulate(bit_stream, bits_to_send, pulses_per_second, fs)
        
    elseif mod_technique == "PSK"
        modulated_signal = PSK.modulate(bit_stream, bits_to_send, pulses_per_second, fs)

    elseif mod_technique == "QAM"
        # QAM(bit_stream, bits_to_send, pulses_per_second) 

    end

    modulated_signal_w_chirp = add_chirp(modulated_signal, fs, pulses_per_second, length(bit_stream))

    PortAudioStream(0, 2) do stream
        write(stream, modulated_signal_w_chirp);
    end

    # demodulated_sig = FSK.demodulate(modulated_signal[1], bits_to_send, pulses_per_second, modulated_signal[2], modulated_signal[3], fs)
    return modulated_signal_w_chirp

end

function receive(mod_technique::String, frequencies, spacing::Int = 2000, fs::Int = 48000, bits_to_send::Int = 1, pulses_per_second::Int = 1)

    filename = "test/audio.wav"; time = (CAP_OFFSET + 30/pulses_per_second)
    # recording = @async begin run(`ffmpeg -y -f alsa -i hw:CARD=PCH,DEV=0 -sample_rate $fs -t $time $filename`); end;
    # run(`ffmpeg -y -f alsa -i hw:CARD=PCH,DEV=0 -sample_rate $fs -t $time $filename`);

    # wait(recording)
    recorded_signal = load(filename)[1][:,1] # recorded audio may be stereo (2-channels), only using one 
    
    synched_signal = match_filter(recorded_signal, fs)

    if mod_technique == "FSK"
        demodulated_sig = FSK.demodulate(synched_signal, bits_to_send, pulses_per_second, spacing, frequencies, fs)
    elseif mod_technique == "PSK"
        demodulated_sig = PSK.demodulate(synched_signal, bits_to_send, pulses_per_second, frequencies, fs)
    elseif mod_technique == "QAM"
        # QAM(bit_stream, bits_to_send, pulses_per_second) 
    end

    return recorded_signal, synched_signal, demodulated_sig

end

function webpage()

    app = dash()

    app.layout = html_div() do
        html_h1("Acoustic Comms")
    end

    run_server(app, "0.0.0.0", 8080)

end
# start of program

println("----- Acoustic Communication thingy -------")

# achirp_fig, synched_signal_fig, transfer_function_fig = check_nulls()

# savefig(achirp_fig, "plots/Plot achirp_fig.html")
# savefig(synched_signal_fig, "plots/Plot synched_signal_fig.html")
# savefig(transfer_function_fig, "plots/Plot transfer_function_fig.html")

println("Transmitting or Receiving?")

choices = "1) Transmitting \n2) Receiving"

println(choices); transceive = readline()

# println("Choose modulation scheme"); 

# choices = "1) FSK \n2) PSK"

# println(choices); 

mod_technique = "FSK" #readline()

# println("Total number of bits?"); 
total_bits = 11 #readline()

# println("Number of bits to send within a pulse?");

bits_to_send = 1
bit_stream = "";
for i in 1:total_bits #==Generate random bit stream==#
    global bit_stream
    bit_stream = bit_stream * string(rand(0:1))
end

# @show bit_stream
n = 15; k = 11
encoded = hamming_encode(bit_stream, n, k)

# @show encoded

# decoded = hamming_decode(encoded);

# @show decoded

n_carriers = 2^bits_to_send; 
frequencies = Array{Int}(undef,0);
f1=2000                # initial carrier value
spacing = 2000         # spacing between frequencies

for n in 1:n_carriers
    append!(frequencies, f1+(n-1)*spacing)
end

if transceive == "1"
    modulated_signal_w_chirp = transmit(bit_stream, mod_technique)
    sleep(5)
    encoded_modulated_signal_w_chirp = transmit(encoded, mod_technique)

    modulated_signal_w_chirp_fig = plot_td([], modulated_signal_w_chirp, "Modulated transmitted signal")
    encoded_modulated_signal_w_chirp_fig = plot_td([], encoded_modulated_signal_w_chirp, "Modulated transmitted signal with encoding")

    savefig(modulated_signal_w_chirp_fig, "plots/Plot modulated_signal_w_chirp.html")
    savefig(encoded_modulated_signal_w_chirp_fig, "plots/Plot encoded_modulated_signal_w_chirp.html")


    @show bit_stream, encoded

elseif transceive == "2"
    recorded_signal, synched_signal, (msg, t_array, signal_array) = receive(mod_technique, frequencies);
    syndrome = hamming_syndrome(parse.(Int, split(msg[1:15], "")))
    decoded_msg = hamming_decode(msg[1:15])
    @show syndrome
    @show decoded_msg

    synched_signal_fig = plot_td([], synched_signal, "Received signal synchronised")
    savefig(synched_signal_fig, "plots/Plot 1.html")

    pl = 2

    # for arr in 1:length(signal_array)
    #     savefig(plot_td(t_array[arr], signal_array[arr]), "Plot " * string(pl))
    #     pl += 1
    # end 
end

#== TODO: fix order of operation
1   add plots to each function
        create plot objects inside each function    DONE
        add variable that will plot only if needed
2   fix PSK for n case 
3   fix QAM

==#

println("------- End -------")