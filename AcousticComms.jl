# acoustic communication script
# EEE4022F
# S Louw, 2022

include("modulation/FSK.jl")
include("modulation/PSK.jl")
include("modulation/ASK.jl")
include("modulation/QAM.jl")
include("encoding/encode.jl")
include("encoding/decode.jl")
include("plotting/plotting.jl")

using Plots, FFTW, WAV, CSV, DataFrames, SampledSignals, LibSndFile
using FileIO: load, save, loadstreaming, savestreaming
plotly()

# constants
CAP_OFFSET = 3; # delay when recording for capacitive effect
fs = 48000; # sampling rate

# inline functions
rect(t) = (abs.(t).<=0.5)*1.0;
chirp(t, T, f₀, K) = rect.(t/T).*cos.(2*pi*(f₀*t .+ 0.5*K*t.^2))

function add_chirp(signal, T = 0.5, f₀ = 10e3, B = 20e3)
    
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

    achirp = add_chirp(zeros(fs*3), 0.5) # three second pause to register echoing

    # start recording
    filename = "test/nulls.wav"; time = (CAP_OFFSET + 20)
    recording = @async begin run(`arecord -r $fs -d $time -t wav $filename`) end;

    # transmit
    sleep(3)
    wavplay(achirp, fs)


    wait(recording)

    # process recorded_signal
    recorded_signal = load(filename)[1][:,1] # recorded audio may be stereo (2-channels), only using one 

    synched_signal = matched_filter(recorded_signal, fs, 0.5, 10e3, 20e3, true)[1:length(achirp)]

    # get transfer function
    transfer_function = fft(synched_signal)./fft(achirp)

    # generate plots
    achirp_fig = plot_fd([], achirp, "Transmitted Chirp Spectrum")
    synched_signal_fig = plot_fd([], synched_signal, "Received Chirp Spectrum")
    transfer_function_fig  = plot_fd([], (ifft(transfer_function)), "Transfer function of Room")

    return achirp_fig, synched_signal_fig, transfer_function_fig

end

function matched_filter(recorded_signal, fs, T = 0.5, f₀ = 10e3, B = 20e3, incl_chirp = false)

    K = B/T #chirp rate
    Rmax = 10 #maximum range 
    c = 343
    t1 = 0:1/fs:(Rmax/c + T);

    chirp_length = length(chirp(t1.-T/2, T, f₀, K))
    padded_chirp = cat(chirp(t1.-T/2, T, f₀, K), zeros(length(recorded_signal) - chirp_length),dims=(1,1))

    H = conj(fft(padded_chirp))
    V = fft(recorded_signal)
    y4 = real(ifft(H .* V))

    matched_index = Int(floor(findmax(y4)[2])) #****************
    
    if incl_chirp
        return recorded_signal[matched_index:end]
    else
        return recorded_signal[matched_index + chirp_length:end]
    end

end

function transmit(bit_stream, mod_technique::String, carriers, pulses_per_second::Int = 1, bits_to_send::Int = 1)

    if mod_technique == "FSK"
        modulated_signal = FSK.modulate(bit_stream, carriers, bits_to_send, pulses_per_second, fs)
        
    elseif mod_technique == "PSK"
        modulated_signal = PSK.modulate(bit_stream, carriers, bits_to_send, pulses_per_second, fs)

    elseif mod_technique == "ASK"
        modulated_signal = ASK.modulate(bit_stream, carriers, bits_to_send, pulses_per_second, fs)

    elseif mod_technique == "QAM"
        modulated_signal = QAM.modulate(bit_stream, carriers, bits_to_send, pulses_per_second, fs)

    end

    modulated_signal_w_chirp = add_chirp(modulated_signal, 1)

    wavplay(modulated_signal_w_chirp, fs)

    return modulated_signal_w_chirp

end

function receive(mod_technique::String, spacing::Int, frequencies, phases, pulses_per_second, ratio::Int, n::Int, bits_to_send::Int, recording = true, filename = "test/audio.wav")


    if recording
        time = 2*Int(CAP_OFFSET + ceil(ratio*n/pulses_per_second))
        run(`arecord -r $fs -d $time -t wav $filename`)
    end

    recorded_signal = load(filename)[1][:,1] # recorded audio may be stereo (2-channels), only using one 
    
    synched_signal = matched_filter(recorded_signal, fs, 1)

    if mod_technique == "FSK"
        demodulated_sig = FSK.demodulate(synched_signal, bits_to_send, pulses_per_second, spacing, frequencies, ratio, n, fs)
    elseif mod_technique == "PSK"
        demodulated_sig = PSK.demodulate(synched_signal, bits_to_send, pulses_per_second, frequencies[1], phases, ratio, n, fs)
    elseif mod_technique == "ASK"

    elseif mod_technique == "QAM"

    end

    return recorded_signal, synched_signal, demodulated_sig

end

function automate(mod_technique, bits_to_send, pulses_per_second, total_bits, encoding, f1, spacing, row_bits, filename)
    
    if encoding == "Hamming"
        col_bits = Int(total_bits/row_bits)

        bit_stream = Int.(zeros(col_bits,row_bits))
        for i in 1:total_bits
            bit_stream[i] = rand(0:1)
        end

        encoded = zeros(col_bits, row_bits)

        for i in 1:Int(total_bits/row_bits)
            encoded[i,1:n] = hamming_encode(bit_stream[i,1:row_bits], n, k)
        end
    elseif encoding == "Majority" 
        bit_stream = [rand(0:1) for i in 1:total_bits]
        encoded = majority_logic_encode(bit_stream, row_bits); 
    end


    col_bits = length(encoded[:,1])

    Δt = 1/fs; t = 0:Δt:(length(encoded)/bits_to_send-Δt)/pulses_per_second;

    n_carriers = 2^bits_to_send; 
    carriers = []; 
    frequencies = Array{Int}(undef,0); 
    phases = Array{Float64}(undef,0);
    amplitudes = Array{Float64}(undef,0);

    if mod_technique == "FSK"

        for n in 1:n_carriers
            append!(carriers, [cos.(2*π*(f1+(n-1)*spacing).*t)])
            append!(frequencies, f1+(n-1)*spacing)
        end

    elseif mod_technique == "PSK"
        
        append!(frequencies, f1)
        for n in 1:n_carriers

            ϕ = 2*π*(n-1)/n_carriers
            append!(carriers, [cos.(2*π*f1.*t .- ϕ)])
            append!(phases, ϕ)
        end

    elseif mod_technique == "ASK"
        
        append!(frequencies, f1)
        for n in 1:n_carriers

            a = (n-1)/n_carriers
            append!(carriers, [a*cos.(2*π*f1.*t)])
            append!(amplitudes, a)
        end
    
    elseif mod_technique == "QAM"

        append!(frequencies, f1)
        for n in 1:n_carriers

            a = n/n_carriers
            append!(carriers, [a*cos.(2*π*f1.*t)])
            append!(carriers, [a*sin.(2*π*f1.*t)])
            append!(amplitudes, a)
        end

    end

    time = 2*Int(CAP_OFFSET + ceil(col_bits*row_bits/pulses_per_second))
    recording = @async begin run(`arecord -r $fs -d $time -t wav test/$filename`) end
    sleep(5)
    encoded_modulated_signal_w_chirp = transmit(transpose(encoded), mod_technique, carriers, pulses_per_second, bits_to_send)
    wait(recording)
    recorded_signal, synched_signal, (msg, t_array, signal_array) = receive(mod_technique, spacing, frequencies, phases, pulses_per_second, col_bits, row_bits, bits_to_send, false, "test/$filename");
    decoded = Int.(majority_logic_decode(transpose(msg)))


    rm("test/$filename")
    # calculate BER using bit_stream and decoded
    ber = 1 - sum(bit_stream .== decoded)/total_bits

    CSV.write("$mod_technique.csv", DataFrame(t = pulses_per_second, BER = ber, tot = total_bits), header = false, append = true)   

    return bit_stream, decoded, ber
end

function main()

    println("----- Acoustic Communication thingy -------")

    println("Transmitting or Receiving?")

    choices = "0) Check frequency nulls\n1) Transmitting \n2) Receiving"

    println(choices); transceive = readline()


    if transceive == "0"

    achirp_fig, synched_signal_fig, transfer_function_fig = check_nulls()

    savefig(achirp_fig, "plots/Plot achirp_fig.html")
    savefig(synched_signal_fig, "plots/Plot synched_signal_fig.html")
    savefig(transfer_function_fig, "plots/Plot transfer_function_fig.html")

    else 

        # println("Choose modulation scheme"); 

        # choices = "1) FSK \n2) PSK"

        # println(choices); 

        mod_technique = "FSK" #readline()

        bits_to_send = 1
        pulses_per_second = 1

        # println("Total number of bits?"); 
        # println("Number of bits to send within a pulse?");
        total_bits = 5 #readline()

        encoding = "2" #readline()

        if encoding == "Hamming"
            row_bits = 4; col_bits = Int(total_bits/row_bits)

            bit_stream = Int.(zeros(col_bits,row_bits))
            for i in 1:total_bits
                bit_stream[i] = rand(0:1)
            end

            n = row_bits; k = 4
            encoded = zeros(Int(total_bits/row_bits), n)

            for i in 1:Int(total_bits/row_bits)
                encoded[i,1:n] = hamming_encode(bit_stream[i,1:row_bits], n, k)
            end
        elseif encoding == "2"  # "Majority Logic" 
            bit_stream = [rand(0:1) for i in 1:total_bits]
            n = 1
            encoded = majority_logic_encode(bit_stream, n); col_bits = length(encoded[:,1])

        end

        Δt = 1/fs; t = 0:Δt:(length(encoded)/bits_to_send-Δt)/pulses_per_second;

        n_carriers = 2^bits_to_send; 
        carriers = []; 
        frequencies = Array{Int}(undef,0); 
        phases = Array{Float64}(undef,0);
        amplitudes = Array{Float64}(undef,0);

        f1=2000                # initial carrier value
        spacing = 1000         # spacing between frequencies

        if mod_technique == "FSK"

            for n in 1:n_carriers
                append!(carriers, [cos.(2*π*(f1+(n-1)*spacing).*t)])
                append!(frequencies, f1+(n-1)*spacing)
            end

        elseif mod_technique == "PSK"
            
            append!(frequencies, f1)
            for n in 1:n_carriers

                ϕ = 2*π*(n-1)/n_carriers
                append!(carriers, [cos.(2*π*f1.*t .- ϕ)])
                append!(phases, ϕ)
            end

        elseif mod_technique == "ASK"
            
            append!(frequencies, f1)
            for n in 1:n_carriers

                a = (n-1)/n_carriers
                append!(carriers, [a*cos.(2*π*f1.*t)])
                append!(amplitudes, a)
            end
        
        elseif mod_technique == "QAM"

            append!(frequencies, f1)
            for n in 1:n_carriers

                a = n/n_carriers
                append!(carriers, [a*cos.(2*π*f1.*t)])
                append!(carriers, [a*sin.(2*π*f1.*t)])
                append!(amplitudes, a)
            end

        end


        if transceive == "1" # transmit
            encoded_modulated_signal_w_chirp = transmit(transpose(encoded), mod_technique, carriers, pulses_per_second, bits_to_send)

            encoded_modulated_signal_w_chirp_fig = plot_td([], encoded_modulated_signal_w_chirp, "Modulated transmitted signal with encoding")
            savefig(encoded_modulated_signal_w_chirp_fig, "plots/Plot encoded_modulated_signal_w_chirp.html")

            @show bit_stream
            @show "----------------"
            @show encoded

        elseif transceive == "2"

            recorded_signal, synched_signal, (msg, t_array, signal_array) = receive(mod_technique, spacing, frequencies, phases, pulses_per_second, col_bits, n, bits_to_send);
            # synched_signal_fig = plot_td([], synched_signal, "Received signal synchronised")
            # savefig(synched_signal_fig, "plots/Plot 1.html")

            decoded = Int.(majority_logic_decode(transpose(msg)))
            @show decoded
            
            # CSV.write("geek.csv", DataFrame(transpose([bit_stream, decoded]), :auto), header = false, append = true)   


        end

    println("------- End -------")

end

end

