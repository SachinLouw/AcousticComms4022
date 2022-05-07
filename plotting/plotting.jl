using Plots

default(size=(800,300)); #Plot canvas size
default(label=""); # Turn off legends
default(ticks=:native);  #Ticks on x-axis are labelled nicely when zooming in.

function plot_td(x_data, y_data, plot_title="", plot_xlabel="", plot_ylabel="", colour = "green", fs = 48000)
    
    Δt = 1/fs

    if !(length(x_data)>0)
        x_data = (0:1:(length(y_data)-1)) * Δt
    end
    fig = plot(x_data,y_data, color = colour)
    xlabel!(plot_xlabel); ylabel!(plot_ylabel); title!(plot_title)
    return fig

end
    
function plot_fd(x_data, y_data, plot_title="", plot_xlabel="", plot_ylabel="", ft = false, colour = "green", fs = 48000)
    
    Δt = 1/fs

    x_data = (0:1:(length(y_data)-1)) * Δt
    N = length(x_data);
    Δf = 1/(N*Δt) 
    f = (0:N/2-1)*Δf;

    if ft == false # check if fft applied

        V = fft(y_data)
        fig = plot(f, abs.(V[1:Int(floor(length(V)/2))]), color = colour)

    else

        fig = plot(f, y_data, color = colour)

    end

    xlabel!(plot_xlabel); ylabel!(plot_ylabel); title!(plot_title)
    
    return fig

end