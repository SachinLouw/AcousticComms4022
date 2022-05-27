using Dates;

include("AcousticComms.jl")

pulses_per_second = [1 2 5 10 20 50 75]

# for p in pulses_per_second

#     for i in 1:5
#         println(automate("PSK", 1, p, 20, "Majority", 1500, 1000, 3, string(Dates.now()) * ".wav"))
#     end

# end

for p in pulses_per_second

    for i in 1:5
        println(automate("PSK", 1, p, 50, "Majority", 1500, 1000, 3, string(Dates.now()) * ".wav"))
    end

end

    # for i in 1:5
    #     println(automate("FSK", 2, 50, 2000, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
    #     println(automate("FSK", 3, 50, 2000, "Majority", 2000, 1000, 3, string(Dates.now()) * ".wav"))
    
    # end

    for i in 1:5
        println(automate("PSK", 1, 50, 3000, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
        println(automate("PSK", 1, 100, 3000, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
    
    end
    # for i in 1:5
    #     println(automate("PSK", 1, 50, 2000, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
    #     println(automate("PSK", 1, 75, 2000, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
    # end

# end

# for p in Int.(pulses_per_second./10)



#     # for i in 1:5
#     #     println(automate("FSK", 1, p, 500, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
#     # end

#     for i in 1:5
#         println(automate("FSK", 1, p*10, 200, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
#     end

#     for i in 1:5
#         println(automate("FSK", 1, p*10, 500, "Majority", 2000, 2000, 3, string(Dates.now()) * ".wav"))
#     end
# end