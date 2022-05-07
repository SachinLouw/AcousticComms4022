# https://www.geeksforgeeks.org/linear-block-code-using-matlab/




using LinearAlgebra

# Given H Matrix

H = [1 0 1 1 0 0    # parity check 
     1 1 1 0 1 0
     0 1 1 0 0 1]

k = 3;
n = 6;

# Generating G Matrix

# Taking the H Matrix Transpose
P = transpose(H);

# Making a copy of H Transpose Matrix
L = P[1:k,:]


# Making a 4 x 7 Matrix
G = [I(k) L]         # generator matrix
# G = []; append!(G, I(k)); append!(G, L)

# Generate U data vector, denoting all information sequences
no = 2 ^ k

u = zeros(no,k)

# Iterate through an Unit-Spaced Vector
for i in 1 : 2^k

    # Iterate through Vector with Specified Increment
    # or in simple words here we are decrementing 4 till we get 1
    

    for j in k : -1 : 1
        if rem(i - 1, 2 ^ (-j + k + 1)) >= 2 ^ (-j + k)
            u[i,j] = 1
        end    
    end

end


# Generate CodeWords
c = rem.(u * G, 2)
# @show c

# Find the min distance
sum_c = []

for i in 2:length(transpose(c)[1,:]) # column width

    # @show sum(transpose(c)[:,i])
    append!(sum_c, sum(transpose(c)[:,i]))

end

w_min = findmin(sum_c) # ?? sum columnwise entries and find minimum


# Given Received codeword
r = [0 1 0 1 1 1];

p = G[:, n-k+2:n];


#Find Syndrome
ht = transpose(H)

s = rem.(r * ht, 2)
@show ht
@show s

i = 1

while i < length(ht[:,1])
    # @show ht[i,1:3] .== s == true
    
    if(vec(ht[i,:])==vec(s))                # ith row, indices 1:3
        r[i] = 1-r[i];
        break;
    end
    global i += 1;
end

println("The Error is in bit:")
@show(i)

println("The Corrected Codeword is :")
@show(r)
