import cmath
import numpy as np

# Function for finding an value that is the nearest power of 2
def nearest_pow_2(x):

    pow_2_val = 1

    while pow_2_val < x:
        pow_2_val*=2
    
    return pow_2_val

# Defining a function for evaluating FFT 
def fft(polynomial_coefficients):

    number_of_coefficients = len(polynomial_coefficients)
    
    # Checking for the base condition
    if number_of_coefficients == 1:
        return polynomial_coefficients
    
    # Handling the case where the number of coefficients is not a power of 2
    pow_2_val = nearest_pow_2(number_of_coefficients)
    if(pow_2_val>number_of_coefficients):
        for i in range(pow_2_val-number_of_coefficients):
            polynomial_coefficients.append(0)
        
        number_of_coefficients = pow_2_val

    # Splitting the polynomials
    pe = []
    po = []
    for i in range(0,len(polynomial_coefficients),2):
        pe.append(polynomial_coefficients[i])
        po.append(polynomial_coefficients[i+1])

    # Evaluating the fundamental frequency
    omega = cmath.exp(2 * cmath.pi * 1j / number_of_coefficients)

    # Recursively calling the fft function
    ye = fft(pe)
    yo = fft(po)

    # Evaluting the polynomial on all the roots of the unity
    y = [0] * number_of_coefficients 
    for j in range(number_of_coefficients//2):
        y[j] = ye[j] + (omega**j)*yo[j]
        y[j + number_of_coefficients//2] = ye[j] - (omega**j)*yo[j]
    
    return y

# Defining a function for evaluating IFFT 
def ifft(frequency_coefficients):

    number_of_coefficients = len(frequency_coefficients)
    
    # Checking for the base condition
    if number_of_coefficients == 1:
        return frequency_coefficients
    
    # Handling the case where the number of coefficients is not a power of 2
    pow_2_val = nearest_pow_2(number_of_coefficients)
    if(pow_2_val>number_of_coefficients):
        for i in range(pow_2_val-number_of_coefficients):
            frequency_coefficients.append(0)
        
        number_of_coefficients = pow_2_val

    # Splitting the polynomials
    pe = []
    po = []
    for i in range(0,len(frequency_coefficients),2):
        pe.append(frequency_coefficients[i])
        po.append(frequency_coefficients[i+1])

    # Evaluating the fundamental frequency
    omega = cmath.exp( - 2 * cmath.pi * 1j / number_of_coefficients)

    # Recursively calling the fft function
    ye = fft(pe)
    yo = fft(po)

    # Evaluting the polynomial on all the roots of the unity
    y = [0] * number_of_coefficients 
    for j in range(number_of_coefficients//2):
        y[j] = (ye[j] + (omega**j)*yo[j])/number_of_coefficients
        y[j + number_of_coefficients//2] = (ye[j] - (omega**j)*yo[j])/number_of_coefficients
    
    return y


# Function for checking the correctness of our implementation of the FFT function
def check_correctness(number_of_coefficients):

    # Generate random coefficients
    # coefficients2 = np.random.rand(number_of_coefficients)
    number_of_coefficients = 5
    coefficients = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
    
    # Compute FFT using your implementation
    fft_result = fft(list(coefficients))
    ifft_result  = ifft(fft_result)

    for each in fft_result:
        print(each)
    print("\n")

    for each in ifft_result:
        print(each)


    # # Compute FFT using NumPy
    # numpy_fft_result = list(np.sort(np.fft.fft(coefficients, n=nearest_pow_2(number_of_coefficients))))

    # # Sorting both the results
    # fft_result.sort()
    # numpy_fft_result.sort()

    # # Compare the results
    # print("Your FFT result:", fft_result)
    # print("NumPy FFT result:", numpy_fft_result)
    # print("Results match:", np.allclose(fft_result, numpy_fft_result, rtol=1e-05))

check_correctness(9)

