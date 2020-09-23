'''
This program is the function and mathematics behind the LMM powering
the six degree of freedom numerical integration

Method: Adams-Bashforth s = 4(max), uses between 0 and 3 previous data points
'''

#These functions do the mathy bits
#https://en.wikipedia.org/wiki/Linear_multistep_method#Families_of_multistep_methods
#ignore the commented out lines for now, they were useful to me earlier
def em(ic, derivative, h):
    return ic+h*(derivative[0])
def ab2(ic,derivative, h):
    #stage1 = em(ic,derivative,h)
    return ic + h*(1.5*derivative[1]-.5*derivative[0])
def ab3(ic,derivative,h):
    #stage2 = ab2(ic,derivative,h)
    return ic + h*(23/12*derivative[2] - 16/12*derivative[1] + 5/12*derivative[0])
def ab4(ic,derivative,h):
    #stage3 = ab3(ic,derivative,h)
    return ic + h*(55/24*derivative[3]-59/24*derivative[2]+37/24*derivative[1]-9/24*derivative[0])

'''
This function handles the switching in the following cases:
The derivative list does not have enough pre-computed derivatives to succeed
in the computation of the desired stage of the multi-step method
A limiter is supplied to improve the computation speed of the method

This function is also important for handling the derivative list, and converting
its format to a format which the above functions understand.
@arguments: ic should be the value of the initial conditions
derivative should be a list containing the previously computed rates of changes
and the 'current' rate of change for the next step
h should be the step size of the simulation
limiter is how many derivatives should be supplied to the numerical method maximum
'''
def adams_bashford(ic,derivative, h, limiter=4):
    #the derivative argument supplied to the function may be any length, but
    #the switching logic and computation will be confused unless derivative is
    #truncated to

    try:
        derivative = derivative[(-1*limiter):-1]
    except:
        #if we hit this except, there is no need to truncate the derivative list
        #before the switcher (number of available derivates in the list < limiter)
        pass
    available_derivatives = len(derivative)
    if (available_derivatives == 1 or limiter ==  1):
        ab_result = em(ic,derivative, h)
    elif (available_derivatives == 2 or limiter == 2):
        ab_result = ab2(ic,derivative,h)
    elif (available_derivatives == 3 or limiter == 3):
        ab_result = ab3(ic,derivative,h)
    elif (available_derivatives == 4 and limiter == 4):
        ab_result = ab4(ic,derivative,h)
    else:
        print('fxn: adams_bashford arguments failed to satisfy any case')
        print('available derivatives: {} limiter: {}').format(available_derivatives,limiter)
        raise

    return ab_result




if __name__ == '__main__':
    h = .5
    position_at_step = [2]
    derivative_list = [3]
    #derivative_list should be [3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    #by the end of execution
    for i in range(20):
        next_step = adams_bashford(position_at_step[-1],derivative_list, h)
        position_at_step.append(next_step)
        derivative_list.append(i)
    print('derivative_list: {}').format(derivative_list)
    print('position_at_step: {}').format(position_at_step)
