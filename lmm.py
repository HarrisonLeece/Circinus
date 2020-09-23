'''
This program is the function and mathematics behind the LMM powering
the six degree of freedom numerical integration

Method: Adams-Bashforth s = 4(max)
'''

#These functions for testing
#https://en.wikipedia.org/wiki/Linear_multistep_method#Families_of_multistep_methods
def em(ic, derivative, h):
    return ic[0]+h*(derivative[0])
def ab2(ic,derivative, h):
    stage1 = em(ic,derivative,h)
    return stage1 + h*(1.5*derivative[1]-.5*derivative[0])
def ab3(ic,derivative,h):
    stage2 = ab2(ic,derivative,h)
    return stage2 + h*(23/12*derivative[2] - 16/12*derivative[1] + 5/12*derivative[0])
def ab4(ic,derivative,h):
    stage3 = ab3(ic,derivative,h)
    return stage3 + h*(55/24*derivative[3]-59/24*derivative[2]+37/24*derivative[1]-9/24*derivative[0])

def adams_bashford(ic,derivative, h, limiter=4):
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
        print('fxn: adams_bashford failed to satisfy any case')
        print('available derivatives: {} limiter: {}').format(available_derivatives,limiter)
        raise

    return ab_result




if __name__ == '__main__':
    h = .001

    for i in range(20):
        adams_bashford()
