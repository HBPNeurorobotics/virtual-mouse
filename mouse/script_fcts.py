#~ filename = "../../common/script_importpaths.py"
#~ exec(compile(open(filename).read(), filename, 'exec'))


import numpy as np

import math


# colorful text (yay!)
color_dictionary = {}
color_dictionary['white']  = '\033[97m'
color_dictionary['blue']   = '\033[94m'
color_dictionary['green']  = '\033[92m'
color_dictionary['red']    = '\033[91m'
color_dictionary['yellow'] = '\033[93m'
color_dictionary['cyan']   = '\033[96m'
color_dictionary['purple'] = '\033[95m'
def printCsX(str_, c = None):
    global color_dictionary
    if c==None:
        print(str(str_))
    else:
        print(color_dictionary[c]+str(str_)+'\033[0m')



def get_str_ids_in_list(str1, list_str):
    ids_ = []
    for i,istr in enumerate(list_str):
        if istr.find(str1)>=0:
            ids_.append(i)
    return ids_

def is_list_in_list(list1, list2):
    list_in = []
    for i1,il1 in enumerate(list1):
        for i2,il2 in enumerate(list2):
            if il1==il2:
                list_in.append(i1)
                break;
    return list_in


# progressbar display
def progressbar(progress, pbnum=40):
    global color_dictionary
    import sys
    loadbar = ""
    for i in range(pbnum):
        if i<=int(progress*float(pbnum)):
            loadbar = loadbar + "="
        else:
            loadbar = loadbar + " "
    #~ sys.stdout.write("\r"+"["+loadbar+"] "+str(int(progress*100.0))+"%")
    if progress < 0.33:
        sys.stdout.write("\r"+"["+color_dictionary['red'   ]+loadbar+'\033[0m'+"]"+" "+str(int(progress*100.0))+"%")
    elif progress < 0.66:
        sys.stdout.write("\r"+"["+color_dictionary['yellow']+loadbar+'\033[0m'+"]"+" "+str(int(progress*100.0))+"%")
    else:
        sys.stdout.write("\r"+"["+color_dictionary['green' ]+loadbar+'\033[0m'+"]"+" "+str(int(progress*100.0))+"%")
    if progress>=1.0:
        print(" "+'\033[0m'+"OK "+'\033[0m')




def bezierPowCsX(t):
	if t<0.0: return 0.0;
	if t>1.0: return 1.0;
	return 3.0*(t)*(t) - 2.0*(t)*(t)*(t);
	

def bezierPowNpCsX(t):
    ret_ = np.zeros(t.shape)
    ret_[np.nonzero(t>=1.0)] = 1.0
    in_ids = np.nonzero((t>=0.0) * (t<=1.0))
    t_in = t[in_ids]
    ret_[in_ids] = 3.0*(t_in)*(t_in) - 2.0*(t_in)*(t_in)*(t_in)
    return ret_



def interp3(x, y, z, v, xi, yi, zi, method='cubic'):
    from scipy.interpolate import interp1d
    """Interpolation on 3-D. x, y, xi, yi should be 1-D
    and z.shape == (len(x), len(y), len(z))"""
    q = (x, y, z)
    qi = (xi, yi, zi)
    for j in [2,1,0]:
        v = interp1d(q[j], v, axis=j, kind=method)(qi[j])
    return v

def linear_repetitive_bezier(x, l):
    rep = x/l
    rep -= l*float(np.floor(rep/l))
    return rep


def deg_to_rad(deg):
    return math.pi*deg/180.0

def rad_to_deg(rad):
    return rad*180.0/math.pi



def add(a,b):
    ret = []
    for i in range(len(a)):
        ret.append(a[i]+b[i])
    return ret


def sub(a,b):
    ret = []
    for i in range(len(a)):
        ret.append(a[i]-b[i])
    return ret


def div(a,b):
    if b==0.0:return a
    else:
        finalV = []
        for i in range(len(a)):
            finalV.append(a[i]/b)
        return finalV


def mult(a,b):
    ret = []
    for i in range(len(a)):
        ret.append(a[i]*b)
    return ret


def newref(vec, mat):
    return [
    vec[0]*mat[0][0] + vec[1]*mat[0][1] + vec[2]*mat[0][2],
    vec[0]*mat[1][0] + vec[1]*mat[1][1] + vec[2]*mat[1][2],
    vec[0]*mat[2][0] + vec[1]*mat[2][1] + vec[2]*mat[2][2]
    ]


def norm(vec):
    norm_ = 0.0
    for i in range(len(vec)):
        norm_ += vec[i]*vec[i]
    return math.sqrt(norm_)


def cross(vec1, vec2):
    return [
    vec1[1]*vec2[2] - vec1[2]*vec2[1],
    vec1[2]*vec2[0] - vec1[0]*vec2[2],
    vec1[0]*vec2[1] - vec1[1]*vec2[0],
    ] 


def dot(vec1, vec2):
    ret = 0.0
    for i in range(len(vec1)):
        ret += vec1[i]*vec2[i]
    return ret


def cross2D(vec1, vec2):
    return vec1[0]*vec2[1] - vec1[1]*vec2[0]


def sign_triangle(v1, v2, v3):
    return (v1[0]-v3[0])*(v2[1]-v3[1]) - (v2[0]-v3[0])*(v1[1]-v3[1])



def point_in_triangle(pt, v1, v2, v3):
    st1 = sign_triangle(pt, v1, v2)
    st2 = sign_triangle(pt, v2, v3)
    st3 = sign_triangle(pt, v3, v1)
    if   st1<=0.0 and st2<=0.0 and st3<=0.0:
        return True;
    elif st1>=0.0 and st2>=0.0 and st3>=0.0:
        return True;
    else:
        return False;



def inverse2x2(Tmat):
    det_ = Tmat[0,0]*Tmat[1,1] - Tmat[0,1]*Tmat[1,0]
    if det_==0.0:return Tmat*0.0;
    return np.array([[Tmat[1,1], -Tmat[0,1]], [-Tmat[1,0], Tmat[0,0]]])/det_
    

def matrix_vec_mult_2x2(mat, vec):
    return np.array([ mat[0,0]*vec[0]+mat[0,1]*vec[1] ,  mat[1,0]*vec[0]+mat[1,1]*vec[1] ])


def simple_gauss3D(vec1, vec2, sigma):
    dist = norm(sub(vec1, vec2))
    return math.exp(-dist*dist/sigma)
    


def sigmoid_CsX(x, Tk, dec):
    return 1.0/(1.0+np.exp(-(x-dec)/Tk))



