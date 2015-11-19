import numpy as np
from scipy.interpolate import interp1d

C = np.zeros((6,6))
E = 69000.0
nu = 0.3
diag = 1 - nu
off = nu
second = 1-2*nu
C[:3,:3] = off
C[0,0] = diag
C[1,1] = diag
C[2,2] = diag
C[3,3] = second
C[4,4] = second
C[5,5] = second
C *= E/((1+nu)*(1-2*nu))

stress_map = [0,50,
72.54489667,
93.22058869,
113.7482966,
134.1905299,
154.5692997,
174.8999821,
195.1930175,
215.4531558,
235.6813211,
255.8795509,
276.0478163,
295.9778533,
315.0527053,
333.1793872,
339.0694504,
344.909689,
350.711548,
356.4805061,
362.217269,
367.9270898,
373.6096429,
379.2628778,
384.8872602,
390.4817816,
396.05096,
401.5984653,
407.1287588,
412.6425772,
418.1389707,
423.6160698,
429.0736518,
434.510607,
439.9277682,
445.3236341,
450.6881221,
456.0100603,
461.2546786,
466.3235345,
471.1609634,
475.7487647,
480.0367866,
483.9800262,
487.621757,
491.0040667]

pl_map = [0,0,
2.48198E-07,
3.6038E-06,
8.74649E-06,
1.4766E-05,
2.13337E-05,
2.82263E-05,
3.52911E-05,
4.24586E-05,
4.97159E-05,
5.70343E-05,
6.44141E-05,
7.48846E-05,
9.76938E-05,
0.000134379,
0.000148847,
0.000164039,
0.000179782,
0.00019599,
0.000212652,
0.000229689,
0.000247104,
0.000264928,
0.000283156,
0.000301804,
0.000320807,
0.000340107,
0.000359636,
0.000379383,
0.000399361,
0.000419602,
0.000440109,
0.000460901,
0.000481967,
0.000503329,
0.000525135,
0.00054754,
0.000571051,
0.000597102,
0.000626521,
0.000659623,
0.000697191,
0.000739914,
0.000787158,
0.000838311]
int_pl = interp1d(stress_map, pl_map, kind='cubic')

def vToM(voigt):
	temp = np.zeros((3,3))
	temp[0,0] = voigt[0]
	temp[1,1] = voigt[1]
	temp[2,2] = voigt[2]
	temp[0,1] = voigt[5]
	temp[1,0] = voigt[5]
	temp[0,2] = voigt[4]
	temp[2,0] = voigt[4]
	temp[1,2] = voigt[3]
	temp[2,1] = voigt[3]
	return temp

def predPl(S):
	dev_S = S
	h = np.sum(dev_S[:3])/3.0
	dev_S[:3] -= h
	temp = vToM(dev_S)
	#print(temp)
	eigs, dirs = np.linalg.eig(temp)
	#print(eigs)
	max_S = np.max(eigs)
	pl = int_pl(max_S)
	return pl*dev_S/max_S
	
def predElFIP(E_tot):
	""" MUST BE VOIGT ORDER 11 22 33 23 13 12 """
	temp = np.reshape(E_tot, (6,1))
	S = np.dot(C, temp)
	pl = predPl(S)
	#print(S)
	#print(pl)
	pl = vToM(pl)
	S = vToM(S)
	pl_p, pl_dirs = np.linalg.eig(pl)
	pl_sort = np.argsort(pl_p)
	pl_shear_dir = (pl_dirs[:,pl_sort[0]]+pl_dirs[:,pl_sort[2]])/2
	pl_shear_dir = pl_shear_dir/np.linalg.norm(pl_shear_dir)
	stress = np.dot(S,pl_shear_dir)
	stress = np.dot(stress,pl_shear_dir)
	#print(stress)
	shear = pl_p[pl_sort[2]]-pl_p[pl_sort[0]]
	#print(pl_p)
	#print(pl_sort)
	#print(shear)
	return shear*(1+0.5*stress/517.0)

def predFIPs(ms, E_tot):
    fips = np.zeros(ms.shape)
    for i in range(ms.shape[0]):
        for j in range(ms.shape[1]):
            for k in range(ms.shape[2]):
                if(ms[i,j,k] != 1):
                    fips[i,j,k] = predElFIP(E_tot[:,i,j,k])
        print("Finished row %d" % i)
    return fips
    
if __name__ == "__main__":
	E_tot = np.asarray([.001,-.0003,-.0003,0,0,0])
	f = predElFIP(E_tot)
	print(f)