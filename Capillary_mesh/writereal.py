from pylab import *
from IPython import embed



art = loadtxt('arteries.txt')
vei = loadtxt('veins.txt')


def check_dupl(A):
    min_dist = 0.01
    n = len(A)
    for i in range(n):
        for j in range(n):
            if norm(A[i] - A[j]) < min_dist:
                if i!=j:
                    return True

    return False
assert(not check_dupl(art) and not check_dupl(vei))



r = 0.015
dx = 0.04
dx_c = dx/2
x0,y0,L,H = loadtxt('origin_L_H.txt')
print(L,H)
#origin = array(x0,y0)   # don't need origin, arteries/veins already scaled to origin

def exclude_boarder_vessels(A, A_r):
    new_vessels = A.copy()
    rerun = True
    while rerun:
        rerun = False
        for i in range(len(new_vessels)):
            p = new_vessels[i]
            if (p[0] > (L - 2*A_r)) or (p[0] < 2*A_r) or (p[1] > (H - 2*A_r)) or (p[1] < 2*A_r):
                print('found boarder-vessel', len(new_vessels),i)
                new_vessels = append(new_vessels[:i], new_vessels[i+1:]).reshape(len(new_vessels)-1,2)
                rerun = True
                break;
    return new_vessels

art = exclude_boarder_vessels(art, r)
vei = exclude_boarder_vessels(vei, r)


if __name__ == '__main__':
    outfile = open('real_output_geo_large.geo', 'w')


    outfile.write('Point(1) = {0, 0, 0, %f};\nPoint(2) = {%f, 0, 0, %f};\nPoint(3) = {%f, %f, 0, %f};\nPoint(4) = {0, %f, 0, %f};\n\n'%(dx,L,dx,L,H,dx,H,dx))

    n_a = len(art)
    n_v = len(vei)

    artery_lines = []

    for i in range(1,n_a+1):
        x,y = art[i-1,0], art[i-1,1]
        p = 9*i
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p, x, y, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+1, x-r, y, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+2, x, y+r, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+3, x+r, y, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+4, x, y-r, dx_c))

        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+5, p+1, p, p+2))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+6, p+2, p, p+3))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+7, p+3, p, p+4))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+8, p+4, p, p+1))
        artery_lines.append([p+5,p+6,p+7,p+8])

    artery_end = p+9

    vein_lines = []
    r = 0.02             # OBSOBS!!!
      
    for i in range(artery_end,n_v+artery_end):
        x,y = vei[i-artery_end,0], vei[i-artery_end,1]
        p = 9*i
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p, x, y, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+1, x-r, y, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+2, x, y+r, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+3, x+r, y, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+4, x, y-r, dx_c))

        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+5, p+1, p, p+2))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+6, p+2, p, p+3))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+7, p+3, p, p+4))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+8, p+4, p, p+1))
        vein_lines.append([p+5,p+6,p+7,p+8])





    N = p+9

    for i in range(1,5):
        outfile.write('//+\n')
        outfile.write('Line(%d) = {%d, %d};\n'%(N+i, i, i%4 + 1))

    N = N+i+1
    outfile.write('Line(%d) = {1, 2};\nLine(%d) = {2, 3};\nLine(%d) = {3, 4};\nLine(%d) = {4, 1};\n'%(N, N+1, N+2, N+3))
    outfile.write('Line Loop(%d) = {%d, %d, %d, %d};\n'%(N+4, N, N+1, N+2, N+3))
    loop_start = N+4
    N = N+5
    for i in range(N,N+len(artery_lines)):
        idx = i - N
        outfile.write('Line Loop(%d) = {%d, %d, %d, %d};\n'%(i,artery_lines[idx][0],artery_lines[idx][1],artery_lines[idx][2],artery_lines[idx][3]))

    N = N+len(artery_lines)
    for i in range(N,N+len(vein_lines)):    
        idx = i - N
        outfile.write('Line Loop(%d) = {%d, %d, %d, %d};\n'%(i,vein_lines[idx][0],vein_lines[idx][1],vein_lines[idx][2],vein_lines[idx][3]))

    loop_end = i
    outfile.write('Plane Surface(1) = {')
    for j in range(loop_start,loop_end):
        outfile.write('%d, '%j)
    outfile.write('%d};\n'%loop_end)

    outfile.close()
