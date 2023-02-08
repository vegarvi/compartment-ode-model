

outfile = open('output_geo.geo', 'w')


n = 4;
d = 0.28;
L = (n)*d;
H = (n)*d;
r = 0.015
dx = 0.04
dx_c = dx/2

outfile.write('//+\nSetFactory("OpenCASCADE");\nPoint(1) = {0, 0, 0, %f};\nPoint(2) = {%f, 0, 0, %f};\nPoint(3) = {%f, %f, 0, %f};\nPoint(4) = {0, %f, 0, %f};\n\n'%(dx,L,dx,L,H,dx,H,dx))

for x in range(1,n+1):
    for y in range(1,n+1):
        p = 9*(x*n + y) - 1
        c = 9*(x*n + y)
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p, x*d - d/2, y*d - d/2, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+1, x*d - d/2 - r,  y*d - d/2, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+2, x*d - d/2 + r,  y*d - d/2, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+3, x*d - d/2,  y*d - d/2 - r, dx_c))
        outfile.write('Point(%d) = {%f, %f, 0, %f};\n'%(p+4, x*d - d/2,  y*d - d/2 + r, dx_c))
        
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+5, p+1, p, p+4))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+6, p+4, p, p+2))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+7, p+2, p, p+3))
        outfile.write('Circle(%d) = {%d, %d, %d};\n'%(p+8, p+3, p, p+1))





N = p+9

for i in range(1,5):
    outfile.write('//+\n')
    outfile.write('Line(%d) = {%d, %d};\n'%(N+i, i, i%4 + 1))






outfile.close() 
