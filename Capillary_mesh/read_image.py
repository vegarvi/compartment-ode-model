import cv2
import matplotlib.pyplot as plt
from numpy import array, savetxt
from IPython import embed


class TestClass():
    def __init__(self, imagename):
        self.fname = imagename
        self.img = cv2.imread(self.fname)

    def getCoord(self, plottitle = 'Mark points'):
        self.point = []
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.imshow(self.img)
        cid = fig.canvas.mpl_connect('button_press_event', self.__onclick__)
        plt.title(plottitle)
        plt.show()
        return self.point

    def __onclick__(self,click):
        self.point.append((click.xdata,click.ydata))
        return self.point



a = TestClass('image_of_vessels.png')
arteries = a.getCoord('Click all arteries, then exit')
veins = a.getCoord('Click all veins, then exit')


ax = array([[0.080297061429064343, 629.46294271765964], [837.53717980133104, 630.65420428627976], [838.1328105856411, 29.662742917459923], [0.080297061429064343, 27.280219780219909], [611.19748176351959, 2.8593576235086857], [825.02893333082045, 5.2418807607488134]])
x_bar = ax[-1,0] - ax[-2,0]      # 1 mm bar corresponds to x_bar in image units
y_bar = - x_bar
origin = ax[0]/x_bar                   # origin coordinates in mm
L = (ax[1,0] - ax[0,0])/x_bar
H = (ax[3,1] - ax[0,1])/y_bar


def rescale(points):
    points = (points-ax[0])
    points[:,0] *= x_bar**-1
    points[:,1] *= y_bar**-1
    return points


arteries = rescale(arteries)
veins = rescale(veins)

savetxt('origin_L_H.txt',[origin[0],origin[1], L, H])
#savetxt('arteries2.txt', arteries)
#savetxt('veins2.txt', veins)




print(arteries, veins)
