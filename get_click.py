import numpy as np
import matplotlib.pyplot as plt

# class TestClass():
#     def __init__(self):
#         self.fname = 'image.jpg'
#         self.img = cv2.imread(self.fname)
#         self.point = ()

def get_click_coord(array):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(array)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return(point)

def onclick(click):
    global point
    point = (click.xdata,click.ydata)
    return(point)

def main():
    image_array = np.zeros((5,5))
    coords = get_click_coord(image_array)
    print(coords)

main()
