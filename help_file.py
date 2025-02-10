import tkinter
root = tkinter.Tk()

canvas = tkinter.Canvas(width = 600, height = 600)
canvas.pack()

def func(event):
    print (event.keysym)
root.bind("<Key>", func)




def pretinaju_sa(A, B, C, D):
    # A, B, C, D su styri vrcholy dvoch useciek
    # funkcia zistuje, ci sa pretinaju

    canvas.create_line(A,B, width=3)
    canvas.create_line(C,D, width = 3)
    canvas.update()

    Ax, Bx, Cx, Dx = A[0], B[0], C[0], D[0]
    Ay, By, Cy, Dy = A[1], B[1], C[1], D[1]

    cond1 = (max(Ax, Bx) < min(Cx, Dx) or min(Ax, Bx) > max(Cx, Dx))
    cond2 = (max(Ay, By) < min(Cy, Dy) or min(Ay, By) > max(Cy, Dy))

    if cond1 or cond2:
        # tak usecky nemozu mat prienik
        return False

    # sklon je nutne urobit v spravnej interpretacii!
    # - v zavislosti od polohy bodov sa meni znamienko.

    def znamienko_sklonu(Ax, Ay, Bx, By):
        if Bx >= Ax and By <= Ay:
            znam = 1
        elif Bx >= Ax and By >= Ay:
            znam = -1
        elif Bx <= Ax and By <= Ay:
            znam = -1
        elif Bx <= Ax and By >= Ay:
            znam = 1
        return znam

    # sklon prvej usecky:
    dx = Ax - Bx
    dy = Ay - By
    if dx == 0:
        k1 = 0
    else:
        k1 = abs(dy/dx) * znamienko_sklonu(Ax,Ay, Bx, By)

    # sklon druhej usecky:
    dx = Cx - Dx
    dy = Cy - Dy
    if dx == 0:
        k2 = 0
    else:
        k2 = abs(dy/dx) * znamienko_sklonu(Cx, Cy, Dx, Dy)
    
    spolocne_x = spolocna_cast(Ax, Bx, Cx, Dx)


    # kvoli otocenej symetrii suradnicovej sustavy ( y rastie, ako ide dolu),
    # treba prehodit znamienka: y --> -y
    # konstantny koeficient linearnej rovnice:
    # -y = ax + b => b = -y -ax
    b1 = - Ay - k1*Ax
    b2 = - Cy - k2*Cx

    # priesecniky: y-suradnica
    # l - at lower bound; g - at greater bound
    l1 = (k1 * spolocne_x[0] + b1) * (-1)
    g1 = (k1 * spolocne_x[1] + b1) * (-1)
    l2 = (k2 * spolocne_x[0] + b2) * (-1)
    g2 = (k2 * spolocne_x[1] + b2) * (-1)

    if spolocne_x[0] == spolocne_x[1]:
        # vdaka predchadzajucej minmax-ovej podmienke
        # mozeme rovno vratit:
        return True

    # print('Prave sa riesi dvojica:')
    # print(A, B, ";" , C, D)
    # print('k1 =', k1)
    # print('k2 =', k2)
    # print('b1 =', b1)
    # print('b2 =', b2)
    # print('l1 =',l1)
    # print('g1 =',g1)
    # print('l2 =',l2)
    # print('g2 =',g2)

    if l1 >= l2 and g1 >= g2:
        # => prva usecka je vzdy pod druhou useckou => ziadny prienik
        return False
    elif l1 <= l2 and g1 <= g2:
        # => prva usecka je vzdy nad druhou useckou => ziadny prienik
        return False
    else:
        return True
    

def spolocna_cast(a, b, c, d):
    # funkcia funguje iba pre body, ktore naisto maju spolocnu cast

    mensi_1 = min(a,b)
    vacsi_1 = max(a,b)
    mensi_2 = min(c,d)
    vacsi_2 = max(c,d)
    
    if mensi_1 <= mensi_2 and vacsi_2 <= vacsi_1:
        spolocna = [mensi_2,vacsi_2]
    elif mensi_1 <= mensi_2 and vacsi_2 >= vacsi_1:
        spolocna = [mensi_2, vacsi_1]
    elif mensi_1 >= mensi_2 and vacsi_1 >= vacsi_2:
        spolocna = [mensi_1, vacsi_2]
    elif mensi_1 >= mensi_2 and vacsi_1 <= vacsi_2:
        spolocna = [mensi_1, vacsi_1]

    return spolocna

def maju_prienik(poly1, poly2):
    # funkcia postupne kontroluje vsetky dvojice hran, ktore by sa mohli pretnut.
    # polygon je v tvare:
    # [(a,b),(c,d),(e,f),(g,h)...]

    for i in range(len(poly1)):
        #print(poly1[i], poly1[i+1])
        if i < len(poly1)-1:
            dvojica = (poly1[i],poly1[i+1])
        else:   # posledny bod s prvym - usecka
            dvojica = (poly1[-1],poly1[0])
        
        for j in range(len(poly2)):
            #print(poly2[j], poly2[j+1])
            if j < len(poly2)-1:
                dvoj = (poly2[j],poly2[j+1])
            else:   # posledny bod s prvym - usecka
                dvoj = (poly2[-1],poly2[0])
            
            #print('\n',dvojica, dvoj,'\n')
            
            if pretinaju_sa(*dvojica, *dvoj):
                return True
    
    return False
            
'''
# test 1. podmienky
print(spolocna_cast(*(1,5),*(2,10)))
# test 2. podmienky
print(spolocna_cast(*(10,15),*(12,20)))
# test 3. podmienky
print(spolocna_cast(*(13,20),*(10,15)))
# test 4. podmienky
print(spolocna_cast(*(15,18),*(11,20)))

print(spolocna_cast(*(100,50),*(75,75)))
'''
'''
print(pretinaju_sa((100,100),(200,200),(100,200),(200,100)))
print(pretinaju_sa((200,200),(100,100),(200,100),(100,200)))

print(pretinaju_sa((200,200),(200,0),(250,100),(100,200)))

print(pretinaju_sa((0,0),(10,10),(0,5),(10,5)))
print(pretinaju_sa((10,10),(30,60),(10,50),(50,50)))


vrcholy2 = [(100,200),(200,200)]    
vrcholy1 = [(200,100),(150,220)] 
canvas.create_line(vrcholy1, fill ='blue')
canvas.create_line(vrcholy2, fill = 'red')

print(maju_prienik(vrcholy1, vrcholy2))



stvorec = [(100,100),(100,0),(0,0),(0,100)]
trojuholnik = [(90,50),(150,0),(150,150)]
canvas.create_polygon(stvorec, fill ='blue')
canvas.create_polygon(trojuholnik, fill = 'red')

print(maju_prienik(stvorec,trojuholnik))

# polygon - vstupne body sa napajaju za sebou
stvorec = [(100,200),(100,300),(400,300),(400,200)]
obdlznik = [(200,50),(300,50),(300,400), (200,400)]
canvas.create_polygon(stvorec, fill ='blue')
canvas.create_polygon(obdlznik, fill = 'red')

print(maju_prienik(stvorec,obdlznik))
'''
from math import sin, cos, pi

# funkcia z hry:
def get_vertices(x, y, a, o, r):
    
    A = (round(x+cos(a+o) * r) , round(y + sin(a+o) * r))
    B = (round(x+cos(a-o) * r) , round(y + sin(a-o) * r)) 
    C = (round(x-cos(a+o) * r) , round(y - sin(a+o) * r)) 
    D = (round(x-cos(a-o) * r) , round(y - sin(a-o) * r))

    return A, B, C, D


prvy = get_vertices(100,100,pi/3, pi/6, 100)
druhy = get_vertices(250,100,4*pi, pi/6, 100)
canvas.create_polygon(prvy, fill ='green')
canvas.create_polygon(druhy, fill = 'yellow')
print(maju_prienik(prvy, druhy))

root.mainloop()


