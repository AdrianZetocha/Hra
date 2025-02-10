import tkinter as tk
from PIL import Image, ImageTk
from math import sin, cos, tan, atan, asin, radians, degrees, pi
from time import time
import json

class GameField:
    timestep = 40   # 40 ms
    field_width = 1500
    field_height = 1000
    tk = tk
    root = tk.Tk()
    canvas = tk.Canvas(width = field_width, height = field_height, bg = f'#000000')
    canvas.pack()
    root.title('Gravity')
    BlackHoleImage = tk.PhotoImage(file = 'BlackHoleResized.png').subsample(4,4)

    def __init__(self):
        # obrazok je 1400 x 700
        self.canvas.create_image(self.field_width/2, 
                self.field_height/2, image = self.BlackHoleImage)
        #r = 40     # event horizon radius
        #xc, yc = field_width/2, field_height/2      # centre coordinates
        # canvas.create_oval(xc-r,yc-r,xc+r,yc+r, outline='white')


class MovingParts(GameField):
    def __init__(self, pos_x, pos_y, angle):
        self.position = [pos_x, pos_y]
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.angle = angle         # uhol je v radianoch!


class Bullet(MovingParts):
    bullet_width = 8
    bullet_length = 40
    bullet_instances = []        # zoznam instancii nabojov
    bullets_on_canvas = []
    radius = (bullet_length**2 + bullet_width**2)**(1/2) / 2
    omega = tan(bullet_length / bullet_width)

    def __init__(self, pos_x, pos_y, angle, vel_x, vel_y, color):
        super().__init__(pos_x, pos_y, angle)
        self.vel_x = vel_x
        self.vel_y = vel_y
        self.velocity = (self.vel_x ** 2 + self.vel_y**2) ** (1/2)
        self.color = color
        self.time_of_launch = time()

        bullet = self.canvas.create_polygon(self.get_vertices(), 
                                            fill = self.color)
        Bullet.bullet_instances.append(self)
        Bullet.bullets_on_canvas.append(bullet)
    
    def get_vertices(self):
        x, y = self.pos_x, self.pos_y                       # suradnice stredu naboja 
        o, a, r = self.omega, self.angle, self.radius       # omega - stredovy uhol
        
        A = (round(x+cos(a+o) * r) , round(y + sin(a+o) * r))
        B = (round(x+cos(a-o) * r) , round(y + sin(a-o) * r)) 
        C = (round(x-cos(a+o) * r) , round(y - sin(a+o) * r)) 
        D = (round(x-cos(a-o) * r) , round(y - sin(a-o) * r))

        return *A, *B, *C, *D
    
    def pretinaju_sa(u1, u2):
        #print(u1, u2)
        # funkcia zistuje, ci sa dve usecky pretinaju
        A, B = u1[0], u1[1]     # vrcholy
        C, D = u2[0], u2[1]
        Ax, Bx, Cx, Dx = A[0], B[0], C[0], D[0]     # suradnice
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
        
        spolocne_x = Bullet.spolocna_cast(Ax, Bx, Cx, Dx)

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

        # print('Tu som')
        # print('poly1',poly1)
        # print(poly2)

        for i in range(len(poly1)):
            if i < len(poly1)-1:
                usecka1 = (poly1[i],poly1[i+1])
            else:   # posledny bod s prvym - usecka
                usecka1 = (poly1[-1],poly1[0])
            
            for j in range(len(poly2)):
                if j < len(poly2)-1:
                    usecka2 = (poly2[j],poly2[j+1])
                else:   # posledny bod s prvym - usecka
                    usecka2 = (poly2[-1],poly2[0])
                
                if Bullet.pretinaju_sa(usecka1, usecka2):
                    return True
        
        return False

    def uprava_suradnic(first_vertices, second_vertices):
        #print(second_vertices)
        x_1 = first_vertices[::2]           # slicing: [start:end:step] 
        y_1 = first_vertices[1::2]
        x_2 = second_vertices[::2]           
        y_2 = second_vertices[1::2]
        res = True

        # obdlzniky nemozu mat priesecnik ked:
        if (max(x_1) < min(x_2) or max(y_1) < min(y_2) or 
            min(x_1) > max(x_2) or min(y_1) > max(y_2)):
            return 0,0, False
        else:
            # zoskupime suradnice do bodov - do vhodneho tvaru.
            fv, sv = first_vertices, second_vertices
            fv = [(fv[i], fv[i+1]) for i in range(0,len(fv),2)]
            sv = [(sv[i], sv[i+1]) for i in range(0,len(sv),2)]
            return fv, sv, res
        

    def collision_checker(p):
        to_explode = []
        to_erase = []
        exit = False
        c = 0 # counter
        cx, cy = GameField.field_width/2, GameField.field_height/2
        w, h = GameField.field_width, GameField.field_height
        bhr = 50    # black hole radius
        m1,n1, m2,n2 = *p.P1_planet.position, *p.P2_planet.position
        R1, R2 = p.P1_planet.radius, p.P2_planet.radius
        winner_is_1, winner_is_2 = False, False

        # suradnice objektov
        baza1 = GameField.canvas.coords(p.P1_planet.delo.basis)
        delo1 = GameField.canvas.coords(p.P1_planet.delo.turret)
        baza2 = GameField.canvas.coords(p.P2_planet.delo.basis)
        delo2 = GameField.canvas.coords(p.P2_planet.delo.turret)

        #print(baza1, delo1, baza2, delo2)


        for b in Bullet.bullet_instances:
            vymazat, explosion = False, False
            # ak je daleko von z hracej plochy:
            if b.pos_x < -500 or b.pos_x > w + 500:
                vymazat = True
            elif b.pos_y < -500 or b.pos_y > h  + 500:
                vymazat = True
            # cierna diera:
            elif (cx-b.pos_x)**2 + (cy-b.pos_y)**2 < bhr**2:
                vymazat = True
            # mozno kolizia s delom 1: kruznica o 40 pixelov sirsi polomer ako planeta
            elif (m1-b.pos_x)**2 + (n1-b.pos_y)**2 < (R1+40)**2:

                # uprava suradnic do potrebneho tvaru:
                naboj = b.get_vertices()
                sb, sn1, res1 = Bullet.uprava_suradnic(baza1, naboj)
                sd, sn2, res2 = Bullet.uprava_suradnic(delo1, naboj)

                if res1 and Bullet.maju_prienik(sb, sn1):         
                    if time() - b.time_of_launch > 0.7:         
                        print('Hra sa skoncila! Vyhral P2.')
                        vymazat, explosion = True, True
                        exit, winner_is_2   = True, True

                elif res2 and Bullet.maju_prienik(sd, sn2):       
                    if time() - b.time_of_launch > 0.7:         
                        print('Hra sa skoncila! Vyhral P2.')
                        vymazat, explosion = True, True
                        exit, winner_is_2   = True, True


                # kolizia s planetou 1:
                elif (m1-b.pos_x)**2 + (n1-b.pos_y)**2 < R1**2:
                    vymazat, explosion = True, True

                    # ak naboj trafil do zakladne:
                    uhol = atan((n1-b.pos_y)/(b.pos_x-m1))
                    alfa = p.P1_planet.angle
                    omega = radians(p.P1_planet.omega)

                    if b.pos_y < n1 and b.pos_x < m1:
                        uhol = uhol + pi
                    elif b.pos_y > n1 and b.pos_x < m1:
                        uhol = uhol - pi
                    
                    # planeta 1 ide po smere hodinovych ruciciek - proti uhlom; alfa vzdy klesa!
                    # alfa zacina od 0; obluk siaha po +pi/2.
                    if (alfa + omega >= uhol >= alfa or
                         alfa + omega + 2*pi >= uhol >= alfa + 2*pi):
                        p.P1_planet.base_hp -= 1
                        p.dialog.update_hp(1, p.P1_planet.base_hp)
                        if p.P1_planet.base_hp == 0:
                            exit, winner_is_2 = True, True
                        if p.shrinking:
                            p.P1_planet.omega -= p.P1_planet.omega / 10
                            p.canvas.itemconfig(p.P1_planet.arc, 
                                            extent = p.P1_planet.omega)
                            # => bude sa zmensovat coraz pomalsie.
                            # ak by som dal namiesto /10 /p.P1_planet.base_hp,
                            # zmensovala by sa linearne.

                        # print('Trafil Zakladnu 1', p.P1_planet.base_hp)
                        # print(b.pos_x-m1, b.pos_y - n1)
                        # print()

            # mozna kolizia s planetou 2:
            elif (m2-b.pos_x)**2 + (n2-b.pos_y)**2 < (R2+40)**2:
                # uprava suradnic do potrebneho tvaru:
                naboj = b.get_vertices()
                sb, sn1, res1 = Bullet.uprava_suradnic(baza2, naboj)
                sd, sn2, res2 = Bullet.uprava_suradnic(delo2, naboj)

                # prienik bazy dela a naboja
                if res1 and Bullet.maju_prienik(sb, sn1):
                    # nema vybuchnut pri pomalom vystreleni...     
                    if time() - b.time_of_launch > 0.7:         
                        print('Hra sa skoncila! Vyhral P1.')                                          
                        vymazat, explosion = True, True
                        exit, winner_is_1  = True, True

                # prienik dela a naboja
                elif res2 and Bullet.maju_prienik(sd, sn2):      
                    if time() - b.time_of_launch > 0.7:
                        print('Hra sa skoncila! Vyhral P1.')
                        vymazat, explosion = True, True
                        exit, winner_is_1  = True, True

                # kolizia s planetou 2:
                elif (m2-b.pos_x)**2 + (n2-b.pos_y)**2 < R2**2:
                    vymazat, explosion = True, True

                    # ak naboj trafil do zakladne:
                    uhol = asin((n2-b.pos_y)/R2)
                    alfa = p.P2_planet.angle 
                    omega = radians(p.P2_planet.omega)
                    

                    if b.pos_y < n2 and b.pos_x > m2:       # 1. kvadrant
                        uhol = uhol
                    elif b.pos_y < n2 and b.pos_x < m2:     # 2. kvadrant
                        uhol = pi - uhol
                    elif b.pos_y > n2 and b.pos_x < m2:     # 3. kvadrant
                        uhol = pi - uhol
                    else:
                        uhol = 2*pi + uhol                  # 4. kvadrant
                
                    # planeta 2 ide proti smeru hodinovych ruciciek - alfa rastie.
                    # alfa zacina od 0; obluk siaha po +pi/2.
                    if (alfa <= uhol <= alfa + omega or 
                        alfa - 2*pi <= uhol <= alfa + omega - 2*pi):
                        p.P2_planet.base_hp -= 1
                        p.dialog.update_hp(2, p.P2_planet.base_hp)
                        if p.P1_planet.base_hp == 0:
                            exit, winner_is_2 = True, True
                        if p.shrinking:
                            p.P2_planet.omega -= p.P2_planet.omega / 10
                            p.canvas.itemconfig(p.P2_planet.arc,
                                             extent = p.P2_planet.omega)
                        # print('Trafil Zakladnu 2', p.P2_planet.base_hp)
                        # print('uhol',uhol,'alfa',alfa,'omega',omega)
                        # print(b.pos_x-m2, b.pos_y - n2)

            else:
                # kolizia s inym nabojom
                for i in range(c+1, len(Bullet.bullet_instances), 1):
                    b2 = Bullet.bullet_instances[i]
                    big_radius = Bullet.radius * 2
                    small_radius = Bullet.bullet_width * 2

                    x_squared = (b.pos_x - b2.pos_x)**2
                    y_squared = (b.pos_y - b2.pos_y)**2

                    if x_squared + y_squared <= big_radius**2:
                        if x_squared + y_squared <= small_radius**2:
                            to_explode.extend([c,i])
                            to_erase.extend([c,i])
                        else:
                            # treba zistit, ci maju prienik
                            first_vertices = b.get_vertices()
                            second_vertices = b2.get_vertices()
                            
                            x_1 = first_vertices[::2]           # slicing: [start:end:step] 
                            y_1 = first_vertices[1::2]
                            x_2 = second_vertices[::2]           
                            y_2 = second_vertices[1::2]

                            # obdlzniky nemozu mat priesecnik ked:
                            if (max(x_1) < min(x_2) or max(y_1) < min(y_2) or 
                                min(x_1) > max(x_2) or min(y_1) > max(y_2)):
                                pass    # (nic sa nedeje)
                            else:
                                # zoskupime suradnice do bodov - do vhodneho tvaru.
                                fv, sv, res = Bullet.uprava_suradnic(
                                        first_vertices, second_vertices)
                                if res and Bullet.maju_prienik(fv, sv):
                                    to_explode.extend([c,i])
                                    to_erase.extend([c,i])
                            
                            # Funkcie zaoberajuce sa s najdenim prieniku funguje v skratke takto:
                            # obdlzniky sa rozdelia na usecky
                            # najde sa x-ova oblast, kde su obidve usecky definovane
                            # najde sa rovnica priamky oboch usecok
                            # nakoniec sa najdu priesecniky s obdlznikovou oblastou a urobit rozsudok.  
            
            if explosion:
                if c not in to_explode:
                    to_explode.append(c)
            if vymazat:
                if c not in to_erase:
                    to_erase.append(c)
            if exit:
                p.pause(draw = False)
                p.dialog.draw_winner(winner_is_1, winner_is_2)
                break
                
            c += 1

        # Sme vonku z cyklu, takze nehrozi ze bude index out of range
        # zoradit indexy v klesajucom poradi
        to_explode = sorted(list(set(to_explode)),reverse=True)    
        to_erase = sorted(list(set(to_erase)), reverse=True)

        for i in to_explode:
            b = Bullet.bullet_instances[i]
            #print(type(b))
            b.explode()
        
        for i in to_erase:
            GameField.canvas.delete(Bullet.bullets_on_canvas[i])
            Bullet.bullet_instances.pop(i)
            Bullet.bullets_on_canvas.pop(i)
            

    
    def explode(self):
        x, y = self.pos_x, self.pos_y
        r = self.radius / 2
        kruh = GameField.canvas.create_oval(x-r,y-r,x+r,y+r,
                                    fill='white', outline='white')
        t = self.timestep 
        c = 0               # bude premennou v kvadratickej funkcii

        def rek1(t, r, c):
            if r > 0:
                GameField.canvas.coords(kruh, x-r,y-r,x+r,y+r)
                GameField.canvas.after(t, rek1, t, 
                                        (r + 5*c - c**2), c+1)          
            else:
                GameField.canvas.delete(kruh)
                
        rek1(t, r, c)


    def fly():
        c = 0   # counter
        for b in Bullet.bullet_instances:
            b.vel_x, b.vel_y, b.velocity, b.angle = Bullet.acceleration(b.pos_x,
                                                    b.pos_y, b.vel_x, b.vel_y)
            dx = b.vel_x * Bullet.timestep / 10
            dy = b.vel_y * Bullet.timestep / 10
            b.pos_x = b.pos_x + dx
            b.pos_y = b.pos_y + dy
            #print(Bullet.bullets_on_canvas[c], dx, dy)
            GameField.canvas.coords(Bullet.bullets_on_canvas[c], 
                                     b.get_vertices())
            GameField.canvas.update()
            c += 1
    
    def acceleration(pos_x, pos_y, vel_x, vel_y, slowing = None):
        # distance from black hole:
        X, Y = GameField.field_width/2, GameField.field_height/2
        r_squared = (pos_x - X)**2 + (pos_y - Y)**2
        mass = 1000
        if slowing is not None:
            mass = mass / slowing
        F = mass / r_squared

        # v tomto smere posobi sila
        direction = abs(atan((pos_y - Y) / (pos_x - X)))
        # treba rozlozit na zlozky sily: pozdlz osi x aj y

        if pos_x > GameField.field_width/2:     
            # napravo od ciernej diery
            # akceleracia do zaporneho smeru (dolava)
            sign_x = -1                        
        else:
            sign_x = 1 
            
        if pos_y > GameField.field_height/2:   
            # pod ciernou dierou 
            # akceleracia do zaporneho smeru (hore)  
            sign_y = -1                          
        else:
            sign_y = 1

        Ax = abs(cos(direction)) * F
        Ay = abs(sin(direction)) * F

        dvx = Ax * GameField.timestep
        dvy = Ay * GameField.timestep
        vel_x += dvx * sign_x
        vel_y += dvy * sign_y
        velocity = (vel_x**2 + vel_y**2)*(1/2)
        # smerovy vektor:
        angle = atan(vel_y/vel_x)

        return vel_x, vel_y, velocity, angle






class Planet(MovingParts):
    def __init__(self, pos_x, pos_y, angle, radius, color, rotation_direction, player, base_hp = None, file_name = None):
        super().__init__(pos_x, pos_y, angle)
        self.smer = rotation_direction
        self.color = color
        self.radius = radius
        self.player = player    # 1 or 2
        self.base_hp = base_hp
        
        if file_name is not None:
            self.pokracuj_v_init(file_name)

    def pokracuj_v_init(self, file_name):

        self.file_name = file_name
        self.num_of_frames = 720
        self.actual_frame = 0
        
        # turret inicialisation: at the beginning in horizontal position
        x, y = self.pos_x, self.pos_y - self.radius
        # preto uhol pi/2 radianov
        a, r, c, s = pi/2, self.radius, self.color, self.smer 
        # uhol sa prekonvertuje na radiany    
        self.delo = Turret(x, y, a, r, c, s, self.player)  

        # following central angle determines base width:
        self.omega = 90     # degrees
        self.base_hp = 10   # also concerning the base 

        # zakladna - obluk
        # (white) square in which the arc will turn:
        self.top_right_corner = [self.pos_x + self.radius - 9,
                                  self.pos_y - self.radius + 9]
        self.bottom_left_corner = [self.pos_x - self.radius + 9,
                                    self.pos_y + self.radius - 9]
        #self.canvas.create_rectangle(top_right_corner,bottom_left_corner,outline='white')


    def initialise_graphic_objects(self):
        
        self.im = Image.open(self.file_name)
        k = 2* self.radius / self.im.width
        self.moon_image = self.im.resize((round(self.im.width * k),
                                           round(self.im.height * k)))
        
        self.frames = []
        for i in range(self.num_of_frames):
            self.moon = self.moon_image.rotate(i*360/self.num_of_frames)
            self.frames.append(ImageTk.PhotoImage(self.moon))
        
        self.planet_image = self.canvas.create_image(self.pos_x, 
            self.pos_y, anchor = 'c', tags = f'planeta{self.player}',
            image = self.frames[self.actual_frame % self.num_of_frames])

        self.arc = self.canvas.create_arc(*self.top_right_corner,
            *self.bottom_left_corner, start = self.angle, 
            extent = self.omega, outline = self.color, style = tk.ARC,
            width = 15, tags = f'zakladna{self.player}')
    

    def move_planet(self):      # and base as well
        # delta alpha - angle change
        da = self.smer * 2*pi / self.num_of_frames      
        self.angle += da
        self.actual_frame += self.smer * 1

        self.canvas.itemconfig(self.arc, 
            start = round(degrees(self.angle)), extent = self.omega)
        self.canvas.itemconfig(self.planet_image, 
            image = self.frames[self.actual_frame % self.num_of_frames])
        self.delo.planet_drag(da)       # uhol sa meni zarovno
        

        if abs(self.angle) > 2*pi:
            self.angle -= 2*pi * self.smer

        # self.canvas.after(self.timestep, self.move_planet)
    
    def circular_orbit(self):
        sx, sy = GameField.field_width/2, GameField.field_height/2
        x, y = self.pos_x, self.pos_y          # 40ms stare suradnice
        self.gama += self.beta
        self.pos_x = sx + self.orbit_radius * cos(self.gama)
        self.pos_y = sy - self.orbit_radius * sin(self.gama)
        dx, dy = self.pos_x - x, self.pos_y - y

        self.many_updates()
        self.vel_x = dx/GameField.timestep
        self.vel_y = dy/GameField.timestep
        # translacna rychlost dosahuje max 0.5


    
    def planet_orbit(self):
        # toto sa robi rovnako ako s nabojmi (/10) takze rychlosti su realne
        self.vel_x, self.vel_y, v,a = Bullet.acceleration(self.pos_x, 
                self.pos_y, self.vel_x, self.vel_y, slowing = self.k)
        self.pos_x += self.vel_x * GameField.timestep / 10
        self.pos_y += self.vel_y * GameField.timestep / 10
        self.many_updates()
        # if self.player == 1:
        #     print('Planet orbit report: vel_x:', self.vel_x, 'vel_y:', self.vel_y, 'player:', self.player)

    
    def many_updates(self):
        # this method updates the position of planets, turrets and bases, if the planets are moving.
        # update planet
        self.position = [self.pos_x, self.pos_y]
        self.canvas.coords(self.planet_image, self.pos_x, self.pos_y)
        # update turret
        self.delo.radius_centre = [self.pos_x, self.pos_y]
        # update base:
        self.top_right_corner = [self.pos_x + self.radius - 7, 
                                 self.pos_y - self.radius + 7]
        self.bottom_left_corner = [self.pos_x - self.radius + 7, 
                                    self.pos_y + self.radius - 7]
        self.canvas.coords(self.arc, *self.top_right_corner, 
                                        *self.bottom_left_corner)

        

class Turret(Planet):
    def __init__(self, pos_x, pos_y, angle, radius, color, rotation_direction, player):
        super().__init__(pos_x, pos_y, angle, radius, color, rotation_direction, player)
        # pos_x, pos_y  - coordinates of the centre of the turret
        # shape of the basis:
        self.h = 10 //2          # halves for easier use
        self.bottom = 100 //2
        self.top = 80 //2
        x, y = self.pos_x, self.pos_y
        self.radius_centre = [self.pos_x, self.pos_y + self.radius]
        self.coordinates = [x+self.bottom, y + self.h, x-self.bottom, 
            y + self.h, x-self.top, y-self.h, x+self.top, y - self.h,]

        # shape of the turret:  just a rectangle
        # turret width
        self.t_w = round(Bullet.bullet_width * 1.2 // 2)        
        # turret height
        self.t_h = round(Bullet.bullet_length * 1.2 // 2)       
        self.t_coords = [x + self.t_w, y, x - self.t_w, y,
            x - self.t_w, y - self.t_h, x + self.t_w, y - self.t_h]

        self.speed = 2*pi/100
        self.tilt = 0

        # Ammo charging:
        self.krok = 1000       # ms
        self.off = False
        self.stop_token = None
        self.obdlzniky = []

    def initialise_graphic_objects(self):
        self.basis = self.canvas.create_polygon(*self.coordinates,
            fill = 'gray', outline='white',tags = f'baza{self.player}')
        self.turret = self.canvas.create_polygon(*self.t_coords, 
            fill = 'gray', outline ='white',tags = f'delo{self.player}')

    def set_easy_mode(self):
        self.bv = -5
        self.v_step = 5
        self.max_bv = 5
        self.min_bv = 5
        self.n = 2
    
    def set_harder_modes(self):
        self.min_bv = 2       # minimal bullet velocity
        self.max_bv = 8       # maximal bullet velocity
        self.n = 4
        self.v_step = (self.max_bv - self.min_bv) / self.n
        # vel. changes over time:
        self.bv = self.min_bv - self.v_step*2         
       
    def get_centre(self):
        m, n = self.radius_centre
        x = m + cos(self.angle) * self.radius
        y = n - sin(self.angle) * self.radius
        return (x, y)
    
    def get_basis_vertices(self):
        x, y = self.get_centre()
        cosine, sine = cos(self.angle), sin(self.angle)

        # bottom points;  at 90 degrees, y1 = y2 ; x1 - bottom half = x2 + bottom half
        ax = x - cosine * self.h - self.bottom * sine
        ay = y - self.bottom * cosine + self.h * sine
        bx = x - cosine * self.h + self.bottom * sine
        by = y + self.bottom * cosine + self.h * sine

        # top points;  at 0 degrees, x1 = x2 ; y1 - bottom half = y2 + bottom half
        cx = x + cosine * self.h + self.top * sine
        cy = y + self.top * cosine - self.h * sine
        dx = x + cosine * self.h - self.top * sine
        dy = y - self.top * cosine - self.h * sine

        return ax, ay, bx, by, cx, cy, dx, dy
    
    def get_turret_vertices(self):
        x, y = self.get_centre()
        cosine = cos(self.angle+self.tilt)
        sine = sin(self.angle+self.tilt)

        ax, ay = x - self.t_w * sine, y - self.t_w * cosine
        bx, by = x + self.t_w * sine, y + self.t_w * cosine
        cx = x + self.t_h * cosine + self.t_w * sine
        cy = y + self.t_w * cosine - self.t_h * sine
        dx = x + self.t_h * cosine - self.t_w * sine
        dy = y - self.t_w * cosine - self.t_h * sine
        return ax, ay, bx, by, cx, cy, dx, dy

    def move_turret(self, smer):
        # this method moves the basis (and the turret) in response to keyboard inputs
        self.angle += smer * self.speed

    def planet_drag(self, da):
        self.angle += da
        # this method moves the basis because the planet is moving
        #print('Am dragged')
        self.canvas.coords(self.basis, self.get_basis_vertices())
        self.canvas.coords(self.turret, self.get_turret_vertices())
    
    def tilt_turret(self, smer):
        # this method tilts the turret in order to aim
        if pi/2 >= self.tilt + pi / 40 * smer >= -pi/2:
            self.tilt += pi / 40 * smer

    def charge(self):
        if self.stop_token is None and self.bv < self.max_bv:
            self.bv += self.v_step
            #print(self.bv)
            if self.bv >= self.min_bv - 0.01:
                self.draw_ammo()
            
            self.off = False
            GameField.canvas.after(self.krok, self.charge)
        else:
            self.off = True
            self.draw_ammo()
    
    def shoot(self, t_vel_x, t_vel_y):
        
        if self.bv >= self.min_bv:

            # this method creates a bullet
            # Bullet init:
            # __init__(self, pos_x, pos_y, angle, vel, color, player):
            ax, ay, bx, by, cx, cy, dx, dy = self.get_turret_vertices()
            x, y = (cx+dx)/2 , (cy+dy)/2
            if self.player == 1:
                color = f'#ffaaaa'
            else:
                color = f'#aaaaff'

            bv_x = cos((self.angle + self.tilt) * (-1)) * self.bv
            bv_y = sin((self.angle + self.tilt) * (-1)) * self.bv
            #print('Shoot report:',self.player, bv_x, bv_y, 'translational vel:', t_vel_x, t_vel_y)        

            Bullet(x, y, (self.angle+self.tilt)*(-1), 
                   bv_x + t_vel_x , bv_y + t_vel_y, color)


            for o in self.obdlzniky:
                GameField.canvas.delete(o)
            self.obdlzniky.clear()

            self.bv = self.min_bv - self.v_step

            if self.off:
                GameField.canvas.after(self.krok, self.charge)
                  
    def draw_ammo(self):
        u = 40 * 5 / (self.n + 1)
        i = (self.bv - self.min_bv) / self.v_step
        if self.player == 1 and self.min_bv <= self.bv:
            o = GameField.canvas.create_rectangle(100 + u * i,
                70, 100 + u*(i+1), 70 + 20, fill = 'orange', 
                outline='black', width = 3)
            self.obdlzniky.append(o)
        elif self.player == 2 and self.min_bv <= self.bv:
            x = GameField.field_width - Dialog.dlzka_listy - 100
            y = 50
            o = GameField.canvas.create_rectangle(x + u * i, 70, 
                x + u*(i+1), 70 + 20, fill = f'#44cccc', 
                                outline='black', width = 3)
            self.obdlzniky.append(o)

    def redraw_ammo(self):
        '''This function redraws ammo when loading a saved game.'''
        ammo = (self.bv - self.min_bv)//self.v_step
        self.bv = self.min_bv
        for i in range(int(ammo)):           # so many rectangles
            self.draw_ammo()
            self.bv += self.v_step
        self.bv -= self.v_step



            


class Dialog(GameField):
    
    dlzka_listy = 200

    def __init__(self, prog):
        # trieda Dialog predsa potrebuje pristup k metodam Program-u!
        self.prog = prog
        self.draw_mainpage()
    
    def draw_mainpage(self):
        # return from pause
        try:
            self.delete_back_to_mainpage_button()
            self.clear_pause()
        except AttributeError:
            pass
        # return from game-over
        try:
            self.clear_winner()
        except AttributeError:
            pass

        w, h = GameField.field_width//2, GameField.field_height//2
        self.panel = GameField.canvas.create_rectangle(w-200,h-150,
                    w+200,h+150, fill = f'#ffaaaa', outline = 'white')
        self.title = GameField.canvas.create_text(w,h-100, 
                            text = 'GRAVITY', font = 'Arial 25 bold')
        self.new_game_button = GameField.tk.Button(text = 'New Game',
                                        command = self.start_new_game)
        self.new_game_button.place(x = w, y = h-25, anchor = 'c')
        self.last_game_button = GameField.tk.Button(
            text = 'Resume last game', command = self.prog.load_game)
        self.last_game_button.place(x = w, y = h+25, anchor = 'c')
        self.intro_button = GameField.tk.Button(text = 'Intro', 
                                                command = self.intro)
        self.intro_button.place(x = w, y = h + 75, anchor = 'c')

    def start_new_game(self):
        w, h = GameField.field_width//2, GameField.field_height//2
        self.intro_button.destroy()
        self.new_game_button.destroy()
        self.last_game_button.destroy()
        # tieto butony sa daju vytvorit aj elegantnejsie - dynamicky
        # a butony nie je nutne zakazdym znicit a vytvarat, ale daju sa schovat a zobrazit.
        self.choose = GameField.canvas.create_text(w, h-50,
                                        text = 'Choose difficulty:')
        self.button_l1 = GameField.tk.Button(text = 'level 1 - easy',
                                        command = self.set_level_1)
        self.button_l2 = GameField.tk.Button(text = 'level 2',
                                        command = self.set_level_2)
        self.button_l3 = GameField.tk.Button(text = 'level 3',
                                        command = self.set_level_3)
        self.button_l4 = GameField.tk.Button(text = 'level 4',
                                        command = self.set_level_4)
        self.button_l5 = GameField.tk.Button(text = 'level 5',
                                        command = self.set_level_5)
        self.button_l6 = GameField.tk.Button(text = 'level 6 - hard',
                                        command = self.set_level_6)
        
        self.button_l1.place(x = w, y = h - 30, anchor='c')
        self.button_l2.place(x = w, y = h + 0,  anchor = 'c',)
        self.button_l3.place(x = w, y = h + 30, anchor = 'c',)
        self.button_l4.place(x = w, y = h + 60, anchor = 'c',)
        self.button_l5.place(x = w, y = h + 90, anchor = 'c',)
        self.button_l6.place(x = w, y = h + 120, anchor = 'c',)

    def intro(self):
        self.intro_button.destroy()
        self.set_level_1()
        self.spravy = [
            ("Press space to shoot", "Press enter to shoot"),
            ("Two projectiles annihilate when they meet."),
            ("Press 'a' and 'd' to move your turret", 
            "Use arrow keys to move your turret"),
            ("Press 't' and 'g' to aim", "Press up and down arrows to aim"),
            ("Press p to pause the game"),
            ("Defend your base, but don't let enemy projectiles strike your turret: a single hit ends the game."),
            ("Otherwise, the first player to strike the opponent's base 10 times, wins."),
            ("May the Force be on your side")]

        def display_text(pair):
            if len(pair) == 2:
                self.oznam1 = GameField.canvas.create_text(350,
                    300, text = pair[0], fill = 'white', 
                                font = 'Arial 16', anchor = 'c')
                self.oznam2 = GameField.canvas.create_text(1100,
                    300, text = pair[1], fill = 'white', 
                                font = 'Arial 16', anchor = 'c')
            else:
                x,y = GameField.field_width/2, GameField.field_height/2
                self.oznam = GameField.canvas.create_text(x, y-200,
                text = pair, fill = 'white', font = 'Arial 16', 
                                                    anchor = 'c')

        def hide_text(pair):
            if len(pair) == 2:
                GameField.canvas.delete(self.oznam1, self.oznam2)
            else:
                GameField.canvas.delete(self.oznam)

        self.wait_times = [0,5,10,15,20,25,33,40, 45]  # in seconds
        self.wait_del = self.wait_times[1:] + [50]
        c = 0
        for pair in self.spravy:
            GameField.canvas.after(self.wait_times[c]*1000,
                                        display_text, pair)
            GameField.canvas.after(self.wait_del[c]*1000,
                                        hide_text, pair)
            c += 1

    # clears all mainpage widgets, including buttons, title and panel.
    def clear_mainpage(self):
        try:
            GameField.canvas.delete('all')
        except AttributeError:
            pass
        
        try:
            for j in (self.button_l1,self.button_l2, self.button_l6,
                self.button_l3, self.button_l4, self.button_l5, 
                self.last_game_button, self.new_game_button,
                self.intro_button):
                j.destroy()
        except AttributeError:
            for j in (self.last_game_button, self.new_game_button,
                       self.intro_button):
                j.destroy()

    # setters: modes 1-6
    def set_level_1(self):
        self.prog.launch(1, True, True)
    def set_level_2(self):
        self.prog.launch(2, True, True)
    def set_level_3(self):
        self.prog.launch(3, True, True)
    def set_level_4(self):
        self.prog.launch(4, True, True)
    def set_level_5(self):
        self.prog.launch(5, True, True)
    def set_level_6(self):
        self.prog.launch(6, True, True)
    
    # drawing the ammo
    def drawing(self, x, x2):

        self.h = self.dlzka_listy/10
        self.s1 = (100,100)
        self.w = GameField.field_width
        self.s2 = (self.w - 100 - self.dlzka_listy, 100)

        # Text: hp
        self.text1 = GameField.canvas.create_text(self.s1[0]-50, 
            self.s1[1] + self.h/2, anchor='w', text='Hp:', 
            font = 'Arial 16 bold', fill = f'#ffaaaa')
        self.text2 = GameField.canvas.create_text(self.s2[0]-50,
            self.s2[1] + self.h/2, anchor='w', text='Hp:', 
            font = 'Arial 16 bold', fill = f'#aaaaff')
        
        # base hp: obdlzniky
        self.frame1 = GameField.canvas.create_rectangle(*self.s1,
            self.s1[0] + self.dlzka_listy, self.s1[1]+self.h, 
            outline= 'white', fill = 'black')
        self.frame2 = GameField.canvas.create_rectangle(*self.s2, 
            self.s2[0] + self.dlzka_listy, self.s1[1]+self.h, 
            outline= 'white', fill = 'black')
        self.p1_base_hp = GameField.canvas.create_rectangle(*self.s1,
            self.s1[0] + x * self.h, self.s1[1]+self.h, 
            fill = f'#ffaaaa', outline = 'red')
        self.p2_base_hp = GameField.canvas.create_rectangle(*self.s2,
            self.s2[0] + x2 * self.h, self.s2[1]+self.h, 
            fill = f'#aaaaff', outline = 'blue')

        # Text: Ammo
        self.text3 = GameField.canvas.create_text(self.s1[0]-85, 
            self.s1[1] - self.h, anchor='w', text='Ammo:', 
            font = 'Arial 16 bold', fill = f'#ffaaaa')

        self.text4 = GameField.canvas.create_text(self.s2[0]-85,
            self.s2[1] - self.h, anchor='w', text='Ammo:', 
            font = 'Arial 16 bold', fill = f'#aaaaff')

    # draws hp1 and calls draw winner if x==0
    def update_hp(self, player, x):
        if player == 1:
            self.canvas.coords(self.p1_base_hp, *self.s1, 
                self.s1[0]+x * self.h, self.s1[1]+self.h)
        elif player == 2:
            self.canvas.coords(self.p2_base_hp, *self.s2,
                 self.s2[0]+x*self.h, self.s2[1]+self.h)
        if x == 0:
            if player == 1:
                self.draw_winner(False, True)
            else:
                self.draw_winner(True, False)


    
    # draws the game-over panel
    def draw_winner(self, w1, w2):      # w je winner_is
        w, h = GameField.field_width/2, GameField.field_height/2
        if w1 and w2:
            bg_color, color = f'#99ff99', 'green'
            napis = 'DRAW'
        elif w2:
            bg_color, color = f'#9999ff', 'blue'
            napis = 'Player 2 won'
        elif w1:
            bg_color, color = f'#ff9999', 'red'
            napis = 'Player 1 won'

        self.rect_winner = GameField.canvas.create_rectangle(w-200,
                                h-150,w+200, h+150, fill=bg_color)
        self.text_winner = GameField.canvas.create_text(w,h, 
                text = napis, font = 'Arial 16 bold', fill=color)
        self.back_to_mainpage()

    # clears the game-over panel
    def clear_winner(self):
        GameField.canvas.delete(self.rect_winner, self.text_winner)
        self.back_button.destroy()


    def draw_pause(self):
        w, h = GameField.field_width//2, GameField.field_height//2
        self.panel = GameField.canvas.create_rectangle(w-200,h-150,
                w+200,h+150, fill = f'#ffaaaa', outline = 'white')
        self.info = GameField.canvas.create_text(w,h-100, 
                text = f"Game paused. Press p to resume.", 
                fill = 'black', font = 'Arial 16 bold', anchor ='c')
        self.save_button = GameField.tk.Button(text = 'Save current game',
                            command = self.prog.save_game)
        self.save_button.place(x = w,y = h, anchor = 'c')
        self.back_to_mainpage()
        # clears pause widgets

    def clear_pause(self):
        GameField.canvas.delete(self.panel, self.info)
        self.save_button.destroy()
        self.back_button.destroy()

    def back_to_mainpage(self):
        w, h = GameField.field_width//2, GameField.field_height//2
        self.back_button = GameField.tk.Button(text = 'Back to mainpage', 
        command = self.draw_mainpage, anchor = 'c')
        self.back_button.place(x = w, y = h +50, anchor='c')
    
    def delete_back_to_mainpage_button(self):
        self.back_button.destroy()




class Program(GameField):

    def __init__(self):
        self.game_n = 0
        self.dialog = Dialog(self)
        # ostatne funkcie sa spustaju cez Dialog!

    # setting of modes is realised by clicking the buttons.
    def launch(self, mode, cast1, cast2):
        if cast1:
            self.game_n += 1              # game number
            # clear garbage from previous game
            self.clear_old()                    
            self.dialog.clear_mainpage()

            self.mode = mode
            self.set_up_default_classes()
            self.set_up_modes()

        if cast2:
            
            self.black_hole = GameField()
            self.graphic_objects()
            self.P1_planet.delo.charge()
            self.P2_planet.delo.charge()
            self.running = True
            self.run()


    def clear_old(self):
        if self.game_n > 1:
            Bullet.bullet_instances.clear()
            Bullet.bullets_on_canvas.clear()
            self.P1_planet.frames.clear()
            self.P2_planet.frames.clear()
            


    def set_up_default_classes(self):
        # this function initialises the core of the program and features related to levels
        # Core: these will be used in all levels
        self.P1_planet = Planet(pos_x= 300, pos_y= self.field_height//2,
            angle= 0, radius= 125, color= 'red', rotation_direction= -1,
            player= 1, file_name= 'Europa-moon.jpg', base_hp = 10)   
                                 
        self.P2_planet = Planet(pos_x= self.field_width - 300, 
            pos_y= self.field_height//2, angle= 0, radius= 125, 
            color= 'blue',rotation_direction= 1, player= 2, 
            file_name='Callisto-moon.jpg', base_hp = 10)
        
        self.controls = Controls(self.P1_planet, self.P2_planet, self.mode)
        GameField.root.bind_all('<p>', self.pause)

    
    def graphic_objects(self):
        for instance in (self.P1_planet, self.P2_planet):
            instance.initialise_graphic_objects()
            instance.delo.initialise_graphic_objects()

            self.dialog.drawing(instance.base_hp, instance.base_hp)

    def set_up_modes(self):
        self.orbiting = False
        self.shrinking = False
        if self.mode == 1:                   # one bullet velocity
            self.P1_planet.delo.set_easy_mode()
            self.P2_planet.delo.set_easy_mode()

        if self.mode >= 2:                   # five bullet velocities
            self.P1_planet.delo.set_harder_modes()
            self.P2_planet.delo.set_harder_modes()
        
        if self.mode >= 3:                   # shrinking bases
            self.shrinking = True

        if self.mode >= 4:                   # smaller turret bases
            k = 1.5
            self.P1_planet.delo.bottom /= k
            self.P2_planet.delo.bottom /= k
            self.P1_planet.delo.top /= k
            self.P2_planet.delo.top /= k

        if self.mode == 5:                   # circular orbit of planets
            self.P1_planet.beta = pi/720              # orbital angular velocity
            self.P1_planet.gama = (self.P1_planet.player-1) * pi    # orbital angle
            self.P1_planet.orbit_radius = self.P1_planet.pos_x - GameField.field_width/2 
            
            self.P2_planet.beta = pi/720              # orbital angular velocity
            self.P2_planet.gama = (self.P2_planet.player-1) * pi    # orbital angle
            self.P2_planet.orbit_radius = self.P1_planet.pos_x - GameField.field_width/2

            self.orbiting = True

        if self.mode == 6:                   # randomised eliptic orbit
            from random import randint, choice
            self.P1_planet.k = 6                    # slowing constant
            self.P2_planet.k = 6
            # translational velocity in the x-direction
            # translational velocity in the y-direction
            self.P1_planet.vel_y = choice((randint(-11,-8),randint(8,11)))  / self.P2_planet.k
            self.P1_planet.vel_x = (2-abs(self.P1_planet.vel_y)) * choice((-1,1))
            self.P2_planet.vel_x = - self.P1_planet.vel_x   
            self.P2_planet.vel_y = - self.P1_planet.vel_y  

            self.orbiting = True


    def pause(self, draw = True, event=None):
        '''Metoda prerusi hlavny cyklus v metode run. 
        Po vypnuti treba run znovu nastartovat.
        V tejto metode je znovu start zahrnuty.'''

        if self.running:                           
            # switch
            self.running = False
            self.P1_planet.delo.stop_token = True
            self.P2_planet.delo.stop_token = True
            print('You paused the game.')

            if draw:
                self.dialog.draw_pause()

        else:   # pause was un-paused
            self.running = True
            self.dialog.clear_pause()
            self.P1_planet.delo.stop_token = None
            self.P2_planet.delo.stop_token = None
            self.P1_planet.delo.charge()
            self.P2_planet.delo.charge()
            self.run()

            print('You resumed the game.')
            

    def run(self):
        # beh programu, jeden krok za jeden timestep
        if self.running:
            self.P1_planet.move_planet()
            self.P2_planet.move_planet()
            self.controls.key_handler()

            if self.orbiting:
                if self.mode == 5:
                    self.P1_planet.circular_orbit()
                    self.P2_planet.circular_orbit()
                elif self.mode == 6:
                    self.P1_planet.planet_orbit()
                    self.P2_planet.planet_orbit()

            Bullet.fly()
            Bullet.collision_checker(self)

            self.canvas.after(self.timestep, self.run)
        
    def save_game(self):
        # Dynamicke zapisanie do suboru a precitanie pomocou json a setattr
        print('Game saved!')
        
        # json-readable attributes grouped in dictionaries:
        program_attr = {'mode':self.mode,'running':self.running, 
            'orbiting': self.orbiting, 'shrinking':self.shrinking}
        planet_1_attr = {}    
        planet_2_attr = {}  
        turret_1_attr = {}        
        turret_2_attr = {}
        
        # we need to filter out non-json-readable attributes.
        for a, val in self.P1_planet.__dict__.items():
            if isinstance(val, (str, int, float)):
                planet_1_attr[a] = val
            # obrazky v zozname frames neslobodno ulozit.
            elif isinstance(val, list) and len(a) > 0:   
                if isinstance(val[0], (str,int,float)):
                    planet_1_attr[a] = val

        for a, val in self.P2_planet.__dict__.items():
            if isinstance(val, (str, int, float)):
                planet_2_attr[a] = val
            elif isinstance(val, list) and len(a) > 0:
                if isinstance(val[0], (str,int,float)):
                    planet_2_attr[a] = val

        for a, val in self.P1_planet.delo.__dict__.items():
            if isinstance(val, (str, int, float, list)):
                turret_1_attr[a] = val
        for a, val in self.P2_planet.delo.__dict__.items():
            if isinstance(val, (str, int, float, list)):
                turret_2_attr[a] = val

        with open('last_game.txt', 'w') as file:
            # vraj jsonovsky subor ma obsahovat len jeden objekt!
            d = [ program_attr,  planet_1_attr, planet_2_attr, 
                turret_1_attr, turret_2_attr,
                [b.__dict__ for b in Bullet.bullet_instances]]
            json.dump(d, file, indent = 3)


    def load_game(self):
        print('Loading game')

        try:
            with open('last_game.txt', 'r') as file:
                data = json.load(file)

        except FileNotFoundError:
            print('No game has been saved yet.')

        program_attributes = data[0]
        planet_1_attributes = data[1]
        planet_2_attributes = data[2]
        turret_1_attributes = data[3]
        turret_2_attributes = data[4]
        bul_ins = data[5]
        # bullets on canvas sa vytvoria spravne pocas behu programu

        # takto vzniknu instancie vsetkych potrebnych tried
        mode = data[0]['mode']
        self.launch(mode, True, False) 
        
        slovniky = [program_attributes, planet_1_attributes, 
        planet_2_attributes, turret_1_attributes, turret_2_attributes]

        instancie = [self, self.P1_planet, self.P2_planet, 
                     self.P1_planet.delo, self.P2_planet.delo]

        for i in range(len(slovniky)):
            for key, value in slovniky[i].items():
                setattr(instancie[i], key, value)
     
        # najprv vytvorim instancie Nabojov s nulovymi hodnotami
        # potom priradim tie spravne, ulozene
        for i in range(len(bul_ins)):
        # farbu treba setnut spravne, lebo je parametrom v canvas objekte
            b = Bullet(0,0,0,0,0, bul_ins[i]['color'])      
            # pocas tohto sa inicializuje aj bullets on canvas.
            for key, value in bul_ins[i].items():
                setattr(b, key, value)


        for i in (self.P1_planet.delo, self.P2_planet.delo):
            i.stop_token = None
            i.redraw_ammo()
        
        print(self.P1_planet.base_hp, self.P2_planet.base_hp)
            
        # az teraz sa vytvoria graficke objekty, budu mat spravne ID
        self.launch(self.mode, False, True)


    


class Controls(Program):
    klavesy = ('a','d','t','g','Left','Right', 'Up', 'Down')
    # (okrem space, return a pause - tie su osobitne)

    def __init__(self, planet1, planet2, mode):
        self.p1 = planet1
        self.p2 = planet2
        self.mode = mode
        self.keyboard = {'a':False, 'd':False, 'Left':False, 
                         'Right':False}
        self.root.bind('<KeyPress>', self.key_pressed_down)
        self.root.bind('<KeyRelease>', self.key_released)
        self.root.bind('<space>', self.vystrel)
        self.root.bind('<Return>', self.vystrel)


    # these two functions keep track of which keys are pressed
    def key_pressed_down(self, event):
        key = event.keysym
        if key in Controls.klavesy:         
            self.keyboard[key] = True
    
    def key_released(self, event):
        key = event.keysym
        if key in Controls.klavesy:
            self.keyboard[key] = False

    def vystrel(self, event):
        if event.keysym == 'space':
            if self.mode < 5:
                # (0,0) - translational velocity of the planet
                self.p1.delo.shoot(0,0)     
            elif self.mode >= 5:
                 self.p1.delo.shoot(self.p1.vel_x, self.p1.vel_y)
        else:
            if self.mode < 5:
                self.p2.delo.shoot(0,0)
            elif self.mode >= 5:
                self.p2.delo.shoot(self.p2.vel_x, self.p2.vel_y)

    def key_handler(self):
        #print(self.keyboard)
        for k in self.keyboard.keys():
            if self.keyboard[k] == True:
                if k == 'a':
                    self.p1.delo.move_turret(1)
                elif k == 'd':
                    self.p1.delo.move_turret(-1)
                elif k == 'Left':
                    self.p2.delo.move_turret(1)
                elif k == 'Right':
                    self.p2.delo.move_turret(-1)
                elif k == 't':
                    self.p1.delo.tilt_turret(1)
                elif k == 'g':
                    self.p1.delo.tilt_turret(-1)
                elif k == 'Up':
                    self.p2.delo.tilt_turret(-1)
                elif k == 'Down':
                    self.p2.delo.tilt_turret(1)


if __name__ == '__main__':
    p = Program()
    GameField.root.mainloop()