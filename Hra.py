import tkinter as tk
from math import sin, cos, tan


class Space:        # Tu su iba triedne atributy, bez init
    timestep = 40   # 40 ms
    sirka = 1200
    vyska = 800
    root = tk.Tk()
    canvas = tk.Canvas(width = sirka, height = vyska, bg = f'#000000')
    canvas.pack()

    # obrazok je 1400 x 700
    BlackHoleImage = tk.PhotoImage(file = 'BlackHole.png').subsample(4,4)
    canvas.create_image(sirka/2, vyska/2, image = BlackHoleImage)


class Projectile(Space):
    proj_list = []

    def __init__(self, x=100, y=100, vel=50, angle=90):

        self.x_pos = x          # center of mass coords
        self.y_pos = y
        self.velocity = vel
        self.angular_orientation = angle

        self.width = 40
        self.lenght = 40

        a = (x+cos(angle) * self.lenght/2 , y + sin(angle) * self.width/2) 
        b = (x+cos(angle) * self.lenght/2 , y - sin(angle) * self.width/2) 
        c = (x-cos(angle) * self.lenght/2 , y + sin(angle) * self.width/2) 
        d = (x-cos(angle) * self.lenght/2 , y - sin(angle) * self.width/2) 

        bullet = self.canvas.create_polygon(a,b,d,c, fill = 'white')
        Projectile.proj_list.append(bullet)





    def gravity(self):
        center_x = Program.sirka /2
        center_y = Program.vyska /2
        M = 100     # mass of the Black Hole
        m = 1       # mass of the Projectile
        r = ((self.x_pos - center_x)**2 + (self.y_pos - center_y)**2)   # distance of projectile from Black Hole
        Fg = (M * m) / (r**2)
        direction_of_force = tan((self.pos_y - center_y) / (self.pos_x - center_x))

        # nasledne treba updatovat parametre projektilu - rychlost a smer 
    
    def black_hole(self):
        # ak sa projektil dostane k event horizon, zanika.
        ...

    def collision(self):
        # ak sa stretne s projektilom supera, vzajomne sa anihiluju
        ...

    def hit(self):
        # ak narazi do pevneho predmetu, zanika
        ...

    def fly(self):
        # iba ine slovo pre move
        for bullet in Projectile.proj_list:
            s = self.velocity * self.timestep / 1000
            a = self.angular_orientation
            self.canvas.move(bullet, s * cos(a), s * sin(a))
            if 0 < self.x_pos < self.sirka and 0 < self.y_pos < self.vyska:
                self.canvas.after(self.timestep, self.fly) 

   



class Controls(Space):

    def __init__(self):
        self.keyboard = {'a':False, 'd':False, 'Left':False, 'Right':False}
        self.root.bind('<KeyPress>', self.key_pressed_down)
        self.root.bind('<KeyRelease>', self.key_released)
        self.root.bind_all('<space>',self.shoot)
        Projectile()


    def key_pressed_down(self, event):
        key = event.keysym            
        self.keyboard[key] = True
    
    def key_released(self, event):
        key = event.keysym
        self.keyboard[key] = False

    #def run(self):
    #    self.root.mainloop()     # musi byt na konci


   
    # spravit pre Player1 aj Player2
    def shoot(self, event):
        Projectile()
        print('Hello, shoot')
    #Space.canvas.bind_all('<space>', hoot)

    def P1_Move_Left(event):
        print('a')
    Space.canvas.bind_all('a',P1_Move_Left)

    def P1_Move_Right(event):
        ...
    Space.canvas.bind_all('a',P1_Move_Right)


    
class Turret:
    def __init__(self, x, y, vel, angle):
        self.x_pos = x
        self.y_pos = y
        self.velocity = vel
        self.angular_orientation = angle

        self.width = 20
        self.length = 70
    
class Base:
    def __init__(self, width):
        self.arc_width = width


class Planet(Space):
    def __init__(self, r, speed, left_edge, right_edge):
        self.radius = r
        self.angular_velocity = speed
        self.canvas.create_oval()
    
    def revolve(self):
        ...


class Program:
    # inicializuj ovladanie:
    c = Controls()
    Space.root.mainloop()
    
p = Program()
print(p.c.keyboard)
