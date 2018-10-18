# -*- coding: utf-8 -*-
from __future__ import division
from sympy import *
import matplotlib
import matplotlib.pyplot as plt
y, t = symbols('y t')

#globals:
c_bashforth =   [
                [1],
                [3/2, -1/2],
                [22/12, -4/3, 5/12],
                [55/24, -59/24, 37/24, -3/8],
                [1901/720, -1387/360, 109/30, -637/360, 251/720],
                [4277/1440, -2641/480, 4991/720, -3649/720, 959/480, -95/288],
                [198721/60480, -18367/2520, 235183/20160, -10754/945, 135713/20160, -5603/2520, 19087/60480],
			    [16083/4480, -1152169/120960, -242653/13440, 296053/13440, 2102243/120960, -115747/13440, 32863/13440, -5257/17280]
                ]
c_multon =  [
                [1],
                [1/2 , 1/2],
                [5/12 , 2/3 , -1/12],
                [3/8 , 19/24 , -5/24 , 1/24],
                [251/720 , 323/360 , -11/30 , 53/360 , -19/720],
                [95/288 , 1427/1440 , -133/240 , 241/720 , -173/1440 , 3/160],
                [19087/60408 , 2713/2520 , -15487/20160 , 586/945 , -6737/20160 , 263/2520 , -863/60408],
                [5257/17280 , 139849/120960 , -4511/4480 , 123133/120960 , -88547/120960 , 1537/4480 , -11351/120960 , 275/24192]
            ]
vx = []
vy = []

def euler(y0, t0, h, n, fun):
    #Output's first line
    print("Metodo de Euler\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    vy.append(yNow)
    vx.append(tNow)
    #Loop for calculation and printing
    for i in range(1 , n + 1):
        F = fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        yNow += h*F
        tNow += h

        vx.append(tNow)
        vy.append(yNow)

        print(str(i) + " " + str(yNow))

    return yNow

def euler_inverso(y0, t0, h, n, fun):
    #Output's first line
    print("Metodo de Euler Inverso\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0  
    tNow = t0
    vx.append(tNow)
    vy.append(yNow)
    #Loop for calculation and printing
    for i in range(1 , n + 1):
        #Using the result from euler to calculate the next step
        yAfter = yNow + h*fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        tAfter = tNow + h
        F = fun.subs([(y,yAfter),(t,tAfter)])
        yNow += h*F
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))

    return

def euler_aprimorado(y0, t0, h, n, fun):
    #Output's first line
    print("Metodo de Euler Aprimorado\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    vx.append(tNow)
    vy.append(yNow)
    #Loop for calculation and printing
    for i in range(1 , n + 1):
        yAfter = vy[i-1] + h*fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        tAfter = vx[i-1] + h
        fAfter = fun.subs([(y,yAfter),(t,tAfter)])
        fNow = fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        yNow += h/2*(fAfter + fNow)
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))
    
    return

def runge_kutta(y0, t0, h, n, fun):
    #Output's first line
    print("Metodo de Runge Kutta\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    #Loop for calculation and printing
    for i in range(1 , n + 1):
        k1 = fun.subs([(y,yNow),(t,tNow)])
        k2 = fun.subs([(y,yNow + 0.5*h*k1),(t,tNow + 0.5*h)])
        k3 = fun.subs([(y,yNow + 0.5*h*k2),(t,tNow + 0.5*h)])
        k4 = fun.subs([(y,yNow + h*k3),(t,tNow + h)])

        yNow += h/6*(k1 + 2*k2 + 2*k3 + k4)
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
    
        print(str(i) + " " + str(yNow))
 
    return
 
def adam_bashforth(y0, t0, h, n, fun, grau):
    #Output's first line
    print("Metodo de Adams Bashforth\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    for i in range(1, grau):
        aux = t0 + h*(i-1)
        vx.append(aux)
    #Loop to show preseted values of y, excluding y0
    for i in range(1, grau):
        print(str(i) + " " + str(vy[i-1]))
    #Loop for calculation and printing
    for i in range(grau , n + 1):
        fAux = 0
        for j in range(0, grau):
            fAux += c_bashforth[grau-1][j]*fun.subs([(y,vy[i-j-2]),(t,vx[i-j-2])])
        yNow = vy[i-2] + h*fAux
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))

    return yNow

def adam_bashforth_by_euler(y0, t0, h, n, fun, grau):
    #Output's first line
    print("Metodo de Adams Bashforth por Euler\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    vy.append(yNow)
    vx.append(tNow)
    #Calculate the past values of y, using the euler method
    for i in range(1, grau):
        F = fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        yNow += h*F
        tNow += h

        vx.append(tNow)
        vy.append(yNow)

        print(str(i) + " " + str(yNow))
    #Loop for calculation and printing
    for i in range(grau , n + 1):
        fAux = 0
        for j in range(0, grau):
            fAux += (c_bashforth[grau - 1][j])*(fun.subs([(y,vy[i-j-1]),(t,vx[i-j-1])]))
        yNow = vy[i-1] + h*fAux
        tNow = vx[i-1] + h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))

    return

def adam_bashforth_by_euler_inverso(y0, t0, h, n, fun, grau):
    #Output's first line
    print("Metodo de Adams Bashforth por Euler Inverso\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    vy.append(yNow)
    vx.append(tNow)
    for i in range(1 , grau):
        #Using the result from euler to calculate the next step
        yAfter = yNow + h*fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        tAfter = tNow + h
        F = fun.subs([(y,yAfter),(t,tAfter)])
        yNow += h*F
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))
    #Loop for calculation and printing
    for i in range(grau , n + 1):
        fAux = 0
        for j in range(0, grau):
            fAux += c_bashforth[grau-1][j]*fun.subs([(y,vy[i-j-1]),(t,vx[i-j-1])])
        yNow = vy[i-1] + h*fAux
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))

    return

def adam_bashforth_by_euler_aprimorado(y0, t0, h, n, fun, grau):
    #Output's first line
    print("Metodo de Adams Bashforth por Euler Aprimorado\ny({0}) = {1}\nh = {2}\n0 {1}".format(t0,y0,h))
    yNow = y0
    tNow = t0
    vy.append(yNow)
    vx.append(tNow)
    #Calculate the past values of y, using the euler method
    for i in range(1 , grau):
        yAfter = vy[i-1] + h*fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        tAfter = vx[i-1] + h
        fAfter = fun.subs([(y,yAfter),(t,tAfter)])
        fNow = fun.subs([(y,vy[i-1]),(t,vx[i-1])])
        yNow += h/2*(fAfter + fNow)
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))
    #Loop for calculation and printing
    for i in range(grau , n + 1):
        fAux = 0
        for j in range(0, grau):
            fAux += c_bashforth[grau-1][j]*fun.subs([(y,vy[i-j-1]),(t,vx[i-j-1])])
        yNow = vy[i-1] + h*fAux
        tNow += h

        vx.append(tNow)
        vy.append(yNow)
        
        print(str(i) + " " + str(yNow))

    return

def main():
    with open('entrada.txt') as arq:
        entrada = arq.read().splitlines()
    for lines in entrada:
        argumentos = lines.split()
        #methods - - - - - - - - - - - - - -
        if(argumentos[0] == "euler"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            euler(y0, t0, h, n, fun)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()
            

        if(argumentos[0] == "euler_inverso"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            euler_inverso(y0, t0, h, n, fun)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()

        if(argumentos[0] == "euler_aprimorado"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            euler_aprimorado(y0, t0, h, n, fun)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()

        if(argumentos[0] == "runge_kutta"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            runge_kutta(y0, t0, h, n, fun)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()

        if(argumentos[0] == "adam_bashforth"):
            grau = int(argumentos[-1])
            for i in range(1,grau):
                aux = float(argumentos[i])
                vy.append(aux)
            t0 = float(argumentos[-5])
            h =  float(argumentos[-4])
            n = int(argumentos[-3])
            fun = sympify(argumentos[-2])
            adam_bashforth(vy[0], t0, h, n, fun, grau)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()

        if(argumentos[0] == "adam_bashforth_by_euler"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            grau = int(argumentos[6])
            adam_bashforth_by_euler(y0, t0, h, n, fun, grau)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()

        if(argumentos[0] == "adam_bashforth_by_euler_inverso"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            grau = int(argumentos[6])
            adam_bashforth_by_euler_inverso(y0, t0, h, n, fun, grau)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()

        if(argumentos[0] == "adam_bashforth_by_euler_aprimorado"):
            y0 = float(argumentos[1])
            t0 = float(argumentos[2])
            h =  float(argumentos[3])
            n = int(argumentos[4])
            fun = sympify(argumentos[5])
            grau = int(argumentos[6])
            adam_bashforth_by_euler_aprimorado(y0, t0, h, n, fun, grau)
            print("\n")
            #Ploting graphic
            plt.xlabel("t")
            plt.ylabel("y")     
            plt.plot(vx, vy, 'go')
            plt.plot(vx, vy, 'k:', color='blue')
            #Clear list
            vx.clear()
            vy.clear()
            #Show graphic
            plt.show()
    
    return

#script mode
if __name__ == "__main__":
    main()