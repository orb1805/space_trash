
import numpy as np
import pprint
import pickle
import FormOfSpaceObjects
import math
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
import os


class PlanetSystem():
    def __init__(self, planets, spaceShip='NotEnoughGold'):
        self.planets = planets
        if spaceShip != 'NotEnoughGold':
            self.spaceShip = spaceShip
        else:
            self.spaceShip = []


    def AddNewPlanet(self, planet):
        self.planets.append(planet)

    def AddSpaceShip(self, spaceShip):
        self.spaceShip = spaceShip

    def ReplaceSystem(self, X, Y, VX, VY, X_Sh=0, Y_Sh=0, VX_Sh=0, VY_Sh=0):
        for planet, x, y, vx, vy in zip(self.planets, X, Y, VX, VY):
            planet.Replace(x, y, vx, vy)
            planet.ReDraw()
        if (self.spaceShip):
            self.spaceShip.Replace(X_Sh, Y_Sh, VX_Sh, VY_Sh)
            self.spaceShip.ReDraw()

    def Draw(self, axes):
        for planet in self.planets:
            planet.Draw(axes)
        if (self.spaceShip):
            self.spaceShip.Draw(axes)

    def GetMoveEquations(self):
        n = len(self.planets)
        _strX = ''
        _strY = ''
        _strVx = ''
        _strVy = ''
        for i in range(n):
            _strX += f'x{i}, '
            _strY += f'y{i}, '
            _strVx += f'Vx{i}, '
            _strVy += f'Vy{i}, '

        X = sp.symbols(_strX)
        Y = sp.symbols(_strY)
        VX = sp.symbols(_strVx)
        VY = sp.symbols(_strVy)

        DX = [Vx for Vx in VX]
        DY = [Vy for Vy in VY]

        DVX = [
            sum([   
                    (planet.m * (x - cur_x))/(sp.sqrt((x - cur_x)**2 + (y - cur_y)**2)**3)
                    for x, y, planet in zip(X, Y, self.planets)
                    if(x != cur_x)
                ])
            for cur_x, cur_y, current_planet in zip(X, Y, self.planets)
        ]

        DVY = [
            sum([   
                    (planet.m * (y - cur_y))/(sp.sqrt((x - cur_x)**2 + (y - cur_y)**2)**3)
                    for x, y, planet in zip(X, Y, self.planets)
                    if(x != cur_x)
                ])
            for cur_x, cur_y, current_planet in zip(X, Y, self.planets)
        ]
        self.SpaceBodyMoveEquations = sp.lambdify([X, Y, VX, VY], [DX, DY, DVX, DVY])

        if (self.spaceShip):
            X_Sh = sp.symbols('x_Sh')
            Y_Sh = sp.symbols('y_Sh')
            VX_Sh = sp.symbols('Vx_Sh')
            VY_Sh = sp.symbols('Vy_Sh')

            F_dv = sp.symbols('f_dv')
            Alpha = sp.symbols('alpha')

            DX_Sh = VX_Sh
            DY_Sh = VY_Sh

            DVX_Sh = sum([
                        (planet.m * (x - X_Sh))/(sp.sqrt((x - X_Sh)**2 + (y - Y_Sh)**2)**3)
                        for x, y, planet in zip(X, Y, self.planets)
                    ]) + F_dv/self.spaceShip.m*sp.cos(Alpha)

            DVY_Sh = sum([
                        (planet.m * (y - Y_Sh))/(sp.sqrt((x - X_Sh)**2 + (y - Y_Sh)**2)**3)
                        for x, y, planet in zip(X, Y, self.planets)
                    ]) + F_dv/self.spaceShip.m*sp.sin(Alpha)
        self.SpaceShipMoveEquations = sp.lambdify([X_Sh, Y_Sh, VX_Sh, VY_Sh, X, Y, VX, VY, F_dv, Alpha],
                                                  [DX_Sh, DY_Sh, DVX_Sh, DVY_Sh])



    def GetStateVectors(self):
        X = np.zeros(len(self.planets))
        Y = np.zeros(len(self.planets))
        VX = np.zeros(len(self.planets))
        VY = np.zeros(len(self.planets))

        for i in range(len(self.planets)):
            X[i] = self.planets[i].x
            Y[i] = self.planets[i].y
            VX[i] = self.planets[i].Vx
            VY[i] = self.planets[i].Vy

        return X, Y, VX, VY

class Planet():
    def __init__(self, x0, y0, Vx0, Vy0, m, R, color):
        self.x0 = x0
        self.y0 = y0
        self.Vx0 = Vx0
        self.Vy0 = Vy0 
        self.m = m
        self.R = R  
        self.color = color

        self.x = x0
        self.y = y0
        self.Vx = Vx0
        self.Vy = Vy0

        phi = np.linspace(0, 6.28, 20)
        self.PlanetX = self.R*np.sin(phi)
        self.PlanetY = self.R*np.cos(phi)

        self.TraceX = np.array([self.x])
        self.TraceY = np.array([self.y])

    def Replace(self, x, y, vx, vy):
        self.x = x
        self.y = y
        self.Vx = vx
        self.Vy = vy

        self.TraceX = np.append(self.TraceX, x)
        self.TraceY = np.append(self.TraceY, y)

    def Draw(self, axes):
        self.DrawedPlanet = axes.plot(self.x + self.PlanetX, self.y + self.PlanetY, color=self.color)[0]
        self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, ':')[0]
        

    def ReDraw(self):
        self.DrawedPlanet.set_data(self.x + self.PlanetX, self.y + self.PlanetY)
        self.DrawedTrace.set_data(self.TraceX, self.TraceY)

class SpaceShip():
    def __init__(self, x0, y0, Vx0, Vy0, m, R, color, flame_color, Phi0, F_dv):
        self.x0 = x0
        self.y0 = y0
        self.Vx0 = Vx0
        self.Vy0 = Vy0
        self.Phi0 = Phi0
        self.m = m
        self.R = R
        self.color = color
        self.flame_color = flame_color
        self.F_dv = F_dv

        self.x = x0
        self.y = y0
        self.Vx = Vx0
        self.Vy = Vy0
        self.Phi = Phi0

        self.SpaceShipX = self.R * np.array([1, 0.5, 0, -0.25, -0.5, -1, -0.6, -0.6, -1, -0.5, -0.25, 0, 0.5, 1])
        self.SpaceShipY = self.R * np.array([0, 0.2, 0.25, 0.23, 0.5, 0.5, 0.2, -0.2, -0.5, -0.5, -0.23, -0.25, -0.2, 0])

        self.SpaceShipFlameX = self.R * np.array([0, -3.5, -3, -4, -3, -3.5, 0])
        self.SpaceShipFlameY = self.R * np.array([0.2, 0.15, 0.1, 0, -0.1, -0.15, -0.2])

        self.TraceX = np.array([self.x])
        self.TraceY = np.array([self.y])

    def Replace(self, x, y, vx, vy):
        self.x = x
        self.y = y
        self.Vx = vx
        self.Vy = vy

        self.TraceX = np.append(self.TraceX, x)
        self.TraceY = np.append(self.TraceY, y)

    def Draw(self, axes):
        RShipX, RShipY = Rot2D(self.SpaceShipX, self.SpaceShipY, self.Phi)
        self.DrawedSpaceShip = axes.plot(self.x + RShipX, self.y + RShipY, color=self.color)[0]
        SmX = self.SpaceShipFlameX-0.6
        RFlameX, RFlameY = Rot2D(SmX, self.SpaceShipFlameY, self.Phi)
        self.DrawedSpaceShipFlame = axes.plot(self.x + RFlameX, self.y + RFlameY, color=self.color)[0]
        self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, ':')[0]

    def ReDraw(self):
        RShipX, RShipY = Rot2D(self.SpaceShipX, self.SpaceShipY, self.Phi)
        self.DrawedSpaceShip.set_data(self.x + RShipX, self.y + RShipY)
        SmX = self.SpaceShipFlameX - 0.6
        RFlameX, RFlameY = Rot2D(SmX, self.SpaceShipFlameY, self.Phi)
        self.DrawedSpaceShipFlame.set_data(self.x + RFlameX, self.y + RFlameY)
        self.DrawedTrace.set_data(self.TraceX, self.TraceY)

class SpaceWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

class SpaceWidget(QMainWindow, FormOfSpaceObjects.Ui_MainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)

        self.setWindowTitle("Творение")

        self.listUniverseObjectsWidget.currentItemChanged.connect(self.ChoseRow)
        self.AddPlanet.clicked.connect(self.AddPlanetFunction)
        self.AddStar.clicked.connect(self.AddStarFunction)
        self.AddShip.clicked.connect(self.AddShipFunction)
        self.ShowSystem.clicked.connect(self.ShowSystemFunction)
        self.saveDataButton.clicked.connect(self.SaveChanges)
        self.deleteButton.clicked.connect(self.DeleteItemFunction)

    def DeleteItemFunction(self):
        currentIndex = self.listUniverseObjectsWidget.currentRow()
        item = self.listUniverseObjectsWidget.item(currentIndex).text()

        FileName = self.Fname_field.text() + '.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']

            if('Космический' in item.split(' ')[1]):
                plSystem.spaceShip = []
            elif('Планета' in item.split(' ')[1]):
                plSystem.planets.pop(currentIndex)
               
            PS_dictionary = {'PS': plSystem}
            with open(FileName, 'wb') as PlanetSystemFile:
                pickle.dump(PS_dictionary, PlanetSystemFile)

            self.ShowSystemFunction()

    def SaveChanges(self):
        currentIndex = self.listUniverseObjectsWidget.currentRow()
        item = self.listUniverseObjectsWidget.item(currentIndex).text()

        FileName = self.Fname_field.text() + '.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']

            print('!!!!!', item.split(' '))
            if('Космический' in item.split(' ')[1]):
                plSystem.spaceShip.x0 = float(self.x0_field.text())
                plSystem.spaceShip.y0 = float(self.y0_field.text())
                plSystem.spaceShip.Vx0 = float(self.Vx0_field.text())
                plSystem.spaceShip.Vy0 = float(self.Vy0_field.text())
                plSystem.spaceShip.m = float(self.M_field.text())
                plSystem.spaceShip.R = float(self.R_field.text())
                plSystem.spaceShip.color = [float(i) for i in self.Color_field.text().split(', ')]
                plSystem.spaceShip.flame_color = [float(i) for i in self.Flamecolor_field.text().split(', ')]

            else:
                plSystem.planets[currentIndex].x0 = float(self.x0_field.text())
                plSystem.planets[currentIndex].y0 = float(self.y0_field.text())
                plSystem.planets[currentIndex].Vx0 = float(self.Vx0_field.text())
                plSystem.planets[currentIndex].Vy0 = float(self.Vy0_field.text())
                plSystem.planets[currentIndex].m = float(self.M_field.text())
                plSystem.planets[currentIndex].R = float(self.R_field.text())
                plSystem.planets[currentIndex].color = [float(i) for i in self.Color_field.text().split(', ')]


            PS_dictionary = {'PS': plSystem}
            with open(FileName, 'wb') as PlanetSystemFile:
                pickle.dump(PS_dictionary, PlanetSystemFile)

            i = 1 
            self.listUniverseObjectsWidget.clear()
            objects = []
            for planet in plSystem.planets:
                str = f' \nПланета {i}: x0={planet.x0}, y0={planet.y0}, Vx0={planet.Vx0}, Vy0={planet.Vy0},  \nm={planet.m}, R={planet.R}, Цвет: {planet.color}\n'
                i += 1
                objects.append(str)
            if plSystem.spaceShip != []:
                str = f' \nКосмический корабль: \nx0={plSystem.spaceShip.x0}, y0={plSystem.spaceShip.y0}, ' \
                          f'Vx0={plSystem.spaceShip.Vx0}, Vy0={plSystem.spaceShip.Vy0},  \nm={plSystem.spaceShip.m}, ' \
                          f'R={plSystem.spaceShip.R}, Цвет: {plSystem.spaceShip.color}\n'
                objects.append(str)

            self.listUniverseObjectsWidget.addItems(objects)

    def ChoseRow(self):
        try:
            currentIndex = self.listUniverseObjectsWidget.currentRow()

            item = self.listUniverseObjectsWidget.item(currentIndex).text()

            data = item.split(':')[1].split(',')
            for i in data:
                if(' x0' in i or '\nx0' in i):
                    x0 = i.split('=')[1]
                if(' y0' in i):
                    y0 = i.split('=')[1]
                if('Vx0' in i):
                    Vx0 = i.split('=')[1]
                if('Vy0' in i):
                    Vy0 = i.split('=')[1]
                if('m' in i): m = i.split('=')[1]
                if('R' in i):
                    R = i.split('=')[1]

            color = item.split(':')[2][2:-2:]

            self.M_field.setText(m)
            self.x0_field.setText(x0)
            self.y0_field.setText(y0)
            self.Vx0_field.setText(Vx0)
            self.Vy0_field.setText(Vy0)
            self.R_field.setText(R)
            self.Color_field.setText(color)

        except: 
            print('Пустой лист')
        
    def AddPlanetFunction(self):

        m = float(self.M_field.text())
        R = float(self.R_field.text())

        Color = [float(n) for n in self.Color_field.text().split(', ')]

        x = float(self.x0_field.text())
        y = float(self.y0_field.text())
        Vx = float(self.Vx0_field.text())
        Vy = float(self.Vy0_field.text())

        P = Planet(x, y, Vx, Vy, m, R, Color)

        FileName = self.Fname_field.text()+'.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']
            plSystem.AddNewPlanet(P)
        else:
            plSystem = PlanetSystem([P])
        PS_dictionary = {'PS': plSystem}
        with open(FileName, 'wb') as PlanetSystemFile:
            pickle.dump(PS_dictionary, PlanetSystemFile)
        self.ShowSystemFunction()

    def AddStarFunction(self):
        FileName = self.Fname_field.text()

    def AddShipFunction(self):

        m = float(self.M_field.text())
        R = float(self.R_field.text())

        Color = [float(n) for n in self.Color_field.text().split(', ')]

        x = float(self.x0_field.text())
        y = float(self.y0_field.text())
        Vx = float(self.Vx0_field.text())
        Vy = float(self.Vy0_field.text())

        Phi = float(self.Phi_field.text())
        F_dv = float(self.F_dv_field.text())
        Flame_color = [float(i) for i in self.Flamecolor_field.text().split(', ')]

        Sh = SpaceShip(x, y, Vx, Vy, m, R, Color, Flame_color, Phi, F_dv)

        FileName = self.Fname_field.text()+'.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']
            plSystem.AddSpaceShip(Sh)
        else:
            plSystem = PlanetSystem([], Sh)
        PS_dictionary = {'PS': plSystem}
        with open(FileName, 'wb') as PlanetSystemFile:
            pickle.dump(PS_dictionary, PlanetSystemFile)
        self.ShowSystemFunction()

    def ShowSystemFunction(self):
        FileName = self.Fname_field.text() + '.universe'
        Filepath = f'{os.path.abspath(os.path.dirname(__file__))}\{FileName}'
        if os.path.exists(Filepath):
            with open(FileName, 'rb') as PlanetSystemFile:
                PS_dictionary = pickle.load(PlanetSystemFile)
            plSystem = PS_dictionary['PS']
            text = 'Планеты:'
            i = 1
            for planet in plSystem.planets:
                str = f' \nПланета {i}: x0={planet.x0}, y0={planet.y0}, Vx0={planet.Vx0}, Vy0={planet.Vy0},  \nm={planet.m}, R={planet.R}, Цвет: {planet.color}\n'
                i += 1
                text = text + str
            if plSystem.spaceShip != []:
                str = f' \nКосмический корабль: \nx0={plSystem.spaceShip.x0}, y0={plSystem.spaceShip.y0}, ' \
                      f'Vx0={plSystem.spaceShip.Vx0}, Vy0={plSystem.spaceShip.Vy0},  \nm={plSystem.spaceShip.m}, ' \
                      f'R={plSystem.spaceShip.R}, Цвет: {plSystem.spaceShip.color}\n'
                text = text + str
            self.UniverseBrowser.setText(text)
        else:
            self.UniverseBrowser.setText('Такого файла не существует')

        i = 1 
        self.listUniverseObjectsWidget.clear()

        objects = []
        for planet in plSystem.planets:
            str = f' \nПланета {i}: x0={planet.x0}, y0={planet.y0}, Vx0={planet.Vx0}, Vy0={planet.Vy0},  \nm={planet.m}, R={planet.R}, Цвет: {planet.color}\n'
            i += 1
            objects.append(str)
        if plSystem.spaceShip != []:
            str = f' \nКосмический корабль: \nx0={plSystem.spaceShip.x0}, y0={plSystem.spaceShip.y0}, ' \
                      f'Vx0={plSystem.spaceShip.Vx0}, Vy0={plSystem.spaceShip.Vy0},  \nm={plSystem.spaceShip.m}, ' \
                      f'R={plSystem.spaceShip.R}, Цвет: {plSystem.spaceShip.color}\n'
            objects.append(str)

        self.listUniverseObjectsWidget.addItems(objects)


app = QApplication([])
window = SpaceWidget()
window.show()
app.exec_()