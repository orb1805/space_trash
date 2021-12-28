from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import FormOfUniverseSimulation
import numpy as np
from matplotlib.animation import FuncAnimation
import sympy as sp
import pickle


class PlanetSystem():
    def __init__(self, planets, spaceShip='NotEnoughGold'):
        self.planets = planets
        if spaceShip != 'NotEnoughGold':
            self.spaceShip = spaceShip

    def AddNewPlanet(self, planet):
        self.planets.append(planet)

    def ReplaceSystem(self, X, Y, VX, VY, X_Sh=0, Y_Sh=0, VX_Sh=0, VY_Sh=0, Phi_Sh=0, F_curr=0):
        for planet, x, y, vx, vy in zip(self.planets, X, Y, VX, VY):
            planet.Replace(x, y, vx, vy)
            planet.ReDraw()
        if (self.spaceShip):
            self.spaceShip.Replace(X_Sh, Y_Sh, VX_Sh, VY_Sh, Phi_Sh, F_curr)
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
                (planet.m * (x - cur_x)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2) ** 3)
                for x, y, planet in zip(X, Y, self.planets)
                if (x != cur_x)
            ])
            for cur_x, cur_y, current_planet in zip(X, Y, self.planets)
        ]
        DVY = [
            sum([
                (planet.m * (y - cur_y)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2) ** 3)
                for x, y, planet in zip(X, Y, self.planets)
                if (x != cur_x)
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
                (planet.m * (x - X_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2) ** 3)
                for x, y, planet in zip(X, Y, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.cos(Alpha)
            DVY_Sh = sum([
                (planet.m * (y - Y_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2) ** 3)
                for x, y, planet in zip(X, Y, self.planets)
            ]) + F_dv / self.spaceShip.m * sp.sin(Alpha)
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
        self.PlanetX = self.R * np.sin(phi)
        self.PlanetY = self.R * np.cos(phi)

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
    def __init__(self, x0, y0, Vx0, Vy0, m, R, color, F_max):
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

        self.phi = 0
        self.F_dv = F_max
        self.F_curr = 0

        self.SpaceShipX = self.R * np.array([1, 0.5, 0, -0.25, -0.5, -1, -0.6, -0.6, -1, -0.5, -0.25, 0, 0.5, 1])
        self.SpaceShipY = self.R * np.array(
            [0, 0.2, 0.25, 0.23, 0.5, 0.5, 0.2, -0.2, -0.5, -0.5, -0.23, -0.25, -0.2, 0])

        self.SpaceShipFlameX = self.R * np.array([0, -3.5, -3, -4, -3, -3.5, 0])
        self.SpaceShipFlameY = self.R * np.array([0.2, 0.15, 0.1, 0, -0.1, -0.15, -0.2])

        self.TraceX = np.array([self.x])
        self.TraceY = np.array([self.y])

    def Replace(self, x, y, vx, vy, phi, F_curr):
        self.x = x
        self.y = y
        self.Vx = vx
        self.Vy = vy
        self.phi = phi
        self.F_curr = F_curr

        self.TraceX = np.append(self.TraceX, x)
        self.TraceY = np.append(self.TraceY, y)

    def Draw(self, axes):
        self.DrawedSpaceShip = axes.plot(self.x + self.SpaceShipX, self.y + self.SpaceShipY, color=self.color)[0]
        self.DrawedSpaceShipFlame = \
            axes.plot(self.x + self.SpaceShipFlameX, self.y + self.SpaceShipFlameY, color=self.flame_color)[0]
        self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, ':')[0]

    def ReDraw(self):
        RotSpaceShipX, RotSpaceShipY = Rot2D(self.SpaceShipX, self.SpaceShipY, self.phi)
        self.DrawedSpaceShip.set_data(self.x + RotSpaceShipX, self.y + RotSpaceShipY)

        PFlameX = self.F_curr / self.F_dv * self.SpaceShipFlameX + self.SpaceShipX[int(len(self.SpaceShipX) / 2)]
        RotSpaceShipFlameX, RotSpaceShipFlameY = Rot2D(PFlameX, self.SpaceShipFlameY, self.phi)
        self.DrawedSpaceShipFlame.set_data(self.x + RotSpaceShipFlameX, self.y + RotSpaceShipFlameY)
        self.DrawedTrace.set_data(self.TraceX, self.TraceY)


def Rot2D(X, Y, phi):
    RotX = X * np.cos(phi) - Y * np.sin(phi)
    RotY = X * np.sin(phi) + Y * np.cos(phi)
    return RotX, RotY


def DrawTheSpace(axes):
    print(4)
    axes.fill([-100, 100, 100, - 100], [-100, - 100, 100, 100], 'black')
    print(5)
    nstars = 5000
    xstars = 200 * np.random.random(nstars) - 100
    ystars = 200 * np.random.random(nstars) - 100
    brstars = 0.5 + np.random.random(nstars) / 2
    sizestars = np.random.random(nstars)
    print(6)
    for i in range(nstars):
        axes.plot(xstars[i], ystars[i], marker='o', markersize=sizestars[i], color=[brstars[i], brstars[i], brstars[i]])


global t


class SpaceWidget(QMainWindow, FormOfUniverseSimulation.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)

        self.setWindowTitle("В одной далёкой галактике")

        self.StartButton.clicked.connect(self.HereAreWeGo)
        # self.StopButton.clicked.connect(self.HereAreWeGo)

        # self.F_Bar.valueChanged.connect(self.F_valuechange)
        # self.Angle_Bar.valueChanged.connect(self.Angle_valuechange)
        self.addToolBar(NavigationToolbar(self.SpWidget.canvas, self))

    def HereAreWeGo(self):
        global t, dt, plSystem, X, Y, VX, VY, Dx, Dy, DVx, DVy, X_Sh, Y_Sh, VX_Sh, VY_Sh, Dx_Sh, Dy_Sh, DVx_Sh, DVy_Sh, F_max, F_dv, Alpha
        print(0)
        #     Параметры массы
        FileName = self.FileName_field.text() + '.universe'
        dt = float(self.TStep_field.text())
        K = float(self.K_field.text())

        print(1)

        with open(FileName, 'rb') as config_dictionary_file:
            # Step 3
            config_dictionary = pickle.load(config_dictionary_file)

        plSystem = config_dictionary['PS']
        print(1.5)
        plSystem.GetMoveEquations()

        print(2)
        print(plSystem.planets)
        print(plSystem.spaceShip)

        X, Y, VX, VY = plSystem.GetStateVectors()
        if (plSystem.spaceShip):
            X_Sh = plSystem.spaceShip.x
            Y_Sh = plSystem.spaceShip.y
            VX_Sh = plSystem.spaceShip.Vx
            VY_Sh = plSystem.spaceShip.Vy
            F_max = plSystem.spaceShip.F_dv
            F_dv = self.F_Bar.value() * F_max / 100
            Alpha = self.Angle_Bar.value() / 360 * 6.28 + 1.57
        else:
            X_Sh = 0
            Y_Sh = 0
            VX_Sh = 0
            VY_Sh = 0
            F_max = 0
            F_dv = 0
            Alpha = 0

        self.SpWidget.canvas.axes.clear()
        # self.SpWidget.canvas.axes.grid(True)
        self.SpWidget.canvas.axes.axis('scaled')
        Side = K * 20

        self.SpWidget.canvas.axes.set(xlim=[-2 * Side, 2 * Side], ylim=[-Side, Side])
        self.SpWidget.canvas.axes.set_title('Это космос')

        t = 0.0
        print(3)
        DrawTheSpace(self.SpWidget.canvas.axes)
        print(7)
        plSystem.Draw(self.SpWidget.canvas.axes)
        self.SpWidget.canvas.show()

        def NewPoints(i):
            global t, dt, plSystem, X, Y, VX, VY, Dx, Dy, DVx, DVy, X_Sh, Y_Sh, VX_Sh, VY_Sh, Dx_Sh, Dy_Sh, DVx_Sh, DVy_Sh, F_max, F_dv, Alpha
            t += dt
            F_dv = self.F_Bar.value() * F_max / 100
            Alpha = -self.Angle_Bar.value() / 360 * 6.28 + 1.57
            # Методом Эйлера
            # Dx, Dy, DVx, DVy = plSystem.SpaceBodyMoveEquations(X, Y, VX, VY)
            # Dx = np.array(Dx)
            # Dy = np.array(Dy)
            # DVx = np.array(DVx)
            # DVy = np.array(DVy)
            # X = X + dt * Dx
            # Y = Y + dt * Dy
            # VX = VX + dt * DVx
            # VY = VY + dt * DVy

            # Методом Рунге - Кутты
            Dx1, Dy1, DVx1, DVy1 = plSystem.SpaceBodyMoveEquations(X, Y, VX, VY)
            Dx1_Sh, Dy1_Sh, DVx1_Sh, DVy1_Sh = plSystem.SpaceShipMoveEquations(X_Sh, Y_Sh, VX_Sh, VY_Sh, X, Y, VX, VY,
                                                                               F_dv,
                                                                               Alpha)

            # print(Dx1_Sh, Dy1_Sh, DVx1_Sh, DVy1_Sh)
            Dx1 = np.array(Dx1)
            Dy1 = np.array(Dy1)
            DVx1 = np.array(DVx1)
            DVy1 = np.array(DVy1)
            Dx1_Sh = np.array(Dx1_Sh)
            Dy1_Sh = np.array(Dy1_Sh)
            DVx1_Sh = np.array(DVx1_Sh)
            DVy1_Sh = np.array(DVy1_Sh)

            Dx2, Dy2, DVx2, DVy2 = plSystem.SpaceBodyMoveEquations(X + Dx1 / 2 * dt, Y + Dy1 / 2 * dt,
                                                                   VX + DVx1 / 2 * dt,
                                                                   VY + DVy1 / 2 * dt)
            Dx2_Sh, Dy2_Sh, DVx2_Sh, DVy2_Sh = plSystem.SpaceShipMoveEquations(
                X_Sh + Dx1_Sh / 2 * dt, Y_Sh + Dy1_Sh / 2 * dt, VX_Sh + DVx1_Sh / 2 * dt, VY_Sh + DVy1_Sh / 2 * dt,
                X + Dx1 / 2 * dt, Y + Dy1 / 2 * dt, VX + DVx1 / 2 * dt, VY + DVy1 / 2 * dt, F_dv, Alpha)

            Dx2 = np.array(Dx2)
            Dy2 = np.array(Dy2)
            DVx2 = np.array(DVx2)
            DVy2 = np.array(DVy2)
            Dx2_Sh = np.array(Dx2_Sh)
            Dy2_Sh = np.array(Dy2_Sh)
            DVx2_Sh = np.array(DVx2_Sh)
            DVy2_Sh = np.array(DVy2_Sh)

            Dx3, Dy3, DVx3, DVy3 = plSystem.SpaceBodyMoveEquations(X + Dx2 / 2 * dt, Y + Dy2 / 2 * dt,
                                                                   VX + DVx2 / 2 * dt,
                                                                   VY + DVy2 / 2 * dt)
            Dx3_Sh, Dy3_Sh, DVx3_Sh, DVy3_Sh = plSystem.SpaceShipMoveEquations(
                X_Sh + Dx2_Sh / 2 * dt, Y_Sh + Dy2_Sh / 2 * dt, VX_Sh + DVx2_Sh / 2 * dt, VY_Sh + DVy2_Sh / 2 * dt,
                X + Dx2 / 2 * dt, Y + Dy2 / 2 * dt, VX + DVx2 / 2 * dt, VY + DVy2 / 2 * dt, F_dv, Alpha)

            Dx3 = np.array(Dx3)
            Dy3 = np.array(Dy3)
            DVx3 = np.array(DVx3)
            DVy3 = np.array(DVy3)
            Dx3_Sh = np.array(Dx3_Sh)
            Dy3_Sh = np.array(Dy3_Sh)
            DVx3_Sh = np.array(DVx3_Sh)
            DVy3_Sh = np.array(DVy3_Sh)

            Dx4, Dy4, DVx4, DVy4 = plSystem.SpaceBodyMoveEquations(X + Dx3 * dt, Y + Dy3 * dt, VX + DVx3 * dt,
                                                                   VY + DVy3 * dt)
            Dx4_Sh, Dy4_Sh, DVx4_Sh, DVy4_Sh = plSystem.SpaceShipMoveEquations(
                X_Sh + Dx3_Sh * dt, Y_Sh + Dy3_Sh * dt, VX_Sh + DVx3_Sh * dt, VY_Sh + DVy3_Sh * dt,
                X + Dx3 * dt, Y + Dy3 * dt, VX + DVx3 * dt, VY + DVy3 * dt, F_dv, Alpha)

            Dx4 = np.array(Dx4)
            Dy4 = np.array(Dy4)
            DVx4 = np.array(DVx4)
            DVy4 = np.array(DVy4)
            Dx4_Sh = np.array(Dx4_Sh)
            Dy4_Sh = np.array(Dy4_Sh)
            DVx4_Sh = np.array(DVx4_Sh)
            DVy4_Sh = np.array(DVy4_Sh)

            X = X + dt / 6 * (Dx1 + 2 * Dx2 + 2 * Dx3 + Dx4)
            Y = Y + dt / 6 * (Dy1 + 2 * Dy2 + 2 * Dy3 + Dy4)
            VX = VX + dt / 6 * (DVx1 + 2 * DVx2 + 2 * DVx3 + DVx4)
            VY = VY + dt / 6 * (DVy1 + 2 * DVy2 + 2 * DVy3 + DVy4)
            X_Sh = X_Sh + dt / 6 * (Dx1_Sh + 2 * Dx2_Sh + 2 * Dx3_Sh + Dx4_Sh)
            Y_Sh = Y_Sh + dt / 6 * (Dy1_Sh + 2 * Dy2_Sh + 2 * Dy3_Sh + Dy4_Sh)
            VX_Sh = VX_Sh + dt / 6 * (DVx1_Sh + 2 * DVx2_Sh + 2 * DVx3_Sh + DVx4_Sh)
            VY_Sh = VY_Sh + dt / 6 * (DVy1_Sh + 2 * DVy2_Sh + 2 * DVy3_Sh + DVy4_Sh)
            # print(X_Sh, Y_Sh, VX_Sh, VY_Sh)

            plSystem.ReplaceSystem(X, Y, VX, VY, X_Sh, Y_Sh, VX_Sh, VY_Sh, Alpha, F_dv)

            drPlanets = [planet.DrawedPlanet for planet in plSystem.planets]
            drTraces = [planet.DrawedTrace for planet in plSystem.planets]
            self.SpWidget.canvas.axes.axis('scaled')
            self.SpWidget.canvas.axes.set(xlim=[-2 * Side + X_Sh, 2 * Side + X_Sh], ylim=[-Side + Y_Sh, Side + Y_Sh])

            return drPlanets + drTraces + [plSystem.spaceShip.DrawedSpaceShip] \
                   + [plSystem.spaceShip.DrawedSpaceShipFlame] + [plSystem.spaceShip.DrawedTrace]

        fig = self.SpWidget.canvas.figure
        animation = FuncAnimation(fig, NewPoints, interval=dt * 1000, blit=True)
        self.SpWidget.canvas.draw()


app = QApplication([])
window = SpaceWidget()
window.show()
app.exec_()
