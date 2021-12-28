import math
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import sympy as sp
import pprint
import math
import time
import scipy.io as io
import pickle


class PlanetSystem():
    def __init__(self, planets, asteroids, spaceShip='NotEnoughGold'):
        self.planets = planets
        self.asteroids = asteroids
        if spaceShip != 'NotEnoughGold':
            self.spaceShip = spaceShip

    def AddNewPlanet(self, planet):
        self.planets.append(planet)

    def ReplaceSystem(self, X, Y, VX, VY, X_Sh=0, Y_Sh=0, VX_Sh=0, VY_Sh=0):
        n = len(self.planets)
        for planet, x, y, vx, vy in zip(self.planets, X[0:n], Y[0:n], VX[0:n], VY[0:n]):
            planet.Replace(x, y, vx, vy)
            planet.ReDraw()
        for asteroid, x, y, vx, vy in zip(self.asteroids, X[n:], Y[n:], VX[n:], VY[n:]):
            asteroid.Replace(x, y, vx, vy)
            asteroid.ReDraw()
        if self.spaceShip:
            self.spaceShip.Replace(X_Sh, Y_Sh, VX_Sh, VY_Sh)
            self.spaceShip.ReDraw()

    def Draw(self, axes):
        for planet in self.planets:
            planet.Draw(axes)
        for asteroid in self.asteroids:
            asteroid.Draw(axes)
        if (self.spaceShip):
            self.spaceShip.Draw(axes)

    def GetMoveEquations(self):
        n = len(self.planets)
        m = len(self.asteroids)
        _strX_planet = ''
        _strY_planet = ''
        _strVx_planet = ''
        _strVy_planet = ''

        for i in range(n):
            _strX_planet += f'x{i}, '
            _strY_planet += f'y{i}, '
            _strVx_planet += f'Vx{i}, '
            _strVy_planet += f'Vy{i}, '

        for i in range(m):
            _strX_planet += f'x{i + n}, '
            _strY_planet += f'y{i + n}, '
            _strVx_planet += f'Vx{i + n}, '
            _strVy_planet += f'Vy{i + n}, '

        x_planet = sp.symbols(_strX_planet)
        y_planet = sp.symbols(_strY_planet)
        vx_planet = sp.symbols(_strVx_planet)
        vy_planet = sp.symbols(_strVy_planet)

        dx_planet = [Vx for Vx in vx_planet]
        dy_planet = [Vy for Vy in vy_planet]

        fullSystem = self.planets + self.asteroids

        dvx_planet = [
            sum([
                (planet.m * (x - cur_x)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2) ** 3)
                for x, y, planet in zip(x_planet, y_planet, self.planets)
                if (x != cur_x)
            ])
            for cur_x, cur_y, current_planet in zip(x_planet, y_planet, fullSystem)
        ]
        dvy_planet = [
            sum([
                (planet.m * (y - cur_y)) / (sp.sqrt((x - cur_x) ** 2 + (y - cur_y) ** 2) ** 3)
                for x, y, planet in zip(x_planet, y_planet, self.planets)
                if (x != cur_x)
            ])
            for cur_x, cur_y, current_planet in zip(x_planet, y_planet, fullSystem)
        ]
        self.SpaceBodyMoveEquations = sp.lambdify([x_planet, y_planet, vx_planet, vy_planet],
                                                  [dx_planet, dy_planet, dvx_planet, dvy_planet])

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
                for x, y, planet in zip(x_planet, y_planet, fullSystem)
            ]) + F_dv / self.spaceShip.m * sp.cos(Alpha)

            DVY_Sh = sum([
                (planet.m * (y - Y_Sh)) / (sp.sqrt((x - X_Sh) ** 2 + (y - Y_Sh) ** 2) ** 3)
                for x, y, planet in zip(x_planet, y_planet, fullSystem)
            ]) + F_dv / self.spaceShip.m * sp.sin(Alpha)
        self.SpaceShipMoveEquations = sp.lambdify(
            [X_Sh, Y_Sh, VX_Sh, VY_Sh, x_planet, y_planet, vx_planet, vy_planet, F_dv, Alpha],
            [DX_Sh, DY_Sh, DVX_Sh, DVY_Sh])

    def GetStateVectors(self):
        X = np.zeros(len(self.planets) + len(self.asteroids))
        Y = np.zeros(len(self.planets) + len(self.asteroids))
        VX = np.zeros(len(self.planets) + len(self.asteroids))
        VY = np.zeros(len(self.planets) + len(self.asteroids))
        n = len(self.planets)
        for i in range(len(self.planets)):
            X[i] = self.planets[i].x
            Y[i] = self.planets[i].y
            VX[i] = self.planets[i].Vx
            VY[i] = self.planets[i].Vy

        for i in range(len(self.asteroids)):
            X[i + n] = self.asteroids[i].x
            Y[i + n] = self.asteroids[i].y
            VX[i + n] = self.asteroids[i].Vx
            VY[i + n] = self.asteroids[i].Vy

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


class Asteroid():
    def __init__(self, x0, y0, Vx0, Vy0, R, color):
        self.x0 = x0
        self.y0 = y0
        self.Vx0 = Vx0
        self.Vy0 = Vy0
        self.m = 0.00001
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

        self.SpaceShipX = self.R * np.array([1, 0.5, 0, -0.25, -0.5, -1, -0.6, -0.6, -1, -0.5, -0.25, 0, 0.5, 1])
        self.SpaceShipY = self.R * np.array(
            [0, 0.2, 0.25, 0.23, 0.5, 0.5, 0.2, -0.2, -0.5, -0.5, -0.23, -0.25, -0.2, 0])

        self.SpaceShipFlameX = self.R * np.array([0, -0.4, 0])
        self.SpaceShipFlameY = self.R * np.array([0.2, 0, -0.2])

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
        self.DrawedSpaceShip = axes.plot(self.x + self.SpaceShipX, self.y + self.SpaceShipY, color=self.color)[0]
        self.DrawedSpaceShipFlame = \
            axes.plot(self.x + self.SpaceShipFlameX, self.y + self.SpaceShipFlameY, color=self.color)[0]
        self.DrawedTrace = axes.plot(self.TraceX, self.TraceY, ':')[0]

    def ReDraw(self):
        self.DrawedSpaceShip.set_data(self.x + self.SpaceShipX, self.y + self.SpaceShipY)
        self.DrawedTrace.set_data(self.TraceX, self.TraceY)


def NewPoints(i):
    global t, dt, plSystem, X, Y, VX, VY, Dx, Dy, DVx, DVy, X_Sh, Y_Sh, VX_Sh, VY_Sh, Dx_Sh, Dy_Sh, DVx_Sh, DVy_Sh, F_dv, Alpha
    t += dt

    # Методом Рунге - Кутты
    Dx1, Dy1, DVx1, DVy1 = plSystem.SpaceBodyMoveEquations(X, Y, VX, VY)
    Dx1_Sh, Dy1_Sh, DVx1_Sh, DVy1_Sh = plSystem.SpaceShipMoveEquations(X_Sh, Y_Sh, VX_Sh, VY_Sh, X, Y, VX, VY, F_dv,
                                                                       Alpha)
    # print(Dx1, Dy1, DVx1, DVy1)
    # print(Dx1_Sh, Dy1_Sh, DVx1_Sh, DVy1_Sh)
    Dx1 = np.array(Dx1)
    Dy1 = np.array(Dy1)
    DVx1 = np.array(DVx1)
    DVy1 = np.array(DVy1)
    Dx1_Sh = np.array(Dx1_Sh)
    Dy1_Sh = np.array(Dy1_Sh)
    DVx1_Sh = np.array(DVx1_Sh)
    DVy1_Sh = np.array(DVy1_Sh)

    Dx2, Dy2, DVx2, DVy2 = plSystem.SpaceBodyMoveEquations(X + Dx1 / 2 * dt, Y + Dy1 / 2 * dt, VX + DVx1 / 2 * dt,
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

    Dx3, Dy3, DVx3, DVy3 = plSystem.SpaceBodyMoveEquations(X + Dx2 / 2 * dt, Y + Dy2 / 2 * dt, VX + DVx2 / 2 * dt,
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

    Dx4, Dy4, DVx4, DVy4 = plSystem.SpaceBodyMoveEquations(X + Dx3 * dt, Y + Dy3 * dt, VX + DVx3 * dt, VY + DVy3 * dt)
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

    plSystem.ReplaceSystem(X, Y, VX, VY, X_Sh, Y_Sh, VX_Sh, VY_Sh)

    drPlanets = [planet.DrawedPlanet for planet in (plSystem.planets + plSystem.asteroids)]
    drTraces = [planet.DrawedTrace for planet in (plSystem.planets + plSystem.asteroids)]

    return drPlanets + drTraces + [plSystem.spaceShip.DrawedSpaceShip] \
           + [plSystem.spaceShip.DrawedSpaceShipFlame] + [plSystem.spaceShip.DrawedTrace]# + [ax.plot(R * np.sin(np.linspace(0, 6.28, 200)), R * np.cos(np.linspace(0, 6.28, 200)))[0]]


if __name__ == '__main__':
    global t, dt, plSystem, X, Y, VX, VY, Dx, Dy, DVx, DVy, X_Sh, Y_Sh, VX_Sh, VY_Sh, Dx_Sh, Dy_Sh, DVx_Sh, DVy_Sh, F_dv, Alpha, ax
    np.seterr(all='warn', over='raise')
    mult = 1
    R = 10 * mult
    m1 = 1000 * mult
    m2 = 1 * mult
    v = 10
    pl1 = Planet(0, 0, 0, 0, m1, 1, 'red')
    pl2 = Planet(R, 0, 0, v, m2, 0.2, 'blue')
    pl3 = Planet(0, -10, 10, 0, 50, 0.3, 'black')
    # pl4 = Planet(-5, 0, 0, -10, 10, 0.4, 'green')
    rnd = np.zeros(20)  # np.random.rand(16)

    Our_Rocket = SpaceShip(100000, -5, 10, 0, 50, 0.3, 'black')

    # plSystem = PlanetSystem([pl1, pl2, pl3, pl4])
    '''as1 = Asteroid(rnd[0], -7+rnd[1], 10 +rnd[2], rnd[3], 0.3, 'black')
    as2 = Asteroid(rnd[4], 7+rnd[5], -10+rnd[6], rnd[7], 0.3, 'black')
    as3 = Asteroid(7+rnd[8], rnd[9], rnd[10], 10+rnd[1], 0.3, 'black')
    as4 = Asteroid(-7+rnd[12], rnd[13], rnd[14], -10+rnd[15], 0.3, 'black')

    plSystem = PlanetSystem([pl1], [as1, as2, as3, as4], Our_Rocket)'''
    alpha = m2 / (m1 + m2)
    x1 = R * (1 - (alpha / 3) ** (1 / 3))
    x2 = R * (1 + (alpha / 3) ** (1 / 3))
    x3 = - R * (1 + 5 / 12 * alpha)
    beta = (m1 - m2) / (m1 + m2)
    x4 = R / 2 * beta
    y4 = R / 2 * 3 ** (1 / 2)
    x5 = R / 2 * beta
    y5 = -R / 2 * 3 ** (1 / 2)
    speed_plus = 0.1
    as1 = Asteroid(rnd[0] + x1, rnd[1], rnd[2], rnd[3] + v * x1 / R + speed_plus, 0.1, 'black')
    as2 = Asteroid(rnd[4] + x2, rnd[5], rnd[6], rnd[7] + v * x2 / R + speed_plus, 0.1, 'black')
    as3 = Asteroid(rnd[8] + x3, rnd[9], rnd[10], rnd[11] + v * x3 / R + speed_plus, 0.1, 'black')
    as4 = Asteroid(x4 + rnd[12], y4 + rnd[13], rnd[14] - 2 * v * x4 / R * np.sqrt(3) / 2, rnd[15] + 2 * v * x4 / (R * 2), 0.1, 'black')
    as5 = Asteroid(x5 + rnd[16], y5 + rnd[17], rnd[18] + 2 * v * x5 * np.sqrt(3) / (2 * R), rnd[19] + 2 * v * x5 / (R * 2), 0.1, 'black')

    plSystem = PlanetSystem([pl1, pl2], [as1, as2, as3, as4, as5], Our_Rocket)

    io.savemat('FirstUniverse.mat', {'PS': plSystem})

    config_dictionary = {'PS': plSystem}

    print(2)

    with open('FirstUniverse', 'wb') as config_dictionary_file:
        pickle.dump(config_dictionary, config_dictionary_file)

    # D = io.loadmat('FirstUniverse.mat')
    # print(D['PS'][0][0])
    # print(type(D['PS'][0][0]))

    # plSystem = D['PS']
    # plSystem.planets[1]

    # with open('FirstUniverse', 'rb') as config_dictionary_file:
    #     # Step 3
    #     config_dictionary = pickle.load(config_dictionary_file)
    #
    #     # After config_dictionary is read from file
    #     print(config_dictionary)
    #
    # plSystem=config_dictionary['PS']
    # print(type(plSystem))

    plSystem.GetMoveEquations()

    X, Y, VX, VY = plSystem.GetStateVectors()
    X_Sh = plSystem.spaceShip.x
    Y_Sh = plSystem.spaceShip.y
    VX_Sh = plSystem.spaceShip.Vx
    VY_Sh = plSystem.spaceShip.Vy

    F_dv = 70
    Alpha = -0.4

    t, dt = 0.0, 0.01

    fig = plt.figure(figsize=[13, 9])
    ax = fig.add_subplot(1, 1, 1)

    # ax = plt.gca()
    xmin = -15 * mult
    xmax = 15 * mult
    ymin = -15 * mult
    ymax = 15 * mult

    ax.axis('equal')
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    plSystem.Draw(ax)

    animation = FuncAnimation(fig, NewPoints, blit=True)

    plt.show()
