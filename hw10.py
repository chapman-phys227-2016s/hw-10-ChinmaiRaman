#! /usr/bin/env python

"""
File: hw10.py
Copyright (c) 2016 Chinmai Raman
License: MIT
Course: PHYS227
Assignment: Homework 10
Date: May 22, 2016
Email: raman105@mail.chapman.edu
Name: Chinmai Raman
Description: Solving ODEs using various methods
"""
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

def euler(u0, v0, T, n):
    t = np.zeros(n + 1)
    u = np.zeros(n + 1)
    v = np.zeros(n + 1)
    u[0] = u0
    v[0] = v0
    dt = T / float(n)
    for k in range(n):
        t[k + 1] = t[k] + dt
        u[k + 1] = u[k] + dt * v[k]
        v[k + 1] = v[k] - dt * u[k]
    return t, u, v

def heun(u0, v0, T, n):
    t = np.zeros(n + 1)
    u = np.zeros(n + 1)
    v = np.zeros(n + 1)
    u[0] = u0
    v[0] = v0
    dt = T / float(n)
    for k in range(n):
        t[k + 1] = t[k] + dt
        u_star = u[k] + dt * v[k]
        v_star = v[k] - dt * u[k]
        u[k + 1] = u[k] + 0.5 * dt * v[k] + 0.5 * dt * v_star
        v[k + 1] = v[k] - 0.5 * dt * u[k] - 0.5 * dt * u_star
    return t, u, v

def rungekutta2(u0, v0, T, n):
    t = np.zeros(n + 1)
    u = np.zeros(n + 1)
    v = np.zeros(n + 1)
    u[0] = u0
    v[0] = v0
    dt = T / float(n)
    for k in range(n):
        t[k + 1] = t[k] + dt
        K1_u = dt *  v[k]
        K1_v = -dt * u[k]
        u[k + 1] = u[k] + dt * (v[k] + 0.5 * K1_u)
        v[k + 1] = v[k] + dt * (-u[k] + 0.5 * K1_v)
    return t, u, v
    
def rungekutta4(u0, v0, T, n):
    t = np.zeros(n + 1)
    u = np.zeros(n + 1)
    v = np.zeros(n + 1)
    u[0] = u0
    v[0] = v0
    dt = T / float(n)
    for k in range(n):
        t[k + 1] = t[k] + dt
        K1_u = dt * v[k]
        K1_v = -dt * u[k]
        K2_u = dt * (v[k] + 0.5 * K1_u)
        K2_v = dt * (-u[k] + 0.5 * K1_v)
        K3_u = dt * (v[k] + 0.5 * K2_u)
        K3_v = dt * (-u[k] + 0.5 * K2_v)
        K4_u = dt * (v[k] + K3_u)
        K4_v = dt * (-u[k] + K3_v)
        u[k + 1] = u[k] + (K1_u + 2 * K2_u + 2 * K3_u + K4_u) / 6.0
        v[k + 1] = v[k] + (K1_v + 2 * K2_v + 2 * K3_v + K4_v) / 6.0
    return t, u, v
        
def plot(f, n):
    t, u, v = eval(f)(1, 0, 10 * np.pi, n)
    fig = plt.figure(1)
    x = np.linspace(0, 10 * np.pi, n)
    y_sin = -1 * np.sin(x)
    y_cos = np.cos(x)
    plt.plot(x, y_sin, 'r--')
    plt.plot(x, y_cos, 'b--')
    plt.plot(t, u, 'g-')
    plt.plot(t, v, 'y-')
    plt.xlabel('t')
    plt.ylabel('f(t)')
    plt.title('f(t)')
    plt.axis([0, 35, -2, 2])
    plt.show()