# Code made for Uy Que Tusa
# 11 March 2021
# License MIT
# Introduction to Control: Python Program HAVOC Assigment

import control as co
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.integrate import solve_ivp

sns.set()

T = np.genfromtxt("Data1", delimiter=';')
V = np.genfromtxt("Data2", delimiter=';')

fig, ax = plt.subplots(2)
ax[0].plot(T[0], T[1], 'r', alpha=0.75)
ax[1].plot(V[0], V[1], 'b', alpha=0.75)
ax[0].set_ylim(20, 30)
ax[1].set_xlabel("Time \t$[s]$ ", fontsize=14)
ax[0].set_ylabel("Temperature \t$[^{\circ}C]$ ", fontsize=14)
ax[1].set_ylabel("Voltage\t$[V]$", fontsize=14)
fig.suptitle('HVAC CONTROL\nTemperature $[^{\circ}C]$ & Voltage $[W]$', fontsize=16)

# Achieve the position of the Steady State
T_postn_1 = np.where(np.isclose(T[1], 24, 1e-5))[0][0]
T_postn_2 = np.where(np.isclose(T[1], 26.1, 1e-5))[0][0]

# Values from the graphic for Temperature and Voltage
T_1_ss = T[1][T_postn_1]
V_1_ss = V[1][T_postn_1]
T_2_ss = T[1][T_postn_2]
V_2_ss = V[1][T_postn_2]

# Values for Density_air,Cp_air,Area_tube,Velocity_air
rho = 1.29  # [Kg/m^3]
Cp = 1.005  # [Kj/Kg*K]
Long = 0.14  # [m]
Diam = 0.045  # [m]
Area = (np.pi * (Diam ** 2)) / 4  # [m^2]
Vel = 1.1  # [m/s]
I_curr = ((T_2_ss - T_1_ss) / (V_2_ss - V_1_ss)) * rho * Cp * Area * Vel

# Phenomenological Parameters
Tau = Long / Vel  # [1/s]
K_u = I_curr / (rho * Cp * Area * Vel)


# Define the ED of the phenomelogical system
def HVAC(t, Y, K_u, Tau, V_2_ss, T_1_ss, Tau_s):
    # Parameters
    Ku = K_u
    tau = Tau
    u1 = V_2_ss
    T_amb = T_1_ss
    K_s = 1

    # ED in the form dy/dt = y +...
    dy_dt = (1 / tau) * (- Y + Ku * u1 + T_amb) + (1 / Tau_s) * K_s
    return dy_dt


# Solve the ED
sol_u = solve_ivp(HVAC, [0, 7 * Tau], [T_1_ss], args=(K_u, Tau, V_2_ss, T_1_ss, 4),
                  t_eval=np.arange(0, 7 * Tau, 0.0001))

# plotting the solve of the ED
plt.figure(2)
plt.plot(sol_u.t, sol_u.y[0])
plt.xlabel('Time $[s]$', fontsize=14)
plt.ylabel('Temperature $[^{\circ}C]$', fontsize=14)
plt.title('Solution E.D of HAVOC System', fontsize=16)

# Define the system in Transfer Functions and Analyze the Open Loop
G_p = co.tf([K_u], [Tau, 1])
G_c = co.tf([1, 0.2], [1, 0])
G_a = co.tf([4], [1])
G_s = co.tf([1], [4, 1])
G_d = co.tf([1], [Tau, 1])

# Open loop
G_u_Sys_ol = co.series(G_p, G_s)
G_d_Sys_ol = co.series(G_d, G_s)

# Here define the time of Analysis
t_i = 20
t_ch = 50
t_f = 80
step_op = 1

t = np.linspace(t_i, t_ch, 1000)
t_1, y1 = co.forced_response(G_u_Sys_ol, t, step_op)
t2 = np.linspace(t_ch, t_f, 1000)
t_2, y2 = co.forced_response(G_d_Sys_ol, t2, step_op)

# close loop
G_u_Sys_cl = co.series(G_c, G_a, G_p, G_s)
G_d_Sys_ol = co.series(G_d, G_s)

G_fbck_u = co.feedback(G_u_Sys_cl)
G_fbck_d = co.feedback(G_d_Sys_ol, G_c)

# Here define the time of Analysis
step_cl_u = 3
step_cl_d = 1
t = np.linspace(t_i, t_ch, 1000)
t_3, y3 = co.forced_response(G_fbck_u, t, step_cl_u)
t2 = np.linspace(t_ch, t_f, 1000)
t_4, y4 = co.forced_response(G_fbck_d, t2, step_cl_d)


# Class of the Figure 4
class figure_4:
    fig4, ax = plt.subplots(3, 2, constrained_layout=True, gridspec_kw={'height_ratios': [2, 1, 1]})
    fig4.suptitle('HVAC CONTROL\nTemperature $[^{\circ}C]$ & Voltage $[W]$', fontsize=16)
    ax[0, 0].plot(t_1, T_1_ss + y1, label='Step Response from u open loop', alpha=0.75)
    ax[0, 0].plot(t_2, T_1_ss + y2 + y1[999], label='Step Response from d open loop')
    ax[0, 0].set_xlim(t_i - 10, t_f)
    ax[0, 1].plot(t_3, T_1_ss + y3, label='Step Response from u closed loop')
    ax[0, 1].plot(t_4, T_1_ss + y4 + y3[999], label='Step Response from d closed loop')
    ax[0, 1].set_xlim(t_i - 10, t_f)
    ax[1, 0].hlines(y=0, xmin=0, xmax=t_i, color='green')
    ax[1, 0].hlines(y=step_op, xmin=t_i, xmax=t_f, color='green')
    ax[1, 0].vlines(x=t_i, ymin=0, ymax=step_op, color='green')
    ax[1, 0].set_xlim(t_i - 10, t_f)
    ax[1, 1].hlines(y=0, xmin=0, xmax=t_i, color='green')
    ax[1, 1].hlines(y=step_cl_u, xmin=t_i, xmax=t_f, color='green')
    ax[1, 1].vlines(x=t_i, ymin=0, ymax=step_cl_u, color='green')
    ax[1, 1].set_xlim(t_i - 10, t_f)
    ax[2, 0].hlines(y=T_1_ss, xmin=0, xmax=t_ch, color='red')
    ax[2, 0].hlines(y=T_1_ss + step_op, xmin=t_ch, xmax=t_f, color='red')
    ax[2, 0].vlines(x=t_ch, ymin=T_1_ss, ymax=T_1_ss + step_op, color='red')
    ax[2, 0].set_xlim(t_i - 10, t_f)
    ax[2, 1].hlines(y=T_1_ss, xmin=0, xmax=t_ch, color='red')
    ax[2, 1].hlines(y=T_1_ss + step_cl_d, xmin=t_ch, xmax=t_f, color='red')
    ax[2, 1].vlines(x=t_ch, ymin=T_1_ss, ymax=T_1_ss + step_cl_d, color='red')
    ax[2, 1].set_xlim(t_i - 10, t_f)
    ax[0, 0].set_ylim(T_1_ss, T_1_ss + step_cl_u + 1)
    ax[0, 1].set_ylim(T_1_ss, T_1_ss + step_cl_u + 1)
    ax[2, 0].set_xlabel("Time \t$[s]$ ", fontsize=14)
    ax[2, 1].set_xlabel("Time \t$[s]$ ", fontsize=14)
    ax[0, 0].set_ylabel("Temperature \t$[^{\circ}C]$ ", fontsize=14)
    ax[1, 0].set_ylabel("Voltage\t$[V]$", fontsize=14)
    ax[2, 0].set_ylabel("$Temperature_{surr}$ \t$[^{\circ}C]$", fontsize=14)
    ax[0, 0].set_title("Open Loop", fontsize=14)
    ax[0, 1].set_title("Closed Loop", fontsize=14)


figure4 = figure_4()

# Roots
plt.figure(5)
roots = co.rlocus(G_u_Sys_ol + G_d_Sys_ol)
plt.figure(6)
roots2 = co.rlocus(G_fbck_u + G_fbck_d)

# print(G_fe_bck)
print("Transfer Function Open Loop:\n", G_u_Sys_ol, '+', G_d_Sys_ol)
print("Transfer Function Closed Loop:\n", G_fbck_u, '+', G_fbck_d)

plt.show()
