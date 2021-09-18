# -*- coding: utf-8 -*-
"""
Created on Tue June 30 13:21:00 2020

@author: Ola Nilsen
"""

import tkinter as tk
from time import time, localtime
from pathlib import Path
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.animation as animation
import math
from collections import deque
from autotune import PIDAutotune

#from PID import PID # PID-controller
"""Ivmech PID Controller is simple implementation of a Proportional-Integral-Derivative (PID) Controller in the Python Programming Language.
More information about PID Controller: http://en.wikipedia.org/wiki/PID_controller
"""



class PID:
    """PID Controller
    """

    def __init__(self, P=0.2, I=0.0, D=0.0, current_time=None):

        self.Kp = P
        self.Ki = I
        self.Kd = D

        self.sample_time = 0.00
#        self.current_time = current_time if current_time is not None else time.time()
        self.current_time = current_time if current_time is not None else time()
        self.last_time = self.current_time

        self.clear()

    def clear(self):
        """Clears PID computations and coefficients"""
        self.SetPoint = 0.0

        self.PTerm = 0.0
        self.ITerm = 0.0
        self.DTerm = 0.0
        self.last_error = 0.0

        # Windup Guard
        self.int_error = 0.0
        self.windup_guard = 20.0

        self.output = 0.0

    def update(self, feedback_value, current_time=None):
        """Calculates PID value for given reference feedback
        .. math::
            u(t) = K_p e(t) + K_i \int_{0}^{t} e(t)dt + K_d {de}/{dt}
        .. figure:: images/pid_1.png
           :align:   center
           Test PID with Kp=1.2, Ki=1, Kd=0.001 (test_pid.py)
        """
        error = self.SetPoint - feedback_value

#        self.current_time = current_time if current_time is not None else time.time()
        self.current_time = current_time if current_time is not None else time()
        delta_time = self.current_time - self.last_time
        delta_error = error - self.last_error

        if (delta_time >= self.sample_time):
            self.PTerm = self.Kp * error
            self.ITerm += error * delta_time

            if (self.ITerm < -self.windup_guard):
                self.ITerm = -self.windup_guard
            elif (self.ITerm > self.windup_guard):
                self.ITerm = self.windup_guard

            self.DTerm = 0.0
            if delta_time > 0:
                self.DTerm = delta_error / delta_time

            # Remember last time and last error for next calculation
            self.last_time = self.current_time
            self.last_error = error

            self.output = self.PTerm + (self.Ki * self.ITerm) + (self.Kd * self.DTerm)

    def setKp(self, proportional_gain):
        """Determines how aggressively the PID reacts to the current error with setting Proportional Gain"""
        self.Kp = proportional_gain

    def setKi(self, integral_gain):
        """Determines how aggressively the PID reacts to the current error with setting Integral Gain"""
        self.Ki = integral_gain

    def setKd(self, derivative_gain):
        """Determines how aggressively the PID reacts to the current error with setting Derivative Gain"""
        self.Kd = derivative_gain

    def setWindup(self, windup):
        """Integral windup, also known as integrator windup or reset windup,
        refers to the situation in a PID feedback controller where
        a large change in setpoint occurs (say a positive change)
        and the integral terms accumulates a significant error
        during the rise (windup), thus overshooting and continuing
        to increase as this accumulated error is unwound
        (offset by errors in the other direction).
        The specific problem is the excess overshooting.
        """
        self.windup_guard = windup

    def setSampleTime(self, sample_time):
        """PID that should be updated at a regular interval.
        Based on a pre-determined sampe time, the PID decides if it should compute or return immediately.
        """
        self.sample_time = sample_time
# Imports for connections
#import Adafruit_GPIO # Part of thermocouple configuration
#from Adafruit_MAX31856 import MAX31856 as MAX31856 # Library for thermocouple-chip.
#from gpiozero import DigitalOutputDevice as DOD # Control for relay

#from picamera import PiCamera

#import for pwm
#import RPi.GPIO as GPIO

#https://stackoverflow.com/questions/2398661/schedule-a-repeating-event-in-python-3
from threading import Timer

class RepeatedTimer(object):
    def __init__(self, interval, function, *args, **kwargs):
        self._timer     = None
        self.function   = function
        self.interval   = interval
        self.args       = args
        self.kwargs     = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self):
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False

#End: https://stackoverflow.com/questions/2398661/schedule-a-repeating-event-in-python-3

class Command:
    Next = 0
    Previous = 0
    line = 0

    def __init__(self, ID, line):
        self.ID = ID
        self.line = line

class com_section(Command):

    counter = 0
    temp_rate = 0
    temp_b = 0

    def __init__(self, ID, line, section_no,Start_temp,End_temp,Time,Capture_rate,Start_time):
        self.ID = int(ID)
        self.line = int(line)
        self.section_no = int(section_no)
        self.Start_temp = float(Start_temp)
        self.End_temp = float(End_temp)
        self.Time = int(Time)*60 #Converting minutes to sectonds for somplicity
        self.Capture_rate = int(Capture_rate)
        self.Start_time = Start_time

        self.temp_rate = (self.End_temp-self.Start_temp)/self.Time
        if self.temp_rate*self.Start_time != 0:
#            self.temp_b = self.Start_temp/(self.temp_rate*self.Start_time)
            self.temp_b = self.Start_temp-(self.temp_rate*self.Start_time)
        else:
            self.temp_b = self.Start_temp

    def hello(self):
        return str("Appended command Section in line: "+str(self.line)+" with Time (s): "+str(self.Time)+"\n")

    def plot_val(self):
        return [self.Start_time,self.Start_temp,float(self.Start_time+self.Time),self.End_temp]

    def set_temp(self,now_time):
        return float(self.temp_rate*now_time+self.temp_b)

class com_end(Command):

    def __init__(self, ID, line, Time):
        self.ID = ID
        self.line = line
        self.Time = Time

    def hello(self):
        return str("Appended command End in line: "+str(self.line)+"\n")


class sens:
    temperature = 0
    init_time = 0

    def __init__(self):
        self.temperature = 20
        self.init_time = time()

    def read_temp_c(self):
        self.temperature = self.oven.temperature()
        return self.temperature

    def setOven(self,oven):
        self.oven = oven

class PulseWM:
    pwm = 0

    def __init__(self):
        self.pwm = 0

    def ChangeDutyCycle(self,val):
        self.pwm = val

    def start(self,val):
        self.pwm = val

class Oven:
    delay = 3 #Delay tid i respons på ovnen
    TimeConst = 60 #Tidskonstant for hvor hurtig ovnen når 63.2% av endringen den er satt til å gjøre
    gain = 0.25 #hvor mye målt verdi endres for PV (er det 100% oppnåelig??)
    oneByTi = 1

    PV_history = []
    OP_history = []

    num_history = 3
    start_time = 0
    frequency = 1

    def __init__(self,D,T,G):
        self.delay = D
        self.TimeConst = T
        self.gain = G
        self.oneByTi = math.exp(-1/self.TimeConst)

    def start(self,start_time,frequency):
        self.start_time = start_time
        self.frequency = frequency
        self.num_history = (self.delay+1)*frequency

        for x in range(self.num_history): #gjør klar ovnen med ingen historie
            self.PV_history.append(20) # setter romtemperatur for alle cellene bakover i tid
            self.OP_history.append(0)

    def update(self,power,time):
        for x in range(len(self.PV_history)-1):
            self.PV_history[x] = self.PV_history[x+1]
            self.OP_history[x] = self.OP_history[x+1]

        self.OP_history[-1] = power
        self.PV_history[-1] = self.OP_history[0]*self.gain*(1-self.oneByTi)+self.PV_history[-2]*self.oneByTi

    def temperature(self):
        return self.PV_history[-1]

class Window(tk.Frame):

    status = 0    # 0 = None loaded, 1 = loaded not started, 2 = started, 3 = pause, 4 = quit, 5 = ended fine

    tuning = False

    commands = []
    output_folder = ""
    run_time = 0

    tune_run_time = 0
    tune_start_time = 0
    tune_status = 0

    start_run_time = 0
    pause_time = 0
    max_time = 0
    set_point = 0

    setT = 20

    #For plotting
    x_profile = [] #punktene som skal inn itemperaturprofil
    y_profile = []

    x_readtemp = []
    y_readtemp = []

    #Default settings for MAX31856 board - to be updated from config file
    SPI_PORT = 0
    SPI_DEVICE = 0

    #Default settings for Relay - to be updated from config file
    RelayPin = 1

    #Default PID settings - to be updated from config file
    P = 0.3371666944639706 #0.616 # 1      0.616
    I = 0.04249983878222571 #0.084 # 0.05
    D = 0.19106186738105319 #1.1218 # 0.05
    Capture_rate = 30

#    camera = PiCamera()


    # Define settings upon initialization. Here you can specify
    def __init__(self, master=None):

        # parameters that you want to send through the Frame class.
        tk.Frame.__init__(self, master)
        #reference to the master widget, which is the tk window
        self.master = master
        #with that, we want to then run init_window, which doesn't yet exist
        self.init_window()

    #Creation of init_window
    def init_window(self):

        # changing the title of our master widget
        self.master.title("PID Tuner v0.8 by Ola Nilsen")
        # allowing the widget to take the full space of the root window
        #self.pack(fill=tk.BOTH, expand=1)
        self.place(relheight = 1, relwidth = 1)

        self.buttonLoad     = tk.Button(self,text="Load profile", command=self.load_runfile)
        self.buttonRun      = tk.Button(self,text="Start", command=self.run_program)
        self.buttonPause    = tk.Button(self,text="Pause", command=self.pause_sequence)
        self.buttonStop     = tk.Button(self,text="Stop", command=self.stop_sequence)
        self.buttonTune     = tk.Button(self,text="Tune",command=self.begin_tuning)
        self.buttonPlot     = tk.Button(self,text="Edit PID",command=self.edit_PID)

        height = 35
        sep = 2
        width = 100

        self.buttonLoad.place(height = height, width = width, x = sep, y = sep)
        self.buttonRun.place(height = height, width = width, x = sep*2+width, y = sep)
        self.buttonPause.place(height = height, width = width, x = sep, y = (height*1+sep*2))
        self.buttonStop.place(height = height, width = width, x = sep*2+width, y = (height*1+sep*2))
        self.buttonTune.place(height = height, width = width, x = sep*3+width*2, y = sep)
        self.buttonPlot.place(height = height, width = width, x = sep*3+width*2, y = (height*1+sep*2))

        self.chkSimVar = tk.IntVar()
        self.chkSimulation  = tk.Checkbutton(self, text = "Simulating oven", variable = self.chkSimVar)
        self.chkSimulation.place(x = sep*5+width*6, y = (height*1+sep*2))
        self.chkSimulation.select()

        self.canvasPlot = tk.Canvas(self,width = 800, height=400)
        self.canvasPlot.place(x = 1, y = 75)

        labelheigth = 16 #24
        labelshift = 80

        self.labelTimeT = tk.Label(self, text = "Time (sec)     :")
        self.labelTimeT.place(x = width*3+sep*4, y = sep)

        self.labelSetTempT = tk.Label(self, text = "Set temp (C)  :")
        self.labelSetTempT.place(x = width*3+sep*4, y = sep+labelheigth)

        self.labelMesTempT = tk.Label(self, text = "Mes temp (C):")
        self.labelMesTempT.place(x = width*3+sep*4, y = sep*2+labelheigth*2)

        self.labelPowerT = tk.Label(self, text = "Power (%)    :")
        self.labelPowerT.place(x = width*3+sep*4, y = sep*2+labelheigth*3)

#

        self.labelTimeText = tk.StringVar()
        self.labelTime = tk.Label(self, textvariable = self.labelTimeText)
        self.labelTime.place(x = labelshift+width*3+sep*4, y = sep)
        self.labelTimeText.set(0)

        self.labelSetTempText = tk.StringVar()
        self.labelSetTemp = tk.Label(self, textvariable = self.labelSetTempText)
        self.labelSetTemp.place(x = labelshift+width*3+sep*4, y = sep+labelheigth)
        self.labelSetTempText.set(30)

        self.labelMesTempText = tk.StringVar()
        self.labelMesTemp = tk.Label(self, textvariable = self.labelMesTempText)
        self.labelMesTemp.place(x = labelshift+width*3+sep*4, y = sep*2+labelheigth*2)
        self.labelMesTempText.set(89)

        self.labelFileNameText = tk.StringVar()
        self.labelFileName = tk.Label(self, textvariable = self.labelFileNameText)
        self.labelFileName.place(x = labelshift+width*4+sep*5, y = sep*1+labelheigth*1)
        self.labelFileNameText.set("")

        self.labelPowerText = tk.StringVar()
        self.labelPower = tk.Label(self, textvariable = self.labelPowerText)
        self.labelPower.place(x = labelshift+width*3+sep*4, y = sep*2+labelheigth*3)
        self.labelPowerText.set(0)

        self.labelEnterSp = tk.Label(self, text = "Setpoint:")
        self.labelEnterSp.place(x = width*5+sep*5, y = sep*1+labelheigth*0)
        self.enterSp = tk.Entry(self)
        self.enterSp.place(x = width*5+sep*5+55, y = sep*1+labelheigth*0)
        self.enterSp.insert(2,"30")

        self.set_point = float(self.enterSp.get())

        self.direct_load_config_file() # Oppstartsprosedyre nr. 2

        self.x_readtemp = deque(maxlen = 600)
        self.y_readtemp = deque(maxlen = 600)
        self.x_profile = deque(maxlen = 600)
        self.y_profile = deque(maxlen = 600)

        # Plotteting
        plt.ion()
        self.fig = plt.Figure() #Plot av lest temperaturprofil // Byttet til stor Figure
        self.ax = self.fig.add_subplot(111)
        #plt.style.use('fivethirtyeight')
        self.line, = self.ax.plot(self.x_readtemp, self.y_readtemp, color='blue',linewidth=1)
        self.line2, = self.ax.plot(self.x_profile, self.y_profile, color='red',linewidth=1)

        self.tempProfile = FigureCanvasTkAgg(self.fig, self.canvasPlot)
        self.tempProfile.get_tk_widget().place(height = 400, width = 800, x = 1, y = 1)

        self.ani = animation.FuncAnimation(self.fig, self.animate, interval = 1000, blit=False)

        self.start_run_time = time()

    def donothing(self):
        return

    def animate(self,i): #Her kommer det to tilstander: Ikke startet, startet.

        #les temperatur
        temp = self.sensor.read_temp_c()
        self.labelMesTempText.set(str('%.2f' % temp)) #Setter variablen på lest temp


        if(self.status==0): # none loaded
            tid = time()-self.start_run_time
            self.x_readtemp.append(tid)
            self.y_readtemp.append(temp)
            self.x_profile.append(tid)
#            self.y_profile.append(float(self.enterSp.get()))
            self.y_profile.append(float(self.labelSetTempText.get()))
            self.donothing()
        elif(self.status==1): #loaded not started
            tid = time()-self.start_run_time
            self.x_readtemp.append(tid)
            self.y_readtemp.append(temp)
            self.x_profile.append(tid)
#            self.y_profile.append(float(self.enterSp.get()))
            self.y_profile.append(float(self.labelSetTempText.get()))
            self.donothing()
        elif(self.status==2): # started
            tid = time()-self.start_run_time
            self.x_readtemp.append(tid)
            self.y_readtemp.append(temp)
            self.x_profile.append(tid)
#            self.y_profile.append(float(self.enterSp.get()))
            self.y_profile.append(float(self.labelSetTempText.get()))
        elif(self.status==3): # pause/halt
            tid = self.run_time
            self.x_readtemp.append(tid)
            self.y_readtemp.append(temp)
            self.x_profile.append(tid)
#            self.y_profile.append(float(self.enterSp.get()))
            self.y_profile.append(float(self.labelSetTempText.get()))
        elif(self.status==4): #quit
            self.donothing()
        elif(self.status==5): #ended fine
            self.donothing()

        self.line.set_xdata(self.x_readtemp)  # update the data
        self.line.set_ydata(self.y_readtemp)  # update the data
        self.line2.set_xdata(self.x_profile)  # update the data
        self.line2.set_ydata(self.y_profile)  # update the data
        self.ax.relim()
        self.ax.autoscale_view()

        return self.line,


    def client_exit(self):
        exit()

    def begin_tuning(self):
        print("Begin tuning")
        self.buttonTune.config(relief = "sunken")
        self.tuning = True
        self.tune_status = 2 #begin tuning
        self.tune_start_time = time()

        maxlen = max(1, round(15 / 1)) #maxlen = max(1, round(delay / args.sampletime))
        self.delayed_temps = deque(maxlen=maxlen)
        self.delayed_temps.extend(maxlen * [float(self.enterSp.get())])



        self.PID_AutoTune = PIDAutotune(self.set_temp(time()),10,1,60) #Rydde opp her!!
#    def __init__(self, setpoint, out_step=10, sampletime=5, lookback=60,
        self.PID_AutoTune.set_initial_output(float(self.labelPowerText.get()))

    def tune_a_time(self):

        self.tune_run_time = time() - self.tune_start_time

        if not self.PID_AutoTune.run(self.sensor.read_temp_c()):

            self.PID_pwm_set(self.PID_AutoTune.get_output())

            self.delayed_temps.append(self.sensor.read_temp_c())
        else:
            self.tuning = False
            #print og rapporter nye PID parametere!
            print("Tuning er ferdig!")
            self.tune_status = 10 # Done tuning

            if self.PID_AutoTune.state == PIDAutotune.STATE_SUCCEEDED:
                for rule in self.PID_AutoTune.tuning_rules:
                    params = self.PID_AutoTune.get_pid_parameters(rule)
                    print('rule: {0}'.format(rule))
                    print('Kp: {0}'.format(params.Kp))
                    print('Ki: {0}'.format(params.Ki))
                    print('Kd: {0}'.format(params.Kd))
                    print()

            #Skriv nye PID til fil
            start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]
            timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
                start_time[1], str(start_time[0])[2:],
                start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.

            tune_out = open("PID_Param_"+self.output_folder+timestamp+".txt","w")

            if self.PID_AutoTune.state == PIDAutotune.STATE_SUCCEEDED:
                for rule in self.PID_AutoTune.tuning_rules:
                    params = self.PID_AutoTune.get_pid_parameters(rule)
                    tune_out.write('rule: {0}\n'.format(rule))
                    tune_out.write('Kp: {0}\n'.format(params.Kp))
                    tune_out.write('Ki: {0}\n'.format(params.Ki))
                    tune_out.write('Kd: {0}\n'.format(params.Kd))
                    tune_out.write("\n")
            tune_out.close()

            self.buttonTune.config(relief = "raised")

            #self.stop_sequence()




    def set_output_folder(self):
        folder = filedialog.askdirectory() #path = Path('folder_of_camera_images')  # Path of folder containing jpg images
        self.output_folder = Path(folder)  # Path of folder containing jpg images
        print(str("Out-folder set to : "+str(self.output_folder)+"\n"))

    def edit_PID(self):
        print("P = "+str(self.P)+"\n")
        print("I = "+str(self.I)+"\n")
        print("D = "+str(self.D)+"\n")
        self.donothing()


    def pause_sequence(self):
        if(self.status == 2): # Program started
            self.status =3 #put it on pause
            self.buttonPause.config(text ="Paused", relief = "sunken")
            self.pause_time = time()
            self.run_time = time()-self.start_run_time
        elif(self.status == 3): # Program already paused
            self.buttonPause.config(text ="Pause", relief = "raised")
            self.status =2 #put it back on track
            t = self.start_run_time + (time() - self.pause_time) #fikser en ny init-tid som tar hensyn til at den har vært pauset
            self.start_run_time = t

    def stop_sequence(self):
        try:
            self.rt.stop()
            self.f_out.close()
        except:
            print("No program was run")
        self.status=5


    def set_temp(self, time):

        if(self.status == 2): #started
            i = 0
            while time > self.commands[i].Time:
                i = i + 1

            temp = self.commands[i].set_temp(time)
#            self.enterSp.insert(2,temp)

        else:
            temp = float(self.enterSp.get()) # Her kommer fint sted å hente verdi fra variabel !!!

        return temp



    def run_sequence(self): #denne kjøres ved en gitt frekvens. Standard = 1 Hz
        t = 0
        if(self.status == 2): #program started
            self.setT = self.set_temp(time()-self.start_run_time)
            self.pid.SetPoint = self.setT
            self.labelSetTempText.set(str('%.2f' % self.setT))
            self.labelTimeText.set(str('%.i' % (time()-self.start_run_time))) #Oppdaterer label på førstesiden
            t = time()-self.start_run_time
            if (self.tuning):
                self.tune_a_time()
        elif(self.status == 3): #Paused
            t = self.run_time
            self.donothing()
        elif(self.status == 0): #nothing loaded, men... da er den jo også ikke startet...
            t = time()-self.start_run_time
            self.setT = float(self.enterSp.get())
            self.pid.SetPoint = self.setT
            self.labelSetTempText.set(str('%.2f' % float(self.setT)))
            if (self.tuning):
                self.tune_a_time()
        elif(self.status == 1): #loaded, but not started
            t = time()-self.start_run_time
            self.setT = float(self.enterSp.get())
            self.pid.SetPoint = self.setT
            self.labelSetTempText.set(str('%.2f' % float(self.setT)))
            if (self.tuning):
                self.tune_a_time()



        #Les-temp:
        temp = self.sensor.read_temp_c()
        if not self.tuning:
            self.pid.update(temp) # Calculates PID-value
            self.PID_pwm_reg(temp)
        #Plot i graf = skjer i animate
        #Plot i datafil
        self.f_out.write(str('%.2f' % t)+"\t"+str('%.2f' % self.setT)+"\t"+str('%.2f' % temp)+"\tNoCapt\n")
        #feil!!! fila var lukket

    def rate(self,val,limit):
        if(val>limit):
            return 100
        elif(val<0):
            return 0
        else:
            return (val/limit)*100


    def PID_pwm_set(self,val):
        self.oven.update(val,time())
        self.labelPowerText.set(str('%.1f' % val))
        #val = 100*max(min(pid.output, 1),0) # Turns PID-value into something between 0 and 10
        #print(("Val: {:.3f}").format(val))
        self.p.ChangeDutyCycle(val)


    def PID_pwm_reg(self,temp):
        self.pid.update(temp) # Calculates PID-value

        val = self.rate(self.pid.output,10)
        self.oven.update(val,time())
        self.labelPowerText.set(str('%.1f' % val))
        #val = 100*max(min(pid.output, 1),0) # Turns PID-value into something between 0 and 10
        #print(("Val: {:.3f}").format(val))
        self.p.ChangeDutyCycle(val)

#        self.pid.SetPoint = float(self.enterSp.get())
#        self.pid.update(temp) # Calculates PID-value


    def run_program(self):

        self.status=2 #Run program
        self.start_run_time = time() #computer start-time (since epoc)
#        start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]
 #       timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
 #           start_time[1], str(start_time[0])[2:],
 #           start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.

#        self.f_out = open(self.output_folder+timestamp+".txt","w")
#        self.f_out.write("Time\tSetTemp\tRealTemp\tCaptureNo\n")

        self.x_readtemp.clear()
        self.y_readtemp.clear()
        self.x_profile.clear()
        self.y_profile.clear()


#        self.rt = RepeatedTimer(1,self.run_sequence) #Controls temp and log
#       Legges nå til å starte ved oppstart

    def load_runfile(self):
        ID = 0
        line = 1
        in_file = tk.filedialog.askopenfilename(title = "Select runfile",filetypes = (("TXT","*.txt"),("all files","*.*")))
        print(in_file)
        print(Path(in_file).stem) # Filnavnet uten suffix
        self.labelFileNameText.set(Path(in_file).stem)
        print(Path(in_file).parent) # filstien uten filnavnet

        self.buttonLoad.config(state ="disabled")
        self.buttonRun.config(state ="normal")

        self.path = Path(in_file)
        with open(in_file,'r') as f:
            read_data = f.readline()
            while read_data != '':
                data = read_data.strip().split(" ")
                if data[0]=='': data[0]=" "

                #Setting default values for a section in case it is lacking from the file
                Start_temp = 20
                End_temp = 20
                Time = 10
                Capture_rate = 30

                if data[0].lower() == "section":
                    if data[1] is not None:
                        if int(data[1]) > 0:
                            section_no = data[1]
                            v = f.readline().strip().split("=")
                            line = line + 1
                            while len(v) == 2:
                                if str(v[0]).lower().strip() == "start_temp":
                                    Start_temp = int(v[1])
                                    print(str("Start_temp set to "+str(Start_temp)+"\n"))
                                elif str(v[0]).lower().strip() == "end_temp":
                                    End_temp = int(v[1])
                                    print(str("End_temp set to "+str(End_temp)+"\n"))
                                elif str(v[0]).lower().strip() == "time":
                                    Time = int(v[1])
                                    print(str("Time set to "+str(Time)+"\n"))
                                elif str(v[0]).lower().strip() == "capture_rate":
                                    Capture_rate = int(v[1])
                                    print(str("Capture_rate set to "+str(Capture_rate)+"\n"))
                                else:
                                    print(str("Syntax error in line: "+str(line)+"\n"))
                                    print(str("Expected something in a section, and found something else.\n"))
                                v = f.readline().strip().split("=")

                            if str(v[0]).lower().strip() == "end_section":
                                self.commands.append(com_section(ID,line,section_no,Start_temp,End_temp,Time,Capture_rate,self.max_time))
                                print(self.commands[ID].hello())
                                self.max_time = self.max_time+Time*60
                            else:
                                print(str("Syntax error in line: "+str(line)+"\n"))
                                print(str("Expected 'End_section' while found: "+str(v[0])+".\n"))

                            ID = ID + 1
                        else:
                            print(str("Syntax error in line: "+str(line)+"\n"))
                            print(str("Found zero or negative value.\n"))
                    else:
                        print(str("Syntax error in line: "+str(line)+"\n"))
                        print(str("Lack value after command: Section.\n"))

                elif data[0].lower() == "end":
                    self.commands.append(com_end(ID,line,self.max_time))
                    print(self.commands[ID].hello())
                    ID = ID + 1


                elif data[0][0].lower() == "#":
                    print("Found #!\n")
                elif data[0].lower() == " ":
                    print("Found space!\n")
                else:
                    print(str("Syntax error in line: "+str(line)+"\n"))
                    print(str("Found \""+str(data[0])+"\" and did not understand it.\n"))

                line = line + 1
                read_data = f.readline()

        f.close()
        print("Precursortester program file was loaded.\n")
        print("Total runtime = "+str('%.3f' % self.max_time)+" s.\n")
        if(len(self.commands)>0): #Måtte legge inn en slik sjekk så den ikke slet med tom liste ved oppstart.
            self.labelSetTempText.set(self.commands[0].set_temp(0)) #Settter variablen på starttemp på programmet. tid = 0
        self.status = 1 # Loaded, not started



    def direct_load_config_file(self):
        with open("Autotune_config.txt",'r') as f:
            data = f.readline().strip().split(" ")
            while data[0] != '':
                if data[0].lower() == "[max31856]":
                    print("[MAX31856]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0].lower()).strip() == "spi_port":
                            SPI_PORT = int(v[1])
                            print(str("SPI_PORT set to "+str(SPI_PORT)+"\n"))
                        elif str(v[0]).lower().strip() == "spi_device":
                            SPI_DEVICE = int(v[1])
                            print(str("SPI_DEVICE set to "+str(SPI_DEVICE)+"\n"))
                        else:
                            print(str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                if data[0].lower() == "[relay]":
                    print("[Relay]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0]).lower().strip() == "relaypin":
                            RelayPin = int(v[1])
                            print(str("RelayPin set to "+str(RelayPin)+"\n"))
                        else:
                            print(str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                if data[0].lower() == "[pid]":
                    print("[PID]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0]).lower().strip() == "p":
                            self.P = float(v[1])
                            print(str("P set to "+str(self.P)+"\n"))
                        elif str(v[0]).lower().strip() == "i":
                            self.I = float(v[1])
                            print(str("I set to "+str(self.I)+"\n"))
                        elif str(v[0]).lower().strip() == "d":
                            self.D = float(v[1])
                            print(str("D set to "+str(self.D)+"\n"))
                        else:
                            print(str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                data = f.readline().strip().split(" ")
        f.close()

#        Raspberry Pi hardware SPI configuration.
        self.SPI_PORT   = SPI_PORT
        self.SPI_DEVICE = SPI_DEVICE
        self.RelayPin = RelayPin
        self.sensor = MAX31856(hardware_spi=Adafruit_GPIO.SPI.SpiDev(self.SPI_PORT, self.SPI_DEVICE), tc_type=MAX31856.MAX31856_K_TYPE) #remove tc_type for non-K thermocouples
#        self.sensor = sens() #Dette er min dummy for å få programmet til å kjøre
#        self.oven = Oven(3,100,1) #Oven(3,10,0.25)
#        self.oven = Oven(3,100,1) #Oven(3,10,0.25)
#        self.oven.start(time(),1)
#        self.sensor.setOven(self.oven)
        # Relay configuration.
        # relay = DOD(0) # Creates a general DigitalOutputDevice using GPIO pin 0

        # PWM configuration
        GPIO.setmode(GPIO.BCM)
        GPIO.setup(RelayPin, GPIO.OUT)
        self.p = GPIO.PWM(RelayPin, 2) #BCM pin 1 set to 0.5 Hz frequency PWM
        self.p.start(0) #Starter PWM med 0 av 100 i duty cycle
##        self.p = PulseWM()


        self.pid = PID(self.P, self.I, self.D)
        self.pid.SetPoint = 20
        self.pid.setSampleTime(1)

        print("Autotune config file was loaded.\n")
        self.status = 0
        self.rt = RepeatedTimer(1,self.run_sequence) #Controls temp and log

        start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]
        timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
            start_time[1], str(start_time[0])[2:],
            start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.

        self.f_out = open(self.output_folder+timestamp+".txt","w")
        self.f_out.write("Time\tSetTemp\tRealTemp\tCaptureNo\n")





# root window created. Here, that would be the only window, but
# you can later have windows within windows.
root = tk.Tk()
root.geometry("800x480")

#creation of an instance
app = Window(root)

#mainloop
root.mainloop()
