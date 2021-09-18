# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 08:58:51 2019

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

    commands = []
    output_folder = ""
    run_time = 0
    start_run_time = 0
    pause_time = 0
    max_time = 0

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
    P = 1
    I = 0.05
    D = 0.05
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
        self.master.title("Precursortester v0.8 by Ola Nilsen")
        # allowing the widget to take the full space of the root window
        #self.pack(fill=tk.BOTH, expand=1)
        self.place(relheight = 1, relwidth = 1)

        # creating a menu instance
        menu = tk.Menu(self.master)
        self.master.config(menu=menu)

        # create the file object)
        file = tk.Menu(menu)
        # adds a command to the menu option, calling it exit, and the
        # command it runs on event is client_exit
        file.add_command(label="Load Temperature program", command=self.load_runfile)
        file.add_command(label="Load Config file", command=self.load_config_file)
        file.add_command(label="Choose Out-folder", command=self.set_output_folder)
        file.add_separator()
        # create the file object)
        file.add_command(label="Plot temperature profile", command=self.plot_temp_profile)
        file.add_separator()
        file.add_command(label="Run program", command=self.run_program)
        file.add_separator()
        file.add_command(label="Stop program", command=self.stop_sequence)
        # create the file object)
        file.add_separator()
        file.add_command(label="Exit", command=self.client_exit)
        #added "file" to our menu
        menu.add_cascade(label="File", menu=file)

        camera = tk.Menu(menu)
        camera.add_command(label="Camera preview on", command=self.client_exit)
        camera.add_command(label="Capture one image", command=self.client_exit)
        camera.add_command(label="Camera preview off", command=self.client_exit)
        camera.add_separator()
        menu.add_cascade(label="Camera Control", menu=camera)



        self.buttonLoad     = tk.Button(self,text="Load",command=self.load_runfile)
        self.buttonRun      = tk.Button(self,text="Run", state = "disabled", command=self.run_program)
        self.buttonPause    = tk.Button(self,text="Pause", state = "disabled",command=self.pause_sequence) # denne prosedyren må lages!
        self.buttonStop     = tk.Button(self,text="Stop", state = "disabled",command=self.stop_sequence)
        self.buttonExit     = tk.Button(self,text="Capture",command=self.client_exit)
        self.buttonPlot     = tk.Button(self,text="Focus",command=self.plot_temp_profile)

        height = 35
        sep = 2
        width = 100

        self.buttonLoad.place(height = height, width = width, x = sep, y = sep)
        self.buttonRun.place(height = height, width = width, x = sep*2+width, y = sep)
        self.buttonPause.place(height = height, width = width, x = sep, y = (height*1+sep*2))
        self.buttonStop.place(height = height, width = width, x = sep*2+width, y = (height*1+sep*2))
        self.buttonExit.place(height = height, width = width, x = sep*3+width*2, y = sep)
        self.buttonPlot.place(height = height, width = width, x = sep*3+width*2, y = (height*1+sep*2))
#        self.buttonPause.place(height = height, width = 100, x = 5, y = (480-height*2-sep*2))
#        self.buttonQuit.place(height = height, width = 100, x = 5, y = (480-height*1-sep*1))

#        self.buttonRectangle = tk.Button(self,text="Confirm rectangle",command=self.getRectangle)
#        self.buttonRectangle.pack(side=tk.TOP)
#        self.buttonRectangle.config(state="disabled")
        
#        self.progressbar = ttk.Progressbar(self.master,orient="horizontal",length=300,mode="determinate")
#        self.progressbar.pack(side=tk.BOTTOM)

        self.scroll = tk.Scrollbar(self)
        self.textbox = tk.Text(self,height=6,width=75) #,height=4,width=50)
#        self.scroll.pack(side=tk.RIGHT, fill=tk.Y)
#        self.textbox.pack(side=tk.LEFT, fill=tk.Y)

        self.textbox.place(height = (480-(height*1+sep*3)), width = 335, x = 440, y = (height*1+sep*2))
        self.scroll.place(height = (480-(height*1+sep*3)), width = 20, x = 775, y = (height*1+sep*2))
        self.scroll.config(command=self.textbox.yview)
        self.textbox.config(yscrollcommand=self.scroll.set)

        self.textbox.insert(tk.END,"Welcome to the Precursortester.\n")
        self.textbox.insert(tk.END,"Your first task is to load the config file.\n")


        self.canvasPlot = tk.Canvas(self,width = 400, height=400)
        self.canvasPlot.place(x = 1, y = 75)
        self.canvasImage = tk.Canvas(self,width = 400, height=400)
        self.canvasImage.place(x = 400, y = 75)

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
        self.labelSetTempText.set(73)

        self.labelMesTempText = tk.StringVar()
        self.labelMesTemp = tk.Label(self, textvariable = self.labelMesTempText)
        self.labelMesTemp.place(x = labelshift+width*3+sep*4, y = sep*2+labelheigth*2)
        self.labelMesTempText.set(89)

        self.labelFileNameText = tk.StringVar()
        self.labelFileName = tk.Label(self, textvariable = self.labelFileNameText)
        self.labelFileName.place(x = labelshift+width*4+sep*5, y = sep*1+labelheigth*0)
        self.labelFileNameText.set("")

        self.labelPowerText = tk.StringVar()
        self.labelPower = tk.Label(self, textvariable = self.labelPowerText)
        self.labelPower.place(x = labelshift+width*3+sep*4, y = sep*2+labelheigth*3)
        self.labelPowerText.set(0)


        self.direct_load_config_file() # Leser inn konfig-fila direkte
 #       self.load_config_file() # Leser inn konfig-fila direkte

        # Plotteting
        plt.ion()
        self.fig = plt.Figure() #Plot av lest temperaturprofil // Byttet til stor Figure
        self.ax = self.fig.add_subplot(111)
        #plt.style.use('fivethirtyeight')
        self.line, = self.ax.plot(self.x_readtemp, self.y_readtemp, color='blue',linewidth=1)        

        self.tempProfile = FigureCanvasTkAgg(self.fig, self.canvasPlot)
        self.tempProfile.get_tk_widget().place(height = 400, width = 400, x = 1, y = 1)

        self.ani = animation.FuncAnimation(self.fig,self.animate,interval = 1000,blit=False)
        self.ax.plot(self.x_profile,self.y_profile,color='red',linewidth=1)        

    def donothing(self):
        return

    def animate(self,i): #Her kommer det to tilstander: Ikke startet, startet.

        #les temperatur
        temp = self.sensor.read_temp_c()
        self.labelMesTempText.set(str('%.2f' % temp)) #Setter variablen på lest temp

        
        if(self.status==0): # none loaded
            self.donothing()
        elif(self.status==1): #loaded not started
            self.donothing()
        elif(self.status==2): # started
            tid = time()-self.start_run_time
#            self.labelTimeText.set(tid)
            self.x_readtemp.append(tid)
            self.y_readtemp.append(temp)
#            self.labelSetTempText.set(self.set_temp(tid))
        elif(self.status==3): # pause/halt
            tid = self.run_time
#            self.labelTimeText.set(tid)
            self.x_readtemp.append(tid)
            self.y_readtemp.append(temp)
        elif(self.status==4): #quit
            self.donothing()
        elif(self.status==5): #ended fine
            self.donothing()

        self.line.set_xdata(self.x_readtemp)  # update the data
        self.line.set_ydata(self.y_readtemp)  # update the data
        
        
        return self.line,
    

    def client_exit(self):
        exit()

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
                                    self.textbox.insert(tk.END,str("Start_temp set to "+str(Start_temp)+"\n"))
                                elif str(v[0]).lower().strip() == "end_temp":
                                    End_temp = int(v[1])
                                    self.textbox.insert(tk.END,str("End_temp set to "+str(End_temp)+"\n"))
                                elif str(v[0]).lower().strip() == "time":
                                    Time = int(v[1])
                                    self.textbox.insert(tk.END,str("Time set to "+str(Time)+"\n"))
                                elif str(v[0]).lower().strip() == "capture_rate":
                                    Capture_rate = int(v[1])
                                    self.textbox.insert(tk.END,str("Capture_rate set to "+str(Capture_rate)+"\n"))
                                else:
                                    self.textbox.insert(tk.END,str("Syntax error in line: "+str(line)+"\n"))
                                    self.textbox.insert(tk.END,str("Expected something in a section, and found something else.\n"))
                                v = f.readline().strip().split("=")

                            if str(v[0]).lower().strip() == "end_section":
                                self.commands.append(com_section(ID,line,section_no,Start_temp,End_temp,Time,Capture_rate,self.max_time))
                                self.textbox.insert(tk.END,self.commands[ID].hello())
                                self.max_time = self.max_time+Time*60
                            else:
                                self.textbox.insert(tk.END,str("Syntax error in line: "+str(line)+"\n"))
                                self.textbox.insert(tk.END,str("Expected 'End_section' while found: "+str(v[0])+".\n"))
                            
                            ID = ID + 1
                        else: 
                            self.textbox.insert(tk.END,str("Syntax error in line: "+str(line)+"\n"))
                            self.textbox.insert(tk.END,str("Found zero or negative value.\n"))
                    else:
                        self.textbox.insert(tk.END,str("Syntax error in line: "+str(line)+"\n"))
                        self.textbox.insert(tk.END,str("Lack value after command: Section.\n"))

                elif data[0].lower() == "end":
                    self.commands.append(com_end(ID,line,self.max_time))
                    self.textbox.insert(tk.END,self.commands[ID].hello())
                    ID = ID + 1


                elif data[0][0].lower() == "#":
                    self.textbox.insert(tk.END,"Found #!\n")
                elif data[0].lower() == " ":
                    self.textbox.insert(tk.END,"Found space!\n")
                else:
                    self.textbox.insert(tk.END,str("Syntax error in line: "+str(line)+"\n"))
                    self.textbox.insert(tk.END,str("Found \""+str(data[0])+"\" and did not understand it.\n"))

                line = line + 1
                read_data = f.readline()

        f.close()
        self.textbox.insert(tk.END,"Precursortester program file was loaded.\n")
        self.textbox.insert(tk.END,"Total runtime = "+str('%.3f' % self.max_time)+" s.\n")
        if(len(self.commands)>0): #Måtte legge inn en slik sjekk så den ikke slet med tom liste ved oppstart.
            self.labelSetTempText.set(self.commands[0].set_temp(0)) #Settter variablen på starttemp på programmet. tid = 0
        self.status = 1 # Loaded, not started
        self.plot_temp_profile()

    def direct_load_config_file(self):
        with open("Precursortester_config.txt",'r') as f:
            data = f.readline().strip().split(" ")
            while data[0] != '': 
                if data[0].lower() == "[max31856]":
                    self.textbox.insert(tk.END,"[MAX31856]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0].lower()).strip() == "spi_port":
                            SPI_PORT = int(v[1])
                            self.textbox.insert(tk.END,str("SPI_PORT set to "+str(SPI_PORT)+"\n"))
                        elif str(v[0]).lower().strip() == "spi_device":
                            SPI_DEVICE = int(v[1])
                            self.textbox.insert(tk.END,str("SPI_DEVICE set to "+str(SPI_DEVICE)+"\n"))
                        else:
                            self.textbox.insert(tk.END,str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                if data[0].lower() == "[relay]":
                    self.textbox.insert(tk.END,"[Relay]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0]).lower().strip() == "relaypin":
                            RelayPin = int(v[1])
                            self.textbox.insert(tk.END,str("RelayPin set to "+str(RelayPin)+"\n"))
                        else:
                            self.textbox.insert(tk.END,str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                if data[0].lower() == "[pid]":
                    self.textbox.insert(tk.END,"[PID]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0]).lower().strip() == "p":
                            self.P = float(v[1])
                            self.textbox.insert(tk.END,str("P set to "+str(self.P)+"\n"))
                        elif str(v[0]).lower().strip() == "i":
                            self.I = float(v[1])
                            self.textbox.insert(tk.END,str("I set to "+str(self.I)+"\n"))
                        elif str(v[0]).lower().strip() == "d":
                            self.D = float(v[1])
                            self.textbox.insert(tk.END,str("D set to "+str(self.D)+"\n"))
                        else:
                            self.textbox.insert(tk.END,str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                data = f.readline().strip().split(" ")
        f.close()
        
#        Raspberry Pi hardware SPI configuration.
##        self.SPI_PORT   = SPI_PORT
##        self.SPI_DEVICE = SPI_DEVICE
##        self.RelayPin = RelayPin 
##        self.sensor = MAX31856(hardware_spi=Adafruit_GPIO.SPI.SpiDev(self.SPI_PORT, self.SPI_DEVICE), tc_type=MAX31856.MAX31856_K_TYPE) #remove tc_type for non-K thermocouples
        self.sensor = sens() #Dette er min dummy for å få programmet til å kjøre
        self.oven = Oven(3,100,1) #Oven(3,10,0.25)
        self.oven.start(time(),1)
        self.sensor.setOven(self.oven)
        # Relay configuration.
        # relay = DOD(0) # Creates a general DigitalOutputDevice using GPIO pin 0
        
        # PWM configuration
##        GPIO.setmode(GPIO.BCM)
##        GPIO.setup(RelayPin, GPIO.OUT)
##        self.p = GPIO.PWM(RelayPin, 2) #BCM pin 1 set to 0.5 Hz frequency PWM
##        self.p.start(0) #Starter PWM med 0 av 100 i duty cycle
        self.p = PulseWM()


        self.pid = PID(self.P, self.I, self.D)
        self.pid.SetPoint = 20
        self.pid.setSampleTime(1)
                
        self.textbox.insert(tk.END,"Precursortester config file was loaded.\n")

    def load_config_file(self):
        in_file = tk.filedialog.askopenfilename(title = "Select config file",filetypes = (("TXT","*.txt"),("all files","*.*")))
        print(in_file)
        with open(in_file,'r') as f:
            data = f.readline().strip().split(" ")
            while data[0] != '': 
                if data[0].lower() == "[max31856]":
                    self.textbox.insert(tk.END,"[MAX31856]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0].lower()).strip() == "spi_port":
                            SPI_PORT = int(v[1])
                            self.textbox.insert(tk.END,str("SPI_PORT set to "+str(SPI_PORT)+"\n"))
                        elif str(v[0]).lower().strip() == "spi_device":
                            SPI_DEVICE = int(v[1])
                            self.textbox.insert(tk.END,str("SPI_DEVICE set to "+str(SPI_DEVICE)+"\n"))
                        else:
                            self.textbox.insert(tk.END,str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                if data[0].lower() == "[relay]":
                    self.textbox.insert(tk.END,"[Relay]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0]).lower().strip() == "relaypin":
                            RelayPin = int(v[1])
                            self.textbox.insert(tk.END,str("RelayPin set to "+str(RelayPin)+"\n"))
                        else:
                            self.textbox.insert(tk.END,str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                if data[0].lower() == "[pid]":
                    self.textbox.insert(tk.END,"[PID]\n")
                    v = f.readline().strip().split("=")
                    while len(v) == 2:
                        if str(v[0]).lower().strip() == "p":
                            self.P = float(v[1])
                            self.textbox.insert(tk.END,str("P set to "+str(self.P)+"\n"))
                        elif str(v[0]).lower().strip() == "i":
                            self.I = float(v[1])
                            self.textbox.insert(tk.END,str("I set to "+str(self.I)+"\n"))
                        elif str(v[0]).lower().strip() == "d":
                            self.D = float(v[1])
                            self.textbox.insert(tk.END,str("D set to "+str(self.D)+"\n"))
                        else:
                            self.textbox.insert(tk.END,str("Did not understand: "+str(v[0]).lower()+"\n"))
                        v = f.readline().strip().split("=")
                data = f.readline().strip().split(" ")
        f.close()
        
#        Raspberry Pi hardware SPI configuration.
        self.SPI_PORT   = SPI_PORT
        self.SPI_DEVICE = SPI_DEVICE
        self.RelayPin = RelayPin 
        self.sensor = MAX31856(hardware_spi=Adafruit_GPIO.SPI.SpiDev(self.SPI_PORT, self.SPI_DEVICE), tc_type=MAX31856.MAX31856_K_TYPE) #remove tc_type for non-K thermocouples
        
        # Relay configuration.
        # relay = DOD(0) # Creates a general DigitalOutputDevice using GPIO pin 0
        
        # PWM configuration
        GPIO.setmode(GPIO.BCM)
        GPIO.setup(RelayPin, GPIO.OUT)
        self.p = GPIO.PWM(RelayPin, 2) #BCM pin 1 set to 0.5 Hz frequency PWM
        self.p.start(0) #Starter PWM med 0 av 100 i duty cycle
        
        self.pid = PID(self.P, self.I, self.D)
        self.pid.SetPoint = 20
        self.pid.setSampleTime(1)
                
        self.textbox.insert(tk.END,"Precursortester config file was loaded.\n")

    def set_output_folder(self):
        folder = filedialog.askdirectory() #path = Path('folder_of_camera_images')  # Path of folder containing jpg images
        self.output_folder = Path(folder)  # Path of folder containing jpg images
        self.textbox.insert(tk.END,str("Out-folder set to : "+str(self.output_folder)+"\n"))

    def plot_temp_profile(self):
        self.textbox.insert(tk.END,"Plotting temperature profile.\n")
        ID = 0
        max_ID = len(self.commands)-1 #Substract one to avoid including the End-command
     
        while ID < max_ID:
            x1,y1,x2,y2 = self.commands[ID].plot_val()
            self.x_profile.append(x1)
            self.y_profile.append(y1)
            self.x_profile.append(x2)
            self.y_profile.append(y2)
            ID = ID+1

        self.ax.plot(self.x_profile,self.y_profile,color='red',linewidth=1)        
        self.textbox.insert(tk.END,"Plotting done.\n")


    def pause_sequence(self):
        if(self.status == 2): # Program started
            self.textbox.insert(tk.END,"Paused the program!\n")
            self.status =3 #put it on pause
            self.buttonPause.config(text ="Paused", relief = "sunken")
            self.pause_time = time()
            self.run_time = time()-self.start_run_time
        elif(self.status == 3): # Program already paused
            self.textbox.insert(tk.END,"Unpaused the program!\n")
            self.buttonPause.config(text ="Pause", relief = "raised")
            self.status =2 #put it back on track
            t = self.start_run_time + (time() - self.pause_time) #fikser en ny init-tid som tar hensyn til at den har vært pauset
            self.start_run_time = t

    def stop_sequence(self):
        self.rt.stop()
##        self.rc.stop()
        self.f_out.close()
        self.textbox.insert(tk.END,"Stopped the program!\n")
        self.status=5

        self.buttonRun.config(state ="normal")
        self.buttonPause.config(state = "normal")
        self.buttonStop.config(state = "disabled")


    def set_temp(self, time):
        i = 0
        while time > self.commands[i].Time:
            i = i + 1
          
        temp = self.commands[i].set_temp(time)

        #Updates the Capture rate when going from one ID to another
        if self.commands[i].Capture_rate != self.Capture_rate:
            self.Capture_rate = self.commands[i].Capture_rate
            self.rc.interval = self.Capture_rate
            
        return temp

    def run_sequence(self): #denne kjøres ved en gitt frekvens. Standard = 1 Hz
        t = 0
        if(self.status == 2): #program started
            self.setT = self.set_temp(time()-self.start_run_time)
            self.pid.SetPoint = self.setT
            self.labelSetTempText.set(str('%.2f' % self.setT))
            self.labelTimeText.set(str('%.i' % (time()-self.start_run_time))) #Oppdaterer label på førstesiden
            self.textbox.insert(tk.END,"Running sequential = "+str(time()-self.start_run_time)+" s\r")
            if (time()-self.start_run_time) > self.max_time:
                self.stop_sequence()
            t = time()-self.start_run_time
        elif(self.status == 3): #Paused
            t = self.run_time
            self.donothing()
      
        #Les-temp:
        temp = self.sensor.read_temp_c()
        self.pid.update(temp) # Calculates PID-value
        self.PID_pwm_reg(temp)        
        #Plot i graf = skjer i animate
        #Plot i datafil
        self.f_out.write(str('%.2f' % t)+"\t"+str('%.2f' % self.setT)+"\t"+str('%.2f' % temp)+"\tNoCapt\n")        
        

    def time_lapse_photo(self):
        this_time = time()-self.start_run_time
        set_temp = self.set_temp(this_time)
        temp = self.sensor.read_temp_c()

        start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]
        timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
            start_time[1], str(start_time[0])[2:],
            start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.        

        self.camera.capture(self.output_folder+"Capture_at_"+timestamp+"_"+str('%.2f' % temp)+".jpg")
        self.textbox.insert(tk.END,"Captured an image at "+str(time()-self.start_run_time)+"s, with temperature at "+str(temp)+" C\n")

        #Plot i datafil
        self.f_out.write(str('%.5f' % this_time)+"\t"+str('%.2f' % set_temp)+"\t"+str('%.2f' % temp)+"\tCapt")        

    def rate(self,val,limit):
        if(val>limit):
            return 100
        elif(val<0):
            return 0
        else:
            return (val/limit)*100


    def PID_pwm_reg(self,temp):
        self.pid.update(temp) # Calculates PID-value
        
        val = self.rate(self.pid.output,10)
        self.oven.update(val,time())
        self.labelPowerText.set(str('%.1f' % val))
        #val = 100*max(min(pid.output, 1),0) # Turns PID-value into something between 0 and 10
        #print(("Val: {:.3f}").format(val))
        self.p.ChangeDutyCycle(val)


    def run_program(self):

        self.buttonRun.config(state ="disabled")
        self.buttonPause.config(state = "normal")
        self.buttonStop.config(state = "normal")

        self.textbox.insert(tk.END,"Beginning program.\n")
        max_ID = len(self.commands)
        self.max_time = self.commands[max_ID-1].Time
        self.textbox.insert(tk.END,"Expected duration = "+str(self.max_time/60)+" minutes\n")
        self.status=2 #Run program
        
        self.start_run_time = time() #computer start-time (since epoc)
        start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]
        timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
            start_time[1], str(start_time[0])[2:],
            start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.        

        self.f_out = open(self.output_folder+timestamp+".txt","w")
        self.f_out.write("Time\tSetTemp\tRealTemp\tCaptureNo\n")

        self.rt = RepeatedTimer(1,self.run_sequence) #Controls temp and log
##        self.rc = RepeatedTimer(self.Capture_rate,self.time_lapse_photo,init_time) #Controls camera


    def camera_preview_on(self):
        self.camera.start_preview()
        return 0
    
    def camera_preview_off(self):
        self.camera.stop_preview()
        return 0

    def camera_capture(self):
        temp = self.sensor.read_temp_c()
        start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]
        timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
            start_time[1], str(start_time[0])[2:],
            start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.        

        self.camera.capture(self.output_folder+"Capture_at_"+timestamp+"_"+str('%.2f' % temp)+".jpg")
        self.textbox.insert(tk.END,"Captured an image at "+str(timestamp)+", with temperature at "+str(temp)+" C\n")


    
    

# root window created. Here, that would be the only window, but
# you can later have windows within windows.
root = tk.Tk()
root.geometry("800x480")

#creation of an instance
app = Window(root)

#mainloop 
root.mainloop()  


