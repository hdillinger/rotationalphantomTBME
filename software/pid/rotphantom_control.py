# Created by: Yu Wang (yu.wang5800@yahoo.com)
# Adapted by: Hannes Dillinger (hannes.dillinger@ptb.de)

# Form layout created with: PyQt5 UI code generator 5.15.10

# General imports
import time
import sys
import numpy as np
import csv

# Tachometer
import socket

# Mass flow controller
from pycomm3 import (CIPDriver, Services, USINT, UINT, REAL, STRING)

# PyQt for GUI
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import QObject, QTimer
from PyQt5.QtWidgets import (QApplication, QDialog, QSizePolicy, 
                             QTabWidget, QWidget, 
                             QFrame, QLayout, QFormLayout, QGridLayout, QHBoxLayout, 
                             QGroupBox, QSpinBox, QDoubleSpinBox,
                             QPushButton, QLineEdit, QPlainTextEdit, 
                             QLCDNumber, QLabel,
                             QMainWindow
                             )
from QLed import QLed
import pyqtgraph as pg


class PIDController(QObject):
    '''
    Class for PID control
    '''
    def __init__(self, kp, ki, kd):
        '''
        Initialize PID parameters and error values.
        Set threshold for Anti-Windup.
        '''
        # PID parameters
        self.kp = kp
        self.ki = ki
        self.kd = kd

        # Errors
        self.prev_error = 0.0
        self.errorsum = 0.0

        # Anti windup (is turned off when set to zero)
        self.windup = 200

    def update_pid(self, error, dt):
        self.errorsum += error
        self.integral = self.errorsum * dt
        self.integralterm = self.ki * self.integral
        derivative = (error - self.prev_error) / dt

        # anti windup for integral term
        self.antiwindup()
        # print(self.integralterm)

        # calculate output
        output_rev = self.kp * error + self.integralterm + self.kd * derivative
        revtoflow = 1 # in L = (L/min)/(1/min)
        output_flow = revtoflow * output_rev # in % * L/min

        # update error variable
        self.prev_error = error # in 1/min
        return output_flow
    
    def antiwindup(self):
        '''
        Anti Wind-Up.
        Method to reduce the effect of saturation to prevent overshoot (and oscillations).
        Only operates on the integral term.
        '''
        # stop the integral term if it is greater than threshold
        if self.windup != 0:
            if self.integralterm > self.windup:
                self.integralterm = self.windup
            elif self.integralterm < -self.windup:
                self.integralterm = -self.windup

class Ui_RotationalPhantomControl(QMainWindow):
    '''
    Main class
    '''
    def setupUi(self, RotationalPhantomControl):
        '''
        Create all required widgets and lay them out on GUI.
        '''
        # App Layout
        RotationalPhantomControl.setObjectName("RotationalPhantomControl")
        RotationalPhantomControl.resize(1143, 1012)
        self.formLayout_8 = QFormLayout(RotationalPhantomControl)
        self.formLayout_8.setObjectName("formLayout_8")

        # Groupbox Input
        self.groupBoxInput = QGroupBox(RotationalPhantomControl)
        self.groupBoxInput.setEnabled(True)
        self.groupBoxInput.setObjectName("groupBoxInput")
        self.gridLayout_5 = QGridLayout(self.groupBoxInput)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.lcdStartMR = QLCDNumber(self.groupBoxInput)
        self.lcdStartMR.setObjectName("lcdStartMR")
        self.gridLayout_5.addWidget(self.lcdStartMR, 2, 0, 1, 1)
        self.pushStartMR = QPushButton(self.groupBoxInput)
        self.pushStartMR.setObjectName("pushStartMR")
        self.gridLayout_5.addWidget(self.pushStartMR, 1, 0, 1, 1)
        self.pushStopMR = QPushButton(self.groupBoxInput)
        self.pushStopMR.setObjectName("pushStopMR")
        self.gridLayout_5.addWidget(self.pushStopMR, 3, 0, 1, 1)

        # Tab Revolutions
        self.tabWidgetInput = QTabWidget(self.groupBoxInput)
        self.tabWidgetInput.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        self.tabWidgetInput.setUsesScrollButtons(False)
        self.tabWidgetInput.setObjectName("tabWidgetInput")
        self.tabRevolutions = QWidget()
        self.tabRevolutions.setObjectName("tabRevolutions")
        self.formLayout_3 = QFormLayout(self.tabRevolutions)
        self.formLayout_3.setObjectName("formLayout_3")
        self.labelRev = QLabel(self.tabRevolutions)
        self.labelRev.setObjectName("labelRev")
        self.formLayout_3.setWidget(0, QFormLayout.LabelRole, self.labelRev)
        self.spinBoxRev = QSpinBox(self.tabRevolutions)
        self.spinBoxRev.setMaximum(900)
        self.spinBoxRev.setSingleStep(50)
        self.spinBoxRev.setValue(200)
        self.spinBoxRev.setObjectName("spinBoxRev")
        self.formLayout_3.setWidget(0, QFormLayout.FieldRole, self.spinBoxRev)
        self.labelSpeed = QLabel(self.tabRevolutions)
        self.labelSpeed.setObjectName("labelSpeed")
        self.formLayout_3.setWidget(1, QFormLayout.LabelRole, self.labelSpeed)
        self.lcdSpeed = QLCDNumber(self.tabRevolutions)
        self.lcdSpeed.setObjectName("lcdSpeed")
        self.formLayout_3.setWidget(1, QFormLayout.FieldRole, self.lcdSpeed)
        self.labelAcc = QLabel(self.tabRevolutions)
        self.labelAcc.setObjectName("labelAcc")
        self.formLayout_3.setWidget(2, QFormLayout.LabelRole, self.labelAcc)
        self.lcdAcc = QLCDNumber(self.tabRevolutions)
        self.lcdAcc.setObjectName("lcdAcc")
        self.formLayout_3.setWidget(2, QFormLayout.FieldRole, self.lcdAcc)
        self.tabWidgetInput.addTab(self.tabRevolutions, "")

        # Tab Speed
        self.tabSpeed = QWidget()
        self.tabSpeed.setObjectName("tabSpeed")
        self.formLayout_2 = QFormLayout(self.tabSpeed)
        self.formLayout_2.setObjectName("formLayout_2")
        self.labelRev_2 = QLabel(self.tabSpeed)
        self.labelRev_2.setObjectName("labelRev_2")
        self.formLayout_2.setWidget(0, QFormLayout.LabelRole, self.labelRev_2)
        self.lcdRev = QLCDNumber(self.tabSpeed)
        self.lcdRev.setSmallDecimalPoint(False)
        self.lcdRev.setDigitCount(5)
        self.lcdRev.setSegmentStyle(QLCDNumber.Filled)
        self.lcdRev.setProperty("value", 100.0)
        self.lcdRev.setProperty("intValue", 100)
        self.lcdRev.setObjectName("lcdRev")
        self.formLayout_2.setWidget(0, QFormLayout.FieldRole, self.lcdRev)
        self.labelSpeed_2 = QLabel(self.tabSpeed)
        self.labelSpeed_2.setObjectName("labelSpeed_2")
        self.formLayout_2.setWidget(1, QFormLayout.LabelRole, self.labelSpeed_2)
        self.labelAcc_2 = QLabel(self.tabSpeed)
        self.labelAcc_2.setObjectName("labelAcc_2")
        self.formLayout_2.setWidget(2, QFormLayout.LabelRole, self.labelAcc_2)
        self.lcdAcc_2 = QLCDNumber(self.tabSpeed)
        self.lcdAcc_2.setObjectName("lcdAcc_2")
        self.formLayout_2.setWidget(2, QFormLayout.FieldRole, self.lcdAcc_2)
        self.dbSpinBoxSpeed = QDoubleSpinBox(self.tabSpeed)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dbSpinBoxSpeed.sizePolicy().hasHeightForWidth())
        self.dbSpinBoxSpeed.setSizePolicy(sizePolicy)
        self.dbSpinBoxSpeed.setMaximum(4.0)
        self.dbSpinBoxSpeed.setSingleStep(1.0)
        self.dbSpinBoxSpeed.setObjectName("dbSpinBoxSpeed")
        self.formLayout_2.setWidget(1, QFormLayout.FieldRole, self.dbSpinBoxSpeed)
        self.tabWidgetInput.addTab(self.tabSpeed, "")

        # Tab Acceleration
        self.tabAcceleration = QWidget()
        self.tabAcceleration.setObjectName("tabAcceleration")
        self.formLayout_4 = QFormLayout(self.tabAcceleration)
        self.formLayout_4.setObjectName("formLayout_4")
        self.labelRev_3 = QLabel(self.tabAcceleration)
        self.labelRev_3.setObjectName("labelRev_3")
        self.formLayout_4.setWidget(0, QFormLayout.LabelRole, self.labelRev_3)
        self.lcdRev_2 = QLCDNumber(self.tabAcceleration)
        self.lcdRev_2.setObjectName("lcdRev_2")
        self.formLayout_4.setWidget(0, QFormLayout.FieldRole, self.lcdRev_2)
        self.labelSpeed_3 = QLabel(self.tabAcceleration)
        self.labelSpeed_3.setObjectName("labelSpeed_3")
        self.formLayout_4.setWidget(1, QFormLayout.LabelRole, self.labelSpeed_3)
        self.lcdSpeed_2 = QLCDNumber(self.tabAcceleration)
        self.lcdSpeed_2.setSmallDecimalPoint(False)
        self.lcdSpeed_2.setDigitCount(5)
        self.lcdSpeed_2.setSegmentStyle(QLCDNumber.Filled)
        self.lcdSpeed_2.setProperty("value", 100.0)
        self.lcdSpeed_2.setProperty("intValue", 100)
        self.lcdSpeed_2.setObjectName("lcdSpeed_2")
        self.formLayout_4.setWidget(1, QFormLayout.FieldRole, self.lcdSpeed_2)
        self.labelAcc_3 = QLabel(self.tabAcceleration)
        self.labelAcc_3.setObjectName("labelAcc_3")
        self.formLayout_4.setWidget(2, QFormLayout.LabelRole, self.labelAcc_3)
        self.dbSpinBoxAcc = QDoubleSpinBox(self.tabAcceleration)
        self.dbSpinBoxAcc.setObjectName("dbSpinBoxAcc")
        self.formLayout_4.setWidget(2, QFormLayout.FieldRole, self.dbSpinBoxAcc)
        self.tabWidgetInput.addTab(self.tabAcceleration, "")
        self.gridLayout_5.addWidget(self.tabWidgetInput, 5, 0, 1, 1)
        self.lcdStopMR = QLCDNumber(self.groupBoxInput)
        self.lcdStopMR.setObjectName("lcdStopMR")
        self.gridLayout_5.addWidget(self.lcdStopMR, 4, 0, 1, 1)
        self.formLayout_8.setWidget(0, QFormLayout.LabelRole, self.groupBoxInput)

        # Groupbox Control Settings
        self.groupBoxFlowControl = QGroupBox(RotationalPhantomControl)
        self.groupBoxFlowControl.setObjectName("groupBoxFlowControl")
        self.formLayout_8.setWidget(0, QFormLayout.FieldRole, self.groupBoxFlowControl)
        self.groupBoxControlSettings = QGroupBox(RotationalPhantomControl)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBoxControlSettings.sizePolicy().hasHeightForWidth())
        self.groupBoxControlSettings.setSizePolicy(sizePolicy)
        self.groupBoxControlSettings.setMinimumSize(QtCore.QSize(270, 0))
        self.groupBoxControlSettings.setMaximumSize(QtCore.QSize(400, 400))
        self.groupBoxControlSettings.setObjectName("groupBoxControlSettings")
        self.formLayout = QFormLayout(self.groupBoxControlSettings)
        self.formLayout.setSizeConstraint(QLayout.SetDefaultConstraint)
        self.formLayout.setContentsMargins(-1, -1, -1, 11)
        self.formLayout.setSpacing(10)
        self.formLayout.setObjectName("formLayout")
        self.labelSamplingTime = QLabel(self.groupBoxControlSettings)
        self.labelSamplingTime.setObjectName("labelSamplingTime")
        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.labelSamplingTime)
        self.spinBoxSamplingTime = QDoubleSpinBox(self.groupBoxControlSettings)
        sizePolicy = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.spinBoxSamplingTime.sizePolicy().hasHeightForWidth())
        self.spinBoxSamplingTime.setSizePolicy(sizePolicy)
        self.spinBoxSamplingTime.setMaximum(1199)
        self.spinBoxSamplingTime.setMinimum(0.1)
        self.spinBoxSamplingTime.setFixedWidth(100)
        self.spinBoxSamplingTime.setSingleStep(0.1)
        self.spinBoxSamplingTime.setValue(0.1)
        self.spinBoxSamplingTime.setObjectName("spinBoxSamplingTime")
        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.spinBoxSamplingTime)
        self.labelPhantomRadius = QLabel(self.groupBoxControlSettings)
        self.labelPhantomRadius.setObjectName("labelPhantomRadius")
        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.labelPhantomRadius)
        self.linePhantomRadius = QLineEdit(self.groupBoxControlSettings)
        self.linePhantomRadius.setFixedWidth(100)
        self.linePhantomRadius.setObjectName("linePhantomRadius")
        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.linePhantomRadius)
        # Kp
        self.label_Kp = QLabel(self.groupBoxControlSettings)
        self.label_Kp.setObjectName("label_Kp")
        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.label_Kp)
        self.dbSpinBox_Kp = QDoubleSpinBox(self.groupBoxControlSettings)
        self.dbSpinBox_Kp.setDecimals(3)
        self.dbSpinBox_Kp.setMinimum(0.0)
        self.dbSpinBox_Kp.setMaximum(20.0)
        self.dbSpinBox_Kp.setValue(3.0)
        self.dbSpinBox_Kp.setObjectName("dbSpinBox_Kp")
        self.dbSpinBox_Kp.setFixedWidth(100)
        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.dbSpinBox_Kp)
        # Ki
        self.label_Ki = QLabel(self.groupBoxControlSettings)
        self.label_Ki.setObjectName("label_Ki")
        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.label_Ki)
        self.dbSpinBox_Ki = QDoubleSpinBox(self.groupBoxControlSettings)
        self.dbSpinBox_Ki.setDecimals(3)
        self.dbSpinBox_Ki.setMinimum(0.0)
        self.dbSpinBox_Ki.setMaximum(10.0)
        self.dbSpinBox_Ki.setValue(0.5)
        self.dbSpinBox_Ki.setObjectName("dbSpinBox_Ki")
        self.dbSpinBox_Ki.setFixedWidth(100)
        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.dbSpinBox_Ki)
        # Kd
        self.label_Kd = QLabel(self.groupBoxControlSettings)
        self.label_Kd.setObjectName("label_Kd")
        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.label_Kd)
        self.dbSpinBox_Kd = QDoubleSpinBox(self.groupBoxControlSettings)
        self.dbSpinBox_Kd.setDecimals(3)
        self.dbSpinBox_Kd.setMinimum(0.0)
        self.dbSpinBox_Kd.setMaximum(10.0)
        self.dbSpinBox_Kd.setObjectName("dbSpinBox_Kd")
        self.dbSpinBox_Kd.setFixedWidth(100)
        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.dbSpinBox_Kd)
        # IP MFC
        self.label_ipAdress_MFC = QLabel(self.groupBoxControlSettings)
        self.label_ipAdress_MFC.setObjectName("label_ipAdress_MFC")
        self.formLayout.setWidget(5, QFormLayout.LabelRole, self.label_ipAdress_MFC)
        self.line_ip_MFC = QLineEdit(self.groupBoxControlSettings)
        self.line_ip_MFC.setObjectName("line_ip_MFC")
        self.line_ip_MFC.setFixedWidth(100)
        self.formLayout.setWidget(5, QFormLayout.FieldRole, self.line_ip_MFC)
        # IP Tacho
        self.label_ipAdress_Tacho = QLabel(self.groupBoxControlSettings)
        self.label_ipAdress_Tacho.setObjectName("label_ipAdress_Tacho")
        self.formLayout.setWidget(6, QFormLayout.LabelRole, self.label_ipAdress_Tacho)
        self.lineEdit_ip_Tacho = QLineEdit(self.groupBoxControlSettings)
        self.lineEdit_ip_Tacho.setObjectName("lineEdit_ip_Tacho")
        self.lineEdit_ip_Tacho.setFixedWidth(100)
        self.formLayout.setWidget(6, QFormLayout.FieldRole, self.lineEdit_ip_Tacho)
        # Port Tacho
        self.label_port_Tacho = QLabel(self.groupBoxControlSettings)
        self.label_port_Tacho.setObjectName("label_ipAdress_Tacho")
        self.formLayout.setWidget(7, QFormLayout.LabelRole, self.label_port_Tacho)
        self.lineEdit_port_Tacho = QLineEdit(self.groupBoxControlSettings)
        self.lineEdit_port_Tacho.setObjectName("lineEdit_port_Tacho")
        self.lineEdit_port_Tacho.setFixedWidth(100)
        self.formLayout.setWidget(7, QFormLayout.FieldRole, self.lineEdit_port_Tacho)
        # Set CW and CCW values
        self.CWbutton = QPushButton(self.groupBoxControlSettings)
        self.CWbutton.setObjectName("Use CW")
        self.CWbutton.setFixedWidth(100)
        self.formLayout.setWidget(8, QFormLayout.LabelRole, self.CWbutton)
        self.CCWbutton = QPushButton(self.groupBoxControlSettings)
        self.CCWbutton.setObjectName("Use CCW")
        self.CCWbutton.setFixedWidth(100)
        self.formLayout.setWidget(8, QFormLayout.FieldRole, self.CCWbutton)
        self.formLayout_8.setWidget(1, QFormLayout.LabelRole, self.groupBoxControlSettings) # Layout
        
        # Groupbox Speed (plot)
        self.groupBoxSpeed = QGroupBox(RotationalPhantomControl)
        self.groupBoxSpeed.setEnabled(True)
        self.groupBoxSpeed.setObjectName("groupBoxSpeed")
        self.formLayout_8.setWidget(1, QFormLayout.FieldRole, self.groupBoxSpeed)

        # Start, Stop, Save, Check buttons and LED
        self.frameButtons = QFrame(RotationalPhantomControl)
        self.frameButtons.setObjectName("frameButtons")
        self.horizontalLayout_2 = QHBoxLayout(self.frameButtons)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.STARTButton = QPushButton(self.frameButtons)
        self.STARTButton.setIconSize(QtCore.QSize(40, 20))
        self.STARTButton.setObjectName("STARTButton")
        self.horizontalLayout_2.addWidget(self.STARTButton)
        self.STOPButton = QPushButton(self.frameButtons)
        self.STOPButton.setObjectName("STOPButton")
        self.horizontalLayout_2.addWidget(self.STOPButton)
        self.line = QFrame(self.frameButtons)
        self.line.setFrameShape(QFrame.VLine)
        self.line.setFrameShadow(QFrame.Sunken)
        self.line.setObjectName("line")
        self.horizontalLayout_2.addWidget(self.line)
        self.lineEdit_savecsv = QLineEdit(self.frameButtons)
        self.lineEdit_savecsv.setObjectName("lineEdit_savecsv")
        self.lineEdit_savecsv.setFixedWidth(200)
        self.horizontalLayout_2.addWidget(self.lineEdit_savecsv)
        self.pushButton = QPushButton(self.frameButtons)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout_2.addWidget(self.pushButton)
        self.checkButtonmfc = QPushButton(self.frameButtons)
        self.line2 = QFrame(self.frameButtons)
        self.line2.setFrameShape(QFrame.VLine)
        self.line2.setFrameShadow(QFrame.Sunken)
        self.line2.setObjectName("line2")
        self.horizontalLayout_2.addWidget(self.line2)
        self.checkButtonmfc.setObjectName("checkButtonmfc")
        self.horizontalLayout_2.addWidget(self.checkButtonmfc)
        self.checkButtontacho = QPushButton(self.frameButtons)
        self.checkButtontacho.setObjectName("checkButtontacho")
        self.horizontalLayout_2.addWidget(self.checkButtontacho)
        self.led = QLed(self.frameButtons, onColour=QLed.Green, shape=QLed.Circle)
        self.led.value = False
        self.horizontalLayout_2.addWidget(self.led)
        self.msg = QLabel("")
        self.horizontalLayout_2.addWidget(self.msg)
        self.formLayout_8.setWidget(2, QFormLayout.SpanningRole, self.frameButtons)

        # Groupbox Control Information
        self.groupBoxControlInfo = QGroupBox(RotationalPhantomControl)
        self.groupBoxControlInfo.setObjectName("groupBoxControlInfo")
        self.gridLayout_9 = QGridLayout(self.groupBoxControlInfo)
        self.gridLayout_9.setObjectName("gridLayout_9")
        self.labelRevOut = QLabel(self.groupBoxControlInfo)
        self.labelRevOut.setObjectName("labelRevOut")
        self.gridLayout_9.addWidget(self.labelRevOut, 1, 0, 1, 1)
        self.lcdRevOut = QLCDNumber(self.groupBoxControlInfo)
        self.lcdRevOut.setObjectName("lcdRevOut")
        self.gridLayout_9.addWidget(self.lcdRevOut, 1, 1, 1, 1)
        self.labelAccOut = QLabel(self.groupBoxControlInfo)
        self.labelAccOut.setObjectName("labelAccOut")
        self.gridLayout_9.addWidget(self.labelAccOut, 1, 2, 1, 1)
        self.lcdAccOut = QLCDNumber(self.groupBoxControlInfo)
        self.lcdAccOut.setObjectName("lcdAccOut")
        self.gridLayout_9.addWidget(self.lcdAccOut, 1, 3, 1, 1)
        self.labelFlow = QLabel(self.groupBoxControlInfo)
        self.labelFlow.setObjectName("labelFlow")
        self.gridLayout_9.addWidget(self.labelFlow, 0, 0, 1, 1)
        self.lcdFlow = QLCDNumber(self.groupBoxControlInfo)
        self.lcdFlow.setObjectName("lcdFlow")
        self.gridLayout_9.addWidget(self.lcdFlow, 0, 1, 1, 1)
        self.labelSpeedOut = QLabel(self.groupBoxControlInfo)
        self.labelSpeedOut.setObjectName("labelSpeedOut")
        self.gridLayout_9.addWidget(self.labelSpeedOut, 0, 2, 1, 1)
        self.lcdSpeedOut = QLCDNumber(self.groupBoxControlInfo)
        self.lcdSpeedOut.setObjectName("lcdSpeedOut")
        self.gridLayout_9.addWidget(self.lcdSpeedOut, 0, 3, 1, 1)
        self.labelOperationTime = QLabel(self.groupBoxControlInfo)
        self.labelOperationTime.setObjectName("labelOperationTime")
        self.gridLayout_9.addWidget(self.labelOperationTime, 2, 2, 1, 1)
        self.lcdOperationTime = QLCDNumber(self.groupBoxControlInfo)
        self.lcdOperationTime.setObjectName("lcdOperationTime")
        self.gridLayout_9.addWidget(self.lcdOperationTime, 2, 3, 1, 1)
        self.labelTempOut = QLabel(self.groupBoxControlInfo)
        self.labelTempOut.setObjectName("labelTempOut")
        self.gridLayout_9.addWidget(self.labelTempOut, 2, 0, 1, 1)
        self.lcdTemp = QLCDNumber(self.groupBoxControlInfo)
        self.lcdTemp.setObjectName("lcdTemp")
        self.gridLayout_9.addWidget(self.lcdTemp, 2, 1, 1, 1)
        self.formLayout_8.setWidget(3, QFormLayout.FieldRole, self.groupBoxControlInfo)

        # Groupbox Step Testing
        self.groupBoxStepResponse = QGroupBox(RotationalPhantomControl)
        self.groupBoxStepResponse.setObjectName("groupBoxStepResponse")
        self.gridLayout_10 = QGridLayout(self.groupBoxStepResponse)
        self.gridLayout_10.setObjectName("gridLayout_10")

        self.StartStepTestbutton = QPushButton(self.frameButtons)
        self.StartStepTestbutton.setObjectName("StartStepTestbutton")
        self.StartStepTestbutton.setFixedWidth(122)
        self.gridLayout_10.addWidget(self.StartStepTestbutton, 0, 0, 1, 1)

        self.StopStepTestbutton = QPushButton(self.frameButtons)
        self.StopStepTestbutton.setObjectName("StopStepTestbutton")
        self.StopStepTestbutton.setFixedWidth(122)
        self.gridLayout_10.addWidget(self.StopStepTestbutton, 0, 1, 1, 1)

        self.StepNumber = QLabel(self.groupBoxStepResponse)
        self.StepNumber.setObjectName("StepNumber")
        self.gridLayout_10.addWidget(self.StepNumber, 1, 0, 1, 1)
        self.spinBoxStepNumber = QSpinBox(self.groupBoxStepResponse)
        self.spinBoxStepNumber.setObjectName("spinBoxStepNumber")
        self.spinBoxStepNumber.setValue(1)
        self.gridLayout_10.addWidget(self.spinBoxStepNumber, 1, 1, 1, 1)
        self.StepAmplitude = QLabel(self.groupBoxStepResponse)
        self.StepAmplitude.setObjectName("StepAmplitude")
        self.gridLayout_10.addWidget(self.StepAmplitude, 2, 0, 1, 1)
        self.spinBoxStepAmplitude = QSpinBox(self.groupBoxStepResponse)
        self.spinBoxStepAmplitude.setObjectName("spinBoxStepAmplitude")
        self.spinBoxStepAmplitude.setMinimum(0)
        self.spinBoxStepAmplitude.setMaximum(80)
        self.spinBoxStepAmplitude.setValue(40)
        self.spinBoxStepAmplitude.setSingleStep(2)
        self.gridLayout_10.addWidget(self.spinBoxStepAmplitude, 2, 1, 1, 1)
        self.formLayout_8.setWidget(3, QFormLayout.LabelRole, self.groupBoxStepResponse)

        self.lineEdit_savecsvstep = QLineEdit(self.frameButtons)
        self.lineEdit_savecsvstep.setObjectName("lineEdit_savecsvstep")
        self.lineEdit_savecsvstep.setFixedWidth(122)
        self.gridLayout_10.addWidget(self.lineEdit_savecsvstep, 3, 0, 1, 1)
        self.SaveStepTestbutton = QPushButton(self.frameButtons)
        self.SaveStepTestbutton.setObjectName("SaveStepTestbutton")
        self.gridLayout_10.addWidget(self.SaveStepTestbutton, 3, 1, 1, 1)

        # Notes
        self.plainNotesStep = QPlainTextEdit(RotationalPhantomControl)
        self.plainNotesStep.setObjectName("plainNotesStep")
        self.formLayout_8.setWidget(4, QFormLayout.LabelRole, self.plainNotesStep)
        self.plainNotesControl = QPlainTextEdit(RotationalPhantomControl)
        self.plainNotesControl.setObjectName("plainNotesControl")
        self.formLayout_8.setWidget(4, QFormLayout.FieldRole, self.plainNotesControl)

        self.retranslateUi(RotationalPhantomControl)
        self.tabWidgetInput.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(RotationalPhantomControl)

        ### MODIFICATIONS for GUI below ###
        # Add tacho and mfc plot widget
        self.start_flag = False

        self.gridLayout_plot = QGridLayout(self.groupBoxSpeed)
        self.plot_graph_tacho = pg.PlotWidget(self.groupBoxSpeed)
        self.gridLayout_plot.addWidget(self.plot_graph_tacho, 0, 0, 1, 1)

        self.gridLayout_plot = QGridLayout(self.groupBoxFlowControl)
        self.plot_graph_mfc = pg.PlotWidget(self.groupBoxFlowControl)
        self.gridLayout_plot.addWidget(self.plot_graph_mfc, 0, 0, 1, 1)

        # Connect Button widgets
        self.STARTButton.clicked.connect(self.START)
        self.STOPButton.clicked.connect(self.STOP)

        self.pushButton.clicked.connect(self.SAVE)
        self.checkButtonmfc.clicked.connect(self.check_mfc)
        self.checkButtontacho.clicked.connect(self.check_tacho)

        self.pushStartMR.setEnabled(False)
        self.pushStartMR.clicked.connect(self.startMR)
        self.pushStopMR.setEnabled(False)
        self.pushStopMR.clicked.connect(self.stopMR)

        self.StartStepTestbutton.clicked.connect(self.StartStepTesting)
        self.StopStepTestbutton.clicked.connect(self.StopStepTesting)
        self.SaveStepTestbutton.clicked.connect(self.SaveSteptest)

        self.CWbutton.clicked.connect(self.clockwise)
        self.CCWbutton.clicked.connect(self.counterclockwise)

        self.pushButton.setEnabled(False)
        self.checkButtonmfc.setEnabled(True)
        self.checkButtontacho.setEnabled(True)
        self.SaveStepTestbutton.setEnabled(False)

    def retranslateUi(self, RotationalPhantomControl):
        '''
        Internalization and localization
        '''
        _translate = QtCore.QCoreApplication.translate
        RotationalPhantomControl.setWindowTitle(_translate("RotationalPhantomControl", "Rotational Phantom Control"))
        self.groupBoxInput.setTitle(_translate("RotationalPhantomControl", "Input"))
        self.pushStartMR.setText(_translate("RotationalPhantomControl", "Set Marker: MRI started"))
        self.pushStopMR.setText(_translate("RotationalPhantomControl", "Set Marker: MRI stopped"))
        self.labelRev.setText(_translate("RotationalPhantomControl", "Number of revolutions [1/min]"))
        self.labelSpeed.setText(_translate("RotationalPhantomControl", "Speed [cm/s]"))
        self.labelAcc.setText(_translate("RotationalPhantomControl", "Acceleration [m/s^2]"))
        self.tabWidgetInput.setTabText(self.tabWidgetInput.indexOf(self.tabRevolutions), _translate("RotationalPhantomControl", "Revolution speed"))
        self.labelRev_2.setText(_translate("RotationalPhantomControl", "Number of revolutions [1/min]"))
        self.labelSpeed_2.setText(_translate("RotationalPhantomControl", "Speed [cm/s]"))
        self.labelAcc_2.setText(_translate("RotationalPhantomControl", "Acceleration [m/s^2]"))
        self.tabWidgetInput.setTabText(self.tabWidgetInput.indexOf(self.tabSpeed), _translate("RotationalPhantomControl", "Speed"))
        self.labelRev_3.setText(_translate("RotationalPhantomControl", "Number of revolutions [1/min]"))
        self.labelSpeed_3.setText(_translate("RotationalPhantomControl", "Speed [cm/s]"))
        self.labelAcc_3.setText(_translate("RotationalPhantomControl", "Acceleration [m/s^2]"))
        self.tabWidgetInput.setTabText(self.tabWidgetInput.indexOf(self.tabAcceleration), _translate("RotationalPhantomControl", "Acceleration"))
        self.groupBoxFlowControl.setTitle(_translate("RotationalPhantomControl", "Flow Control"))
        self.groupBoxControlSettings.setTitle(_translate("RotationalPhantomControl", "Control Settings"))
        self.labelSamplingTime.setText(_translate("RotationalPhantomControl", "Sampling time [s]"))
        self.labelPhantomRadius.setText(_translate("RotationalPhantomControl", "Phantom radius [m]"))
        self.linePhantomRadius.setText(_translate("RotationalPhantomControl", "0.0675"))
        self.label_Kp.setText(_translate("RotationalPhantomControl", "Proportional Gain Kp"))
        self.label_Ki.setText(_translate("RotationalPhantomControl", "Integral Gain Ki"))
        self.label_Kd.setText(_translate("RotationalPhantomControl", "Derivative Gain Kd"))
        # self.label_Ti.setText(_translate("RotationalPhantomControl", "Time constant Ti [s]"))
        self.label_ipAdress_MFC.setText(_translate("RotationalPhantomControl", "IP Adress Bürkert MFC 8741"))
        self.line_ip_MFC.setText(_translate("RotationalPhantomControl", "192.168.1.100"))
        self.label_ipAdress_Tacho.setText(_translate("RotationalPhantomControl", "IP Adress Monarch ACT-3X"))
        self.lineEdit_ip_Tacho.setText(_translate("RotationalPhantomControl", "192.168.1.102"))
        self.label_port_Tacho.setText(_translate("RotationalPhantomControl", "Port Monarch ACT-3X"))
        self.lineEdit_port_Tacho.setText(_translate("RotationalPhantomControl", "5000"))
        self.CWbutton.setText(_translate("RotationalPhantomControl", "Use CW"))
        self.CCWbutton.setText(_translate("RotationalPhantomControl", "Use CCW"))
        self.groupBoxSpeed.setTitle(_translate("RotationalPhantomControl", "Revolution Speed/Speed/Acceleration"))
        self.STARTButton.setText(_translate("RotationalPhantomControl", "START"))
        self.STOPButton.setText(_translate("RotationalPhantomControl", "STOP"))
        self.lineEdit_savecsv.setText(_translate("RotationalPhantomControl", "20240821_Nr01_CCW_200Umin"))
        self.pushButton.setText(_translate("RotationalPhantomControl", "Save to .csv file"))
        self.checkButtonmfc.setText(_translate("RotationalPhantomControl", "Check mfc"))
        self.checkButtontacho.setText(_translate("RotationalPhantomControl", "Check tacho"))
        self.groupBoxControlInfo.setTitle(_translate("RotationalPhantomControl", "Control Information"))
        self.labelRevOut.setText(_translate("RotationalPhantomControl", "Number of revolutions [1/min]"))
        self.labelAccOut.setText(_translate("RotationalPhantomControl", "Acceleration [m/s^2]"))
        self.labelFlow.setText(_translate("RotationalPhantomControl", "Flow [L/min]"))
        self.labelSpeedOut.setText(_translate("RotationalPhantomControl", "Speed [cm/s]"))
        self.labelOperationTime.setText(_translate("RotationalPhantomControl", "Operational time [s]"))
        self.labelTempOut.setText(_translate("RotationalPhantomControl", "Temperature [°C]"))
        self.groupBoxStepResponse.setTitle(_translate("RotationalPhantomControl", "Step Response Testing"))
        self.StartStepTestbutton.setText(_translate("RotationalPhantomControl", "Start Test"))
        self.SaveStepTestbutton.setText(_translate("RotationalPhantomControl", "Save Test Data"))
        self.StopStepTestbutton.setText(_translate("RotationalPhantomControl", "Stop Test"))
        self.lineEdit_savecsvstep.setText(_translate("RotationalPhantomControl", "CCW_240821"))
        self.StepNumber.setText(_translate("RotationalPhantomControl", "Test Number"))
        self.StepAmplitude.setText(_translate("RotationalPhantomControl", "Amplitude in Nl/min"))
        self.plainNotesStep.setPlainText(_translate("RotationalPhantomControl", "Workflow Step Response Testing:\n"
        "- Set Amplitude (Revolution Speed)\n"
        "- Set Test Number (+1)\n"
        "- Write CW or CCW\n"
        "- Click Start Test button\n"
        "- Experiment runs until Stop Test button clicked\n"
        "\n"
        "- Click Save Test Data to save to .csv file\n"
        "- Name of file is generated automatically with Format: Steptest_AmpXX_NrXX_CW.csv\n"
        ))
        self.plainNotesControl.setPlainText(_translate("RotationalPhantomControl", "Workflow Control:\n"
        "- Check IP and port addresses of devices\n"
        "- Check Phantom radius\n"
        "- Set Control Settings (PID control parameters and Sampling time)\n"
        "- Set Speed Input\n"
        "- Click Check tacho and Check mfc to check if they are connected properly\n"
        "- Click START\n"
        "- Click Set Marker for when MRI starts and stops\n"
        "- Click STOP to stop the Program\n"
        "\n"
        "- Write filename\n"
        "- Click Save to .csv file to save to .csv file\n"
        ))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.MRstarted = 0
        self.MRended = 0
        
        # Timer for LED - turns off after two seconds
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.turnoff_led)
    
    ## Tachometer functions
    def connect_tacho(self):
        '''
        Establishes Connection to Tachometer
        A Socket Object has to be created and opened.
        '''
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM) # create socket object
        self.sock.connect((self.address)) # establish connection

        command = "CH_1,TYPE,RPM,GATE,STD\r" # command byte string that is sent
        self.sock.sendto(command.encode('ascii'), self.address)

        time.sleep(0.1) # wait to receive answer
        data = self.sock.recv(1024) # set buffer size to get whole byte string
        
        print(data)

        print ("Connection to Tachometer successful.")
        print("****************************************")

    def disconnect_tacho(self):
        '''
        Closes Connection to Tachometer
        Socket Object created before has to be closed.
        '''
        self.sock.close()
        print ("Connection to Tachometer closed.")
        print("****************************************")

    def send_display_data(self):
        '''
        Sends the Command @D1 to Tachometer to
        send current displayed value once.
        '''
        # command = "CH_1,TYPE=RPM\r@D1\r" # command byte string that is sent
        command = "@D3\r" # command byte string that is sent
        self.sock.sendto(command.encode('ascii'), self.address)
        # time.sleep(0.1) # wait to receive answer
        data = self.sock.recv(1024) # set buffer size to get whole byte string
        print(data.decode('ascii').split('\r'))
        data_split = data.decode('ascii').split('\r')
        # data_decoded = float(data[-14:])/100  # last digits to get speed only, 8 stickers
        for data_point in data_split:
            try:
                # data_decoded = float(data_point)/60
                data_decoded = float(data_point)/20
                break
            except ValueError:
                data_decoded = None
        # data_decoded = float(data)  # last digits to get speed only, 8 stickers
        print(data_decoded)
        return data_decoded
    
    def set_display_rate(self, rate):
        '''
        NOT USED.
        Sends the Command DSPLR to Tachometer to
        Set the maximum display update rate to 0.5, 1 or 1.5 seconds.

        :param float rate: display rate (0.5, 1.0 or 1.5)
        '''
        if rate == 0.5:
            command3 = "DSPLR=HALF\r" #byte string
        elif rate == 1:
            command3 = "DSPLR=1_SEC\r"
        else:
            command3 = "DSPLR=1.5_S\r"
        self.sock.sendto(command3.encode('ascii'), self.address)
        time.sleep(0.1) # wait to receive answer
        self.data = self.sock.recv(32)

    ## MFC functions
    def read_acyclic(self, input):
        '''
        Read Acyclic data from MFC via IP address.
        Class, Instance and Attribute are passed.

        :param str input: Object to be read
        '''
        drive_path = self.mfc_ip
        data_type=UINT
        if input=="name":
            instance=b'\x08'
            data_type=STRING
            name="Device Name"
        elif input=="identnr":
            instance=b'\x09'
            name="Ident number"
        elif input=="serialnr":
            instance=b'\x0E'
            name="Serial Number"
        elif input=="temp":
            instance=b'\x12'
            data_type=REAL
            name="Device Temperature"
        elif input=="time":
            instance=b'\x14'
            name="Operation Time"

        with CIPDriver(drive_path) as drive:
            acyclic = drive.generic_message(
                service=Services.get_attribute_single, # b"\x0E"
                class_code=b'\xC7',
                instance=instance, 
                attribute=b'\x03',
                data_type=data_type,
                connected=False,
                #unconnected_send=True, # if connected=False and you have a route path
                route_path=True,
                name=name
            )

        value=acyclic.value
        if input=="temp": # Kelvin to Celsius
            value=acyclic.value-273.15
        return value
    
    def read_cyclic(self, input):
        '''
        Read Cyclic data from MFC via IP address.
        Class, Instance and Attribute are passed.

        :param str input: Object to be read
        '''
        drive_path = self.mfc_ip
        data_type=REAL
        class_code=b'\x64'
        if input=="status":
            instance=b'\x01'
            data_type=USINT
            name="NAMUR Status"
        elif input=="flow":
            instance=b'\x02'
            name="Actual Flow"
        elif input=="mtemp":
            instance=b'\x03'
            data_type=UINT
            name="Medium Temperature"
        elif input=="ctotalizer":
            instance=b'\x04'
            name="Current Totalizer"
        elif input=="csetpoint":
            instance=b'\x05'
            name="Current Setpoint"
        elif input=="setpoint":
            instance=b'\x01'
            class_code=b'\x65'
            data_type=UINT
            name="Setpoint"
        elif input=="cduty":
            instance=b'\x06'
            name="Current Dutycycle"

        with CIPDriver(drive_path) as drive:
            cyclic = drive.generic_message(
                service=Services.get_attribute_single, # b"\x0E"
                class_code=class_code,
                instance=instance, 
                attribute=b'\x03',
                data_type=data_type,
                connected=False,
                #unconnected_send=True, # if connected=False and you have a route path
                route_path=True,
                name=name
            )

        value=cyclic.value
        if input=="mtemp": # Kelvin to Celsius
            value=0.05*cyclic.value-273.
        return value
    
    def set_setpoint(self, input, setting):
        '''
        Sets Setpoint of MFC via IP address.
        Class, Instance and Attribute are passed.

        :param str input: Object to be set
        :param int/float setting: Value of Object
        '''
        drive_path = self.mfc_ip
        if input=="SetSPsource":
            instance=b'\x03'
            request_data=UINT.encode(setting)
            name='Set SP Source'
        elif input=="SetSPpercent":
            instance=b'\x04'
            request_data=REAL.encode(setting)
            name='Set SP Percent'

        with CIPDriver(drive_path) as drive:
            rw_controller = drive.generic_message(
                service=Services.set_attribute_single,
                class_code=b'\xC4',
                instance=instance, 
                attribute=b'\x03',
                request_data=request_data,
                # data_type=data_type,
                connected=True,
                #unconnected_send=True, #only needed if connected=False and you have a route path
                route_path=True,
                name=name
            )        
        value=rw_controller.value
        return value

    def read_setpoint(self, input):
        '''
        Reads Setpoint of MFC via IP address.
        Class, Instance and Attribute are passed.

        :param str input: Object to be read
        '''
        drive_path = self.mfc_ip
        if input=="SPsource":
            instance=b'\x03'
            data_type=UINT
            name='SP Source'
        elif input=="SPpercent":
            instance=b'\x04'
            data_type=REAL
            name='SP Percent'

        with CIPDriver(drive_path) as drive:
            rw_controller = drive.generic_message(
                service=Services.get_attribute_single,
                class_code=b'\xC4',
                instance=instance, 
                attribute=b'\x03',
                data_type=data_type,
                connected=True,
                #unconnected_send=True, # if connected=False and you have a route path
                route_path=True,
                name=name
            )        
        value=rw_controller.value
        return value

    ## Check connections
    def check_mfc(self):
        '''
        Check if MFC is properly connected and displays error if not.
        '''
        # Read IP Adress and Port
        self.mfc_ip = self.line_ip_MFC.text() # ip on computer: 192.168.1.103
        # self.mfc_port = 502
        self.checkButtonmfc.setEnabled(False)
        self.checkButtontacho.setEnabled(False)
        try:
            list = CIPDriver.list_identity(self.mfc_ip)
            print("****************************************")
            print("Connection to MFC successful.")
            # LED
            self.led.setOnColour(QLed.Green)
            self.led.value = True
            self.timer.start(1000)
            self.msg.setText("Successful.")
        except Exception as e:
            print(f"Connection to MFC failed: {e}")
            # LED
            self.led.setOnColour(QLed.Red)
            self.led.value = True
            self.timer.start(1000)
            self.msg.setText("Failed.")
        print("****************************************")

    def check_tacho(self):
        '''
        Check if tachometer is properly connected and displays error if not.
        '''
        # Read IP Adress and Port of Tachometer
        self.tacho_ip = self.lineEdit_ip_Tacho.text() # ip on computer: 192.168.1.101
        self.tacho_port = int(self.lineEdit_port_Tacho.text())
        self.address = (self.tacho_ip, self.tacho_port)
        self.checkButtonmfc.setEnabled(False)
        self.checkButtontacho.setEnabled(False)
        try:
            # (socket.AF_INET) specifies the IP address family for IPv4
            # (socket.SOCK_STREAM) indicates that we are using a TCP socket
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock: 
                # establish connection to a remote server
                sock.connect(self.address) 
            print("****************************************")
            print(f"Connection to tachometer successful.")
            # LED
            self.led.setOnColour(QLed.Green)
            self.led.value = True
            self.timer.start(1000)
            self.msg.setText("Successful.")
        except socket.error as e:
            print(f"Connection to tachometer failed: {e}")
            print("Make sure PmRemote Software is disconnected!")
            # LED
            self.led.setOnColour(QLed.Red)
            self.led.value = True
            self.timer.start(1000)
            self.msg.setText("Failed.")
        print("****************************************")

    def turnoff_led(self):
        self.led.value = False
        self.timer.stop()
        self.msg.setText("")
        self.checkButtonmfc.setEnabled(True)
        self.checkButtontacho.setEnabled(True)

    ## Main program
    def START(self):
        '''
        Establish Connection to all devices.
        Set Setpoint Source of MFC.
        Start Program including PID Control and Plotting.
        '''
        # while loop variable
        self.start_flag = True

        # disable Buttons/Boxes
        self.spinBoxRev.setEnabled(False)
        self.CWbutton.setEnabled(False)
        self.CCWbutton.setEnabled(False)
        self.spinBoxSamplingTime.setEnabled(False)
        self.dbSpinBox_Kp.setEnabled(False)
        self.dbSpinBox_Ki.setEnabled(False)
        self.dbSpinBox_Kd.setEnabled(False)
        self.pushButton.setEnabled(False)
        self.checkButtonmfc.setEnabled(False)
        self.checkButtontacho.setEnabled(False)
        self.SaveStepTestbutton.setEnabled(False)
        # enable Buttons/Boxes
        self.pushStartMR.setEnabled(True)
        self.pushStopMR.setEnabled(True)

        # Read phantom radius and convert Input to display in lcd
        self.radius = float(self.linePhantomRadius.text())
        revtospeed = 100*2*np.pi*self.radius/60 #cm/s
        revtoacc = 4*np.pi*np.pi*self.radius/(60*60) #m/s
        self.lcdSpeed.display(self.spinBoxRev.value()*revtospeed)
        self.lcdAcc.display(self.spinBoxRev.value()*self.spinBoxRev.value()*revtoacc)

        # Read IP Adress and Port of Tachometer
        self.tacho_ip = self.lineEdit_ip_Tacho.text() # ip on computer: 192.168.1.101
        self.tacho_port = int(self.lineEdit_port_Tacho.text())
        self.address = (self.tacho_ip, self.tacho_port)
        # Read IP Adress and Port of MFC
        self.mfc_ip = self.line_ip_MFC.text() # ip on computer: 192.168.1.103
        # self.mfc_port = 502

        # Sampling time is 0.1 without pausing or else subtract 0.1
        self.dt = self.spinBoxSamplingTime.value()
        self.pausedt = self.dt - 0.1 # tacho pausing time

        # Clear plots and values
        self.plot_graph_tacho.clear()
        self.plot_graph_mfc.clear()
        self.lcdStartMR.display(0.0)
        self.lcdStopMR.display(0.0)

        # Modify tacho and mfc plot visuals
        self.plot_graph_tacho.setBackground("k")
        self.plot_graph_tacho.addLegend()
        self.plot_graph_tacho.enableAutoRange()
        plotItem_tacho = self.plot_graph_tacho.getPlotItem()
        plotItem_tacho.showGrid(x=True, y=True, alpha=0.5)

        self.plot_graph_mfc.setBackground("k")
        self.plot_graph_mfc.addLegend()
        self.plot_graph_mfc.enableAutoRange()
        plotItem_mfc = self.plot_graph_mfc.getPlotItem()
        plotItem_mfc.showGrid(x=True, y=True, alpha=0.5)

        # Read control parameter values
        self.pid = PIDController(self.dbSpinBox_Kp.value(), self.dbSpinBox_Ki.value(), self.dbSpinBox_Kd.value())

        # Set Setpoint Source to Stored SP Percent (2)
        self.set_setpoint('SetSPsource', 2)
        # Establish connection to tachometer (socket object)
        self.connect_tacho()
        # Start plots
        self.display_plots()
        

    def STOP(self):
        '''
        Stop Program including Control and Plots.
        Close Tachometer connection and close MFC to zero.
        '''
        self.start_flag = False
        self.start_step = False
        # enable Buttons/Boxes
        self.spinBoxRev.setEnabled(True)
        self.CWbutton.setEnabled(True)
        self.CCWbutton.setEnabled(True)
        self.spinBoxSamplingTime.setEnabled(True)
        self.dbSpinBox_Kp.setEnabled(True)
        self.dbSpinBox_Ki.setEnabled(True)
        self.dbSpinBox_Kd.setEnabled(True)
        self.pushButton.setEnabled(True)
        self.checkButtonmfc.setEnabled(True)
        self.checkButtontacho.setEnabled(True)
        self.SaveStepTestbutton.setEnabled(True)
        # disable Buttons/Boxes
        self.pushStartMR.setEnabled(False)
        self.pushStopMR.setEnabled(False)

        # Disconnect Tachometer
        self.disconnect_tacho()

        # Close MFC
        self.set_setpoint('SetSPpercent', 0)

    def SAVE(self):
        '''
        Save data from both plots to csv file.
        Name of file taken from user input.
        '''
        firstrow = ['timestamps', 'Target Speed', 'Current Speed', 'Setpoint Flow', 'Current Flow', 'MR Start time', 'MR End time']
        rows = [self.timestamps,
                self.data_points_tacho_target,
                self.data_points_tacho,
                self.data_points_mfc,
                self.data_points_mfc_flow]
        lastrow = [self.MRstarted, self.MRended]

        with open("Software/GUI/textfiles/" + self.lineEdit_savecsv.text() + ".csv", "w") as file:
            write = csv.writer(file)
            write.writerow(firstrow)
            write.writerows(rows)
            write.writerow(lastrow)


    ## Plotting
    def display_plots(self):
        '''
        Open MFC to 100.
        Plot speed information of the Rotational Phantom received from tachometer and air flow through the MFC.
        '''
        # Initialize empty lists and values
        self.timestamps = []
        self.data_points_tacho = []
        self.data_points_tacho_target = []
        self.data_points_mfc = []
        self.data_points_mfc_flow = []

        self.timestamps.append(0)
        self.data_points_tacho.append(0)
        self.data_points_tacho_target.append(self.spinBoxRev.value())
        self.data_points_mfc.append(0)
        self.data_points_mfc_flow.append(0)

        # Open mfc to Start Phantom rotation
        self.set_setpoint('SetSPpercent', 100)

        # Plot Legend/Labels
        tacho1 = self.plot_graph_tacho.plot(self.timestamps[-100:], self.data_points_tacho_target[-100:], pen='r', name="Target Rotational Speed in RPM")
        tacho2 = self.plot_graph_tacho.plot(self.timestamps[-100:], self.data_points_tacho[-100:], name="Current Rotational Speed in RPM")
        mfc1 = self.plot_graph_mfc.plot(self.timestamps[-100:], self.data_points_mfc[-100:], pen='r', name="Setpoint Flow in Nl/min")
        mfc2 = self.plot_graph_mfc.plot(self.timestamps[-100:], self.data_points_mfc_flow[-100:], name="Current Flow in Nl/min")
        # Start time
        self.start_time = time.time()

        data_diff = []
        while self.start_flag:
            data_decoded = self.send_display_data()
            #data_decoded = self.data_points_tacho[-1]

            # # update only if last value is valid
            if data_decoded is not None:
                # data_decoded = np.mean(np.array(self.data_points_tacho[-1:-4:-1] + [data_decoded]))
                # data_decoded = self.data_points_tacho[-1]
            #     # Check for unrealistic speed data/steps
            #     # data_diff.append(abs(data_decoded-self.data_points_tacho[-1]))
                data_diff = abs(data_decoded-self.data_points_tacho[-1])
                # if len(data_diff) > 5:
                # if data_diff > 10:
                if data_diff > 20 and len(self.data_points_tacho) > 50:
                #     # data_decoded = self.data_points_tacho[-1]
                #     # data_diff = []
                    data_decoded = np.mean(np.array(self.data_points_tacho[-1:-2:-1]))
                # else:
                #     pass
            
                

                # Record timestamp and received data
                self.timestamps.append(time.time() - self.start_time)
                self.data_points_tacho.append(data_decoded)
                self.data_points_tacho_target.append(self.spinBoxRev.value())
                self.data_points_mfc.append(self.read_setpoint('SPpercent')/100*80)
                self.data_points_mfc_flow.append(self.read_cyclic('flow'))

            # Plot
            #self.plot_graph_tacho.plot(self.timestamps[-10:], self.data_points_tacho_target[-10:], pen='r')
            #self.plot_graph_tacho.plot(self.timestamps[-10:], self.data_points_tacho[-10:])
            #self.plot_graph_mfc.plot(self.timestamps[-10:], self.data_points_mfc[-10:], pen='r')
            #self.plot_graph_mfc.plot(self.timestamps[-10:], self.data_points_mfc_flow[-10:])
            
            tacho1.setData(self.timestamps[-100:], self.data_points_tacho_target[-100:])
            tacho2.setData(self.timestamps[-100:], self.data_points_tacho[-100:])
            mfc1.setData(self.timestamps[-100:], self.data_points_mfc[-100:])
            mfc2.setData(self.timestamps[-100:], self.data_points_mfc_flow[-100:])

            # Control Output
            lastFlow = self.data_points_mfc_flow[-1] #L/min
            lastRev = self.data_points_tacho[-1] #1/min
            lastTemp = self.read_acyclic('temp') #°C
            lastSpeed = 100*lastRev*(2*np.pi*float(self.linePhantomRadius.text()))/60 #cm/s
            lastAcc = lastRev*lastRev*4*np.pi*np.pi*float(self.linePhantomRadius.text())/(60*60)
            OpTime = self.timestamps[-1]

            self.lcdFlow.display(lastFlow)
            self.lcdSpeedOut.display(lastSpeed)
            self.lcdRevOut.display(lastRev)
            self.lcdTemp.display(lastTemp)
            self.lcdAccOut.display(lastAcc)
            self.lcdOperationTime.display(OpTime)

            # Update PID
            self.update_mfc_sp()

            # Pause for sampling time minus 0.1 (time that tacho sleeps until answer received)
            if self.pausedt != 0:
                time.sleep(self.pausedt)
            else:
                pass

            # Update plot
            QApplication.processEvents()

    ## PID Control
    def clockwise(self):
        '''
        Set the PID parameters to standard values for clockwise driven phantom.
        The values were determined through step testing.
        '''
        self.dbSpinBox_Kp.setValue(3)
        self.dbSpinBox_Ki.setValue(0.5)
        self.dbSpinBox_Kd.setValue(0)
    
    def counterclockwise(self):
        '''
        Set the PID parameters to standard values for counterclockwise driven phantom.
        The values were determined through step testing.
        '''
        self.dbSpinBox_Kp.setValue(2)
        self.dbSpinBox_Ki.setValue(0.4)
        self.dbSpinBox_Kd.setValue(0)
    
    def update_mfc_sp(self):
        '''
        Update MFC Setpoint value using PID control.
        '''
        # Read input of PID parameters and sampling time
        self.pid.kp = self.dbSpinBox_Kp.value()
        self.pid.ki = self.dbSpinBox_Ki.value()
        self.pid.kd = self.dbSpinBox_Kd.value()
        dt = self.dt

        # Read desired speed
        self.target_rev = self.spinBoxRev.value()

        # Get last speed value
        # self.actual_rev = self.data_points_tacho[-1]
        self.actual_rev = np.mean(self.data_points_tacho[-1:-4:-1])
        

        # Calculate speed error and get clipped control output for MFC
        error = (self.target_rev - self.actual_rev)
        control_signal = self.pid.update_pid(error, dt)
        control_signal_clipped = np.clip(control_signal, 0, 100)

        # Open MFC according to control output
        self.set_setpoint('SetSPpercent', control_signal_clipped)
        
    ## MR Marker
    def displayMRstart(self):
        '''
        Mark MRI start time by clicking Start MRI button.
        Display the timestamp and plot a vertical line.
        '''
        self.MRstarted = self.timestamps[-1]
        self.lcdStartMR.display(self.MRstarted)
        self.plot_graph_tacho.addLine(x=self.MRstarted)
        self.plot_graph_mfc.addLine(x=self.MRstarted)

    def startMR(self):
        self.displayMRstart()

    def displayMRstop(self):
        '''
        Mark MRI end time by clicking Stop MRI button.
        Display the timestamp and plot a vertical line.
        '''
        self.MRended = self.timestamps[-1]
        self.lcdStopMR.display(self.MRended)
        self.plot_graph_tacho.addLine(x=self.MRended)
        self.plot_graph_mfc.addLine(x=self.MRended)

    def stopMR(self):
        self.displayMRstop()

    ## Step testing
    def StartStepTesting(self):
        '''
        Starts Step Response Testing.
        '''
        # while loop variable
        self.start_step = True

        # disable Buttons/Boxes
        self.spinBoxRev.setEnabled(False)
        self.spinBoxStepAmplitude.setEnabled(False)
        self.pushStartMR.setEnabled(False)
        self.pushStopMR.setEnabled(False)
        self.SaveStepTestbutton.setEnabled(False)
        self.checkButtonmfc.setEnabled(False)
        self.checkButtontacho.setEnabled(False)
        # enable Buttons/Boxes

        # Read phantom radius
        self.radius = float(self.linePhantomRadius.text())

        # Read IP addresses and ports
        self.tacho_ip = self.lineEdit_ip_Tacho.text() # ip on computer: 192.168.1.101
        self.tacho_port = int(self.lineEdit_port_Tacho.text())
        self.address = (self.tacho_ip, self.tacho_port)
        self.mfc_ip = self.line_ip_MFC.text() # ip on computer: 192.168.1.103
        # self.mfc_port = 502

        # Clear plots and values
        self.plot_graph_tacho.clear()
        self.plot_graph_mfc.clear()
        self.lcdStartMR.display(0.0)
        self.lcdStopMR.display(0.0)

        # Modify plot visuals
        self.plot_graph_tacho.setBackground("k")
        self.plot_graph_tacho.addLegend()
        plotItem_tacho = self.plot_graph_tacho.getPlotItem()
        plotItem_tacho.showGrid(x=True, y=True, alpha=0.5)

        self.plot_graph_mfc.setBackground("k")
        self.plot_graph_mfc.addLegend()
        plotItem_mfc = self.plot_graph_mfc.getPlotItem()
        plotItem_mfc.showGrid(x=True, y=True, alpha=0.5)

        # Set Setpoint Source to Stored SP Percent (2)
        self.set_setpoint('SetSPsource', 2)
        # Establish connection to tachometer
        self.connect_tacho()

        # Read Step Amplitude and open MFC (Maximum 80L/min,...)
        self.amplitude = 100*self.spinBoxStepAmplitude.value()/80 # in %
        self.set_setpoint('SetSPpercent', self.amplitude)

        # Start plots
        self.display_StepTesting()

    def display_StepTesting(self):
        '''
        Plot Step Response Testing Data.
        '''
        # Initialize empty lists and values
        self.steptest_timestamps = []
        self.steptest_tacho = []
        self.steptest_mfc = []
        self.steptest_mfc_flow = []

        self.steptest_timestamps.append(0)
        self.steptest_tacho.append(0)
        self.steptest_mfc.append(self.read_setpoint('SPpercent')/100*80)
        self.steptest_mfc_flow.append(0)

        # Plot Legend/Labels
        self.plot_graph_tacho.plot(self.steptest_timestamps, self.steptest_tacho, name="Current Rotational Speed in RPM")
        self.plot_graph_mfc.plot(self.steptest_timestamps, self.steptest_mfc, pen='r', name="Setpoint Flow in Nl/min")
        self.plot_graph_mfc.plot(self.steptest_timestamps, self.steptest_mfc_flow, name="Current Flow in Nl/min")

        # Start time
        self.start_time = time.time()

        while self.start_step:
            data_decoded = self.send_display_data()
            
            
            if data_decoded is not None:

                # Record timestamp and received data
                self.steptest_timestamps.append(time.time() - self.start_time)
                self.steptest_tacho.append(data_decoded)
                self.steptest_mfc.append(self.read_setpoint('SPpercent')/100*80)
                self.steptest_mfc_flow.append(self.read_cyclic('flow'))

                # Plot
                self.plot_graph_tacho.plot(self.steptest_timestamps, self.steptest_tacho)
                self.plot_graph_mfc.plot(self.steptest_timestamps, self.steptest_mfc, pen='r')
                self.plot_graph_mfc.plot(self.steptest_timestamps, self.steptest_mfc_flow)
                
                # Control Output
                lastFlow = self.steptest_mfc_flow[-1] #L/min
                lastRev = self.steptest_tacho[-1] #1/min
                lastTemp = self.read_acyclic('temp') #°C
                lastSpeed = 100*lastRev*(2*np.pi*self.radius)/60 #cm/s
                lastAcc = 0
                OpTime = self.steptest_timestamps[-1]

                self.lcdFlow.display(lastFlow)
                self.lcdSpeedOut.display(lastSpeed)
                self.lcdRevOut.display(lastRev)
                self.lcdTemp.display(lastTemp)
                self.lcdAccOut.display(lastAcc)
                self.lcdOperationTime.display(OpTime)

            QApplication.processEvents() # update plot

    def StopStepTesting(self):
        '''
        Stop Step Testing including Plots.
        Close Tachometer connection and close MFC to zero.
        '''
        self.start_flag = False
        self.start_step = False

        # enable Buttons/Boxes
        self.spinBoxStepAmplitude.setEnabled(True)
        self.spinBoxRev.setEnabled(True)
        self.SaveStepTestbutton.setEnabled(True)
        self.checkButtonmfc.setEnabled(True)
        self.checkButtontacho.setEnabled(True)
        # disable Buttons/Boxes

        # Disconnect
        self.disconnect_tacho()

        # Close MFC - Setpoint to zero percent
        self.set_setpoint('SetSPpercent', 0)

    def SaveSteptest(self):
        '''
        Save data from both plots to csv file.
        Name of file partially taken from user input.
        '''
        firstrow_step = ['timestamps', 'Current Speed', 'Setpoint Flow', 'Current Flow']
        rows_step = [self.steptest_timestamps,
                    self.steptest_tacho,
                    self.steptest_mfc,
                    self.steptest_mfc_flow]

        with open("Software/GUI/steptests/Steptest_Amp" + str(self.spinBoxStepAmplitude.value()) + "_Nr" + str(self.spinBoxStepNumber.value()) + self.lineEdit_savecsvstep.text() + ".csv", "w") as file_step:
            write = csv.writer(file_step)
            write.writerow(firstrow_step)
            write.writerows(rows_step)


##########################################
# Execute Application
if __name__ == '__main__':

    pg.setConfigOption('leftButtonPan', False)

    # Create object of QApplication class
    app = QApplication(sys.argv) # sys.argv List of command line parameters
    Dialog = QDialog()

    # Set up font size to 14pt for all widgets
    custom_font = QtGui.QFont()
    custom_font.setPointSize(14)

    # Create window
    ui = Ui_RotationalPhantomControl()
    ui.setupUi(Dialog)
    Dialog.show()

    # Run application until user closes it
    sys.exit(app.exec_())