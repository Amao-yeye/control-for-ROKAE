#!/usr/bin/env python
'''
This module consists of the method of Transform, namely:
    Rotation_Z,Rotation_Y,Rotation_X
    Fixed_RotationXYZ,Inverse_Fixed_RotationXYZ
    Euler_RotationXYZ,Euler_RotationZYZ,Inverse_Euler_RotationZYZ
    TransformMatrix,DH_Transform
'''
import numpy as np
from math import atan2, cos, sin, sqrt, pi

def Rotation_Z(angle):
    R = np.matrix([
        [cos(angle), -sin(angle),           0],
        [sin(angle),  cos(angle),           0],
        [         0,           0,           1]
    ])
    return R

def Rotation_Y(angle):
    R = np.matrix([
        [ cos(angle),           0,  sin(angle)],
        [          0,           1,           0],
        [-sin(angle),           0,  cos(angle)]
    ])
    return R

def Rotation_X(angle):
    R = np.matrix([
        [         1,           0,           0],
        [         0,  cos(angle), -sin(angle)],
        [         0,  sin(angle),  cos(angle)]
    ])
    return R

def Fixed_RotationXYZ(tx,ty,tz):
    return Rotation_Z(tz)*Rotation_Y(ty)*Rotation_X(tx)

def Inverse_Fixed_RotationXYZ(R):
    ty = atan2(-R[2,0],sqrt(R[0,0]**2+R[1,0]**2))
    if ty!= pi/2 and ty!= -pi/2:
        tx = atan2(R[2,1]/cos(ty),R[2,2]/cos(ty))
        tz = atan2(R[1,0]/cos(ty),R[0,0]/cos(ty))
    elif ty == pi/2:
        tz = 0
        tx = atan2(R[0,1],R[1,1])
    else:
        tz = 0
        tx = -atan2(R[0,1],R[1,1])
    return tx,ty,tz


def Euler_RotationXYZ(tx,ty,tz):
    return Rotation_X(tx)*Rotation_Y(ty)*Rotation_Z(tz)

def Euler_RotationZYZ(tz1,ty,tz2):
    return Rotation_Z(tz1)*Rotation_Y(ty)*Rotation_Z(tz2)

def Inverse_Euler_RotationZYZ(R):
    ty = atan2(sqrt(R[2,0]**2+R[2,1]**2),R[2,2])
    if ty!= 0 and ty!= pi:
        tz1 = atan2(R[1,2]/sin(ty), R[0,2]/sin(ty))
        tz2 = atan2(R[2,1]/sin(ty),-R[2,0]/sin(ty))
    elif ty == 0:
        tz1 = 0
        tz2 = atan2(-R[0,1],R[1,1])
    else:
        tz1 = 0
        tz2 = atan2(R[0,1],-R[1,1])
    #print(tz1*180/pi,ty*180/pi,tz2*180/pi)
    return tz1,ty,tz2

def TransformMatrix(x,y,z,tx,ty,tz):
    R = Fixed_RotationXYZ(tx,ty,tz)
    t = np.matrix([[x],[y],[z]])
    T_up   = np.hstack((R,t))
    T_down = np.matrix([0,0,0,1])
    T = np.vstack((T_up,T_down))
    return T

def Inverse_TransformMatrix(T):
    trans = [T[0,3],T[1,3],T[2,3]]
    tx,ty,tz = Inverse_Fixed_RotationXYZ(T[0:3,0:3])
    rot = [tx,ty,tz]
    return trans,rot

def MergeTransformMatrix(trans,rot):
    T = TransformMatrix(trans[0],trans[1],trans[2],rot[0],rot[1],rot[2])
    return T

def DH_Transform(DH):
    alpha = DH[0]
    a = DH[1]
    theta = DH[2]
    d = DH[3]
    T = np.matrix([
        [           cos(theta),           -sin(theta),           0,              a],
        [sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha),  -sin(alpha)*d],
        [sin(theta)*sin(alpha), cos(theta)*sin(alpha),  cos(alpha),   cos(alpha)*d],
        [                    0,                     0,           0,              1]
    ])
    return T    

