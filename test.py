# 几何法仿真验证 -*
from os import write
from ROBOT import ROBOT
from ik import GeneralKinematic
import numpy as np
import sim
import math
import time
import csv
import until
from matplotlib import pyplot as plt


def main():
    #定义坐标轴
    # fig = plt.figure()
    # ax1 = plt.axes(projection = '3d')
    # x = []
    # y = []
    # z = []
    # for t in range(1000):
    #     time = t * 0.01
    #     k = np.array([2*math.pi, 100, 100])
    #     initPoint = np.array([563, 0, 150.2])
    #     curvePostion = until.trajPlanning(time, k, initPoint)
    #     x.append(curvePostion[0])
    #     y.append(curvePostion[1])
    #     z.append(curvePostion[2])
    
    # ax1.scatter3D(x,y,z,cmap = 'Blues')
    # ax1.plot3D(x,y,z,'gray')
    # plt.show()    
    robot = ROBOT()
    jointStates = np.array([0, 30, 0, 60, 0, 90, 0])*math.pi/180
    robot.conditionUpdate(jointStates)
    trans,rot = robot.EndPose(robot.JointStates)
    transN,rotN = robot.Pose(robot.JointStates,7)
    print(trans)
    print(transN)





if __name__ == '__main__':
    main()