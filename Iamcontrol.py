# 几何法仿真验证 -*
from os import write
from ROBOT import ROBOT
from ik import GeneralKinematic
import numpy as np
import sim
import math
import time
import csv
from matplotlib import pyplot as plt
import until


q =  np.array([0, 0, 0, 0, 0, 0, 0])*math.pi/180

#关节极限
q_min_armc = np.array([-170, -120, -170, -120, -170, -120, -180])*math.pi/180
q_max_armc = np.array([170, 120, 170, 120, 170, 120, 180])*math.pi/180

#用于连接
sim.simxFinish(-1) # just in case, close all opened connections
clientID=sim.simxStart('127.0.0.1',19999,True,True,5000,5) # Connect to CoppeliaSim

#备注：
#矩阵乘法：np.dot   点乘: *    叉乘：np.cross

def main():
    #record
   # csvfile = open('simulatedata','wb')
   # writer = csv.writer(csvfile)
   # fig = plt.figure()
   # ax1 = plt.axes(projection = '3d')
    #connect parameter
    sim.simxFinish(-1) # just in case, close all opened connections
    clientID=sim.simxStart('127.0.0.1',19999,True,True,5000,5) # Connect to CoppeliaSim
    #===初始化工作==================
    # init visual robot
    visualRobot = ROBOT()
    visualRobot.JointStates = np.array([0, 30, 0, 60, 0, 90, 0])*math.pi/180
    visualRobot.p_RCM_goal = np.array([563, 0, 200])
    visualRobot.RCM_Set_Flag = True
    visualRobot.conditionUpdate(visualRobot.JointStates)

    transRecordx = []
    transRecordy = []
    transRecordz = []
    if clientID !=-1:
        print ('Connected to remote API server')
        #------用户给定-------------------------------------------------------------------------
        myrcm     = np.array([563, 0, 200]) #用户设定的RCM位置
        mytest    = np.array([0, 1220,340]) #用户设定的起始点
        #------获取仿真中的对象并进行相应的设置---------------------------------------------------------------------
        retu1, my_test  = sim.simxGetObjectHandle(clientID,'test',sim.simx_opmode_blocking)
        retu3, my_rcm   = sim.simxGetObjectHandle(clientID,'RCM',sim.simx_opmode_blocking)
        sim.simxSetObjectPosition(clientID, my_test, -1, mytest/100, sim.simx_opmode_blocking)
        sim.simxSetObjectPosition(clientID, my_rcm,   -1,  myrcm/100, sim.simx_opmode_blocking)
        jointHandle = []
        for i in range(visualRobot.JointStatesNum):
            joinName = 'joint'+ str(i+1)
            ret1, joint    = sim.simxGetObjectHandle(clientID, joinName,sim.simx_opmode_blocking)
            jointHandle.append(joint)
        ## q为关节角，先摆到初始位置

        t = 0
        dt = 0.01
        initPoint = np.array([563, 0, 150.2])
        k = np.array([2*math.pi, 10, 10])
        while(True):
            ## 时间更新
            t = t + dt
            print(t)
            if(t > 2):
                break
            goal = until.trajPlanning(t, k, initPoint)
            joint = visualRobot.RCM_PathPlan(visualRobot.JointStates,goal)
            visualRobot.conditionUpdate(joint)
            ## 求解并且给定求解出来的角度于仿真中的机械臂
            for i in range(visualRobot.JointStatesNum):
                sim.simxSetJointPosition(clientID, jointHandle[i], visualRobot.JointStates[i], sim.simx_opmode_streaming)  
            ## 记录机械臂末端的位置
                transRecord,rotRecord = visualRobot.EndPose(visualRobot.JointStates)
                transRecordx.append(transRecord[0])
                transRecordy.append(transRecord[1])
                transRecordz.append(transRecord[2])
        
      #  ax1.scatter3D(transRecordx,transRecordy,transRecordz,cmap = 'Blues')
       # ax1.plot3D(transRecordx,transRecordy,transRecordz,'gray')
        plt.show()    

    else:
        print ('Failed connecting to remote API server')



if __name__ == '__main__':
    main()