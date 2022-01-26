#!/usr/bin/env python
# license removed for brevity
from time import sleep
import xml.dom.minidom

from numpy.core.fromnumeric import squeeze
import Transform
import numpy as np
import math
import until

# from urdf_parser_py.urdf import URDF
# from pykdl_utils.kdl_kinematics import KDLKinematics
# from pykdl_utils.kdl_parser import kdl_tree_from_urdf_model

class ROBOT:
    def __init__(self,Lambda = 0.5):
        # Count the number of continuous joint, in order to know the dimension of joint states
        self.JointStatesNum = 7
        # DH parameter
        self.DH = np.array([[0, 0, 0, 341.5],
                    [-math.pi/2, 0, 0, 0],
                    [ math.pi/2, 0, 0, 394],
                    [-math.pi/2, 0, 0, 0],
                    [math.pi/2, 0, 0, 366],
                    [-math.pi/2, 0, 0, 0],
                    [math.pi/2, 0, 0, 250.3]])
        self.endEffector = np.array([0,0,0,282.2])
        # Robot JointState
        self.JointStates = [0,0,0,0,0,0,0]
        # Get the basic transfome matirx for robot, and store in a list 
        self.BasicTransformMatirx = self.GetBasicBasicTransformMatirx()
        # Get the endEffector transfome matirx for robot
        self.endEffectorMatirx = self.GetEndeffectorTransformMatirx()
        ## eps, used for numerically differentiation
        self.eps = 0.001
        ## tolerence, used for determining if a path is accurate enough
        self.tolerance = 0.02
        ## limit velocity of the robot
        self.VelocityLimit = 0.1
        ## Lambda for RCM setting, default set to be 2
        self.Lambda = Lambda
        ## Max iterative for PathPlanner
        self.max_iter = 1000
        ## iterative recorder, used for online PathPlan
        self.iter = 0
        ## Flag to mark if RCM has been set
        self.RCM_Set_Flag = 0
        ## Max intergration, prevent Overflow
        self.Inte_max = 1
        ## PID Control, recording error is necessary
        self.ErrorSum = None
        self.LastError = None
        # RCM position in grobal
        self.p_RCM_goal = np.zeros(3)

        self.TrackingOver = 0
    
    def conditionUpdate(self,JointStates):
        self.JointStates = JointStates
        self.BasicTransformMatirx = self.GetBasicBasicTransformMatirx()

    def GetBasicBasicTransformMatirx(self):
        BasicTransformMatirxList = []
        for joint in range(self.JointStatesNum):
            TM = Transform.DH_Transform(self.DH[joint]+[0, 0, self.JointStates[joint], 0])
            BasicTransformMatirxList.append(TM)
        return BasicTransformMatirxList
    
    def GetEndeffectorTransformMatirx(self):
        EndeffectorTransformMatirx = Transform.DH_Transform(self.endEffector)
        return EndeffectorTransformMatirx
    
    def Pose(self,JointStates,FrameIndex):#LI：返回命令的坐标系的位姿
        if len(JointStates)!=self.JointStatesNum:
            print("dismatch JointStatesNum!!!")
            return []
        pose = np.eye(4)
        self.BasicTransformMatirx = self.GetBasicBasicTransformMatirx()
        for i in range(FrameIndex):
            pose = np.dot(pose,self.BasicTransformMatirx[i])
        trans,rot = Transform.Inverse_TransformMatrix(pose)
        return trans,rot
        #return pose 4X4 matrix
    
    def EndPose(self,JointStates):
        trans,rot = self.Pose(JointStates,self.JointStatesNum)
        T = Transform.MergeTransformMatrix(trans,rot)##LI：根据旋转和平移得到齐次变换矩阵
        pose = np.dot(T,self.endEffectorMatirx)##LI：左乘T可以得到变换后的位姿
        transE,rotE = Transform.Inverse_TransformMatrix(pose) ##LI：最终返回的是平移和欧拉角
        return transE,rotE
            
    def Jocobian(self,JointStates,FrameIndex):
        Jn = np.zeros([6,FrameIndex])
        transN,rotN = self.Pose(JointStates,FrameIndex)
        transformN = Transform.MergeTransformMatrix(transN,rotN)
        for i in range(FrameIndex):
            transI,rotI = self.Pose(JointStates,i+1)
            transformI = Transform.MergeTransformMatrix(transI,rotI)
            Jv = np.cross(transformI[0:3,2].reshape(1,3),(transformN[0:3,3].reshape(1,3) - transformI[0:3,3].reshape(1,3)))
            Jw = transformI[0:3,2].reshape(3,)
            Jn[0:3,i] = Jv
            Jn[3:6,i] = Jw
        return Jn
    
    def JocobianEnd(self,JointStates):
        Jn = np.zeros([6,self.JointStatesNum])
        transN,rotN = self.EndPose(JointStates)
        transformN = Transform.MergeTransformMatrix(transN,rotN)
        for i in range(self.JointStatesNum):
            transI,rotI = self.Pose(JointStates,i+1)
            transformI = Transform.MergeTransformMatrix(transI,rotI)
            Jv = np.cross(transformI[0:3,2].reshape(1,3),(transformN[0:3,3].reshape(1,3) - transformI[0:3,3].reshape(1,3)))
            Jw = transformI[0:3,2].reshape(3,)
            Jn[0:3,i] = Jv
            Jn[3:6,i] = Jw
        return Jn

    #Inverse solver TCPpose(x,y,z,rx,ry,rz)
    def IKSolver(self,TCPPose,efs = pow(10,-10),iter = 1000):
        joint = self.JointStates
        deltaE = 1
        count = 0

        while(deltaE > efs):
            transEnd,rotEnd = self.EndPose(joint)
            dTCPEnd = np.zero(6)
            dTCPEnd[0:3] = TCPPose[0:3] - transEnd
            dTCPEnd[3:6] = TCPPose[3:6] - rotEnd

            jocobianE = self.JocobianEnd(joint)
            djoint = np.dot(np.linalg.pinv(jocobianE),dTCPEnd)
            joint = joint + djoint
            deltaE = np.linalg.norm(djoint)
            count += 1
            if (count > iter):
                print("Solution wouldn't converge")
                return self.JointStates
        
        ikjoint = until.qq_choose(joint)
        return ikjoint
        


    def IK_PathPlan(self,JointStates,trans_goal,rot_goal):
        # print JointStates,trans_goal,rot_goal
        Path = [JointStates]
        while True: ## Planning
            if len(Path)>self.max_iter:
                print("Plan failed! Please Make sure the point can be reached.")
                return []
            Current_JointStates = Path[-1]
            trans,rot = self.Pose(Current_JointStates,len(self.Joints))
            error = (np.array(trans_goal+rot_goal)-np.array(trans+rot)).reshape(-1,1)
            # error = (np.array(trans_goal)-np.array(trans)).reshape(-1,1)
            square_error = np.dot(error.reshape(1,-1),error)[0]
            #print square_error
            if square_error[0]<self.tolerance**2:
                print("Over!")
                print(error)
                print(square_error[0])
                return Path
            
            J = self.Jocobian(Current_JointStates,len(self.Joints))
            
            ### pseudoinverse of J
            # J_pinv = np.linalg.pinv(J[0:3,:])
            J_pinv = np.linalg.pinv(J)
            # # print J_verse

            ### Proportional Control
            Kp = 20*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
            # Kp = 10*np.diag(np.array([0.01,0.01,0.01]))            
            velocity = np.dot(np.dot(J_pinv, Kp),error)

            New_JointStates = []
            for i in range(self.JointStatesNum):
                if velocity[i,0]>self.VelocityLimit:
                    v = self.VelocityLimit
                elif velocity[i,0]<-self.VelocityLimit:
                    v = -self.VelocityLimit
                else:
                    v = velocity[i,0]
                New_JointStates.append(Current_JointStates[i] + v)
            Path.append(New_JointStates)

    def RCM_Set(self,JointStates,Lambda = None):
        ## Setting RCM Goal Point
        if not Lambda is None:
            self.Lambda = Lambda
        # print self.Lambda
        trans1, rot = self.Pose(JointStates,len(self.Joints))
        trans2, rot = self.Pose(JointStates,len(self.Joints)-1)
        self.p_RCM_goal = np.array(trans2) + self.Lambda * (np.array(trans1)-np.array(trans2))
        self.p_RCM_goal = self.p_RCM_goal.reshape(3,1)

    ###用的是这个函数
    def RCM_PathPlan(self,JointStates,trans_goal):
        ## We assume that the task is about the position of the last Frame. Will Modify later

        if self.RCM_Set_Flag == 0: 
        ## Setting RCM Goal Point
            self.RCM_Set(JointStates)
            self.RCM_Set_Flag = 1

        ## init the Path and LambdaList
        Path = [JointStates]#LI:Path里面存放的是为了完成task的一系列角度值
        LambdaList = [self.Lambda]## Modify RCM point by changing initial Lambda LI：Lambda就是RCM和Pi之间的距离

        while True: ## Planning
            if len(Path)>self.max_iter:
                print("Plan failed! Please Make sure the point can be reached.")
                return []            
            ## refreshing current jointstats
            Current_JointStates = Path[-1] ##LI:Path数组中的最后一位
            self.JointStates = Current_JointStates
            Current_Lambda = LambdaList[-1]

            ## task trans and rot of the Last Frame
            transEnd,rotEnd = self.EndPose(Current_JointStates)#LI：得到最后一个坐标系的平移和欧拉角
            transN,rotN = self.Pose(Current_JointStates,self.JointStatesNum)
            ## vector form of the position of the last two frames
            p1 = np.array(transEnd).reshape(3,1) 
            p2 = np.array(transN).reshape(3,1)
            ## get Jocobian Matrix of the last two frames
            JocobianEnd = self.JocobianEnd(Current_JointStates)
            JocobianN = self.Jocobian(Current_JointStates,self.JointStatesNum)
            ## Currently don't care about orientation goal, use only the first three rows
            J1 = JocobianEnd[0:3,:]
            J2 = JocobianN[0:3,:]
            
            p_RCM = p2 + Current_Lambda * (p1 - p2) #LI: RCM的位置
            J_RCM = J2 + Current_Lambda * (J1 - J2)
            J_RCM = np.hstack((J_RCM,p1-p2))# 竖直方向上堆叠，这个才是J_RCM

            ## Extended-Task Error
            ## Calculate Error, which is a column vector, defined as the difference of goal and current state
            error = (np.array(trans_goal)-np.array(transEnd)).reshape(3,1)
            ## Calculate RCM Error
            error_RCM = self.p_RCM_goal.reshape(3,1) - p_RCM
            ## Extended-Error
            Extended_Error = np.vstack((error,error_RCM))#LI：论文中的e_t
            
            ## Use square error to judge if stop the planning, by comparing squre error and tollerance
            square_error = np.dot(Extended_Error.reshape(1,-1),Extended_Error)[0]
            #print(square_error[0])
            if square_error[0]<self.tolerance**2:
                # print("Planning Over!")
                # print("Final Error:")
                # print(Extended_Error)
                # print("Final square Error:")
                # print(square_error[0])
                # return Path, LambdaList
                return Path[-1]

            ## Error still too large, continue planning
            ## Task Jocobian
            J_t = J1
            ## Extended-Task Jocobian
            J_E = np.hstack((J_t,np.zeros((3,1))))
            ## Planning Jocobian
            J = np.vstack((J_E,J_RCM))
            ### pseudoinverse of J ##LI：J的伪逆
            J_pinv = np.linalg.pinv(J)

            ### Proportional Control
            Kp = 20*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
            
            velocity = np.dot(np.dot(J_pinv, Kp),Extended_Error)
            ### Get new JointStates
            New_JointStates = []
            for i in range(self.JointStatesNum):
                if velocity[i,0]>self.VelocityLimit:
                    v = self.VelocityLimit
                elif velocity[i,0]<-self.VelocityLimit:
                    v = -self.VelocityLimit
                else:
                    v = velocity[i,0]
                New_JointStates.append( Current_JointStates[i] + v)
            Path.append(New_JointStates)
            ## when VelocityLimit happens, Lambda refresh should be different, modify later.
            LambdaList.append(Current_Lambda+velocity[self.JointStatesNum,0])

    def RCM_PathPlan_PID(self,JointStates,trans_goal):

        ## We assume that the task is about the position of the last Frame. Will Modify later

        ## Setting RCM Goal Point
        if self.RCM_Set_Flag == 0:
            self.RCM_Set(JointStates)
            self.RCM_Set_Flag = 1

        ## init the Path and LambdaList
        Path = [JointStates]
        LambdaList = [self.Lambda]## Modify RCM point by changing initial Lambda
        ## PID Control, recording error is necessary
        ErrorSum = None
        LastError = None

        while True: ## Planning
            # sleep(0.5)
            if len(Path)>self.max_iter:
                print("Plan failed! Please Make sure the point can be reached.")
                return []            
            ## refreshing current jointstats
            Current_JointStates = Path[-1]
            Current_Lambda = LambdaList[-1]

            ## task trans and rot of the Last Frame
            trans1,rot1 = self.Pose(Current_JointStates,len(self.Joints))
            trans2,rot2 = self.Pose(Current_JointStates,len(self.Joints)-1)
            ## vector form of the position of the last two frames
            p1 = np.array(trans1).reshape(3,1)
            p2 = np.array(trans2).reshape(3,1)
            ## get Jocobian Matrix of the last two frames
            Jocobian1 = self.Jocobian(Current_JointStates,len(self.Joints))
            Jocobian2 = self.Jocobian(Current_JointStates,len(self.Joints)-1)
            ## Currently don't care about orientation goal, use only the first three rows
            J1 = Jocobian1[0:3,:]
            J2 = Jocobian2[0:3,:]
            
            p_RCM = p2 + Current_Lambda * (p1 - p2)
            J_RCM = J2 + Current_Lambda * (J1 - J2)
            J_RCM = np.hstack((J_RCM,p1-p2))

            ## Extended-Task Error
            ## Calculate Error, which is a column vector, defined as the difference of goal and current state
            error = (np.array(trans_goal)-np.array(trans1)).reshape(3,1)
            ## Calculate RCM Error
            error_RCM = self.p_RCM_goal-p_RCM
            ## Extended-Error
            Extended_Error = np.vstack((error,error_RCM))

            if ErrorSum is None:
                ErrorSum = Extended_Error
            else:
                ErrorSum = ErrorSum + Extended_Error
                ErrorSum = np.multiply(np.multiply(ErrorSum,ErrorSum<self.Inte_max),ErrorSum>-self.Inte_max) +\
                                        self.Inte_max * (ErrorSum>self.Inte_max) +\
                                        (-self.Inte_max) * (ErrorSum<-self.Inte_max)


            if LastError is None:
                LastError = Extended_Error
                ErrorDiff = 0
            else:
                ErrorDiff = Extended_Error- LastError
                LastError = Extended_Error

            # print "Error"
            # print Extended_Error
            # print "Sum"
            # print ErrorSum
            # print "Diff"
            # print ErrorDiff
            # print "  "
            
            ## Use square error to judge if stop the planning, by comparing squre error and tollerance
            square_error = np.dot(Extended_Error.reshape(1,-1),Extended_Error)[0]
            if square_error[0]<self.tolerance**2:
                print("Planning Over!")
                print("Final Error:")
                print(Extended_Error)
                print("Final square Error:")
                print(square_error[0])
                print("Iterated Times:")
                print(len(Path))
                # return Path, LambdaList
                return Path

            ## Error still too large, continue planning
            ## Task Jocobian
            J_t = J1
            ## Extended-Task Jocobian
            J_E = np.hstack((J_t,np.zeros((3,1))))
            ## Planning Jocobian
            J = np.vstack((J_E,J_RCM))
            ### pseudoinverse of J
            J_pinv = np.linalg.pinv(J)

            ### PID Control
            Kp = 10*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
            # KI = 10*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
            # KD = 10*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
            KI = 0.05*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
            KD = 11*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))            

            PID_Error = np.dot(Kp,Extended_Error) +\
                                     np.dot(KI,ErrorSum) +\
                                     np.dot(KD,ErrorDiff) 
            velocity = np.dot(J_pinv,PID_Error)
            ### Get new JointStates
            New_JointStates = []
            for i in range(self.JointStatesNum):
                if velocity[i,0]>self.VelocityLimit:
                    v = self.VelocityLimit
                elif velocity[i,0]<-self.VelocityLimit:
                    v = -self.VelocityLimit
                else:
                    v = velocity[i,0]
                New_JointStates.append( Current_JointStates[i] + v)
            Path.append(New_JointStates)
            ## when VelocityLimit happens, Lambda refresh should be different, modify later.
            LambdaList.append(Current_Lambda+velocity[self.JointStatesNum,0])

    def RCM_PathPlan_PID_Online(self,JointStates,Lambda,trans_goal):
        if self.TrackingOver == 1:
            return JointStates,Lambda
        ## We assume that the task is about the position of the last Frame. Will Modify later

        ## Setting RCM Goal Point
        if self.RCM_Set_Flag == 0:
            self.RCM_Set(JointStates)
            self.RCM_Set_Flag = 1
        
        # while True: ## Planning
        # sleep(0.5)

        ## if iterating time overflows, give up tracking
        if self.iter>self.max_iter:
            print("Track failed! Please Make sure the point can be reached.")
            ## Initiate parameters, waiting for next plan
            self.ErrorSum = None
            self.LastError = None
            self.TrackingOver = 1
            return JointStates,Lambda

        
        ## refreshing current jointstats
        Current_JointStates = JointStates
        Current_Lambda = Lambda

        ## task trans and rot of the Last Frame
        trans1,rot1 = self.Pose(Current_JointStates,len(self.Joints))
        trans2,rot2 = self.Pose(Current_JointStates,len(self.Joints)-1)
        ## vector form of the position of the last two frames
        p1 = np.array(trans1).reshape(3,1)
        p2 = np.array(trans2).reshape(3,1)
        ## get Jocobian Matrix of the last two frames
        Jocobian1 = self.Jocobian(Current_JointStates,len(self.Joints))
        Jocobian2 = self.Jocobian(Current_JointStates,len(self.Joints)-1)
        ## Currently don't care about orientation goal, use only the first three rows
        J1 = Jocobian1[0:3,:]
        J2 = Jocobian2[0:3,:]
        
        p_RCM = p2 + Current_Lambda * (p1 - p2)
        J_RCM = J2 + Current_Lambda * (J1 - J2)
        J_RCM = np.hstack((J_RCM,p1-p2))

        ## Extended-Task Error
        ## Calculate Error, which is a column vector, defined as the difference of goal and current state
        error = (np.array(trans_goal)-np.array(trans1)).reshape(3,1)
        ## Calculate RCM Error
        error_RCM = self.p_RCM_goal-p_RCM
        print(self.p_RCM_goal)
        ## Extended-Error
        Extended_Error = np.vstack((error,error_RCM))

        if self.ErrorSum is None:
            self.ErrorSum = Extended_Error
        else:
            self.ErrorSum = self.ErrorSum + Extended_Error
            self.ErrorSum = np.multiply(np.multiply(self.ErrorSum,self.ErrorSum<self.Inte_max),self.ErrorSum>-self.Inte_max) +\
                                    self.Inte_max * (self.ErrorSum>self.Inte_max) +\
                                    (-self.Inte_max) * (self.ErrorSum<-self.Inte_max)


        if self.LastError is None:
            self.LastError = Extended_Error
            ErrorDiff = Extended_Error
        else:
            ErrorDiff = Extended_Error- self.LastError
            self.LastError = Extended_Error

        print("Error")
        print(Extended_Error)
        
        ## Use square error to judge if stop the planning, by comparing squre error and tollerance
        square_error = np.dot(Extended_Error.reshape(1,-1),Extended_Error)[0]
        if square_error[0]<self.tolerance**2:
            print("Tracking Over!")
            self.TrackingOver = 1
            print("Final Error:")
            print(Extended_Error)
            print("Final square Error:")
            print(square_error[0])
            print("Iterated Times:")
            print(self.iter)
            # return Path, LambdaList
            self.ErrorSum = None
            self.LastError = None
            return JointStates,Lambda

        self.iter += 1
        ## Error still too large, continue tracking
        ## Task Jocobian
        J_t = J1
        ## Extended-Task Jocobian
        J_E = np.hstack((J_t,np.zeros((3,1))))
        ## Planning Jocobian
        J = np.vstack((J_E,J_RCM))
        ### pseudoinverse of J
        J_pinv = np.linalg.pinv(J)

        ### PID Control
        Kp = 3*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
        # KI = 10*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
        # KD = 10*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
        KI = 0*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))
        KD = 3*np.diag(np.array([0.01,0.01,0.01,0.05,0.05,0.05]))            

        PID_Error = np.dot(Kp,Extended_Error) +\
                                    np.dot(KI,self.ErrorSum) +\
                                    np.dot(KD,ErrorDiff) 
        velocity = np.dot(J_pinv,PID_Error)
        print("velocity")
        print(velocity)
        ### Get new JointStates
        New_JointStates = []
        for i in range(self.JointStatesNum):
            if velocity[i,0]>self.VelocityLimit:
                v = self.VelocityLimit
            elif velocity[i,0]<-self.VelocityLimit:
                v = -self.VelocityLimit
            else:
                v = velocity[i,0]
            New_JointStates.append( Current_JointStates[i] + v)
            New_Lambda = Current_Lambda+velocity[self.JointStatesNum,0]
        return New_JointStates, New_Lambda
        ## when VelocityLimit happens, Lambda refresh should be different, modify later.
        # LambdaList.append(Current_Lambda+velocity[self.JointStatesNum,0])

# robot = ROBOT('robot_description.urdf')
# robot.IK_PathPlan([0,0,0,0,0,0,0],[ 0.0598489865661,-0.418889224529,0.941015541553],[0,0,0])

