import math
import numpy as np


# 将关节角计算到正负pi
def qq_choose(self, qq):
    '''
        本函数用于选着关节角范围
        input:qq为计算出的关节角
        output:q关节角范围[-pi,pi]
    '''
    q = np.copy(qq)
    for i in range(self.n):
        while (q[i] > math.pi):
            q[i] = q[i] - 2 * math.pi
        while (q[i] < - math.pi):
            q[i] = q[i] + 2 * math.pi
    return q

def trajPlanning(time,k,initPoint):
    '''
        本函数用于阿基米德螺旋线生成
        input：
            time:规划时间点
            k:1X3数组 旋转参数 半径参数 高度参数
            initPoint:螺旋线起始点
        output：1X3数组 三维坐标
    '''
    theta = k[0] * time
    r  = k[1] * time
    high = k[2] * time

    x = r * math.cos(theta)
    y = r * math.sin(theta)
    z = -high 
    
    goal = np.zeros(3)
    goal[0] = initPoint[0] + x
    goal[1] = initPoint[1] + y
    goal[2] = initPoint[2] + z
    return goal