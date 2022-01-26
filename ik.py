#==========================通用运动学类======================#
import math
import numpy as np


class GeneralKinematic(object):
	'''
	函数依赖math和numpy
	'''
	def __init__(self, DH_0,q_min, q_max):
		self.DH_0 = DH_0
		self.theta = DH_0[:, 0]
		self.alpha = DH_0[:, 1]
		self.a = DH_0[:, 2]
		self.d = DH_0[:, 3]
		self.q_min = q_min
		self.q_max = q_max

		self.n = len(self.theta)

	#相邻关节传递矩阵
	def trans(self, theta, alpha, a, d):
		T = np.array([[math.cos(theta), -math.sin(theta) * math.cos(alpha),
					   math.sin(theta) * math.sin(alpha), a * math.cos(theta)],
					  [math.sin(theta), math.cos(theta) * math.cos(alpha),
					   -math.cos(theta) * math.sin(alpha), a * math.sin(theta)],
					  [0, math.sin(alpha), math.cos(alpha), d],
					  [0, 0, 0, 1]])
		return T

	# ZYX欧拉角转变为旋转矩阵
	def euler_zyx2rot(self,phi):
		'''
			ZYX欧拉角转变为旋转矩阵
			input:欧拉角
			output:旋转矩阵
		'''
		R = np.array([[np.cos(phi[0]) * np.cos(phi[1]),np.cos(phi[0]) * np.sin(phi[1]) * np.sin(phi[2]) - np.sin(phi[0]) * np.cos(phi[2]),
					   np.cos(phi[0]) * np.sin(phi[1]) * np.cos(phi[2]) + np.sin(phi[0]) * np.sin(phi[2])],
					  [np.sin(phi[0]) * np.cos(phi[1]),np.sin(phi[0]) * np.sin(phi[1]) * np.sin(phi[2]) + np.cos(phi[0]) * np.cos(phi[2]),
					   np.sin(phi[0]) * np.sin(phi[1]) * np.cos(phi[2]) - np.cos(phi[0]) * np.sin(phi[2])],
					  [-np.sin(phi[0]), np.cos(phi[1]) * np.sin(phi[2]), np.cos(phi[1]) * np.cos(phi[2])]])
		return R

	# 旋转矩阵转变为ZYX欧拉角
	def rot2euler_zyx(self, Re):
		'''
			ZYX欧拉角速度变为姿态角速度转化矩阵
			input:旋转矩阵
			output:欧拉角[alpha,beta,gamma]
		'''
		euler_zyx = np.zeros(3)
		if(abs(abs(Re[2, 0]) - 1) < math.pow(10, -6)):
			if(Re[2,0] < 0):
				beta = math.pi/2
				alpha = np.arctan2(-Re[1,2],Re[1,1])
				gamma = 0
			else:
				beta = -math.pi/2
				alpha = -np.arctan2(-Re[1, 2], Re[1, 1])
				gamma = 0
		else:
			p_beta = math.asin(-Re[2,0])
			cb = np.cos(p_beta)
			alpha = math.atan2(Re[1,0]*cb,Re[0,0]*cb)
			gamma = math.atan2(Re[2,1]*cb,Re[2,2]*cb)
			if((math.sin(gamma)*Re[2,1]) < 0):
				beta = math.pi - p_beta
			else:
				beta = p_beta
		euler_zyx[0] = alpha
		euler_zyx[1] = beta
		euler_zyx[2] = gamma
		for i in range(3):
			if(euler_zyx[i]>=3.14 or euler_zyx[i]<=-3.14):
				euler_zyx[i] = 0.0
		return euler_zyx

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

	#正运动学,返回齐次矩阵
	def fkine(self, qr):
		An = np.eye(4)
		for i in range(self.n):
			T = self.trans(self.theta[i] + qr[i], self.alpha[i], self.a[i], self.d[i])
			An = np.dot(An, T)  # 末端到惯性坐标系传递矩阵
		return An

	#正运动学,输出六维末端位姿,姿态用zyx欧拉角表示
	def fkine_euler(self, qr):
		xe = np.zeros(6)
		An = np.eye(4)
		for i in range(self.n):
			T = self.trans(self.theta[i] + qr[i], self.alpha[i], self.a[i], self.d[i])
			An = np.dot(An, T)  # 末端到惯性坐标系传递矩阵
		xe[0:3] = An[0:3, 3]
		xe[3:6] = self.rot2euler_zyx(An[0:3, 0:3])
		return xe
		
	#求取雅克比矩阵
	def jeco(self, qr):
		# 计算雅克比矩阵
		U = np.eye(4)
		Jn = np.zeros([6, self.n])
		T = np.zeros([4, 4, self.n])
		for i in range(self.n):
			i = self.n - i - 1
			T[:, :, i] = self.trans(self.theta[i] + qr[i], self.alpha[i], self.a[i], self.d[i])
			U = np.dot(T[:, :, i], U)
			dd = np.array([-U[0, 0] * U[1, 3] + U[1, 0] * U[0, 3],
						   -U[0, 1] * U[1, 3] + U[1, 1] * U[0, 3],
						   -U[0, 2] * U[1, 3] + U[1, 2] * U[0, 3]])
			Jn[0:3, i] = dd
			Jn[3:6, i] = U[2, 0:3]

		An = self.fkine(qr)
		R = An[0:3, 0:3]
		J_R = np.zeros([6, 6])
		J_R[0:3, 0:3] = R
		J_R[3:6, 3:6] = R

		J0 = np.dot(J_R, Jn)
		return J0

	# ***基于雅克比矩阵迭代求解逆运动学***#
	def iterate_ikine(self, q_guess, Te, efs=pow(10, -12), i_max=1000):
		'''
			本函数基于雅克比迭代求解n自由度机械臂逆运动学方程
				 q_ready是上一时刻的位置,单位:弧度;
			     T0e为DH坐标系确定的DH{0}坐标系与DH{6}之间的关系(目标矩阵);
			     efs求解误差阀值，默认值10^(-10)
				 i_limit迭代最大次数,默认值1000
			output:qq为相对与DH_q0的转动角度,单位:弧度;已处理到[-pi, pi] 之间
		'''
		# 建立初时刻迭代初值
		q_r =self.theta + q_guess

		# 计数及标签
		deltaQ = 1
		temp_count = 0

		# 迭代循环求解
		while (deltaQ > efs):

			# 求解正运动学
			An = np.eye(4)
			T = np.zeros([4, 4, self.n])

			for i in range(self.n):
				T[:, :, i] = self.trans(q_r[i], self.alpha[i], self.a[i], self.d[i])
				An = np.dot(An, T[:, :, i])

			# 计算末端误差
			dA = np.zeros(6)
			dA[0:3] = Te[0:3, 3] - An[0:3, 3]
			dA[3:6] = 0.5 * (np.cross(An[0:3, 0], Te[0:3, 0]) + np.cross(An[0:3, 1], Te[0:3, 1])
							 + np.cross(An[0:3, 2], Te[0:3, 2]))

			# 计算雅克比矩阵
			U = np.eye(4)
			Jn = np.zeros([6, self.n])
			for i in range(self.n):
				i = self.n - i - 1
				U = np.dot(T[:, :, i], U)

				dd = np.array([-U[0, 0] * U[1, 3] + U[1, 0] * U[0, 3],
							   -U[0, 1] * U[1, 3] + U[1, 1] * U[0, 3],
							   -U[0, 2] * U[1, 3] + U[1, 2] * U[0, 3]])
				Jn[0:3, i] = dd
				Jn[3:6, i] = U[2, 0:3]

			R = An[0:3, 0:3]
			J_R = np.zeros([6, 6])
			J_R[0:3, 0:3] = R
			J_R[3:6, 3:6] = R

			J0 = np.dot(J_R, Jn)
			# 求取关节角关节角度偏差值
			dq = np.dot(np.linalg.pinv(J0), dA)
			q_r = q_r + dq
			deltaQ = np.linalg.norm(dq)
			temp_count = temp_count + 1
			if (temp_count > i_max):
				print("Solution wouldn't converge")
				return q_guess

		q_tmp = q_r - self.theta
		q = self.qq_choose(q_tmp)
		return q

	#机械臂关节极限判断，返回值为0或1
	def exceed_joint_limit(self, qq, q_min, q_max):
		'''
			判断关节角是否超出限制
			input:关节角，关节角范围
			outpu：0,未超出，1超出
		'''
		n = len(qq)
		limit = False
		for i in range(n):
			if((qq[i] < q_min[i]) or (qq[i] > q_max[i])):
				print ("第", i+1, "关节超出极限:", qq[i]*180/np.pi)
				limit = True
				break
		return limit
	#带关节限制
	def iterate_ikine_limit_xyz(self, q_guess, Xe):
		Te = np.eye(4)
		Te[0:3, 0:3] = self.euler_zyx2rot(Xe[3:])
		Te[0:3, 3] = Xe[:3]
		print("Te:", Te)
		qr = self.iterate_ikine(q_guess, Te)
		flag = self.exceed_joint_limit(qr ,self.q_min, self.q_max)
		if(flag):
			#print "flag:", flag
			qr = np.copy(q_guess)
		return qr

	# 带关节限制
	def iterate_ikine_limit(self, q_guess, Te):
		qr = self.iterate_ikine(q_guess, Te)
		flag = self.exceed_joint_limit(qr, self.q_min, self.q_max)
		if (flag):
			# print "flag:", flag
			qr = np.copy(q_guess)
		return qr
