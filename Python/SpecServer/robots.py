def e2z(e):
    '''
    Replace with code to convert energy to int z position
    '''

    return int(e * 11.3 + 2)


class Robot:
    # Add this to Robot class

    robots = {}

    @classmethod
    def SLAP_ALL_API(cls, **kwargs):
        for r in cls.robots.values():
            r.SLAP(r.TPP, float(kwargs['PR_pos']), int(kwargs['loops']), int(kwargs['tests']))

    @classmethod
    def RunRowland_API(cls, **kwargs):
        ei = float(kwargs['ei'])
        ef = float(kwargs['ef'])

        zi = e2z(ei)
        zf = e2z(ef)
        cls.RunRowland(zi, zf, int(kwargs['points']))

    @classmethod
    def RunRowland(cls, zi, zf, points):
        '''
        Move the runrowland code here:
        '''
        print(f"RunRowland {zi} {zf} {points}")

    @classmethod
    def SetTarget_API(cls, **kwargs):
        pr_pos = kwargs['PR_pos']
        cls.SetTarget(pr_pos)

    @classmethod
    def SetTarget(cls, pr_pos):
        '''
        Move the set target code here
        '''

        print(f"SetTarget {pr_pos:.3f}")

    # def __init__(self, webName, Name, IP, Tool, TPP, APP, LEDT, Home, ESP):

    def __init__(self, webName, IP):
        self.webName = webName
        #self.Name = RDK.Item(Name)
        self.IP = IP
        #self.Tool = RDK.Item(Tool)
        self.TPP = webName
        #self.TPP = RDK.Item(TPP)
        #self.APP = RDK.Item(APP)
        #self.LEDT = RDK.Item(LEDT)
        #self.Home = RDK.Item(Home)
        #Error = np.array(Pose_2_KUKA(self.TPP.Pose())) - np.array(Pose_2_KUKA(self.APP.Pose()))
        #self.Error = Error
        #self.Ex = Error[0]
        #self.Ey = Error[1]
        #self.Ez = Error[2]
        #self.ESP = ESP
        #self.JPos = self.Name.Joints()

        Robot.robots[webName] = self    # add to dictionary of all robots

    def SafeJM(self, target, trans_coord):
        # SAFE collsion free joint movement based on robodk station
        # Input robodk item as target and coordintate offset

        print(f"SafeJM {self.webName}, {target}, {trans_coord}")

        #Final_Pose = target.Pose() * transl(trans_coord)
        #self.Name.setSpeed(-1, 50)
        # if self.Name.MoveL_Test(self.JPos, Final_Pose) == 0:         # linear move to the approach position

            # self.Name.MoveL(Final_Pose)
            #self.JPos = self.Name.Joints()
        # else:
            #print("collision avoided at " + str(trans_coord))

    def SafeLM(self, target, trans_coord):
            # SAFE collsion free linear movement based on robodk station
        print(f"SafeLM {self.webName}, {target}, {trans_coord}")

            # Input robodk item as target and coordintate offset
            #Final_Pose = target.Pose() * transl(trans_coord)
            #self.Name.setSpeed(-1, 50)
            # if self.Name.MoveL_Test(self.JPos, Final_Pose) == 0:         # linear move to the approach position

                # self.Name.MoveL(Final_Pose)
                #self.JPos = self.Name.Joints()
            # else:
                #print("collision avoided at " + str(trans_coord))
    def SLAP(self, target, pr_pos, loops, tests):
        print(f"SLAP {self.webName}, {target}, pr_pos {pr_pos}, loops {loops}, tests {tests}")

    def HomeAPI(self, **kwargs):
        print(f"{self.webName} HomeAPI {kwargs}")
        # self.Connect_and_Home()

    def JMAPI(self, **kwargs):
        self.SafeJM(self.TPP, (float(kwargs['x']), float(kwargs['y']), float(kwargs['z'])))

    def LMAPI(self, **kwargs):
        self.SafeLM(self.TPP, (float(kwargs['x']), float(kwargs['y']), float(kwargs['z'])))

    def SLAPAPI(self, **kwargs):
        self.SLAP(self.TPP, float(kwargs['PR_pos']), int(kwargs['loops']), int(kwargs['tests']))


R1 = Robot('R1', '192.168.0.100')
R2 = Robot('R2', '192.168.0.101')
