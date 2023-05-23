
from wsgiref.simple_server import make_server
import falcon
import json
#from robots import Robot
import ROXS
from ROXS import Robot


'''
API

#  Control of individual robots:
=================================
# /Home
R1 is the target robot object instance 
# Web API:    /home robotId

curl -H "Content-Type: application/json" --request POST  localhost:8000/Home -d "{\"id\": \"R1\"}"

R1.Connect_and_Home()

# /LM and /JM
# x, y, z are the target coordinates (floats)
# WebAPI:  /?M {'id': robotId,  'x': x, 'y': y, 'z': z} 


curl -H "Content-Type: application/json" --request POST  localhost:8000/LM -d "{\"id\": \"R1\", \"x\": 0, \"y\": 0,  \"z\": 0}"
curl -H "Content-Type: application/json" --request POST  localhost:8000/JM -d "{\"id\": \"R2\", \"x\": 0, \"y\": 0,  \"z\": 0}"


R1.SafeLM(R1.TPP, (x, y, z))
R1.SafeJM(R1.TPP, (x, y, z))

# WebAPI: /SLAP {'id':robotId, 'PR_POS' : pos, 'loops' : loops, 'tests': tests}

curl -H "Content-Type: application/json" --request POST  localhost:8000/SLAP -d "{\"id\": \"R1\", \"PR_pos\": 13.1, \"loops\": 3,  \"tests\": 2}"

R1.SLAP(PR_pos, loops, tests)

# Run SLAP on all robots  
# WebAPI:  /SLAP_ALL, {'PR_POS' : pos, 'loops' : loops, 'tests': tests}

curl -H "Content-Type: application/json" --request POST  localhost:8000/SLAP_ALL -d "{\"PR_pos\": 9.4, \"loops\": 2,  \"tests\": 4}"

robots.SLAP_ALL(pos, loops, tests)


# Set the position of the photo resistor   
# Parameter PR_POS float
# WebAPI:  /SetTarget {'PR_POS': pos}

curl -H "Content-Type: application/json" --request POST  localhost:8000/SetTarget -d "{\"PR_pos\": 7.7}"

robots.SetTargets(PR_pos)



#  Multi Robot Commands:
# 
# Parameters?

# ei:  float
# ef:  float
# points:  integer


curl -H "Content-Type: application/json" --request POST  localhost:8000/RunRowland -d "{\"ei\": 7.7, \"ef\": 14.9,  \"points\": 21}"

#WebAPI  /RunRowland {'ei' : ei, 'ef' : ef, 'points' : points} 

# Function to convert energy to position 
zi=e2z(ei)
zf=e2z(ef)

robots.RunRowland(zi, zf, points)

'''


class robotCommand:

    def __init__(self, commandString, *args):
        self.args = args
        self.commandString = commandString

    # def on_get(self, req, resp):
        # Build a json string from robot instance data with required status. Example uses static string
        # resp.body = json.dumps(....)

    def on_post(self, req, resp):

        kwargs = req.media

        id = kwargs.pop("id")

        missing = set(self.args) - set(kwargs.keys())
        extra = set(kwargs.keys()) - set(self.args)

        if not id in Robot.robots:
            resp.status = falcon.HTTP_404
            resp.text = json.dumps({"success": False})
            print(f"Invalid robot id {id}")

        elif missing != set():
            txt = ','.join(missing)
            resp.status = falcon.HTTP_400
            resp.text = json.dumps({"success": False})
            print(f"Missing parameters: {txt}")

        elif extra != set():
            txt = ','.join(extra)
            resp.status = falcon.HTTP_400
            resp.text = json.dumps({"success": False})
            print(f"Extra parameters: {txt}")
        else:

            print(f"{self.commandString+'API'}, Robot {id}, {kwargs}")
            resp.status = falcon.HTTP_200
            # Call function to trigger requested action here. (Probably should run in a separate thread)

            self.r = Robot.robots[id]

            func = getattr(self.r, self.commandString + "API")

            func(**kwargs)

            #resp.body = f"Robot {id}: {command}"
            resp.text = json.dumps({"success": True})


class HomeCommand(robotCommand):

    '''
    Add this function to robot class
        
    def HomeAPI(self, **kwargs):
        self.Connect_and_Home()
    
    '''

    def __init__(self, app):
        args = []
        super().__init__('Home', *args)
        app.add_route('/Home', self)


class JMCommand(robotCommand):

    '''
    Add this function to robot class
        
    def JMAPI(self, **kwargs):
        self.SafeJM(self.TPP,(float(kwargs['x']), float(kwargs['y']), float(kwargs['z'])))
    
    '''

    def __init__(self, app):
        args = ['x', 'y', 'z']
        super().__init__('JM', *args)
        app.add_route('/JM', self)


class LMCommand(robotCommand):
    '''
    Add this function to robot class
        
    def JMAPI(self, **kwargs):
        self.SafeLM(self.TPP, (float(kwargs['x']), float(kwargs['y']), float(kwargs['z'])))
    
    '''

    def __init__(self, app):
        args = ['x', 'y', 'z']
        super().__init__('LM', *args)
        app.add_route('/LM', self)


class SLAPCommand(robotCommand):
    '''
    Add this function to robot class
        
    def SLAPAPI(self, **kwargs):
        self.SLAP(self.TPP, float(kwargs['PR_pos']), int(kwargs['loops']), int(kwargs['tests'] ))   
    '''

    def __init__(self, app):
        args = ['PR_pos', 'loops', 'tests']
        super().__init__('SLAP', *args)
        app.add_route('/SLAP', self)


class systemCommand:

    def __init__(self, commandString, *args):
        self.args = args
        self.commandString = commandString

    # def on_get(self, req, resp):
        # Build a json string from system status
        #resp.body = json.dumps()

    def on_post(self, req, resp):

        kwargs = req.media

        missing = set(self.args) - set(kwargs.keys())
        extra = set(kwargs.keys()) - set(self.args)

        if missing != set(): 
            txt = ','.join(missing)
            resp.status = falcon.HTTP_400
            resp.text = json.dumps({"success": False})
            print(f"Missing parameters: {txt}")

        elif extra != set():
            txt = ','.join(extra)
            resp.status = falcon.HTTP_400
            resp.text = json.dumps({"success": False})
            print(f"Extra parameters: {txt}")
        else:

            print(f"{self.commandString+'API'}, Robot {id}, {kwargs}")
            resp.status = falcon.HTTP_200
            # Call function to trigger requested action here. (Probably should run in a separate thread)

            func = getattr(Robot, self.commandString + "_API")

            func(**kwargs)

            resp.text = json.dumps({"success": True})


class SetTargetCommand(systemCommand):
    '''
    Add this classmethod Robot class
        
    def Set_Target_API(cls, **kwargs):
        Robot.SetTargets(float(kwargs['PR_pos']))   
    '''

    def __init__(self, app):
        args = ['PR_pos']
        super().__init__('SetTarget', *args)
        app.add_route('/SetTarget', self)
        
class SetRowPos_mmCommand(systemCommand):
    '''
    Add this classmethod Robot class
        
    def Set_Target_API(cls, **kwargs):
        Robot.SetTargets(float(kwargs['PR_pos']))   
    '''

    def __init__(self, app):
        args = ['pos']
        super().__init__('SetRowPos_mm', *args)
        app.add_route('/SetRowPos_mm', self)
        
class RSetRowPos_Command(systemCommand):
    '''
    Add this classmethod Robot class
        
    def Set_Target_API(cls, **kwargs):
        Robot.SetTargets(float(kwargs['PR_pos']))   
    '''

    def __init__(self, app):
        args = ['hv']
        super().__init__('RSetRowPos', *args)
        app.add_route('/RSetRowPos', self)
        
class RSetRowPos_mm_Command(systemCommand):
    '''
    Add this classmethod Robot class
        
    def Set_Target_API(cls, **kwargs):
        Robot.SetTargets(float(kwargs['PR_pos']))   
    '''

    def __init__(self, app):
        args = ['pos']
        super().__init__('RSetRowPos_mm', *args)
        app.add_route('/RSetRowPos_mm', self)
        
class SetRowPosCommand(systemCommand):
    '''
    Add this classmethod Robot class
        
    def Set_Target_API(cls, **kwargs):
        Robot.SetTargets(float(kwargs['PR_pos']))   
    '''

    def __init__(self, app):
        args = ['hv']
        super().__init__('SetRowPos', *args)
        app.add_route('/SetRowPos', self)
        
        
        


        


class SLAP_ALLCommand(systemCommand):
    '''
    Add this classmethod Robot class
        
    def SLAP_ALL_API(cls, **kwargs):
        for r in cls.robots.values():
             r.SLAP(r.TPP, float(kwargs['PR_pos']), int(kwargs['loops']), int(kwargs['tests'] ))   
    '''

    def __init__(self, app):
        args = ['PR_pos', 'loops', 'tests']
        super().__init__('SLAP_ALL', *args)
        app.add_route('/SLAP_ALL', self)
        
class SLAP_ALL_I_mmCommand(systemCommand):
    '''
    Add this classmethod Robot class
        
    def SLAP_ALL_API(cls, **kwargs):
        for r in cls.robots.values():
             r.SLAP(r.TPP, float(kwargs['PR_pos']), int(kwargs['loops']), int(kwargs['tests'] ))   
    '''

    def __init__(self, app):
        args = ['zi','zf', 'loops', 'tests']
        super().__init__('SLAP_ALL_I_mm', *args)
        app.add_route('/SLAP_ALL_I_mm', self)
        
class SLAP_ALL_ICommand(systemCommand):
    '''
    Add this classmethod Robot class
        
    def SLAP_ALL_API(cls, **kwargs):
        for r in cls.robots.values():
             r.SLAP(r.TPP, float(kwargs['PR_pos']), int(kwargs['loops']), int(kwargs['tests'] ))   
    '''

    def __init__(self, app):
        args = ['ei','ef', 'loops', 'tests']
        super().__init__('SLAP_ALL_I', *args)
        app.add_route('/SLAP_ALL_I', self)

class setdCommand(systemCommand):
    '''
    Add these classmethods to the Robot class
        
    @classmethod
    def RunRowland_API(cls, **kwargs):
        ei = float(kwargs['ei'])
        ef = float(kwargs['ef'])

        zi = e2f(ei)
        zf = e2f(ef)
        cls.RunRowland(zi, zf, int(kwargs['points']))

    @classmethod
    def RunRowland(cls, zi, zf, points):
        """
        Move your runrowland code here:
        """

        pass
    '''

    def __init__(self, app):
        args = ['material','h','k','l']
        super().__init__('setd', *args)
        app.add_route('/setd', self)
        
class RunRowlandCommand(systemCommand):
    '''
    Add these classmethods to the Robot class
        
    @classmethod
    def RunRowland_API(cls, **kwargs):
        ei = float(kwargs['ei'])
        ef = float(kwargs['ef'])

        zi = e2f(ei)
        zf = e2f(ef)
        cls.RunRowland(zi, zf, int(kwargs['points']))

    @classmethod
    def RunRowland(cls, zi, zf, points):
        """
        Move your runrowland code here:
        """

        pass
    '''

    def __init__(self, app):
        args = ['ei', 'ef', 'points']
        super().__init__('RunRowland', *args)
        app.add_route('/RunRowland', self)
        
class RunRowland_mmCommand(systemCommand):
    '''
    Add these classmethods to the Robot class
        
    @classmethod
    def RunRowland_API(cls, **kwargs):
        ei = float(kwargs['ei'])
        ef = float(kwargs['ef'])

        zi = e2f(ei)
        zf = e2f(ef)
        cls.RunRowland(zi, zf, int(kwargs['points']))

    @classmethod
    def RunRowland(cls, zi, zf, points):
        """
        Move your runrowland code here:
        """

        pass
    '''

    def __init__(self, app):
        args = ['zi', 'zf', 'points']
        super().__init__('RunRowland_mm', *args)
        app.add_route('/RunRowland_mm', self)


if __name__ == '__main__':

    app = falcon.App()

    h = HomeCommand(app)
    jm = JMCommand(app)
    lm = LMCommand(app)
    slap = SLAPCommand(app)
    stgt = SetTargetCommand(app)
    slapall = SLAP_ALLCommand(app)
    slapalli = SLAP_ALL_ICommand(app)
    slapallimm= SLAP_ALL_I_mmCommand(app)
    rwlnd = RunRowlandCommand(app)
    rwlndmm = RunRowland_mmCommand(app)
    setrowposmm = SetRowPos_mmCommand(app)
    Rsetrowposmm = RSetRowPos_mm_Command(app)
    Rsetrowpos = RSetRowPos_Command(app)
    setrowpos=SetRowPosCommand(app)
    setd=setdCommand(app)
    

    with make_server('', 8000, app) as httpd:
        print('Serving on port 8000...')

        # Serve until process is killed
        httpd.serve_forever()
