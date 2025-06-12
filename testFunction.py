
import gdsfactory as gf
import cmath
import math
from kfactory.kcell import show
import os
from matplotlib import pyplot as plt
os.environ['GIT_PYTHON_REFRESH'] = 'quiet'

def ArraySAS(focal, Ldev, Lctrl, dL, Ltap, Ls, slab_angle, dch, da, Nch, Na, SW=0):
    print("Str-Arc-Str Array Scheme")   
    '''
    focal : FPR Length
    Ldev : horizontal distance of full device [ from center point of input FPR to center point of output fpr)
    Lctrl : initial value to optimize path parameters to make fixed path difference.
    dL : target path difference
    Ltap : taper length
    Ls : straight length
    slab_algle : FPR angle in degrees
    dch : radial spacing of input, output waveguids
    da : radial spacing of array waveguides at slab junction
    Nch : output channel count
    Na : array count
    SW : it uses root finding method, such as bisection method to make optimal layout [ 1 : optimize on, but not included , 0 : off ]
    '''
    
    
    c = gf.Component()
        
    portArAngles = []
    refAngle = math.radians(slab_angle)
    csSM = gf.cross_section.strip(width=4.5)
    csTAP = gf.cross_section.strip(width=6.0)
    
    dAr = da / focal
    Lref = focal * 2 + (Ltap + Ls) * 2 
    res = Lctrl
    paths = gf.Path()

    Eki = 0 # straight-to-arc offset coefficient

    if SW ==1:
        print("optimize")
        res = bisec(funcSAS,Lctrl-50000,Lctrl+50000,*[Lctrl, focal, Ldev, Ltap, Ls, slab_angle, da, Na])
        print(f"Root: {res}")
        print("optimize end")
    
    L0 = Lref + res    
    for i in range(0,Na):
        portArAngles.append(refAngle - dAr * (Na-1) * 0.5 + dAr * i)

    
    I_ain = []
    S_ain = []
    R_ain = []
    H = []
    for i in range(0,Na):
        TA = portArAngles[i]
        #print(TA)
        I_ain.append(L0 + dL * i)
        S_ain.append(0.5 * (I_ain[i] - (Ldev * TA / math.sin(TA))) / (1 - TA * math.cos(TA) / math.sin(TA)))
        R_ain.append((0.5 * Ldev - (S_ain[i] * math.cos(TA))) / math.sin(TA))
        # Left Straight offset points region #
        
        p1 = [cmath.rect(focal,TA).real,cmath.rect(focal,TA).imag]
        p2 = [p1[0] + cmath.rect(Ltap,TA).real, p1[1] + cmath.rect(Ltap,TA).imag]
        p3 = [p2[0] + cmath.rect(Ls,TA).real, p2[1] + cmath.rect(Ls,TA).imag]
        p4 = [p3[0] + cmath.rect(S_ain[i] - focal -Ltap- Ls ,TA).real, p3[1] + cmath.rect(S_ain[i] - focal -Ltap- Ls ,TA).imag]
        if i==0 or i==(Na-1):
            print(p1,p2,p3,p4)

        offsetR = Eki / R_ain[i]
        
        # Left Taper draw #
        pathTap = gf.path.transition(cross_section1 = csTAP, cross_section2 = csSM, width_type="linear")
        taper = gf.path.straight(length=Ltap,npoints=2)
        taper.drotate(math.degrees(TA))
        taper.dmove(p1)
        dTaper = gf.path.extrude_transition(taper, pathTap)

        pathStr = gf.path.straight(length=S_ain[i]-focal-Ltap)
        pathStr.drotate(math.degrees(TA))
        pathStr.dmove(p2)
        dStr = gf.path.extrude(pathStr,csSM)
        
        narcpts = int(math.degrees(TA)/0.05)        
        pathArc = gf.path.arc(radius=R_ain[i]-(Eki/R_ain[i]),angle = -math.degrees(TA), npoints = narcpts)
        pathArc.drotate(math.degrees(TA))
        pathArc.dmove(p4)
        pathArc.dmovex(Eki/R_ain[i] * math.sin(TA))
        pathArc.dmovey(-Eki/R_ain[i] * math.cos(TA))
        dArc = gf.path.extrude(pathArc,csSM)
        
        c<<dTaper
        c<<dStr
        c<<dArc
        
        H.append(S_ain[i]*math.sin(TA)+R_ain[i]*(1-math.cos(TA)))
    

    #plt.plot(R_ain)
    
    return c

    
def Slab(focal, slab_angle = 64.5, dch = 25.5, da= 8.5, Nch = 96, Na = 528, dAng = 0.01,ratio_in=3.0, ratio_out=1.5, preFix ="N"):
    
    c = gf.Component()
    
    print(f"Slab arc resolution : {dAng}")

    dPhi = float( dch / focal )
    dAr = float( da / focal )
    angleIn = dPhi * ( Nch - 1 ) * ratio_in
    angleOut = dAr * ( Na - 1 ) * ratio_out
    point_num_in = int( angleIn / math.radians(dAng))
    point_num_out = int( angleOut / math.radians(dAng))
    dangle_in = angleIn / point_num_in
    dangle_out = angleOut / point_num_out
    refAngle =  math.radians(slab_angle)

    angleList = []

    portArAngles = []
    portArPoints = []
    
    for i in range(0,Na):
        portArAngles.append(math.degrees(refAngle - dAr * (Na-1) * 0.5 + dAr * i))
        portArPoints.append([cmath.rect(focal, refAngle - dAr * (Na-1) * 0.5 + dAr * i ).real,
                      cmath.rect(focal, refAngle - dAr * (Na-1) * 0.5 + dAr * i ).imag])
                            

    points = []

    for i in range(0,point_num_in):
        points.append( [cmath.rect(focal, refAngle - angleIn * 0.5 + dangle_in * i ).real,
                      cmath.rect(focal, refAngle - angleIn * 0.5 + dangle_in * i ).imag] 
                     )

    mid_point = [cmath.rect(0.5 * focal, math.radians(slab_angle)).real,
                 cmath.rect(0.5 * focal, math.radians(slab_angle)).imag]

    print(mid_point)
                 
    for i in range(0,point_num_out):
        points.append( [mid_point[0] + cmath.rect(0.5 * focal, math.pi + refAngle - angleOut * 0.5 + dangle_out * i ).real,
                      mid_point[1] + cmath.rect(0.5 * focal, math.pi + refAngle - angleOut * 0.5 + dangle_out * i ).imag] 
                     )
    #print(points)
    points.append([cmath.rect(focal, refAngle - angleIn * 0.5 ).real,
                      cmath.rect(focal, refAngle - angleIn * 0.5 ).imag])
    c.add_polygon(points,layer = LayerSlab)
    for i in range(0,Na):
        c.add_port(name=f"AR {preFix} {i+1}",
                   width = da * 0.8,
                   orientation = portArAngles[i],
                   center=(portArPoints[i][0],portArPoints[i][1]),layer = LayerPortNum)
        
    c.draw_ports()
    c.pprint_ports()
    return c


LayerSlab = [1,0]
LayerArray = [2,0]
LayerInTap = [3,0]
LayerOutTap = [4,0]
LayerSeg = [5,0]
LayerArTap = [6,0]
LayerPortNum = [999,0]

c = gf.Component()
refPnt = [0,0]

SlabIn = c<<Slab(26751.49415357,64.5,preFix="IN")
SlabOut = c<<Slab(26751.49415357,180-64.5,preFix="OUT")
SlabOut.movex(36500)

c2 = ArraySAS(26751.49415357,36500,10257.764912353196,30.842574502627887,500,500,64.5,25.5,8.5,96,528,SW=0)
c<<c2
c.show()

