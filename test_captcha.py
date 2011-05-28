from utils import get_points_from_img,get_elements
from SC import SC
import cv
from numpy import *
import sys

if __name__ == '__main__':
    import sys
    def make_graph(P1,P2,COST,LINES):
        from matplotlib import pylab
               
        x = []
        y = []
        al = P1+P2
        for i in xrange(len(al)):
            x.append(al[i][0])
            y.append(al[i][1])
        
        ax = pylab.subplot(111)
        pylab.grid(True)
        
        pylab.plot(P1[0],P1[1],'go',P2[0],P2[1],'ro')
        
        ax.set_title('Total cost: %s' % COST)       
        
        for l in LINES:     
            pylab.plot((l[0][0],l[1][0]),(l[0][1],l[1][1]), 'k-')
        
        pylab.show()

    a = SC(r_outer=2)
    sampls = 30

    imgs = get_elements('test.png')

    points1,t1 = get_points_from_img('9M.png',simpleto=sampls)
    P = a.compute(points1)
    x1 = [p[0] for p in points1]
    y1 = [400-p[1] for p in points1]

    for i,img in enumerate(imgs):
        cv.ShowImage('Preview - %s' % i,img)
        cv.WaitKey() 
    sys.exit()
    
    for img in imgs:
    
        cv.ShowImage('Preview2',img)
        cv.WaitKey()   
    
        points2,t2 = get_points_from_img(img,simpleto=sampls)
    
        if points2:

            Q = a.compute(points2)
            x2 = [p[0] for p in points2]
            y2 = [400-p[1] for p in points2]
                
            COST,indexes = a.diff(Q,P)
            
            # getting correspoding points arrays for interpolation
            pp = []
            qp = []
            for i,k in indexes:
                qp.append(points2[i])
                pp.append(points1[k])
                

            fx,fy,diff,affcost = a.interpolate(qp,pp)
            x3 = [fx(p[0]) for p in points2]
            y3 = [400-fy(p[1]) for p in points2] 
            LINES = []
            for q1,p1 in indexes:
                LINES.append([[points1[p1][0],400-points1[p1][1]],[points2[q1][0],400-points2[q1][1]]])
            
            polarity_flag = 1
            ori_weight = 0.1

            costmat_shape = a.cost(Q,P)
            theta_diff = kron(ones((1,sampls)),t1) - kron(ones((sampls,1)),t2.H)
            if polarity_flag:
                # use edge polarity
                costmat_theta=0.5*(1-cos(theta_diff))
            else:
                # ignore edge polarity
                costmat_theta=0.5*(1-cos(2*theta_diff))
                
                
            print costmat_shape.shape,costmat_theta.shape
            costmat=(1-ori_weight)*costmat_shape+ori_weight*costmat_theta;
            
            a1=costmat.min(0)
            a2=costmat.min(1)
            sc_cost=max(mean(a1),mean(a2));

            print "Shape cost: %s\nBending energy: %s\nAffine Cost: %s\n" % (sc_cost,diff,affcost)

            TOTAL = 0.1*diff+sc_cost+0.3*affcost
            print 'TOTAL MATCH:',TOTAL

            make_graph((x1,y1),(x2,y2),TOTAL,LINES)  

            #make_graph((x1,y1),(x2,y2),COST,LINES)  
            #make_graph((x1,y1),(x3,y3),diff)  
                    
            
    
