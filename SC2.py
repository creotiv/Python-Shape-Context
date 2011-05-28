from numpy import *
import math
# Hungurian algorithm implementation
import munkres
from utils import get_points_from_img,get_elements
import time
import heapq
import cv

def logspace(d1, d2, n):
    sp =  [( 10 **(d1 + k * (d2-d1)/(n - 1)))   for k in xrange(0, n -1)]
    sp.append(10 ** d2)
    return sp
    
def euclid_distance(p1,p2):
    return math.sqrt( ( p2[0] - p1[0] ) ** 2 + ( p2[1] - p1[1] ) ** 2 )
    
    
def get_angle(p1,p2):
    """Return angle in radians"""
    return math.atan2((p2[1] - p1[1]),(p2[0] - p1[0]))
    
    
class SC(object):

    HUNGURIAN = 1

    def __init__(self,nbins_r=5,nbins_theta=12,r_inner=0.1250,r_outer=2.0):
        self.nbins_r        = nbins_r
        self.nbins_theta    = nbins_theta
        self.r_inner        = r_inner
        self.r_outer        = r_outer
        self.nbins          = nbins_theta*nbins_r


    def _dist2(self, x, c):
        result = zeros((len(x), len(c)))
        for i in xrange(len(x)):
            for j in xrange(len(c)):
                result[i,j] = euclid_distance(x[i],c[j])
        return result
        
        
    def _get_angles(self, x):
        result = zeros((len(x), len(x)))
        for i in xrange(len(x)):
            for j in xrange(len(x)):
                result[i,j] = get_angle(x[i],x[j])
        return result
        
    
    def get_mean(self,matrix):
        """ This is not working. Should delete this and make something better"""
        h,w = matrix.shape
        mean_vector = matrix.mean(1)
        mean = mean_vector.mean()
        
        return mean

        
    def compute(self,points,r=None):
        
        t = time.time()
        r_array = self._dist2(points,points)
        mean_dist = r_array.mean()
        r_array_n = r_array / mean_dist
        
        r_bin_edges = logspace(log10(self.r_inner),log10(self.r_outer),self.nbins_r)  

        r_array_q = zeros((len(points),len(points)), dtype=int)
        for m in xrange(self.nbins_r):
           r_array_q +=  (r_array_n < r_bin_edges[m])

        fz = r_array_q > 0
        
        theta_array = self._get_angles(points)
        # 2Pi shifted
        theta_array_2 = theta_array + 2*math.pi * (theta_array < 0)
        #theta_array_q = 1 + floor(theta_array_2 /(2 * math.pi / self.nbins_theta))
        # norming by mass(mean) angle v.0.1 ############################################
        # By Andrey Nikishaev
        theta_array_delta = theta_array - theta_array.mean()
        theta_array_delta_2 = theta_array_delta + 2*math.pi * (theta_array_delta < 0)
        theta_array_q = 1 + floor(theta_array_delta_2 /(2 * math.pi / self.nbins_theta))
        ################################################################################

        
        BH = zeros((len(points),self.nbins))
        for i in xrange(len(points)):
            sn = zeros((self.nbins_r, self.nbins_theta))
            for j in xrange(len(points)):
                if (fz[i, j]):
                    sn[r_array_q[i, j] - 1, theta_array_q[i, j] - 1] += 1
            BH[i] = sn.reshape(self.nbins)
            
        print 'PROFILE: ' + str(time.time()-t)     
            
        return BH,theta_array_2        
        
        
    def _cost(self,hi,hj):
        cost = 0
        for k in xrange(self.nbins):
            if (hi[k] + hj[k]):
                cost += ( (hi[k] - hj[k])**2 ) / ( hi[k] + hj[k] )
            
        return cost*0.5
        
    
    def cost(self,P,Q):
        p,_ = P.shape
        p2,_ = Q.shape
        C = zeros((p,p2))
        for i in xrange(p):
            for j in xrange(p2):
                C[i,j] = self._cost(Q[j]/p,P[i]/p2)    
        
        return C
        
    def __hungurian_method(self,C):
        t = time.time()
        m = munkres.Munkres()
        indexes = m.compute(C.tolist())
        total = 0
        for row, column in indexes:
            value = C[row][column]
            total += value
        print 'PROFILE2: ' + str(time.time()-t)     

        return total,indexes

    def quick_diff(self,P,Qs,method=HUNGURIAN):
        res = []
        
        p,_ = P.shape
        q,_ = Qs.shape
        for i in xrange(p):
            for j in xrange(q):
                heapq.heappush(res,(self._cost(P[i],Qs[j]),i) )
        
        data = zeros((q,self.nbins))
        for i in xrange(q):
            data[i] = P[heapq.heappop(res)[1]]
       
        return self.diff(data,Qs)

    
    def cost_angles(self,Pa,Qa):
        pass

    def diff(self,P,Q,beta=0.1,method=HUNGURIAN):
        result = None
        C = self.cost(P[0],Q[0])*(1-beta) + beta*self.cost_angles(P[1],Q[1])

        if method == self.HUNGURIAN:
            result = self.__hungurian_method(C)
        else:
            raise Exception('No such optimization method.')
            
        return result
            

    def get_contextes(self,BH,r=5):
        res = zeros((r,self.nbins))
        used = []
        sums = []
        
        # get r shape contexts with maximum number of connected elements
        # this gives same result for same query
        for i in xrange(len(BH)):
            heapq.heappush(sums,(-BH[i].sum(),i))
            
        for i in xrange(r):
            _,l = heapq.heappop(sums)
            res[i] = BH[l]
            used.append(l)
            
        del sums     
        
        """
        # get random r shape contexts
        # this not good because gives different result for same query
        while len(used) < r:
            i = random.randint(0,len(BH))
            if i not in used:
                res[len(used)] = BH[i].reshape(self.nbins)
                used.append(i)
        """
        
        return res,used


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

    a = SC()
    sampls = 100

    imgs = get_elements('test.png')

    points = get_points_from_img('9M2.png',simpleto=sampls)
    P = a.compute(points)
    x1 = [p[0] for p in points]
    y1 = [400-p[1] for p in points]
        
    for img in imgs:
    
        points2 = get_points_from_img(img,simpleto=sampls)
        
        if points2:
            Q = a.compute(points2)
            x2 = [p[0] for p in points2]
            y2 = [400-p[1] for p in points2]
            
            # get rendom r shape contexts from query shape
            Qs,points_ids = a.get_contextes(Q,5)
            points2s = [points2[i] for i in points_ids]
            COST,indexes = a.quick_diff(P,Qs)
            
            LINES = []
            for p1,q1 in indexes:
                LINES.append([[points[p1][0],400-points[p1][1]],[points2s[q1][0],400-points2s[q1][1]]])

            make_graph((x1,y1),(x2,y2),COST,LINES)    
    """
    COST,indexes = a.diff(P,Q)
    
    LINES = []
    for p1,q1 in indexes:
        LINES.append([[points[p1][0],400-points[p1][1]],[points2[q1][0],400-points2[q1][1]]])

    make_graph((x1,y1),(x2,y2),COST,LINES)  
    
    """
    
    
