import sys
from math import sin, cos, sqrt, pi
import cv
import urllib2
import time
import math
from numpy import *
from scipy.spatial.distance import euclidean
import heapq as hq

CANNY       = 1

def get_elements(filename,treshold=50,minheight=15,minarea=200,elements=6):
    src = cv.LoadImage(filename, cv.CV_LOAD_IMAGE_GRAYSCALE)
    test = cv.CreateImage(cv.GetSize(src),32,3)  
    dst = cv.CreateImage(cv.GetSize(src), 8, 1)
    storage = cv.CreateMemStorage(0)
    cv.Canny(src, dst, treshold, treshold*3, 3)

    storage = cv.CreateMemStorage(0)
    seqs = cv.FindContours(dst, storage,cv.CV_RETR_TREE, cv.CV_CHAIN_APPROX_NONE, (0,0))

    res = []

    c = seqs.h_next()
    while True:
        if not c:
            break
        box = cv.BoundingRect(c)
        area = box[2]*box[3]
        #and (area > minarea)
        if (box[3] > minheight):
            res.append(box) 
        c = c.h_next()
    
    if len(res) < elements:
        while len(res) < elements:
            m = 0
            c = 0
            for i,e in enumerate(res):
                if e[3] > m:
                    m = e[3]
                    c = i
                    
            big = res.pop(c)
            res.append((big[0],big[1],int(big[2]*1.0/2),big[3]))        
            res.append((big[0]+int(big[2]*1.0/2),big[1],int(big[2]*1.0/2),big[3])) 
    
    #for box in res:
    #    cv.Rectangle(dst, (box[0],box[1]), (box[0]+box[2],box[1]+box[3]), cv.RGB(255,255,255))     
        
    #cv.ShowImage('Preview2',dst)
    #cv.WaitKey()   
        
    imgs = []
    print len(res)
    for box in res:
        cv.SetImageROI(src, box);
     
        tmp = cv.CreateImage((box[2],box[3]),8,1)
         
        cv.Copy(src, tmp);
        hq.heappush(imgs,(box[0],tmp))
        
        cv.ResetImageROI(src);
        
    res = [hq.heappop(imgs)[1] for i in xrange(len(res))]  
    return res   
    

def euclid_distance(p1,p2):
    return math.sqrt( ( p2[0] - p1[0] ) ** 2 + ( p2[1] - p1[1] ) ** 2 )
    
    
def get_points_from_img(src,treshold=50,simpleto=100,t=CANNY):
    ts = time.time()
    if isinstance(src,str):
        src = cv.LoadImage(src, cv.CV_LOAD_IMAGE_GRAYSCALE)
    test = cv.CreateImage(cv.GetSize(src),8,1)  
    if t == CANNY:
        dst = cv.CreateImage(cv.GetSize(src), 8, 1)
        storage = cv.CreateMemStorage(0)
        cv.Canny(src, dst, treshold, treshold*3, 3)
        
    A = zeros((cv.GetSize(src)[1],cv.GetSize(src)[0]))
    for y in xrange(cv.GetSize(src)[1]):
        for x in xrange(cv.GetSize(src)[0]):
            A[y,x] = src[y,x]    
    
    px,py = gradient(A)
    points = []
    w,h = cv.GetSize(src)
    for y in xrange(h):
        for x in xrange(w):
            try:
                c = dst[y,x]
            except:
                print x,y
            if c == 255:
                points.append((x,y))
    
    r = 2
    while len(points) > simpleto:
        newpoints = points
        xr = range(0,w,r)
        yr = range(0,h,r)
        for p in points:
            if p[0] not in xr and p[1] not in yr:
                newpoints.remove(p)
                if len(points) <= simpleto:
                    T = zeros((simpleto,1)) 
                    for i,(x,y) in enumerate(points):
                        T[i] = math.atan2(py[y,x],px[y,x])+pi/2;    
                    return points,asmatrix(T)
        r += 1
    T = zeros((simpleto,1)) 
    for i,(x,y) in enumerate(points):
        T[i] = math.atan2(py[y,x],px[y,x])+pi/2;    
        
    return points,asmatrix(T)
    
    
def dist2(x,c):
    """
        Euclidian distance matrix
    """
    ncentres = c.shape[0]
    ndata = x.shape[0]
    return (ones((ncentres, 1)) * (((power(x,2)).H)).sum(axis=0)).H + ones((ndata, 1)) * ((power(c,2)).H).sum(axis=0) - multiply(2,(x*(c.H)));
    
def bookenstain(X,Y,beta):
    """
        Bookstein PAMI89
    
        Article: Principal Warps: Thin-Plate Splines and the Decomposition of Deformations

    """
    X = asmatrix(X)
    Y = asmatrix(Y)

    N = X.shape[0]
    r2 = dist2(X,X)
    K = multiply(r2,log(r2+eye(N,N)))
    P = concatenate((ones((N,1)),X),1)
    L = bmat([[K, P], [P.H, zeros((3,3))]])
    V = concatenate((Y.H,zeros((2,3))),1)

    L[0:N,0:N] = L[0:N,0:N] + beta * eye(N,N)

    invL = linalg.inv(L)

    # L^-1 * v^T = (W | a_1 a_x a_y)^T
    c = invL*(V.H)
    cx = c[:,0]
    cy = c[:,1]
    
    Q = (c[0:N,:].H) * K * c[0:N,:]
    E = mean(diag(Q))

    n_good = 10

    A=concatenate((cx[n_good+2:n_good+3,:],cy[n_good+2:n_good+3,:]),1);
    s=linalg.svd(A);
    aff_cost=log(s[0]/s[1])

    return cx,cy,E,aff_cost,L
    
def gauss_kernel(N):
    """
        Gaussian kernel
    """
    g=2**(1-N)*diag(fliplr(pascal(N)));
    W=g*g.H;
   
   
def pascal(n, k = 0):
    """
        Pascal matrix
    """
    p = diag( (-1)**arange(n) )
    p[:, 0] = ones(n)

    #  Generate the Pascal Cholesky factor (up to signs).
    for j in range(1, n - 1):
        for i in range(j+1, n):
            p[i, j] = p[i-1, j] - p[i-1, j-1]

    if k == 0:
        p = matrix(p) * matrix(p.T)

    elif k == 2:
        p = rot90(p, 3)
        if n/2 == round(n/2):
            p = -p

    return p
