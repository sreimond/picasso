# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 09:11:40 2018

@author: sreimond
"""

import numpy as np
from scipy import special

def colorNames():
    return {'aliceblue':'#F0F8FF',
            'antiquewhite':'#FAEBD7',
            'aqua':'#00FFFF',
            'aquamarine':'#7FFFD4',
            'azure':'#F0FFFF',
            'beige':'#F5F5DC',
            'bisque':'#FFE4C4',
            'black':'#000000',
            'blanchedalmond':'#FFEBCD',
            'blue':'#0000FF',
            'blueviolet':'#8A2BE2',
            'brown':'#A52A2A',
            'burlywood':'#DEB887',
            'cadetblue':'#5F9EA0',
            'chartreuse':'#7FFF00',
            'chocolate':'#D2691E',
            'coral':'#FF7F50',
            'cornflowerblue':'#6495ED',
            'cornsilk':'#FFF8DC',
            'crimson':'#DC143C',
            'cyan':'#00FFFF',
            'darkblue':'#00008B',
            'darkcyan':'#008B8B',
            'darkgoldenrod':'#B8860B',
            'darkgray':'#A9A9A9',
            'darkgrey':'#A9A9A9',
            'darkgreen':'#006400',
            'darkkhaki':'#BDB76B',
            'darkmagenta':'#8B008B',
            'darkolivegreen':'#556B2F',
            'darkorange':'#FF8C00',
            'darkorchid':'#9932CC',
            'darkred':'#8B0000',
            'darksalmon':'#E9967A',
            'darkseagreen':'#8FBC8F',
            'darkslateblue':'#483D8B',
            'darkslategray':'#2F4F4F',
            'darkslategrey':'#2F4F4F',
            'darkturquoise':'#00CED1',
            'darkviolet':'#9400D3',
            'deeppink':'#FF1493',
            'deepskyblue':'#00BFFF',
            'dimgray':'#696969',
            'dimgrey':'#696969',
            'dodgerblue':'#1E90FF',
            'firebrick':'#B22222',
            'floralwhite':'#FFFAF0',
            'forestgreen':'#228B22',
            'fuchsia':'#FF00FF',
            'gainsboro':'#DCDCDC',
            'ghostwhite':'#F8F8FF',
            'gold':'#FFD700',
            'goldenrod':'#DAA520',
            'gray':'#808080',
            'grey':'#808080',
            'green':'#008000',
            'greenyellow':'#ADFF2F',
            'honeydew':'#F0FFF0',
            'hotpink':'#FF69B4',
            'indianred':'#CD5C5C',
            'indigo':'#4B0082',
            'ivory':'#FFFFF0',
            'khaki':'#F0E68C',
            'lavender':'#E6E6FA',
            'lavenderblush':'#FFF0F5',
            'lawngreen':'#7CFC00',
            'lemonchiffon':'#FFFACD',
            'lightblue':'#ADD8E6',
            'lightcoral':'#F08080',
            'lightcyan':'#E0FFFF',
            'lightgoldenrodyellow':'#FAFAD2',
            'lightgray':'#D3D3D3',
            'lightgrey':'#D3D3D3',
            'lightgreen':'#90EE90',
            'lightpink':'#FFB6C1',
            'lightsalmon':'#FFA07A',
            'lightseagreen':'#20B2AA',
            'lightskyblue':'#87CEFA',
            'lightslategray':'#778899',
            'lightslategrey':'#778899',
            'lightsteelblue':'#B0C4DE',
            'lightyellow':'#FFFFE0',
            'lime':'#00FF00',
            'limegreen':'#32CD32',
            'linen':'#FAF0E6',
            'magenta':'#FF00FF',
            'maroon':'#800000',
            'mediumaquamarine':'#66CDAA',
            'mediumblue':'#0000CD',
            'mediumorchid':'#BA55D3',
            'mediumpurple':'#9370DB',
            'mediumseagreen':'#3CB371',
            'mediumslateblue':'#7B68EE',
            'mediumspringgreen':'#00FA9A',
            'mediumturquoise':'#48D1CC',
            'mediumvioletred':'#C71585',
            'midnightblue':'#191970',
            'mintcream':'#F5FFFA',
            'mistyrose':'#FFE4E1',
            'moccasin':'#FFE4B5',
            'navajowhite':'#FFDEAD',
            'navy':'#000080',
            'oldlace':'#FDF5E6',
            'olive':'#808000',
            'olivedrab':'#6B8E23',
            'orange':'#FFA500',
            'orangered':'#FF4500',
            'orchid':'#DA70D6',
            'palegoldenrod':'#EEE8AA',
            'palegreen':'#98FB98',
            'paleturquoise':'#AFEEEE',
            'palevioletred':'#DB7093',
            'papayawhip':'#FFEFD5',
            'peachpuff':'#FFDAB9',
            'peru':'#CD853F',
            'pink':'#FFC0CB',
            'plum':'#DDA0DD',
            'powderblue':'#B0E0E6',
            'purple':'#800080',
            'rebeccapurple':'#663399',
            'red':'#FF0000',
            'rosybrown':'#BC8F8F',
            'royalblue':'#4169E1',
            'saddlebrown':'#8B4513',
            'salmon':'#FA8072',
            'sandybrown':'#F4A460',
            'seagreen':'#2E8B57',
            'seashell':'#FFF5EE',
            'sienna':'#A0522D',
            'silver':'#C0C0C0',
            'skyblue':'#87CEEB',
            'slateblue':'#6A5ACD',
            'slategray':'#708090',
            'slategrey':'#708090',
            'snow':'#FFFAFA',
            'springgreen':'#00FF7F',
            'steelblue':'#4682B4',
            'tan':'#D2B48C',
            'teal':'#008080',
            'thistle':'#D8BFD8',
            'tomato':'#FF6347',
            'turquoise':'#40E0D0',
            'violet':'#EE82EE',
            'wheat':'#F5DEB3',
            'white':'#FFFFFF',
            'whitesmoke':'#F5F5F5',
            'yellow':'#FFFF00',
            'yellowgreen':'#9ACD32'}

def rgb2hsv(rgb):
    rgb0 = np.array(rgb)/255.0
    r = rgb0[0]
    g = rgb0[1]
    b = rgb0[2]
    cmax = np.amax(rgb0)
    cmin = np.amin(rgb0)
    delta = cmax-cmin
    v = cmax
    if delta==0:
        h = 0.0
        s = 0.0
    else:
        s = delta/cmax
        if cmax==r:
            h = 60.0 * (((g-b)/delta) % 6.0)
        elif cmax==g:
            h = 60.0 * (((b-r)/delta) + 2.0)
        elif cmax==b:
            h = 60.0 * (((r-g)/delta) + 4.0)
    return np.array([h,s,v])

def hsv2rgb(hsv):
    h = hsv[0]
    s = hsv[1]
    v = hsv[2]
    c = v * s
    x = c * (1.0-np.fabs(((h/60.0) % 2.0) - 1.0))
    m = v - c
    if h<60.0:
        r = c
        g = x
        b = 0
    elif h>=60.0 and h<120.0:
        r = x
        g = c
        b = 0
    elif h>=120.0 and h<180.0:
        r = 0
        g = c
        b = x
    elif h>=180.0 and h<240.0:
        r = 0
        g = x
        b = c
    elif h>=240.0 and h<300.0:
        r = x
        g = 0
        b = c
    elif h>=300.0:
        r = c
        g = 0
        b = x
    return np.array([(r+m),(g+m),(b+m)]) * 255.0

def rgb2yuv(rgb):
    R = np.ones((3,3))    
    R[0,:] = [0.299,0.587,0.114]
    R[1,:] = [-0.14713,-0.28889,0.436]
    R[2,:] = [0.615,-0.51499,-0.10001]
    rgbv = np.ones((3,1))
    rgbv[0] = rgb[0]
    rgbv[1] = rgb[1]
    rgbv[2] = rgb[2]
    return np.dot(R,rgbv).ravel()

def yuv2rgb(yuv):
    R = np.ones((3,3))    
    R[0,:] = [1.0,0.0,1.13983]
    R[1,:] = [1.0,-0.39465,-0.58060]
    R[2,:] = [1.0,2.03211,0.0]
    yuvv = np.ones((3,1))
    yuvv[0] = yuv[0]
    yuvv[1] = yuv[1]
    yuvv[2] = yuv[2]
    return np.dot(R,yuv).ravel()
    
def rgb2hex(rgb):
    rgb = np.array(np.round(rgb),dtype=int)
    return '#%02x%02x%02x' % (rgb[0],rgb[1],rgb[2])

def hex2rgb(hexcode):
    r = int(hexcode[1:3],16)
    g = int(hexcode[3:5],16)
    b = int(hexcode[5:7],16)
    return np.array([r,g,b])

def hex2hsv(hexcode):
    return rgb2hsv(hex2rgb(hexcode))
    
def hsv2hex(hsv):
    return rgb2hex(hsv2rgb(hsv))

def hsv2yuv(hsv):
    return rgb2yuv(hsv2rgb(hsv))
    
def yuv2hsv(yuv):
    return rgb2hsv(yuv2rgb(yuv))

def hex2yuv(hexcode):
    return rgb2yuv(hex2rgb(hexcode))

def yuv2hex(yuv):
    return rgb2hex(yuv2rgb(yuv))
    
def name2rgb(name):
    return hex2rgb(colorNames()[name.lower()])

def name2hsv(name):
    return hex2hsv(colorNames()[name.lower()])

def name2yuv(name):
    return rgb2yuv(name2rgb(name))    

def rgb2name(rgb):
    # compares rgb to list and finds 'nearest neighbour' from list
    yuvi = rgb2yuv(rgb)
    d = np.inf
    for key,value in colorNames().items():
        yuvj = name2yuv(key)
        dij = ((yuvi[0]-yuvj[0])**2.0 + 
               (yuvi[1]-yuvj[1])**2.0 +  
               (yuvi[2]-yuvj[2])**2.0) ** 0.5
        if dij < d:
            d = dij
            name = key
    return name

def hsv2name(hsv):
    return rgb2name(hsv2rgb(hsv))

def hex2name(hexcode):
    return rgb2name(hex2rgb(hexcode))
    
def yuv2name(yuv):
    return rgb2name(yuv2rgb(yuv))
    
def color_gradient(colors,color_count=100,color_representation='rgb',interp='rgb'):
    if interp=='rgb':
        return color_gradient_rgb(colors,
                                  color_count=color_count,
                                  color_representation=color_representation)
    elif interp=='hsv':
        return color_gradient_hsv(colors,
                                  color_count=color_count,
                                  color_representation=color_representation)
    elif interp=='yuv':
        return color_gradient_yuv(colors,
                                  color_count=color_count,
                                  color_representation=color_representation)

def color_gradient_rgb(colors,color_count=100,color_representation='rgb'):
    nc = len(colors)
    xp = range(nc)
    xi = np.linspace(0,nc-1,num=color_count,endpoint=True)
    fpr = []
    fpg = []
    fpb = []
    gradient = ()
    if color_representation.lower()=='rgb':
        for rgb in colors:
            fpr.append(rgb[0])
            fpg.append(rgb[1])
            fpb.append(rgb[2])
    elif color_representation.lower()=='hsv':
        for hsv in colors:
            rgb = hsv2rgb(hsv)
            fpr.append(rgb[0])
            fpg.append(rgb[1])
            fpb.append(rgb[2])
    elif color_representation.lower()=='hex':
        for hexcode in colors:
            rgb = hex2rgb(hexcode)
            fpr.append(rgb[0])
            fpg.append(rgb[1])
            fpb.append(rgb[2])
    elif color_representation.lower()=='yuv':
        for yuv in colors:
            rgb = yuv2rgb(yuv)
            fpr.append(rgb[0])
            fpg.append(rgb[1])
            fpb.append(rgb[2])
    elif color_representation.lower()=='name':
        for name in colors:
            rgb = name2rgb(name)
            fpr.append(rgb[0])
            fpg.append(rgb[1])
            fpb.append(rgb[2])
    fir = np.interp(xi,xp,fpr)
    fig = np.interp(xi,xp,fpg)
    fib = np.interp(xi,xp,fpb)
    for i in range(color_count):
        gradient += ([fir[i],fig[i],fib[i]],)
    return gradient    
    
def color_gradient_hsv(colors,color_count=100,color_representation='rgb'):
    nc = len(colors)
    xp = range(nc)
    xi = np.linspace(0,nc-1,num=color_count,endpoint=True)
    fph = []
    fps = []
    fpv = []
    gradient = ()
    if color_representation.lower()=='rgb':
        for rgb in colors:
            hsv = rgb2hsv(rgb)
            fph.append(hsv[0])
            fps.append(hsv[1])
            fpv.append(hsv[2])
    elif color_representation.lower()=='hsv':
        for hsv in colors:
            fph.append(hsv[0])
            fps.append(hsv[1])
            fpv.append(hsv[2])
    elif color_representation.lower()=='yuv':
        for yuv in colors:
            hsv = yuv2hsv(yuv)
            fph.append(hsv[0])
            fps.append(hsv[1])
            fpv.append(hsv[2])
    elif color_representation.lower()=='hex':
        for hexcode in colors:
            hsv = hex2hsv(hexcode)
            fph.append(hsv[0])
            fps.append(hsv[1])
            fpv.append(hsv[2])
    elif color_representation.lower()=='name':
        for name in colors:
            hsv = name2hsv(name)
            fph.append(hsv[0])
            fps.append(hsv[1])
            fpv.append(hsv[2])       
    fih = np.interp(xi,xp,fph)
    fis = np.interp(xi,xp,fps)
    fiv = np.interp(xi,xp,fpv)
    for i in range(color_count):
        gradient += (hsv2rgb([fih[i],fis[i],fiv[i]]),)
    return gradient

def color_gradient_yuv(colors,color_count=100,color_representation='rgb'):
    nc = len(colors)
    xp = range(nc)
    xi = np.linspace(0,nc-1,num=color_count,endpoint=True)
    fpy = []
    fpu = []
    fpv = []
    gradient = ()
    if color_representation.lower()=='rgb':
        for rgb in colors:
            yuv = rgb2yuv(rgb)
            fpy.append(yuv[0])
            fpu.append(yuv[1])
            fpv.append(yuv[2])
    elif color_representation.lower()=='hsv':
        for hsv in colors:
            yuv = hsv2yuv(hsv)
            fpy.append(yuv[0])
            fpu.append(yuv[1])
            fpv.append(yuv[2])
    elif color_representation.lower()=='yuv':
        for yuv in colors:
            fpy.append(yuv[0])
            fpu.append(yuv[1])
            fpv.append(yuv[2])
    elif color_representation.lower()=='hex':
        for hexcode in colors:
            yuv = hex2yuv(hexcode)
            fpy.append(yuv[0])
            fpu.append(yuv[1])
            fpv.append(yuv[2])
    elif color_representation.lower()=='name':
        for name in colors:
            yuv = name2yuv(name)
            fpy.append(yuv[0])
            fpu.append(yuv[1])
            fpv.append(yuv[2])    
    fiy = np.interp(xi,xp,fpy)
    fiu = np.interp(xi,xp,fpu)
    fiv = np.interp(xi,xp,fpv)
    for i in range(color_count):
        gradient += (yuv2rgb([fiy[i],fiu[i],fiv[i]]),)
    return gradient

def color_schemes_tol(n,nature):
    # https://personal.sron.nl/~pault/
    x = np.linspace(0,1,n)
    if nature=='qualitative':
        cmap = qualitative_table_tol(n)
    elif nature=='sequential':
        r = 1.0 - 0.392 * (1.0+special.erf((x-0.869)/(0.255)))
        g = 1.021 - 0.456 * (1.0+special.erf((x-0.527)/(0.376)))
        b = 1.0 - 0.493 * (1.0+special.erf((x-0.272)/(0.309)))
        cmap = ()
        for i in range(n):
            cmap += ([r[i]*255,g[i]*255,b[i]*255],)
    elif nature=='diverging':        
        r = 0.237 - 2.13*x + 26.92*x**2.0 - 65.5*x**3.0 + 63.5*x**4.0 - 22.36*x**5.0
        g = ((0.572+1.524*x-1.811*x**2.0) / (1.0 - 0.291*x + 0.1574*x**2.0))**2.0
        b = 1.0 / (1.579-4.03*x + 12.92*x**2.0 - 31.4*x**3.0 + 48.6*x**4.0 - 23.36*x**5.0)
        cmap = ()
        for i in range(n):
            cmap += ([r[i]*255,g[i]*255,b[i]*255],)
    elif nature=='rainbow':        
        r = (0.472 - 0.567*x + 4.05*x**2.0) / (1.0 + 8.72*x - 19.17*x**2.0 + 14.1*x**2.0)
        g = 0.108932 - 1.22635*x + 27.284*x**2.0 - 98.577*x**3.0 + 163.3*x**4.0 - 131.395*x**5.0 + 40.634*x**6.0
        b = 1.0 / (1.97 + 3.54*x - 68.5*x**2.0 + 243*x**3.0 - 297*x**4 + 125*x**5.0)
        cmap = ()
        for i in range(n):
            cmap += ([r[i]*255,g[i]*255,b[i]*255],)
    elif nature=='hues':
        cmap = hues_table_tol(n);
    return cmap
        
def qualitative_table_tol(n):
    if n==1:        
        r = [68]
        g = [119]
        b = [170]
    elif n==2:
        r = [68, 204]
        g = [119, 102]
        b = [170, 119]
    elif n==3:
        r = [68, 221, 204]
        g = [119, 204, 102]
        b = [170, 119, 119]
    elif n==4:
        r = [68, 17, 221, 204]
        g = [119, 119, 204, 102]
        b = [170, 51, 119, 119]
    elif n==5:
        r = [51, 136, 17, 221, 204]
        g = [34, 204, 119, 204, 102]
        b = [136, 238, 51, 119, 119]
    elif n==6:
        r = [51, 136, 17, 221, 204, 170]
        g = [34, 204, 119, 204, 102, 68]
        b = [136, 238, 51, 119, 119, 153]
    elif n==7:
        r = [51, 136, 68, 17, 221, 204, 170]
        g = [34, 204, 170, 119, 204, 102, 68]
        b = [136, 238, 153, 51, 119, 119, 153]
    elif n==8:
        r = [51, 136, 68, 17, 153, 221, 204, 170]
        g = [34, 204, 170, 119, 153, 204, 102, 68]
        b = [136, 238, 153, 51, 51, 119, 119, 153]
    elif n==9:
        r = [51, 136, 68, 17, 153, 221, 204, 136, 170]
        g = [34, 204, 170, 119, 153, 204, 102, 34, 68]
        b = [136, 238, 153, 51, 51, 119, 119, 85, 153]
    elif n==10:
        r = [51, 136, 68, 17, 153, 221, 102, 204, 136, 170]
        g = [34, 204, 170, 119, 153, 204, 17, 102, 34, 68]
        b = [136, 238, 153, 51, 51, 119, 0, 119, 85, 153]    
    elif n==11:
        r = [51, 102, 136, 68, 17, 153, 221, 102, 204, 136, 170]
        g = [34, 153, 204, 170, 119, 153, 204, 17, 102, 34, 68]
        b = [136, 204, 238, 153, 51, 51, 119, 0, 119, 85, 153]
    elif n==12:
        r = [51, 102, 136, 68, 17, 153, 221, 102, 204, 170, 136, 170]
        g = [34, 153, 204, 170, 119, 153, 204, 17, 102, 68, 34, 68]
        b = [136, 204, 238, 153, 51, 51, 119, 0, 119, 102, 85, 153]
    else:
        r = np.random.randint(0,256,n)
        g = np.random.randint(0,256,n)
        b = np.random.randint(0,256,n)
    cmap = ()
    for i in range(n):
        cmap += ([r[i],g[i],b[i]],)
    return cmap
        
def hues_table_tol(n):
    if n==1:
        r = [119, 68, 17]
        g = [170, 119, 68]
        b = [221, 170, 119]    
    elif n==2:
        r = [119, 68, 17, 119, 68, 17]
        g = [170, 119, 68, 204, 170, 119]
        b = [221, 170, 119, 204, 170, 119]        
    elif n==3:
        r = [119, 68, 17, 119, 68, 17, 136, 68, 17]
        g = [170, 119, 68, 204, 170, 119, 204, 170, 119]
        b = [221, 170, 119, 204, 170, 119, 170, 119, 68]    
    elif n==4:
        r = [119, 68, 17, 119, 68, 17, 136, 68, 17, 221, 170, 119]
        g = [170, 119, 68, 204, 170, 119, 204, 170, 119, 221, 170, 119]
        b = [221, 170, 119, 204, 170, 119, 170, 119, 68, 119, 68, 17]
    elif n==5:        
        r = [119, 68, 17, 119, 68, 17, 136, 68, 17, 221, 170, 119, 221, 170, 119]
        g = [170, 119, 68, 204, 170, 119, 204, 170, 119, 221, 170, 119, 170, 119, 68]
        b = [221, 170, 119, 204, 170, 119, 170, 119, 68, 119, 68, 17, 119, 68, 17]
    elif n==6:
        r = [119, 68, 17, 119, 68, 17, 136, 68, 17, 221, 170, 119, 221, 170, 119, 221, 170, 119]
        g = [170, 119, 68, 204, 170, 119, 204, 170, 119, 221, 170, 119, 170, 119, 68, 119, 68, 17]
        b = [221, 170, 119, 204, 170, 119, 170, 119, 68, 119, 68, 17, 119, 68, 17, 136, 85, 34]
    elif n==7:
        r = [119, 68, 17, 119, 68, 17, 136, 68, 17, 221, 170, 119, 221, 170, 119, 221, 170, 119, 221, 170, 119]
        g = [170, 119, 68, 204, 170, 119, 204, 170, 119, 221, 170, 119, 170, 119, 68, 119, 68, 17, 119, 68, 17]
        b = [221, 170, 119, 204, 170, 119, 170, 119, 68, 119, 68, 17, 119, 68, 17, 136, 85, 34, 136, 85, 34]
    elif n==8:
        r = [119, 68, 17, 119, 68, 17, 136, 68, 17, 221, 170, 119, 221, 170, 119, 221, 170, 119, 221, 170, 119, 204, 170, 119]
        g = [170, 119, 68, 204, 170, 119, 204, 170, 119, 221, 170, 119, 170, 119, 68, 119, 68, 17, 119, 68, 17, 153, 68, 17]
        b = [221, 170, 119, 204, 170, 119, 170, 119, 68, 119, 68, 17, 119, 68, 17, 136, 85, 34, 136, 85, 34, 187, 136, 85]
    else:
        r = np.random.randint(0,256,n)
        g = np.random.randint(0,256,n)
        b = np.random.randint(0,256,n)
    cmap = ()
    for i in range(n):
        cmap += ([r[i],g[i],b[i]],)
    return cmap
