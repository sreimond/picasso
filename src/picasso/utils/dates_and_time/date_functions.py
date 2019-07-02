# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 13:45:27 2018

@author: sreimond
"""

import numpy as np
from datetime import datetime

def is_leap( years, cal='auto' ):
    """
    The `is_leap` function enables array input.
    Documentation see the `_is_leap` function.
    """
    years = np.array(years,ndmin=1)
    years_count = np.size(years)
    ret = np.zeros(years_count,dtype=np.bool_)
    for ix in range(years_count):
        try:
            ret[ix] = _is_leap( years[ix], cal=cal )
        except:
            ret[ix] = np.nan
    return ret
    
def is_julian( y, m, d ):
    """
    The `is_julian` function enables array input.
    Documentation see the `_is_julian` function.
    """
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    days = np.array(d,ndmin=1)
    years_count = np.size(years)
    dim_check = ((years_count == np.size(months)) 
                 and (years_count == np.size(days)))
    if not dim_check:
        raise ValueError('dimension mismatch')
    ret = np.zeros(years_count,dtype=np.bool_)
    for ix in range(years_count):
        try:
            ret[ix] = _is_julian( years[ix], months[ix], days[ix] )
        except:
            ret[ix] = np.nan            
    return ret

def is_gregorian( y, m, d ):
    """
    The `is_gregorian` function enables array input.
    Documentation see the `_is_julian` function.
    """
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    days = np.array(d,ndmin=1)
    years_count = np.size(years)
    dim_check = ((years_count == np.size(months)) 
                 and (years_count == np.size(days)))
    if not dim_check:
        raise ValueError('dimension mismatch')
    ret = np.zeros(years_count,dtype=np.bool_)
    for ix in range(years_count):
        try:
            ret[ix] = _is_gregorian( years[ix], months[ix], days[ix] )
        except:
            ret[ix] = np.nan
    return ret    

def ymd2jd( y, m, d, cal='auto' ): 
    """
    The `ymd2jd` function enables array input.
    Documentation see the `_ymd2jd` function.
    """
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    days = np.array(d,ndmin=1)
    years_count = np.size(years)
    dim_check = ((years_count == np.size(months)) 
                 and (years_count == np.size(days)))
    if not dim_check:
        raise ValueError('dimension mismatch')
    ret = np.zeros(years_count,dtype=np.float_)
    for ix in range(years_count):
        try:
            ret[ix] = _ymd2jd( years[ix], months[ix], days[ix], cal=cal )
        except:
            ret[ix] = np.nan
    return ret
    
def jd2ymd( jd, cal='auto' ):
    """
    The `jd2ymd` function enables array input.
    Documentation see the `_jd2ymd` function.
    """
    jd = np.array(jd,ndmin=1)
    jd_count = np.size(jd)
    years = np.zeros(jd_count,dtype=np.int_)
    months = np.zeros(jd_count,dtype=np.int_)
    days = np.zeros(jd_count,dtype=np.float_)
    for ix in range(jd_count):
        try:
            years[ix], months[ix], days[ix] = _jd2ymd( jd[ix], cal=cal )
        except:
            years[ix], months[ix], days[ix] = np.nan, np.nan, np.nan
    return years, months, days
        
def ymd2mjd( y, m, d, cal='auto' ): 
    """
    The `ymd2mjd` function enables array input.
    Documentation see the `_ymd2mjd` function.
    """
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    days = np.array(d,ndmin=1)
    years_count = np.size(years)
    dim_check = ((years_count == np.size(months)) 
                 and (years_count == np.size(days)))
    if not dim_check:
        raise ValueError('dimension mismatch')
    ret = np.zeros(years_count,dtype=np.float_)
    for ix in range(years_count):
        try:
            ret[ix] = _ymd2mjd( years[ix], months[ix], days[ix], cal=cal )
        except:
            ret[ix] = np.nan
    return ret

def mjd2ymd( mjd, cal='auto' ):
    """
    The `mjd2ymd` function enables array input.
    Documentation see the `_mjd2ymd` function.
    """
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)
    years = np.zeros(mjd_count,dtype=np.int_)
    months = np.zeros(mjd_count,dtype=np.int_)
    days = np.zeros(mjd_count,dtype=np.float_)
    for ix in range(mjd_count):
        try:
            years[ix], months[ix], days[ix] = _mjd2ymd( mjd[ix], cal=cal )
        except:
            years[ix], months[ix], days[ix] = np.nan, np.nan, np.nan
    return years, months, days

def jd2dow( jd ):
    """
    The `jd2dow` function enables array input.
    Documentation see the `_jd2dow` function.
    """
    jd = np.array(jd,ndmin=1)
    jd_count = np.size(jd)
    days_number = np.zeros(jd_count,dtype=np.int_)
    days_name = np.zeros(jd_count,dtype='|S3')    
    for ix in range(jd_count):
        try:
            days_number[ix], days_name[ix] = _jd2dow( jd[ix] )
        except:
            days_number[ix], days_name[ix] = np.nan, 'nan'
    return days_number, days_name
    
def mjd2dow( mjd ):
    """
    The `mjd2dow` function enables array input.
    Documentation see the `_mjd2dow` function.
    """
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)
    days_number = np.zeros(mjd_count,dtype=np.int_)
    days_name = np.zeros(mjd_count,dtype='|S3')    
    for ix in range(mjd_count):
        try:
            days_number[ix], days_name[ix] = _mjd2dow( mjd[ix] )
        except:
            days_number[ix], days_name[ix] = np.nan, 'nan'
    return days_number, days_name

def dhms2day( d, h, m, s ):
    """
    The `dhms2day` function enables array input.
    Documentation see the `_dhms2day` function.
    """
    days = np.array(d,ndmin=1)
    hours = np.array(h,ndmin=1)
    minutes = np.array(m,ndmin=1)
    seconds = np.array(s,ndmin=1)
    days_count = np.size(days)
    dim_check = ((days_count == np.size(hours)) 
                 and (days_count == np.size(minutes))
                 and (days_count == np.size(seconds)))
    if not dim_check:
        raise ValueError('dimension mismatch')
    ret = np.zeros(days_count,dtype=np.float_)
    for ix in range(days_count):
        try:
            ret[ix] = _dhms2day( days[ix], hours[ix], minutes[ix], seconds[ix] )
        except:
            ret[ix] = np.nan
    return ret

def day2dhms( day ):
    """
    The `day2dhms` function enables array input.
    Documentation see the `_day2dhms` function.
    """
    day = np.array(day,ndmin=1)
    day_count = np.size(day)
    days = np.zeros(day_count,dtype=np.int_)
    hours = np.zeros(day_count,dtype=np.int_)
    minutes = np.zeros(day_count,dtype=np.int_)
    seconds = np.zeros(day_count,dtype=np.float_)
    for ix in range(day_count):
        try:
            days[ix], hours[ix], minutes[ix], seconds[ix] = _day2dhms( day[ix] )
        except:
            days[ix], hours[ix], minutes[ix], seconds[ix] = np.nan, np.nan, np.nan, np.nan
    return days, hours, minutes, seconds
    
def ymd2doy( y, m, d ):
    """
    The `ymd2doy` function enables array input.
    Documentation see the `_ymd2doy` function.
    """
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    days = np.array(d,ndmin=1)
    years_count = np.size(years)
    dim_check = ((years_count == np.size(months)) 
                 and (years_count == np.size(days)))
    if not dim_check:
        raise ValueError('dimension mismatch')
    ret = np.zeros(years_count,dtype=np.int_)
    for ix in range(years_count):
        try:
            ret[ix] = _ymd2doy( years[ix], months[ix], days[ix] )
        except:
            ret[ix] = np.nan
    return ret
    
def jd2doy( jd ):
    """
    The `jd2doy` function enables array input.
    Documentation see the `_jd2doy` function.
    """
    jd = np.array(jd,ndmin=1)
    jd_count = np.size(jd)
    doys = np.zeros(jd_count,dtype=np.int_)
    for ix in range(jd_count):
        try:
            doys[ix] = _jd2doy( jd[ix] )
        except:
            doys[ix] = np.nan
    return doys
    
def mjd2doy( mjd ):
    """
    The `mjd2doy` function enables array input.
    Documentation see the `_mjd2doy` function.
    """
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)
    doys = np.zeros(mjd_count,dtype=np.int_)
    for ix in range(mjd_count):
        try:
            doys[ix] = _mjd2doy( mjd[ix] )
        except:
            doys[ix] = np.nan
    return doys
    
def doy2ymd( y, doy ):
    """
    The `doy2ymd` function enables array input.
    Documentation see the `_doy2ymd` function.
    """
    ys = np.array(y,ndmin=1)
    doys = np.array(doy,ndmin=1)
    ys_count = np.size(ys)
    dim_check = (ys_count == np.size(doys))
    if not dim_check:
        raise ValueError('dimension mismatch')        
    years = np.zeros(ys_count,dtype=np.int_)
    months = np.zeros(ys_count,dtype=np.int_)
    days = np.zeros(ys_count,dtype=np.float_)
    for ix in range(ys_count):
        try:
            years[ix], months[ix], days[ix] = _doy2ymd( ys[ix], doys[ix] )
        except:
            years[ix], months[ix], days[ix] = np.nan, np.nan, np.nan
    return years, months, days

def doy2jd( y, doy ):
    """
    The `doy2jd` function enables array input.
    Documentation see the `_doy2jd` function.
    """
    ys = np.array(y,ndmin=1)
    doys = np.array(doy,ndmin=1)
    ys_count = np.size(ys)    
    ret = np.zeros(ys_count,dtype=np.float_)
    for ix in range(ys_count):
        try:
            ret[ix] = _doy2jd( ys[ix], doys[ix] )
        except:
            ret[ix] = np.nan
    return ret
    
def doy2mjd( y, doy ):
    """
    The `doy2mjd` function enables array input.
    Documentation see the `_doy2mjd` function.
    """
    ys = np.array(y,ndmin=1)
    doys = np.array(doy,ndmin=1)
    ys_count = np.size(ys)    
    ret = np.zeros(ys_count,dtype=np.float_)
    for ix in range(ys_count):
        try:
            ret[ix] = _doy2mjd( ys[ix], doys[ix] )
        except:
            ret[ix] = np.nan
    return ret

def mjd2jd( mjd ):
    """
    The `mjd2jd` function enables array input.
    Documentation see the `_mjd2jd` function.
    """
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)    
    jd = np.zeros(mjd_count,dtype=np.float_)
    for ix in range(mjd_count):
        try:
            jd[ix] = _mjd2jd( mjd[ix] )
        except:
            jd[ix] = np.nan
    return jd

def jd2mjd( jd ):
    """
    The `jd2mjd` function enables array input.
    Documentation see the `_jd2mjd` function.
    """
    jd = np.array(jd,ndmin=1)
    jd_count = np.size(jd)    
    mjd = np.zeros(jd_count,dtype=np.float_)
    for ix in range(jd_count):
        try:
            mjd[ix] = _jd2mjd( jd[ix] )
        except:
            mjd[ix] = np.nan
    return mjd

def mjd2datetimenumber( mjd ):
    """
    The `mjd2datetimenumber` function enables array input.
    Documentation see the `_mjd2datetimenumber` function.
    """
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)    
    ret = np.zeros(mjd_count,dtype=datetime)
    for ix in range(mjd_count):
        try:
            ret[ix] = _mjd2datetimenumber( mjd[ix] )
        except:
            ret[ix] = np.nan
    return ret

def datetimenumber2mjd( t ):
    """
    The `datetimenumber2mjd` function enables array input.
    Documentation see the `_datetimenumber2mjd` function.
    """
    t = np.array(t,ndmin=1)
    t_count = np.size(t)    
    mjd = np.zeros(t_count,dtype=np.float_)
    for ix in range(t_count):
        try:
            mjd[ix] = _datetimenumber2mjd( t[ix] )
        except:
            mjd[ix] = np.nan
    return mjd

def ym2jdmrange( y, m ):
    """
    The `ym2jdmrange` function enables array input.
    Documentation see the `_ym2jdmrange` function.
    """    
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    years_count = np.size(years)
    dim_check = (years_count == np.size(months)) 
    if not dim_check:
        raise ValueError('dimension mismatch')
    jd1 = np.zeros(years_count,dtype=np.float_)
    jd2 = np.zeros(years_count,dtype=np.float_)    
    for ix in range(years_count):
        try:
            jd1[ix], jd2[ix] = _ym2jdmrange( years[ix], months[ix] )
        except:
            jd1[ix], jd2[ix] = np.nan, np.nan
    return jd1, jd2
    
def ym2mjdmrange( y, m ):
    """
    The `ym2mjdmrange` function enables array input.
    Documentation see the `_ym2mjdmrange` function.
    """    
    years = np.array(y,ndmin=1)
    months = np.array(m,ndmin=1)
    years_count = np.size(years)
    dim_check = (years_count == np.size(months)) 
    if not dim_check:
        raise ValueError('dimension mismatch')
    mjd1 = np.zeros(years_count,dtype=np.float_)
    mjd2 = np.zeros(years_count,dtype=np.float_)    
    for ix in range(years_count):
        try:
            mjd1[ix], mjd2[ix] = _ym2mjdmrange( years[ix], months[ix] )
        except:
            mjd1[ix], mjd2[ix] = np.nan, np.nan
    return mjd1, mjd2

def jd2jdmrange( jd ):
    """
    The `jd2jdmrange` function enables array input.
    Documentation see the `_jd2jdmrange` function.
    """    
    jd = np.array(jd,ndmin=1)
    jd_count = np.size(jd)    
    jd1 = np.zeros(jd_count,dtype=np.float_)
    jd2 = np.zeros(jd_count,dtype=np.float_)    
    for ix in range(jd_count):
        try:
            jd1[ix], jd2[ix] = _jd2jdmrange( jd[ix] )
        except:
            jd1[ix], jd2[ix] = np.nan, np.nan
    return jd1, jd2

def mjd2mjdmrange( mjd ):
    """
    The `mjd2mjdmrange` function enables array input.
    Documentation see the `_mjd2mjdmrange` function.
    """       
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)    
    mjd1 = np.zeros(mjd_count,dtype=np.float_)
    mjd2 = np.zeros(mjd_count,dtype=np.float_)    
    for ix in range(mjd_count):
        try:
            mjd1[ix], mjd2[ix] = _mjd2mjdmrange( mjd[ix] )
        except:
            mjd1[ix], mjd2[ix] = np.nan, np.nan
    return mjd1, mjd2
    
def ymstring2jd( ymstr ):
    """
    The `ymstring2jd` function enables array input.
    Documentation see the `_ymstring2jd` function.
    """    
    ymstr = np.array(ymstr,ndmin=1)
    ymstr_count = np.size(ymstr)    
    jd = np.zeros(ymstr_count,dtype=np.float_)
    for ix in range(ymstr_count):
        try:
            jd[ix] = _ymstring2jd( ymstr[ix] )
        except:
            jd[ix] = np.nan
    return jd   

def ymstring2mjd( ymstr ):
    """
    The `ymstring2mjd` function enables array input.
    Documentation see the `_ymstring2mjd` function.
    """    
    ymstr = np.array(ymstr,ndmin=1)
    ymstr_count = np.size(ymstr)    
    mjd = np.zeros(ymstr_count,dtype=np.float_)
    for ix in range(ymstr_count):
        try:
            mjd[ix] = _ymstring2mjd( ymstr[ix] )
        except:
            mjd[ix] = np.nan
    return mjd   

def jd2ymstring( jd, ymstring ):
    """
    The `jd2ymstring` function enables array input.
    Documentation see the `_jd2ymstring` function.
    """    
    jd = np.array(jd,ndmin=1)
    jd_count = np.size(jd)
    dummy = datetime.now().strftime(ymstring)
    type_string = '|S%d' % len(dummy)
    ret = np.zeros(jd_count,dtype=type_string)    
    for ix in range(jd_count):
        try:
            ret[ix] = _jd2ymstring( jd[ix], ymstring )
        except:
            ret[ix] = np.nan
    return ret   

def mjd2ymstring( mjd, ymstring ):
    """
    The `mjd2ymstring` function enables array input.
    Documentation see the `_mjd2ymstring` function.
    """    
    mjd = np.array(mjd,ndmin=1)
    mjd_count = np.size(mjd)
    dummy = datetime.now().strftime(ymstring)
    type_string = '|S%d' % len(dummy)
    ret = np.zeros(mjd_count,dtype=type_string)    
    for ix in range(mjd_count):
        try:
            ret[ix] = _mjd2ymstring( mjd[ix], ymstring )
        except:
            ret[ix] = 'nan'
    return ret   
    
#.............................................................................#

def _is_leap( year, cal='auto' ):    
    """
    The `_is_leap` function determines whether the input year is a leap year or 
    not.
    
    Returns a boolean.    
    
    The result depends on the choice of the calender type:
    cal='julian': all years divisible by 4 are leap years.
    cal='gregorian': same rule except centurial years that are not divisible by
                     400.
    cal='auto': automatic choice of calender.
    
    Requires the `_is_julian` function.
    """    
    ij  = _is_julian( year, 1, 1 )
    c1a = ij and (cal=='auto' or cal=='julian')
    c1b = not ij and cal=='julian'
    c2a = not ij and (cal=='auto' or cal=='gregorian')
    c2b = ij and cal=='gregorian'
    ijl = year % 4 == 0
    igl = ijl and ((year % 100 != 0) or (year % 400 == 0))
    if c1a or c1b:
        il = np.copy(ijl)
    elif c2a or c2b:
        il = np.copy(igl)
    return il

def _is_julian( y, m, d ):      
    """
    The `_is_julian` function determines whether the input date is in the range 
    of the Julian date calender or not.
    
    Returns a boolean.
    """
    ret = False
    if y<1582:
        ret = True
    if y==1582 and m<10:
        ret = True
    if y==1582 and m==10 and d<=4:
        ret = True
    return ret
    
def _is_gregorian( y, m, d):
    """
    The `_is_gregorian` function determines whether the input date is in the 
    range of the Gregorian date calender or not.
    
    Returns a boolean.
    """
    ret = False
    if y>1582:
        ret = True
    if y==1582 and m>10:
        ret = True
    if y==1582 and m==10 and d>=15:
        ret = True
    return ret

def _is_valid( y, m, d, cal='auto' ):
    """
    The `_is_valid` function determines whether the input date is valid in the 
    following sense: the day after 1582-10-04 (Julian calendar) was followed by
    1582-10-15 (Gregorian calendar), any dates inbetween are invalid if the 
    input variable cal is set to 'gregorian'. Otherwise, Julian calendar is 
    assumed.
    
    Raises error if not valid.

    Requires the functions: `_is_julian`, `_is_gregorian`.
    
    """
    ij  = _is_julian( y, m, d )
    ig  = _is_gregorian( y, m, d )
    cig = (cal=='gregorian')
    if not (ij or ig) and cig:
        raise ValueError('the date %04d-%02d-%4.2f did not exist in the '
                         'Gregorian calendar.' % (y,m,d))        
    
def _ymd2jd( y, m, d, cal='auto' ): 
    """
    The `_ymd2jd` function converts a calender date into a Julian date 
    (serial day number).
    
    Returns a float.
    
    Requires the `_is_julian` function.
    
    Reference: Meeus 1991 (ISBN:0943396352)
    """
    _is_valid( y, m, d, cal=cal )
    if m <= 2:
        a = y - 1.0
        b = m + 12.0
    else:
        a = np.copy( y )
        b = np.copy( m )
    ij  = _is_julian( y, m, d )
    c1a = ij and (cal=='auto' or cal=='julian')
    c1b = not ij and cal=='julian'
    c2a = not ij and (cal=='auto' or cal=='gregorian')
    c2b = ij and cal=='gregorian'
    if c1a or c1b:
        B = 0
    elif c2a or c2b:
        A = np.floor( a/100.0 )
        B = 2.0 - A + np.floor( A/4.0 )
    JD = (np.floor(365.25 * (a+4716.0)) 
          + np.floor(30.6001*(b+1.0)) 
          + d 
          + B 
          - 1524.5)
    return JD
    
def _jd2ymd( jd, cal='auto' ):
    """
    The `_jd2ymd` function converts a Julian date (serial day number) into a 
    calender date (year, month, day).
    
    Returns three values: int (year), int (month), float (day).
        
    Reference: Meeus 1991 (ISBN:0943396352)
    """
    J = jd + 0.5
    Z = np.floor( J )
    F = J - Z
    jg  = 2299161
    c1a = (Z < jg)  and (cal=='auto' or cal=='julian')
    c1b = (Z >= jg) and cal=='julian'
    c2a = (Z >= jg) and (cal=='auto' or cal=='gregorian')
    c2b = (Z < jg)  and cal=='gregorian'
    if c1a or c1b:
        A = np.copy(Z)    
    elif c2a or c2b:
        a = np.floor((Z - 1867216.25)/36524.25)
        A = Z + 1.0 + a - np.floor(a/4.0)
    B = A + 1524.0
    C = np.floor((B-122.1)/365.25)
    D = np.floor(365.25*C)
    E = np.floor((B-D)/30.6001)
    day   = B - D - np.floor(30.6001*E) + F
    if E<14:
        month = E - 1
    else:
        month = E - 13
    if month>2:
        year = C - 4716
    else:
        year = C - 4715
    y = np.int_(year)
    m = np.int_(month)
    d = day
    _is_valid( y, m, d, cal=cal )
    return y, m, d

def _ymd2mjd( y, m, d, cal='auto' ):
    """
    The `_ymd2mjd` function converts a calender date into a Modified Julian date 
    (serial day number).
    
    Returns a float.
    
    Requires the functions: `_jd2mjd`, `_ymd2jd`, `_is_julian`.
    
    Reference: Meeus 1991 (ISBN:0943396352)
    """
    return _jd2mjd(_ymd2jd(y, m, d, cal))
    
def _mjd2ymd( mjd, cal='auto' ):
    """
    The `_mjd2ymd` function converts a Modified Julian date (serial day number) 
    into a calender date (year, month, day).
    
    Returns three values: int (year), int (month), float (day).
        
    Requires the functions: `_mjd2jd`, `_jd2ymd`.
    
    Reference: Meeus 1991 (ISBN:0943396352)
    """
    return _jd2ymd(_mjd2jd(mjd),cal)
    
def _jd2dow( jd ):
    """
    The `_jd2dow` function determines the day of the week for a given a Julian 
    date (serial day number).
    
    Returns two values: int (day number), str (day name).
    """
    # JD=0 was a Monday (if it would have existed already)
    days_nam = ['Sun','Mon','Tue','Wed','Thu','Fri','Sat'] 
    days_num = [7,1,2,3,4,5,6] 
    J = np.floor( jd )
    A = jd - J
    if A >= 0.5:
        B = J + 2.0
    else:
        B = J + 1.0
    n = np.int_(B % 7)
    return days_num[n], days_nam[n]    
        
def _mjd2dow( mjd ): 
    """
    The `_mjd2dow` function determines the day of the week for a given a 
    Modified Julian date (serial day number).
    
    Returns two values: int (day number), str (day name).
    
    Requires the `_mjd2jd` function.
    """ 
    # November 17, 1858 was a Wednesday
    return _jd2dow( _mjd2jd( mjd ) )

def _dhms2day( d, h, m, s ):
    """
    The `_dhms2day` function converts hour, minute and second format for a given 
    day (decimal or integer) into a decimal day.
    
    Returns a float.
    """ 
    return d + h/24.0 + m/1440.0 + s/86400.0

def _day2dhms( day ):
    """
    The `_day2dhms` function converts a given decimal day into the day, hour, 
    minute, second format.
    
    Returns four values: int (day), int (hour), int (minute), float (second).
    """ 
    d = np.floor( day )
    h = np.floor( (day-d) * 24.0 )
    m = np.floor( ((day-d) * 24.0 - h ) * 60.0 )
    s = ((((day-d) * 24.0 - h ) * 60.0 ) - m ) * 60.0
    return np.int_(d), np.int_(h), np.int_(m), s

def _ymd2doy( y, m, d ):
    """
    The `_ymd2doy` function determines the day of the year for a given date.
    
    Returns an integer.
    
    Requires the functions: `_is_leap`, `_is_julian`.
    """ 
    months = [31, 28, 31, 30, 31, 30, 
              31, 31, 30, 31, 30, 31]
    if _is_leap( y ):
        months[1] = 29
    return np.sum( months[0:m-1] ) + d

def _jd2doy( jd ):
    """
    The `_jd2doy` function determines the day of the year for a given Julian 
    date (serial day number).
    
    Returns an integer.
    
    Requires the functions: `_jd2ymd`, `_ymd2doy`, `_is_leap`, `_is_julian`.
    """ 
    y, m, d = _jd2ymd( jd )
    return _ymd2doy( y, m, d )

def _mjd2doy( mjd ):
    """
    The `_mjd2doy` function determines the day of the year for a given Modified
    Julian date (serial day number).
    
    Returns an integer.
    
    Requires the functions: `_mjd2ymd`, `_mjd2jd`, `_jd2ymd`, `_ymd2doy`, 
                            `_is_leap`, `_is_julian`.
    """ 
    y, m, d = _mjd2ymd( mjd )
    return _ymd2doy( y, m, d )

def _doy2ymd( y, doy ):
    """
    The `_doy2ymd` function determines the date (year, month, day) for a given 
    day of a given year.
    
    Returns three values: int (year), int (month), float (day).
    
    Requires the functions: `_is_leap`, `_is_julian`.
    """ 
    months = [31, 28, 31, 30, 31, 30, 
              31, 31, 30, 31, 30, 31]
    if _is_leap( y ):
        months[1] = 29
    months_cs = np.cumsum( months )
    m = 1
    while (months_cs[m-1] - doy) < 0:
        m += 1
    d = months[m-1] - (months_cs[m-1] - doy)
    return y, m, d

def _doy2jd( y, doy ):
    """
    The `_doy2jd` function determines the Julian date (serial day number) for a 
    given day of a given year.
    
    Returns a float.
    
    Requires the functions: `_doy2ymd`, `_ymd2jd`, `_is_leap`, `_is_julian`.
    """ 
    y, m, d = _doy2ymd( y, doy )
    return _ymd2jd( y, m, d )

def _doy2mjd( y, doy ):
    """
    The `_doy2mjd` function determines the Modified Julian date (serial day 
    number) for a given day of a given year.
    
    Returns a float.
    
    Requires the functions: `_doy2ymd`, `_ymd2mjd`, `_is_leap`, `_is_julian`, 
                            `_jd2mjd`, `_ymd2jd`.
    """ 
    y, m, d = _doy2ymd( y, doy )
    return _ymd2mjd( y, m, d )

def _mjd2jd( mjd ):
    """
    The `_mjd2jd` function converts the Modified Julian date (serial day number)
    into a Julian date (serial day number).
    
    Returns a float.
    """ 
    return mjd + 2400000.5
    
def _jd2mjd( jd ):
    """
    The `_jd2mjd` function converts the Julian date (serial day number)into a 
    Modified Julian date (serial day number).
    
    Returns a float.
    """ 
    return jd - 2400000.5  

def _mjd2datetimenumber( mjd ):
    """
    The `_mjd2datetimenumber` function converts the Modified Julian date (serial 
    day number) into the Python datetime format.
    
    Returns a datetime.datetime object.
    
    Requires the functions: `_mjd2ymd`, `_day2dhms`, `_mjd2jd`, `_jd2ymd`.
    """ 
    y,m,d = _mjd2ymd( mjd )
    _,h,mins,s = _day2dhms( d )    
    ms = (s-np.floor(s)) * 1e6
    return datetime(y,m,int(d),h,mins,int(np.floor(s)),int(ms))

def _datetimenumber2mjd( t ):
    """
    The `_datetimenumber2mjd` function converts the Python datetime format into 
    a Modified Julian date (serial day number).
    
    Returns a float.
    
    Requires the functions: `_ymd2mjd`, _jd2mjd`, `_ymd2jd`, `_is_julian`.
    """ 
    y = t.year
    m = t.month
    d = (t.day
        +t.hour/24.
        +t.minute/3600.
        +t.second/86400.
        +t.microsecond/86400000000.)
    return _ymd2mjd(y,m,d)
    
def _ym2jdmrange( y, m ):
    """
    The `_ym2jdmrange` function determines the first and last day of a given 
    month in a given year and converts these two dates into the Julian Date 
    (serial day number).
    
    Returns two values: float (first day), float (last day).
    
    Requires the functions: `_ymd2jd`, `_is_julian`.
    """ 
    ti = _ymd2jd( y, m, 1 )
    if m==12:
        tk = _ymd2jd( y+1, 1, 1 )
    else:
        tk = _ymd2jd( y, m+1, 1 )  
    tj = tk - 1
    return ti, tj    
    
def _ym2mjdmrange( y, m ):
    """
    The `_ym2mjdmrange` function determines the first and last day of a given 
    month in a given year and converts these two dates into the Modified Julian
    Date (serial day number).
    
    Returns two values: float (first day), float (last day).
    
    Requires the functions: `_ymd2mjd`,`_jd2mjd`, `_ymd2jd`, `_is_julian`.
    """ 
    ti = _ymd2mjd( y, m, 1 )
    if m==12:
        tk = _ymd2mjd( y+1, 1, 1 )
    else:
        tk = _ymd2mjd( y, m+1, 1 )  
    tj = tk - 1
    return ti, tj    
        
def _jd2jdmrange( jd ):
    """
    The `_jd2jdmrange` function determines the first and last day of the month
    from a given Julian date and returns these dates in the same format.
    
    Returns two values: float (first day), float (last day).
    
    Requires the functions: `_ym2jdmrange`, `_ymd2jd`, `_is_julian`, `_jd2ymd`.
    """ 
    y,m,d = _jd2ymd( jd )
    return _ym2jdmrange( y, m )
    
def _mjd2mjdmrange( mjd ):
    """
    The `_mjd2mjdmrange` function determines the first and last day of the month
    from a given Modified Julian date and returns these dates in the same 
    format.
    
    Returns two values: float (first day), float (last day).
    
    Requires the functions: `_ym2mjdmrange`, `_ymd2mjd`,`_jd2mjd`, `_ymd2jd`, 
                            `_is_julian`, `_mjd2jd`, `_jd2ymd`.
    """ 
    y,m,d = _mjd2ymd( mjd )
    return _ym2mjdmrange( y, m )
    
def _ymstring2jd( ymstr ):
    """
    The `_ymstring2jd` function determines the date from a date string given in 
    the format 'yyyy-mm' or 'yyyy-mm-dd' and converts this date into the 
    Julian date format.
    
    Returns a float.
    
    Requires the functions: `_ymd2jd`, `_is_julian`, `_dhms2day`.
    """ 
    if len(ymstr)==10:
        ymd = [int(ymstr[0:4]),int(ymstr[5:7]),int(ymstr[8:10])]        
    elif len(ymstr)==7:
        ymd = [int(ymstr[0:4]),int(ymstr[5:7]),1]
    elif len(ymstr)==19:
        ymd = [int(ymstr[0:4]),int(ymstr[5:7]),int(ymstr[8:10])] 
        hms = [int(ymstr[11:13]),int(ymstr[14:16]),int(ymstr[17:19])] 
        ymd[2] = _dhms2day(ymd[2],hms[0],hms[1],hms[2])
    else:
        raise ValueError('unknown format')
    return _ymd2jd(ymd[0],ymd[1],ymd[2])    
    
def _ymstring2mjd( ymstr ):
    """
    The `_ymstring2mjd` function determines the date from a date string given in 
    the format 'yyyy-mm' or 'yyyy-mm-dd' and converts this date into the 
    Modified Julian date format.
    
    Returns a float.
    
    Requires the functions: `_ymd2mjd`, `_jd2mjd`, `_ymd2jd`, `_is_julian`,
                            `_dhms2day`.
    """ 
    if len(ymstr)==10:
        ymd = [int(ymstr[0:4]),int(ymstr[5:7]),int(ymstr[8:10])]        
    elif len(ymstr)==7:
        ymd = [int(ymstr[0:4]),int(ymstr[5:7]),1]
    elif len(ymstr)==19:
        ymd = [int(ymstr[0:4]),int(ymstr[5:7]),int(ymstr[8:10])] 
        hms = [int(ymstr[11:13]),int(ymstr[14:16]),int(ymstr[17:19])] 
        ymd[2] = _dhms2day(ymd[2],hms[0],hms[1],hms[2])
    else:
        raise ValueError('unknown format')
    return _ymd2mjd(ymd[0],ymd[1],ymd[2])    
    
def _jd2ymstring( jd, ymstring ):
    """
    The `_jd2ymstring` function converts a date given in the Julian date format 
    into a string specified by the ymstring format flag (Python's datetime 
    notation).
    
    Returns a string.
    
    Requires the functions: `_mjd2datetimenumber`, `_jd2mjd`, `_mjd2ymd`, 
                            `_day2dhms`, `_mjd2jd`, `_jd2ymd`.
    """ 
    t = _mjd2datetimenumber(_jd2mjd(jd))
    return t.strftime(ymstring)
    
def _mjd2ymstring( mjd, ymstring ):
    """
    The `_mjd2ymstring` function converts a date given in the Modified Julian 
    date format into a string specified by the ymstring format flag (Python's 
    datetime notation).
    
    Returns a string.
    
    Requires the functions: `_mjd2datetimenumber`, `_mjd2ymd`, 
                            `_day2dhms`, `_mjd2jd`, `_jd2ymd`.
    """ 
    t = _mjd2datetimenumber(mjd)
    return t.strftime(ymstring)
    
    
    
    
    
    
    
    
    