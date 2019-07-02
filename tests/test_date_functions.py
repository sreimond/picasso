# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.dates_and_time import date_functions
import numpy as np
from datetime import datetime


def test_leap_years_1899to1904_and_1999to2004():
	"""
	Test if the function `is_leap` correctly determines whether the years 
	around the turns of the last two centuries are leap years or not.        
	"""
	years = [1899,1900,1901,1902,1903,1904,1999,2000,2001,2002,2003,2004]
	leaps = date_functions.is_leap( years )
	assert not leaps[0]  # 1899 was not a leap year
	assert not leaps[1]  # 1900 was not a leap year
	assert not leaps[2]  # 1901 was not a leap year
	assert not leaps[3]  # 1902 was not a leap year
	assert not leaps[4]  # 1903 was not a leap year
	assert leaps[5]      # 1904 was a leap year
	assert not leaps[6]  # 1999 was not a leap year
	assert leaps[7]      # 2000 was a leap year
	assert not leaps[8]  # 2001 was not a leap year
	assert not leaps[9]  # 2002 was not a leap year
	assert not leaps[10] # 2003 was not a leap year
	assert leaps[11]     # 2004 was a leap year
    
    
def test_julian_dates_10_1000_1582_2000():
	"""
	Test if the function `is_leap` correctly determines whether the years 
	around the turns of the last two centuries are leap years or not.        
	"""
	years = [10,1000,1582,1582,2000]
	months = [5,12,10,10,1]
	days = [30,20,4,5,1]
	julians = date_functions.is_julian( years, months, days )
	assert julians[0] # 10-05-30 was a Julian date
	assert julians[1] # 1000-12-20 was a Julian date
	assert julians[2] # 1582-10-04 was a Julian date
	assert not julians[3] # 1582-10-05 was not a Julian date
	assert not julians[4] # 2000-01-01 was not a Julian date
        
def test_gregorian_dates_10_1000_1582_2000():
	"""
	Test if the function `is_leap` correctly determines whether the years 
	around the turns of the last two centuries are leap years or not.        
	"""
	years = [10,1000,1582,1582,2000]
	months = [5,12,10,10,1]
	days = [30,20,4,5,1]
	gregorians = date_functions.is_gregorian( years, months, days )
	assert not gregorians[0] # 10-05-30 was not a Gegorian date
	assert not gregorians[1] # 1000-12-20 was not a Gregorian date
	assert not gregorians[2] # 1582-10-04 was not a Gregorian date
	assert not gregorians[3] # 1582-10-05 was not a Gregorian date
	assert gregorians[4] # 2000-01-01 was a Gregorian date
    
def test_ymd2jd():
	"""
	Test if the function `ymd2jd` correctly converts dates into Julian date
	day numbers or not.        
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	years = [2000,1999,1987,1987,1988,1988,1900,1600,1600,837,
			 -123,-122,-1000,-1000,-1001,-4712]
	months = [1,1,1,6,1,6,1,1,12,4,12,1,7,2,8,1]
	days = [1.5,1.0,27.0,19.5,27.0,19.5,1.0,1.0,31.0,10.3,31.0,1.0,12.5,
			29.0,17.9,1.5]
	jd_true = [2451545.0,2451179.5,2446822.5,2446966.0,2447187.5,2447332.0,
			   2415020.5,2305447.5,2305812.5,2026871.8,1676496.5,1676497.5,
			   1356001.0,1355866.5,1355671.4,0.0]
	jd_test = date_functions.ymd2jd( years, months, days )
	assert jd_test[0] == pytest.approx(jd_true[0],1e-12)
	assert jd_test[1] == pytest.approx(jd_true[1],1e-12)
	assert jd_test[2] == pytest.approx(jd_true[2],1e-12)
	assert jd_test[3] == pytest.approx(jd_true[3],1e-12)
	assert jd_test[4] == pytest.approx(jd_true[4],1e-12)
	assert jd_test[5] == pytest.approx(jd_true[5],1e-12)
	assert jd_test[6] == pytest.approx(jd_true[6],1e-12)
	assert jd_test[7] == pytest.approx(jd_true[7],1e-12)
	assert jd_test[8] == pytest.approx(jd_true[8],1e-12)
	assert jd_test[9] == pytest.approx(jd_true[9],1e-12)
	assert jd_test[10] == pytest.approx(jd_true[10],1e-12)
	assert jd_test[11] == pytest.approx(jd_true[11],1e-12)
	assert jd_test[12] == pytest.approx(jd_true[12],1e-12)
	assert jd_test[13] == pytest.approx(jd_true[13],1e-12)
	assert jd_test[14] == pytest.approx(jd_true[14],1e-12)
	assert jd_test[15] == pytest.approx(jd_true[15],1e-12)        

def test_jd2ymd():
	"""
	Test if the function `jd2ymd` correctly converts Julian dates into date
	(years, months, days) or not.        
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	years_true = [2000,1999,1987,1987,1988,1988,1900,1600,1600,837,
			 -123,-122,-1000,-1000,-1001,-4712,1957,333,-584]
	months_true = [1,1,1,6,1,6,1,1,12,4,12,1,7,2,8,1,10,1,5]
	days_true = [1.5,1.0,27.0,19.5,27.0,19.5,1.0,1.0,31.0,10.3,31.0,1.0,
				 12.5,29.0,17.9,1.5,4.81,27.5,28.63]
	jd = [2451545.0,2451179.5,2446822.5,2446966.0,2447187.5,2447332.0,
		  2415020.5,2305447.5,2305812.5,2026871.8,1676496.5,1676497.5,
		  1356001.0,1355866.5,1355671.4,0.0,2436116.31,1842713.0,
		  1507900.13]        
	years_test, months_test, days_test = date_functions.jd2ymd( jd )
	assert years_test[0] == years_true[0]
	assert years_test[1] == years_true[1]
	assert years_test[2] == years_true[2]
	assert years_test[3] == years_true[3]
	assert years_test[4] == years_true[4]
	assert years_test[5] == years_true[5]
	assert years_test[6] == years_true[6]
	assert years_test[7] == years_true[7]
	assert years_test[8] == years_true[8]
	assert years_test[9] == years_true[9]
	assert years_test[10] == years_true[10]
	assert years_test[11] == years_true[11]
	assert years_test[12] == years_true[12]
	assert years_test[13] == years_true[13]
	assert years_test[14] == years_true[14]
	assert years_test[15] == years_true[15]
	assert years_test[16] == years_true[16]
	assert years_test[17] == years_true[17]
	assert years_test[18] == years_true[18]
	
	assert months_test[0] == months_true[0]
	assert months_test[1] == months_true[1]
	assert months_test[2] == months_true[2]
	assert months_test[3] == months_true[3]
	assert months_test[4] == months_true[4]
	assert months_test[5] == months_true[5]
	assert months_test[6] == months_true[6]
	assert months_test[7] == months_true[7]
	assert months_test[8] == months_true[8]
	assert months_test[9] == months_true[9]
	assert months_test[10] == months_true[10]
	assert months_test[11] == months_true[11]
	assert months_test[12] == months_true[12]
	assert months_test[13] == months_true[13]
	assert months_test[14] == months_true[14]
	assert months_test[15] == months_true[15]
	assert months_test[16] == months_true[16]
	assert months_test[17] == months_true[17]
	assert months_test[18] == months_true[18]
	
	assert days_test[0] == pytest.approx(days_true[0],1e-5)
	assert days_test[1] == pytest.approx(days_true[1],1e-5)
	assert days_test[2] == pytest.approx(days_true[2],1e-5)
	assert days_test[3] == pytest.approx(days_true[3],1e-5)
	assert days_test[4] == pytest.approx(days_true[4],1e-5)
	assert days_test[5] == pytest.approx(days_true[5],1e-5)
	assert days_test[6] == pytest.approx(days_true[6],1e-5)
	assert days_test[7] == pytest.approx(days_true[7],1e-5)
	assert days_test[8] == pytest.approx(days_true[8],1e-5)
	assert days_test[9] == pytest.approx(days_true[9],1e-5)
	assert days_test[10] == pytest.approx(days_true[10],1e-5)
	assert days_test[11] == pytest.approx(days_true[11],1e-5)
	assert days_test[12] == pytest.approx(days_true[12],1e-5)
	assert days_test[13] == pytest.approx(days_true[13],1e-5)
	assert days_test[14] == pytest.approx(days_true[14],1e-5)
	assert days_test[15] == pytest.approx(days_true[15],1e-5)
	assert days_test[16] == pytest.approx(days_true[16],1e-5)
	assert days_test[17] == pytest.approx(days_true[17],1e-5)
	assert days_test[18] == pytest.approx(days_true[18],1e-5)
                
def test_ymd2mjd():
	"""
	Test if the function `ymd2jd` correctly converts dates into Modified 
	Julian date day numbers or not.        
	
	Reference: Meeus 1991 (ISBN:0943396352)
			   https://bowie.gsfc.nasa.gov/time/
			   https://de.mathworks.com/help/aerotbx/ug/mjuliandate.html
	"""
	years = [1858,1980,2006]
	months = [11,1,12]
	days = [17.0,30.0,19.0]
	mjd_true = [0.0,44268.0,54088.0]
	mjd_test = date_functions.ymd2mjd( years, months, days )
	assert mjd_test[0] == pytest.approx(mjd_true[0],1e-12)
	assert mjd_test[1] == pytest.approx(mjd_true[1],1e-12)
	assert mjd_test[2] == pytest.approx(mjd_true[2],1e-12)
        
def test_mjd2ymd():
	"""
	Test if the function `mjd2ymd` correctly converts Modified Julian 
	dates into date (years, months, days) or not.        
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	years_true = [1858,1980,2006]
	months_true = [11,1,12]
	days_true = [17.0,30.0,19.0]
	mjd = [0.0,44268.0,54088.0]
	years_test, months_test, days_test = date_functions.mjd2ymd( mjd )
	assert years_test[0] == years_true[0]
	assert years_test[1] == years_true[1]
	assert years_test[2] == years_true[2]
	
	assert months_test[0] == months_true[0]
	assert months_test[1] == months_true[1]
	assert months_test[2] == months_true[2]
	
	assert days_test[0] == pytest.approx(days_true[0],1e-5)
	assert days_test[1] == pytest.approx(days_true[1],1e-5)
	assert days_test[2] == pytest.approx(days_true[2],1e-5)
        
def test_jd2dow():
	"""
	Test if the function `jd2dow` correctly determines the day of the week
	from a given Julian date.        
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	jd = [2434923.5,2458130.5]
	dnum_true = [3,5]
	dnam_true = np.array(['Wed','Fri'],dtype='|S3')
	dnum_test, dnam_test = date_functions.jd2dow( jd )
	
	assert dnum_test[0] == dnum_true[0]
	assert dnum_test[1] == dnum_true[1]
	assert dnam_test[0] == dnam_true[0]
	assert dnam_test[1] == dnam_true[1]
        
def test_mjd2dow():
	"""
	Test if the function `mjd2dow` correctly determines the day of the week
	from a given Modified Julian date.        
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	mjd = [34923.0,58130.5]
	dnum_true = [3,5]
	dnam_true = np.array(['Wed','Fri'],dtype='|S3')
	dnum_test, dnam_test = date_functions.mjd2dow( mjd )
	
	assert dnum_test[0] == dnum_true[0]
	assert dnum_test[1] == dnum_true[1]
	assert dnam_test[0] == dnam_true[0]
	assert dnam_test[1] == dnam_true[1]
            
def test_dhms2day():
	"""
	Test if the function `dhms2day` correctly converts from the day, hours,
	minutes, seconds format to the decimal day format. 
	
	https://en.wikipedia.org/wiki/Decimal_time
	MATLAB
	"""
	d = [0,737072]
	h = [7,16]
	m = [12,43]
	s = [0.0,25.12]
	d_true = [0.3,737072.696818518569]
	d_test = date_functions.dhms2day(d,h,m,s)
	
	assert d_test[0] == pytest.approx(d_true[0],1e-8)
	assert d_test[1] == pytest.approx(d_true[1],1e-8)
        
def test_day2dhms():
	"""
	Test if the function `day2dhms` correctly converts from the the decimal
	day format to the day, hours, minutes, seconds format. 
	
	https://en.wikipedia.org/wiki/Decimal_time
	MATLAB
	"""
	d = [0.300520833333,737072.696818518569]
	d_true = [0,737072]
	h_true = [7,16]
	m_true = [12,43]
	s_true = [45.0,25.12]
	
	d_test, h_test, m_test, s_test = date_functions.day2dhms( d )
	
	assert d_test[0] == d_true[0]
	assert d_test[1] == d_true[1]
	
	assert h_test[0] == h_true[0]
	assert h_test[1] == h_true[1]
	
	assert m_test[0] == m_true[0]
	assert m_test[1] == m_true[1]
	
	assert s_test[0] == pytest.approx(s_true[0],1e-5)
	assert s_test[1] == pytest.approx(s_true[1],1e-5)
	  
def test_ymd2doy():
	"""
	Test if the function `ymd2doy` correctly converts from the year, month,
	day format to the day of the year. 
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	years = [1978,1988]
	months = [11,4]
	days = [14,22]
	dn_true = [318,113]
	dn_test = date_functions.ymd2doy(years,months,days)
	
	assert dn_test[0] == dn_true[0]
	assert dn_test[1] == dn_true[1]

def test_jd2doy():
	"""
	Test if the function `jd2doy` correctly converts from the Julian date 
	format to the day of the year. 
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	jd = [2443826.5,2447273.5]
	dn_true = [318,113]
	dn_test = date_functions.jd2doy( jd )
	
	assert dn_test[0] == dn_true[0]
	assert dn_test[1] == dn_true[1]  
        
def test_mjd2doy():
	"""
	Test if the function `mjd2doy` correctly converts from the Modified 
	Julian date format to the day of the year. 
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	mjd = [43826.0,47273.0]
	dn_true = [318,113]
	dn_test = date_functions.mjd2doy( mjd )
	
	assert dn_test[0] == dn_true[0]
	assert dn_test[1] == dn_true[1]
	
def test_doy2ymd():
	"""
	Test if the function `doy2ymd` correctly converts from the day of the 
	year format to the year, month, day format.
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	dn = [318,113]
	years_true = [1978,1988]
	months_true = [11,4]
	days_true = [14,22]
	years_test, months_test, days_test = date_functions.doy2ymd(years_true,
																dn)        

	assert years_test[0] == years_true[0]
	assert years_test[1] == years_true[1]
	
	assert months_test[0] == months_true[0]
	assert months_test[1] == months_true[1]
	
	assert days_test[0] == days_true[0]
	assert days_test[1] == days_true[1]        
	
def test_doy2jd():
	"""
	Test if the function `doy2jd` correctly converts from the day of the 
	year format to the Julian date format.
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	dn = [318,113]
	years_true = [1978,1988]
	jd_true = [2443826.5,2447273.5]
	jd_test = date_functions.doy2jd(years_true,dn)

	assert jd_test[0] == jd_true[0]
	assert jd_test[1] == jd_true[1]
        
def test_doy2mjd():
	"""
	Test if the function `doy2mjd` correctly converts from the day of the 
	year format to the Modified Julian date format.
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	dn = [318,113]
	years_true = [1978,1988]
	mjd_true = [43826,47273]
	mjd_test = date_functions.doy2mjd(years_true,dn)

	assert mjd_test[0] == mjd_true[0]
	assert mjd_test[1] == mjd_true[1]
        
def test_mjd2jd():
	"""
	Test if the function `mjd2jd` correctly converts from the Modified 
	Julian date format to the Julian date format.
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	mjd = [43826,47273]
	jd_true = [2443826.5,2447273.5]
	jd_test = date_functions.mjd2jd( mjd )

	assert jd_test[0] == jd_true[0]
	assert jd_test[1] == jd_true[1]
	
def test_jd2mjd():
	"""
	Test if the function `mjd2jd` correctly converts from the Julian date 
	format to the Modified Julian date format.
	
	Reference: Meeus 1991 (ISBN:0943396352)
	"""
	
	jd = [2443826.5,2447273.5]
	mjd_true = [43826,47273]
	mjd_test = date_functions.jd2mjd( jd )

	assert mjd_test[0] == mjd_true[0]
	assert mjd_test[1] == mjd_true[1]
    
    
def test_mjd2datetimenumber():
	"""
	Test if the function `mjd2datetimenumber` correctly converts from the 
	Modified Julian date format to the Python datetime format.        
	"""
	
	mjd = [43826,47273]
	dt_true = [datetime(1978,11,14),datetime(1988,4,22)]
	dt_test = date_functions.mjd2datetimenumber( mjd )

	assert dt_test[0] == dt_true[0]
	assert dt_test[1] == dt_true[1]
        
def test_datetimenumber2mjd():
	"""
	Test if the function `datetimenumber2mjd` correctly converts from the 
	Python datetime format to the Modified Julian date format.
	"""
	dt = [datetime(1978,11,14),datetime(1988,4,22)]
	mjd_true = [43826,47273]
	mjd_test = date_functions.datetimenumber2mjd( dt )

	assert mjd_test[0] == mjd_true[0]
	assert mjd_test[1] == mjd_true[1]
        
def test_ym2jdmrange():
	"""
	Test if the function `ym2jdmrange` correctly determines the first and
	last day of the given month from the year, month format in the Julian 
	date format.
	"""
	years = [1899,2000,2020]
	months = [12,1,2]
	jd1_true = [2414989.5,2451544.5,2458880.5]
	jd2_true = [2415019.5,2451574.5,2458908.5]
	jd1_test, jd2_test = date_functions.ym2jdmrange( years, months )

	assert jd1_test[0] == jd1_true[0]
	assert jd1_test[1] == jd1_true[1]
	assert jd1_test[2] == jd1_true[2]
	assert jd2_test[0] == jd2_true[0]
	assert jd2_test[1] == jd2_true[1]
	assert jd2_test[2] == jd2_true[2]


def test_ym2mjdmrange():
	"""
	Test if the function `ym2mjdmrange` correctly determines the first and
	last day of the given month from the year, month format in the Modified
	Julian date format.
	"""
	years = [1899,2000,2020]
	months = [12,1,2]
	mjd1_true = [14989,51544,58880]
	mjd2_true = [15019,51574,58908]
	mjd1_test, mjd2_test = date_functions.ym2mjdmrange( years, months )

	assert mjd1_test[0] == mjd1_true[0]
	assert mjd1_test[1] == mjd1_true[1]
	assert mjd1_test[2] == mjd1_true[2]
	assert mjd2_test[0] == mjd2_true[0]
	assert mjd2_test[1] == mjd2_true[1]
	assert mjd2_test[2] == mjd2_true[2]
        
def test_jd2jdmrange():
	"""
	Test if the function `jd2jdmrange` correctly determines the first and
	last day of the given month from the Julian Date format in the Julian 
	date format.
	"""
	jd = [2414999.4,2451573.1,2458907.9]
	jd1_true = [2414989.5,2451544.5,2458880.5]
	jd2_true = [2415019.5,2451574.5,2458908.5]
	jd1_test, jd2_test = date_functions.jd2jdmrange( jd )

	assert jd1_test[0] == jd1_true[0]
	assert jd1_test[1] == jd1_true[1]
	assert jd1_test[2] == jd1_true[2]
	assert jd2_test[0] == jd2_true[0]
	assert jd2_test[1] == jd2_true[1]
	assert jd2_test[2] == jd2_true[2]
        
def test_mjd2mjdmrange():
	"""
	Test if the function `mjd2mjdmrange` correctly determines the first and
	last day of the given month from the Modifed Julian Date format in the 
	Modifed Julian date format.
	"""
	mjd = [14999.1,51560.6,58907.9]
	mjd1_true = [14989,51544,58880]
	mjd2_true = [15019,51574,58908]
	mjd1_test, mjd2_test = date_functions.mjd2mjdmrange( mjd )

	assert mjd1_test[0] == mjd1_true[0]
	assert mjd1_test[1] == mjd1_true[1]
	assert mjd1_test[2] == mjd1_true[2]
	assert mjd2_test[0] == mjd2_true[0]
	assert mjd2_test[1] == mjd2_true[1]
	assert mjd2_test[2] == mjd2_true[2]

def test_ymstring2jd():
	"""
	Test if the function `ymstring2jd` correctly converts the date given in
	the "yyyy-mm" or "yyyy-mm-dd" or "yyyy-mm-dd HH:MM:SS" string format to 
	the Julian date format.
	"""
	ymstrings = ['1952-05','2018-01-16 09:54:22','2000-12-15']
	jd_true = [2434133.5,2458134.9127546297,2451893.5]
	
	jd_test = date_functions.ymstring2jd( ymstrings )

	assert jd_test[0] == jd_true[0]
	assert jd_test[1] == pytest.approx(jd_true[1],1e-8)
	assert jd_test[2] == jd_true[2]
        
def test_ymstring2mjd():
	"""
	Test if the function `ymstring2mjd` correctly converts the date given in
	the "yyyy-mm" or "yyyy-mm-dd" or "yyyy-mm-dd HH:MM:SS" string format to 
	the Modified Julian date format.
	"""
	ymstrings = ['1952-05','2018-01-16 09:54:22','2000-12-15']
	mjd_true = [34133,58134.412754629629,51893]
	
	mjd_test = date_functions.ymstring2mjd( ymstrings )

	assert mjd_test[0] == mjd_true[0]
	assert mjd_test[1] == pytest.approx(mjd_true[1],1e-8)
	assert mjd_test[2] == mjd_true[2]

def test_jd2ymstring():
	"""
	Test if the function `jd2ymstring` correctly converts the date given in
	Julian date format to the Python date string format.
	"""
	jd = [2434133.5,2458134.9127546297,2451893.5]
	ymstrings_true = np.zeros(3,dtype='|S19') 
	ymstrings_true[0] = '1952-05-01 00:00:00'
	ymstrings_true[1] = '2018-01-16 09:54:22'
	ymstrings_true[2] = '2000-12-15 00:00:00'       
	ymstrings_test = date_functions.jd2ymstring( jd, '%Y-%m-%d %H:%M:%S' )

	assert ymstrings_test[0] == ymstrings_true[0]
	assert ymstrings_test[1] == ymstrings_true[1]
	assert ymstrings_test[2] == ymstrings_true[2]
        

def test_mjd2ymstring():
	"""
	Test if the function `mjd2ymstring` correctly converts the date given 
	in the Modfied Julian date format to the Python date string format.
	"""
	mjd = [34133,58134.412754629629,51893]
	ymstrings_true = np.zeros(3,dtype='|S19') 
	ymstrings_true[0] = '1952-05-01 00:00:00'
	ymstrings_true[1] = '2018-01-16 09:54:22'
	ymstrings_true[2] = '2000-12-15 00:00:00'           
	ymstrings_test = date_functions.mjd2ymstring(mjd,'%Y-%m-%d %H:%M:%S')

	assert ymstrings_test[0] == ymstrings_true[0]
	assert ymstrings_test[1] == ymstrings_true[1]
	assert ymstrings_test[2] == ymstrings_true[2]

