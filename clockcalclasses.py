# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 16:39:30 2016

@author: sweel_rafelski
"""

class Clock(object):
    #kwargs will hold all args not used in __init__ for Clock
    def __init__(self,hours=0, minutes=0, seconds=0, **kwargs):
        self.__hours = hours
        self.__minutes = minutes
        self.__seconds = seconds

    def set(self,hours, minutes, seconds=0):
        self.__hours = hours
        self.__minutes = minutes
        self.__seconds = seconds

    def tick(self):
        """ Time will be advanced by one second """
        if self.__seconds == 59:
            self.__seconds = 0
            if (self.__minutes == 59):
                self.__minutes = 0
                self.__hours = 0 if self.__hours==23  else self.__hours+1
	    else:
		self.__minutes += 1;
	else:
            self.__seconds += 1;

    def display(self):
        print("%d:%d:%d" % (self.__hours, self.__minutes, self.__seconds))

    def __str__(self):
        return "%2d:%2d:%2d" % (self.__hours, self.__minutes, self.__seconds)

#x = Clock()
#print(x)
#for i in xrange(10000):
#    x.tick()
#print(x)


#The Calendar Class

class Calendar(object):
    months = (31,28,31,30,31,30,31,31,30,31,30,31)
    #kwargs will hold all args not used in __init__ for Calendar
    def __init__(self, day=1, month=1, year=1900):  
        self.__day = day
        self.__month = month
        self.__year = year

    def leapyear(self,y):
	if y % 4:
	   # not a leap year
	   return 0;
	else:
	   if y % 100:
	     return 1;
           else:
	     if y % 400:
                return 0
	     else:
		return 1;

    def set(self, day, month, year):
        self.__day = day
        self.__month = month
        self.__year = year


    def get():
	return (self, self.__day, self.__month, self.__year)
    def advance(self):
        months = Calendar.months
	max_days = months[self.__month-1]
	if self.__month == 2:
		max_days += self.leapyear(self.__year)
	if self.__day == max_days:
		self.__day = 1
		if (self.__month == 12):
			self.__month = 1
			self.__year += 1
		else:
			self.__month += 1
	else:
		self.__day += 1


    def __str__(self):
       return str(self.__day)+"/"+ str(self.__month)+ "/"+ str(self.__year)

#if __name__ == "__main__":
#   x = Calendar()
#   print(x)
#   x.advance()
#   print(x)


#==============================================================================
# #The Calendar-Clock Class
# 
# from clock import Clock
# from calendar import Calendar
#==============================================================================

class CalendarClock(Clock, Calendar):

   def __init__(self, H=12, Y =3000, **kwargs):   
        #self.mem = x  #note that kwargs is needed here for argument unpacking
        Clock.__init__(self, hours=H, **kwargs)  #kwargs will store args for Calendar
        Calendar.__init__(self, year=Y, **kwargs)  #kwargs will store args for Clock

     

   def __str__(self):
       return Calendar.__str__(self) + ", " + Clock.__str__(self)


if __name__  == "__main__":
   x = CalendarClock(hours=1)
   print(x)
   #print x.mem
#   for i in range(1000):
#      x.tick()
#   for i in range(1000):
#      x.advance()
#   print(x)