# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 20:56:30 2016

@author: sweel_Rafelski
"""
import requests

for i in range(1, 168):
    url = 'http://www.slow-chinese.com/podcasts/Slow_Chinese_{}.mp3'.format(i)
    print "getting on {}".format(url)
    r = requests.get(url, stream=True)
    with open("Slow_Chinese_{:03d}.mp3".format(i), 'wb') as fd:
        for chunk in r.iter_content(1024):
            fd.write(chunk)
