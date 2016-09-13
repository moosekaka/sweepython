# -*- coding: utf-8 -*-
"""
Created on Fri Sep 09 16:44:59 2016

@author: sweel_Rafelski
"""
import io
import requests
import bs4

#with open('text.html', 'w') as f:
with io.open('text3.html', 'w', encoding='utf8') as f:
    for i in range(113, 168):

        url = 'http://www.slow-chinese.com/podcasts/{}'.format(i)
        print "Scrapping {}".format(url)
        r = requests.get(url)
        soup = bs4.BeautifulSoup(r.content.decode('utf8'))
        #print title
        for i in soup.findAll('h1'):
            f.write(''.join(i.findAll(text=True)))
            f.write(u'<br>')
        #print body
        temp = soup.findAll('div', attrs={'id': '-0'})
        if len(temp):
            for node in temp:
                f.write(''.join(node.findAll(text=True)))
                f.write(u'<br><br>')
        else:
            for node in soup.findAll('article'):
                for subnode in node.findAll('p', attrs={'class': None}):
                    f.write(''.join(subnode.findAll(text=True)))
                    f.write(u'<br>')