#!/usr/bin/env python
# -*- coding: utf-8 -*-

from time import sleep
import sys
import bisect
import math
import collections

from scipy import stats
import pylab

import random
random.seed()


DEBUG = False

Statistics = collections.namedtuple('Statistics', ['initial_length',
                                                   'additions',
                                                   'deceased_transplants',
                                                   'living_transplants',
                                                   'removals'])

#CPMC
statistics = {
    'CPMC': {
    #   year             len   add  dtx  ltx  rem
        2014: Statistics(2047, 521, 141,  46, 431),
        2013: Statistics(2127, 573, 112,  55, 486),
        2012: Statistics(2177, 521, 143,  58, 370),
        2011: Statistics(2316, 492, 123,  66, 442),
        2010: Statistics(2311, 541, 109,  73, 349),
        2009: Statistics(2135, 634, 114,  57, 287),
        2008: Statistics(1975, 649, 111,  89, 284),
        2007: Statistics(1728, 704, 130,  54, 272),
        2006: Statistics(1611, 565, 129,  52, 267),
        2005: Statistics(1533, 509, 115,  53, 254),
        2004: Statistics(1475, 493, 107,  55, 273),
        2003: Statistics(1387, 444, 104,  38, 206),
        2002: Statistics(1312, 422, 102,  45, 200),
        2001: Statistics(1196, 478, 102,  40, 219),
        2000: Statistics(1048, 430,  78,  35, 169)
    },
    'UCLA': {
    #   year             len   add  dtx  ltx  rem
        2014: Statistics(2178, 574, 175, 126, 257),
        2013: Statistics(2172, 542, 166, 121, 249),
        2012: Statistics(2067, 615, 176, 107, 226),
        2011: Statistics(2080, 595, 212, 100, 297),
        2010: Statistics(2022, 633, 164, 132, 279),
        2009: Statistics(1973, 650, 135, 134, 332),
        2008: Statistics(2307, 463, 167, 127, 498),
        2007: Statistics(2511, 694, 186, 121, 591),
        2006: Statistics(2138, 944, 180, 104, 252),
        2005: Statistics(1900, 781, 170, 106, 267),
        2004: Statistics(1523, 849, 205, 102, 160),
        2003: Statistics(1304, 677, 208,  91, 159),
        2002: Statistics(1259, 460, 177,  79, 142),
        2001: Statistics(1148, 479, 164,  67, 107),
        2000: Statistics( 964, 541, 154,  95,  58),
        1999: Statistics( 809, 452, 164,  52,  81)
    },
    'UCD': {
    #   year             len   add  dtx  ltx  rem
        2014: Statistics(1066, 627, 257,  44, 8+40+25+30+4+47),
        2013: Statistics( 995, 574, 290,  47, 7+29+29+56+2+43),
        2012: Statistics(1045, 373, 222,  48, 12+32+32+41+1+35),
        2011: Statistics( 944, 519, 208,  53, 7+30+35+5+0+80),
        2010: Statistics( 846, 442, 126,  57, 9+27+35+3+0+87),
        2009: Statistics( 809, 314,  95,  45, 13+29+37+3+0+55),
        2008: Statistics( 714, 321,  79,  33, 9+21+28+5+0+51), 
        2007: Statistics( 624, 303,  73,  58, 10+21+26+3+0+22),
        2006: Statistics( 250, 511,  68,  23, 4+11+17+0+0+14),
        2005: Statistics( 282, 145,  44,  32, 61+10+15+1+0+14),
        2004: Statistics( 239, 150,  39,  20, 2+9+16+0+0+20),
        2003: Statistics( 198, 148,  47,  10, 2+7+19+0+0+22),
        2002: Statistics( 152, 116,  28,  10, 1+8+0+0+0+23),
        2001: Statistics( 148,  69,  23,  14, 2+10+0+0+0+16),
        2000: Statistics( 116,  97,  28,   9, 0+3+0+1+0+24),
        1999: Statistics( 108,  68,  24,   2, 0+3+0+5+0+26)
    }
}

for center in statistics:
    statistics[center]["avg"] = Statistics(
        sum(s[0] for s in statistics[center].values())/len(statistics[center]),
        sum(s[1] for s in statistics[center].values())/len(statistics[center]),
        sum(s[2] for s in statistics[center].values())/len(statistics[center]),
        sum(s[3] for s in statistics[center].values())/len(statistics[center]),
        sum(s[4] for s in statistics[center].values())/len(statistics[center])
    )
    last_four = [statistics[center][2014], statistics[center][2013],
                 statistics[center][2012], statistics[center][2011]]
    statistics[center]["avg4"] = Statistics(
        sum(s[0] for s in last_four)/4,
        sum(s[1] for s in last_four)/4,
        sum(s[2] for s in last_four)/4,
        sum(s[3] for s in last_four)/4,
        sum(s[4] for s in last_four)/4
    )

# Possible modes:
# transplant - simple printout of year tracked patient is transplanted
# statistics - printout of median years on list and on dialysis for all patients,
#              years waited and on dialysis for tracked patient
# histogram - histogram of patients in each sorting year in the waiting list
# results_histogram - histograms of years on the list and years on dialysis
mode = ['transplant']
start_year = {
    'CPMC': 2000,
    'UCLA': 1999,
    'UCD':  1999
}
insert_year = {
    'CPMC': 2012.0164,
    'UCLA': 2011.600,
    'UCD':  2015.8 # expected
}
dialysis_year = 2008.800
change_year = {
    'CPMC': 2015,
    'UCLA': 2006,
    'UCD':  2015
}
end_year = 2050

acceptable_ratio = {
    'CPMC': 0.801,#0.390,
    'UCLA': 0.794,#0.434,
    'UCD':  0.300 #0.467
}

perfect_match_ratio = {
    'CPMC': 0.064*.8, #... why the .8? just experimenting, or...?
    'UCLA': 0.114*.8,
    'UCD':  0.054*.8
}

if len(sys.argv) < 2 or sys.argv[1] not in statistics:
    print "Specify center."
    exit(1)
center = sys.argv[1]

current_year = start_year[center]

class Patient:

    def __init__(self, yl, yd, epts=None, tracked=False):
        self.listed = yl
        self.dialysis = yd
        if epts == None:
            self.epts = random.uniform(0,1)
        else:
            self.epts = epts
        self.tracked = tracked

    def key(self):
        if current_year < change_year[center]:
            return self.listed
        else:
            return self.dialysis

    def __cmp__(self, other):
        return cmp(other.key(), self.key())

    def __hash__(self):
        return hash((self.listed, self.dialysis, self.tracked))

    def __str__(self):
        if self.tracked:
            return "Patient(*%.3f, %.3f*)" % (self.listed, self.dialysis)
        else:
            return "Patient(%.3f, %.3f)" % (self.listed, self.dialysis)

    def __int__(self):
        return int(self.key())

    def __float__(self):
        return float(self.key())



waiting_list = []
years_waited = []
years_on_dialysis = []
tracked_patient = Patient(insert_year[center], dialysis_year, True)

def print_list(n=None):
    global waiting_list
    print '['
    print "  ",
    count = 0
    for patient in waiting_list:
        if patient.tracked:
            print "(*%.2f, %.2f*)" % (patient.listed, patient.dialysis),
        else:
            print "( %.2f, %.2f )" % (patient.listed, patient.dialysis), 
        count += 1
        if count % 3 == 0:
            print
            print "  ",
        if n != None and count >= n:
            break
    if count % 3 != 0:
        print
    print ']'

def print_stats():
    median_list = years_waited[len(years_waited)/2]
    median_dialysis = years_on_dialysis[len(years_on_dialysis)/2]
    patient_list = current_year+1 - insert_year[center]
    patient_dialysis = current_year+1 - dialysis_year
    print "Years waited: %.2f" % (patient_list)
    print "Years on dialysis: %.2f" % (patient_dialysis)
    print "Median years waited: %.2f" % (median_list)
    print "Median years on dialysis: %.2f" % (median_dialysis)

def calc_histogram():
    global waiting_list
    hist = {}
    for patient in waiting_list:
        hist[int(patient.key())] = hist.get(int(patient.key()), 0) + 1
    return hist

def show_histogram(l):
    bins = abs(int(max(l)) - int(min(l)))
    pylab.hist([float(p) for p in l], bins=bins)
    pylab.show()


def dialysis_dist():
    return random.gauss(1.5, 3.0)

def deceased_tx_dist():
    return int(min(len(waiting_list)-1,max(0,random.gauss(.95*len(waiting_list),.03*len(waiting_list)))))

def perfect_match_dist():
    if center == 'CPMC':
        return random.randint(int(len(waiting_list)*0.7), len(waiting_list)-1)
    else:
        return random.randint(0,len(waiting_list)-1)

def living_tx_dist():
    return int(min(len(waiting_list)-1,max(0,random.gauss(.1*len(waiting_list),.1*len(waiting_list)))))

def removed_dist():
    return random.randint(0,len(waiting_list)-1)

def fill_unique(n, pdf, allow_tracked=False):
    dist = set()
    while len(dist) < n:
        dist |= set(pdf() for i in xrange(n-len(dist)))
        if not allow_tracked and current_year >= int(insert_year[center]):
            try:
                dist.remove(waiting_list.index(tracked_patient))
            except KeyError:
                pass
            except ValueError:
                pass
    return sorted(dist, reverse=True)

def fill_deceased(n):
    if current_year < int(insert_year[center]):
        return fill_unique(n, deceased_tx_dist)
    else:
        dist = set()
        while len(dist) < int(n*perfect_match_ratio[center]):
            dist |= set(perfect_match_dist() for i in xrange(int(n*perfect_match_ratio[center]-len(dist))))
        while len(dist) < n:
            dist |= set(deceased_tx_dist() for i in xrange(n-len(dist)))
        return sorted(dist, reverse=True)



def check(action, year, candidates):
    if DEBUG:
        print action, len(candidates)
    if action == 'deceased_transplant':
        for p in candidates:
            years_waited.append(current_year-p.listed+1)
            years_on_dialysis.append(current_year-p.dialysis+1)
    if tracked_patient in candidates:
        if action == 'deceased_transplant':
            if (not (int(tracked_patient.key()) == current_year and
                random.uniform(0,1) < (tracked_patient.key() - current_year)) and
                random.uniform(0,1) < acceptable_ratio[center]):
                    done(action, year)
            else:
                bisect.insort(waiting_list, tracked_patient)
        else:
            # reset everything
            #print action, "in", year, "- resetting..."
            init()
            return True
    return False

def init():
    global waiting_list, current_year, start_year, statistics
    current_year = start_year[center]
    waiting_list = []
    years_waited = []
    years_on_dialysis = []
    if DEBUG:
        for i in xrange(200):
            year = start_year[center]+random.uniform(-6,1)
            waiting_list.append(Patient(year, year-dialysis_dist()))
    else:
        for i in xrange(statistics[center][start_year[center]][0]):
            year = start_year[center]+random.uniform(-6,1)
            waiting_list.append(Patient(year, year-dialysis_dist()))
    waiting_list.sort()

def done(result, year):
    if 'transplant' in mode:
        if result == "deceased_transplant":
            print "Transplanted in", year
        elif result == "nothing":
            print "Still on the list in", year

    if 'statistics' in mode:
        print_stats()
        #print sorted(years_waited)
    if 'histogram' in mode:
        show_histogram(waiting_list)
        #print_list()
    if 'results_histogram' in mode:
        show_histogram(years_waited)
        show_histogram(years_on_dialysis)

    exit(0)



init()

#print_list()
while current_year < end_year:
    current_year += 1

    try:
        new_listings = statistics[center][current_year].additions
        deceased_transplants = statistics[center][current_year].deceased_transplants
        living_transplants = statistics[center][current_year].living_transplants
        removed = statistics[center][current_year].removals
    except KeyError:
        new_listings = statistics[center]["avg4"].additions
        deceased_transplants = statistics[center]["avg4"].deceased_transplants
        living_transplants = statistics[center]["avg4"].living_transplants
        removed = statistics[center]["avg4"].removals

    try:
        target_length = statistics[center][current_year+1].initial_length
    except KeyError:
        target_length = None


    #print current_year, len(waiting_list), new_listings, deceased_transplants, living_transplants, removed

    if DEBUG:
        new_listings = 60
        deceased_transplants = 14
        living_transplants = 5
        removed = 33

    if current_year == int(insert_year[center]):
        bisect.insort(waiting_list, tracked_patient)
        #print "Inserting patient at %d/%d in %d" % (waiting_list.index(tracked_patient), len(waiting_list), current_year)
        #print_list()

    for i in xrange(new_listings):
        year = current_year+random.uniform(0,1)
        bisect.insort(waiting_list, Patient(year, year-dialysis_dist()))

    if current_year == change_year[center]:
        #print "Sorting by dialysis year"
        #if current_year >= insert_year[center]:
        #    old_i = waiting_list.index(tracked_patient)
        #    print "Was at %d/%d" % (old_i, len(waiting_list))
        waiting_list.sort() #sort by dialysis vintage
        #if current_year >= insert_year[center]:
        #    new_i = waiting_list.index(tracked_patient)
        #    print "Now at %d/%d" % (new_i, len(waiting_list))
        #    print "Jumped ahead by %d%%" % (int(float(new_i-old_i)/len(waiting_list)*100))
        #print_list()

    check('deceased_transplant', current_year, [waiting_list.pop(i) for i in fill_deceased(deceased_transplants)])
    check('living_transplant', current_year, [waiting_list.pop(i) for i in fill_unique(living_transplants, living_tx_dist)])
    check('removed', current_year, [waiting_list.pop(i) for i in fill_unique(removed, removed_dist)])
    while target_length != None and len(waiting_list) > target_length:
        check('removed', current_year, [waiting_list.pop(i) for i in fill_unique(len(waiting_list)-target_length, removed_dist)])
    if target_length != None and target_length < len(waiting_list):
        for i in xrange(target_length - len(waiting_list)):
            year = current_year+random.uniform(0,1)
            bisect.insort(waiting_list, Patient(year, year-dialysis_dist()))

    #if current_year >= insert_year and tracked_patient not in waiting_list:
        #print "Sneakily removed patient in", current_year
        #sleep(3)

    if DEBUG:
        print "============"
        print current_year
        print "----"
        print_list()

try:
    print "Still on list in", current_year, ",", waiting_list.index(tracked_patient), "of", len(waiting_list)
except ValueError:
    print len(waiting_list)
    print waiting_list[0], waiting_list[len(waiting_list)-1]
