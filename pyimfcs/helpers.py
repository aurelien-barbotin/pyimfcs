#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 14:17:00 2022

@author: aurelienb
"""

def time_str_2_min(creation_date):
    creation_time = creation_date[-8:].split(':')
    crt_min = int(creation_time[0])*60+int(creation_time[1])+float(creation_time[2])/60
    return crt_min