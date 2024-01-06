#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:48:38 2023

@author: bloodmark
"""

delta_file_path = '/home/bloodmark/workarea/nucmer/tcas6_final_all.delta'

with open(delta_file_path, 'r') as delta_file:
    for line in delta_file:
        if line.startswith('>'):
            fields = line.split()
            if fields[0] == '>1':
                if fields[1] == 'I':
                    print(f"Insertion at position {fields[2]} in sequence 2, length {fields[3]}")
                elif fields[1] == 'D':
                    print(f"Deletion at position {fields[2]} in sequence 1, length {fields[3]}")
                elif fields[1] == 'G':
                    print(f"Gap at position {fields[2]} in both sequences, length {fields[3]}")
