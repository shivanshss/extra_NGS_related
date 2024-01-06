#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:23:58 2023

@author: bloodmark
"""

from pyGenomeTracks import Genome, HorizontalChromosome, LineTrack, BedTrack

# Function to create an ideogram using pyGenomeTracks
def create_ideogram(input_bed, output_pdf):
    # Create a Genome object
    genome = Genome(chromosomes=['your_chromosome_name'])

    # Add a horizontal chromosome track
    genome.add_track(HorizontalChromosome())

    # Add a BED track for your mutations
    genome.add_track(
        BedTrack(
            name='Mutations',
            bed_file=input_bed,
            height=2,
            color='red',
            labels=True
        )
    )

    # Plot the ideogram and save to a PDF file
    genome.plot(output_pdf, from_pos=0, to_pos=1000000)

if __name__ == "__main__":
    input_bed = "output.bed"  # Replace with the path to your BED file
    output_pdf = "ideogram.pdf"

    create_ideogram(input_bed, output_pdf)
