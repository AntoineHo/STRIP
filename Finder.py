#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import re
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

class Region :
    def __init__(self, start, stop, id) :
        self.start = start
        self.stop = stop
        self.id = id

        self.length = self.stop - self.start

    def __str__(self) :
        return "{}:{}-{} ({}bp)".format(self.id, self.start, self.stop, self.length)

    def __eq__(self, region) :
        return self.length == region.length and self.id == region.id

    def Overlap(self, region) :
        # Check if same contig and overlap
        if self.id != region.id :
            return False
        #       start--------------------------------stop
        #   start ------------stop
        #   start ----------------------------stop
        #                                       start----stop
        # ==============A                  B===============
        #
        A = region.start
        B = region.stop
        if (self.start >= A and self.start <= B) or (self.stop >= A and self.stop <= B) or (self.start <= A and self.stop >= B) :
            return True
        else :
            return False

class RestrictionProfile :
    def __init__(self, fasta, site, bed_to_ignore=None) :

        self.fasta = fasta
        if os.path.exists(self.fasta) :
            self.fasta = os.path.abspath(self.fasta)
        else :
            raise Exception("ERROR: file does not exist!")

        self.site = site
        """ UNWORKING FOR NOW """
        self.bed_to_ignore = bed_to_ignore # List of region objects

        self.regions_to_ignore = []
        if self.bed_to_ignore != None and os.path.exists(self.bed_to_ignore) :
            self.bed_to_ignore = os.path.abspath(self.bed_to_ignore)
            f = open(self.bed_to_ignore, "r")
            for line in f :
                sl = line.strip().split("\t")
                current_region = Region(int(sl[1]), int(sl[2]), sl[0])
                self.regions_to_ignore.append(current_region)

        self.sequences = []
        for record in SeqIO.parse(self.fasta, "fasta") :
            self.sequences.append([record.id, str(record.seq)])

        self.profile = None

    def __str__(self) :
        if self.profile == None :
            return "Restriction is not done yet!"
        else :
            return "{} fragments found in perfect restriction situation!".format(len(self.profile))

    def SimulateDeletion(self, minsize=100000, maxsize=500000, mincontigsize=10000000, ndel=1) :

        # 0. Get all sequences longer than size :
        to_be_deleted = []
        for id, seq in self.sequences :
            if len(seq) > maxsize :
                to_be_deleted.append([id, seq])

        will_be_deleted = []
        already_chosen_contig_ids = []
        # 1. Choose randomly n_del element in the list
        for i in range(ndel) :
            contig = random.choice(to_be_deleted)
            while len(contig[1]) < mincontigsize or contig[0] in already_chosen_contig_ids:
                contig = random.choice(to_be_deleted)
            already_chosen_contig_ids.append(contig[0])
            will_be_deleted.append(contig)

        # 2. delete part of the sequence
        for id, seq in will_be_deleted :
            print("{} = {}".format(id, len(seq)))

        deleted = []
        for id, seq in will_be_deleted :
            undeleted_length = len(seq) # Gets undeleted length of the fragment
            size = random.randint(minsize, maxsize) # Choose a size to delete
            start = random.randint(0, undeleted_length-size-1) # Choose a start point on the contig
            end = start + size # Computes the endpoint of the deletion

            deleted.append([id, seq[0:start]+seq[end:]])  # Modifies the sequence

            print("Deletion on {} from {} to {}\n\tDeletion size: {}\n\tSize before: {}\n\tSize after: {}".format(id, start, end, size, undeleted_length, undeleted_length-size))

        for id, seq in deleted :
            print("{} = {}".format(id, len(seq)))

        # 3. Completes the deleted seq list
        deleted_ids = [id for id, seq in deleted]
        for id, seq in self.sequences :
            if id not in deleted_ids :
                deleted.append([id, seq])

        return self.Restrict(deleted)

    def Restrict(self, sequences, save_profile=False) :
        print("Restricting...")
        all_regions = []
        # For each contig in the list
        for id, seq in sequences :
            #print(record.id)
            # Find all coordinates of the restriction site
            all_coords = [m.start() for m in re.finditer(self.site, seq)]
            all_coords += [m.start() for m in re.finditer(self.complement(self.site), seq)]
            all_coords = sorted(list(all_coords))
            #print(len(all_coords))
            # Loops over all coordinates and get regions of fragmentation
            last = 0
            current = 0
            # For each coord : get fragmentation then check if fragment overlaps with regions to exclude if so discard fragment (telomerics) ==> UNDERESTIMATION
            for coord in all_coords :
                current = coord
                current_region = Region(last, current, id)
                #print(current_region)
                if len(self.regions_to_ignore) != 0 :
                    for region in self.regions_to_ignore :
                        if not current_region.Overlap(region) : # If region does not overlap with the excluded regions add it to the region to consider list
                            print("not_overlapping!")
                            all_regions.append(current_region)
                else :
                    all_regions.append(current_region)
                last = coord

        print("Done!")
        if save_profile :
            self.profile = all_regions

        return all_regions

    def complement(self, seq) :
        base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        letters = list(seq)
        letters = [base_complement[base] for base in letters]
        return ''.join(letters)

    def revcom(self, seq) :
        return self.complement(seq[::-1])


    def PlotProfile(self, profile, miny=200000, maxy=2000000) :

        plt.style.use("default")
        fig, ax = plt.subplots(figsize=(3,20))

        x_pos = [0.5 for r in profile]
        y_pos = [r.length for r in profile]

        pcl = []
        for x, y in zip(x_pos, y_pos) :
            pcl.append(mpatches.Rectangle(xy=[x,y], width=1, height=5000))

        pc = PatchCollection(pcl, facecolor="k", alpha=0.6)
        ax.add_collection(pc)

        ax.set_yticks(y_pos)

        ax.set_xlim(0,2)
        ax.set_ylim(miny, maxy)
        ax.get_xaxis().set_visible(False)

        plt.show()

        return fig, ax

    def PlotMultiProfile(self, profile_list, miny=200000, maxy=2000000, labels=None) :

        unique_regions, missing_regions = self.CompareProfiles(profile_list)

        plt.style.use("default")
        np = len(profile_list)
        fig, ax = plt.subplots(nrows=1,ncols=np,figsize=(3*np,20))

        for n, profile in enumerate(profile_list) :
            #print(len(unique_regions))
            #print(len(profile))

            x_pos = [0.5 for r in profile]
            y_pos = [r.length for r in profile]

            uniq = []
            for r in profile :
                if r in unique_regions :
                    uniq.append(True)
                else :
                    uniq.append(False)

            pclUniq = []
            pcl = []
            for x, y, isUniq in zip(x_pos, y_pos, uniq) :
                #print(isUniq)
                if isUniq :
                    pclUniq.append(mpatches.Rectangle(xy=[x,y], width=1, height=5000))
                else :
                    pcl.append(mpatches.Rectangle(xy=[x,y], width=1, height=5000))


            pclMissing = []
            for r in missing_regions :
                if r not in profile :
                    pclMissing.append(mpatches.Rectangle(xy=[0.5,r.length], width=1, height=5000))

            pcMissing = PatchCollection(pclUniq, facecolor="g", alpha=0.4)
            pcUniq = PatchCollection(pclUniq, facecolor="r", alpha=0.9)
            pc = PatchCollection(pcl, facecolor="k", alpha=0.6)
            ax[n].add_collection(pcMissing)
            ax[n].add_collection(pcUniq)
            ax[n].add_collection(pc)

            ax[n].set_yticks(y_pos)

            ax[n].set_xlim(0,2)
            ax[n].set_ylim(miny, maxy)
            ax[n].get_xaxis().set_visible(False)

        if labels != None :
            if len(labels) == len(profile_list) :
                for n, l in enumerate(labels) :
                    ax[n].set_title(l, fontsize=16)
            else :
                print("WARNING: Unmatching label list length: no labels added!")

        plt.subplots_adjust(wspace=0.5, hspace=None)
        plt.show()

        return fig, ax

    def CompareProfiles(self, profile_list) :

        all_fragments = []
        for profile in profile_list :
            for region in profile :
                all_fragments.append(region)

        missing_regions = []
        unique_regions = []
        regions_list = []
        for region in all_fragments :
            regions_list.append([region, all_fragments.count(region)])

        for r, c in regions_list :
            if c == 1 :
                unique_regions.append(r)
            if c < len(profile_list) :
                missing_regions.append(r)

        """ For debugging
        for region in unique_regions :
            print(region)
        """

        return unique_regions, missing_regions

"""
def main():
    ### Main program
    # Code goes over here.

    return 0

if __name__ == "__main__":
    main()
"""
