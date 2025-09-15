from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import sys
import csv
from tkinter import Tk, filedialog
import os
import matplotlib.pyplot as plt
import numpy as np

def getbasepeaks(record):
    abifraw = record.annotations['abif_raw']
    traceA = abifraw['DATA10']
    traceC = abifraw['DATA12']
    traceG = abifraw['DATA9']
    traceT = abifraw['DATA11']
    positions = abifraw['PLOC2']
    seq = str(record.seq)

    results = []
    for i, base in enumerate(seq):
        pos = positions[i]
        results.append((base, pos, traceA[pos], traceC[pos], traceG[pos], traceT[pos]))
    return results

def printpeaks(peaks):
    print("Base\tPos\tA\tC\tG\tT")
    baseseq = ""
    for base, pos, a, c, g, t in peaks:
        baseseq += base
        print(f"{base}\t{pos}\t{a}\t{c}\t{g}\t{t}")

def printalignmentchunks(alignment, width=60):
    seqA, seqB, score, start, end = alignment
    matchline = []
    
    for a, b in zip(seqA, seqB):
        if a == b and a != '-':
            matchline.append('|')
        else:
            matchline.append(' ')
    matchline = ''.join(matchline)
    
    for i in range(0, len(seqA), width):
        print(seqA[i:i+width])
        print(matchline[i:i+width])
        print(seqB[i:i+width])
        print()  
    print(f"Score={score}")

def mapwindow(ab1seq, refseq, windowstart, windowend):
    alignments = pairwise2.align.localms(refseq, ab1seq, 2, -1, -1, -1)
    if not alignments:
        raise ValueError("No alignment found between reference and AB1 sequence")
    
    bestalignment = alignments[0]
    alignedref, alignedab1, score, start, end = bestalignment

    printalignmentchunks(bestalignment, width=120)

    refpositions = []
    ab1positions = []
    currentref = 0
    currentab1 = 0
    
    for refchar, ab1char in zip(alignedref, alignedab1):
        if refchar == '-':
            refpositions.append(None)
        else:
            refpositions.append(currentref)
            currentref += 1
        
        if ab1char == '-':
            ab1positions.append(None)
        else:
            ab1positions.append(currentab1)
            currentab1 += 1
    
    segmentindices = []
    for idx, pos in enumerate(refpositions):
        if pos is not None and windowstart <= pos <= windowend:
            segmentindices.append(idx)
    
    if not segmentindices:
        raise ValueError("Reference window not found in alignment")
    
    ab1segmentpositions = []
    for idx in segmentindices:
        if ab1positions[idx] is not None:
            ab1segmentpositions.append(ab1positions[idx])
    
    if not ab1segmentpositions:
        raise ValueError("Window not found in aligned AB1 sequence")
    
    ab1start = min(ab1segmentpositions)
    ab1end = max(ab1segmentpositions)
    ab1windowseq = ab1seq[ab1start:ab1end+1]
    
    return ab1start, ab1end, ab1windowseq

def getextended(ab1seq, start, end, seq):
    newstart = end + 1
    newend = (end + 30)  
    
    newab1sequence = ab1seq[newstart:newend + 1]
    
    return newstart, newend, newab1sequence

def lookat12bpdeletion(ab1seq, startext, endext, seqext, peaks, similaritysequenceCR, similaritysequenceCS):
    firstdeletionsite = ab1seq[startext:startext + 12]
    print(f"First 12 bp of extended sequence: {firstdeletionsite}")

    print("\nBase-by-base CR vs CS comparison:")
    print("Index | CR Base | CS Base | %CR | %CS | Difference")
    print("-" * 60)

    results = []
    for i in range(12):
        crbase = similaritysequenceCR[i] if i < len(similaritysequenceCR) else "N/A"
        csbase = similaritysequenceCS[i] if i < len(similaritysequenceCS) else "N/A"

        if startext + i < len(peaks):
            actualbase, pos, A, C, G, T = peaks[startext + i]
            peakdict = {'A': A, 'C': C, 'G': G, 'T': T}
            totalpeaks = sum(peakdict.values())

            percentageCR = (peakdict.get(crbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0
            percentageCS = (peakdict.get(csbase, 0) / totalpeaks * 100) if totalpeaks > 0 else 0.0

            if percentageCR + percentageCS > 0:
                difference = ((percentageCR - percentageCS) / (percentageCR + percentageCS)) * 100
                if difference < 0:
                    difference = 0.0
            else:
                difference = 0.0

        else:
            percentageCR = 0.0
            percentageCS = 0.0
            difference = 0.0

        print(f"{i+1:5} | {crbase:7} | {csbase:7} | {percentageCR:6.1f} | {percentageCS:6.1f} | {difference:8.1f}")
        results.append({
            'CRbase': crbase,
            'CSbase': csbase,
            'percentageCR': percentageCR,
            'percentageCS': percentageCS,
            'difference': difference
        })

    return results

def main():
    root = Tk()
    root.withdraw()
    
    print("Select AB1 files: ")
    filenames = filedialog.askopenfilenames(
        title="Select AB1 Files",
        filetypes=[("AB1 files", "*.ab1"), ("All files", "*.*")]
    )
    if not filenames:
        print("No files selected")
        sys.exit(1)
    
    print("Select loaction: ")
    csvpath = filedialog.asksaveasfilename(
        title="Save CSV As",
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    
    if not csvpath:
        print("not saved.")
        sys.exit(1)

    allresults = []
    samplenames = []
    similaritysequenceCR = "GCTGAATGCCCT"
    similaritysequenceCS = "TGCCTTTCGATC"
    
    for filename in filenames:
        print(f"\nProcessing {os.path.basename(filename)}")
        samplenames.append(os.path.basename(filename))
        
        try:
            record = SeqIO.read(filename, "abi")
            peaks = getbasepeaks(record)
            
            ab1seq = str(record.seq)
            refsequence = "TTGCGGCGTGGCCTATCCGGGCGAACTTTTGGCCGTGATGGGCAGTTCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCGCAGGGCATCCAAGTATCGCCATCCGGGATGCGACTGCTCAATGGCCAACCTGTGGACGCCAAGGAGATGCAGGCCAGGTGCGCCTATGTCCAGCAGGATGACCTCTTTATCGGCTCCCTAACGGCCAGGGAACACCTGATTTTCCAAGCCATGGTGCGGATGCCACGACATCTGACCTATCGGCAGCGAGTGGCCCGCGTGGATCAGGTGATCCAGGAGCTTTCGCTCAGCAAATGTCAGCACACGATCATCGGTGTGCCCGGCAGGGTGAAAGGTCTGTCCGGCGGAGAAAGG"
            findpattern = "GTTCCGGTGCCGGAAAGACGACCCT"
            
            windowstart = refsequence.find(findpattern)
            if windowstart == -1:
                print(f"Pattern '{findpattern}' not found.")
                continue
                
            windowend = windowstart + len(findpattern) - 1
            
            start, end, seq = mapwindow(ab1seq, refsequence, windowstart, windowend)
            startext, endext, seqext = getextended(ab1seq, start, end, seq)
            percentages = lookat12bpdeletion(ab1seq, startext, endext, seqext, peaks, similaritysequenceCR, similaritysequenceCS)            
            allresults.append(percentages)
            
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            allresults.append([{'CRbase': similaritysequenceCR[i] if i < len(similaritysequenceCR) else "N/A", 'CSbase': similaritysequenceCS[i] if i < len(similaritysequenceCS) else "N/A", 'percentageCR': 0.0, 'percentageCS': 0.0, 'difference': 0.0} for i in range(12)])
            continue
        
    print("\nAll samples to write:", samplenames)
    print("All results size:", len(allresults), "x", len(allresults[0]) if allresults else 0)

    with open(csvpath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        for samplename, sampleresult in zip(samplenames, allresults):
            writer.writerow([samplename])

            header = ['Index', 'CR Base', 'CS Base', 'Percentage CR', 'Percentage CS', 'Difference']
            writer.writerow(header)

            for i in range(12):
                row = [i + 1, sampleresult[i]['CRbase'], sampleresult[i]['CSbase'],
                    sampleresult[i]['percentageCR'], sampleresult[i]['percentageCS'], sampleresult[i]['difference']]
                writer.writerow(row)

            writer.writerow([])
    
    print(f"\nResults saved @ {csvpath}")
    
    input("done")

if __name__ == "__main__":
    main()