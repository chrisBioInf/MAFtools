#!/usr/bin/env python3
"""
MAF Alignment Extractor
Extracts MAF alignment blocks based on BED coordinates
"""

import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
from typing import List, Tuple, Dict, Optional

import sys
import os
from optparse import OptionParser


__author__ = "Christopher Klapproth"
__institution__= "University Leipzig"
__credits__ = []
__license__ = "GPLv2"
__version__="1.0.0"
__maintainer__ = "Christopher Klapproth"
__email__ = "christopher@bioinf.uni-leipzig.de"
__status__ = "Development"


THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(THIS_FILE)


usage_statement = "Usage: MAFtools [program] [options] [MAF file], where `program` is one of: `describe`, `blockdescribe`, `merge`, `filter`, `window`, `select`, `toBed`"


strand_dict = {
    1 : '+',
    -1 : '-'
}

#
#
#   Define Helper Functions
#
#

def sortRecords(records: list) -> None:
    n = len(records)
    if n < 3:
        return
    swapped = False
    
    for i in range(1, n-1):
        for j in range(1, n-i-1):
            if records[j].id > records[j + 1].id:
                swapped = True
                records[j], records[j + 1] = records[j + 1], records[j]
         
        if not swapped:
            return
        
        
def eliminate_consensus_gaps(records: list) -> list:
    ungapped_seqs = []
    seq_matrix = np.array([list(record.seq) for record in records])
    for i in range(0, len(records)):
        seq_ = ""
        for c in range(0, len(seq_matrix[i])):
            if seq_matrix[i][c] != "-":
                seq_ = seq_ + seq_matrix[i][c]
                continue
            for j in range(0, len(seq_matrix)):
                if seq_matrix[j][c] != "-":
                    seq_ = seq_ + seq_matrix[i][c]
                    break
        ungapped_seqs.append(seq_)
    
    for i in range(0, len(records)):
        records[i].seq = Seq(ungapped_seqs[i])
    
    return records


def max_gap_seqs(records: list, max_gaps: int =0, reference: bool =True) -> list:
    start_index = 0
    if reference == True:
        start_index = 1
    length = len(records)
    
    if start_index >= length:
        return []
    
    columns = len(str(records[0].seq))
    drop_indices = set()

    if columns == 0:
        return []
    
    for i in range(start_index, length):
        gaps = records[i].seq.count("-")
        gap_fraction = gaps / columns
        if (gap_fraction > max_gaps) or (gaps == columns):
            drop_indices.add(i)

    return [records[n] for n in range(0, length) if n not in drop_indices]

def coordinate_distance(end1, start2):
    return start2 - end1


def pairwise_sequence_identity(seq1, seq2):
    count = 0
    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]:
            count += 1
    return count / len(seq1)
 
 
def local_species_consensus(block1, block2):
    block1_ids = set([str(seq.id) for seq in block1])
    block2_ids = set([str(seq.id) for seq in block2])
    consensus_species = len(block1_ids.intersection(block2_ids)) + len(block2_ids.intersection(block1_ids)) 
    consensus_score  = consensus_species / (len(block1_ids) + len(block2_ids))
    return consensus_species, consensus_score


def concat_with_bridge(seq1, seq2, offset, max_offset):
    n_gaps = max_offset - offset
    seq_ = Seq(str(seq1) + "N"*offset + '-'*n_gaps + str(seq2)) 
    return seq_


def filter_by_seq_length(records, reference=True):
    length_dict = {str(record.id) : len(record.seq.replace('-', '')) for record in records}
    record_dict = {str(record.id) : record for record in records}
    mu = np.mean([x for x in length_dict.values()])
    std = np.std([x for x in length_dict.values()])
    
    count = 0
    
    for record in records:
        count += 1
        if (count == 1) and reference:
            continue
        
        x = str(record.id)
        if std > 0:
            z = abs(length_dict.get(x) - mu) / std
        else:
            z = 0
        if z > 2:
            record_dict.pop(x)
    
    return [record for record in record_dict.values()]


def merge_blocks(block1, block2, reference=True, offset_threshold=0):
    block1_ids = [str(record.id) for record in block1]
    block2_ids = [str(record.id) for record in block2]
    
    if len(block1_ids) >= len(block2_ids):
        consensus_blocks = set(block1_ids).intersection(block2_ids)
    else:
        consensus_blocks = set(block2_ids).intersection(block1_ids)
    
    block1_records = [record for record in block1 if record.id in consensus_blocks]
    block2_records = [record for record in block2 if record.id in consensus_blocks]
    sortRecords(block1_records)
    sortRecords(block2_records)

    merged_records = []
    offsets = [coordinate_distance(block1_records[i].annotations["start"] + block1_records[i].annotations["size"], block2_records[i].annotations["start"]) for i in range(0, len(block1_records))]
    block1_merge = []
    block2_merge = []
     
    for i in range(0, len(block1_records)):
        offset = offsets[i]
        record_, record2 = block1_records[i], block2_records[i]
        strand1, strand2 = record_.annotations["strand"], record2.annotations["strand"]
        start1, end1 = record_.annotations["start"], record_.annotations["start"] + record_.annotations["size"]
        start2, end2 = record2.annotations["start"], record2.annotations["start"] + record2.annotations["size"]
        if (offset < 0):
            continue
        if (offset > offset_threshold) or (strand1 != strand2):
            continue
        block1_merge.append(record_)
        block2_merge.append(record2)

    offsets = [coordinate_distance(block1_merge[i].annotations["start"] + block1_merge[i].annotations["size"], block2_merge[i].annotations["start"]) for i in range(0, len(block1_merge))]


    if len(offsets) < 2:
        # print_maf_alignment(block1)
        return block2

    max_offset = max(offsets)
    index = 0

    for (record_, record2) in zip(block1_merge, block2_merge):
        offset = offsets[index]
        record_.seq = concat_with_bridge(record_.seq, record2.seq, offset, max_offset)
        record_.annotations["size"] = (end1 - start1) + offset + (end2 - start2)
        merged_records.append(record_)
        index += 1
    
    return merged_records


def check_block_viability(block1, block2, species_consensus_threshold, block_distance_threshold, block_length_threshold, min_seqs=2):
    merge_flag = True
    reference1 = block1[0]
    reference2 = block2[0]
    end1 = reference1.annotations["start"] + reference1.annotations["size"]
    start2 = reference2.annotations["start"]
    consensus_species, consensus_score = local_species_consensus(block1, block2)

    def check_distance_viability(block1, block2):
        end1_dict = {r.id : (r.annotations.get('start') + r.annotations.get('size')) for r in block1}
        start2_dict = {r.id : r.annotations.get('start') for r in block2}
        viable_records = []

        for name in end1_dict.keys():
            distance = start2_dict.get(name, 0) - end1_dict.get(name)
            if (distance >= 0) and (distance <= block_distance_threshold):
                viable_records.append(name)

        return len(viable_records) / len(block1)

    
    if coordinate_distance(end1, start2) > block_distance_threshold:
        merge_flag = False
    elif reference1.annotations["strand"] != reference2.annotations["strand"]:
        merge_flag = False
    elif (len(reference1.seq) + len(reference2.seq)) > block_length_threshold:
        merge_flag = False 
    elif consensus_score < species_consensus_threshold:
        merge_flag = False
    elif consensus_species < min_seqs:
        merge_flag = False
    elif reference1.id != reference2.id:
        merge_flag = False
    # elif check_distance_viability(block1, block2) < species_consensus_threshold:
    #    merge_flag = False
         
    return merge_flag



class MAFObject:
    def __init__(self, maf_file: str, bed_file: str, reference_species: str):
        """
        Initialize the MAF extractor.
        
        Args:
            maf_file: Path to MAF alignment file
            bed_file: Path to BED annotation file
            reference_species: Name of the reference species in MAF
        """
        self.maf_file = maf_file
        self.bed_file = bed_file
        self.reference_species = reference_species
        
    def parse_bed(self) -> List[Tuple[str, int, int]]:
        """Parse BED file and return list of (chrom, start, end) tuples."""
        intervals = []
        with open(self.bed_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        intervals.append((chrom, start, end))
        
        return intervals
    
    def extract_alignment_region(self, alignment, ref_start: int, ref_size: int, 
                                bed_start: int, bed_end: int) -> Optional[MultipleSeqAlignment]:
        """
        Extract a region from an alignment based on BED coordinates.
        
        Args:
            alignment: Bio.Align.MultipleSeqAlignment object
            ref_start: Start position of reference in MAF block
            ref_size: Size of aligned reference sequence (excluding gaps)
            bed_start: Start position from BED file
            bed_end: End position from BED file
            
        Returns:
            New alignment for the extracted region or None if no overlap
        """
        # Calculate overlap
        ref_end = ref_start + ref_size
        overlap_start = max(ref_start, bed_start)
        overlap_end = min(ref_end, bed_end)

        print("#################################")
        print(bed_start, bed_end)
        print(ref_start, ref_end)

        if overlap_start >= overlap_end:
            print("-> No overlap.")
            return None  # No overlap
        
        print("Overlap: ", overlap_start, " ", overlap_end)
        # Find reference sequence
        ref_record = None
        ref_idx = None
        for i, record in enumerate(alignment):
            if record.id.startswith(self.reference_species):
                ref_record = record
                ref_idx = i
                break
                
        if ref_record is None:
            return None
            
        # Calculate positions in the gapped alignment
        ref_seq = str(ref_record.seq)
        ungapped_pos = ref_start
        start_col = None
        end_col = None
        
        for col, base in enumerate(ref_seq):
            if base != '-':
                if ungapped_pos == overlap_start:
                    start_col = col
                if ungapped_pos == overlap_end:
                    end_col = col
                    break
                ungapped_pos += 1 
        
        # if start_col is None or end_col is None:
        #    return None
        if end_col is None:
            end_col = len(ref_seq)-1
        if start_col is None:
            start_col = 0

        print("Columns:", start_col, " ", end_col)   
            
        # Extract alignment columns
        new_records = []
        for record in alignment:
            new_seq = str(record.seq)[start_col:end_col]
            new_record = SeqRecord(
                Seq(new_seq),
                id=record.id,
                description=record.description,
                annotations=record.annotations.copy() if hasattr(record, 'annotations') else {}
            )
            new_records.append(new_record)
            
        return MultipleSeqAlignment(new_records)
    
    def update_maf_coordinates(self, alignment, ref_start: int, ref_size: int,
                             bed_start: int, bed_end: int) -> Dict[str, Dict]:
        """
        Update MAF coordinates for extracted alignment.
        
        Returns dictionary with updated coordinates for each sequence.
        """
        coords = {}
        
        # Find reference sequence
        ref_record = alignment[0]
        for record in alignment:
            if record.id.startswith(self.reference_species):
                ref_record = record
                break
                
        if ref_record is None:
            return coords
            
        # Calculate new start position for reference
        ref_end = ref_start + ref_size
        overlap_start = max(ref_start, bed_start)
        overlap_end = min(ref_end, bed_end)
        new_ref_start = overlap_start
        new_ref_size = overlap_end - overlap_start
        
        # For each sequence, calculate new coordinates
        for record in alignment:
            seq_id = record.id
            parts = seq_id.split('.')
            species = parts[0] if parts else seq_id
            
            # Get original MAF info from annotations if available
            if hasattr(record, 'annotations'):
                orig_start = record.annotations.get('start', 0)
                orig_strand = record.annotations.get('strand', '+')
                orig_srcsize = record.annotations.get('srcSize', 0)
            else:
                orig_start = 0
                orig_strand = '+'
                orig_srcsize = 0
                
            # Count non-gap bases in extracted sequence
            new_size = sum(1 for base in str(record.seq) if base != '-')
            
            coords[seq_id] = {
                'start': new_ref_start if species == self.reference_species else orig_start,
                'size': new_size,
                'strand': strand_dict.get(orig_strand, orig_strand),
                'srcSize': orig_srcsize
            }
            
        return coords
    
    def format_maf_block(self, alignment, coords: Dict[str, Dict], score: float = 0.0) -> str:
        """Format an alignment block in MAF format."""
        lines = []
        
        # Add 'a' line
        lines.append(f"a score={score}")
        
        # Add 's' lines
        for record in alignment:
            seq_id = record.id
            if seq_id in coords:
                c = coords[seq_id]
                line = f"s {seq_id} {c['start']} {c['size']} {c['strand']} {c['srcSize']} {str(record.seq)}"
                lines.append(line)
                
        return '\n'.join(lines)
    
    def extract_blocks(self, output_file: str):
        """Main method to extract MAF blocks based on BED intervals."""
        # Parse BED intervals
        bed_intervals = self.parse_bed()
        print(bed_intervals)
        
        # Open output file
        with open(output_file, 'w') as out:
            # Write MAF header
            out.write("##maf version=1\n\n")
            
            # Process each MAF alignment block
            alignments = AlignIO.parse(self.maf_file, "maf")
            
            for alignment in alignments:
                # Get reference info from first sequence (assumed to be reference)
                ref_record = None
                for record in alignment:
                    if record.id.startswith(self.reference_species):
                        ref_record = record
                        break
                        
                if ref_record is None:
                    continue
                    
                # Get MAF coordinates
                ref_start = ref_record.annotations.get('start', 0)
                ref_size = ref_record.annotations.get('size', 0)
                
                # Check each BED interval
                for chrom, bed_start, bed_end in bed_intervals:
                    # Determine reference sequence
                    self.reference_species = chrom

                    # Extract if there's overlap
                    extracted = self.extract_alignment_region(
                        alignment, ref_start, ref_size, bed_start, bed_end
                    )
                    
                    if extracted:
                        # Update coordinates
                        coords = self.update_maf_coordinates(
                            extracted, ref_start, ref_size, bed_start, bed_end
                        )
                        
                        # Format and write
                        maf_block = self.format_maf_block(extracted, coords)
                        out.write(maf_block + "\n\n")


    def merge(self, options):
        args = args[1:]

        output_file = options.out_file
        
        alignments = AlignIO.parse(self.maf_file, "maf")
        block1 = next(alignments)
        merged_blocks_n = []
        local_merges = 1
        
        while block1:
            block2 = next(alignments, None)
            
            if not block2:
                records = eliminate_consensus_gaps(block1)
                records = max_gap_seqs(records, options.max_gaps, options.reference)
                merged_blocks_n.append(local_merges)
                
                # Append to output then end loop, as there are no further blocks:
                AlignIO.write(MultipleSeqAlignment(records), handle=open(options.out_file, 'a'), format="maf")
                break

            merge_flag = check_block_viability(block1, block2, 
                                            options.species_consensus,
                                            options.distance,
                                            options.length)
            if merge_flag == True: 
                block1 = merge_blocks(block1, block2, options.reference, 
                                    offset_threshold=options.distance)
                local_merges += 1
            else:
                records = eliminate_consensus_gaps(block1)
                records = max_gap_seqs(records, options.max_gaps, options.reference)
                merged_blocks_n.append(local_merges)
                local_merges = 1
                block1 = block2

                # Append to output:
                AlignIO.write(MultipleSeqAlignment(records), handle=open(options.out_file, 'a'), format="maf")

                        
    def toBed(self, output_bed_file: str):
        """
        Convert MAF alignment blocks to BED format.
        
        Creates a BED file where each line represents the coordinates of one 
        alignment block with respect to the reference sequence.
        
        Args:
            output_bed_file: Path to output BED file
        """
        bed_entries = []
        
        # Process each MAF alignment block
        alignments = AlignIO.parse(self.maf_file, "maf")
        
        for alignment in alignments:
            # Find reference sequence in this alignment block
            ref_record = None
            for record in alignment:
                if record.id.startswith(self.reference_species):
                    ref_record = record
                    break
            
            if ref_record is None:
                continue  # Skip blocks without reference sequence
            
            # Get MAF coordinates for reference
            ref_seq_name = ref_record.id
            ref_start = ref_record.annotations.get('start', 0)
            ref_size = ref_record.annotations.get('size', 0)
            
            # Calculate end position (BED uses 0-based, half-open intervals like MAF)
            ref_end = ref_start + ref_size
            
            # Add to BED entries
            bed_entries.append((ref_seq_name, ref_start, ref_end))
        
        # Write BED file
        with open(output_bed_file, 'w') as out:
            for seq_name, start, end in bed_entries:
                out.write(f"{seq_name}\t{start}\t{end}\n")


    def describe(self):
        alignments = AlignIO.parse(self.maf_file, "maf")

        data = {
            "Identifiers": [],
            "Sequences": [], 
            "Number of blocks": [],
            "Columns": [],  
            "Aligned columns": [], 
            "Aligned refseq nucleotides": [],
            "Average block length": [], 
        }
        sequence_set = set()
        identifier_set = set()
        n_blocks = 0
        n_columns = 0
        nucleotides = 0
        n_sequences_per_block = []

        for alignment in alignments:
            n_blocks += 1
            length = len(alignment[0].seq)
            n_columns += length
            nucleotides += length - str(alignment[0].seq).count('-')
            n_sequences = 0

            for record in alignment:
                identifier_set.add(str(record.id).split('.')[0])
                sequence_set.add(record.id)
                n_sequences += 1

            n_sequences_per_block.append(n_sequences)

        data = {
                "Identifiers": len(identifier_set),
                "Sequences": len(sequence_set), 
                "Number of blocks": n_blocks,
                "Columns": n_columns,  
                "Average sequences per block": round(sum(n_sequences_per_block) / n_blocks, 2), 
                "Aligned reference nucleotides": nucleotides,
                "Average block length": round(n_columns / n_blocks, 2), 
            }
        for (k, v) in data.items():
            print("%s: %s" % (k, v))

    def _calculate_pairwise_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate pairwise similarity between two aligned sequences.
        
        Similarity = matching positions / (total positions - gap-only positions)
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of equal length")
        
        matches = 0
        compared_positions = 0
        
        for c1, c2 in zip(seq1, seq2):
            # Skip positions where both sequences have gaps
            if c1 == '-' and c2 == '-':
                continue
                
            compared_positions += 1
            
            # Count matches (ignoring case)
            if c1 == c2 and c1 != '-':
                matches += 1
        
        if compared_positions == 0:
            return 0.0
        
        return matches / compared_positions
    
    
    def blockdescribe(self, output_tsv_file: str):
        """
        Generate a tab-delimited summary file with statistics for each alignment block.
        
        Creates a TSV file with one row per alignment block containing:
        - refseq: Reference sequence name
        - start: Start coordinate
        - size: Size of reference sequence (non-gap bases)
        - length: Number of columns in alignment
        - mean_pairwise_identity: Average pairwise identity to reference
        - gaps_fraction: Average gap fraction across all sequences
        
        Args:
            output_tsv_file: Path to output TSV file
        """
        # Open output file and write header
        with open(output_tsv_file, 'w') as out:
            # Write header
            out.write("refseq\tstart\tsize\tlength\tmean_pairwise_identity\tgaps_fraction\n")
            
            # Process each alignment block
            alignments = AlignIO.parse(self.maf_file, "maf")
            block_num = 0
            
            for alignment in alignments:
                block_num += 1
                
                # Find reference sequence (first one)
                if len(alignment) == 0:
                    continue
                    
                ref_record = alignment[0]
                ref_seq = str(ref_record.seq).upper()
                
                # Get reference metadata
                refseq_name = ref_record.id
                start = ref_record.annotations.get('start', 0)
                size = ref_record.annotations.get('size', 0)
                
                # Get alignment length (number of columns)
                alignment_length = alignment.get_alignment_length()
                
                # Calculate mean pairwise identity to reference
                pairwise_identities = []
                if len(alignment) > 1:
                    for i in range(1, len(alignment)):
                        seq = str(alignment[i].seq).upper()
                        identity = self._calculate_pairwise_similarity(ref_seq, seq)
                        pairwise_identities.append(identity)
                
                mean_pairwise_identity = sum(pairwise_identities) / len(pairwise_identities) if pairwise_identities else 1.0
                
                # Calculate average gap fraction
                gap_fractions = []
                for record in alignment:
                    seq_str = str(record.seq)
                    gap_count = seq_str.count('-')
                    gap_fraction = gap_count / len(seq_str) if len(seq_str) > 0 else 0
                    gap_fractions.append(gap_fraction)
                
                avg_gap_fraction = sum(gap_fractions) / len(gap_fractions) if gap_fractions else 0.0
                
                # Write row
                out.write(f"{refseq_name}\t{start}\t{size}\t{alignment_length}\t"
                        f"{mean_pairwise_identity:.4f}\t{avg_gap_fraction:.4f}\n")
        
        print(f"Block description complete: {block_num} blocks processed")
        print(f"Output written to: {output_tsv_file}")

            
    def filter(self, output_maf_file: str, 
            maximum_average_gaps: float = 0.25,
            max_pairwise_similarity: float = 1.0,
            min_pairwise_similarity: float = 0.5,
            min_number_seqs: int = 3,
            block_length: Optional[int] = None):
        """
        Filter MAF alignment blocks based on various criteria.
        
        Args:
            output_maf_file: Path to output MAF file for filtered blocks
            maximum_average_gaps: Maximum average gap fraction allowed (default: 0.25)
            max_pairwise_similarity: Maximum pairwise similarity to reference (default: 1.0)
            min_pairwise_similarity: Minimum pairwise similarity to reference (default: 0.5)
            min_number_seqs: Minimum number of sequences required (default: 3)
            block_length: Required alignment length; None means no filter (default: None)
        """
        blocks_total = 0
        blocks_passed = 0
        
        with open(output_maf_file, 'w') as out:
            # Write MAF header
            out.write("##maf version=1\n\n")
            
            # Process each alignment block
            alignments = AlignIO.parse(self.maf_file, "maf")
            
            for alignment in alignments:
                blocks_total += 1
                
                # Check if block passes all filters
                if self._check_filters(alignment, maximum_average_gaps, 
                                    max_pairwise_similarity, min_pairwise_similarity,
                                    min_number_seqs, block_length):
                    blocks_passed += 1
                    
                    # Write the block to output
                    self._write_alignment_block(out, alignment)
        
        print(f"Filtering complete: {blocks_passed}/{blocks_total} blocks passed filters")
        print(f"Output written to: {output_maf_file}")

    def _check_filters(self, alignment, max_avg_gaps: float, max_sim: float, 
                    min_sim: float, min_seqs: int, block_len: Optional[int]) -> bool:
        """Check if alignment block passes all filter criteria."""
        
        # 1. Check minimum number of sequences
        if len(alignment) < min_seqs:
            return False
        
        # 2. Check block length if specified
        if block_len is not None and alignment.get_alignment_length() != block_len:
            return False
        
        # Find reference sequence (first one matching reference species)
        ref_record = None
        ref_idx = None
        for i, record in enumerate(alignment):
            if record.id.startswith(self.reference_species):
                ref_record = record
                ref_idx = i
                break
        
        if ref_record is None:
            return False  # No reference found
        
        # 3. Calculate average gap fraction
        gap_fractions = []
        for record in alignment:
            seq_str = str(record.seq)
            gap_count = seq_str.count('-')
            gap_fraction = gap_count / len(seq_str) if len(seq_str) > 0 else 0
            gap_fractions.append(gap_fraction)
        
        avg_gap_fraction = sum(gap_fractions) / len(gap_fractions)
        if avg_gap_fraction > max_avg_gaps:
            return False
        
        # 4. Calculate pairwise similarities to reference
        ref_seq = str(ref_record.seq).upper()
        
        for i, record in enumerate(alignment):
            if i == ref_idx:  # Skip self-comparison
                continue
                
            seq = str(record.seq).upper()
            similarity = self._calculate_pairwise_similarity(ref_seq, seq)
            
            if similarity > max_sim or similarity < min_sim:
                return False
        
        return True  # Passed all filters


    def _write_alignment_block(self, file_handle, alignment):
        """Write a single alignment block in MAF format."""
        # Get score from alignment annotations if available
        score = alignment.annotations.get('score', 0.0) if hasattr(alignment, 'annotations') else 0.0
        
        # Write 'a' line
        file_handle.write(f"a score={score}\n")
        
        # Write 's' lines for each sequence
        for record in alignment:
            # Get sequence info from annotations
            start = record.annotations.get('start', 0)
            size = record.annotations.get('size', 0)
            strand = record.annotations.get('strand', '+')
            src_size = record.annotations.get('srcSize', 0)
            
            # Write sequence line
            file_handle.write(f"s {record.id} {start} {size} {strand} {src_size} {str(record.seq)}\n")
        
        # Add blank line between blocks
        file_handle.write("\n")

    def windows(self, output_maf_file: str, length: int, step: int):
        """
        Apply sliding window approach to alignment blocks.
        
        Blocks longer than 'length' are cut into overlapping smaller blocks.
        Blocks already at or below 'length' are kept as-is.
        
        Args:
            output_maf_file: Path to output MAF file
            length: Length of resulting smaller alignment blocks
            step: Step size for overlapping windows
        """
        if step <= 0 or step > length:
            raise ValueError("Step must be positive and not exceed window length")
        
        blocks_total = 0
        windows_created = 0
        
        with open(output_maf_file, 'w') as out:
            # Write MAF header
            out.write("##maf version=1\n\n")
            
            # Process each alignment block
            alignments = AlignIO.parse(self.maf_file, "maf")
            
            for alignment in alignments:
                blocks_total += 1
                alignment_length = alignment.get_alignment_length()
                
                # If block is already small enough, keep as-is
                if alignment_length <= length:
                    windows_created += 1
                    self._write_alignment_block(out, alignment)
                else:
                    # Apply sliding window
                    window_starts = range(0, alignment_length - length + 1, step)
                    
                    # Handle last window if needed
                    if window_starts[-1] + length < alignment_length:
                        window_starts = list(window_starts)
                        # Add a final window that ends at the alignment end
                        final_start = alignment_length - length
                        if final_start > window_starts[-1]:
                            window_starts.append(final_start)
                    
                    for window_start in window_starts:
                        window_end = min(window_start + length, alignment_length)
                        window_alignment = self._extract_window(alignment, window_start, window_end)
                        
                        if window_alignment:
                            windows_created += 1
                            self._write_alignment_block(out, window_alignment)
        
        print(f"Windowing complete: {blocks_total} blocks processed")
        print(f"Created {windows_created} windows/blocks total")
        print(f"Output written to: {output_maf_file}")

    def _extract_window(self, alignment, window_start: int, window_end: int) -> MultipleSeqAlignment:
        """
        Extract a window from an alignment block.
        
        Args:
            alignment: Original alignment block
            window_start: Start column position (0-based)
            window_end: End column position (exclusive)
        
        Returns:
            New alignment for the window with updated coordinates
        """
        new_records = []
        
        # Find reference sequence for coordinate calculations
        ref_record = None
        ref_idx = None
        for i, record in enumerate(alignment):
            if record.id.startswith(self.reference_species):
                ref_record = record
                ref_idx = i
                break
        
        # Process each sequence in the alignment
        for i, record in enumerate(alignment):
            # Extract sequence window
            window_seq = str(record.seq)[window_start:window_end]
            
            # Calculate new start position
            if i == ref_idx and ref_record is not None:
                # For reference: count non-gap bases before window
                ref_seq_before = str(record.seq)[:window_start]
                bases_before = sum(1 for base in ref_seq_before if base != '-')
                new_start = record.annotations.get('start', 0) + bases_before
            else:
                # For non-reference sequences: need to calculate based on their gaps
                seq_before = str(record.seq)[:window_start]
                bases_before = sum(1 for base in seq_before if base != '-')
                new_start = record.annotations.get('start', 0) + bases_before
            
            # Count non-gap bases in window
            new_size = sum(1 for base in window_seq if base != '-')
            
            # Create new record with updated annotations
            new_annotations = record.annotations.copy() if hasattr(record, 'annotations') else {}
            new_annotations['start'] = new_start
            new_annotations['size'] = new_size
            
            new_record = SeqRecord(
                Seq(window_seq),
                id=record.id,
                description=record.description,
                annotations=new_annotations
            )
            new_records.append(new_record)
        
        # Create new alignment with updated annotations
        new_alignment = MultipleSeqAlignment(new_records)
        
        # Copy alignment-level annotations if present
        if hasattr(alignment, 'annotations'):
            new_alignment.annotations = alignment.annotations.copy()
        
        return new_alignment

        
def main():
    parser = OptionParser(usage=usage_statement, version="__version__")
    args = sys.argv

    if len(args) < 3:
        print(usage_statement)
        sys.exit()

    maffile = args[-1]

    mafobject = MAFObject(maffile, "", "")
    
    if len(args) < 2:
        print(usage_statement)
        sys.exit()

    if args[1] == "describe":
        mafobject.describe()

    elif args[1] == "blockdescribe":
        parser.add_option("-o", "--output", action="store", default="mafblocks.tsv", type="string", dest="output_tsv", help="Output TSV file (default: mafblocks.tsv)")
        options, args = parser.parse_args()

        mafobject.blockdescribe(options.output_tsv)

    elif args[1] == "merge":
        parser.add_option("-o", "--output", action="store", type="string", dest="out_file", default="output.maf", help="MAF file to write to. If empty, results alignments are redirected to output.maf.")
        parser.add_option("-s", "--species-consensus", action="store", type="float", default=0.75, dest="species_consensus", help="Minimal consensus between neighboring blocks for merging (Default: 0.75).")
        parser.add_option("-r", "--no-reference", action="store_false", default=True, dest="reference", help="Set this flag if the first sequence should NOT be considered as reference.")
        parser.add_option("-d", "--max-distance", action="store", default=0, type="int", dest="distance", help="Maximum distance between genomic coordinates of sequences for merging of neighboring blocks (Default: 0).")
        parser.add_option("-l", "--max-length", action="store", default=1000, type="int", dest="length", help="Merged alignment blocks will not be extended past this block length (Default: 1000).")
        parser.add_option("-g", "--max-gaps", action="store", default=0.5, type="float", dest="max_gaps", help="All sequences with a larger gap fraction than this value will be dropped (Default: 0.5).")
        parser.add_option("-m", "--min-seqs", action="store", default=2, type="int", dest="min_seqs", help="No merging will happen, if the blocks share this few or less sequences (Default: 2).")
        options, args = parser.parse_args()

        mafobject.merge(parser)

    elif args[1] == "select":
        parser.add_option("-b", "--bed", action="store", default="", type="string", dest="bed", help="The bed file with coordinated to extract.")
        parser.add_option("-r", "--right-flank", action="store", default=0, type="int", dest="right_flank", help="Add this many nucleotides to the right side of each coordinate to extract. (Default: 0)")
        parser.add_option("-l", "--left-flank", action="store", default=0, type="int", dest="left_flank", help="Add this many nucleotides to the left side of each coordinate to extract. (Default: 0)")
        parser.add_option("-o", "--output", action="store", default="ouput.maf", type="string", dest="out_maf", help="The MAF file to write to. If none, will direct to output.maf.")
        options, args = parser.parse_args()

        if len(options.bed) == 0:
            print("A bed file with coordinates is required.")
            sys.exit()

        mafobject.bed_file = options.bed
        mafobject.extract_blocks(options.out_maf)

    elif args[1] == "window":
        # output_maf_file: str, length: int, step: int
        parser.add_option("-l", "--length", action="store", default=120, type="int", dest="length", help="Target length of the cut alignment blocks (Default: 120).")
        parser.add_option("-s", "--step", action="store", default=40, type="int", dest="step", help="Step size for the sliding window (Default: 40)")
        parser.add_option("-o", "--output", action="store", default="ouput.maf", type="string", dest="out_maf", help="The MAF file to write to. If none, will direct to output.maf.")
        options, args = parser.parse_args()

        mafobject.windows(output_maf_file=options.out_maf, length=options.length, step=options.step)


    elif args[1] == "filter":
        parser.add_option("-g", "--max-gaps", action="store", default=0.25, type="float", dest="max_gaps", help="Maximum fraction of gaps. If higher, alignmnent block is discarded (Default: 0.25).")
        parser.add_option("-a", "--max-similarity", action="store", default=1.0, type="float", dest="max_similarity", help="Maximum pairwise similarity to the reference for a block to be kept (Default: 1.0).")
        parser.add_option("-m", "--min-similarity", action="store", default=0.5, type="float", dest="min_similarity", help="Minimum pairwise similarity to the reference for a block to be kept (Default: 0.5).")
        parser.add_option("-l", "--min-length", action="store", default=40, type="int", dest="min_length", help="Minimum block length. Shorter blocks are discarded (Default: 40).")
        parser.add_option("-n", "--min-seqs", action="store", default=3, type="int", dest="min_seqs", help="Minimum number of sequences to keep block (Default: 3)")
        parser.add_option("-o", "--output", action="store", default="ouput.maf", type="string", dest="out_maf", help="The MAF file to write to. If none, will direct to output.maf.")
        options, args = parser.parse_args()

        mafobject.filter(output_maf_file=options.out_maf,
                         maximum_average_gaps=options.max_gaps,
                         max_pairwise_similarity=options.max_similarity,
                         min_pairwise_similarity=options.min_similarity,
                         min_number_seqs=options.min_seqs,
                         block_length=options.min_length)

    elif args[1] == "toBed":
        parser.add_option("-o", "--output_bed", action="store", default="ouput.bed", type="string", dest="output_bed", help="The bed file to write to. If none, will direct to output.bed.")
        options, args = parser.parse_args()
        mafobject.toBed(options.output_bed)

    else:
        print(usage_statement)
        sys.exit()
    

if __name__=='__main__':
    main()
