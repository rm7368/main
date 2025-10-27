#!/usr/bin/env python3
"""
RNA-FISH Probe Designer

Design and validate RNA-FISH probe pairs with specificity analysis.
Supports both direct sequence input and NCBI RefSeq IDs with barcode system.
"""

from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import re
import subprocess
import time
import math
import pandas as pd
import argparse
import sys
import os
from collections import Counter

# Global parameters
DEFAULT_POLY_AT = 7
DEFAULT_POLY_CG = 6
DEFAULT_POLY_TRIPLET = 3
DEFAULT_FORBIDDEN_TRINUCS = ["gtc", "gac", "gcc", "ggg", "ccc"]

# Amplifier and spacer combinations
AMPLIFIER_CONFIGS = {
    'B1': {
        'amplifier': 'GAGGAGGGCAGCAAACGGGAAGAGTCTTCCTTTACG',
        'spacer1': 'aa',  # for 5' probe
        'spacer2': 'ta'   # for 3' probe
    },
    'B2': {
        'amplifier': 'CCTCGTAAATCCTCATCAATCATCCAGTAAACCGCC',
        'spacer1': 'aa',
        'spacer2': 'aa'
    },
    'B3': {
        'amplifier': 'GTCCCTGCCTCTATATCTCCACTCAACTTTAACCCG',
        'spacer1': 'tt',
        'spacer2': 'tt'
    },
    'B4': {
        'amplifier': "CCTCAACCTACCTCCAACTCTCACCATATTCGCTTC",
        'spacer1': 'aa',
        'spacer2': 'at'
    },
    'B5': {
        'amplifier': "CTCACTCCCAATCTCTATCTACCCTACAAATCCAAT",
        'spacer1': 'aa',
        'spacer2': 'aa'
    }
}

def validate_input_sequence(sequence):
    """Validate input DNA sequence"""
    if not isinstance(sequence, str):
        raise TypeError("Sequence must be a string")
    
    if len(sequence) < 100:
        raise ValueError(f"Sequence too short ({len(sequence)} bp). Need at least 100 bp for probe design")
    
    # Check for valid nucleotides
    valid_nucs = set('ATCGN')
    sequence_nucs = set(sequence.upper())
    invalid_nucs = sequence_nucs - valid_nucs
    
    if invalid_nucs:
        raise ValueError(f"Invalid nucleotides found: {invalid_nucs}. Only A, T, C, G, N allowed")
    
    return sequence.upper()

def fetch_sequence_from_refseq(refseq_id, email):
    """Fetch sequence from NCBI using RefSeq ID"""
    if not isinstance(refseq_id, str):
        raise TypeError("RefSeq ID must be a string")
    
    if not email or "@" not in email:
        raise ValueError("Valid email address required for NCBI queries")
    
    Entrez.email = email
    
    try:
        print(f"Fetching sequence for {refseq_id}...")
        handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        
        sequence = str(record.seq).upper()
        print(f"Retrieved sequence: {len(sequence)} bp")
        return validate_input_sequence(sequence), record.description
        
    except Exception as e:
        raise ValueError(f"Could not fetch sequence for {refseq_id}: {e}")

def get_gene_name_from_description(description):
    """Extract gene name from NCBI description"""
    if not isinstance(description, str):
        return "Unknown_Gene"
    
    # Common patterns in NCBI descriptions
    if "(" in description and ")" in description:
        start = description.find("(") + 1
        end = description.find(")", start)
        gene_name = description[start:end]
        # Clean up gene name
        return re.sub(r'[^\w\-_]', '_', gene_name)[:20]
    else:
        # Fallback: use first word
        first_word = description.split()[0] if description.split() else "Unknown_Gene"
        return re.sub(r'[^\w\-_]', '_', first_word)[:20]

def validate_barcode_config(barcode_id):
    """Validate barcode ID and return config"""
    if barcode_id not in AMPLIFIER_CONFIGS:
        raise ValueError(f"Invalid barcode {barcode_id}. Must be one of: {list(AMPLIFIER_CONFIGS.keys())}")
    
    config = AMPLIFIER_CONFIGS[barcode_id]
    
    # Validate amplifier length
    if len(config['amplifier']) != 36:
        raise ValueError(f"Amplifier for {barcode_id} must be 36 nucleotides")
    
    return config

def everyprobe(sequence, amplifier, spacer1, spacer2):
    """Generate all possible probe pairs from RNA sequence"""
    if not isinstance(sequence, str):
        raise TypeError("Sequence must be a string")
    
    if len(sequence) < 52:
        raise ValueError(f"Sequence too short ({len(sequence)} bp). Need at least 52 bp")
    
    i = 0
    seqlist = []
    while 0 <= i <= (len(sequence) - 52):
        p1 = Seq(sequence[i:i+25])
        probe1binding = str(p1.reverse_complement())
        p2 = Seq(sequence[i+27:i+52])
        probe2binding = str(p2.reverse_complement())
        
        # Use specific spacers for each probe
        current = (amplifier[:18] + spacer1 + probe1binding, probe2binding + spacer2 + amplifier[18:])
        seqlist.append(current)
        i += 1
    return seqlist

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    if not isinstance(sequence, str) or not sequence:
        return 0.0
    gc_count = sequence.lower().count("g") + sequence.lower().count("c")
    return (gc_count / len(sequence)) * 100

def calculate_tm(sequence):
    """Calculate melting temperature using accurate formula"""
    if not isinstance(sequence, str) or not sequence:
        return 0.0
    
    sequence = sequence.lower()
    g_count = sequence.count("g")
    c_count = sequence.count("c")
    total_length = len(sequence)
    gc_count = g_count + c_count
    tm = 64.9 + 41 * ((gc_count - 16.4) / total_length)
    return tm

def validate_probe(probe, gc_min=45, gc_max=60, tm_min=60, tm_max=75, 
                  amplifier=None, spacer1=None, spacer2=None,
                  polyAT=DEFAULT_POLY_AT, polyCG=DEFAULT_POLY_CG, 
                  polytriplet=DEFAULT_POLY_TRIPLET, badtris=None):
    """Validate individual probe against criteria"""
    if badtris is None:
        badtris = DEFAULT_FORBIDDEN_TRINUCS
    
    if not isinstance(probe, str):
        raise TypeError("Probe must be a string")
    
    results = {}
    
    # Count mononucleotide repeats
    mono_repeats = (len(re.findall(f'A{{{polyAT},}}', probe)) +
                    len(re.findall(f'T{{{polyAT},}}', probe)) +
                    len(re.findall(f'C{{{polyCG},}}', probe)) +
                    len(re.findall(f'G{{{polyCG},}}', probe)))
    
    # Count trinucleotide repeats
    tri_repeats = 0
    for tri in badtris:
        pattern = f'({tri}){{{polytriplet},}}'
        tri_repeats += len(re.findall(pattern, probe.lower()))
    
    # Calculate binding length
    total_binding_length = 0
    if amplifier and spacer1 and probe[:18] == amplifier[:18]:
        # Check for spacer1 presence
        spacer_start = 18
        spacer_end = 18 + len(spacer1)
        if probe[spacer_start:spacer_end] == spacer1:
            total_binding_length = len(probe[spacer_end:])
    elif amplifier and spacer2:
        # Check for spacer2 at end
        amp_end_len = len(amplifier[18:])
        spacer_len = len(spacer2)
        if (len(probe) >= amp_end_len + spacer_len and 
            probe[-(amp_end_len):] == amplifier[18:] and
            probe[-(amp_end_len + spacer_len):-(amp_end_len)] == spacer2):
            total_binding_length = len(probe[:-(amp_end_len + spacer_len)])
    
    total_penalty_score = mono_repeats + tri_repeats
    gc_content = calculate_gc_content(probe)
    tm = calculate_tm(probe)
    gc_pass = gc_min <= gc_content <= gc_max
    tm_pass = tm_min <= tm <= tm_max
    passes = gc_pass and tm_pass
    
    results = {
        "mono_repeats": mono_repeats,
        "tri_repeats": tri_repeats,
        "total_binding_length": total_binding_length,
        "total_penalty_score": total_penalty_score,
        "gc_content": gc_content,
        "tm": tm,
        "passes_filters": passes
    }
    return results

def filter_passing_probes(probe_list, **kwargs):
    """Filter probe pairs where both probes pass validation"""
    if not isinstance(probe_list, list):
        raise TypeError("Probe list must be a list")
    
    passing_probes = []
    for probe_pair in probe_list:
        if not isinstance(probe_pair, (tuple, list)) or len(probe_pair) != 2:
            continue
        
        result_p1 = validate_probe(probe_pair[0], **kwargs)
        result_p2 = validate_probe(probe_pair[1], **kwargs)
        if result_p1["passes_filters"] and result_p2["passes_filters"]:
            passing_probes.append(probe_pair)
    return passing_probes

def penalty_filter(probe_list, top_n=100, penalty_threshold=3, **kwargs):
    """Filter probes by penalty score with distribution awareness"""
    if not isinstance(probe_list, list):
        raise TypeError("Probe list must be a list")
    
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    
    scored_probes = []
    for probe_pair in probe_list:
        penalty = ((validate_probe(probe_pair[0], **kwargs))["total_penalty_score"] + 
                   (validate_probe(probe_pair[1], **kwargs))["total_penalty_score"]) / 2
        if penalty <= penalty_threshold:
            scored_probes.append(probe_pair)
    
    if len(scored_probes) <= top_n:
        return scored_probes
    
    step_size = len(scored_probes) / top_n
    selected_indices = [int(i * step_size) for i in range(top_n)]
    return [scored_probes[i] for i in selected_indices]

def get_target_rna_region(probe_index, sequence):
    """Get target RNA region for specific probe position"""
    if not isinstance(probe_index, int) or probe_index < 0:
        raise ValueError("Probe index must be a non-negative integer")
    
    target_rna = sequence.replace('T', 'U').upper()
    region_start = probe_index
    region_end = probe_index + 52
    
    if region_end > len(target_rna):
        raise ValueError(f"Probe index {probe_index} exceeds sequence length")
    
    return target_rna[region_start:region_end]

def extract_binding_regions(probe_pair, amplifier, spacer1, spacer2):
    """Extract just the binding regions from full probe sequences"""
    if not isinstance(probe_pair, (tuple, list)) or len(probe_pair) != 2:
        raise TypeError("Probe pair must be a tuple or list of 2 sequences")
    
    probe1, probe2 = probe_pair
    
    # Extract probe 1 binding region (has spacer1)
    if probe1[:18] == amplifier[:18]:
        spacer_end = 18 + len(spacer1)
        probe1_binding = probe1[spacer_end:]
    else:
        probe1_binding = probe1[:25]  # Fallback
        
    # Extract probe 2 binding region (has spacer2)
    amp_end_len = len(amplifier[18:])
    spacer_len = len(spacer2)
    if (len(probe2) >= amp_end_len + spacer_len and 
        probe2[-(amp_end_len):] == amplifier[18:]):
        probe2_binding = probe2[:-(amp_end_len + spacer_len)]
    else:
        probe2_binding = probe2[:25]  # Fallback
    
    return probe1_binding, probe2_binding

def parse_energy_from_output(stdout):
    """Parse energy from RNAcofold output"""
    if not isinstance(stdout, str):
        return None
    
    lines = stdout.strip().split('\n')
    for line in lines:
        if '(' in line and ')' in line:
            energy_match = re.search(r'\(\s*([-+]?\d*\.?\d+)\s*\)', line)
            if energy_match:
                return float(energy_match.group(1))
    return None

def calculate_individual_binding_energies(probe_pair, target_rna_region, amplifier, spacer1, spacer2):
    """Calculate binding energy for each probe individually"""
    probe1_binding, probe2_binding = extract_binding_regions(probe_pair, amplifier, spacer1, spacer2)
    
    results = {}
    
    # Analyze probe 1 binding to target region 1 (first 25bp)
    target_region1 = target_rna_region[:25]
    sequences1 = f"{target_region1}&{probe1_binding}"
    
    # Analyze probe 2 binding to target region 2 (last 25bp)
    target_region2 = target_rna_region[27:52]
    sequences2 = f"{target_region2}&{probe2_binding}"
    
    # Calculate binding energy for probe 1
    try:
        result1 = subprocess.run(['RNAcofold', '--noPS'], 
                                input=sequences1, 
                                text=True, 
                                capture_output=True, 
                                timeout=30)
        
        if result1.returncode == 0:
            energy1 = parse_energy_from_output(result1.stdout)
            results['probe1_energy'] = energy1
        else:
            results['probe1_energy'] = None
            
    except Exception:
        results['probe1_energy'] = None
    
    # Calculate binding energy for probe 2
    try:
        result2 = subprocess.run(['RNAcofold', '--noPS'], 
                                input=sequences2, 
                                text=True, 
                                capture_output=True, 
                                timeout=30)
        
        if result2.returncode == 0:
            energy2 = parse_energy_from_output(result2.stdout)
            results['probe2_energy'] = energy2
        else:
            results['probe2_energy'] = None
            
    except Exception:
        results['probe2_energy'] = None
    
    # Combined energy
    if results['probe1_energy'] is not None and results['probe2_energy'] is not None:
        results['combined_energy'] = results['probe1_energy'] + results['probe2_energy']
    else:
        results['combined_energy'] = None
    
    return results

def run_full_vienna_analysis(probe_candidates, sequence, amplifier, spacer1, spacer2):
    """Run ViennaRNA analysis on all probe candidates"""
    if not isinstance(probe_candidates, list):
        raise TypeError("Probe candidates must be a list")
    
    # Generate all probes for index lookup
    all_probes = everyprobe(sequence, amplifier, spacer1, spacer2)
    
    vienna_results = []
    print(f"Running ViennaRNA analysis on {len(probe_candidates)} probe pairs...")
    
    for i, probe_pair in enumerate(probe_candidates):
        try:
            original_index = all_probes.index(probe_pair)
            target_region = get_target_rna_region(original_index, sequence)
            
            binding_result = calculate_individual_binding_energies(probe_pair, target_region, amplifier, spacer1, spacer2)
            
            if binding_result['combined_energy'] is not None:
                result = {
                    'probe_pair': probe_pair,
                    'original_position': original_index,
                    'probe1_energy': binding_result['probe1_energy'],
                    'probe2_energy': binding_result['probe2_energy'],
                    'combined_energy': binding_result['combined_energy']
                }
                vienna_results.append(result)
            
            if (i + 1) % 10 == 0:
                print(f"Processed {i + 1}/{len(probe_candidates)} probes...")
        
        except Exception as e:
            print(f"Warning: Error processing probe {i + 1}: {e}")
            continue
    
    return vienna_results

def fixed_select_distributed_candidates(vienna_results, target_count=45, sequence_length=None):
    """Select candidates with optimal distribution, avoiding duplicates"""
    if not isinstance(vienna_results, list):
        raise TypeError("Vienna results must be a list")
    
    if target_count <= 0:
        raise ValueError("target_count must be positive")
    
    transcript_length = sequence_length or len(vienna_results)
    
    # Create bins across the transcript
    num_bins = target_count
    bin_size = transcript_length / num_bins
    
    # Group probes by bins
    binned_probes = {}
    for result in vienna_results:
        bin_idx = int(result['original_position'] / bin_size)
        if bin_idx not in binned_probes:
            binned_probes[bin_idx] = []
        binned_probes[bin_idx].append(result)
    
    # Select best probe from each bin, avoiding duplicates
    selected_probes = []
    used_probe_ids = set()
    
    for bin_idx in range(num_bins):
        candidates = []
        
        # Look in current bin and adjacent bins
        for search_bin in [bin_idx-1, bin_idx, bin_idx+1]:
            if search_bin in binned_probes:
                candidates.extend(binned_probes[search_bin])
        
        if candidates:
            # Sort by binding energy (best first)
            candidates.sort(key=lambda x: x['combined_energy'])
            
            # Find the best candidate not already selected
            for candidate in candidates:
                probe_id = (candidate['probe_pair'][0], candidate['probe_pair'][1])
                if probe_id not in used_probe_ids:
                    selected_probes.append(candidate)
                    used_probe_ids.add(probe_id)
                    break
    
    # Fill remaining slots with best available candidates
    while len(selected_probes) < target_count and len(selected_probes) < len(vienna_results):
        remaining = [r for r in vienna_results 
                    if (r['probe_pair'][0], r['probe_pair'][1]) not in used_probe_ids]
        
        if remaining:
            remaining.sort(key=lambda x: x['combined_energy'])
            best_remaining = remaining[0]
            selected_probes.append(best_remaining)
            probe_id = (best_remaining['probe_pair'][0], best_remaining['probe_pair'][1])
            used_probe_ids.add(probe_id)
        else:
            break
    
    return selected_probes[:target_count]

def quick_specificity_analysis(blast1_records, blast2_records, candidate, target_gene="Unknown"):
    """Quick analysis focused on counting concerning hits"""
    
    def count_concerning_hits(blast_records):
        concerning = 0
        total = 0
        
        if blast_records and blast_records[0].alignments:
            for alignment in blast_records[0].alignments:
                total += 1
                # Skip if it's the target gene
                title_lower = alignment.title.lower()
                if target_gene.lower() in title_lower:
                    continue
                    
                for hsp in alignment.hsps:
                    identity = hsp.identities / hsp.align_length
                    if identity > 0.8:  # >80% identity to non-target
                        concerning += 1
                        break
        
        return concerning, total
    
    p1_concerning, p1_total = count_concerning_hits(blast1_records)
    p2_concerning, p2_total = count_concerning_hits(blast2_records)
    
    total_concerning = p1_concerning + p2_concerning
    specificity_score = max(0, 100 - (total_concerning * 20))
    
    return {
        'candidate': candidate,
        'concerning_hits': total_concerning,
        'total_hits': p1_total + p2_total,
        'specificity_score': specificity_score
    }

def run_blast_analysis(blast_candidates, gene_name="Unknown", skip_blast=False, amplifier=None, spacer1=None, spacer2=None):
    """Run BLAST analysis on binding regions only"""
    if skip_blast:
        print("Skipping BLAST analysis...")
        return [{
            'candidate': candidate,
            'full_probe_pair': candidate['probe_pair'],
            'specificity_score': 100,
            'concerning_hits': 0
        } for candidate in blast_candidates]
    
    if not isinstance(blast_candidates, list):
        raise TypeError("BLAST candidates must be a list")
    
    if not amplifier or not spacer1 or not spacer2:
        raise ValueError("Amplifier and spacer sequences required for BLAST analysis")
    
    results = []
    print(f"Running BLAST analysis on {len(blast_candidates)} candidates...")
    
    for i, candidate in enumerate(blast_candidates):
        try:
            binding_probe1, binding_probe2 = extract_binding_regions(candidate['probe_pair'], amplifier, spacer1, spacer2)
            
            if (i + 1) % 5 == 0:
                print(f"Progress: {i+1}/{len(blast_candidates)} candidates...")
            
            # BLAST only the binding regions
            result1 = NCBIWWW.qblast("blastn", "refseq_rna", binding_probe1, 
                                   hitlist_size=10, expect=0.1, word_size=7,
                                   entrez_query="Mus musculus[Organism]")
            
            result2 = NCBIWWW.qblast("blastn", "refseq_rna", binding_probe2,
                                   hitlist_size=10, expect=0.1, word_size=7,
                                   entrez_query="Mus musculus[Organism]")
            
            blast1_records = list(NCBIXML.parse(result1))
            blast2_records = list(NCBIXML.parse(result2))
            
            analysis = quick_specificity_analysis(blast1_records, blast2_records, candidate, gene_name)
            analysis['full_probe_pair'] = candidate['probe_pair']
            analysis['binding_regions'] = (binding_probe1, binding_probe2)
            
            results.append(analysis)
            time.sleep(1)
            
        except Exception as e:
            print(f"Warning: Error with BLAST candidate {i+1}: {e}")
            # Add placeholder with poor score
            results.append({
                'candidate': candidate,
                'full_probe_pair': candidate['probe_pair'],
                'specificity_score': 0,
                'concerning_hits': 999,
                'error': str(e)
            })
    
    return results

def create_combined_idt_opool_excel(all_gene_results, output_file):
    """Create combined IDT oPool format Excel file for multiple genes"""
    if not isinstance(all_gene_results, list):
        raise TypeError("Gene results must be a list")
    
    if not all_gene_results:
        raise ValueError("No gene results to save")
    
    pool_data = []
    
    for gene_result in all_gene_results:
        gene_name = gene_result['gene_name']
        final_probes = gene_result['probes']
        
        for i, result in enumerate(final_probes):
            # Handle different result formats
            if 'full_probe_pair' in result:
                probe1, probe2 = result['full_probe_pair']
            elif 'probe_pair' in result:
                probe1, probe2 = result['probe_pair']
            else:
                print(f"Warning: Unexpected result format for {gene_name} probe {i+1}")
                continue
            
            pool_data.append({'Pool Name': gene_name, 'Sequence': probe1})
            pool_data.append({'Pool Name': gene_name, 'Sequence': probe2})
    
    if not pool_data:
        raise ValueError("No valid probe data to save")
    
    df = pd.DataFrame(pool_data)
    
    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='oPool_Order', index=False)
        
        print(f"Combined IDT oPool order file saved: {output_file}")
        print(f"Total genes: {len(all_gene_results)}")
        print(f"Total sequences: {len(pool_data)}")
        
        # Summary by gene
        gene_counts = df['Pool Name'].value_counts()
        print("\nSequences per gene:")
        for gene, count in gene_counts.items():
            print(f"  {gene}: {count} sequences ({count//2} probe pairs)")
        
        return output_file
    
    except Exception as e:
        raise RuntimeError(f"Failed to save Excel file: {e}")

def parse_gene_file(gene_file):
    """Parse file containing RefSeq,Barcode pairs"""
    if not os.path.isfile(gene_file):
        raise FileNotFoundError(f"Gene file not found: {gene_file}")
    
    gene_data = []
    
    with open(gene_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse RefSeq,Barcode
            parts = line.split(',')
            if len(parts) != 2:
                raise ValueError(f"Line {line_num}: Expected format 'RefSeq,Barcode', got: {line}")
            
            refseq_id = parts[0].strip()
            barcode_id = parts[1].strip()
            
            if not refseq_id:
                raise ValueError(f"Line {line_num}: RefSeq ID cannot be empty")
            
            if not barcode_id:
                raise ValueError(f"Line {line_num}: Barcode ID cannot be empty")
            
            # Validate barcode and get config
            try:
                config = validate_barcode_config(barcode_id)
            except ValueError as e:
                raise ValueError(f"Line {line_num}: {e}")
            
            gene_data.append({
                'refseq_id': refseq_id,
                'barcode_id': barcode_id,
                'amplifier': config['amplifier'],
                'spacer1': config['spacer1'],
                'spacer2': config['spacer2']
            })
    
    if not gene_data:
        raise ValueError("No valid gene data found in file")
    
    return gene_data

def process_multiple_genes(gene_data, email, **kwargs):
    """Process multiple genes with their individual amplifiers and spacers"""
    all_gene_results = []
    
    for i, gene_info in enumerate(gene_data):
        try:
            refseq_id = gene_info['refseq_id']
            barcode_id = gene_info['barcode_id']
            amplifier = gene_info['amplifier']
            spacer1 = gene_info['spacer1']
            spacer2 = gene_info['spacer2']
            
            print(f"\n{'='*60}")
            print(f"Processing gene {i+1}/{len(gene_data)}: {refseq_id}")
            print(f"Barcode: {barcode_id}")
            print(f"Amplifier: {amplifier}")
            print(f"Spacers: {spacer1}, {spacer2}")
            print(f"{'='*60}")
            
            # Fetch sequence and gene name
            sequence, description = fetch_sequence_from_refseq(refseq_id, email)
            gene_name = get_gene_name_from_description(description)
            
            print(f"Gene name: {gene_name}")
            
            # Run probe design pipeline with this gene's barcode config
            final_probes = design_probes_pipeline(
                sequence=sequence,
                gene_name=gene_name,
                amplifier=amplifier,
                spacer1=spacer1,
                spacer2=spacer2,
                **kwargs
            )
            
            all_gene_results.append({
                'refseq_id': refseq_id,
                'gene_name': gene_name,
                'barcode_id': barcode_id,
                'amplifier': amplifier,
                'spacer1': spacer1,
                'spacer2': spacer2,
                'probes': final_probes,
                'sequence_length': len(sequence)
            })
            
            print(f"Completed {gene_name}: {len(final_probes)} probe pairs selected")
            
        except Exception as e:
            print(f"Error processing {gene_info['refseq_id']}: {e}")
            print("Continuing with next gene...")
            continue
    
    return all_gene_results

def create_idt_opool_excel(final_probes, gene_name, output_file):
    """Create IDT oPool format Excel file for single gene"""
    if not isinstance(final_probes, list):
        raise TypeError("Final probes must be a list")
    
    if not isinstance(gene_name, str) or not gene_name:
        raise ValueError("Gene name must be a non-empty string")
    
    pool_data = []
    
    for i, result in enumerate(final_probes):
        # Handle different result formats
        if 'full_probe_pair' in result:
            probe1, probe2 = result['full_probe_pair']
        elif 'probe_pair' in result:
            probe1, probe2 = result['probe_pair']
        else:
            print(f"Warning: Unexpected result format for probe {i+1}")
            continue
        
        pool_data.append({'Pool Name': gene_name, 'Sequence': probe1})
        pool_data.append({'Pool Name': gene_name, 'Sequence': probe2})
    
    if not pool_data:
        raise ValueError("No valid probe data to save")
    
    df = pd.DataFrame(pool_data)
    
    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='oPool_Order', index=False)
        
        print(f"IDT oPool order file saved: {output_file}")
        print(f"Total sequences: {len(pool_data)}")
        return output_file
    
    except Exception as e:
        raise RuntimeError(f"Failed to save Excel file: {e}")

def design_probes_pipeline(sequence, gene_name, amplifier, spacer1, spacer2, final_count=30, 
                          gc_min=45, gc_max=60, tm_min=60, tm_max=75, skip_blast=False):
    """Main pipeline function"""
    try:
        # Validate inputs
        sequence = validate_input_sequence(sequence)
        
        if not isinstance(gene_name, str) or not gene_name:
            raise ValueError("Gene name must be a non-empty string")
        
        if final_count <= 0:
            raise ValueError("final_count must be positive")
        
        print("Starting probe design pipeline...")
        
        # Step 1: Generate all possible probes
        print("Step 1: Generating probe pairs...")
        all_probes = everyprobe(sequence, amplifier, spacer1, spacer2)
        print(f"Generated {len(all_probes)} probe pairs")
        
        # Step 2: GC/Tm filtering
        print("Step 2: Applying GC/Tm filters...")
        filtered_probes = filter_passing_probes(all_probes, gc_min=gc_min, gc_max=gc_max, 
                                               tm_min=tm_min, tm_max=tm_max, 
                                               amplifier=amplifier, spacer1=spacer1, spacer2=spacer2)
        print(f"Passed GC/Tm filters: {len(filtered_probes)} probe pairs")
        
        if len(filtered_probes) == 0:
            raise ValueError("No probes passed GC/Tm filters. Consider relaxing parameters.")
        
        # Step 3: Penalty filtering
        print("Step 3: Penalty filtering...")
        vienna_candidates = penalty_filter(filtered_probes, amplifier=amplifier, 
                                         spacer1=spacer1, spacer2=spacer2)
        print(f"Selected for ViennaRNA: {len(vienna_candidates)} probe pairs")
        
        # Step 4: ViennaRNA analysis
        print("Step 4: ViennaRNA binding energy analysis...")
        vienna_results = run_full_vienna_analysis(vienna_candidates, sequence, amplifier, spacer1, spacer2)
        print(f"ViennaRNA analysis complete: {len(vienna_results)} probe pairs")
        
        if len(vienna_results) == 0:
            raise ValueError("No probes survived ViennaRNA analysis")
        
        # Step 5: Coverage distribution
        print("Step 5: Selecting distributed candidates...")
        blast_candidates = fixed_select_distributed_candidates(vienna_results, 
                                                              target_count=min(final_count*2, len(vienna_results)),
                                                              sequence_length=len(sequence))
        print(f"Selected for BLAST: {len(blast_candidates)} probe pairs")
        
        # Step 6: BLAST analysis (optional)
        print("Step 6: BLAST specificity analysis...")
        blast_results = run_blast_analysis(blast_candidates, gene_name, skip_blast, 
                                         amplifier, spacer1, spacer2)
        
        # Step 7: Final selection
        print("Step 7: Final probe selection...")
        # Sort by specificity score, then by binding energy
        sorted_results = sorted(blast_results, 
                               key=lambda x: (x.get('specificity_score', 0), 
                                             -x['candidate'].get('combined_energy', 0)), 
                               reverse=True)
        
        final_probes = sorted_results[:final_count]
        
        print(f"Pipeline complete! Selected {len(final_probes)} final probe pairs")
        return final_probes
        
    except Exception as e:
        print(f"Pipeline error: {e}")
        raise

def main():
    """Main function for command line execution"""
    parser = argparse.ArgumentParser(
        description="Design RNA-FISH probes with specificity validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python probe_designer.py --sequence ATCGATCG... --gene-name GKN3 --barcode B1
  python probe_designer.py --refseq NM_001346162 --email your@email.com --barcode B2
  python probe_designer.py --gene-file genes.txt --email your@email.com

Gene file format (genes.txt):
  NM_010275,B1
  NM_007393,B2
  # Comments start with #

Available barcodes: B1, B2, B3, B4, B5
        """
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--sequence', type=str, help='Input DNA sequence (5\' to 3\')')
    input_group.add_argument('--refseq', type=str, help='NCBI RefSeq ID (e.g., NM_001346162)')
    input_group.add_argument('--gene-file', type=str, help='File with RefSeq,Barcode pairs (one per line)')
    
    # Required/Optional parameters
    parser.add_argument('--gene-name', type=str, help='Gene name for output (auto-detected if using RefSeq)')
    parser.add_argument('--email', type=str, help='Email for NCBI queries (required for RefSeq)')
    parser.add_argument('--output', type=str, default="probe_order.xlsx", help='Output Excel file (default: probe_order.xlsx)')
    parser.add_argument('--final-count', type=int, default=30, help='Number of final probe pairs (default: 30)')
    
    # Validation parameters
    parser.add_argument('--gc-min', type=float, default=45, help='Minimum GC content %% (default: 45)')
    parser.add_argument('--gc-max', type=float, default=60, help='Maximum GC content %% (default: 60)')
    parser.add_argument('--tm-min', type=float, default=60, help='Minimum Tm °C (default: 60)')
    parser.add_argument('--tm-max', type=float, default=75, help='Maximum Tm °C (default: 75)')
    parser.add_argument('--skip-blast', action='store_true', help='Skip BLAST analysis (faster)')
    
    # Barcode for single gene mode
    parser.add_argument('--barcode', type=str, choices=['B1', 'B2', 'B3', 'B4', 'B5'], 
                       help='Barcode ID (B1-B5) for single gene')
    
    args = parser.parse_args()
    
    # Validate arguments
    try:
        if args.refseq and not args.email:
            parser.error("--email is required when using --refseq")
        
        if args.refseq and not args.barcode:
            parser.error("--barcode is required when using --refseq")
        
        if args.gene_file and not args.email:
            parser.error("--email is required when using --gene-file")
        
        if args.sequence and not args.barcode:
            parser.error("--barcode is required when using --sequence")
        
        if args.gc_min >= args.gc_max:
            parser.error("--gc-min must be less than --gc-max")
        
        if args.tm_min >= args.tm_max:
            parser.error("--tm-min must be less than --tm-max")
        
        if args.final_count <= 0:
            parser.error("--final-count must be positive")
        
        # Handle different input types
        if args.sequence:
            # Single sequence with barcode
            sequence = args.sequence
            config = validate_barcode_config(args.barcode)
            gene_name = args.gene_name or "Unknown_Gene"
            print(f"Using provided sequence: {len(sequence)} bp")
            print(f"Using barcode: {args.barcode}")
            
            final_probes = design_probes_pipeline(
                sequence=sequence,
                gene_name=gene_name,
                amplifier=config['amplifier'],
                spacer1=config['spacer1'],
                spacer2=config['spacer2'],
                final_count=args.final_count,
                gc_min=args.gc_min,
                gc_max=args.gc_max,
                tm_min=args.tm_min,
                tm_max=args.tm_max,
                skip_blast=args.skip_blast
            )
            
            # Single gene - use original function
            output_file = create_idt_opool_excel(final_probes, gene_name, args.output)
            
        elif args.refseq:
            # Single RefSeq with barcode
            sequence, description = fetch_sequence_from_refseq(args.refseq, args.email)
            config = validate_barcode_config(args.barcode)
            gene_name = args.gene_name or get_gene_name_from_description(description)
            print(f"Gene name: {gene_name}")
            print(f"Using barcode: {args.barcode}")
            
            final_probes = design_probes_pipeline(
                sequence=sequence,
                gene_name=gene_name,
                amplifier=config['amplifier'],
                spacer1=config['spacer1'],
                spacer2=config['spacer2'],
                final_count=args.final_count,
                gc_min=args.gc_min,
                gc_max=args.gc_max,
                tm_min=args.tm_min,
                tm_max=args.tm_max,
                skip_blast=args.skip_blast
            )
            
            # Single gene - use original function  
            output_file = create_idt_opool_excel(final_probes, gene_name, args.output)
            
        elif args.gene_file:
            # Multiple genes from file
            gene_data = parse_gene_file(args.gene_file)
            print(f"Processing {len(gene_data)} genes from file: {args.gene_file}")
            
            # Process all genes
            all_gene_results = process_multiple_genes(
                gene_data=gene_data,
                email=args.email,
                final_count=args.final_count,
                gc_min=args.gc_min,
                gc_max=args.gc_max,
                tm_min=args.tm_min,
                tm_max=args.tm_max,
                skip_blast=args.skip_blast
            )
            
            if not all_gene_results:
                raise ValueError("No genes were successfully processed")
            
            # Create combined Excel file
            output_file = create_combined_idt_opool_excel(all_gene_results, args.output)
        
        print(f"\nProbe design complete!")
        print(f"Output file: {output_file}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

# Allow both command line and module import
if __name__ == "__main__":
    main()